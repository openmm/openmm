#include <pybind11/pybind11.h>
#include <torch/extension.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

#define private public
#include "openmm/Context.h"
#undef private

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "CudaArray.h"
#include "CudaContext.h"
#include "CudaPlatform.h"

#include <cctype>
#include <cstdint>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace py = pybind11;
using OpenMM::Context;
using OpenMM::ContextImpl;
using OpenMM::CudaArray;
using OpenMM::CudaContext;
using OpenMM::CudaPlatform;

extern "C" void openmm_bridge_extract_xyz(
        std::uintptr_t src_ptr,
        std::uintptr_t correction_ptr,
        int element_size,
        double* dst,
        int n,
        double extract_scale,
        cudaStream_t stream);

extern "C" void openmm_bridge_scatter_xyz(
        std::uintptr_t dst_ptr,
        std::uintptr_t correction_ptr,
        int element_size,
        const double* src,
        int n,
        double scatter_scale,
        cudaStream_t stream);

namespace {

std::uintptr_t parse_hex_pointer_from_repr(const std::string& repr) {
    const std::string marker = "0x";
    const std::size_t start = repr.find(marker);
    if (start == std::string::npos)
        throw std::runtime_error("could not find a SWIG pointer address in object repr: " + repr);
    std::size_t end = start + marker.size();
    while (end < repr.size() && std::isxdigit(static_cast<unsigned char>(repr[end])))
        ++end;
    std::uintptr_t value = 0;
    std::stringstream stream;
    stream << std::hex << repr.substr(start + marker.size(), end - start - marker.size());
    stream >> value;
    if (value == 0)
        throw std::runtime_error("parsed a null SWIG pointer from object repr: " + repr);
    return value;
}

std::uintptr_t pointer_from_swig_object(const py::object& object) {
    py::object raw = py::hasattr(object, "this") ? object.attr("this") : object;
    try {
        py::object as_int = py::module_::import("builtins").attr("int")(raw);
        return as_int.cast<std::uintptr_t>();
    }
    catch (const py::error_already_set&) {
        PyErr_Clear();
    }
    return parse_hex_pointer_from_repr(py::repr(raw).cast<std::string>());
}

Context& unwrap_context(const py::object& context_object) {
    auto* context = reinterpret_cast<Context*>(pointer_from_swig_object(context_object));
    if (context == nullptr)
        throw std::runtime_error("received a null OpenMM Context pointer");
    return *context;
}

CudaContext& cuda_context_from_openmm_context(const py::object& context_object) {
    Context& context = unwrap_context(context_object);
    ContextImpl& impl = context.getImpl();
    void* platform_data = impl.getPlatformData();
    if (platform_data == nullptr)
        throw std::runtime_error("OpenMM Context has no CUDA platform data");
    auto* data = static_cast<CudaPlatform::PlatformData*>(platform_data);
    if (data->contexts.empty() || data->contexts[0] == nullptr)
        throw std::runtime_error("OpenMM CUDA platform data has no CudaContext");
    return *data->contexts[0];
}

std::uintptr_t device_pointer(CudaArray& array) {
    return static_cast<std::uintptr_t>(array.getDevicePointer());
}

std::uintptr_t correction_pointer(CudaArray& correction, int expected_size) {
    if (!correction.isInitialized())
        return 0;
    if (correction.getSize() == 0)
        return 0;
    if (correction.getElementSize() != expected_size)
        return 0;
    return device_pointer(correction);
}

torch::Tensor allocate_xyz_tensor(CudaContext& cuda, int atoms) {
    const int device_index = cuda.getDeviceIndex();
    torch::Device device(torch::kCUDA, device_index);
    return torch::empty({atoms, 3}, torch::TensorOptions().dtype(torch::kFloat64).device(device));
}

void ensure_cuda_tensor(const torch::Tensor& tensor, int atoms, int device_index) {
    if (!tensor.defined())
        throw std::runtime_error("tensor is undefined");
    if (!tensor.is_cuda())
        throw std::runtime_error("tensor must be a CUDA tensor");
    if (tensor.scalar_type() != torch::kFloat64)
        throw std::runtime_error("tensor must have dtype float64");
    if (tensor.dim() != 2 || tensor.size(0) != atoms || tensor.size(1) != 3)
        throw std::runtime_error("tensor must have shape (num_atoms, 3)");
    if (!tensor.is_contiguous())
        throw std::runtime_error("tensor must be contiguous");
    if (tensor.get_device() != device_index)
        throw std::runtime_error("tensor device does not match the OpenMM CUDA device");
}

// posq.xyz is nm; Boltz/genai stack uses Å for coordinates. Velm.xyz is nm/ps (no Å conversion).
constexpr double kBridgeAngstromToInternalPos = 0.1;
constexpr double kBridgeInternalPosToAngstrom = 10.0;

torch::Tensor array_to_tensor(CudaContext& cuda, CudaArray& array, CudaArray* correction, double extract_scale) {
    const int atoms = static_cast<int>(array.getSize());
    torch::Tensor tensor = allocate_xyz_tensor(cuda, atoms);
    CUcontext previous = nullptr;
    cuCtxGetCurrent(&previous);
    cuCtxSetCurrent(cuda.getContext());
    auto stream = reinterpret_cast<cudaStream_t>(cuda.getCurrentStream());
    openmm_bridge_extract_xyz(
            device_pointer(array),
            correction == nullptr ? 0 : correction_pointer(*correction, array.getElementSize()),
            array.getElementSize(),
            tensor.data_ptr<double>(),
            atoms,
            extract_scale,
            stream);
    cudaStreamSynchronize(stream);
    if (previous != nullptr)
        cuCtxSetCurrent(previous);
    return tensor;
}

void tensor_to_array(
        CudaContext& cuda,
        CudaArray& array,
        CudaArray* correction,
        torch::Tensor tensor,
        double scatter_scale) {
    const int atoms = static_cast<int>(array.getSize());
    if (!tensor.is_contiguous())
        tensor = tensor.contiguous();
    ensure_cuda_tensor(tensor, atoms, cuda.getDeviceIndex());
    CUcontext previous = nullptr;
    cuCtxGetCurrent(&previous);
    cuCtxSetCurrent(cuda.getContext());
    auto stream = reinterpret_cast<cudaStream_t>(cuda.getCurrentStream());
    openmm_bridge_scatter_xyz(
            device_pointer(array),
            correction == nullptr ? 0 : correction_pointer(*correction, array.getElementSize()),
            array.getElementSize(),
            tensor.data_ptr<double>(),
            atoms,
            scatter_scale,
            stream);
    cudaStreamSynchronize(stream);
    if (previous != nullptr)
        cuCtxSetCurrent(previous);
}

} // namespace

torch::Tensor positions_to_tensor(const py::object& context_object) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    return array_to_tensor(cuda, cuda.getPosq(), &cuda.getPosqCorrection(), kBridgeInternalPosToAngstrom);
}

void tensor_to_positions(const py::object& context_object, torch::Tensor positions) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    tensor_to_array(cuda, cuda.getPosq(), &cuda.getPosqCorrection(), positions, kBridgeAngstromToInternalPos);
}

torch::Tensor velocities_to_tensor(const py::object& context_object) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    return array_to_tensor(cuda, cuda.getVelm(), nullptr, 1.0);
}

void tensor_to_velocities(const py::object& context_object, torch::Tensor velocities) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    tensor_to_array(cuda, cuda.getVelm(), nullptr, velocities, 1.0);
}

namespace {

constexpr double kInvFixedPointScale = 1.0 / static_cast<double>(0x100000000LL);
// OpenMM forces: kJ/mol/nm  → eV/Å (matches AmberGBCudaTeacher numpy conversion)
constexpr double kKjPerNmToEvPerAng = 1.0 / (96.4853321233100184 * 10.0);

torch::Tensor forces_fixed_to_torch(CudaContext& cuda, int padded_atoms) {
    CudaArray& force = cuda.getForce();
    const int n_elem = static_cast<int>(force.getSize());
    TORCH_CHECK(n_elem == padded_atoms * 3, "CudaContext.force size mismatch");
    const int device_index = cuda.getDeviceIndex();
    torch::Device device(torch::kCUDA, device_index);
    torch::Tensor flat =
            torch::empty({n_elem}, torch::TensorOptions().dtype(torch::kInt64).device(device));
    CUcontext previous = nullptr;
    cuCtxGetCurrent(&previous);
    cuCtxSetCurrent(cuda.getContext());
    auto stream = reinterpret_cast<cudaStream_t>(cuda.getCurrentStream());
    cudaError_t st = cudaMemcpyAsync(
            flat.data_ptr<int64_t>(),
            reinterpret_cast<const void*>(device_pointer(force)),
            static_cast<size_t>(n_elem) * sizeof(int64_t),
            cudaMemcpyDeviceToDevice,
            stream);
    if (st != cudaSuccess)
        throw std::runtime_error(std::string("cudaMemcpyAsync(force): ") + cudaGetErrorString(st));
    cudaStreamSynchronize(stream);
    if (previous != nullptr)
        cuCtxSetCurrent(previous);
    // OpenMM CUDA layout: x components for all atoms, then y, then z (see nonbonded.cu).
    torch::Tensor flat_fp64 = flat.to(torch::kFloat64);
    torch::Tensor fx = flat_fp64.narrow(0, 0, padded_atoms);
    torch::Tensor fy = flat_fp64.narrow(0, padded_atoms, padded_atoms);
    torch::Tensor fz = flat_fp64.narrow(0, 2 * padded_atoms, padded_atoms);
    torch::Tensor out = torch::stack({fx, fy, fz}, /*dim=*/1);
    out.mul_(kInvFixedPointScale * kKjPerNmToEvPerAng);
    return out;
}

struct CudaBridge {
    py::object context_object;
    int groups;

    CudaBridge(py::object ctx, int groups_in)
            : context_object(std::move(ctx)), groups(groups_in) {}

    torch::Tensor evaluate(torch::Tensor positions) {
        if (!positions.is_cuda())
            throw std::runtime_error("positions must be a CUDA tensor");
        torch::Tensor pos_fp64 =
                positions.scalar_type() == torch::kFloat64 ? positions : positions.to(torch::kFloat64);
        if (!pos_fp64.is_contiguous())
            pos_fp64 = pos_fp64.contiguous();
        Context& ctx = unwrap_context(context_object);
        ContextImpl& impl = ctx.getImpl();
        CudaContext& cuda = cuda_context_from_openmm_context(context_object);
        const int padded_atoms = static_cast<int>(cuda.getPosq().getSize());
        const int n_in = static_cast<int>(pos_fp64.size(0));
        TORCH_CHECK(pos_fp64.size(1) == 3, "positions must have shape (N, 3)");
        TORCH_CHECK(n_in <= padded_atoms, "positions has more atoms than the OpenMM CUDA context");

        torch::Tensor pos_upload = pos_fp64;
        if (n_in < padded_atoms) {
            torch::Tensor tail = torch::zeros(
                    {padded_atoms - n_in, 3},
                    torch::TensorOptions().dtype(torch::kFloat64).device(pos_fp64.device()));
            pos_upload = torch::cat({pos_fp64, tail}, /*dim=*/0);
        }

        // Positions tensor is Angstrom per genai/OpenMM-host convention; CUDA posq holds nm (scale in tensor_to_positions).
        tensor_to_positions(context_object, pos_upload);
        impl.calcForcesAndEnergy(true, false, groups);
        torch::Tensor forces_fp64 = forces_fixed_to_torch(cuda, padded_atoms).slice(/*dim=*/0, /*start=*/0, /*end=*/n_in);
        if (positions.scalar_type() != torch::kFloat64)
            return forces_fp64.to(positions.scalar_type());
        return forces_fp64;
    }
};

} // namespace

PYBIND11_MODULE(openmm_cuda_bridge, module) {
    module.doc() =
            "GPU-resident OpenMM CUDA Context <-> PyTorch tensor bridge. "
            "Positions: Angstrom in/out (matches host setPositions); OpenMM posq uses nm internally. "
            "Forces: eV/Å. Velocities via velocities_to_tensor/tensor_to_velocities: nm/ps (PyTorch doubles; scatter/extract scale 1.0).";
    module.def("positions_to_tensor", &positions_to_tensor, py::arg("context"));
    module.def("tensor_to_positions", &tensor_to_positions, py::arg("context"), py::arg("positions"));
    module.def("velocities_to_tensor", &velocities_to_tensor, py::arg("context"));
    module.def("tensor_to_velocities", &tensor_to_velocities, py::arg("context"), py::arg("velocities"));

    py::class_<CudaBridge>(module, "CudaBridge")
            .def(
                    py::init([](py::object ctx, int groups_in) {
                        return std::make_unique<CudaBridge>(std::move(ctx), groups_in);
                    }),
                    py::arg("context"),
                    py::arg("groups") = static_cast<int>(0xFFFFFFFF))
            .def("evaluate", &CudaBridge::evaluate, py::arg("positions"),
                 "CUDA Context forces as a tensor after GPU-side evaluation.");
}
