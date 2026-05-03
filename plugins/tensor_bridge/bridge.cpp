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
        cudaStream_t stream);

extern "C" void openmm_bridge_scatter_xyz(
        std::uintptr_t dst_ptr,
        std::uintptr_t correction_ptr,
        int element_size,
        const double* src,
        int n,
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

torch::Tensor array_to_tensor(CudaContext& cuda, CudaArray& array, CudaArray* correction) {
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
            stream);
    cudaStreamSynchronize(stream);
    if (previous != nullptr)
        cuCtxSetCurrent(previous);
    return tensor;
}

void tensor_to_array(CudaContext& cuda, CudaArray& array, CudaArray* correction, torch::Tensor tensor) {
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
            stream);
    cudaStreamSynchronize(stream);
    if (previous != nullptr)
        cuCtxSetCurrent(previous);
}

} // namespace

torch::Tensor positions_to_tensor(const py::object& context_object) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    return array_to_tensor(cuda, cuda.getPosq(), &cuda.getPosqCorrection());
}

void tensor_to_positions(const py::object& context_object, torch::Tensor positions) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    tensor_to_array(cuda, cuda.getPosq(), &cuda.getPosqCorrection(), positions);
}

torch::Tensor velocities_to_tensor(const py::object& context_object) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    return array_to_tensor(cuda, cuda.getVelm(), nullptr);
}

void tensor_to_velocities(const py::object& context_object, torch::Tensor velocities) {
    CudaContext& cuda = cuda_context_from_openmm_context(context_object);
    tensor_to_array(cuda, cuda.getVelm(), nullptr, velocities);
}

PYBIND11_MODULE(openmm_cuda_bridge, module) {
    module.doc() = "GPU-resident OpenMM CUDA Context <-> PyTorch tensor bridge";
    module.def("positions_to_tensor", &positions_to_tensor, py::arg("context"));
    module.def("tensor_to_positions", &tensor_to_positions, py::arg("context"), py::arg("positions"));
    module.def("velocities_to_tensor", &velocities_to_tensor, py::arg("context"));
    module.def("tensor_to_velocities", &tensor_to_velocities, py::arg("context"), py::arg("velocities"));
}
