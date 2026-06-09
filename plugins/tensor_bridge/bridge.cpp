#include <pybind11/pybind11.h>

#ifdef GENAI_TPS_STABLE_TORCH_ABI
#include "bridge_interop.hpp"
#include "bridge_stable_ops.hpp"
#else
#include <torch/extension.h>
#endif

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
#ifdef GENAI_TPS_STABLE_TORCH_ABI
    return openmm_bridge_interop::aten_from_stable(
            openmm_bridge_stable::allocate_xyz_tensor(device_index, atoms));
#else
    torch::Device device(torch::kCUDA, device_index);
    return torch::empty({atoms, 3}, torch::TensorOptions().dtype(torch::kFloat64).device(device));
#endif
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
    // Sync required before OpenMM reads posq on the same stream; deferring this
    // requires profile evidence (see scripts/profile/run_openmm_bridge_baseline.sh).
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
// OpenMM energies: kJ/mol → eV (matches the eV/Å force scaling above).
constexpr double kKjPerMolToEv = 1.0 / 96.4853321233100184;

torch::Tensor forces_fixed_to_torch(CudaContext& cuda, int padded_atoms) {
    CudaArray& force = cuda.getForce();
    const int n_elem = static_cast<int>(force.getSize());
#ifdef GENAI_TPS_STABLE_TORCH_ABI
    if (n_elem != padded_atoms * 3) {
        throw std::runtime_error("CudaContext.force size mismatch");
    }
    CUcontext previous = nullptr;
    cuCtxGetCurrent(&previous);
    cuCtxSetCurrent(cuda.getContext());
    auto stream = reinterpret_cast<cudaStream_t>(cuda.getCurrentStream());
    torch::stable::Tensor out = openmm_bridge_stable::forces_fixed_to_torch(
            cuda.getDeviceIndex(),
            device_pointer(force),
            padded_atoms,
            stream);
    if (previous != nullptr) {
        cuCtxSetCurrent(previous);
    }
    return openmm_bridge_interop::aten_from_stable(out);
#else
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
    if (st != cudaSuccess) {
        throw std::runtime_error(std::string("cudaMemcpyAsync(force): ") + cudaGetErrorString(st));
    }
    cudaStreamSynchronize(stream);
    if (previous != nullptr) {
        cuCtxSetCurrent(previous);
    }
    torch::Tensor flat_fp64 = flat.to(torch::kFloat64);
    torch::Tensor fx = flat_fp64.narrow(0, 0, padded_atoms);
    torch::Tensor fy = flat_fp64.narrow(0, padded_atoms, padded_atoms);
    torch::Tensor fz = flat_fp64.narrow(0, 2 * padded_atoms, padded_atoms);
    torch::Tensor out = torch::stack({fx, fy, fz}, /*dim=*/1);
    out.mul_(kInvFixedPointScale * kKjPerNmToEvPerAng);
    return out;
#endif
}

// Read the genai-tps per-atom energy buffer (long long fixed point, paddedNumAtoms)
// as a 1-D tensor of length padded_atoms in eV. Same fixed-point scale as forces,
// energy unit conversion kJ/mol → eV.
torch::Tensor atom_energy_fixed_to_torch(CudaContext& cuda, int padded_atoms) {
    CudaArray& atomEnergy = cuda.getAtomEnergyBuffer();
    const int n_elem = static_cast<int>(atomEnergy.getSize());
#ifdef GENAI_TPS_STABLE_TORCH_ABI
    if (n_elem != padded_atoms) {
        throw std::runtime_error("CudaContext.atomEnergyBuffer size mismatch");
    }
    CUcontext previous = nullptr;
    cuCtxGetCurrent(&previous);
    cuCtxSetCurrent(cuda.getContext());
    auto stream = reinterpret_cast<cudaStream_t>(cuda.getCurrentStream());
    torch::stable::Tensor out = openmm_bridge_stable::atom_energies_fixed_to_torch(
            cuda.getDeviceIndex(),
            device_pointer(atomEnergy),
            padded_atoms,
            stream);
    if (previous != nullptr) {
        cuCtxSetCurrent(previous);
    }
    return openmm_bridge_interop::aten_from_stable(out);
#else
    TORCH_CHECK(n_elem == padded_atoms, "CudaContext.atomEnergyBuffer size mismatch");
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
            reinterpret_cast<const void*>(device_pointer(atomEnergy)),
            static_cast<size_t>(n_elem) * sizeof(int64_t),
            cudaMemcpyDeviceToDevice,
            stream);
    if (st != cudaSuccess) {
        throw std::runtime_error(std::string("cudaMemcpyAsync(atomEnergy): ") + cudaGetErrorString(st));
    }
    cudaStreamSynchronize(stream);
    if (previous != nullptr) {
        cuCtxSetCurrent(previous);
    }
    torch::Tensor out = flat.to(torch::kFloat64);
    out.mul_(kInvFixedPointScale * kKjPerMolToEv);
    return out;
#endif
}

struct CudaBridge {
    py::object context_object;
    int groups;

    CudaBridge(py::object ctx, int groups_in)
            : context_object(std::move(ctx)), groups(groups_in) {}

    // Upload positions, run calcForcesAndEnergy, and return forces (n_in, 3) in
    // eV/Å. When ``include_energy`` is true, the group's total potential energy
    // (eV) from the *same* evaluation is written to ``energy_ev_out``. When
    // ``atom_energy_out`` is non-null (requires include_energy), the genai-tps
    // per-atom energy buffer (eV, sliced to n_in) is written there.
    torch::Tensor run_forces(
            torch::Tensor positions,
            bool include_energy,
            double& energy_ev_out,
            torch::Tensor* atom_energy_out = nullptr) {
        if (!positions.is_cuda())
            throw std::runtime_error("positions must be a CUDA tensor");
        Context& ctx = unwrap_context(context_object);
        ContextImpl& impl = ctx.getImpl();
        CudaContext& cuda = cuda_context_from_openmm_context(context_object);
        const int padded_atoms = static_cast<int>(cuda.getPosq().getSize());
#ifdef GENAI_TPS_STABLE_TORCH_ABI
        torch::stable::Tensor pos_stable = openmm_bridge_interop::stable_from_aten(positions);
        torch::stable::Tensor pos_fp64 = openmm_bridge_stable::ensure_fp64_contiguous(pos_stable);
        const int n_in = static_cast<int>(pos_fp64.size(0));
        if (pos_fp64.size(1) != 3)
            throw std::runtime_error("positions must have shape (N, 3)");
        if (n_in > padded_atoms)
            throw std::runtime_error("positions has more atoms than the OpenMM CUDA context");

        torch::stable::Tensor pos_upload = openmm_bridge_stable::pad_positions(pos_fp64, padded_atoms);
        tensor_to_positions(context_object, openmm_bridge_interop::aten_from_stable(pos_upload));
        const double energy_kj = impl.calcForcesAndEnergy(true, include_energy, groups);
        if (include_energy)
            energy_ev_out = energy_kj * kKjPerMolToEv;
        torch::Tensor forces_fp64 = forces_fixed_to_torch(cuda, padded_atoms).slice(/*dim=*/0, /*start=*/0, /*end=*/n_in);
        if (atom_energy_out != nullptr) {
            torch::Tensor ae = atom_energy_fixed_to_torch(cuda, padded_atoms).slice(/*dim=*/0, /*start=*/0, /*end=*/n_in);
            *atom_energy_out = (positions.scalar_type() != torch::kFloat64) ? ae.to(positions.scalar_type()) : ae;
        }
        if (positions.scalar_type() != torch::kFloat64)
            return forces_fp64.to(positions.scalar_type());
        return forces_fp64;
#else
        torch::Tensor pos_fp64 =
                positions.scalar_type() == torch::kFloat64 ? positions : positions.to(torch::kFloat64);
        if (!pos_fp64.is_contiguous())
            pos_fp64 = pos_fp64.contiguous();
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
        const double energy_kj = impl.calcForcesAndEnergy(true, include_energy, groups);
        if (include_energy)
            energy_ev_out = energy_kj * kKjPerMolToEv;
        torch::Tensor forces_fp64 = forces_fixed_to_torch(cuda, padded_atoms).slice(/*dim=*/0, /*start=*/0, /*end=*/n_in);
        if (atom_energy_out != nullptr) {
            torch::Tensor ae = atom_energy_fixed_to_torch(cuda, padded_atoms).slice(/*dim=*/0, /*start=*/0, /*end=*/n_in);
            *atom_energy_out = (positions.scalar_type() != torch::kFloat64) ? ae.to(positions.scalar_type()) : ae;
        }
        if (positions.scalar_type() != torch::kFloat64)
            return forces_fp64.to(positions.scalar_type());
        return forces_fp64;
#endif
    }

    torch::Tensor evaluate(torch::Tensor positions) {
        double energy_ev = 0.0;  // discarded
        return run_forces(std::move(positions), /*include_energy=*/false, energy_ev);
    }

    // Force-and-energy variant: returns (forces (n_in, 3) eV/Å, scalar energy eV)
    // for the bound force group(s) from a single calcForcesAndEnergy evaluation.
    // This is the Tier-1 per-group *scalar* energy used to validate the Tier-2
    // torch bonded per-atom decomposition (and the Tier-3 per-atom kernels).
    std::pair<torch::Tensor, double> evaluate_with_energy(torch::Tensor positions) {
        double energy_ev = 0.0;
        torch::Tensor forces = run_forces(std::move(positions), /*include_energy=*/true, energy_ev);
        return std::make_pair(std::move(forces), energy_ev);
    }

    // Tier-3 variant: returns (forces (n_in, 3) eV/Å, per-atom energy (n_in,) eV)
    // for the bound force group(s). The per-atom energy comes from the genai-tps
    // atomEnergyBuffer populated by the forked nonbonded/GB kernels; summing it
    // over atoms recovers the scalar group energy from evaluate_with_energy.
    std::pair<torch::Tensor, torch::Tensor> evaluate_with_atom_energy(torch::Tensor positions) {
        double energy_ev = 0.0;
        torch::Tensor atom_energy;
        torch::Tensor forces = run_forces(
                std::move(positions), /*include_energy=*/true, energy_ev, &atom_energy);
        return std::make_pair(std::move(forces), std::move(atom_energy));
    }

    torch::Tensor evaluate_batch(torch::Tensor positions) {
        TORCH_CHECK(positions.dim() == 3, "positions must have shape (B, N, 3)");
        TORCH_CHECK(positions.size(2) == 3, "positions must have shape (B, N, 3)");
        const int64_t batch_size = positions.size(0);
        std::vector<torch::Tensor> force_rows;
        force_rows.reserve(static_cast<size_t>(batch_size));
        for (int64_t i = 0; i < batch_size; ++i) {
            force_rows.push_back(evaluate(positions.select(0, i).contiguous()));
        }
        return torch::stack(force_rows, /*dim=*/0);
    }

    // Batch Tier-3: returns (forces (B, n_in, 3), per-atom energy (B, n_in)).
    // Serial loop over the batch, matching evaluate_batch semantics.
    std::pair<torch::Tensor, torch::Tensor> evaluate_batch_with_atom_energy(torch::Tensor positions) {
        TORCH_CHECK(positions.dim() == 3, "positions must have shape (B, N, 3)");
        TORCH_CHECK(positions.size(2) == 3, "positions must have shape (B, N, 3)");
        const int64_t batch_size = positions.size(0);
        std::vector<torch::Tensor> force_rows;
        std::vector<torch::Tensor> energy_rows;
        force_rows.reserve(static_cast<size_t>(batch_size));
        energy_rows.reserve(static_cast<size_t>(batch_size));
        for (int64_t i = 0; i < batch_size; ++i) {
            auto fe = evaluate_with_atom_energy(positions.select(0, i).contiguous());
            force_rows.push_back(std::move(fe.first));
            energy_rows.push_back(std::move(fe.second));
        }
        return std::make_pair(
                torch::stack(force_rows, /*dim=*/0),
                torch::stack(energy_rows, /*dim=*/0));
    }

    // Batch super-system Tier-3: pack (B, N, 3) with per-replica X offsets, one
    // calcForcesAndEnergy on the replicated Context, return (B, N, 3) / (B, N).
    std::pair<torch::Tensor, torch::Tensor> evaluate_super_batch_with_atom_energy(
            torch::Tensor positions,
            int64_t n_per_replica,
            double offset_angstrom) {
        TORCH_CHECK(positions.dim() == 3, "positions must have shape (B, N, 3)");
        TORCH_CHECK(positions.size(2) == 3, "positions must have shape (B, N, 3)");
        const int64_t batch_size = positions.size(0);
        const int64_t n = positions.size(1);
        TORCH_CHECK(n == n_per_replica, "n_per_replica must match positions.size(1)");

        torch::Tensor pos64 =
                positions.scalar_type() == torch::kFloat64 ? positions : positions.to(torch::kFloat64);
        if (!pos64.is_contiguous())
            pos64 = pos64.contiguous();

        torch::Tensor offsets =
                torch::arange(batch_size, pos64.options().dtype(torch::kFloat64)).view({-1, 1, 1}) *
                offset_angstrom;
        torch::Tensor x = pos64.select(/*dim=*/2, /*index=*/0).unsqueeze(/*dim=*/2) + offsets;
        torch::Tensor super_pos =
                torch::cat({x, pos64.slice(/*dim=*/2, /*start=*/1, /*end=*/3)}, /*dim=*/2)
                        .reshape({batch_size * n, 3});

        auto fe = evaluate_with_atom_energy(super_pos.contiguous());
        torch::Tensor forces = fe.first.view({batch_size, n, 3});
        torch::Tensor energies = fe.second.view({batch_size, n});
        return std::make_pair(std::move(forces), std::move(energies));
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
                 "CUDA Context forces as a tensor after GPU-side evaluation.")
            .def("evaluate_with_energy", &CudaBridge::evaluate_with_energy, py::arg("positions"),
                 "(forces (N,3) eV/Å, scalar potential energy eV) for the bound force "
                 "group(s) from a single calcForcesAndEnergy evaluation.")
            .def("evaluate_batch", &CudaBridge::evaluate_batch, py::arg("positions"),
                 "Batch forces for shape (B, N, 3). Current implementation preserves "
                 "single-Context semantics by evaluating each row in sequence.")
            .def("evaluate_with_atom_energy", &CudaBridge::evaluate_with_atom_energy, py::arg("positions"),
                 "(forces (N,3) eV/Å, per-atom energy (N,) eV) for the bound force group(s); "
                 "per-atom energy sums to the scalar group energy.")
            .def("evaluate_batch_with_atom_energy", &CudaBridge::evaluate_batch_with_atom_energy,
                 py::arg("positions"),
                 "Batch (forces (B,N,3), per-atom energy (B,N)) via serial per-row evaluation.")
            .def("evaluate_super_batch_with_atom_energy",
                 &CudaBridge::evaluate_super_batch_with_atom_energy,
                 py::arg("positions"),
                 py::arg("n_per_replica"),
                 py::arg("offset_angstrom"),
                 "Replicated super-system batch: one calcForcesAndEnergy over B replicas.");
}
