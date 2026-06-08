// Stable LibTorch ABI tensor ops for openmm_cuda_bridge (PyTorch >= 2.10).
#include "bridge_stable_ops.hpp"

#include <cuda_runtime_api.h>

#include <torch/csrc/inductor/aoti_torch/c/shim.h>
#include <torch/csrc/inductor/aoti_torch/generated/c_shim_aten.h>
#include <torch/csrc/stable/library.h>
#include <torch/csrc/stable/ops.h>
#include <torch/csrc/stable/stableivalue_conversions.h>
#include <torch/csrc/stable/version.h>
#include <torch/headeronly/core/ScalarType.h>
#include <torch/headeronly/util/HeaderOnlyArrayRef.h>

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace openmm_bridge_stable {
namespace {

constexpr double kInvFixedPointScale = 1.0 / static_cast<double>(0x100000000LL);
constexpr double kKjPerNmToEvPerAng = 1.0 / (96.4853321233100184 * 10.0);

torch::stable::Device cuda_device(int device_index) {
    return torch::stable::Device(torch::stable::DeviceType::CUDA, device_index);
}

torch::stable::Tensor stack_tensors(
        const std::vector<torch::stable::Tensor>& tensors,
        int64_t dim) {
    std::vector<AtenTensorHandle> handles;
    handles.reserve(tensors.size());
    for (const auto& t : tensors) {
        handles.push_back(t.get());
    }
    AtenTensorHandle ret = nullptr;
    TORCH_ERROR_CODE_CHECK(
        aoti_torch_aten_stack(handles.data(), static_cast<int64_t>(handles.size()), dim, &ret));
    return torch::stable::Tensor(ret);
}

torch::stable::Tensor mul_scalar(torch::stable::Tensor self, double value) {
    AtenTensorHandle ret = nullptr;
    TORCH_ERROR_CODE_CHECK(aoti_torch_aten_mul_Scalar(self.get(), value, &ret));
    return torch::stable::Tensor(ret);
}

} // namespace

torch::stable::Tensor allocate_xyz_tensor(int device_index, int atoms) {
    std::array<int64_t, 2> size{atoms, 3};
    return torch::stable::empty(
            torch::headeronly::IntHeaderOnlyArrayRef(size.data(), size.size()),
            torch::headeronly::ScalarType::Double,
            std::nullopt,
            cuda_device(device_index));
}

torch::stable::Tensor forces_fixed_to_torch(
        int device_index,
        std::uintptr_t force_ptr,
        int padded_atoms,
        cudaStream_t stream) {
    const int n_elem = padded_atoms * 3;
    std::array<int64_t, 1> size{static_cast<int64_t>(n_elem)};

    torch::stable::Tensor flat = torch::stable::empty(
            torch::headeronly::IntHeaderOnlyArrayRef(size.data(), size.size()),
            torch::headeronly::ScalarType::Long,
            std::nullopt,
            cuda_device(device_index));

    cudaError_t st = cudaMemcpyAsync(
            flat.mutable_data_ptr<int64_t>(),
            reinterpret_cast<const void*>(force_ptr),
            static_cast<size_t>(n_elem) * sizeof(int64_t),
            cudaMemcpyDeviceToDevice,
            stream);
    if (st != cudaSuccess) {
        throw std::runtime_error(std::string("cudaMemcpyAsync(force): ") + cudaGetErrorString(st));
    }
    cudaStreamSynchronize(stream);

    torch::stable::Tensor flat_fp64 = torch::stable::to(
            flat,
            torch::headeronly::ScalarType::Double,
            std::nullopt,
            std::nullopt,
            std::nullopt,
            false,
            false,
            std::nullopt);

    torch::stable::Tensor fx = torch::stable::narrow(flat_fp64, 0, 0, padded_atoms);
    torch::stable::Tensor fy = torch::stable::narrow(flat_fp64, 0, padded_atoms, padded_atoms);
    torch::stable::Tensor fz = torch::stable::narrow(flat_fp64, 0, 2 * padded_atoms, padded_atoms);
    torch::stable::Tensor out = stack_tensors({fx, fy, fz}, 1);
    return mul_scalar(out, kInvFixedPointScale * kKjPerNmToEvPerAng);
}

torch::stable::Tensor atom_energies_fixed_to_torch(
        int device_index,
        std::uintptr_t atom_energy_ptr,
        int padded_atoms,
        cudaStream_t stream) {
    std::array<int64_t, 1> size{static_cast<int64_t>(padded_atoms)};

    torch::stable::Tensor flat = torch::stable::empty(
            torch::headeronly::IntHeaderOnlyArrayRef(size.data(), size.size()),
            torch::headeronly::ScalarType::Long,
            std::nullopt,
            cuda_device(device_index));

    cudaError_t st = cudaMemcpyAsync(
            flat.mutable_data_ptr<int64_t>(),
            reinterpret_cast<const void*>(atom_energy_ptr),
            static_cast<size_t>(padded_atoms) * sizeof(int64_t),
            cudaMemcpyDeviceToDevice,
            stream);
    if (st != cudaSuccess) {
        throw std::runtime_error(std::string("cudaMemcpyAsync(atomEnergy): ") + cudaGetErrorString(st));
    }
    cudaStreamSynchronize(stream);

    torch::stable::Tensor flat_fp64 = torch::stable::to(
            flat,
            torch::headeronly::ScalarType::Double,
            std::nullopt,
            std::nullopt,
            std::nullopt,
            false,
            false,
            std::nullopt);
    // kJ/mol → eV, same fixed-point scale as forces.
    return mul_scalar(flat_fp64, kInvFixedPointScale * (1.0 / 96.4853321233100184));
}

torch::stable::Tensor ensure_fp64_contiguous(const torch::stable::Tensor& tensor) {
    torch::stable::Tensor out = torch::stable::to(
            tensor,
            torch::headeronly::ScalarType::Double,
            std::nullopt,
            std::nullopt,
            std::nullopt,
            false,
            false,
            std::nullopt);
    return torch::stable::contiguous(out);
}

torch::stable::Tensor pad_positions(
        const torch::stable::Tensor& positions_fp64,
        int padded_atoms) {
    const int64_t n_in = positions_fp64.size(0);
    if (n_in >= padded_atoms) {
        return positions_fp64;
    }
    std::array<int64_t, 2> tail_size{padded_atoms - n_in, 3};
    torch::stable::Tensor tail = torch::stable::zeros(
            torch::headeronly::IntHeaderOnlyArrayRef(tail_size.data(), tail_size.size()),
            torch::headeronly::ScalarType::Double,
            std::nullopt,
            positions_fp64.device());
    return torch::stable::cat({positions_fp64, tail}, /*dim=*/0);
}

torch::stable::Tensor slice_rows(const torch::stable::Tensor& tensor, int rows) {
    return torch::stable::narrow(tensor, 0, 0, rows);
}

torch::stable::Tensor cast_dtype(
        const torch::stable::Tensor& tensor,
        torch::headeronly::ScalarType dtype) {
    return torch::stable::to(
            tensor,
            dtype,
            std::nullopt,
            std::nullopt,
            std::nullopt,
            false,
            false,
            std::nullopt);
}

// Registered stable custom ops (scaffold for future expansion).
STABLE_TORCH_LIBRARY(openmm_bridge, m) {
    m.def("allocate_xyz(int device_index, int atoms) -> Tensor");
    m.def("forces_fixed_to_torch(Tensor force_buf, int padded_atoms) -> Tensor");
}

} // namespace openmm_bridge_stable
