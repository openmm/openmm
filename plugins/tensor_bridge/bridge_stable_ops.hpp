#pragma once

#include <cuda_runtime_api.h>

#include <torch/csrc/stable/device.h>
#include <torch/csrc/stable/ops.h>
#include <torch/csrc/stable/tensor.h>
#include <torch/headeronly/core/ScalarType.h>

namespace openmm_bridge_stable {

torch::stable::Tensor allocate_xyz_tensor(int device_index, int atoms);

torch::stable::Tensor ensure_fp64_contiguous(const torch::stable::Tensor& tensor);

torch::stable::Tensor pad_positions(
        const torch::stable::Tensor& positions_fp64,
        int padded_atoms);

torch::stable::Tensor slice_rows(const torch::stable::Tensor& tensor, int rows);

torch::stable::Tensor cast_dtype(
        const torch::stable::Tensor& tensor,
        torch::headeronly::ScalarType dtype);

torch::stable::Tensor forces_fixed_to_torch(
        int device_index,
        std::uintptr_t force_ptr,
        int padded_atoms,
        cudaStream_t stream);

torch::stable::Tensor atom_energies_fixed_to_torch(
        int device_index,
        std::uintptr_t atom_energy_ptr,
        int padded_atoms,
        cudaStream_t stream);

} // namespace openmm_bridge_stable
