#pragma once

#include <torch/csrc/inductor/aoti_torch/utils.h>
#include <torch/csrc/stable/tensor.h>
#include <torch/extension.h>

namespace openmm_bridge_interop {

inline torch::stable::Tensor stable_from_aten(torch::Tensor tensor) {
    at::Tensor moved = std::move(tensor);
    return torch::stable::Tensor(torch::aot_inductor::new_tensor_handle(std::move(moved)));
}

inline torch::Tensor aten_from_stable(const torch::stable::Tensor& tensor) {
    at::Tensor* ptr = torch::aot_inductor::tensor_handle_to_tensor_pointer(tensor.get());
    if (ptr == nullptr) {
        throw std::runtime_error("null stable tensor handle");
    }
    return *ptr;
}

} // namespace openmm_bridge_interop
