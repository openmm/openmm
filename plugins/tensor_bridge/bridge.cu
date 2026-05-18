// OpenMM CUDA posq stores xyz in nanometers. This PyTorch bridge uses Angstrom for
// positions read/written by positions_to_tensor / tensor_to_positions / CudaBridge.evaluate.
// Velocities use getVelm(): nm/ps internally; tensor_to_velocities uses scale 1.0 (no Å conversion).
#include <cuda_runtime.h>

#include <cstdint>
#include <stdexcept>
#include <string>

namespace {

constexpr int kBlockSize = 256;

inline void check_cuda(cudaError_t status, const char* label) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(label) + ": " + cudaGetErrorString(status));
}

__global__ void extract_xyz_f4_kernel(
        const float4* src,
        const float4* correction,
        double* dst,
        int n,
        double extract_scale) {
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom >= n)
        return;
    const float4 p = src[atom];
    double x = static_cast<double>(p.x);
    double y = static_cast<double>(p.y);
    double z = static_cast<double>(p.z);
    if (correction != nullptr) {
        const float4 c = correction[atom];
        x += static_cast<double>(c.x);
        y += static_cast<double>(c.y);
        z += static_cast<double>(c.z);
    }
    dst[3 * atom] = x * extract_scale;
    dst[3 * atom + 1] = y * extract_scale;
    dst[3 * atom + 2] = z * extract_scale;
}

__global__ void scatter_xyz_f4_kernel(
        float4* dst,
        float4* correction,
        const double* src,
        int n,
        double scatter_scale) {
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom >= n)
        return;
    float4 p = dst[atom];
    const double sx = src[3 * atom] * scatter_scale;
    const double sy = src[3 * atom + 1] * scatter_scale;
    const double sz = src[3 * atom + 2] * scatter_scale;
    const float x = static_cast<float>(sx);
    const float y = static_cast<float>(sy);
    const float z = static_cast<float>(sz);
    p.x = x;
    p.y = y;
    p.z = z;
    dst[atom] = p;
    if (correction != nullptr) {
        float4 c = correction[atom];
        c.x = static_cast<float>(sx - static_cast<double>(x));
        c.y = static_cast<float>(sy - static_cast<double>(y));
        c.z = static_cast<float>(sz - static_cast<double>(z));
        correction[atom] = c;
    }
}

__global__ void extract_xyz_d4_kernel(const double4* src, double* dst, int n, double extract_scale) {
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom >= n)
        return;
    const double4 p = src[atom];
    dst[3 * atom] = p.x * extract_scale;
    dst[3 * atom + 1] = p.y * extract_scale;
    dst[3 * atom + 2] = p.z * extract_scale;
}

__global__ void scatter_xyz_d4_kernel(double4* dst, const double* src, int n, double scatter_scale) {
    const int atom = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom >= n)
        return;
    double4 p = dst[atom];
    p.x = src[3 * atom] * scatter_scale;
    p.y = src[3 * atom + 1] * scatter_scale;
    p.z = src[3 * atom + 2] * scatter_scale;
    dst[atom] = p;
}

} // namespace

extern "C" void openmm_bridge_extract_xyz(
        std::uintptr_t src_ptr,
        std::uintptr_t correction_ptr,
        int element_size,
        double* dst,
        int n,
        double extract_scale,
        cudaStream_t stream) {
    const int blocks = (n + kBlockSize - 1) / kBlockSize;
    if (element_size == static_cast<int>(sizeof(float4))) {
        extract_xyz_f4_kernel<<<blocks, kBlockSize, 0, stream>>>(
                reinterpret_cast<const float4*>(src_ptr),
                correction_ptr == 0 ? nullptr : reinterpret_cast<const float4*>(correction_ptr),
                dst,
                n,
                extract_scale);
    }
    else if (element_size == static_cast<int>(sizeof(double4))) {
        extract_xyz_d4_kernel<<<blocks, kBlockSize, 0, stream>>>(
                reinterpret_cast<const double4*>(src_ptr),
                dst,
                n,
                extract_scale);
    }
    else {
        throw std::runtime_error("unsupported OpenMM CudaArray element size for xyz extraction");
    }
    check_cuda(cudaGetLastError(), "openmm_bridge_extract_xyz launch");
}

extern "C" void openmm_bridge_scatter_xyz(
        std::uintptr_t dst_ptr,
        std::uintptr_t correction_ptr,
        int element_size,
        const double* src,
        int n,
        double scatter_scale,
        cudaStream_t stream) {
    const int blocks = (n + kBlockSize - 1) / kBlockSize;
    if (element_size == static_cast<int>(sizeof(float4))) {
        scatter_xyz_f4_kernel<<<blocks, kBlockSize, 0, stream>>>(
                reinterpret_cast<float4*>(dst_ptr),
                correction_ptr == 0 ? nullptr : reinterpret_cast<float4*>(correction_ptr),
                src,
                n,
                scatter_scale);
    }
    else if (element_size == static_cast<int>(sizeof(double4))) {
        scatter_xyz_d4_kernel<<<blocks, kBlockSize, 0, stream>>>(
                reinterpret_cast<double4*>(dst_ptr),
                src,
                n,
                scatter_scale);
    }
    else {
        throw std::runtime_error("unsupported OpenMM CudaArray element size for xyz scatter");
    }
    check_cuda(cudaGetLastError(), "openmm_bridge_scatter_xyz launch");
}
