#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ContextImpl.h"
#include "CudaContext.h"
#include "CudaPlatform.h"

#include <cuda.h>
#include <pybind11/pybind11.h>
#include <torch/extension.h>

#include <cstdint>
#include <map>
#include <sstream>
#include <string>

using namespace OpenMM;

namespace py = pybind11;

#define CHECK_CUDA_RESULT(result, cu, prefix)                                    \
    if ((result) != CUDA_SUCCESS) {                                             \
        std::stringstream msg;                                                  \
        msg << (prefix) << ": " << (cu).getErrorString(result) << " ("        \
            << (result) << ") at " << __FILE__ << ":" << __LINE__;            \
        throw OpenMMException(msg.str());                                       \
    }

static const char* BRIDGE_KERNEL_SOURCE = R"(
extern "C" __global__
void writePositions(const real* __restrict__ posTensor,
                    real4* __restrict__ posq,
                    int* __restrict__ atomIndex,
                    int numAtoms) {
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < numAtoms; atom += blockDim.x*gridDim.x) {
        int index = atomIndex[atom];
        real4 p = posq[atom];
        p.x = posTensor[3*index] * real(0.1);
        p.y = posTensor[3*index+1] * real(0.1);
        p.z = posTensor[3*index+2] * real(0.1);
        posq[atom] = p;
    }
}

extern "C" __global__
void readForces(const long long* __restrict__ forceBuffers,
                real* __restrict__ forceTensor,
                int* __restrict__ atomIndex,
                int numAtoms,
                int paddedNumAtoms) {
    const real scale = real(1.0 / (4294967296.0 * 96.4853321233100184 * 10.0));
    for (int atom = blockIdx.x*blockDim.x+threadIdx.x; atom < numAtoms; atom += blockDim.x*gridDim.x) {
        int index = atomIndex[atom];
        forceTensor[3*index]   = real(forceBuffers[atom]) * scale;
        forceTensor[3*index+1] = real(forceBuffers[atom+paddedNumAtoms]) * scale;
        forceTensor[3*index+2] = real(forceBuffers[atom+2*paddedNumAtoms]) * scale;
    }
}
)";

static CudaContext& get_cuda_context(ContextImpl& impl) {
    auto* data = static_cast<CudaPlatform::PlatformData*>(impl.getPlatformData());
    if (data == nullptr || data->contexts.empty()) {
        throw OpenMMException("Context is not backed by an initialized CUDA PlatformData");
    }
    return *data->contexts[0];
}

static void* tensor_data_pointer(CudaContext& cu, torch::Tensor& tensor) {
    if (cu.getUseDoublePrecision()) {
        return tensor.data_ptr<double>();
    }
    return tensor.data_ptr<float>();
}

class CudaBridge {
public:
    explicit CudaBridge(std::uintptr_t context_ptr, int groups = 0xFFFFFFFF)
        : context_(reinterpret_cast<Context*>(context_ptr)),
          impl_(&context_->getImplementation()),
          cu_(get_cuda_context(*impl_)),
          groups_(groups) {
        CHECK_CUDA_RESULT(cuDevicePrimaryCtxRetain(&primary_context_, cu_.getDevice()), cu_, "Failed to retain CUDA primary context");
        ContextSelector selector(cu_);
        module_ = cu_.createModule(BRIDGE_KERNEL_SOURCE, std::map<std::string, std::string>());
        write_positions_kernel_ = cu_.getKernel(module_, "writePositions");
        read_forces_kernel_ = cu_.getKernel(module_, "readForces");
    }

    ~CudaBridge() {
        cuDevicePrimaryCtxRelease(cu_.getDevice());
    }

    int num_atoms() const {
        return cu_.getNumAtoms();
    }

    int padded_num_atoms() const {
        return cu_.getPaddedNumAtoms();
    }

    int device_index() const {
        return cu_.getDeviceIndex();
    }

    void set_positions(torch::Tensor positions_angstrom) {
        validate_positions(positions_angstrom);
        torch::Tensor pos = canonical_tensor(positions_angstrom);
        void* pos_data = tensor_data_pointer(cu_, pos);
        int num_atoms = cu_.getNumAtoms();

        CHECK_CUDA_RESULT(cuCtxSynchronize(), cu_, "Failed to synchronize before switching to OpenMM context");
        {
            ContextSelector selector(cu_);
            void* args[] = {
                &pos_data,
                &cu_.getPosq().getDevicePointer(),
                &cu_.getAtomIndexArray().getDevicePointer(),
                &num_atoms,
            };
            cu_.executeKernel(write_positions_kernel_, args, num_atoms);
            CHECK_CUDA_RESULT(cuCtxSynchronize(), cu_, "Failed to synchronize after writing positions");
        }
    }

    torch::Tensor compute_forces() {
        const auto dtype = cu_.getUseDoublePrecision() ? torch::kFloat64 : torch::kFloat32;
        torch::Tensor forces = torch::empty(
            {cu_.getNumAtoms(), 3},
            torch::TensorOptions().device(torch::kCUDA, cu_.getDeviceIndex()).dtype(dtype)
        );
        void* force_data = tensor_data_pointer(cu_, forces);
        int num_atoms = cu_.getNumAtoms();
        int padded_num_atoms = cu_.getPaddedNumAtoms();

        {
            ContextSelector selector(cu_);
            impl_->computeVirtualSites();
            impl_->calcForcesAndEnergy(true, false, groups_);
            void* args[] = {
                &cu_.getForce().getDevicePointer(),
                &force_data,
                &cu_.getAtomIndexArray().getDevicePointer(),
                &num_atoms,
                &padded_num_atoms,
            };
            cu_.executeKernel(read_forces_kernel_, args, num_atoms);
            CHECK_CUDA_RESULT(cuCtxSynchronize(), cu_, "Failed to synchronize after reading forces");
        }
        CHECK_CUDA_RESULT(cuCtxPushCurrent(primary_context_), cu_, "Failed to restore CUDA primary context");
        CUcontext popped;
        CHECK_CUDA_RESULT(cuCtxPopCurrent(&popped), cu_, "Failed to pop restored CUDA primary context");
        return forces;
    }

    torch::Tensor evaluate(torch::Tensor positions_angstrom) {
        set_positions(positions_angstrom);
        return compute_forces();
    }

private:
    void validate_positions(const torch::Tensor& positions_angstrom) const {
        if (!positions_angstrom.is_cuda()) {
            throw OpenMMException("positions_angstrom must be a CUDA tensor");
        }
        if (positions_angstrom.dim() != 2 || positions_angstrom.size(1) != 3) {
            throw OpenMMException("positions_angstrom must have shape (N, 3)");
        }
        if (positions_angstrom.size(0) != cu_.getNumAtoms()) {
            std::stringstream msg;
            msg << "positions_angstrom has " << positions_angstrom.size(0)
                << " atoms, but OpenMM context has " << cu_.getNumAtoms();
            throw OpenMMException(msg.str());
        }
        if (positions_angstrom.get_device() != cu_.getDeviceIndex()) {
            std::stringstream msg;
            msg << "positions tensor is on CUDA device " << positions_angstrom.get_device()
                << ", but OpenMM context is on CUDA device " << cu_.getDeviceIndex();
            throw OpenMMException(msg.str());
        }
    }

    torch::Tensor canonical_tensor(torch::Tensor tensor) const {
        const auto dtype = cu_.getUseDoublePrecision() ? torch::kFloat64 : torch::kFloat32;
        return tensor.to(torch::TensorOptions().device(torch::kCUDA, cu_.getDeviceIndex()).dtype(dtype)).contiguous();
    }

    Context* context_;
    ContextImpl* impl_;
    CudaContext& cu_;
    int groups_;
    CUcontext primary_context_;
    CUmodule module_;
    CUfunction write_positions_kernel_;
    CUfunction read_forces_kernel_;
};

PYBIND11_MODULE(_openmm_cuda_bridge, m) {
    py::class_<CudaBridge>(m, "CudaBridge")
        .def(py::init<std::uintptr_t, int>(), py::arg("context_ptr"), py::arg("groups") = 0xFFFFFFFF)
        .def("num_atoms", &CudaBridge::num_atoms)
        .def("padded_num_atoms", &CudaBridge::padded_num_atoms)
        .def("device_index", &CudaBridge::device_index)
        .def("set_positions", &CudaBridge::set_positions)
        .def("compute_forces", &CudaBridge::compute_forces)
        .def("evaluate", &CudaBridge::evaluate);
}
