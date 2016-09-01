/* -------------------------------------------------------------------------- *
 *                           OpenMMCudaCompiler                               *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "CudaCompilerKernelFactory.h"
#include "CudaCompilerKernels.h"
#include "internal/windowsExportCudaCompiler.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

#ifdef OPENMM_CUDACOMPILER_BUILDING_STATIC_LIBRARY
static void registerKernelFactories() {
#else
extern "C" OPENMM_EXPORT_CUDACOMPILER void registerKernelFactories() {
#endif
    try {
        // Make sure this is at least CUDA 7.0.
        
        int driverVersion;
        cuDriverGetVersion(&driverVersion);
        if (driverVersion >= 7000) {
            Platform& platform = Platform::getPlatformByName("CUDA");
            CudaCompilerKernelFactory* factory = new CudaCompilerKernelFactory();
            platform.registerKernelFactory(CudaCompilerKernel::Name(), factory);
        }
    }
    catch (std::exception ex) {
        // Ignore
    }
}

#ifdef OPENMM_CUDACOMPILER_BUILDING_STATIC_LIBRARY
extern "C" void registerCudaCompilerKernelFactories() {
    registerKernelFactories();
}
#else
extern "C" OPENMM_EXPORT_CUDACOMPILER void registerCudaCompilerKernelFactories() {
    registerKernelFactories();
}
extern "C" OPENMM_EXPORT_CUDACOMPILER void registerPlatforms() {
}
#endif

KernelImpl* CudaCompilerKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    if (name == CudaCompilerKernel::Name())
        return new CudaRuntimeCompilerKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}