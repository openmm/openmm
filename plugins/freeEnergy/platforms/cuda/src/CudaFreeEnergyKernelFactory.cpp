/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "CudaFreeEnergyKernelFactory.h"
#include "CudaFreeEnergyKernels.h"
#include "openmm/freeEnergyKernels.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

#if defined(OPENMM_BUILDING_SHARED_LIBRARY)
    #if defined(WIN32)
      #include <windows.h>
        extern "C" void initOpenMMCudaFreeEnergyPlugin();
        BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
            if (ul_reason_for_call == DLL_PROCESS_ATTACH)
                initOpenMMCudaFreeEnergyPlugin();
            return TRUE;
        }   
    #else
        extern "C" void __attribute__((constructor)) initOpenMMCudaFreeEnergyPlugin();
    #endif
#endif

using namespace OpenMM;

extern "C" void initOpenMMCudaFreeEnergyPlugin() {

//(void) fprintf( stderr, "initOpenMMCudaFreeEnergyPlugin called ");
    if ( gpuIsAvailable() ){
        CudaPlatform* cudaPlatform            = new CudaPlatform();
        CudaFreeEnergyKernelFactory* factory  = new CudaFreeEnergyKernelFactory();
    
//(void) fprintf( stderr, "gpu is available platform=%p", cudaPlatform);
        cudaPlatform->registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
        cudaPlatform->registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
        cudaPlatform->registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
    
        Platform::registerPlatform(cudaPlatform);
    }
//(void) fprintf( stderr, "\n");
}

KernelImpl* CudaFreeEnergyKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {

    CudaPlatform::PlatformData& data = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());

    if (name == CalcNonbondedSoftcoreForceKernel::Name())
        return new CudaFreeEnergyCalcNonbondedSoftcoreForceKernel(name, platform, data, context.getSystem());

    if (name == CalcGBSAOBCSoftcoreForceKernel::Name())
        return new CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel(name, platform, data);

    if (name == CalcGBVISoftcoreForceKernel::Name())
        return new CudaFreeEnergyCalcGBVISoftcoreForceKernel(name, platform, data);

    throw OpenMMException( (std::string("Tried to create kernel with illegal kernel name '") + name + "'").c_str() );
}
