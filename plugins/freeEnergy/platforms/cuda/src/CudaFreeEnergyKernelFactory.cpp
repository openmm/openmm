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
#include "kernels/GpuFreeEnergyCudaKernels.h"

using namespace OpenMM;

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

    // (void) fprintf( stderr, "initOpenMMCudaFreeEnergyPlugin called\n");
    if ( gpuIsAvailableSoftcore() ){
        for( int ii = 0; ii < Platform::getNumPlatforms(); ii++ ){
            Platform& platform = Platform::getPlatform(ii);
            if( platform.getName().compare( "Cuda" ) == 0 ){
                 CudaFreeEnergyKernelFactory* factory = new CudaFreeEnergyKernelFactory();
                 platform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
                 platform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
                 platform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
            }
        }
    }   
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
