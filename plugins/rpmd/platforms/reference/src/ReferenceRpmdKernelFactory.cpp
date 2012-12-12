/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
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

#include "ReferenceRpmdKernelFactory.h"
#include "ReferenceRpmdKernels.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

#if defined(WIN32)
    #include <windows.h>
    extern "C" void initRpmdReferenceKernels();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            initRpmdReferenceKernels();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) initRpmdReferenceKernels();
#endif

extern "C" void initRpmdReferenceKernels() {
    Platform& platform = Platform::getPlatformByName("Reference");
    ReferenceRpmdKernelFactory* factory = new ReferenceRpmdKernelFactory();
    platform.registerKernelFactory(IntegrateRPMDStepKernel::Name(), factory);
}

KernelImpl* ReferenceRpmdKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    if (name == IntegrateRPMDStepKernel::Name())
        return new ReferenceIntegrateRPMDStepKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
