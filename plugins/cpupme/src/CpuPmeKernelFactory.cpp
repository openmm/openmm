/* -------------------------------------------------------------------------- *
 *                              OpenMMCpuPme                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

#include "CpuPmeKernelFactory.h"
#include "CpuPmeKernels.h"
#include "internal/windowsExportPme.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" OPENMM_EXPORT_PME void registerKernelFactories() {
    if (CpuCalcPmeReciprocalForceKernel::isProcessorSupported()) {
        CpuPmeKernelFactory* factory = new CpuPmeKernelFactory();
        for (int i = 0; i < Platform::getNumPlatforms(); i++) {
            Platform::getPlatform(i).registerKernelFactory(CalcPmeReciprocalForceKernel::Name(), factory);
            Platform::getPlatform(i).registerKernelFactory(CalcDispersionPmeReciprocalForceKernel::Name(), factory);
        }
    }
}

#ifdef OPENMM_PME_BUILDING_STATIC_LIBRARY
extern "C" void registerCpuPmeKernelFactories() {
    registerKernelFactories();
}
#else
extern "C" OPENMM_EXPORT_PME void registerCpuPmeKernelFactories() {
    registerKernelFactories();
}
extern "C" OPENMM_EXPORT_PME void registerPlatforms() {
}
#endif

KernelImpl* CpuPmeKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    if (name == CalcPmeReciprocalForceKernel::Name())
        return new CpuCalcPmeReciprocalForceKernel(name, platform);
    if (name == CalcDispersionPmeReciprocalForceKernel::Name())
        return new CpuCalcDispersionPmeReciprocalForceKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
