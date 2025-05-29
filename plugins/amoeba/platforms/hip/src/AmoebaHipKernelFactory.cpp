/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
 * Portions copyright (c) 2021 Advanced Micro Devices, Inc.                   *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "AmoebaHipKernelFactory.h"
#include "AmoebaHipKernels.h"
#include "AmoebaCommonKernels.h"
#include "HipPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"

using namespace OpenMM;

#ifdef OPENMM_BUILDING_STATIC_LIBRARY
static void registerPlatforms() {
#else
extern "C" OPENMM_EXPORT void registerPlatforms() {
#endif
}

#ifdef OPENMM_BUILDING_STATIC_LIBRARY
static void registerKernelFactories() {
#else
extern "C" OPENMM_EXPORT void registerKernelFactories() {
#endif
    try {
        Platform& platform = Platform::getPlatformByName("HIP");
        AmoebaHipKernelFactory* factory = new AmoebaHipKernelFactory();
        platform.registerKernelFactory(CalcAmoebaTorsionTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaMultipoleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaGeneralizedKirkwoodForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaVdwForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaWcaDispersionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcHippoNonbondedForceKernel::Name(), factory);
    }
    catch (...) {
        // Ignore.  The HIP platform isn't available.
    }
}

extern "C" OPENMM_EXPORT void registerAmoebaHipKernelFactories() {
    try {
        Platform::getPlatformByName("HIP");
    }
    catch (...) {
        Platform::registerPlatform(new HipPlatform());
    }
    registerKernelFactories();
}

KernelImpl* AmoebaHipKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    HipPlatform::PlatformData& data = *static_cast<HipPlatform::PlatformData*>(context.getPlatformData());
    HipContext& cu = *data.contexts[0];

    if (name == CalcAmoebaTorsionTorsionForceKernel::Name())
        return new CommonCalcAmoebaTorsionTorsionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaMultipoleForceKernel::Name())
        return new HipCalcAmoebaMultipoleForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaGeneralizedKirkwoodForceKernel::Name())
        return new CommonCalcAmoebaGeneralizedKirkwoodForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaVdwForceKernel::Name())
        return new CommonCalcAmoebaVdwForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaWcaDispersionForceKernel::Name())
        return new CommonCalcAmoebaWcaDispersionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcHippoNonbondedForceKernel::Name())
        return new HipCalcHippoNonbondedForceKernel(name, platform, cu, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
