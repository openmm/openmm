/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#include "AmoebaReferenceKernelFactory.h"
#include "AmoebaReferenceKernels.h"
#include "ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

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
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
             AmoebaReferenceKernelFactory* factory = new AmoebaReferenceKernelFactory();
             platform.registerKernelFactory(CalcAmoebaBondForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaAngleForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaInPlaneAngleForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaPiTorsionForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaStretchBendForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaOutOfPlaneBendForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaTorsionTorsionForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaVdwForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaMultipoleForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaGeneralizedKirkwoodForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcAmoebaWcaDispersionForceKernel::Name(), factory);
             platform.registerKernelFactory(CalcHippoNonbondedForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerAmoebaReferenceKernelFactories() {
    registerKernelFactories();
}

KernelImpl* AmoebaReferenceKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& referencePlatformData = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());

    // create AmoebaReferenceData object if contextToAmoebaDataMap does not contain
    // key equal to current context
    if (name == CalcAmoebaBondForceKernel::Name())
        return new ReferenceCalcAmoebaBondForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaAngleForceKernel::Name())
        return new ReferenceCalcAmoebaAngleForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaInPlaneAngleForceKernel::Name())
        return new ReferenceCalcAmoebaInPlaneAngleForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaPiTorsionForceKernel::Name())
        return new ReferenceCalcAmoebaPiTorsionForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaStretchBendForceKernel::Name())
        return new ReferenceCalcAmoebaStretchBendForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaOutOfPlaneBendForceKernel::Name())
        return new ReferenceCalcAmoebaOutOfPlaneBendForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaTorsionTorsionForceKernel::Name())
        return new ReferenceCalcAmoebaTorsionTorsionForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaVdwForceKernel::Name())
        return new ReferenceCalcAmoebaVdwForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaMultipoleForceKernel::Name())
        return new ReferenceCalcAmoebaMultipoleForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaGeneralizedKirkwoodForceKernel::Name())
        return new ReferenceCalcAmoebaGeneralizedKirkwoodForceKernel(name, platform, context.getSystem());

    if (name == CalcAmoebaWcaDispersionForceKernel::Name())
        return new ReferenceCalcAmoebaWcaDispersionForceKernel(name, platform, context.getSystem());

    if (name == CalcHippoNonbondedForceKernel::Name())
        return new ReferenceCalcHippoNonbondedForceKernel(name, platform, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}