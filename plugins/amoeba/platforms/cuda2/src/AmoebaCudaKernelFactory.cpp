/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "AmoebaCudaKernelFactory.h"
#include "AmoebaCudaKernels.h"
#include "CudaPlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"

using namespace OpenMM;

extern "C" void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        AmoebaCudaKernelFactory* factory = new AmoebaCudaKernelFactory();
        platform.registerKernelFactory(CalcAmoebaHarmonicBondForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaUreyBradleyForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaHarmonicAngleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaHarmonicInPlaneAngleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaPiTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaStretchBendForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaOutOfPlaneBendForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaTorsionTorsionForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaMultipoleForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaGeneralizedKirkwoodForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaVdwForceKernel::Name(), factory);
        platform.registerKernelFactory(CalcAmoebaWcaDispersionForceKernel::Name(), factory);
    }
    catch (...) {
        // Ignore.  The CUDA platform isn't available.
    }
}

extern "C" OPENMM_EXPORT void registerAmoebaCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}

KernelImpl* AmoebaCudaKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaPlatform::PlatformData& data = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    CudaContext& cu = *data.contexts[0];

    if (name == CalcAmoebaHarmonicBondForceKernel::Name())
        return new CudaCalcAmoebaHarmonicBondForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaHarmonicAngleForceKernel::Name())
        return new CudaCalcAmoebaHarmonicAngleForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaHarmonicInPlaneAngleForceKernel::Name())
        return new CudaCalcAmoebaHarmonicInPlaneAngleForceKernel(name, platform, cu, context.getSystem());

//    if (name == CalcAmoebaTorsionForceKernel::Name())
//        return new CudaCalcAmoebaTorsionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaPiTorsionForceKernel::Name())
        return new CudaCalcAmoebaPiTorsionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaStretchBendForceKernel::Name())
        return new CudaCalcAmoebaStretchBendForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaOutOfPlaneBendForceKernel::Name())
        return new CudaCalcAmoebaOutOfPlaneBendForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaTorsionTorsionForceKernel::Name())
        return new CudaCalcAmoebaTorsionTorsionForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaMultipoleForceKernel::Name())
        return new CudaCalcAmoebaMultipoleForceKernel(name, platform, cu, context.getSystem());

//    if (name == CalcAmoebaGeneralizedKirkwoodForceKernel::Name())
//        return new CudaCalcAmoebaGeneralizedKirkwoodForceKernel(name, platform, cu, context.getSystem());

    if (name == CalcAmoebaVdwForceKernel::Name())
        return new CudaCalcAmoebaVdwForceKernel(name, platform, cu, context.getSystem());

//    if (name == CalcAmoebaWcaDispersionForceKernel::Name())
//        return new CudaCalcAmoebaWcaDispersionForceKernel(name, platform, cu, context.getSystem());
//
//    if (name == CalcAmoebaUreyBradleyForceKernel::Name())
//        return new CudaCalcAmoebaUreyBradleyForceKernel(name, platform, cu, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
