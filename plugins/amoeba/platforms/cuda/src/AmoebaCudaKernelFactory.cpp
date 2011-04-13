/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors:                                                                   *
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
#include "AmoebaCudaData.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" void registerPlatforms() {
}

extern "C" OPENMMCUDA_EXPORT void registerKernelFactories() {
//fprintf( stderr,"In registerKernelFactories AmoebaCudaKernelFactory\n" ); fflush( stderr );
    for( int ii = 0; ii < Platform::getNumPlatforms(); ii++ ){
        Platform& platform = Platform::getPlatform(ii);
        if( platform.getName() == "Cuda" ){

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
    }
}

static std::map<ContextImpl*, AmoebaCudaData*> contextToAmoebaDataMap;

// look up AmoebaCudaData for input contextImpl in contextToAmoebaDataMap

extern "C" void* getAmoebaCudaData( ContextImpl& context ) {
    std::map<ContextImpl*, AmoebaCudaData*>::const_iterator mapIterator  = contextToAmoebaDataMap.find(&context);
    if( mapIterator == contextToAmoebaDataMap.end() ){
        return NULL;
    } else {
        return static_cast<void*>(mapIterator->second);
    }
}

// remove AmoebaCudaData from contextToAmoebaDataMap

extern "C" void removeAmoebaCudaDataFromContextMap( void* inputContext ) {
    ContextImpl* context = static_cast<ContextImpl*>(inputContext);
    contextToAmoebaDataMap.erase( context );
    return;
}

KernelImpl* AmoebaCudaKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaPlatform::PlatformData& cudaPlatformData = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());

    // create AmoebaCudaData object if contextToAmoebaDataMap does not contain
    // key equal to current context

    AmoebaCudaData* amoebaCudaData;
    std::map<ContextImpl*, AmoebaCudaData*>::const_iterator mapIterator  = contextToAmoebaDataMap.find(&context);
    if( mapIterator == contextToAmoebaDataMap.end() ){
        amoebaCudaData                         = new AmoebaCudaData( cudaPlatformData );
        contextToAmoebaDataMap[&context]       = amoebaCudaData;
        //amoebaCudaData->setLog( stderr );
        amoebaCudaData->setContextImpl( static_cast<void*>(&context) );
    } else {
        amoebaCudaData                         = mapIterator->second;
    }

    if (name == CalcAmoebaHarmonicBondForceKernel::Name())
        return new CudaCalcAmoebaHarmonicBondForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaHarmonicAngleForceKernel::Name())
        return new CudaCalcAmoebaHarmonicAngleForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaHarmonicInPlaneAngleForceKernel::Name())
        return new CudaCalcAmoebaHarmonicInPlaneAngleForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaTorsionForceKernel::Name())
        return new CudaCalcAmoebaTorsionForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaPiTorsionForceKernel::Name())
        return new CudaCalcAmoebaPiTorsionForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaStretchBendForceKernel::Name())
        return new CudaCalcAmoebaStretchBendForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaOutOfPlaneBendForceKernel::Name())
        return new CudaCalcAmoebaOutOfPlaneBendForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaTorsionTorsionForceKernel::Name())
        return new CudaCalcAmoebaTorsionTorsionForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaMultipoleForceKernel::Name())
        return new CudaCalcAmoebaMultipoleForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaGeneralizedKirkwoodForceKernel::Name())
        return new CudaCalcAmoebaGeneralizedKirkwoodForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaVdwForceKernel::Name())
        return new CudaCalcAmoebaVdwForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaWcaDispersionForceKernel::Name())
        return new CudaCalcAmoebaWcaDispersionForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    if (name == CalcAmoebaUreyBradleyForceKernel::Name())
        return new CudaCalcAmoebaUreyBradleyForceKernel(name, platform, *amoebaCudaData, context.getSystem());

    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
