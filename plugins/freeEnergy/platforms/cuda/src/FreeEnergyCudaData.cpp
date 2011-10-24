/* -------------------------------------------------------------------------- *
 *                               OpenMMFreeEnergy                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "FreeEnergyCudaData.h"
#include "openmm/OpenMMException.h"
#include <sstream>

extern "C" void removeFreeEnergyCudaDataFromContextMap( void* context ); 

namespace OpenMM {

FreeEnergyCudaData::FreeEnergyCudaData( CudaPlatform::PlatformData& data ) : cudaPlatformData(data) {
        
    kernelCount                   = 0;
    freeEnergyGpu                 = freeEnergyGpuInit( cudaPlatformData.gpu );

    log                           = NULL;
    contextImpl                   = NULL;
    gpuInitialized                = false;

    boxDimensions[0]              = 0.0;
    boxDimensions[1]              = 0.0;
    boxDimensions[2]              = 0.0;
}   

FreeEnergyCudaData::~FreeEnergyCudaData() {
    if( getLog() ){
        (void) fprintf( getLog(), "~FreeEnergyCudaData called kernelCount=%d\n", kernelCount );
        (void) fflush( getLog() );
    }   
    freeEnergyGpuShutDown( freeEnergyGpu );
}

void FreeEnergyCudaData::decrementKernelCount( void ) {

    kernelCount--;
    if( getLog() ){
        (void) fprintf( getLog(), "~reeEnergyCudaData decrementKernelCount called. %d\n", kernelCount );
        (void) fflush( getLog() );
    }   
    if( kernelCount == 0 && contextImpl != NULL ){
        removeFreeEnergyCudaDataFromContextMap( contextImpl );
        freeEnergyGpuShutDown( freeEnergyGpu );
    }
}

void FreeEnergyCudaData::incrementKernelCount( void ) {
    kernelCount++;
}

freeEnergyGpuContext FreeEnergyCudaData::getFreeEnergyGpu( void ) const {
    return freeEnergyGpu;
}

void FreeEnergyCudaData::setLog( FILE* inputLog ) {
    log            = inputLog;
    freeEnergyGpu->log = inputLog;
}

FILE* FreeEnergyCudaData::getLog( void ) const {
    return log;
}

void FreeEnergyCudaData::setContextImpl( void* inputContextImpl ) {
    contextImpl = inputContextImpl;
}

void FreeEnergyCudaData::initializeGpu( void ) {

    if( !gpuInitialized ){

        gpuContext gpu = freeEnergyGpu->gpuContext;

        if( freeEnergyGpu->freeEnergySim.nonbondedCutoff != gpu->sim.nonbondedCutoff ){
            std::stringstream msg;
            msg << "The softcore non-bonded cutoff=" << freeEnergyGpu->freeEnergySim.nonbondedCutoff;
            msg << "does not agree with the non-softcore cutoff= " << gpu->sim.nonbondedCutoff;
            throw OpenMM::OpenMMException( msg.str() );
        }
/*
        freeEnergyGpuBuildOutputBuffers( freeEnergyGpu, getHasFreeEnergyGeneralizedKirkwood() );
        freeEnergyGpuBuildThreadBlockWorkList( freeEnergyGpu );

        boxDimensions[0] = freeEnergyGpu->gpuContext->sim.periodicBoxSizeX;
        boxDimensions[1] = freeEnergyGpu->gpuContext->sim.periodicBoxSizeY;
        boxDimensions[2] = freeEnergyGpu->gpuContext->sim.periodicBoxSizeZ;
*/
        gpuBuildExclusionList( gpu );
        gpuSetConstants( gpu );
        freeEnergyGpuSetConstants( freeEnergyGpu );
        gpuInitialized   = true;

        if( log ){
            //gpuPrintCudaFreeEnergyGmxSimulation( freeEnergyGpu, getLog() );
            (void) fprintf( log, "FreeEnergyCudaGpu initialized kernelCount=%d\n", kernelCount );
            (void) fflush( log );
        }

    } else {
/*
        if( boxDimensions[0] != freeEnergyGpu->gpuContext->sim.periodicBoxSizeX ||
            boxDimensions[1] != freeEnergyGpu->gpuContext->sim.periodicBoxSizeY ||
            boxDimensions[2] != freeEnergyGpu->gpuContext->sim.periodicBoxSizeZ ){
            freeEnergyGpuSetConstants( freeEnergyGpu, 1 );
            
            boxDimensions[0] = freeEnergyGpu->gpuContext->sim.periodicBoxSizeX;
            boxDimensions[1] = freeEnergyGpu->gpuContext->sim.periodicBoxSizeY;
            boxDimensions[2] = freeEnergyGpu->gpuContext->sim.periodicBoxSizeZ;
        }
*/

    }

    return;
}

}

