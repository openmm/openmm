/* -------------------------------------------------------------------------- *
 *                               OpenMMAmoeba                                 *
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

#include "AmoebaCudaData.h"
#include "openmm/OpenMMException.h"

extern "C" void removeAmoebaCudaDataFromContextMap( void* context ); 

namespace OpenMM {

AmoebaCudaData::AmoebaCudaData( CudaPlatform::PlatformData& data ) : cudaPlatformData(data) {
        
    kernelCount                   = 0;
    hasAmoebaBonds                = false;
    hasAmoebaGeneralizedKirkwood  = false;
    hasAmoebaMultipole            = false;
    amoebaGpu                     = amoebaGpuInit( cudaPlatformData.gpu );
    localForceKernel              = NULL;
    log                           = NULL;
    contextImpl                   = NULL;
    gpuInitialized                = false;
    applyMultipoleCutoff          = 0;
    useVdwNeighborList            = 0;
    multipoleForceCount           = 0;
}   

AmoebaCudaData::~AmoebaCudaData() {
    amoebaGpuShutDown( amoebaGpu );
}

void AmoebaCudaData::decrementKernelCount( void ) {
    kernelCount--;
    if( kernelCount == 0 && contextImpl != NULL ){
        removeAmoebaCudaDataFromContextMap( contextImpl );
    }
}

void AmoebaCudaData::incrementKernelCount( void ) {
    kernelCount++;
}

void AmoebaCudaData::setHasAmoebaBonds( bool inputHasAmoebaBonds ) {
    hasAmoebaBonds = inputHasAmoebaBonds;
}

void AmoebaCudaData::setHasAmoebaMultipole( bool inputHasAmoebaMultipole ) {
    hasAmoebaMultipole = inputHasAmoebaMultipole;
}

bool AmoebaCudaData::getHasAmoebaMultipole( void ) const {
    return hasAmoebaMultipole;
}

void AmoebaCudaData::setHasAmoebaGeneralizedKirkwood( bool inputHasAmoebaGeneralizedKirkwood ) {
    hasAmoebaGeneralizedKirkwood = inputHasAmoebaGeneralizedKirkwood;
}

bool AmoebaCudaData::getHasAmoebaGeneralizedKirkwood( void ) const {
    return hasAmoebaGeneralizedKirkwood;
}

amoebaGpuContext AmoebaCudaData::getAmoebaGpu( void ) const {
    return amoebaGpu;
}

void AmoebaCudaData::setAmoebaLocalForcesKernel( KernelImpl* inputLocalForceKernel ){
    localForceKernel = inputLocalForceKernel;
}

KernelImpl* AmoebaCudaData::getAmoebaLocalForcesKernel( void ) const {
    return localForceKernel;
}

void AmoebaCudaData::setLog( FILE* inputLog ) {
    log            = inputLog;
    amoebaGpu->log = inputLog;
}

FILE* AmoebaCudaData::getLog( void ) const {
    return log;
}

void AmoebaCudaData::setContextImpl( void* inputContextImpl ) {
    contextImpl = inputContextImpl;
}

void AmoebaCudaData::initializeGpu( void ) {
    if( !gpuInitialized ){
        if( getHasAmoebaGeneralizedKirkwood() && !getHasAmoebaMultipole() ){
            throw OpenMMException("GK force requires Multipole force\n");
        }
        amoebaGpuSetConstants( amoebaGpu );
        gpuInitialized = true;
        if( log ){
            gpuPrintCudaAmoebaGmxSimulation( amoebaGpu, getLog() );
            (void) fprintf( log, "Gpu initialized\n" );
            (void) fflush( log );
        }
    }
    return;
}

void AmoebaCudaData::incrementMultipoleForceCount( void ) {
    multipoleForceCount++;
}

int AmoebaCudaData::getMultipoleForceCount( void ) const {
    return multipoleForceCount;
}

void AmoebaCudaData::setApplyMultipoleCutoff( int inputApplyMultipoleCutoff ) {
    applyMultipoleCutoff = inputApplyMultipoleCutoff;
}

int AmoebaCudaData::getApplyMultipoleCutoff( void ) const {
    return applyMultipoleCutoff;
}

void AmoebaCudaData::setUseVdwNeighborList( int inputUseVdwNeighborList ) {
    useVdwNeighborList = inputUseVdwNeighborList;
}

int AmoebaCudaData::getUseVdwNeighborList( void ) const {
    return useVdwNeighborList;
}

}

