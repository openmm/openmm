/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Scott Le Grand, Peter Eastman                                     *
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

using namespace std;

#include "GpuNonbondedSoftcore.h"
#include "GpuFreeEnergyCudaKernels.h"

// GpuNonbondedSoftcore constructor

GpuNonbondedSoftcore::GpuNonbondedSoftcore( ){
    _softcoreLJLambda     = 1.0f;
    _psSoftcoreLJLambda   = NULL;
}

GpuNonbondedSoftcore::~GpuNonbondedSoftcore( ){
    delete _psSoftcoreLJLambda;
}

// set global softCoreLJLambda

int GpuNonbondedSoftcore::setSoftCoreLJLambda( float softCoreLJLambda ){
    _softcoreLJLambda  = softCoreLJLambda;
    return 0;
}

// get global softCoreLJLambda

float GpuNonbondedSoftcore::getSoftCoreLJLambda( void ) const {
    return _softcoreLJLambda;
}

// initialize SoftCoreLJLambda particle array

int GpuNonbondedSoftcore::initializeParticleSoftCoreLJLambda( unsigned int numberOfParticles ){
    _psSoftcoreLJLambda = new CUDAStream<float>( numberOfParticles, 1, "SoftcoreLJLambda");
    for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
        (*_psSoftcoreLJLambda)[ii] = 1.0f;
    }
    return 0;
}

// set entry in SoftCoreLJLambda particle array

int GpuNonbondedSoftcore::setParticleSoftCoreLJLambda( unsigned int particleIndex, float softCoreLJLambda ){
    (*_psSoftcoreLJLambda)[particleIndex] = softCoreLJLambda;
    return 0;
}

// upload SoftCoreLJLambda array

int GpuNonbondedSoftcore::upload( gpuContext gpu ){
    _psSoftcoreLJLambda->Upload();
    SetCalculateCDLJSoftcoreSupplementarySim( getGpuParticleSoftCoreLJLambda() );
    SetCalculateCDLJObcGbsaSoftcoreSupplementary1Sim( getGpuParticleSoftCoreLJLambda() );
    return 0;
}

// get address for SoftCoreLJLambda particle array on board

float* GpuNonbondedSoftcore::getGpuParticleSoftCoreLJLambda( void ) const {
    return _psSoftcoreLJLambda->_pDevStream[0];
}
