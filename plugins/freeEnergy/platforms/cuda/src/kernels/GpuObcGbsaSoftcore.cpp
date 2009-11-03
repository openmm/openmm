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

#include "GpuObcGbsaSoftcore.h"
#include "GpuFreeEnergyCudaKernels.h"

// GpuObcGbsaSoftcore constructor

GpuObcGbsaSoftcore::GpuObcGbsaSoftcore( ){
    _psNonPolarScalingFactors   = NULL;
}

GpuObcGbsaSoftcore::~GpuObcGbsaSoftcore( ){
    delete _psNonPolarScalingFactors;
}

// initialize NonPolarScalingFactors array

int GpuObcGbsaSoftcore::initializeNonPolarScalingFactors( unsigned int numberOfParticles ){
    _psNonPolarScalingFactors = new CUDAStream<float>( numberOfParticles, 1, "ObcSoftcoreNonPolarScalingFactors");
    for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
        (*_psNonPolarScalingFactors)[ii] = 1.0f;
    }
    return 0;
}

// set entry in NonPolarScalingFactors array

int GpuObcGbsaSoftcore::setNonPolarScalingFactors( unsigned int particleIndex, float nonPolarScalingFactor ){
    (*_psNonPolarScalingFactors)[particleIndex] = nonPolarScalingFactor;
    return 0;
}

// upload NonPolarScalingFactors array

int GpuObcGbsaSoftcore::upload( gpuContext gpu ){
    _psNonPolarScalingFactors->Upload();
    SetCalculateObcGbsaSoftcoreNonPolarScalingFactorsSim( getGpuNonPolarScalingFactors() );
    return 0;
}

// get address for NonPolarScalingFactors array on board

float* GpuObcGbsaSoftcore::getGpuNonPolarScalingFactors( void ) const {
    return _psNonPolarScalingFactors->_pDevStream[0];
}
