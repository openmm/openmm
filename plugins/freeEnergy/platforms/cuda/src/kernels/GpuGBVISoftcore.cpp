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

#include "GpuGBVISoftcore.h"
#include "GpuFreeEnergyCudaKernels.h"

// GpuGBVISoftcore constructor

GpuGBVISoftcore::GpuGBVISoftcore( ){
    _bornRadiiScalingMethod   = 0;
    _quinticLowerLimitFactor  = 0.8;
    _quinticUpperLimit        = 0.008;
    _psSwitchDerivative       = NULL;
}

// GpuGBVISoftcore destructor

GpuGBVISoftcore::~GpuGBVISoftcore( ){
    delete _psSwitchDerivative;
}

// set quintic lower limit factor value

int GpuGBVISoftcore::setQuinticLowerLimitFactor( float inputQuinticLowerLimitFactor ){
    _quinticLowerLimitFactor = inputQuinticLowerLimitFactor;
    return 0;
}

// get quintic lower limit factor value

float GpuGBVISoftcore::getQuinticLowerLimitFactor( void ) const {
    return _quinticLowerLimitFactor;
}

// set quintic upper limit value

int GpuGBVISoftcore::setQuinticUpperLimit( float inputQuinticUpperLimit ){
    _quinticUpperLimit = inputQuinticUpperLimit;
    return 0;
}

// get quintic upper limit value

float GpuGBVISoftcore::getQuinticUpperLimit( void ) const {
    return _quinticUpperLimit;
}

// get Born radii scaling method

int GpuGBVISoftcore::getBornRadiiScalingMethod( void ) const {
    return _bornRadiiScalingMethod;
}

// set Born radii scaling method

int GpuGBVISoftcore::setBornRadiiScalingMethod( int inputBornRadiiScalingMethod ){
    _bornRadiiScalingMethod = inputBornRadiiScalingMethod;
    return 0;
}

// get address for SwitchDerivative array on board

float* GpuGBVISoftcore::getGpuSwitchDerivative( void ) const {
    return _psSwitchDerivative->_pDevStream[0];
}

// get SwitchDerivative array 

CUDAStream<float>* GpuGBVISoftcore::getSwitchDerivative( void ) const {
    return _psSwitchDerivative;
}

// initialize SwitchDerivative array

int GpuGBVISoftcore::initializeGpuSwitchDerivative( unsigned int numberOfParticles ){
    _psSwitchDerivative = new CUDAStream<float>( numberOfParticles, 1, "SwitchDerivative");
    for( unsigned int ii = 0; ii < numberOfParticles; ii++ ){
        (*_psSwitchDerivative)[ii] = 1.0f;
    }   
    return 0;
}

// upload SoftCoreLambda array

int GpuGBVISoftcore::upload( gpuContext gpu ){
    if( getBornRadiiScalingMethod() > 0 ){
        SetCalculateGBVISoftcoreSupplementarySim( this );
    }
    return 0;
}
