#ifndef OPENMM_FREE_ENERGY_GPU_GBVI_SOFTCORE_
#define OPENMM_FREE_ENERGY_GPU_GBVI_SOFTCORE_
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

#include "gputypes.h"
#include "cudatypes.h"
//#include "cudaKernels.h"
#include "openmm/OpenMMException.h"

// info related to nonbonded softcore

class GpuGBVISoftcore {

    public:

   /** 
     * This is an enumeration of the different methods that may be used for scaling of the Born radii.
     */
        /**
         * No scaling method is applied.
         */
        static const int NoScaling          = 0;
        /**
         * Use the method outlined in Proteins 55, 383-394 (2004), Eq. 6
         */
        static const int Tanh               = 1;
        /**
         * Use quintic spline scaling function
         */
        static const int QuinticSpline      = 2;

        GpuGBVISoftcore();
        ~GpuGBVISoftcore();

        /** 
         * Set softcore value
         */

        int setSoftCoreLambda( float softCoreLambda );

        /** 
         * Get softcore value
         */

        float getSoftCoreLambda( void ) const;

        /** 
         * Set quintic lower limit factor value
         */

        int setQuinticLowerLimitFactor( float quinticLowerLimitFactor );

        /** 
         * Get quintic lower limit factor value
         */

        float getQuinticLowerLimitFactor( void ) const;

        /** 
         * Set quintic upper limit value
         */

        int setQuinticUpperLimit( float quinticUpperLimit );

        /** 
         * Get quintic upper limit value
         */

        float getQuinticUpperLimit( void ) const;

        /** 
         * Get Born radii scaling method
         */

        int getBornRadiiScalingMethod( void ) const;

        /** 
         * Set Born radii scaling method
         */

        int setBornRadiiScalingMethod( int bornRadiiScalingMethod );

        // initialize SoftCoreLJLambda particle array

        int initializeGpuSwitchDerivative( unsigned int numberOfParticles );

       /** 
         * Get address for switch derivative array
         * 
         * @return address
         */

         float* getGpuSwitchDerivative( void ) const;

       /** 
         * Get switch derivative array
         * 
         * @return address
         */

         CUDAStream<float>* getSwitchDerivative( void ) const;

        /** 
         * Upload data
         * 
         * @return 0 always
         */

        int upload( gpuContext gpu );

    private:
       float _quinticLowerLimitFactor;
       float _quinticUpperLimit;
       unsigned int _bornRadiiScalingMethod;
       CUDAStream<float>*  _psSwitchDerivative;
};

#endif // OPENMM_FREE_ENERGY_GPU_GBVI_SOFTCORE_
