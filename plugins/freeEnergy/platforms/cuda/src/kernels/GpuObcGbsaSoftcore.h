#ifndef OPENMM_FREE_ENERGY_GPU_OBC_GBSA_SOFTCORE_
#define OPENMM_FREE_ENERGY_GPU_OBC_GBSA_SOFTCORE_

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

// info related to nonbonded softcore

class GpuObcGbsaSoftcore {

    public:

        GpuObcGbsaSoftcore();
        ~GpuObcGbsaSoftcore();

        /** 
         * Initialize NonPolarScalingFactors array
         * 
         * @param numberOfParticles number of particles
         *
         * @return 0 always
         */

        int initializeNonPolarScalingFactors( unsigned int numberOfParticles );

        /** 
         * Upload data
         * 
         * @param implicitSolvent set if implicit solvent is included in system
         *
         * @return 0 always
         */

        int upload( gpuContext gpu );

        /** 
         * Set nonPolarScalingFactor entry
         * 
         * @param particleIndex             index of particle
         * @param nonPolarScalingFactor     nonPolarScalingFactor value
         *
         * @return 0 always
         */

        int setNonPolarScalingFactors( unsigned int particleIndex, float nonPolarScalingFactor );

        /** 
         * Get address for NonPolarScalingFactors array on board
         * 
         * @return address
         */

         float* getGpuNonPolarScalingFactors( void ) const;

    private:

       CUDAStream<float>*  _psNonPolarScalingFactors;
    
};

#endif // OPENMM_FREE_ENERGY_GPU_OBC_GBSA_SOFTCORE_
