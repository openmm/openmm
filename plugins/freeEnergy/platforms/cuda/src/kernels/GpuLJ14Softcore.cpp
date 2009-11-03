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

#include "GpuLJ14Softcore.h"
#include "GpuFreeEnergyCudaKernels.h"

// GpuLJ14Softcore constructor

GpuLJ14Softcore::GpuLJ14Softcore( ){
    psLJ14SoftcoreID          = NULL;
    psLJ14SoftcoreParameter   = NULL;
}

// GpuLJ14Softcore destructor

GpuLJ14Softcore::~GpuLJ14Softcore( ){
    delete psLJ14SoftcoreID;
    delete psLJ14SoftcoreParameter;
}

int GpuLJ14Softcore::flipStrides( gpuContext gpu ){
    int flip = gpu->sim.outputBuffers - 1;
    for (unsigned int ii = 0; ii < psLJ14SoftcoreID->_stride; ii++)
    {
        (*psLJ14SoftcoreID)[ii].z = flip - (*psLJ14SoftcoreID)[ii].z;
        (*psLJ14SoftcoreID)[ii].w = flip - (*psLJ14SoftcoreID)[ii].w;
    }   
    psLJ14SoftcoreID->Upload();

    return 0;
}

