#ifndef FREE_ENERGY_CUDA_TYPES_H
#define FREE_ENERGY_CUDA_TYPES_H

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

#include <kernels/cudatypes.h>

#include <stdarg.h>
#include <limits>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>
#include <builtin_types.h>
#include <vector_functions.h>

enum CudaFreeEnergyNonbondedMethod
{
    FREE_ENERGY_NO_CUTOFF,
    FREE_ENERGY_CUTOFF,
    FREE_ENERGY_PERIODIC
};

struct cudaFreeEnergyGmxSimulation {

    // Constants

    unsigned int   LJ14_count;                      // LJ count
    int4*          pLJ14ID;                         // LJ 14 particles ids
    float4*        pLJ14Parameter;                  // LJ 14 parameters 

    float           epsfac;                         // Epsilon factor for CDLJ calculations
    CudaFreeEnergyNonbondedMethod nonbondedMethod;  // How to handle nonbonded interactions
    float           nonbondedCutoff;                // Cutoff distance for nonbonded interactions
    float           nonbondedCutoffSqr;             // Square of the cutoff distance for nonbonded interactions

    float           periodicBoxSizeX;               // The X dimension of the periodic box
    float           periodicBoxSizeY;               // The Y dimension of the periodic box
    float           periodicBoxSizeZ;               // The Z dimension of the periodic box

    float           invPeriodicBoxSizeX;            // The 1 over the X dimension of the periodic box
    float           invPeriodicBoxSizeY;            // The 1 over the Y dimension of the periodic box
    float           invPeriodicBoxSizeZ;            // The 1 over the Z dimension of the periodic box

    float           recipBoxSizeX;                  // The X dimension of the reciprocal box for Ewald summation
    float           recipBoxSizeY;                  // The Y dimension of the reciprocal box for Ewald summation
    float           recipBoxSizeZ;                  // The Z dimension of the reciprocal box for Ewald summation

    float           cellVolume;                     // Ewald parameter alpha (a.k.a. kappa)

    float           reactionFieldK;                 // Constant for reaction field correction
    float           reactionFieldC;                 // Constant for reaction field correction
    float4*         pSigEps4;                       // sigma, eps, lambda. charge

    int             bornRadiiScalingMethod;         // flag for method to use scaling radii (0=none,1=quintic spline)
    float           quinticLowerLimitFactor;        // lower limit factor for quintic spline
    float           quinticUpperLimit;              // upper limit for quintic spline
    float*          pSwitchDerivative;              // switch deriviatives for quintic spline

    float*          pNonPolarScalingFactors;        // non-polar scaling factors

};

#endif // FREE_ENERGY_CUDA_TYPES_H
