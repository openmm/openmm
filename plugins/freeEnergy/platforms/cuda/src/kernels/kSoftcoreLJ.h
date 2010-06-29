#ifndef _K_SOFTCORE_LJ__H__
#define _K_SOFTCORE_LJ__H__

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

/**
 * This file contains kernel for calculating softcore LJ force prefactor
 */

#ifdef USE_SOFTCORE_LJ
static __device__ float getSoftCoreLJ( float r2, float sig, float  eps, float lambdaI, float lambdaJ, float* energy)
{

   float r                         = sqrt(r2);
   float lambda                    = lambdaI < lambdaJ ? lambdaI : lambdaJ;
   eps                            *= lambda;


    // (r/sig)
    float sig2                     = r/sig;
          sig2                    *= sig2;
    float sig6                     = sig2*sig2*sig2;

    float softcoreLJTerm           = 0.5f*( 1.0f -  lambda) + sig6;
    float softcoreLJInv            = 1.0f/softcoreLJTerm;
    float softcoreLJInv2           = softcoreLJInv*softcoreLJInv;
    *energy                        = eps*(softcoreLJInv2 - softcoreLJInv);

    return eps*softcoreLJInv2*( 12.0f*softcoreLJInv - 6.0f )*sig6;
    
}
#endif

#endif
