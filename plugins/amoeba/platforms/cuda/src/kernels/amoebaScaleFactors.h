#ifndef __AMOEBA_SCALE_FACTORS_H__
#define __AMOEBA_SCALE_FACTORS_H__

/* -------------------------------------------------------------------------- *
 *                          AmoebaOpenMM                                      *
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

static __constant__ float mpoleScale[5]   = { 0.0f, 0.0f, 0.0f, 0.4f, 0.8f };
static __constant__ float polarScale[5]   = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f };
static __constant__ float directScale[5]  = { 0.0f, 1.0f, 1.0f, 1.0f, 1.0f };
//float mutualScale[5]  = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
           
// must be explicitly initialized!

//static __constant__ float mScale[4]       = { 0.0f, 0.4f, 0.8f, 1.0f };
//static __constant__ float pScale[4]       = { 1.0f, 0.5f, 0.0f, -2.0f };
//static __constant__ float dScale[2]       = { 0.0f, 1.0f };
//static __constant__ float uScale[5]       = { 1.0f, 1.0f, 1.0f, 1.0f, 1.0f };
           
// subroutine to get masked scale factors

__device__ static void getMaskedDScaleFactor( unsigned int gridIndex, int scaleMask, float* dScale )
{
    unsigned int mask             = 1 << gridIndex;   
    *dScale                       = (scaleMask & mask) ? 0.0f : 1.0f;

}

__device__ static void getMaskedPScaleFactor( unsigned int gridIndex, int2 scaleMask, float* pScale )
{
    unsigned int mask             = 1 << gridIndex;   
    *pScale                       = (scaleMask.x & mask) ? 0.5f : 1.0f;
    *pScale                      *= (scaleMask.y & mask) ? 0.0f : 1.0f;

}

__device__ static void getMaskedMScaleFactor( unsigned int gridIndex, int2 scaleMask, float* mScale )
{
    unsigned int mask             = 1 << gridIndex;   

    // 0 0 -> 1 -> 1   -> 1.0
    // 1 0 -> 1 -> 0.4 -> 0.4
    // 0 1 -> 1 -> 0.8 -> 0.8
    // 1 1 -> 0 ->   0 -> 0.0

    *mScale                       = (scaleMask.x & mask) && (scaleMask.y & mask) ? 0.0f : 1.0f;
    *mScale                      *= (scaleMask.x & mask) ? 0.8f : 1.0f;
    *mScale                      *= (scaleMask.y & mask) ? 0.4f : 1.0f;

}

// subroutine to get cell coordinates

__device__ static void decodeCell( unsigned int cellId, unsigned int* x, unsigned int* y, bool* exclusions )
{
    *x          = cellId;
    *y          = ((*x >> 2) & 0x7fff) << GRIDBITS;

    *exclusions = (*x & 0x1);
    *x          = (*x >> 17) << GRIDBITS;

}

__device__ static void load3dArrayBufferPerWarp( unsigned int offset, float* forceSum, float* outputForce )
{

    float of; 
    of                                  = outputForce[offset];
    of                                 += forceSum[0];
    outputForce[offset]                 = of;  

    of                                  = outputForce[offset+1];
    of                                 += forceSum[1];
    outputForce[offset+1]               = of;  

    of                                  = outputForce[offset+2];
    of                                 += forceSum[2];
    outputForce[offset+2]               = of;  

}

__device__ static void load3dArray( unsigned int offset, float* forceSum, float* outputForce )
{

    outputForce[offset]                 = forceSum[0];  
    outputForce[offset+1]               = forceSum[1];  
    outputForce[offset+2]               = forceSum[2];  

}

__device__ static void scale3dArray( float scaleFactor, float* force )
{

    force[0]  *= scaleFactor;  
    force[1]  *= scaleFactor;  
    force[2]  *= scaleFactor;  

}

#endif //__AMOEBA_SCALE_FACTORS_H__
