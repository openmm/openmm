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
 * This file contains the kernels for Verlet integration.  It is included
 * several times in kVerletUpdate.cu with different #defines to generate
 * different versions of the kernels.
 */

#ifdef REMOVE_CM
__global__ void kVerletUpdatePart1CM_kernel()
#else
__global__ void kVerletUpdatePart1_kernel()
#endif
{
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
#ifdef REMOVE_CM
    extern __shared__ float3 sCM[];
    float3 CM           = { 0.0f, 0.0f, 0.0f};
    float4 CM1          = { 0.0f, 0.0f, 0.0f, 0.0f };

    // Read CM outputs from previous step
    unsigned int cpos = threadIdx.x;
    while (cpos < gridDim.x)
    {
        CM1             = cSim.pLinearMomentum[cpos];
        CM.x           += CM1.x;
        CM.y           += CM1.y;
        CM.z           += CM1.z;
        cpos           += blockDim.x;
    }
    sCM[threadIdx.x].x  = CM.x;
    sCM[threadIdx.x].y  = CM.y;
    sCM[threadIdx.x].z  = CM.z;
    __syncthreads();

    // Reduce CM
    unsigned int offset = 1;
    unsigned int mask   = 1;
    while (offset < blockDim.x)
    {
        if (((threadIdx.x & mask) == 0) && (threadIdx.x + offset < blockDim.x))
        {
            sCM[threadIdx.x].x += sCM[threadIdx.x + offset].x;
            sCM[threadIdx.x].y += sCM[threadIdx.x + offset].y;
            sCM[threadIdx.x].z += sCM[threadIdx.x + offset].z;
        }
        mask = 2 * mask + 1;
        offset *= 2;
        __syncthreads();
    }
#endif
    while (pos < cSim.atoms)
    {
        float4 apos             = cSim.pPosq[pos];
        float4 velocity         = cSim.pVelm4[pos];
        float4 force            = cSim.pForce4[pos];
        float dtOverMass        = cSim.deltaT*velocity.w;

        cSim.pOldPosq[pos]      = apos;
        velocity.x             += dtOverMass*force.x;
        velocity.y             += dtOverMass*force.y;
        velocity.z             += dtOverMass*force.z;
#ifdef REMOVE_CM
        velocity.x             -= sCM[0].x;
        velocity.y             -= sCM[0].y;
        velocity.z             -= sCM[0].z;
#endif

        apos.x                  = velocity.x*cSim.deltaT;
        apos.y                  = velocity.y*cSim.deltaT;
        apos.z                  = velocity.z*cSim.deltaT;

        cSim.pPosqP[pos]        = apos;
        cSim.pVelm4[pos]        = velocity;
        pos                    += blockDim.x * gridDim.x;
    }
}

#ifdef REMOVE_CM
__global__ void kVerletUpdatePart2CM_kernel()
#else
__global__ void kVerletUpdatePart2_kernel()
#endif
{
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
#ifdef REMOVE_CM
    extern __shared__ float3 sCM[];
    float3 CM                   = {0.0f, 0.0f, 0.0f};
#endif

    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
        float4 apos             = cSim.pPosq[pos];
        float4 xPrime           = cSim.pPosqP[pos];

        velocity.x              = cSim.oneOverDeltaT*(xPrime.x);
        velocity.y              = cSim.oneOverDeltaT*(xPrime.y);
        velocity.z              = cSim.oneOverDeltaT*(xPrime.z);

        xPrime.x               += apos.x;
        xPrime.y               += apos.y;
        xPrime.z               += apos.z;

#ifdef REMOVE_CM
        float mass              = 1.0f / velocity.w;
        CM.x                   += mass * velocity.x;
        CM.y                   += mass * velocity.y;
        CM.z                   += mass * velocity.z;
#endif
        cSim.pPosq[pos]         = xPrime;
        cSim.pVelm4[pos]        = velocity;

        pos                    += blockDim.x * gridDim.x;
    }

#ifdef REMOVE_CM
    // Scale CM
    CM.x *= cSim.inverseTotalMass;
    CM.y *= cSim.inverseTotalMass;
    CM.z *= cSim.inverseTotalMass;
    sCM[threadIdx.x] = CM;
    __syncthreads();

    // Reduce CM for CTA
    unsigned int offset = 1;
    unsigned int mask   = 1;
    while (offset < blockDim.x)
    {
        if (((threadIdx.x & mask) == 0) && (threadIdx.x + offset < blockDim.x))
        {
            sCM[threadIdx.x].x += sCM[threadIdx.x + offset].x;
            sCM[threadIdx.x].y += sCM[threadIdx.x + offset].y;
            sCM[threadIdx.x].z += sCM[threadIdx.x + offset].z;
        }
        mask = 2 * mask + 1;
        offset *= 2;
        __syncthreads();
    }
    if (threadIdx.x == 0)
    {
        float4 CM;
        CM.x                                = sCM[0].x;
        CM.y                                = sCM[0].y;
        CM.z                                = sCM[0].z;
        CM.w                                = 0.0f;
        cSim.pLinearMomentum[blockIdx.x]    = CM;
    }
#endif
}
