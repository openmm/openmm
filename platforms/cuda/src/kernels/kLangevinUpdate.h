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
 * This file contains the kernels for Langevin integration.  It is included
 * several times in kLangevinUpdate.cu with different #defines to generate
 * different versions of the kernels.
 */

#ifdef REMOVE_CM
__global__ void kLangevinUpdatePart1CM_kernel()
#else
__global__ void kLangevinUpdatePart1_kernel()
#endif
{
    __shared__ float params[MaxParams];
    if (threadIdx.x < MaxParams)
        params[threadIdx.x] = cSim.pLangevinParameters[threadIdx.x];
    __syncthreads();
    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos   = cSim.pRandomPosition[blockIdx.x];
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
        float4 velocity         = cSim.pVelm4[pos];
        float4 xVector          = cSim.pxVector4[pos];
        float4 random4a         = cSim.pRandom4a[rpos + pos];
        float2 random2a         = cSim.pRandom2a[rpos + pos];
        float4 apos             = cSim.pPosq[pos];
        float4 force            = cSim.pForce4[pos];

        float3 Vmh;
        float sqrtInvMass       = sqrt(velocity.w);
        Vmh.x                   = xVector.x * params[DOverTauC] + sqrtInvMass * params[Yv] * random4a.x;
        Vmh.y                   = xVector.y * params[DOverTauC] + sqrtInvMass * params[Yv] * random4a.y;
        Vmh.z                   = xVector.z * params[DOverTauC] + sqrtInvMass * params[Yv] * random4a.z;
        float4 vVector;
        vVector.x               = sqrtInvMass * params[V] * random4a.w;
        vVector.y               = sqrtInvMass * params[V] * random2a.x;
        vVector.z               = sqrtInvMass * params[V] * random2a.y;
        vVector.w               = 0.0f;
        cSim.pvVector4[pos]     = vVector;
        velocity.x              = velocity.x * params[EM_V] +
                                  velocity.w * force.x * params[TauOneMinusEM_V] +
                                  vVector.x -
                                  params[EM] * Vmh.x;
        velocity.y              = velocity.y * params[EM_V] +
                                  velocity.w * force.y * params[TauOneMinusEM_V] +
                                  vVector.y -
                                  params[EM] * Vmh.y;
        velocity.z              = velocity.z * params[EM_V] +
                                  velocity.w * force.z * params[TauOneMinusEM_V] +
                                  vVector.z -
                                  params[EM] * Vmh.z;
#ifdef REMOVE_CM
        velocity.x             -= sCM[0].x;
        velocity.y             -= sCM[0].y;
        velocity.z             -= sCM[0].z;
#endif
        cSim.pOldPosq[pos]      = apos;
        apos.x                  = velocity.x * params[Fix1];
        apos.y                  = velocity.y * params[Fix1];
        apos.z                  = velocity.z * params[Fix1];
        cSim.pPosqP[pos]        = apos;
        cSim.pVelm4[pos]        = velocity;
        pos                    += blockDim.x * gridDim.x;
    }
}

#ifdef REMOVE_CM
__global__ void kLangevinUpdatePart2CM_kernel()
#else
__global__ void kLangevinUpdatePart2_kernel()
#endif
{
    __shared__ float params[MaxParams];
    if (threadIdx.x < MaxParams)
        params[threadIdx.x] = cSim.pLangevinParameters[threadIdx.x];
    __syncthreads();
    unsigned int pos            = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int rpos           = cSim.pRandomPosition[blockIdx.x];
#ifdef REMOVE_CM
    extern __shared__ float3 sCM[];
    float3 CM                   = {0.0f, 0.0f, 0.0f};
    __syncthreads();
#endif

    while (pos < cSim.atoms)
    {
        float4 velocity         = cSim.pVelm4[pos];
        float4 xPrime           = cSim.pPosqP[pos];
        float4 vVector          = cSim.pvVector4[pos];
        float4 xVector;
        float4 random4b         = cSim.pRandom4b[rpos + pos];
        float2 random2b         = cSim.pRandom2b[rpos + pos];
        float3 Xmh;
        float sqrtInvMass       = sqrt(velocity.w);
        velocity.x              = xPrime.x * params[OneOverFix1];
        velocity.y              = xPrime.y * params[OneOverFix1];
        velocity.z              = xPrime.z * params[OneOverFix1];
#ifdef REMOVE_CM
        float mass              = 1.0f / velocity.w;
        CM.x                   += mass * velocity.x;
        CM.y                   += mass * velocity.y;
        CM.z                   += mass * velocity.z;
#endif

        Xmh.x                   = vVector.x * params[TauDOverEMMinusOne] +
                                  sqrtInvMass * params[Yx] * random4b.x;
        Xmh.y                   = vVector.y * params[TauDOverEMMinusOne] +
                                  sqrtInvMass * params[Yx] * random4b.y;
        Xmh.z                   = vVector.z * params[TauDOverEMMinusOne] +
                                  sqrtInvMass * params[Yx] * random4b.z;
        xVector.x               = sqrtInvMass * params[X] * random4b.w;
        xVector.y               = sqrtInvMass * params[X] * random2b.x;
        xVector.z               = sqrtInvMass * params[X] * random2b.y;
        xPrime.x               += xVector.x - Xmh.x;
        xPrime.y               += xVector.y - Xmh.y;
        xPrime.z               += xVector.z - Xmh.z;


        cSim.pPosq[pos]         = xPrime;
        cSim.pVelm4[pos]        = velocity;
        cSim.pxVector4[pos]     = xVector;

        pos                    += blockDim.x * gridDim.x;
    }

    // Update random position pointer
    if (threadIdx.x == 0)
    {
        rpos                   += cSim.paddedNumberOfAtoms;
        if (rpos > cSim.randoms)
            rpos               -= cSim.randoms;
        cSim.pRandomPosition[blockIdx.x] = rpos;
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

