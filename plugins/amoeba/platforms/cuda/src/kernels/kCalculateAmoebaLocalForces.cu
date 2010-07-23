/* -------------------------------------------------------------------------- *
 *                                  AmoebaOpenMM                              *
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

#include <stdio.h>
#include <cuda.h>
#include <vector_functions.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include "amoebaGpuTypes.h"

static __constant__ cudaGmxSimulation cSim;
static __constant__ cudaAmoebaGmxSimulation cAmoebaSim;

/* Cuda compiler on Windows does not recognized "static const float" values */
#define LOCAL_HACK_PI        3.1415926535897932384626433832795f
#define LOCAL_HACK_RADIAN   57.29577951308232088f
#define LOCAL_HACK_RADIAN_D 57.29577951308232088

#define DOT3(v1, v2) (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)

#define CROSS_PRODUCT(v1, v2, c) \
    c.x = v1.y * v2.z - v1.z * v2.y; \
    c.y = v1.z * v2.x - v1.x * v2.z; \
    c.z = v1.x * v2.y - v1.y * v2.x;


void SetCalculateAmoebaLocalForcesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status         = cudaMemcpyToSymbol(cSim, &gpu->sim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetCalculateAmoebaLocalForcesSim copy to cSim failed");
    status         = cudaMemcpyToSymbol(cAmoebaSim, &amoebaGpu->amoebaSim, sizeof(cudaAmoebaGmxSimulation));     
    RTERROR(status, "cudaMemcpyToSymbol: SetSim copy to cAmoebaSim failed");
}

void GetCalculateAmoebaLocalForcesSim(amoebaGpuContext amoebaGpu)
{
    cudaError_t status;
    gpuContext gpu = amoebaGpu->gpuContext;
    status = cudaMemcpyFromSymbol(&gpu->sim, cSim, sizeof(cudaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: GetCalculateAmoebaLocalForcesSim copy from cSim failed");
    status = cudaMemcpyFromSymbol(&amoebaGpu->amoebaSim, cAmoebaSim, sizeof(cudaAmoebaGmxSimulation));     
    RTERROR(status, "cudaMemcpyFromSymbol: GetCalculateAmoebaLocalForcesSim copy from cAmoebaSim failed");
}

// bicubic spline

__device__ void bicubic( float4 y, float4 y1i, float4 y2i, float4 y12i, float x1, float x1l, float x1u,
                         float x2, float x2l, float x2u, float* energyOut, float* dang1Out, float* dang2Out )
//                         float4* c0, float4* c1, float4* c2, float4* c3 )
{

   // c[0][j] = cl[0-3]
   // c[1][j] = cl[4-7]
   // c[2][j] = cl[8-11]
   // c[4][j] = cl[12-15]

   float c[4][4];
   float d1   = x1u - x1l;
   float d2   = x2u - x2l;
   float d12  = d1*d2;

   float4 y1;
   y1.x       = d1*y1i.x;
   y1.y       = d1*y1i.y;
   y1.z       = d1*y1i.z;
   y1.w       = d1*y1i.w;

   float4 y2;
   y2.x       = d2*y2i.x;
   y2.y       = d2*y2i.y;
   y2.z       = d2*y2i.z;
   y2.w       = d2*y2i.w;

   float4 y12;
   y12.x      = d12*y12i.x;
   y12.y      = d12*y12i.y;
   y12.z      = d12*y12i.z;
   y12.w      = d12*y12i.w;

   // 1    1    1.000 [0][0]
   c[0][0]    = y.x;

   // 2    9    1.000 [0][1]
   c[0][1]    = y2.x;

   // 3    1   -3.000 [0][2]
   // 3    4    3.000
   // 3    9   -2.000
   // 3   12   -1.000

   c[0][2]    =  3.0f*(y.w - y.x) - (2.0f*y2.x + y2.w);
   
   // 4    1    2.000 [0][3]
   // 4    4   -2.000
   // 4    9    1.000
   // 4   12    1.000
 
   c[0][3]    =  2.0f*(y.x - y.w) + y2.x + y2.w;

   // 5    5    1.000 [1][0]
   c[1][0]    =  y1.x;

   // 6   13    1.000 [1][1]
   c[1][1]    =  y12.x;

   // 7    5   -3.000 [1][2]
   // 7    8    3.000
   // 7   13   -2.000
   // 7   16   -1.000

   c[1][2]    =  3.0f*(y1.w - y1.x) - (2.0f*y12.x + y12.w);

   // 8    5    2.000 [1][3]
   // 8    8   -2.000
   // 8   13    1.000
   // 8   16    1.000

   c[1][3]    =  2.0f*(y1.x - y1.w) + y12.x + y12.w;

   // 9    1   -3.000 [2][0]
   // 9    2    3.000
   // 9    5   -2.000
   // 9    6   -1.000

   c[2][0]    =  3.0f*(y.y - y.x) - (2.0f*y1.x + y1.y);

   // 10    9   -3.000 [2][1]
   // 10   10    3.000
   // 10   13   -2.000
   // 10   14   -1.000
   c[2][1]    =  3.0f*(y2.y - y2.x) - (2.0f*y12.x + y12.y);

   // 11    1    9.000 [2][2]
   // 11    2   -9.000
   // 11    3    9.000
   // 11    4   -9.000

   // 11    5    6.000
   // 11    6    3.000
   // 11    7   -3.000
   // 11    8   -6.000

   // 11    9    6.000
   // 11   10   -6.000
   // 11   11   -3.000
   // 11   12    3.000

   // 11   13    4.000
   // 11   14    2.000
   // 11   15    1.000
   // 11   16    2.000
   c[2][2]    =  9.0f*(y.x - y.y + y.z - y.w) + 6.0f* y1.x  +  3.0f* y1.y  -  3.0f* y1.z  -  6.0f* y1.w +
                                                6.0f* y2.x  -  6.0f* y2.y  -  3.0f* y2.z  +  3.0f* y2.w +
                                                4.0f*y12.x  +  2.0f*y12.y  +       y12.z  +  2.0f*y12.w;

   // 12    1   -6.000 [2][3]
   // 12    2    6.000
   // 12    3   -6.000
   // 12    4    6.000

   // 12    5   -4.000
   // 12    6   -2.000
   // 12    7    2.000
   // 12    8    4.000

   // 12    9   -3.000
   // 12   10    3.000
   // 12   11    3.000
   // 12   12   -3.000

   // 12   13   -2.000
   // 12   14   -1.000
   // 12   15   -1.000
   // 12   16   -2.000

   c[2][3]    =  6.0f*(y.y - y.x + y.w - y.z) + -4.0f* y1.x  -  2.0f* y1.y  +  2.0f* y1.z  +  4.0f* y1.w +
                                                -3.0f* y2.x  +  3.0f* y2.y  +  3.0f* y2.z  -  3.0f* y2.w
                                                -2.0f*y12.x  -       y12.y  -       y12.z  -  2.0f*y12.w;

   // 13    1    2.000 [3][0]
   // 13    2   -2.000
   // 13    5    1.000
   // 13    6    1.000
   c[3][0]    =  2.0f*(y.x - y.y) + y1.x + y1.y;

   // 14    9    2.000 [3][1]
   // 14   10   -2.000
   // 14   13    1.000
   // 14   14    1.000
   c[3][1]    =  2.0f*(y2.x - y2.y) + y12.x + y12.y;

   // 15    1   -6.000 [3][2]
   // 15    2    6.000
   // 15    3   -6.000
   // 15    4    6.000

   // 15    5   -3.000
   // 15    6   -3.000
   // 15    7    3.000
   // 15    8    3.000

   // 15    9   -4.000
   // 15   10    4.000
   // 15   11    2.000
   // 15   12   -2.000

   // 15   13   -2.000
   // 15   14   -2.000
   // 15   15   -1.000
   // 15   16   -1.000

   c[3][2]    =  6.0f*( y.y -  y.x +  y.w -  y.z) +
                 3.0f*(y1.z + y1.w - y1.x - y1.y) +
                 2.0f*( 2.0f*(y2.y - y2.x) + y2.z - y2.w) +
                 -2.0f*(y12.x + y12.y) - y12.z - y12.w;

   // 16    1    4.000 [3][3]
   // 16    2   -4.000
   // 16    3    4.000
   // 16    4   -4.000

   // 16    5    2.000
   // 16    6    2.000
   // 16    7   -2.000
   // 16    8   -2.000

   // 16    9    2.000
   // 16   10   -2.000
   // 16   11   -2.000
   // 16   12    2.000

   // 16   13    1.000
   // 16   14    1.000
   // 16   15    1.000
   // 16   16    1.000

   c[3][3]    =  4.0f*(  y.x -    y.y  +   y.z -  y.w)  +
                 2.0f*( y1.x +   y1.y  -  y1.z -  y1.w) +
                 2.0f*( y2.x -   y2.y  -  y2.z +  y2.w) +
                       y12.x +  y12.y  + y12.z + y12.w;

   float      t = (x1-x1l) / (x1u-x1l);
   float      u = (x2-x2l) / (x2u-x2l);

   float energy =            ((c[3][3]*u + c[3][2])*u + c[3][1])*u + c[3][0];
         energy = t*energy + ((c[2][3]*u + c[2][2])*u + c[2][1])*u + c[2][0];
         energy = t*energy + ((c[1][3]*u + c[1][2])*u + c[1][1])*u + c[1][0];
         energy = t*energy + ((c[0][3]*u + c[0][2])*u + c[0][1])*u + c[0][0];

   float dang1  =           (3.0f*c[3][3]*t + 2.0f*c[2][3])*t + c[1][3];
         dang1  = u*dang1 + (3.0f*c[3][2]*t + 2.0f*c[2][2])*t + c[1][2];
         dang1  = u*dang1 + (3.0f*c[3][1]*t + 2.0f*c[2][1])*t + c[1][1];
         dang1  = u*dang1 + (3.0f*c[3][0]*t + 2.0f*c[2][0])*t + c[1][0];

   float dang2  =           (3.0f*c[3][3]*u + 2.0f*c[3][2])*u + c[3][1];
         dang2  = t*dang2 + (3.0f*c[2][3]*u + 2.0f*c[2][2])*u + c[2][1];
         dang2  = t*dang2 + (3.0f*c[1][3]*u + 2.0f*c[1][2])*u + c[1][1];
         dang2  = t*dang2 + (3.0f*c[0][3]*u + 2.0f*c[0][2])*u + c[0][1];

      dang1     = dang1 / (x1u-x1l);
      dang2     = dang2 / (x2u-x2l);

   *energyOut   = energy;
   *dang1Out    = dang1;
   *dang2Out    = dang2;

}
    
__global__ 
#if (__CUDA_ARCH__ >= 200)
__launch_bounds__(GF1XX_LOCALFORCES_THREADS_PER_BLOCK, 1)
#elif (__CUDA_ARCH__ >= 130)
__launch_bounds__(GT2XX_LOCALFORCES_THREADS_PER_BLOCK, 1)
#else
__launch_bounds__(G8X_LOCALFORCES_THREADS_PER_BLOCK, 1)
#endif
void kCalculateAmoebaLocalForces_kernel()
{
    unsigned int pos = blockIdx.x * blockDim.x + threadIdx.x;
    //Vectors* A = &sV[threadIdx.x];

    float energy = 0.0f;

    while (pos < cAmoebaSim.amoebaBond_offset)
    {
        if (pos < cAmoebaSim.amoebaBonds)
        {
            int4   atom         = cAmoebaSim.pAmoebaBondID[pos];
            float4 atomA        = cSim.pPosq[atom.x];
            float4 atomB        = cSim.pPosq[atom.y];
            float4 bond         = cAmoebaSim.pAmoebaBondParameter[pos];
            float dx            = atomB.x - atomA.x;
            float dy            = atomB.y - atomA.y;
            float dz            = atomB.z - atomA.z;
            float r2            = dx * dx + dy * dy + dz * dz;
            float r             = sqrt(r2);
            float deltaIdeal    = r - bond.x;
#if defined OLD
            energy             += 0.5f * bond.y * deltaIdeal * deltaIdeal;
            float dEdR          = bond.y * deltaIdeal;
            dEdR                = (r > 0.0f) ? (dEdR / r) : 0.0f;
#else
            float deltaIdeal2   = deltaIdeal*deltaIdeal;
            energy             += bond.y * deltaIdeal2*( 1.0f + bond.z*deltaIdeal + bond.w*deltaIdeal2 );
            float dEdR          = 2.0f*bond.y * deltaIdeal*( 1.0f + 1.5f*bond.z*deltaIdeal + 2.0f*bond.w*deltaIdeal2 );
            dEdR                = (r > 0.0f) ? (dEdR / r) : 0.0f;
#endif

            dx                 *= dEdR;
            dy                 *= dEdR;
            dz                 *= dEdR;
            unsigned int offsetA                = atom.x + atom.z * cSim.stride;
            unsigned int offsetB                = atom.y + atom.w * cSim.stride;
            float4 forceA                       = cSim.pForce4[offsetA];
            float4 forceB                       = cSim.pForce4[offsetB];
            forceA.x                           += dx;
            forceA.y                           += dy;
            forceA.z                           += dz;
            forceB.x                           -= dx;
            forceB.y                           -= dy;
            forceB.z                           -= dz;
            cSim.pForce4[offsetA]               = forceA;
            cSim.pForce4[offsetB]               = forceB;    
        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaAngle_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaBond_offset;
        if (pos1 < cAmoebaSim.amoebaAngles )
        {
            int4   atom1            = cAmoebaSim.pAmoebaAngleID1[pos1];  

            /*
               bond_angle1.x        ideal
               bond_angle1.y        k
            */

            float2 bond_angle1      = cAmoebaSim.pAmoebaAngleParameter[pos1];

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];
            float4 v0;
            float4 v1;
            v0.x                    = a2.x - a1.x;
            v0.y                    = a2.y - a1.y;
            v0.z                    = a2.z - a1.z;
            v1.x                    = a2.x - a3.x;
            v1.y                    = a2.y - a3.y;
            v1.z                    = a2.z - a3.z;
            float3 cp;
            CROSS_PRODUCT(v0, v1, cp);
            float rp                = DOT3(cp, cp); //cx * cx + cy * cy + cz * cz;
            rp                      = max(sqrt(rp), 1.0e-06f);
            float r21               = DOT3(v0, v0); // dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
            float r23               = DOT3(v1, v1); // dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
            float dot               = DOT3(v0, v1); // dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
            float cosine            = dot / sqrt(r21 * r23);

// e     = angunit * force * dt2 * (1.0d0+cang*dt+qang*dt2+pang*dt3+sang*dt4)
// deddt = angunit * force * dt * radian * (2.0d0 + 3.0d0*cang*dt + 4.0d0*qang*dt2 + 5.0d0*pang*dt3 + 6.0d0*sang*dt4)

            float angle             = acos(cosine);
            float deltaIdeal        = angle*(180.0f/LOCAL_HACK_PI) - bond_angle1.x;
            float deltaIdeal2       = deltaIdeal *deltaIdeal;
            float deltaIdeal3       = deltaIdeal *deltaIdeal2;
            float deltaIdeal4       = deltaIdeal2*deltaIdeal2;
            energy                 += bond_angle1.y*deltaIdeal2*( 1.0f + cAmoebaSim.amoebaAngleCubicK*deltaIdeal               +
                                                                         cAmoebaSim.amoebaAngleQuarticK*deltaIdeal2            +
                                                                         cAmoebaSim.amoebaAnglePenticK*deltaIdeal3             +
                                                                         cAmoebaSim.amoebaAngleSexticK*deltaIdeal4 );

            float dEdR              = bond_angle1.y*deltaIdeal*( 2.0f + 3.0f*cAmoebaSim.amoebaAngleCubicK*deltaIdeal           +
                                                                        4.0f*cAmoebaSim.amoebaAngleQuarticK*deltaIdeal2        +
                                                                        5.0f*cAmoebaSim.amoebaAnglePenticK*deltaIdeal3         + 
                                                                        6.0f*cAmoebaSim.amoebaAngleSexticK*deltaIdeal4 );

            dEdR                   *= LOCAL_HACK_RADIAN;

            float termA             =  dEdR / (r21 * rp);
            float termC             = -dEdR / (r23 * rp);
            float3 c21;
            float3 c23;
            CROSS_PRODUCT(v0, cp, c21);
            CROSS_PRODUCT(v1, cp, c23);
            c21.x                  *= termA;
            c21.y                  *= termA;
            c21.z                  *= termA;
            c23.x                  *= termC;
            c23.y                  *= termC;
            c23.z                  *= termC;
            int2 atom2              = cAmoebaSim.pAmoebaAngleID2[pos1];
            unsigned int offset     = atom1.x + atom1.w * cSim.stride;
            float4 force            = cSim.pForce4[offset]; 
            force.x                += c21.x;
            force.y                += c21.y;
            force.z                += c21.z;
            cSim.pForce4[offset]    = force;
            offset                  = atom1.y + atom2.x * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= (c21.x + c23.x);
            force.y                -= (c21.y + c23.y);
            force.z                -= (c21.z + c23.z);
            cSim.pForce4[offset]    = force;
            offset                  = atom1.z + atom2.y * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                += c23.x;
            force.y                += c23.y;
            force.z                += c23.z;
            cSim.pForce4[offset]    = force;
        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaInPlaneAngle_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaAngle_offset;
        if (pos1 < cAmoebaSim.amoebaInPlaneAngles )
        {
            int4   atom1            = cAmoebaSim.pAmoebaInPlaneAngleID1[pos1];  

            /*
               bond_angle1.x        ideal
               bond_angle1.y        k
            */

            float2 bond_angle1      = cAmoebaSim.pAmoebaInPlaneAngleParameter[pos1];

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];
            float4 a4               = cSim.pPosq[atom1.w];

            float xad               = a1.x - a4.x;
            float yad               = a1.y - a4.y;
            float zad               = a1.z - a4.z;

            float xbd               = a2.x - a4.x;
            float ybd               = a2.y - a4.y;
            float zbd               = a2.z - a4.z;

            float xcd               = a3.x - a4.x;
            float ycd               = a3.y - a4.y;
            float zcd               = a3.z - a4.z;

            float xt                = yad*zcd - zad*ycd;
            float yt                = zad*xcd - xad*zcd;
            float zt                = xad*ycd - yad*xcd;

            float rt2               = xt*xt + yt*yt + zt*zt;

            float delta             = -(xt*xbd + yt*ybd + zt*zbd) / rt2;

            float xip               = a2.x + xt*delta;
            float yip               = a2.y + yt*delta;
            float zip               = a2.z + zt*delta;

            float xap               = a1.x - xip;
            float yap               = a1.y - yip;
            float zap               = a1.z - zip;

            float xcp               = a3.x - xip;
            float ycp               = a3.y - yip;
            float zcp               = a3.z - zip;

            float rap2              = xap*xap + yap*yap + zap*zap;
            float rcp2              = xcp*xcp + ycp*ycp + zcp*zcp;

            float xm                = ycp*zap - zcp*yap;
            float ym                = zcp*xap - xcp*zap;
            float zm                = xcp*yap - ycp*xap;

            float rm                = sqrtf(xm*xm + ym*ym + zm*zm);
                  rm                = rm > 0.000001f ? rm : 0.000001f;
            float dot               = xap*xcp + yap*ycp + zap*zcp;
            float product           = sqrtf(rap2*rcp2);
            float cosine            = product > 0.0f ? (dot/product) : 0.0f;
                  cosine            = cosine >  1.0f ?  1.0f : cosine;
                  cosine            = cosine < -1.0f ? -1.0f : cosine;
            float angle             = acos(cosine);

            // if product == 0, set force/energy to 0

            float deltaIdeal        = product > 0.0f ? (angle*(180.0f/LOCAL_HACK_PI) - bond_angle1.x) : 0.0f;
            float deltaIdeal2       = deltaIdeal *deltaIdeal;
            float deltaIdeal3       = deltaIdeal *deltaIdeal2;
            float deltaIdeal4       = deltaIdeal2*deltaIdeal2;

            energy                 += bond_angle1.y*deltaIdeal2*( 1.0f +  cAmoebaSim.amoebaInPlaneAngleCubicK*deltaIdeal          +
                                                                          cAmoebaSim.amoebaInPlaneAngleQuarticK*deltaIdeal2       +
                                                                          cAmoebaSim.amoebaInPlaneAnglePenticK*deltaIdeal3        +
                                                                          cAmoebaSim.amoebaInPlaneAngleSexticK*deltaIdeal4 );


            float dEdR              = bond_angle1.y*deltaIdeal*( 2.0f + 3.0f*cAmoebaSim.amoebaInPlaneAngleCubicK*deltaIdeal       + 
                                                                        4.0f*cAmoebaSim.amoebaInPlaneAngleQuarticK*deltaIdeal2    +
                                                                        5.0f*cAmoebaSim.amoebaInPlaneAnglePenticK*deltaIdeal3     + 
                                                                        6.0f*cAmoebaSim.amoebaInPlaneAngleSexticK*deltaIdeal4 );

            dEdR                   *= LOCAL_HACK_RADIAN;

            float terma             = -dEdR/ (rap2*rm);
            float termc             = dEdR / (rcp2*rm);

            float dedxia            = terma * (yap*zm-zap*ym);
            float dedyia            = terma * (zap*xm-xap*zm);
            float dedzia            = terma * (xap*ym-yap*xm);

            float dedxic            = termc * (ycp*zm-zcp*ym);
            float dedyic            = termc * (zcp*xm-xcp*zm);
            float dedzic            = termc * (xcp*ym-ycp*xm);

            float dedxip            = -dedxia - dedxic;
            float dedyip            = -dedyia - dedyic;
            float dedzip            = -dedzia - dedzic;

            float delta2            = 2.0f * delta;
            float ptrt2             = (dedxip*xt + dedyip*yt + dedzip*zt) / rt2;

            float term              = (zcd*ybd-ycd*zbd) + delta2*(yt*zcd-zt*ycd);
            float dpdxia            = delta*(ycd*dedzip-zcd*dedyip) + term*ptrt2;

                  term              = (xcd*zbd-zcd*xbd) + delta2*(zt*xcd-xt*zcd);
            float dpdyia            = delta*(zcd*dedxip-xcd*dedzip) + term*ptrt2;

                  term              = (ycd*xbd-xcd*ybd) + delta2*(xt*ycd-yt*xcd);
            float dpdzia            = delta*(xcd*dedyip-ycd*dedxip) + term*ptrt2;

                  term              = (yad*zbd-zad*ybd) + delta2*(zt*yad-yt*zad);
            float dpdxic            = delta*(zad*dedyip-yad*dedzip) + term*ptrt2;

                  term              = (zad*xbd-xad*zbd) + delta2*(xt*zad-zt*xad);
            float dpdyic            = delta*(xad*dedzip-zad*dedxip) + term*ptrt2;

                  term              = (xad*ybd-yad*xbd) + delta2*(yt*xad-xt*yad);
            float dpdzic            = delta*(yad*dedxip-xad*dedyip) + term*ptrt2;

                  dedxia            = dedxia + dpdxia;
                  dedyia            = dedyia + dpdyia;
                  dedzia            = dedzia + dpdzia;

            float dedxib            = dedxip;
            float dedyib            = dedyip;
            float dedzib            = dedzip;

                  dedxic            = dedxic + dpdxic;
                  dedyic            = dedyic + dpdyic;
                  dedzic            = dedzic + dpdzic;

            float dedxid            = -dedxia - dedxib - dedxic;
            float dedyid            = -dedyia - dedyib - dedyic;
            float dedzid            = -dedzia - dedzib - dedzic;

            int4 atom2              = cAmoebaSim.pAmoebaInPlaneAngleID2[pos1];
            unsigned int offset     = atom1.x + atom2.x * cSim.stride;
            float4 force            = cSim.pForce4[offset]; 
            force.x                -= dedxia;
            force.y                -= dedyia;
            force.z                -= dedzia;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.y + atom2.y * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxib;
            force.y                -= dedyib;
            force.z                -= dedzib;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.z + atom2.z * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxic;
            force.y                -= dedyic;
            force.z                -= dedzic;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.w + atom2.w * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxid;
            force.y                -= dedyid;
            force.z                -= dedzid;
            cSim.pForce4[offset]    = force;
        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaTorsion_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaInPlaneAngle_offset;
        if (pos1 < cAmoebaSim.amoebaTorsions )
        {
            int4   atom1            = cAmoebaSim.pAmoebaTorsionID1[pos1];  

            /*
               torsionParam1.x      amplitude(1)
               torsionParam1.y      phase(1)
               torsionParam1.z      amplitude(2)
               torsionParam1.w      phase(2)
               torsionParam2.x      amplitude(3)
               torsionParam2.y      phase(3)
            */

            float4 torsionParam1    = cAmoebaSim.pAmoebaTorsionParameter1[pos1];
            float2 torsionParam2    = cAmoebaSim.pAmoebaTorsionParameter2[pos1];

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];
            float4 a4               = cSim.pPosq[atom1.w];

            float xba               = a2.x - a1.x;
            float yba               = a2.y - a1.y;
            float zba               = a2.z - a1.z;

            float xcb               = a3.x - a2.x;
            float ycb               = a3.y - a2.y;
            float zcb               = a3.z - a2.z;

            float xdc               = a4.x - a3.x;
            float ydc               = a4.y - a3.y;
            float zdc               = a4.z - a3.z;

            float xt                = yba*zcb - ycb*zba;
            float yt                = zba*xcb - zcb*xba;
            float zt                = xba*ycb - xcb*yba;
            float xu                = ycb*zdc - ydc*zcb;
            float yu                = zcb*xdc - zdc*xcb;
            float zu                = xcb*ydc - xdc*ycb;

            float xtu               = yt*zu - yu*zt;
            float ytu               = zt*xu - zu*xt;
            float ztu               = xt*yu - xu*yt;

            float rt2               = xt*xt + yt*yt + zt*zt;
            float ru2               = xu*xu + yu*yu + zu*zu;

            float rtru              = sqrtf(rt2 * ru2);

            float rcb               = sqrtf(xcb*xcb + ycb*ycb + zcb*zcb);
            float cosine            = rtru > 0.0f ? ( (xt*xu + yt*yu + zt*zu) / rtru) : 0.0f;
            float sine              = rtru > 0.0f ? ( (xcb*xtu + ycb*ytu + zcb*ztu) / (rcb*rtru) ) : 0.0f;

            // set the torsional parameters for this angle;

            float    v1             = torsionParam1.x;
            float    angle          = torsionParam1.y;
            float    c1             = cos( angle );
            float    s1             = sin( angle );

            float    v2             = torsionParam1.z;
                     angle          = torsionParam1.w;
            float    c2             = cos( angle );
            float    s2             = sin( angle );


            float    v3             = torsionParam2.x;
                     angle          = torsionParam2.y;
            float    c3             = cos( angle );
            float    s3             = sin( angle );

            // compute the multiple angle trigonometry and the phase terms

            float    cosine2          = cosine*cosine - sine*sine;
            float    sine2            = 2.0f * cosine * sine;
            float    cosine3          = cosine*cosine2 - sine*sine2;
            float    sine3            = cosine*sine2 + sine*cosine2;

            // not deleted since may be needed in future
/*
            float    cosine4          = cosine*cosine3 - sine*sine3;
            float    sine4            = cosine*sine3 + sine*cosine3;
            float    cosine5          = cosine*cosine4 - sine*sine4;
            float    sine5            = cosine*sine4 + sine*cosine4;
            float    cosine6          = cosine*cosine5 - sine*sine5;
            float    sine6            = cosine*sine5 + sine*cosine5;
*/
            float    phi1             = 1.0f + (cosine*c1 + sine*s1);
            float    phi2             = 1.0f + (cosine2*c2 + sine2*s2);
            float    phi3             = 1.0f + (cosine3*c3 + sine3*s3);
/*
            float    phi4             = 1.0f + (cosine4*c4 + sine4*s4);
            float    phi5             = 1.0f + (cosine5*c5 + sine5*s5);
            float    phi6             = 1.0f + (cosine6*c6 + sine6*s6);
*/
            float    dphi1            = (cosine*s1 - sine*c1);
            float    dphi2            = 2.0f * (cosine2*s2 - sine2*c2);
            float    dphi3            = 3.0f * (cosine3*s3 - sine3*c3);
/*
            float    dphi4            = 4.0f * (cosine4*s4 - sine4*c4);
            float    dphi5            = 5.0f * (cosine5*s5 - sine5*c5);
            float    dphi6            = 6.0f * (cosine6*s6 - sine6*c6);
*/

            // calculate torsional energy and master chain rule term

            energy                   += v1*phi1 + v2*phi2 + v3*phi3;
     //                                 + v4*phi4 + v5*phi5 + v6*phi6;

            float    dedphi           = v1*dphi1 + v2*dphi2 + v3*dphi3;
     //                                 + v4*dphi4 + v5*dphi5 + v6*dphi6;

            // chain rule terms for first derivative components

            float    xca              = a3.x - a1.x;
            float    yca              = a3.y - a1.y;
            float    zca              = a3.z - a1.z;

            float    xdb              = a4.x - a2.x;
            float    ydb              = a4.y - a2.y;
            float    zdb              = a4.z - a2.z;

            float    dedxt            =  dedphi * (yt*zcb - ycb*zt) / (rt2*rcb);
            float    dedyt            =  dedphi * (zt*xcb - zcb*xt) / (rt2*rcb);
            float    dedzt            =  dedphi * (xt*ycb - xcb*yt) / (rt2*rcb);
            float    dedxu            = -dedphi * (yu*zcb - ycb*zu) / (ru2*rcb);
            float    dedyu            = -dedphi * (zu*xcb - zcb*xu) / (ru2*rcb);
            float    dedzu            = -dedphi * (xu*ycb - xcb*yu) / (ru2*rcb);

            // compute first derivative components for this angle

            float    dedxia           = zcb*dedyt - ycb*dedzt;
            float    dedyia           = xcb*dedzt - zcb*dedxt;
            float    dedzia           = ycb*dedxt - xcb*dedyt;
            float    dedxib           = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
            float    dedyib           = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
            float    dedzib           = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;
            float    dedxic           = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
            float    dedyic           = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
            float    dedzic           = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;
            float    dedxid           = zcb*dedyu - ycb*dedzu;
            float    dedyid           = xcb*dedzu - zcb*dedxu;
            float    dedzid           = ycb*dedxu - xcb*dedyu;

            int4 atom2                = cAmoebaSim.pAmoebaTorsionID2[pos1];
            unsigned int offset       = atom1.x + atom2.x * cSim.stride;
            float4 force              = cSim.pForce4[offset]; 
            force.x                  -= dedxia;
            force.y                  -= dedyia;
            force.z                  -= dedzia;
            cSim.pForce4[offset]      = force;

            offset                    = atom1.y + atom2.y * cSim.stride;
            force                     = cSim.pForce4[offset];
            force.x                  -= dedxib;
            force.y                  -= dedyib;
            force.z                  -= dedzib;
            cSim.pForce4[offset]      = force;

            offset                    = atom1.z + atom2.z * cSim.stride;
            force                     = cSim.pForce4[offset];
            force.x                  -= dedxic;
            force.y                  -= dedyic;
            force.z                  -= dedzic;
            cSim.pForce4[offset]      = force;

            offset                    = atom1.w + atom2.w * cSim.stride;
            force                     = cSim.pForce4[offset];
            force.x                  -= dedxid;
            force.y                  -= dedyid;
            force.z                  -= dedzid;
            cSim.pForce4[offset]      = force;
        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaPiTorsion_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaTorsion_offset;
        if (pos1 < cAmoebaSim.amoebaPiTorsions )
        {
            int4   atom1            = cAmoebaSim.pAmoebaPiTorsionID1[pos1];  
            int4   atom2            = cAmoebaSim.pAmoebaPiTorsionID2[pos1];  

            /*
               torsionParam1.x      amplitude(1)
               torsionParam1.y      phase(1)
               torsionParam1.z      amplitude(2)
               torsionParam1.w      phase(2)
               torsionParam2.x      amplitude(3)
               torsionParam2.y      phase(3)
            */

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];
            float4 a4               = cSim.pPosq[atom1.w];
            float4 a5               = cSim.pPosq[atom2.x];
            float4 a6               = cSim.pPosq[atom2.y];

            // compute the value of the pi-orbital torsion angle

            float xad               = a1.x - a4.x;
            float yad               = a1.y - a4.y;
            float zad               = a1.z - a4.z;

            float xbd               = a2.x - a4.x;
            float ybd               = a2.y - a4.y;
            float zbd               = a2.z - a4.z;

            float xec               = a5.x - a3.x;
            float yec               = a5.y - a3.y;
            float zec               = a5.z - a3.z;

            float xgc               = a6.x - a3.x;
            float ygc               = a6.y - a3.y;
            float zgc               = a6.z - a3.z;

            float xip               = yad*zbd - ybd*zad + a3.x;
            float yip               = zad*xbd - zbd*xad + a3.y;
            float zip               = xad*ybd - xbd*yad + a3.z;

            float xiq               = yec*zgc - ygc*zec + a4.x;
            float yiq               = zec*xgc - zgc*xec + a4.y;
            float ziq               = xec*ygc - xgc*yec + a4.z;

            float xcp               = a3.x - xip;
            float ycp               = a3.y - yip;
            float zcp               = a3.z - zip;

            float xdc               = a4.x - a3.x;
            float ydc               = a4.y - a3.y;
            float zdc               = a4.z - a3.z;

            float xqd               = xiq - a4.x;
            float yqd               = yiq - a4.y;
            float zqd               = ziq - a4.z;

            float xt                = ycp*zdc - ydc*zcp;
            float yt                = zcp*xdc - zdc*xcp;
            float zt                = xcp*ydc - xdc*ycp;

            float xu                = ydc*zqd - yqd*zdc;
            float yu                = zdc*xqd - zqd*xdc;
            float zu                = xdc*yqd - xqd*ydc;

            float xtu               = yt*zu - yu*zt;
            float ytu               = zt*xu - zu*xt;
            float ztu               = xt*yu - xu*yt;

            float rt2               = xt*xt + yt*yt + zt*zt;
            float ru2               = xu*xu + yu*yu + zu*zu;

            float rtru              = sqrtf(rt2 * ru2);
            float rdc               = sqrt(xdc*xdc + ydc*ydc + zdc*zdc);

            float cosine            = rtru > 0.0f ? (xt*xu + yt*yu + zt*zu) / rtru : 0.0f;
            float sine              = (rtru*rdc) > 0.0f ? (xdc*xtu + ydc*ytu + zdc*ztu) / (rdc*rtru) : 0.0f;

            // zero energy/force if rtru == 0

            float v2                = cAmoebaSim.pAmoebaPiTorsionParameter[pos1];
                  v2                = rtru > 0.0f ? v2 : 0.0f;

            // compute the multiple angle trigonometry and the phase terms

            float cosine2           = cosine*cosine - sine*sine;
            float sine2             = 2.0f * cosine * sine;
            float phi2              = 1.0f - cosine2;
            float dphi2             = 2.0f * sine2;

            // calculate pi-orbital torsion energy and master chain rule term

            energy                 += v2 * phi2;
            float dedphi            = v2 * dphi2;

            // chain rule terms for first derivative components

            float xdp               = a4.x - xip;
            float ydp               = a4.y - yip;
            float zdp               = a4.z - zip;

            float xqc               = xiq - a3.x;
            float yqc               = yiq - a3.y;
            float zqc               = ziq - a3.z;

            float dedxt             =  dedphi * (yt*zdc - ydc*zt) / (rt2*rdc);
            float dedyt             =  dedphi * (zt*xdc - zdc*xt) / (rt2*rdc);
            float dedzt             =  dedphi * (xt*ydc - xdc*yt) / (rt2*rdc);

            float dedxu             = -dedphi * (yu*zdc - ydc*zu) / (ru2*rdc);
            float dedyu             = -dedphi * (zu*xdc - zdc*xu) / (ru2*rdc);
            float dedzu             = -dedphi * (xu*ydc - xdc*yu) / (ru2*rdc);

            // compute first derivative components for pi-orbital angle

            float dedxip            = zdc*dedyt - ydc*dedzt;
            float dedyip            = xdc*dedzt - zdc*dedxt;
            float dedzip            = ydc*dedxt - xdc*dedyt;

            float dedxic            = ydp*dedzt - zdp*dedyt + zqd*dedyu - yqd*dedzu;
            float dedyic            = zdp*dedxt - xdp*dedzt + xqd*dedzu - zqd*dedxu;
            float dedzic            = xdp*dedyt - ydp*dedxt + yqd*dedxu - xqd*dedyu;

            float dedxid            = zcp*dedyt - ycp*dedzt + yqc*dedzu - zqc*dedyu;
            float dedyid            = xcp*dedzt - zcp*dedxt + zqc*dedxu - xqc*dedzu;
            float dedzid            = ycp*dedxt - xcp*dedyt + xqc*dedyu - yqc*dedxu;

            float dedxiq            = zdc*dedyu - ydc*dedzu;
            float dedyiq            = xdc*dedzu - zdc*dedxu;
            float dedziq            = ydc*dedxu - xdc*dedyu;

            // compute first derivative components for individual atoms

            float dedxia            = ybd*dedzip - zbd*dedyip;
            float dedyia            = zbd*dedxip - xbd*dedzip;
            float dedzia            = xbd*dedyip - ybd*dedxip;

            float dedxib            = zad*dedyip - yad*dedzip;
            float dedyib            = xad*dedzip - zad*dedxip;
            float dedzib            = yad*dedxip - xad*dedyip;

            float dedxie            = ygc*dedziq - zgc*dedyiq;
            float dedyie            = zgc*dedxiq - xgc*dedziq;
            float dedzie            = xgc*dedyiq - ygc*dedxiq;

            float dedxig            = zec*dedyiq - yec*dedziq;
            float dedyig            = xec*dedziq - zec*dedxiq;
            float dedzig            = yec*dedxiq - xec*dedyiq;

                  dedxic            = dedxic + dedxip - dedxie - dedxig;
                  dedyic            = dedyic + dedyip - dedyie - dedyig;
                  dedzic            = dedzic + dedzip - dedzie - dedzig;
                  dedxid            = dedxid + dedxiq - dedxia - dedxib;
                  dedyid            = dedyid + dedyiq - dedyia - dedyib;
                  dedzid            = dedzid + dedziq - dedzia - dedzib;

            int4 atom3              = cAmoebaSim.pAmoebaPiTorsionID3[pos1];
            unsigned int offset     = atom1.x + atom2.z * cSim.stride;
            float4 force            = cSim.pForce4[offset]; 
            force.x                -= dedxia;
            force.y                -= dedyia;
            force.z                -= dedzia;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.y + atom2.w * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxib;
            force.y                -= dedyib;
            force.z                -= dedzib;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.z + atom3.x * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxic;
            force.y                -= dedyic;
            force.z                -= dedzic;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.w + atom3.y * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxid;
            force.y                -= dedyid;
            force.z                -= dedzid;
            cSim.pForce4[offset]    = force;

            offset                  = atom2.x + atom3.z * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxie;
            force.y                -= dedyie;
            force.z                -= dedzie;
            cSim.pForce4[offset]    = force;

            offset                  = atom2.y + atom3.w * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxig;
            force.y                -= dedyig;
            force.z                -= dedzig;
            cSim.pForce4[offset]    = force;

        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaStretchBend_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaPiTorsion_offset;
        if (pos1 < cAmoebaSim.amoebaStretchBends )
        {
            int4   atom1            = cAmoebaSim.pAmoebaStretchBendID1[pos1];  

            /*
               parameters.x        length AB
               parameters.y        length CB
               parameters.z        angle (in radians)
               parameters.w        k
            */

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];

            // compute the value of the bond angle

            float xab               = a1.x - a2.x;
            float yab               = a1.y - a2.y;
            float zab               = a1.z - a2.z;

            float xcb               = a3.x - a2.x;
            float ycb               = a3.y - a2.y;
            float zcb               = a3.z - a2.z;

            float rab               = sqrtf(xab*xab + yab*yab + zab*zab);
            float rcb               = sqrtf(xcb*xcb + ycb*ycb + zcb*zcb);

            float xp                = ycb*zab - zcb*yab;
            float yp                = zcb*xab - xcb*zab;
            float zp                = xcb*yab - ycb*xab;

            float rp                = sqrt(xp*xp + yp*yp + zp*zp);

            float    dot            = xab*xcb + yab*ycb + zab*zcb;
            float    cosine         = rab*rcb > 0.0f ? (dot / (rab*rcb)) : 1.0f;
                     cosine         = cosine >  1.0f ?  1.0f : cosine;
                     cosine         = cosine < -1.0f ? -1.0f : cosine;
            float    angle          = acos(cosine);

            // find chain rule terms for the bond angle deviation

            float4 parameters       = cAmoebaSim.pAmoebaStretchBendParameter[pos1];

            float    dt             = LOCAL_HACK_RADIAN*(angle - parameters.z);
            float    terma          = rab*rp != 0.0f ? (-LOCAL_HACK_RADIAN/(rab*rab*rp)) : 0.0f;
            float    termc          = rcb*rp != 0.0f ? ( LOCAL_HACK_RADIAN/(rcb*rcb*rp)) : 0.0f;

            float    ddtdxia        = terma * (yab*zp-zab*yp);
            float    ddtdyia        = terma * (zab*xp-xab*zp);
            float    ddtdzia        = terma * (xab*yp-yab*xp);

            float    ddtdxic        = termc * (ycb*zp-zcb*yp);
            float    ddtdyic        = termc * (zcb*xp-xcb*zp);
            float    ddtdzic        = termc * (xcb*yp-ycb*xp);

            // find chain rule terms for the bond length deviations

            float    dr             = parameters.x > 0.0f ? (rab - parameters.x) : 0.0f;
                     terma          = parameters.x > 0.0f ? (1.0f/rab) : 0.0f;

                     dr            += parameters.y > 0.0f ? (rcb - parameters.y) : 0.0f;
                     termc          = parameters.y > 0.0f ? (1.0f/rcb) : 0.0f;

            float    ddrdxia        = terma * xab;
            float    ddrdyia        = terma * yab;
            float    ddrdzia        = terma * zab;

            float    ddrdxic        = termc * xcb;
            float    ddrdyic        = termc * ycb;
            float    ddrdzic        = termc * zcb;

            // get the energy and master chain rule terms for derivatives

            float    term           = (rp != 0.0f) ? parameters.w : 0.0f;
                     energy        += term*dt*dr;

            float    dedxia         = term * (dt*ddrdxia+ddtdxia*dr);
            float    dedyia         = term * (dt*ddrdyia+ddtdyia*dr);
            float    dedzia         = term * (dt*ddrdzia+ddtdzia*dr);

            float    dedxic         = term * (dt*ddrdxic+ddtdxic*dr);
            float    dedyic         = term * (dt*ddrdyic+ddtdyic*dr);
            float    dedzic         = term * (dt*ddrdzic+ddtdzic*dr);

            float    dedxib         = -dedxia - dedxic;
            float    dedyib         = -dedyia - dedyic;
            float    dedzib         = -dedzia - dedzic;

            // increment the total stretch-bend energy and derivatives

            unsigned int offset     = atom1.x + atom1.w * cSim.stride;
            float4 force            = cSim.pForce4[offset]; 
            force.x                -= dedxia;
            force.y                -= dedyia;
            force.z                -= dedzia;
            cSim.pForce4[offset]    = force;

            int2 atom2              = cAmoebaSim.pAmoebaStretchBendID2[pos1];
            offset                  = atom1.y + atom2.x * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxib;
            force.y                -= dedyib;
            force.z                -= dedzib;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.z + atom2.y * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxic;
            force.y                -= dedyic;
            force.z                -= dedzic;
            cSim.pForce4[offset]    = force;

        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaOutOfPlaneBend_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaStretchBend_offset;
        if (pos1 < cAmoebaSim.amoebaOutOfPlaneBends )
        {
            int4   atom1            = cAmoebaSim.pAmoebaOutOfPlaneBendID1[pos1];  

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];
            float4 a4               = cSim.pPosq[atom1.w];

            // compute the value of the bond angle

            float xab               = a1.x - a2.x;
            float yab               = a1.y - a2.y;
            float zab               = a1.z - a2.z;

            float xcb               = a3.x - a2.x;
            float ycb               = a3.y - a2.y;
            float zcb               = a3.z - a2.z;

            // compute the out-of-plane bending angle

            float xdb               = a4.x - a2.x;
            float ydb               = a4.y - a2.y;
            float zdb               = a4.z - a2.z;

            float xad               = a1.x - a4.x;
            float yad               = a1.y - a4.y;
            float zad               = a1.z - a4.z;

            float xcd               = a3.x - a4.x;
            float ycd               = a3.y - a4.y;
            float zcd               = a3.z - a4.z;

            float rdb2              = xdb*xdb + ydb*ydb + zdb*zdb;
            float rad2              = xad*xad + yad*yad + zad*zad;
            float rcd2              = xcd*xcd + ycd*ycd + zcd*zcd;

            float ee                = xab*(ycb*zdb-zcb*ydb) +  yab*(zcb*xdb-xcb*zdb) + zab*(xcb*ydb-ycb*xdb);

            float dot               = xad*xcd + yad*ycd + zad*zcd;
            float cc                = rad2*rcd2 - dot*dot;
            float bkk2              = (cc != 0.0f) ? (rdb2 - ee*ee/cc) : 0.0f;

            float cosine            = (rdb2  !=  0.0f) ? sqrtf(bkk2/rdb2) : 0.0f;
                  cosine            = (cosine >  1.0f) ?  1.0f            : cosine;
                  cosine            = (cosine < -1.0f) ? -1.0f            : cosine;

/*
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
            if (opbtyp .eq. 'W-D-C') then
               rab2 = xab*xab + yab*yab + zab*zab
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               dot = xab*xcb+yab*ycb+zab*zcb
               cc = rab2*rcb2 - dot*dot
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
            else if (opbtyp .eq. 'ALLINGER') then
               rad2 = xad*xad + yad*yad + zad*zad
               rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
               dot = xad*xcd + yad*ycd + zad*zcd
               cc = rad2*rcd2 - dot*dot
            end if
c
c     find the out-of-plane angle bending energy
c
            ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &              + zdb*(xab*ycb-yab*xcb)
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            if (rdb2.ne.0.0d0 .and. cc.ne.0.0d0) then
               bkk2 = rdb2 - ee*ee/cc
               cosine = sqrt(bkk2/rdb2)
               cosine = min(1.0d0,max(-1.0d0,cosine))
               angle = radian * acos(cosine)

*/
            // find the out-of-plane energy and master chain rule terms

            float    dt             = LOCAL_HACK_RADIAN_D*acos(cosine);
            float    dt2            = dt  * dt;
            float    dt3            = dt2 * dt;
            float    dt4            = dt2 * dt2;
            float      k            = (rdb2 != 0.0f && cc != 0.0f) ? cAmoebaSim.pAmoebaOutOfPlaneBendParameter[pos1] : 0.0f;
                    
                   energy          += k*dt2*(1.0f + (cAmoebaSim.amoebaOutOfPlaneBendCubicK*  dt ) +
                                                    (cAmoebaSim.amoebaOutOfPlaneBendQuarticK*dt2) +
                                                    (cAmoebaSim.amoebaOutOfPlaneBendPenticK* dt3) +
                                                    (cAmoebaSim.amoebaOutOfPlaneBendSexticK* dt4) );

            float    deddt          = k*dt*LOCAL_HACK_RADIAN*(2.0f                                           + 
                                                             (3.0f*cAmoebaSim.amoebaOutOfPlaneBendCubicK*  dt ) +
                                                             (4.0f*cAmoebaSim.amoebaOutOfPlaneBendQuarticK*dt2) +
                                                             (5.0f*cAmoebaSim.amoebaOutOfPlaneBendPenticK* dt3) +
                                                             (6.0f*cAmoebaSim.amoebaOutOfPlaneBendSexticK* dt4) );

            float    eeSign         = ee >= 0.0f ? 1.0f : -1.0f;
            float    dedcos         = -deddt*eeSign/sqrtf(cc*bkk2);

            // chain rule terms for first derivative components

            float    term            = ee / cc;

            float    dccdxia         = (xad*rcd2-xcd*dot) * term;
            float    dccdyia         = (yad*rcd2-ycd*dot) * term;
            float    dccdzia         = (zad*rcd2-zcd*dot) * term;

            float    dccdxic         = (xcd*rad2-xad*dot) * term;
            float    dccdyic         = (ycd*rad2-yad*dot) * term;
            float    dccdzic         = (zcd*rad2-zad*dot) * term;

            float    dccdxid         = -dccdxia - dccdxic;
            float    dccdyid         = -dccdyia - dccdyic;
            float    dccdzid         = -dccdzia - dccdzic;

                     term            = ee / rdb2;

            float    deedxia         = ydb*zcb - zdb*ycb;
            float    deedyia         = zdb*xcb - xdb*zcb;
            float    deedzia         = xdb*ycb - ydb*xcb;

            float    deedxic         = yab*zdb - zab*ydb;
            float    deedyic         = zab*xdb - xab*zdb;
            float    deedzic         = xab*ydb - yab*xdb;

            float    deedxid         = ycb*zab - zcb*yab + xdb*term;
            float    deedyid         = zcb*xab - xcb*zab + ydb*term;
            float    deedzid         = xcb*yab - ycb*xab + zdb*term;

            // compute first derivative components for this angle

            float    dedxia          = dedcos * (dccdxia+deedxia);
            float    dedyia          = dedcos * (dccdyia+deedyia);
            float    dedzia          = dedcos * (dccdzia+deedzia);

            float    dedxic          = dedcos * (dccdxic+deedxic);
            float    dedyic          = dedcos * (dccdyic+deedyic);
            float    dedzic          = dedcos * (dccdzic+deedzic);

            float    dedxid          = dedcos * (dccdxid+deedxid);
            float    dedyid          = dedcos * (dccdyid+deedyid);
            float    dedzid          = dedcos * (dccdzid+deedzid);

            float    dedxib          = -dedxia - dedxic - dedxid;
            float    dedyib          = -dedyia - dedyic - dedyid;
            float    dedzib          = -dedzia - dedzic - dedzid;

            // increment the out-of-plane bending energy and gradient;

            int4   atom2            = cAmoebaSim.pAmoebaOutOfPlaneBendID2[pos1];  

            unsigned int offset     = atom1.x + atom2.x * cSim.stride;
            float4 force            = cSim.pForce4[offset]; 
            force.x                -= dedxia;
            force.y                -= dedyia;
            force.z                -= dedzia;

force.x = bkk2;
force.y = rdb2;
force.z = cosine;
force.w = dt;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.y + atom2.y * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxib;
            force.y                -= dedyib;
            force.z                -= dedzib;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.z + atom2.z * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxic;
            force.y                -= dedyic;
            force.z                -= dedzic;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.w + atom2.w * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxid;
            force.y                -= dedyid;
            force.z                -= dedzid;
            cSim.pForce4[offset]    = force;

        }
        pos += blockDim.x * gridDim.x;
    }

    while (pos < cAmoebaSim.amoebaTorsionTorsion_offset)
    {
        unsigned int pos1   = pos - cAmoebaSim.amoebaOutOfPlaneBend_offset;
        if (pos1 < cAmoebaSim.amoebaTorsionTorsions )
        {
            int4   atom1            = cAmoebaSim.pAmoebaTorsionTorsionID1[pos1];  
            int4   atom2            = cAmoebaSim.pAmoebaTorsionTorsionID2[pos1];  

            // atom2.y = chiral atom index
            // atom2.z = grid index

            float4 a1               = cSim.pPosq[atom1.x];
            float4 a2               = cSim.pPosq[atom1.y];
            float4 a3               = cSim.pPosq[atom1.z];
            float4 a4               = cSim.pPosq[atom1.w];
            float4 a5               = cSim.pPosq[atom2.x];

            float xba               = a2.x - a1.x;
            float yba               = a2.y - a1.y;
            float zba               = a2.z - a1.z;

            float xcb               = a3.x - a2.x;
            float ycb               = a3.y - a2.y;
            float zcb               = a3.z - a2.z;

            float xdc               = a4.x - a3.x;
            float ydc               = a4.y - a3.y;
            float zdc               = a4.z - a3.z;

            float xed               = a5.x - a4.x;
            float yed               = a5.y - a4.y;
            float zed               = a5.z - a4.z;

            float xt                = yba*zcb - ycb*zba;
            float yt                = zba*xcb - zcb*xba;
            float zt                = xba*ycb - xcb*yba;

            float xu                = ycb*zdc - ydc*zcb;
            float yu                = zcb*xdc - zdc*xcb;
            float zu                = xcb*ydc - xdc*ycb;

            float rt2               = xt*xt + yt*yt + zt*zt;
            float ru2               = xu*xu + yu*yu + zu*zu;

            float rtru              = sqrtf(rt2 * ru2);

            float xv                = ydc*zed - yed*zdc;
            float yv                = zdc*xed - zed*xdc;
            float zv                = xdc*yed - xed*ydc;

            float rv2               = xv*xv + yv*yv + zv*zv;
            float rurv              = sqrtf(ru2 * rv2);

            float    rcb            = sqrtf(xcb*xcb + ycb*ycb + zcb*zcb);
            float    cosine1        = rtru != 0.0f ? ((xt*xu + yt*yu + zt*zu) / rtru) : 0.0f;
                     cosine1        = cosine1 >  1.0f ?  1.0f : cosine1;
                     cosine1        = cosine1 < -1.0f ? -1.0f : cosine1;

            float    angle1         = LOCAL_HACK_RADIAN * acos(cosine1);
            float    sign           = xba*xu + yba*yu + zba*zu;
                     angle1         = sign < 0.0f ? -angle1 : angle1;
            float    value1         = angle1;

            float    rdc            = sqrtf(xdc*xdc + ydc*ydc + zdc*zdc);
            float    cosine2        = (xu*xv + yu*yv + zu*zv) / rurv;
                     cosine2        = cosine2 >  1.0f ?  1.0f : cosine2;
                     cosine2        = cosine2 < -1.0f ? -1.0f : cosine2;
            float    angle2         = LOCAL_HACK_RADIAN * acos(cosine2);

                     sign           = xcb*xv + ycb*yv + zcb*zv;
                     angle2         = sign < 0.0f ? -angle2 : angle2;
            float    value2         = angle2;

            // check for inverted chirality at the central atom

            // if atom2.y < 0, then no chiral check required
            // sign is set to 1.0 in this case
            // user atom2.x for the atom index to avoid warp divergence

            int chiralAtomIndex     = (atom2.y  > -1) ? atom2.y : atom2.x;
            float4 a6               = cSim.pPosq[chiralAtomIndex];

            float xac               = a6.x - a3.x;
            float yac               = a6.y - a3.y;
            float zac               = a6.z - a3.z;

            float xbc               = a2.x - a3.x;
            float ybc               = a2.y - a3.y;
            float zbc               = a2.z - a3.z;

            // xdc, ydc, zdc appear above

            float xdc1              = a4.x - a3.x;
            float ydc1              = a4.y - a3.y;
            float zdc1              = a4.z - a3.z;

            float c1                = ybc*zdc1 - zbc*ydc1;
            float c2                = ydc1*zac - zdc1*yac;
            float c3                = yac*zbc  - zac*ybc;
            float vol               = xac*c1 + xbc*c2 + xdc1*c3;
                  sign              = vol > 0.0f ? 1.0f : -1.0f;
                  sign              = atom2.y < 0 ? 1.0f : sign;
                  value1           *= sign;
                  value2           *= sign;

           // use bicubic interpolation to compute spline values
           // compute indices into grid based on angles

           int index1               = (int) ((value1 - cAmoebaSim.amoebaTorTorGridBegin[atom2.z])/cAmoebaSim.amoebaTorTorGridDelta[atom2.z] + 1.0e-05f);
           float fIndex             = (float) index1;
           float x1l                = cAmoebaSim.amoebaTorTorGridDelta[atom2.z]*fIndex + cAmoebaSim.amoebaTorTorGridBegin[atom2.z];
           float x1u                = x1l + cAmoebaSim.amoebaTorTorGridDelta[atom2.z];

           int index2               = (int) ((value2 - cAmoebaSim.amoebaTorTorGridBegin[atom2.z])/cAmoebaSim.amoebaTorTorGridDelta[atom2.z] + 1.0e-05f);
                 fIndex             = (float) index2;
           float x2l                = cAmoebaSim.amoebaTorTorGridDelta[atom2.z]*fIndex + cAmoebaSim.amoebaTorTorGridBegin[atom2.z];
           float x2u                = x2l + cAmoebaSim.amoebaTorTorGridDelta[atom2.z];

           int posIndex1            = index2 + index1*cAmoebaSim.amoebaTorTorGridNy[atom2.z];
               posIndex1           += cAmoebaSim.amoebaTorTorGridOffset[atom2.z];

           int posIndex2            = index2 + (index1+1)*cAmoebaSim.amoebaTorTorGridNy[atom2.z];
               posIndex2           += cAmoebaSim.amoebaTorTorGridOffset[atom2.z];

            // load grid points surrounding angle

            float4 y;
            float4 y1;
            float4 y2;
            float4 y12;

            int localIndex         = posIndex1;
                   y.x             = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].x;
                   y1.x            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].y;
                   y2.x            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].z;
                   y12.x           = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].w;

                   localIndex      = posIndex2;
                   y.y             = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].x;
                   y1.y            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].y;
                   y2.y            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].z;
                   y12.y           = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].w;

                   localIndex      = posIndex2 + 1;
                   y.z             = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].x;
                   y1.z            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].y;
                   y2.z            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].z;
                   y12.z           = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].w;

                   localIndex      = posIndex1 + 1;
                   y.w             = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].x;
                   y1.w            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].y;
                   y2.w            = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].z;
                   y12.w           = cAmoebaSim.pAmoebaTorsionTorsionGrids[localIndex].w;

            // perform interpolation

            float    e;
            float    dedang1;
            float    dedang2;
            //float4   cx0,cx1,cx2,cx3;

            bicubic( y, y1, y2, y12, value1, x1l, x1u, value2,  x2l,  x2u, &e, &dedang1, &dedang2 );
            //bicubic( y, y1, y2, y12, value1, x1l, x1u, value2,  x2l,  x2u, &e, &dedang1, &dedang2, &cx0, &cx1, &cx2, &cx3 );
            energy                    += e;
            dedang1                   *=  sign * LOCAL_HACK_RADIAN;
            dedang2                   *=  sign * LOCAL_HACK_RADIAN;

            // chain rule terms  for first angle derivative components

            float    xca               = a3.x - a1.x;
            float    yca               = a3.y - a1.y;
            float    zca               = a3.z - a1.z;

            float    xdb               = a4.x - a2.x;
            float    ydb               = a4.y - a2.y;
            float    zdb               = a4.z - a2.z;

            float    dedxt             =  dedang1 * (yt*zcb - ycb*zt) / (rt2*rcb);
            float    dedyt             =  dedang1 * (zt*xcb - zcb*xt) / (rt2*rcb);
            float    dedzt             =  dedang1 * (xt*ycb - xcb*yt) / (rt2*rcb);
            float    dedxu             = -dedang1 * (yu*zcb - ycb*zu) / (ru2*rcb);
            float    dedyu             = -dedang1 * (zu*xcb - zcb*xu) / (ru2*rcb);
            float    dedzu             = -dedang1 * (xu*ycb - xcb*yu) / (ru2*rcb);

            // compute  first derivative components for first angle

            float    dedxia            = zcb*dedyt - ycb*dedzt;
            float    dedyia            = xcb*dedzt - zcb*dedxt;
            float    dedzia            = ycb*dedxt - xcb*dedyt;

            float    dedxib            = yca*dedzt - zca*dedyt + zdc*dedyu - ydc*dedzu;
            float    dedyib            = zca*dedxt - xca*dedzt + xdc*dedzu - zdc*dedxu;
            float    dedzib            = xca*dedyt - yca*dedxt + ydc*dedxu - xdc*dedyu;

            float    dedxic            = zba*dedyt - yba*dedzt + ydb*dedzu - zdb*dedyu;
            float    dedyic            = xba*dedzt - zba*dedxt + zdb*dedxu - xdb*dedzu;
            float    dedzic            = yba*dedxt - xba*dedyt + xdb*dedyu - ydb*dedxu;

            float    dedxid            = zcb*dedyu - ycb*dedzu;
            float    dedyid            = xcb*dedzu - zcb*dedxu;
            float    dedzid            = ycb*dedxu - xcb*dedyu;

            // chain rule terms  for second angle derivative components

            float    xec               = a5.x - a3.x;
            float    yec               = a5.y - a3.y;
            float    zec               = a5.z - a3.z;

            float    dedxu2            =  dedang2 * (yu*zdc - ydc*zu) / (ru2*rdc);
            float    dedyu2            =  dedang2 * (zu*xdc - zdc*xu) / (ru2*rdc);
            float    dedzu2            =  dedang2 * (xu*ydc - xdc*yu) / (ru2*rdc);
            float    dedxv2            = -dedang2 * (yv*zdc - ydc*zv) / (rv2*rdc);
            float    dedyv2            = -dedang2 * (zv*xdc - zdc*xv) / (rv2*rdc);
            float    dedzv2            = -dedang2 * (xv*ydc - xdc*yv) / (rv2*rdc);

            // compute  first derivative components for second angle

            float    dedxib2           = zdc*dedyu2 - ydc*dedzu2;
            float    dedyib2           = xdc*dedzu2 - zdc*dedxu2;
            float    dedzib2           = ydc*dedxu2 - xdc*dedyu2;
            float    dedxic2           = ydb*dedzu2 - zdb*dedyu2 + zed*dedyv2 - yed*dedzv2;
            float    dedyic2           = zdb*dedxu2 - xdb*dedzu2 + xed*dedzv2 - zed*dedxv2;
            float    dedzic2           = xdb*dedyu2 - ydb*dedxu2 + yed*dedxv2 - xed*dedyv2;
            float    dedxid2           = zcb*dedyu2 - ycb*dedzu2 + yec*dedzv2 - zec*dedyv2;
            float    dedyid2           = xcb*dedzu2 - zcb*dedxu2 + zec*dedxv2 - xec*dedzv2;
            float    dedzid2           = ycb*dedxu2 - xcb*dedyu2 + xec*dedyv2 - yec*dedxv2;
            float    dedxie2           = zdc*dedyv2 - ydc*dedzv2;
            float    dedyie2           = xdc*dedzv2 - zdc*dedxv2;
            float    dedzie2           = ydc*dedxv2 - xdc*dedyv2;

            // increment the torsion-torsion energy and gradient

/*
            float    ett = ett + e
            float    dett(1,ia) = dett(1,ia) + dedxia 
            float    dett(2,ia) = dett(2,ia) + dedyia 
            float    dett(3,ia) = dett(3,ia) + dedzia
            float    dett(1,ib) = dett(1,ib) + dedxib + dedxib2
            float    dett(2,ib) = dett(2,ib) + dedyib + dedyib2
            float    dett(3,ib) = dett(3,ib) + dedzib + dedzib2
            float    dett(1,ic) = dett(1,ic) + dedxic + dedxic2
            float    dett(2,ic) = dett(2,ic) + dedyic + dedyic2 
            float    dett(3,ic) = dett(3,ic) + dedzic + dedzic2
            float    dett(1,id) = dett(1,id) + dedxid + dedxid2
            float    dett(2,id) = dett(2,id) + dedyid + dedyid2
            float    dett(3,id) = dett(3,id) + dedzid + dedzid2
            float    dett(1,ie) = dett(1,ie) + dedxie2 
            float    dett(2,ie) = dett(2,ie) + dedyie2 
            float    dett(3,ie) = dett(3,ie) + dedzie2
*/
            // increment the torsion-torsion gradient;
            int4   atom3            = cAmoebaSim.pAmoebaTorsionTorsionID3[pos1];  

            unsigned int offset     = atom1.x + atom2.w * cSim.stride;
            float4 force            = cSim.pForce4[offset]; 
            force.x                -= dedxia;
            force.y                -= dedyia;
            force.z                -= dedzia;
            cSim.pForce4[offset]    = force;

            offset                  = atom1.y + atom3.x * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= (dedxib + dedxib2);
            force.y                -= (dedyib + dedyib2);
            force.z                -= (dedzib + dedzib2);
            cSim.pForce4[offset]    = force;

            offset                  = atom1.z + atom3.y * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= (dedxic + dedxic2);
            force.y                -= (dedyic + dedyic2);
            force.z                -= (dedzic + dedzic2);
            cSim.pForce4[offset]    = force;

            offset                  = atom1.w + atom3.z * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= (dedxid + dedxid2);
            force.y                -= (dedyid + dedyid2);
            force.z                -= (dedzid + dedzid2);
            cSim.pForce4[offset]    = force;

            offset                  = atom2.x + atom3.w * cSim.stride;
            force                   = cSim.pForce4[offset];
            force.x                -= dedxie2;
            force.y                -= dedyie2;
            force.z                -= dedzie2;
            cSim.pForce4[offset]    = force;

        }
        pos += blockDim.x * gridDim.x;
    }
    cSim.pEnergy[blockIdx.x * blockDim.x + threadIdx.x] += energy;
}


void kCalculateAmoebaLocalForces(amoebaGpuContext gpu)
{
   
    if( gpu->log ){
        static int call = 0;
        if( call == 0 ){
            (void) fprintf( gpu->log,"kCalculateAmoebaLocalForces: blks=%u thrds/blk=%u\n",
                            gpu->gpuContext->sim.blocks, gpu->gpuContext->sim.localForces_threads_per_block); fflush( gpu->log );
            call++;
        }
    }

    kCalculateAmoebaLocalForces_kernel<<<gpu->gpuContext->sim.blocks, gpu->gpuContext->sim.localForces_threads_per_block>>>();
    LAUNCHERROR("kCalculateAmoebaLocalForces");

/*    
    fprintf( stderr, "Done kCalculateAmoebaLocalForces %p \n", gpu->gpuContext->psForce4); fflush( stderr );
    gpu->gpuContext->psForce4->Print( stderr ); fflush( stderr );
    gpu->gpuContext->psEnergy->Print( stderr ); fflush( stderr );
    GetCalculateAmoebaLocalForcesSim( gpu );
    gpu->gpuContext->psForce4->Print( stderr ); fflush( stderr );
    gpu->gpuContext->psEnergy->Print( stderr ); fflush( stderr );
    fprintf( stderr, "Ez %p\n", (float* ) gpu->gpuContext->sim.pEnergy );
*/
    
}

