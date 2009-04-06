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
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/**
 * This file contains the kernel for evalauating nonbonded forces
 * using the Ewald summation method.
 */

#include <cuComplex.h>

/* Define complex multiply operations */

__device__ cuComplex ComplexMul(cuComplex a, cuComplex b)
{
     cuComplex c;
     c.x = a.x * b.x - a.y * b.y;
     c.y = a.x * b.y + a.y * b.x;
     return c;

}

__device__ cuComplex ComplexConjMul(cuComplex a, cuComplex b)
{
     cuComplex c;
     c.x = a.x*b.x + a.y*b.y;
     c.y = a.y*b.x - a.x*b.y;
     return c;

}

__device__ cuComplex FloatComplexMul(float r, cuComplex a)
{
     cuComplex b;
     b.x = r*a.x;
     b.y = r*a.y;
     return b;

}


/* This kernel is under development */


__global__ void  kCalculateCDLJEwaldForces_kernel(unsigned int* workUnit, int numWorkUnits)
{
/*

    extern __shared__ Atom sA[];
    unsigned int totalWarps = cSim.nonbond_blocks*cSim.nonbond_threads_per_block/GRID;
    unsigned int warp = (blockIdx.x*blockDim.x+threadIdx.x)/GRID;
    int pos = warp*numWorkUnits/totalWarps;
    int end = (warp+1)*numWorkUnits/totalWarps;

//###########################################################################
// EWALD RECIP SPACE
//###########################################################################

// *******************************************************************
    float alphaEwald       =  3.123413f;
    float factorEwald = -1.0 / (4*alphaEwald*alphaEwald);
    float PI = 3.14159265358979323846f;


    float SQRT_PI  = sqrt(PI);
    float TWO_PI   = 2.0f * PI;
    float epsilon  =  1.0f;

/////##############################################################################

    float recipCoeff = 4.0*PI/(cSim.periodicBoxSizeX * cSim.periodicBoxSizeY * cSim.periodicBoxSizeZ) /epsilon;

// setup reciprocal box

    float recipBoxSizeX = TWO_PI / cSim.periodicBoxSizeX;
    float recipBoxSizeY = TWO_PI / cSim.periodicBoxSizeY;
    float recipBoxSizeZ = TWO_PI / cSim.periodicBoxSizeZ;

// setup K-vectors

    unsigned int numRx = 60+1;
    unsigned int numRy = 60+1;
    unsigned int numRz = 60+1;
    const int kmax = 61;

    cuComplex eir[kmax][numberOfAtoms][3];
    cuComplex tab_xy[numberOfAtoms];
    cuComplex tab_qxyz[numberOfAtoms];


  for(unsigned int i = 0; i < numberOfAtoms; i++) {
    for(unsigned int m = 0; (m < 3); m++) {
      eir[0][i][m].x = 1;
      eir[0][i][m].y = 0;
    }


      eir[1][i][0].x = cos(apos.x*recipBoxSizeX);
      eir[1][i][0].y = sin(apos.x*recipBoxSizeX);

      eir[1][i][1].x = cos(apos.y*recipBoxSizeY);
      eir[1][i][1].y = sin(apos.y*recipBoxSizeY);

      eir[1][i][2].x = cos(apos.z*recipBoxSizeZ);
      eir[1][i][2].y = sin(apos.z*recipBoxSizeZ);

    for(unsigned int j=2; (j<kmax); j++) {
      for(unsigned int m=0; (m<3); m++) {
        eir[j][i][m] = ComplexMul (eir[j-1][i][m] , eir[1][i][m]);
      }
    }
  }

// *******************************************************************


    int lowry = 0;
    int lowrz = 1;

    for(int rx = 0; rx < numRx; rx++) {

      float kx = rx * recipBoxSizeX;

      for(int ry = lowry; ry < numRy; ry++) {

        float ky = ry * recipBoxSizeY;

        if(ry >= 0) {
          for(int n = 0; n < numberOfAtoms; n++) {
            tab_xy[n] = ComplexMul (eir[rx][n][0] , eir[ry][n][1]);
          }
        }

        else {
          for(int n = 0; n < numberOfAtoms; n++) {
            tab_xy[n]= ComplexConjMul (eir[rx][n][0] , (eir[-ry][n][1]));
          }
        }

        for (int rz = lowrz; rz < numRz; rz++) {

          float kz = rz * recipBoxSizeZ;
          float k2 = kx * kx + ky * ky + kz * kz;
          float ak = exp(k2*factorEwald) / k2;
          float akv = 2.0 * ak * (1.0/k2 - factorEwald);

          if( rz >= 0) {

           for( int n = 0; n < numberOfAtoms; n++) {
             tab_qxyz[n] = FloatComplexMul ( Charge[n] * ComplexMul (tab_xy[n] , eir[rz][n][2]));
           }
          }

          else {

            for( int n = 0; n < numberOfAtoms; n++) {
              tab_qxyz[n] = FloatComplexMul( Charge[n] * ComplexConjMul (tab_xy[n] , conj(eir[-rz][n][2])) );
            }
          }

          float cs = 0.0f;
          float ss = 0.0f;

          for( int n = 0; n < numberOfAtoms; n++) {
            cs += tab_qxyz[n].x;
            ss += tab_qxyz[n].y;
          }

          recipEnergy += ak * ( cs * cs + ss * ss);
          float vir =  akv * ( cs * cs + ss * ss);

          for(int n = 0; n < numberOfAtoms; n++) {
            float force = ak * (cs * tab_qxyz[n].y - ss * tab_qxyz[n].x);
            forces[n][0] += 2.0 * recipCoeff * force * kx ;
            forces[n][1] += 2.0 * recipCoeff * force * ky ;
            forces[n][2] += 2.0 * recipCoeff * force * kz ;
          } 

          lowrz = 1 - numRz;
        }
        lowry = 1 - numRy;
      }
    }


//###########################################################################


//###########################################################################
// END EWALD RECIP SPACE
//###########################################################################


    while (pos < end)
    {

        // Extract cell coordinates from appropriate work unit
        unsigned int x = workUnit[pos];
        unsigned int y = ((x >> 2) & 0x7fff) << GRIDBITS;
        bool bExclusionFlag = (x & 0x1);
        x = (x >> 17) << GRIDBITS;
        float4      apos;   // Local atom x, y, z, q
        float3      af;     // Local atom fx, fy, fz
        float dx;
        float dy;
        float dz;
        float r2;
        float r;
        float invR;
        float sig;
        float sig2;
        float sig6;
        float eps;
        float dEdR;
        unsigned int tgx = threadIdx.x & (GRID - 1);
        unsigned int tbx = threadIdx.x - tgx;
        int tj = tgx;
        Atom* psA = &sA[tbx];
        unsigned int i      = x + tgx;
        apos                = cSim.pPosq[i];
        float2 a            = cSim.pAttr[i];
        af.x                = 0.0f;
        af.y                = 0.0f;
        af.z                = 0.0f;
        
                int j                   = y + tgx;
                float4 temp             = cSim.pPosq[j];
                float2 temp1            = cSim.pAttr[j];
                sA[threadIdx.x].x       = temp.x;
                sA[threadIdx.x].y       = temp.y;
                sA[threadIdx.x].z       = temp.z;
                sA[threadIdx.x].q       = temp.w;
                sA[threadIdx.x].sig     = temp1.x;
                sA[threadIdx.x].eps     = temp1.y;

            sA[threadIdx.x].fx      = 0.0f;
            sA[threadIdx.x].fy      = 0.0f;
            sA[threadIdx.x].fz      = 0.0f;
            apos.w                 *= cSim.epsfac;

                // Read fixed atom data into registers and GRF
                unsigned int excl       = cSim.pExclusion[x * cSim.exclusionStride + y + tgx];
                excl                    = (excl >> tgx) | (excl << (GRID - tgx));

                for (unsigned int j = 0; j < GRID; j++)
                {
                    dx              = psA[tj].x - apos.x;
                    dy              = psA[tj].y - apos.y;
                    dz              = psA[tj].z - apos.z;

                    dx -= floor(dx/cSim.periodicBoxSizeX+0.5f)*cSim.periodicBoxSizeX;
                    dy -= floor(dy/cSim.periodicBoxSizeY+0.5f)*cSim.periodicBoxSizeY;
                    dz -= floor(dz/cSim.periodicBoxSizeZ+0.5f)*cSim.periodicBoxSizeZ;

                    r2              = dx * dx + dy * dy + dz * dz;
                    r               = sqrt(r2);
                    invR            = 1.0f / sqrt(r2);
                    sig             = a.x + psA[tj].sig;
                    sig2            = invR * sig;
                    sig2           *= sig2;
                    sig6            = sig2 * sig2 * sig2;
                    eps             = a.y * psA[tj].eps;

// ##### LJ
                    dEdR            = eps * (12.0f * sig6 - 6.0f) * sig6;

// ##### SHORT RANGE EWALD
                    float alphaR    = alphaEwald * r;
                    dEdR += apos.w * psA[tj].q * invR * invR * invR * (erfc(alphaR) + 2.0 * alphaR * exp ( - alphaR * alphaR) / SQRT_PI );


                    if (!(excl & 0x1) || r2 > cSim.nonbondedCutoffSqr)
                    {
                        dEdR = 0.0f;
                    }


                    dx             *= dEdR;
                    dy             *= dEdR;
                    dz             *= dEdR;
                    af.x           -= dx;
                    af.y           -= dy;
                    af.z           -= dz;
                    psA[tj].fx     += dx;
                    psA[tj].fy     += dy;
                    psA[tj].fz     += dz;
                    excl          >>= 1;
                    tj              = (tj + 1) & (GRID - 1);
                }

            // Write results
            float4 of;

            of.x                                = af.x;
            of.y                                = af.y;
            of.z                                = af.z;
            of.w                                = 0.0f;
            int offset                          = x + tgx + (y >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = of;
            of.x                                = sA[threadIdx.x].fx;
            of.y                                = sA[threadIdx.x].fy;
            of.z                                = sA[threadIdx.x].fz;
            offset                              = y + tgx + (x >> GRIDBITS) * cSim.stride;
            cSim.pForce4a[offset]               = of;

        pos++;
    }
*/
}
