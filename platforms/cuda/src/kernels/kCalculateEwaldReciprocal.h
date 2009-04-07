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

    unsigned int pos    = threadIdx.x + blockIdx.x * blockDim.x;

    cuComplex eir[kmax][cSim.atoms][3];
    cuComplex tab_xy[cSim.atoms];
    cuComplex tab_qxyz[cSim.atoms];


  while (pos < cSim.atoms)
  {

      float4 apos             = cSim.pPosq[pos];

      for(unsigned int m = 0; (m < 3); m++) {
        eir[0][pos][m].x = 1;
        eir[0][pos][m].y = 0;
      }


      eir[1][pos][0].x = cos(apos.x*recipBoxSizeX);
      eir[1][pos][0].y = sin(apos.x*recipBoxSizeX);

      eir[1][pos][1].x = cos(apos.y*recipBoxSizeY);
      eir[1][pos][1].y = sin(apos.y*recipBoxSizeY);

      eir[1][pos][2].x = cos(apos.z*recipBoxSizeZ);
      eir[1][pos][2].y = sin(apos.z*recipBoxSizeZ);

      for(unsigned int j=2; (j<kmax); j++) {
        for(unsigned int m=0; (m<3); m++) {
          eir[j][pos][m] = ComplexMul (eir[j-1][pos][m] , eir[1][pos][m]);
        }
      }

      pos                    += blockDim.x * gridDim.x;
  }

// *******************************************************************


    int lowry = 0;
    int lowrz = 1;

    for(int rx = 0; rx < numRx; rx++) {

      float kx = rx * recipBoxSizeX;

      for(int ry = lowry; ry < numRy; ry++) {

        float ky = ry * recipBoxSizeY;

        if(ry >= 0) {
          while (pos < cSim.atoms)
          {
            tab_xy[pos] = ComplexMul (eir[rx][pos][0] , eir[ry][pos][1]);
            pos                    += blockDim.x * gridDim.x;
          }
        }

        else {
          while (pos < cSim.atoms)
          {
            tab_xy[pos]= ComplexConjMul (eir[rx][pos][0] , (eir[-ry][pos][1]));
            pos                    += blockDim.x * gridDim.x;
          }
        }

        for (int rz = lowrz; rz < numRz; rz++) {

          float kz = rz * recipBoxSizeZ;
          float k2 = kx * kx + ky * ky + kz * kz;
          float ak = exp(k2*factorEwald) / k2;
          float akv = 2.0 * ak * (1.0/k2 - factorEwald);

          if( rz >= 0) {

           while (pos < cSim.atoms)
           {
             float4 apos         = cSim.pPosq[pos];
             apos.w             *= cSim.epsfac;
             tab_qxyz[pos] = FloatComplexMul ( apos.w * ComplexMul (tab_xy[pos] , eir[rz][pos][2]));
             pos                    += blockDim.x * gridDim.x;
           }
          }

          else {

            while (pos < cSim.atoms)
            {
              float4 apos         = cSim.pPosq[pos];
              apos.w             *= cSim.epsfac;
              tab_qxyz[pos] = FloatComplexMul( apos.w * ComplexConjMul (tab_xy[pos] , conj(eir[-rz][pos][2])) );
              pos                    += blockDim.x * gridDim.x;
            }
          }

          float cs = 0.0f;
          float ss = 0.0f;

          while (pos < cSim.atoms)
          {
            cs += tab_qxyz[pos].x;
            ss += tab_qxyz[pos].y;
            pos                    += blockDim.x * gridDim.x;
          }

          recipEnergy += ak * ( cs * cs + ss * ss);
          float vir =  akv * ( cs * cs + ss * ss);

          while (pos < cSim.atoms)
          {
            float4 force            = cSim.pForce4[pos];

            float dEdR = ak * (cs * tab_qxyz[pos].y - ss * tab_qxyz[pos].x);
            force.x += 2.0 * recipCoeff * dEdR * kx ;
            force.y += 2.0 * recipCoeff * dEdR * ky ;
            force.z += 2.0 * recipCoeff * dEdR * kz ;
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

}
