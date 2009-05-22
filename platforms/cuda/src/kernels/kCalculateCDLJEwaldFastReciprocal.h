/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Rossen P. Apostolov, Peter Eastman                                     *
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

  /* Define multiply operations for floats */
  
  __device__ float2 MultofFloat2(float2 a, float2 b)
  {
       float2 c;
       c.x = a.x * b.x - a.y * b.y;
       c.y = a.x * b.y + a.y * b.x;
       return c;
  
  }
  
  __device__ float2 ConjMultofFloat2(float2 a, float2 b)
  {
       float2 c;
       c.x = a.x*b.x + a.y*b.y;
       c.y = a.y*b.x - a.x*b.y;
       return c;
  
  }
  
  __device__ float2 FloatMultFloat2(float r, float2 a)
  {
       float2 b;
       b.x = r*a.x;
       b.y = r*a.y;
       return b;
  
  }
   
  __device__ float2 FloatMultConjFloat2(float r, float2 a)
  {
       float2 b;
       b.x = r*a.y;
       b.y = r*a.x;
       return b;
  
  }
  
__global__ void kCalculateEwaldFastEikr_kernel()
{

    int kmax               = cSim.kmax;
    float4 apos;

    unsigned int atom    = threadIdx.x + blockIdx.x * blockDim.x;

    while (atom < cSim.atoms)
    {

          apos                         = cSim.pPosq[atom];

//generic form of the array
// pEikr[ atomID*kmax*3 + k*3 + m]


// k = 0, explicitly
        for(unsigned int m = 0; (m < 3); m++) {
          cSim.pEikr[atom*kmax*3 + 0 + m].x = 1;
          cSim.pEikr[atom*kmax*3 + 0 + m].y = 0;
        }
// k = 1, explicitly
          cSim.pEikr[atom*kmax*3 + 3 + 0].x = cos(apos.x*cSim.recipBoxSizeX);
          cSim.pEikr[atom*kmax*3 + 3 + 0].y = sin(apos.x*cSim.recipBoxSizeX);

          cSim.pEikr[atom*kmax*3 + 3 + 1].x = cos(apos.y*cSim.recipBoxSizeY);
          cSim.pEikr[atom*kmax*3 + 3 + 1].y = sin(apos.y*cSim.recipBoxSizeY);

          cSim.pEikr[atom*kmax*3 + 3 + 2].x = cos(apos.z*cSim.recipBoxSizeZ);
          cSim.pEikr[atom*kmax*3 + 3 + 2].y = sin(apos.z*cSim.recipBoxSizeZ);

// k > 1, by recursion
        for(unsigned int k=2; (k<kmax); k++) {
          for(unsigned int m=0; (m<3); m++) {
            cSim.pEikr[atom*kmax*3 + k*3 + m] = MultofFloat2 (cSim.pEikr[atom*kmax*3 + (k-1)*3 + m] , cSim.pEikr[atom*kmax*3 + 3 + m]);
          }
        }

        atom                            += blockDim.x * gridDim.x;
    }
}

__global__ void kCalculateEwaldFastStructureFactors_kernel()
{

// hard-coded maximum k-vectors, no interface yet
    int kmax               = cSim.kmax;

    float4 apos;
    int lowry = 0;
    int lowrz = 1;
    int numRx = 20+1;
    int numRy = 20+1;
    int numRz = 20+1;
    unsigned int totalK  = (numRx*2 - 1) * (numRy*2 - 1) * (numRz*2 - 1);

    float2 tab_xy;
    int index;

    unsigned int atom    = threadIdx.x + blockIdx.x * blockDim.x;

    while (atom < cSim.atoms)
    {

           apos                         = cSim.pPosq[atom];

           // cSim.pEikr[atom*kmax*3 + k*3 + m] 

           for(int rx = 0; rx < numRx; rx++) 
           {

             for(int ry = lowry; ry < numRy; ry++) 
             {

               if(ry >= 0) 
               {
                   tab_xy = MultofFloat2 (cSim.pEikr[atom*kmax*3 + rx*3 + 0] , cSim.pEikr[atom*kmax*3 + ry*3 + 1]);
               }

               else 
               {
                   tab_xy = ConjMultofFloat2 (cSim.pEikr[atom*kmax*3 + rx*3 + 0] , cSim.pEikr[atom*kmax*3 - ry*3 + 1]);
               }

               for (int rz = lowrz; rz < numRz; rz++) 
               {

                   index = rx * (numRy*2 - 1 ) * (numRz*2 - 1) + (ry + numRy -1 ) * (numRz * 2 - 1) + (rz + numRz -1 );

                   if( rz >= 0) 
                   {
                      cSim.pStructureFactor[atom*totalK + index] = FloatMultFloat2 ( (apos.w ), MultofFloat2 (tab_xy ,cSim.pEikr[atom*kmax*3 + rz*3 + 2] ));
                   }

                   else 
                   {
                      cSim.pStructureFactor[atom*totalK + index] = FloatMultFloat2 ( ( apos.w ), ConjMultofFloat2 (tab_xy ,cSim.pEikr[atom*kmax*3 - rz*3 + 2] ));
                   }

                     cSim.pCosSinSum[index].x = 0.0;
                     cSim.pCosSinSum[index].y = 0.0;
       
                   lowrz = 1 - numRz;
               }

               lowry = 1 - numRy;
             }

           }

        atom                            += blockDim.x * gridDim.x;
    }
}


__global__ void kCalculateEwaldFastCosSinSums_kernel()
{

//    float2 eikr;
    int lowry = 0;
    int lowrz = 1;
    int numRx = 20+1;
    int numRy = 20+1;
    int numRz = 20+1;
    unsigned int totalK  = (numRx*2 - 1) * (numRy*2 - 1) * (numRz*2 - 1);

    int index;

    unsigned int rx    = threadIdx.x + blockIdx.x * blockDim.x;

    while (rx < numRx)
    {
// **********************************************************************

// cSim.pEikr[atom*kmax*3 + k*3 + m] 

//    for(int rx = 0; rx < numRx; rx++) {

      for(int ry = lowry; ry < numRy; ry++) {

        for (int rz = lowrz; rz < numRz; rz++) {

          index = rx * (numRy*2 - 1 ) * (numRz*2 - 1) + (ry + numRy -1 ) * (numRz * 2 - 1) + (rz + numRz -1 );

             for (int atom = 0; atom < cSim.atoms; atom++)
             {
                 cSim.pCosSinSum[index].x += cSim.pStructureFactor[atom*totalK + index].x;
                 cSim.pCosSinSum[index].y += cSim.pStructureFactor[atom*totalK + index].y;
             }

          lowrz = 1 - numRz;
        }
        lowry = 1 - numRy;
      }

        rx                            += blockDim.x * gridDim.x;
    }

}

__global__ void kCalculateEwaldFastForces_kernel()
{

    float PI               = 3.14159265358979323846f;
    const float epsilon     =  1.0;
    float recipCoeff = (4*PI/cSim.V/epsilon);

    int lowry = 0;
    int lowrz = 1;
    int numRx = 20+1;
    int numRy = 20+1;
    int numRz = 20+1;
    unsigned int totalK  = (numRx*2 - 1) * (numRy*2 - 1) * (numRz*2 - 1);
    int index;

    unsigned int atom    = threadIdx.x + blockIdx.x * blockDim.x;

    while (atom < cSim.atoms)
    {

    for(int rx = 0; rx < numRx; rx++) {

      float kx = rx * cSim.recipBoxSizeX;

      for(int ry = lowry; ry < numRy; ry++) {

      float ky = ry * cSim.recipBoxSizeY;

        for (int rz = lowrz; rz < numRz; rz++) {

        float kz = rz * cSim.recipBoxSizeZ;

       // next one is scary!
          index = rx * (numRy*2 - 1 ) * (numRz*2 - 1) + (ry + numRy -1 ) * (numRz * 2 - 1) + (rz + numRz -1 );

          float k2 = kx * kx + ky * ky + kz * kz;
          float ak = exp(k2*cSim.factorEwald) / k2;

          float dEdR = ak * (cSim.pCosSinSum[index].x  * cSim.pStructureFactor[atom*totalK + index].y - cSim.pCosSinSum[index].y * cSim.pStructureFactor[atom*totalK + index].x);
 
          cSim.pForce4[atom].x += 2 * recipCoeff * dEdR * kx ;
          cSim.pForce4[atom].y += 2 * recipCoeff * dEdR * ky ;
          cSim.pForce4[atom].z += 2 * recipCoeff * dEdR * kz ;

          lowrz = 1 - numRz;
        }
        lowry = 1 - numRy;
      }
    }



        atom                            += blockDim.x * gridDim.x;
    }
}
