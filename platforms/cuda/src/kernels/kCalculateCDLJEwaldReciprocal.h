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

/**
 * This file contains the kernel for evaluating nonbonded forces using the
 * Ewald summation method (Reciprocal space summation).
 */


__global__ void kCalculateCDLJEwaldReciprocalForces_kernel()
{
    float alphaEwald       = cSim.alphaEwald;
    float eps0             = 5.72765E-4;

    int numRx              = 20+1;
    int numRy              = 20+1;
    int numRz              = 20+1;

    int lowry, lowrz;

    float kx, ky, kz, k2, ek;
    float Qi, Qj, SinI, SinJ, CosI, CosJ;

    float V = cSim.cellVolume;

    float4 apos1, apos2 ;

    unsigned int atomID1    = threadIdx.x + blockIdx.x * blockDim.x;

    while (atomID1 < cSim.stride * cSim.outputBuffers)
    {
        apos1             = cSim.pPosq[atomID1];

        unsigned int atomID2    = 0;
        while (atomID2 < cSim.atoms)
        {
                apos2             = cSim.pPosq[atomID2];

                lowry = 0;
                lowrz = 1;

                for(int rx = 0; rx < numRx; rx++) {

                  kx = rx * cSim.recipBoxSizeX;

                  for(int ry = lowry; ry < numRy; ry++) {

                    ky = ry * cSim.recipBoxSizeY;

                    for (int rz = lowrz; rz < numRz; rz++) {

                      kz = rz * cSim.recipBoxSizeZ;

                      k2 = kx * kx + ky * ky + kz * kz;
                      ek = exp (  k2 * cSim.factorEwald);

                      Qi = apos1.w ;
                      Qj = apos2.w ;

                      SinI = sin ( kx * apos1.x + ky * apos1.y + kz * apos1.z );
                      SinJ = sin ( kx * apos2.x + ky * apos2.y + kz * apos2.z );
                      CosI = cos ( kx * apos1.x + ky * apos1.y + kz * apos1.z );
                      CosJ = cos ( kx * apos2.x + ky * apos2.y + kz * apos2.z );

                      cSim.pForce4[atomID1].x -= (2.0 / (V * eps0 ))  * Qi * ( kx/k2) * ek * ( - SinI * Qj * CosJ + CosI * Qj * SinJ);
                      cSim.pForce4[atomID1].y -= (2.0 / (V * eps0 ))  * Qi * ( ky/k2) * ek * ( - SinI * Qj * CosJ + CosI * Qj * SinJ);
                      cSim.pForce4[atomID1].z -= (2.0 / (V * eps0 ))  * Qi * ( kz/k2) * ek * ( - SinI * Qj * CosJ + CosI * Qj * SinJ);

                      lowrz = 1 - numRz;
                    }

                    lowry = 1 - numRy;
                  }
                }

                atomID2++;
       }

       atomID1                            += blockDim.x * gridDim.x;

    }
}
