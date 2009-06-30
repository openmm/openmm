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
    const float eps0 = 1.0f/(4.0f*3.1415926535f*cSim.epsfac);

    unsigned int atomID1    = threadIdx.x + blockIdx.x * blockDim.x;

    while (atomID1 < cSim.atoms)
    {
        float4 apos1 = cSim.pPosq[atomID1];
        float4 af    = cSim.pForce4[atomID1];
        unsigned int atomID2 = 0;
        while (atomID2 < cSim.atoms)
        {
            float4 apos2 = cSim.pPosq[atomID2];
            float scale = 2.0f*apos1.w*apos2.w/(cSim.cellVolume*eps0);

            int lowry = 0;
            int lowrz = 1;

            for(int rx = 0; rx < cSim.kmax; rx++)
            {
                float kx = rx*cSim.recipBoxSizeX;
                for(int ry = lowry; ry < cSim.kmax; ry++)
                {
                    float ky = ry*cSim.recipBoxSizeY;
                    for (int rz = lowrz; rz < cSim.kmax; rz++)
                    {
                        float kz = rz*cSim.recipBoxSizeZ;
                        float k2 = kx*kx + ky*ky + kz*kz;
                        float ek = exp(k2*cSim.factorEwald);

                        float arg1 = kx*apos1.x + ky*apos1.y + kz*apos1.z;
                        float arg2 = kx*apos2.x + ky*apos2.y + kz*apos2.z;
                        float sinI = sin(arg1);
                        float sinJ = sin(arg2);
                        float cosI = cos(arg1);
                        float cosJ = cos(arg2);

                        float f = scale * ek * (-sinI*cosJ + cosI*sinJ) / k2;
                        af.x -= kx*f;
                        af.y -= ky*f;
                        af.z -= kz*f;

                        lowrz = 1 - cSim.kmax;
                    }
                    lowry = 1 - cSim.kmax;
                }
            }
            atomID2++;
       }
       cSim.pForce4[atomID1] = af;
       atomID1 += blockDim.x * gridDim.x;

    }
}
