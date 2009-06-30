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

#define EIR(x, y, z) cSim.pEwaldEikr[(x)+(y)*cSim.kmax+(z)*cSim.kmax*3]

__global__ void kCalculateEwaldFastEikr_kernel()
{

    int kmax = cSim.kmax;
    float4 apos;

    unsigned int atom    = threadIdx.x + blockIdx.x * blockDim.x;

    while (atom < cSim.atoms)
    {

          apos                         = cSim.pPosq[atom];

//generic form of the array
// pEikr[ atomID*kmax*3 + k*3 + m]


// k = 0, explicitly
        for(unsigned int m = 0; (m < 3); m++) {
            EIR(atom, 0, m).x = 1;
            EIR(atom, 0, m).y = 0;
        }
// k = 1, explicitly
          EIR(atom, 1, 0).x = cos(apos.x*cSim.recipBoxSizeX);
          EIR(atom, 1, 0).y = sin(apos.x*cSim.recipBoxSizeX);

          EIR(atom, 1, 1).x = cos(apos.y*cSim.recipBoxSizeY);
          EIR(atom, 1, 1).y = sin(apos.y*cSim.recipBoxSizeY);

          EIR(atom, 1, 2).x = cos(apos.z*cSim.recipBoxSizeZ);
          EIR(atom, 1, 2).y = sin(apos.z*cSim.recipBoxSizeZ);

// k > 1, by recursion
        for(unsigned int k=2; (k<kmax); k++) {
          for(unsigned int m=0; (m<3); m++) {
            EIR(atom, k, m) = MultofFloat2(EIR(atom, k-1, m), EIR(atom, 1, m));
          }
        }

        atom                            += blockDim.x * gridDim.x;
    }
}

__global__ void kCalculateEwaldFastCosSinSums_kernel()
{
    const unsigned int ksize = 2*cSim.kmax-1;
    const unsigned int totalK = ksize*ksize*ksize;
    unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;
    while (index < totalK)
    {
        int rx = index/(ksize*ksize);
        int remainder = index - rx*ksize*ksize;
        int ry = remainder/ksize;
        int rz = remainder - ry*ksize - cSim.kmax + 1;
        ry += -cSim.kmax + 1;
        float2 sum = make_float2(0.0f, 0.0f);
        for (int atom = 0; atom < cSim.atoms; atom++)
        {
            float2 tab_xy = (ry >= 0 ? MultofFloat2(EIR(atom, rx, 0), EIR(atom, ry, 1)) :
                                       ConjMultofFloat2(EIR(atom, rx, 0), EIR(atom, -ry, 1)));
            float charge = cSim.pPosq[atom].w;
            float2 structureFactor = (rz >= 0 ? MultofFloat2(tab_xy, EIR(atom, rz, 2)) :
                                                ConjMultofFloat2(tab_xy, EIR(atom, -rz, 2)));
            sum.x += charge*structureFactor.x;
            sum.y += charge*structureFactor.y;
        }
        cSim.pEwaldCosSinSum[index] = sum;
        index += blockDim.x * gridDim.x;
    }
}

__global__ void kCalculateEwaldFastForces_kernel()
{

    float PI = 3.14159265358979323846f;
    const float epsilon =  1.0;
    float recipCoeff = cSim.epsfac*(4*PI/cSim.cellVolume/epsilon);

    int lowry = 0;
    int lowrz = 1;
    const int numRx = cSim.kmax;
    const int numRy = cSim.kmax;
    const int numRz = cSim.kmax;

    unsigned int atom = threadIdx.x + blockIdx.x * blockDim.x;

    while (atom < cSim.atoms)
    {
        float4 force = cSim.pForce4[atom];
        float charge = cSim.pPosq[atom].w;
        for (int rx = 0; rx < numRx; rx++) {
            float kx = rx * cSim.recipBoxSizeX;
            for (int ry = lowry; ry < numRy; ry++) {
                float ky = ry * cSim.recipBoxSizeY;
                float2 tab_xy = (ry >= 0 ? MultofFloat2(EIR(atom, rx, 0), EIR(atom, ry, 1)) :
                                           ConjMultofFloat2(EIR(atom, rx, 0), EIR(atom, -ry, 1)));
                for (int rz = lowrz; rz < numRz; rz++) {
                    float kz = rz * cSim.recipBoxSizeZ;
                    int index = rx*(numRy*2-1)*(numRz*2-1) + (ry+numRy-1)*(numRz*2-1) + (rz+numRz-1);
                    float k2 = kx*kx + ky*ky + kz*kz;
                    float ak = exp(k2*cSim.factorEwald) / k2;
                    float2 structureFactor = (rz >= 0 ? MultofFloat2(tab_xy, EIR(atom, rz, 2)) :
                                                        ConjMultofFloat2(tab_xy, EIR(atom, -rz, 2)));
                    float2 cosSinSum = cSim.pEwaldCosSinSum[index];
                    float dEdR = ak * charge * (cosSinSum.x*structureFactor.y - cosSinSum.y*structureFactor.x);
                    force.x += 2 * recipCoeff * dEdR * kx;
                    force.y += 2 * recipCoeff * dEdR * ky;
                    force.z += 2 * recipCoeff * dEdR * kz;
                    lowrz = 1 - numRz;
                }
                lowry = 1 - numRy;
            }
        }
        cSim.pForce4[atom] = force;
        atom += blockDim.x * gridDim.x;
    }
}
