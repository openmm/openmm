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

/* Cuda compiler on Windows does not recognized "static const float" values */
#define LOCAL_HACK_PI 3.1415926535897932384626433832795

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

/**
 * Precompute the cosine and sine sums which appear in each force term.
 */

__global__ void kCalculateEwaldFastCosSinSums_kernel()
{
    const float epsilon =  1.0;
    const float recipCoeff = cSim.epsfac*(4*LOCAL_HACK_PI/cSim.cellVolume/epsilon);
    const unsigned int ksizex = 2*cSim.kmaxX-1;
    const unsigned int ksizey = 2*cSim.kmaxY-1;
    const unsigned int ksizez = 2*cSim.kmaxZ-1;
    const unsigned int totalK = ksizex*ksizey*ksizez;
    unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;
    float energy = 0.0f;
    while (index < (cSim.kmaxY-1)*ksizez+cSim.kmaxZ)
        index += blockDim.x * gridDim.x;
    while (index < totalK)
    {
        // Find the wave vector (kx, ky, kz) this index corresponds to.

        int rx = index/(ksizey*ksizez);
        int remainder = index - rx*ksizey*ksizez;
        int ry = remainder/ksizez;
        int rz = remainder - ry*ksizez - cSim.kmaxZ + 1;
        ry += -cSim.kmaxY + 1;
        float kx = rx*cSim.recipBoxSizeX;
        float ky = ry*cSim.recipBoxSizeY;
        float kz = rz*cSim.recipBoxSizeZ;

        // Compute the sum for this wave vector.

        float2 sum = make_float2(0.0f, 0.0f);
        for (int atom = 0; atom < cSim.atoms; atom++)
        {
            float4 apos = cSim.pPosq[atom];
            float phase = apos.x*kx;
            float2 structureFactor = make_float2(cos(phase), sin(phase));
            phase = apos.y*ky;
            structureFactor = MultofFloat2(structureFactor, make_float2(cos(phase), sin(phase)));
            phase = apos.z*kz;
            structureFactor = MultofFloat2(structureFactor, make_float2(cos(phase), sin(phase)));
            sum.x += apos.w*structureFactor.x;
            sum.y += apos.w*structureFactor.y;
        }
        cSim.pEwaldCosSinSum[index] = sum;

        // Compute the contribution to the energy.

        float k2 = kx*kx + ky*ky + kz*kz;
        float ak = exp(k2*cSim.factorEwald) / k2;
        energy += recipCoeff*ak*(sum.x*sum.x + sum.y*sum.y);
        index += blockDim.x * gridDim.x;
    }
    cSim.pEnergy[blockIdx.x*blockDim.x+threadIdx.x] += energy;
}

/**
 * Compute the reciprocal space part of the Ewald force, using the precomputed sums from the
 * previous routine.
 */

__global__ void kCalculateEwaldFastForces_kernel()
{
    const float epsilon =  1.0;
    float recipCoeff = cSim.epsfac*(4*LOCAL_HACK_PI/cSim.cellVolume/epsilon);

    unsigned int atom = threadIdx.x + blockIdx.x * blockDim.x;

    while (atom < cSim.atoms)
    {
        float4 force = cSim.pForce4[atom];
        float4 apos = cSim.pPosq[atom];

        // Loop over all wave vectors.

        int lowry = 0;
        int lowrz = 1;
        for (int rx = 0; rx < cSim.kmaxX; rx++) {
            float kx = rx * cSim.recipBoxSizeX;
            for (int ry = lowry; ry < cSim.kmaxY; ry++) {
                float ky = ry * cSim.recipBoxSizeY;
                float phase = apos.x*kx;
                float2 tab_xy = make_float2(cos(phase), sin(phase));
                phase = apos.y*ky;
                tab_xy = MultofFloat2(tab_xy, make_float2(cos(phase), sin(phase)));
                for (int rz = lowrz; rz < cSim.kmaxZ; rz++) {
                    float kz = rz * cSim.recipBoxSizeZ;

                    // Compute the force contribution of this wave vector.

                    int index = rx*(cSim.kmaxY*2-1)*(cSim.kmaxZ*2-1) + (ry+cSim.kmaxY-1)*(cSim.kmaxZ*2-1) + (rz+cSim.kmaxZ-1);
                    float k2 = kx*kx + ky*ky + kz*kz;
                    float ak = exp(k2*cSim.factorEwald) / k2;
                    phase = apos.z*kz;
                    float2 structureFactor = MultofFloat2(tab_xy, make_float2(cos(phase), sin(phase)));
                    float2 cosSinSum = cSim.pEwaldCosSinSum[index];
                    float dEdR = ak * apos.w * (cosSinSum.x*structureFactor.y - cosSinSum.y*structureFactor.x);
                    force.x += 2 * recipCoeff * dEdR * kx;
                    force.y += 2 * recipCoeff * dEdR * ky;
                    force.z += 2 * recipCoeff * dEdR * kz;
                    lowrz = 1 - cSim.kmaxZ;
                }
                lowry = 1 - cSim.kmaxY;
            }
        }

        // Record the force on the atom.

        cSim.pForce4[atom] = force;
        atom += blockDim.x * gridDim.x;
    }
}
