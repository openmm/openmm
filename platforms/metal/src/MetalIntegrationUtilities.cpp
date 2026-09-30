/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform; long-only virtual-site scheduling
 * follows platforms/cuda/src/CudaIntegrationUtilities.cpp.
 * Source: platforms/opencl/src/OpenCLIntegrationUtilities.cpp
 *
 * Original OpenCL Platform code:
 * Portions copyright (c) 2009-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *
 * Metal Platform code:
 * Portions copyright (c) 2026 Chun-Chi Hung.
 * Authors: Chun-Chi Hung
 *
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

#include "MetalIntegrationUtilities.h"
#include "MetalContext.h"

using namespace OpenMM;
using namespace std;

MetalIntegrationUtilities::MetalIntegrationUtilities(MetalContext& context, const System& system) : IntegrationUtilities(context, system) {
    ccmaConvergedHostBuffer.initialize<int>(context, 1, "ccmaConvergedHostBuffer");
}

MetalArray& MetalIntegrationUtilities::getPosDelta() {
    return dynamic_cast<MetalContext&>(context).unwrap(posDelta);
}

MetalArray& MetalIntegrationUtilities::getRandom() {
    return dynamic_cast<MetalContext&>(context).unwrap(random);
}

MetalArray& MetalIntegrationUtilities::getStepSize() {
    return dynamic_cast<MetalContext&>(context).unwrap(stepSize);
}

void MetalIntegrationUtilities::applyConstraintsImpl(bool constrainVelocities, double tol) {
    ComputeKernel settleKernel, shakeKernel, ccmaForceKernel;
    if (constrainVelocities) {
        settleKernel = settleVelKernel;
        shakeKernel = shakeVelKernel;
        ccmaForceKernel = ccmaVelForceKernel;
    }
    else {
        settleKernel = settlePosKernel;
        shakeKernel = shakePosKernel;
        ccmaForceKernel = ccmaPosForceKernel;
    }
    if (settleAtoms.isInitialized()) {
        settleKernel->setArg(1, (float) tol);
        settleKernel->execute(settleAtoms.getSize());
    }
    if (shakeAtoms.isInitialized()) {
        shakeKernel->setArg(1, (float) tol);
        shakeKernel->execute(shakeAtoms.getSize());
    }
    if (ccmaConstraintAtoms.isInitialized()) {
        if (ccmaConstraintAtoms.getSize() <= 1024) {
            // Use the version of CCMA that runs in a single kernel with one workgroup.
            ccmaFullKernel->setArg(0, (int) constrainVelocities);
            ccmaFullKernel->setArg(14, (float) tol);
            ccmaFullKernel->execute(128, 128);
        }
        else {
            // Use the version of CCMA that uses multiple kernels.
            ccmaForceKernel->setArg(6, ccmaConvergedHostBuffer);
            ccmaForceKernel->setArg(7, (float) tol);
            ccmaDirectionsKernel->execute(ccmaConstraintAtoms.getSize());
            const int checkInterval = 4;
            int converged[2] = {0, 0};
            int zero = 0;
            ccmaConvergedHostBuffer.upload(&zero);
            ccmaUpdateKernel->setArg(4, constrainVelocities ? context.getVelm() : posDelta);
            for (int i = 0; i < 150; i++) {
                ccmaForceKernel->setArg(8, i);
                ccmaForceKernel->execute(ccmaConstraintAtoms.getSize());
                ccmaMultiplyKernel->setArg(5, i);
                ccmaMultiplyKernel->execute(ccmaConstraintAtoms.getSize());
                ccmaUpdateKernel->setArg(9, i);
                ccmaUpdateKernel->execute(context.getNumAtoms());
                if ((i+1)%checkInterval == 0) {
                    // Read only the convergence flags; all constraint updates execute on the GPU.
                    ccmaConverged.download(converged);
                    if (converged[i%2])
                        break;
                }
            }
        }
    }
}

void MetalIntegrationUtilities::distributeForcesFromVirtualSites() {
    if (numVsites > 0) {
        Vec3 boxVectors[3];
        context.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        mm_double4 recipBoxVectorsDouble[3];
        context.computeReciprocalBoxVectors(recipBoxVectorsDouble);
        mm_float4 boxVectorsFloat[3], recipBoxVectorsFloat[3];
        for (int i = 0; i < 3; i++) {
            boxVectorsFloat[i] = mm_float4((float) boxVectors[i][0], (float) boxVectors[i][1], (float) boxVectors[i][2], 0);
            recipBoxVectorsFloat[i] = mm_float4((float) recipBoxVectorsDouble[i].x, (float) recipBoxVectorsDouble[i].y, (float) recipBoxVectorsDouble[i].z, 0);
        }
        vsiteForceKernel->setArg(18, boxVectorsFloat[0]);
        vsiteForceKernel->setArg(19, boxVectorsFloat[1]);
        vsiteForceKernel->setArg(20, boxVectorsFloat[2]);
        vsiteForceKernel->setArg(21, recipBoxVectorsFloat[0]);
        vsiteForceKernel->setArg(22, recipBoxVectorsFloat[1]);
        vsiteForceKernel->setArg(23, recipBoxVectorsFloat[2]);
        for (int i = numVsiteStages-1; i >= 0; i--) {
            vsiteForceKernel->setArg(2, context.getLongForceBuffer());
            vsiteForceKernel->setArg(25, i);
            vsiteForceKernel->execute(numVsites);
        }
    }
}
