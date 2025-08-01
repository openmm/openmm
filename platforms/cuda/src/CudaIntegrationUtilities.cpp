/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2021 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "CudaIntegrationUtilities.h"
#include "CudaContext.h"
#include "openmm/common/ContextSelector.h"

using namespace OpenMM;
using namespace std;

#define CHECK_RESULT(result) CHECK_RESULT2(result, errorMessage);
#define CHECK_RESULT2(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<dynamic_cast<CudaContext&>(context).getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
    }

CudaIntegrationUtilities::CudaIntegrationUtilities(CudaContext& context, const System& system) : IntegrationUtilities(context, system),
        ccmaConvergedMemory(NULL) {
        CHECK_RESULT2(cuEventCreate(&ccmaEvent, context.getEventFlags()), "Error creating event for CCMA");
        CHECK_RESULT2(cuMemHostAlloc((void**) &ccmaConvergedMemory, sizeof(int), CU_MEMHOSTALLOC_DEVICEMAP), "Error allocating pinned memory");
        CHECK_RESULT2(cuMemHostGetDevicePointer(&ccmaConvergedDeviceMemory, ccmaConvergedMemory, 0), "Error getting device address for pinned memory");
}

CudaIntegrationUtilities::~CudaIntegrationUtilities() {
    ContextSelector selector(context);
    if (ccmaConvergedMemory != NULL) {
        cuMemFreeHost(ccmaConvergedMemory);
        cuEventDestroy(ccmaEvent);
    }
}

CudaArray& CudaIntegrationUtilities::getPosDelta() {
    return dynamic_cast<CudaContext&>(context).unwrap(posDelta);
}

CudaArray& CudaIntegrationUtilities::getRandom() {
    return dynamic_cast<CudaContext&>(context).unwrap(random);
}

CudaArray& CudaIntegrationUtilities::getStepSize() {
    return dynamic_cast<CudaContext&>(context).unwrap(stepSize);
}

void CudaIntegrationUtilities::applyConstraintsImpl(bool constrainVelocities, double tol) {
    ContextSelector selector(context);
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
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            settleKernel->setArg(1, tol);
        else
            settleKernel->setArg(1, (float) tol);
        settleKernel->execute(settleAtoms.getSize());
    }
    if (shakeAtoms.isInitialized()) {
        if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
            shakeKernel->setArg(1, tol);
        else
            shakeKernel->setArg(1, (float) tol);
        shakeKernel->execute(shakeAtoms.getSize());
    }
    if (ccmaConstraintAtoms.isInitialized()) {
        if (ccmaConstraintAtoms.getSize() <= 1024) {
            // Use the version of CCMA that runs in a single kernel with one workgroup.
            ccmaFullKernel->setArg(0, (int) constrainVelocities);
            if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
                ccmaFullKernel->setArg(14, tol);
            else
                ccmaFullKernel->setArg(14, (float) tol);
            ccmaFullKernel->execute(128, 128);
        }
        else {
            ccmaForceKernel->setArg(6, ccmaConvergedDeviceMemory);
            if (context.getUseDoublePrecision() || context.getUseMixedPrecision())
                ccmaForceKernel->setArg(7, tol);
            else
                ccmaForceKernel->setArg(7, (float) tol);
            ccmaDirectionsKernel->execute(ccmaConstraintAtoms.getSize());
            const int checkInterval = 4;
            ccmaConvergedMemory[0] = 0;
            ccmaUpdateKernel->setArg(4, constrainVelocities ? context.getVelm() : posDelta);
            for (int i = 0; i < 150; i++) {
                ccmaForceKernel->setArg(8, i);
                ccmaForceKernel->execute(ccmaConstraintAtoms.getSize());
                if ((i+1)%checkInterval == 0)
                    CHECK_RESULT2(cuEventRecord(ccmaEvent, 0), "Error recording event for CCMA");
                ccmaMultiplyKernel->setArg(5, i);
                ccmaMultiplyKernel->execute(ccmaConstraintAtoms.getSize());
                ccmaUpdateKernel->setArg(9, i);
                ccmaUpdateKernel->execute(context.getNumAtoms());
                if ((i+1)%checkInterval == 0) {
                    CHECK_RESULT2(cuEventSynchronize(ccmaEvent), "Error synchronizing on event for CCMA");
                    if (ccmaConvergedMemory[0])
                        break;
                }
            }
        }
    }
}

void CudaIntegrationUtilities::distributeForcesFromVirtualSites() {
    ContextSelector selector(context);
    if (numVsites > 0) {
        Vec3 boxVectors[3];
        context.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        mm_double4 recipBoxVectorsDouble[3];
        context.computeReciprocalBoxVectors(recipBoxVectorsDouble);
        if (context.getUseDoublePrecision()) {
            mm_double4 boxVectorsDouble[3];
            for (int i = 0; i < 3; i++)
                boxVectorsDouble[i] = mm_double4(boxVectors[i][0], boxVectors[i][1], boxVectors[i][2], 0);
            vsiteForceKernel->setArg(18, boxVectorsDouble[0]);
            vsiteForceKernel->setArg(19, boxVectorsDouble[1]);
            vsiteForceKernel->setArg(20, boxVectorsDouble[2]);
            vsiteForceKernel->setArg(21, recipBoxVectorsDouble[0]);
            vsiteForceKernel->setArg(22, recipBoxVectorsDouble[1]);
            vsiteForceKernel->setArg(23, recipBoxVectorsDouble[2]);
        }
        else {
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
        }
        for (int i = numVsiteStages-1; i >= 0; i--) {
            vsiteForceKernel->setArg(2, context.getLongForceBuffer());
            vsiteForceKernel->setArg(25, i);
            vsiteForceKernel->execute(numVsites);
        }
    }
}
