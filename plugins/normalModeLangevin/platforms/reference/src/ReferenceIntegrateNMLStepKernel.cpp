/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Chris Sweet                                                       *
 * Contributors: Christopher Bruns                                            *
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

#include "ReferenceIntegrateNMLStepKernel.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKReference/ReferenceCCMAAlgorithm.h"
#include <vector>

using namespace std;

namespace OpenMM {

// These methods lifted from ReferenceKernels.cpp, so NormalModeLangevin could be a plugin
//   allocateIntArray
//   disposeIntArray
//   extractPositions
//   extractVelocities
//   extractForces
//   findAnglesForCCMA

static int** allocateIntArray(int length, int width) {
    int** array = new int*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new int[width];
    return array;
}

static void disposeIntArray(int** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static RealOpenMM** extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->positions;
}

static RealOpenMM** extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->velocities;
}

static RealOpenMM** extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->forces;
}

static void findAnglesForCCMA(const System& system, vector<ReferenceCCMAAlgorithm::AngleInfo>& angles) {
    for (int i = 0; i < system.getNumForces(); i++) {
        const HarmonicAngleForce* force = dynamic_cast<const HarmonicAngleForce*>(&system.getForce(i));
        if (force != NULL) {
            for (int j = 0; j < force->getNumAngles(); j++) {
                int atom1, atom2, atom3;
                double angle, k;
                force->getAngleParameters(j, atom1, atom2, atom3, angle, k);
                angles.push_back(ReferenceCCMAAlgorithm::AngleInfo(atom1, atom2, atom3, (RealOpenMM)angle));
            }
        }
    }
}


ReferenceIntegrateNMLStepKernel::~ReferenceIntegrateNMLStepKernel() {
    if (dynamics)
        delete dynamics;
    if (constraints)
        delete constraints;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (constraintDistances)
        delete[] constraintDistances;
    if (projectionVectors)
        delete projectionVectors;
}

void ReferenceIntegrateNMLStepKernel::initialize(const System& system, const NMLIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    constraintDistances = new RealOpenMM[numConstraints];
    for (int i = 0; i < numConstraints; ++i) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraintIndices[i][0] = particle1;
        constraintIndices[i][1] = particle2;
        constraintDistances[i] = static_cast<RealOpenMM>(distance);
    }
    //SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateNMLStepKernel::execute(ContextImpl& context, const NMLIntegrator& integrator, const double currentPE, const int stepType) {

    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    double* dProjectionVectors = integrator.getProjectionVectors();
    unsigned int numProjectionVectors = integrator.getNumProjectionVectors();
    bool projVecChanged = integrator.getProjVecChanged();
    double minimumLimit = integrator.getMinimumLimit();
    double maxEig = integrator.getMaxEigenvalue();


    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** velData = extractVelocities(context);
    RealOpenMM** forceData = extractForces(context);

    //projection vectors
    if (projectionVectors == 0 || projVecChanged) {
        std::cout << "Setting vectors " << projVecChanged << std::endl;
        unsigned int arraySz = numProjectionVectors * context.getSystem().getNumParticles() * 3;
        if(projectionVectors == 0) {
            projectionVectors = new RealOpenMM[arraySz];
        }
        for(unsigned int i=0; i<arraySz; i++) {
            projectionVectors[i] = static_cast<RealOpenMM>(dProjectionVectors[i]);
        }

    }


    if (dynamics == 0 ){ //|| temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.

        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );

        dynamics = new ReferenceNMLDynamics(context.getSystem().getNumParticles(),
                                            static_cast<RealOpenMM>(stepSize),
                                            static_cast<RealOpenMM>(tau),

                                            static_cast<RealOpenMM>(temperature),
                                            projectionVectors,
                                            numProjectionVectors,
                                            static_cast<RealOpenMM>(minimumLimit),
                                            static_cast<RealOpenMM>(maxEig));

        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(context.getSystem(), angles);
        constraints = new ReferenceCCMAAlgorithm(context.getSystem().getNumParticles(), numConstraints, constraintIndices, constraintDistances, masses, angles, (RealOpenMM)integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    //for minimization quadratic fit we need the potential energy
    RealOpenMM currentPE2 = static_cast<RealOpenMM>(currentPE);
    //
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses, currentPE2, stepType );
    //update at dynamic step 2
    if(stepType == 2){
        data.time += stepSize;
        data.stepCount++;
    }

}

} /* namespace OpenMM */

