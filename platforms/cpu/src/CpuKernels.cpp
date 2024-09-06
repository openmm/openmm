/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2024 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

#include "CpuKernels.h"
#include "ReferenceAngleBondIxn.h"
#include "ReferenceBondForce.h"
#include "ReferenceConstraints.h"
#include "ReferenceKernelFactory.h"
#include "ReferenceKernels.h"
#include "ReferenceLJCoulomb14.h"
#include "ReferenceProperDihedralBond.h"
#include "ReferenceRbDihedralBond.h"
#include "ReferenceTabulatedFunction.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/Vec3.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/internal/vectorize.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include <iostream>
#include "lepton/ParsedExpression.h"

using namespace OpenMM;
using namespace std;

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->positions;
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->velocities;
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->forces;
}

static Vec3& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->periodicBoxSize;
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return data->periodicBoxVectors;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->constraints;
}

static const ReferenceVirtualSites& extractVirtualSites(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->virtualSites;
}

static map<string, double>& extractEnergyParameterDerivatives(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->energyParameterDerivatives;
}

/**
 * Make sure an expression doesn't use any undefined variables.
 */
static void validateVariables(const Lepton::ExpressionTreeNode& node, const set<string>& variables) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("Unknown variable in expression: "+op.getName());
    for (auto& child : node.getChildren())
        validateVariables(child, variables);
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.

    vector<Vec3> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (masses[i] > 0)
            shiftedVel[i] = velData[i]+forceData[i]*(timeShift/masses[i]);
        else
            shiftedVel[i] = velData[i];
    }
    
    // Apply constraints to them.

    if (timeShift != 0) {
        vector<double> inverseMasses(numParticles);
        for (int i = 0; i < numParticles; i++)
            inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
        extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    }
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}

/**
 * Copy particle charges into the fourth element of the posq array.
 */
static void copyChargesToPosq(ContextImpl& context, const vector<float>& charges, int index) {
    CpuPlatform::PlatformData& data = CpuPlatform::getPlatformData(context);
    if (index == data.currentPosqIndex)
        return;
    data.currentPosqIndex = index;
    AlignedArray<float>& posq = data.posq;
    for (int i = 0; i < charges.size(); i++)
        posq[4*i+3] = charges[i];
}

CpuCalcForcesAndEnergyKernel::CpuCalcForcesAndEnergyKernel(std::string name, const Platform& platform, CpuPlatform::PlatformData& data, ContextImpl& context) :
        CalcForcesAndEnergyKernel(name, platform), data(data) {
    // Create a Reference platform version of this kernel.
    
    ReferenceKernelFactory referenceFactory;
    referenceKernel = Kernel(referenceFactory.createKernelImpl(name, platform, context));
}

void CpuCalcForcesAndEnergyKernel::initialize(const System& system) {
    referenceKernel.getAs<ReferenceCalcForcesAndEnergyKernel>().initialize(system);
    lastPositions.resize(system.getNumParticles(), Vec3(1e10, 1e10, 1e10));
}

void CpuCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    referenceKernel.getAs<ReferenceCalcForcesAndEnergyKernel>().beginComputation(context, includeForce, includeEnergy, groups);
    
    // Convert positions to single precision and clear the forces.

    int numParticles = context.getSystem().getNumParticles();
    bool positionsValid = true;
    data.threads.execute([&] (ThreadPool& threads, int threadIndex) {
        // Convert the positions to single precision and apply periodic boundary conditions

        AlignedArray<float>& posq = data.posq;
        vector<Vec3>& posData = extractPositions(context);
        Vec3* boxVectors = extractBoxVectors(context);
        double boxSize[3] = {boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]};
        double invBoxSize[3] = {1/boxVectors[0][0], 1/boxVectors[1][1], 1/boxVectors[2][2]};
        bool triclinic = (boxVectors[0][1] != 0 || boxVectors[0][2] != 0 || boxVectors[1][0] != 0 || boxVectors[1][2] != 0 || boxVectors[2][0] != 0 || boxVectors[2][1] != 0);
        int numParticles = context.getSystem().getNumParticles();
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (data.isPeriodic) {
            if (triclinic) {
                for (int i = start; i < end; i++) {
                    Vec3 pos = posData[i];
                    pos -= boxVectors[2]*floor(pos[2]*invBoxSize[2]);
                    pos -= boxVectors[1]*floor(pos[1]*invBoxSize[1]);
                    pos -= boxVectors[0]*floor(pos[0]*invBoxSize[0]);
                    posq[4*i] = (float) pos[0];
                    posq[4*i+1] = (float) pos[1];
                    posq[4*i+2] = (float) pos[2];
                }
            }
            else {
                for (int i = start; i < end; i++) {
                    for (int j = 0; j < 3; j++) {
                        double x = posData[i][j];
                        double base = floor(x*invBoxSize[j])*boxSize[j];
                        posq[4*i+j] = (float) (x-base);
                    }
                }
            }
        }
        else
            for (int i = start; i < end; i++) {
                posq[4*i] = (float) posData[i][0];
                posq[4*i+1] = (float) posData[i][1];
                posq[4*i+2] = (float) posData[i][2];
            }
        
        // Check for invalid positions.
        
        for (int i = 4*start; i < 4*end; i += 4)
            if (posq[i] != posq[i] || posq[i+1] != posq[i+1] || posq[i+2] != posq[i+2])
                positionsValid = false;

        // Clear the forces.

        fvec4 zero(0.0f);
        for (int j = 0; j < numParticles; j++)
            zero.store(&data.threadForce[threadIndex][j*4]);
    });
    data.threads.waitForThreads();
    if (!positionsValid)
        throw OpenMMException("Particle coordinate is NaN.  For more information, see https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan");

    // Determine whether we need to recompute the neighbor list.
        
    if (data.neighborList != NULL && data.cutoff > 0.0) {
        double padding = data.paddedCutoff-data.cutoff;;
        bool needRecompute = false;
        double closeCutoff2 = 0.25*padding*padding;
        double farCutoff2 = 0.5*padding*padding;
        int maxNumMoved = numParticles/10;
        vector<int> moved;
        vector<Vec3>& posData = extractPositions(context);
        for (int i = 0; i < numParticles; i++) {
            Vec3 delta = posData[i]-lastPositions[i];
            double dist2 = delta.dot(delta);
            if (dist2 > closeCutoff2) {
                moved.push_back(i);
                if (dist2 > farCutoff2 || moved.size() > maxNumMoved) {
                    needRecompute = true;
                    break;
                }
            }
        }
        if (!needRecompute && moved.size() > 0) {
            // Some particles have moved further than half the padding distance.  Look for pairs
            // that are missing from the neighbor list.

            int numMoved = moved.size();
            double cutoff2 = data.cutoff*data.cutoff;
            double paddedCutoff2 = data.paddedCutoff*data.paddedCutoff;
            for (int i = 1; i < numMoved && !needRecompute; i++)
                for (int j = 0; j < i; j++) {
                    Vec3 delta = posData[moved[i]]-posData[moved[j]];
                    if (delta.dot(delta) < cutoff2) {
                        // These particles should interact.  See if they are in the neighbor list.
                        
                        Vec3 oldDelta = lastPositions[moved[i]]-lastPositions[moved[j]];
                        if (oldDelta.dot(oldDelta) > paddedCutoff2) {
                            needRecompute = true;
                            break;
                        }
                    }
                }
        }
        if (needRecompute) {
            data.neighborList->computeNeighborList(numParticles, data.posq, data.exclusions, extractBoxVectors(context), data.isPeriodic, data.paddedCutoff, data.threads);
            lastPositions = posData;
        }
    }
}

double CpuCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups, bool& valid) {
    // Sum the forces from all the threads.
    
    data.threads.execute([&] (ThreadPool& threads, int threadIndex) {
        // Sum the contributions to forces that have been calculated by different threads.
        
        int numParticles = context.getSystem().getNumParticles();
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        vector<Vec3>& forceData = extractForces(context);
        for (int i = start; i < end; i++) {
            fvec4 f(0.0f);
            for (int j = 0; j < numThreads; j++)
                f += fvec4(&data.threadForce[j][4*i]);
            forceData[i][0] += f[0];
            forceData[i][1] += f[1];
            forceData[i][2] += f[2];
        }
    });
    data.threads.waitForThreads();
    return referenceKernel.getAs<ReferenceCalcForcesAndEnergyKernel>().finishComputation(context, includeForce, includeEnergy, groups, valid);
}

void CpuCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    angleIndexArray.resize(numAngles, vector<int>(3));
    angleParamArray.resize(numAngles, vector<double>(2));
    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        angleParamArray[i][0] = angle;
        angleParamArray[i][1] = k;
    }
    bondForce.initialize(system.getNumParticles(), numAngles, 3, angleIndexArray, data.threads);
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double CpuCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceAngleBondIxn angleBond;
    if (usePeriodic)
        angleBond.setPeriodic(extractBoxVectors(context));
    bondForce.calculateForce(posData, angleParamArray, forceData, includeEnergy ? &energy : NULL, angleBond);
    return energy;
}

void CpuCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    for (int i = firstAngle; i <= lastAngle; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        if (particle1 != angleIndexArray[i][0] || particle2 != angleIndexArray[i][1] || particle3 != angleIndexArray[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        angleParamArray[i][0] = angle;
        angleParamArray[i][1] = k;
    }
}

void CpuCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray.resize(numTorsions, vector<int>(4));
    torsionParamArray.resize(numTorsions, vector<double>(3));
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = k;
        torsionParamArray[i][1] = phase;
        torsionParamArray[i][2] = periodicity;
    }
    bondForce.initialize(system.getNumParticles(), numTorsions, 4, torsionIndexArray, data.threads);
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double CpuCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceProperDihedralBond periodicTorsionBond;
    if (usePeriodic)
        periodicTorsionBond.setPeriodic(extractBoxVectors(context));
    bondForce.calculateForce(posData, torsionParamArray, forceData, includeEnergy ? &energy : NULL, periodicTorsionBond);
    return energy;
}

void CpuCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = firstTorsion; i <= lastTorsion; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        torsionParamArray[i][0] = k;
        torsionParamArray[i][1] = phase;
        torsionParamArray[i][2] = periodicity;
    }
}

void CpuCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray.resize(numTorsions, vector<int>(4));
    torsionParamArray.resize(numTorsions, vector<double>(6));
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = c0;
        torsionParamArray[i][1] = c1;
        torsionParamArray[i][2] = c2;
        torsionParamArray[i][3] = c3;
        torsionParamArray[i][4] = c4;
        torsionParamArray[i][5] = c5;
    }
    bondForce.initialize(system.getNumParticles(), numTorsions, 4, torsionIndexArray, data.threads);
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double CpuCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceRbDihedralBond rbTorsionBond;
    if (usePeriodic)
        rbTorsionBond.setPeriodic(extractBoxVectors(context));
    bondForce.calculateForce(posData, torsionParamArray, forceData, includeEnergy ? &energy : NULL, rbTorsionBond);
    return energy;
}

void CpuCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        torsionParamArray[i][0] = c0;
        torsionParamArray[i][1] = c1;
        torsionParamArray[i][2] = c2;
        torsionParamArray[i][3] = c3;
        torsionParamArray[i][4] = c4;
        torsionParamArray[i][5] = c5;
    }
}

class CpuCalcNonbondedForceKernel::PmeIO : public CalcPmeReciprocalForceKernel::IO {
public:
    PmeIO(float* posq, float* force, int numParticles) : posq(posq), force(force), numParticles(numParticles) {
    }
    float* getPosq() {
        return posq;
    }
    void setForce(float* f) {
        for (int i = 0; i < numParticles; i++) {
            force[4*i] += f[4*i];
            force[4*i+1] += f[4*i+1];
            force[4*i+2] += f[4*i+2];
        }
    }
private:
    float* posq;
    float* force;
    int numParticles;
};

CpuNonbondedForce* createCpuNonbondedForceVec(const CpuNeighborList& neighbors);

CpuCalcNonbondedForceKernel::CpuCalcNonbondedForceKernel(string name, const Platform& platform, CpuPlatform::PlatformData& data) : CalcNonbondedForceKernel(name, platform),
        data(data), hasInitializedPme(false), hasInitializedDispersionPme(false), nonbonded(NULL) {
}

CpuCalcNonbondedForceKernel::~CpuCalcNonbondedForceKernel() {
    if (nonbonded != NULL)
        delete nonbonded;
}

void CpuCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    chargePosqIndex = data.requestPosqIndex();
    ljPosqIndex = data.requestPosqIndex();

    // Identify which exceptions are 1-4 interactions.

    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    vector<int> nb14s;
    map<int, int> nb14Index;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end()) {
            nb14Index[i] = nb14s.size();
            nb14s.push_back(i);
        }
    }

    // Record the particle parameters.

    num14 = nb14s.size();
    bonded14IndexArray.resize(num14, vector<int>(2));
    bonded14ParamArray.resize(num14, vector<double>(3));
    particleParams.resize(numParticles);
    charges.resize(numParticles);
    C6params.resize(numParticles);
    baseParticleParams.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
       force.getParticleParameters(i, baseParticleParams[i][0], baseParticleParams[i][1], baseParticleParams[i][2]);
    
    // Record exception parameters.
    
    baseExceptionParams.resize(num14);
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        force.getExceptionParameters(nb14s[i], particle1, particle2, baseExceptionParams[i][0], baseExceptionParams[i][1], baseExceptionParams[i][2]);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
    }
    bondForce.initialize(system.getNumParticles(), num14, 2, bonded14IndexArray, data.threads);
    
    // Record information about parameter offsets.
    
    hasParticleOffsets = (force.getNumParticleParameterOffsets() > 0);
    hasExceptionOffsets = (force.getNumExceptionParameterOffsets() > 0);
    particleParamOffsets.resize(force.getNumParticles());
    exceptionParamOffsets.resize(force.getNumExceptions());
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string param;
        int particle;
        double charge, sigma, epsilon;
        force.getParticleParameterOffset(i, param, particle, charge, sigma, epsilon);
        auto paramPos = find(paramNames.begin(), paramNames.end(), param);
        int paramIndex;
        if (paramPos == paramNames.end()) {
            paramIndex = paramNames.size();
            paramNames.push_back(param);
        }
        else
            paramIndex = paramPos-paramNames.begin();
        particleParamOffsets[particle].push_back(make_tuple(charge, sigma, epsilon, paramIndex));
    }
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        auto paramPos = find(paramNames.begin(), paramNames.end(), param);
        int paramIndex;
        if (paramPos == paramNames.end()) {
            paramIndex = paramNames.size();
            paramNames.push_back(param);
        }
        else
            paramIndex = paramPos-paramNames.begin();
        exceptionParamOffsets[nb14Index[exception]].push_back(make_tuple(charge, sigma, epsilon, paramIndex));
    }
    paramValues.resize(paramNames.size(), 0.0);

    // Record other parameters.
    
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff) {
        data.requestNeighborList(0.0, 0.0, true, exclusions);
        useSwitchingFunction = false;
    }
    else {
        data.requestNeighborList(nonbondedCutoff, 0.25*nonbondedCutoff, true, exclusions);
        useSwitchingFunction = force.getUseSwitchingFunction();
        switchingDistance = force.getSwitchingDistance();
    }
    if (nonbondedMethod == Ewald) {
        double alpha;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmax[0], kmax[1], kmax[2]);
        ewaldAlpha = alpha;
    }
    else if (nonbondedMethod == PME) {
        double alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSize[0], gridSize[1], gridSize[2], false);
        ewaldAlpha = alpha;
    }
    else if (nonbondedMethod == LJPME) {
        double alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSize[0], gridSize[1], gridSize[2], false);
        ewaldAlpha = alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, dispersionGridSize[0], dispersionGridSize[1], dispersionGridSize[2], true);
        ewaldDispersionAlpha = alpha;
        useSwitchingFunction = false;
    }
    if (nonbondedMethod == NoCutoff || nonbondedMethod == CutoffNonPeriodic)
        exceptionsArePeriodic = false;
    else
        exceptionsArePeriodic = force.getExceptionsUsePeriodicBoundaryConditions();
    rfDielectric = force.getReactionFieldDielectric();
    if (force.getUseDispersionCorrection())
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
    data.isPeriodic |= (nonbondedMethod == CutoffPeriodic || nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME);
    nonbonded = createCpuNonbondedForceVec(*data.neighborList);
}

double CpuCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    if (!hasInitializedPme) {
        hasInitializedPme = true;
        useOptimizedPme = false;
        computeParameters(context, false);
        if (nonbondedMethod == PME) {
            // If available, use the optimized PME implementation.

            vector<string> kernelNames;
            kernelNames.push_back("CalcPmeReciprocalForce");
            useOptimizedPme = getPlatform().supportsKernels(kernelNames);
            if (useOptimizedPme) {
                optimizedPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), context);
                optimizedPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSize[0], gridSize[1], gridSize[2], numParticles, ewaldAlpha, data.deterministicForces);
            }
        }
        if (nonbondedMethod == LJPME) {
            // If available, use the optimized PME implementation.

            vector<string> kernelNames;
            kernelNames.push_back("CalcPmeReciprocalForce");
            kernelNames.push_back("CalcDispersionPmeReciprocalForce");
            useOptimizedPme = getPlatform().supportsKernels(kernelNames);
            if (useOptimizedPme) {
                optimizedPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), context);
                optimizedPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSize[0], gridSize[1], gridSize[2], numParticles, ewaldAlpha, data.deterministicForces);
                optimizedDispersionPme = getPlatform().createKernel(CalcDispersionPmeReciprocalForceKernel::Name(), context);
                optimizedDispersionPme.getAs<CalcDispersionPmeReciprocalForceKernel>().initialize(dispersionGridSize[0], dispersionGridSize[1],
                                                                                                  dispersionGridSize[2], numParticles, ewaldDispersionAlpha, data.deterministicForces);
            }
        }
    }
    computeParameters(context, true);
    copyChargesToPosq(context, charges, chargePosqIndex);
    AlignedArray<float>& posq = data.posq;
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    Vec3* boxVectors = extractBoxVectors(context);
    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    bool ljpme = (nonbondedMethod == LJPME);
    if (nonbondedMethod != NoCutoff)
        nonbonded->setUseCutoff(nonbondedCutoff, rfDielectric);
    if (data.isPeriodic) {
        Vec3* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 1.999999*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded->setPeriodic(boxVectors);
        nonbonded->setPeriodicExceptions(exceptionsArePeriodic);
    }
    if (ewald)
        nonbonded->setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        nonbonded->setUsePME(ewaldAlpha, gridSize);
    if (useSwitchingFunction)
        nonbonded->setUseSwitchingFunction(switchingDistance);
    if (ljpme){
        nonbonded->setUsePME(ewaldAlpha, gridSize);
        nonbonded->setUseLJPME(ewaldDispersionAlpha, dispersionGridSize);
    }
    double nonbondedEnergy = 0;
    if (includeDirect)
        nonbonded->calculateDirectIxn(numParticles, &posq[0], posData, particleParams, C6params, exclusions, data.threadForce, includeEnergy ? &nonbondedEnergy : NULL, data.threads);
    if (includeReciprocal) {
        if (useOptimizedPme) {
            PmeIO io(&posq[0], &data.threadForce[0][0], numParticles);
            Vec3 periodicBoxVectors[3] = {boxVectors[0], boxVectors[1], boxVectors[2]};
            optimizedPme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, periodicBoxVectors, includeEnergy);
            nonbondedEnergy += optimizedPme.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
            if (nonbondedMethod == LJPME) {
                copyChargesToPosq(context, C6params, ljPosqIndex);
                optimizedDispersionPme.getAs<CalcDispersionPmeReciprocalForceKernel>().beginComputation(io, periodicBoxVectors, includeEnergy);
                nonbondedEnergy += optimizedDispersionPme.getAs<CalcDispersionPmeReciprocalForceKernel>().finishComputation(io);
            }
        }
        else
            nonbonded->calculateReciprocalIxn(numParticles, &posq[0], posData, particleParams, C6params, exclusions, forceData, includeEnergy ? &nonbondedEnergy : NULL);
    }
    energy += nonbondedEnergy;
    if (includeDirect) {
        ReferenceLJCoulomb14 nonbonded14;
        if (exceptionsArePeriodic) {
            Vec3* boxVectors = extractBoxVectors(context);
            nonbonded14.setPeriodic(boxVectors);
        }
        bondForce.calculateForce(posData, bonded14ParamArray, forceData, includeEnergy ? &energy : NULL, nonbonded14);
        if (data.isPeriodic && nonbondedMethod != LJPME)
            energy += dispersionCoefficient/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]);
    }
    return energy;
}

void CpuCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException) {
    if (force.getNumParticles() != numParticles)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Identify which exceptions are 1-4 interactions.

    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end())
            nb14s.push_back(i);
    }
    if (nb14s.size() != num14)
        throw OpenMMException("updateParametersInContext: The number of non-excluded exceptions has changed");

    // Record the values.

    for (int i = firstParticle; i <= lastParticle; ++i)
       force.getParticleParameters(i, baseParticleParams[i][0], baseParticleParams[i][1], baseParticleParams[i][2]);
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        force.getExceptionParameters(nb14s[i], particle1, particle2, baseExceptionParams[i][0], baseExceptionParams[i][1], baseExceptionParams[i][2]);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
    }
    computeParameters(context, false);
    
    // Recompute the coefficient for the dispersion correction.

    NonbondedForce::NonbondedMethod method = force.getNonbondedMethod();
    if (force.getUseDispersionCorrection() && (method == NonbondedForce::CutoffPeriodic || method == NonbondedForce::Ewald || method == NonbondedForce::PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
}

void CpuCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME && nonbondedMethod != LJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (useOptimizedPme)
        optimizedPme.getAs<const CalcPmeReciprocalForceKernel>().getPMEParameters(alpha, nx, ny, nz);
    else {
        alpha = ewaldAlpha;
        nx = gridSize[0];
        ny = gridSize[1];
        nz = gridSize[2];
    }
}

void CpuCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != LJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (useOptimizedPme)
        optimizedDispersionPme.getAs<const CalcPmeReciprocalForceKernel>().getPMEParameters(alpha, nx, ny, nz);
    else {
        alpha = ewaldDispersionAlpha;
        nx = dispersionGridSize[0];
        ny = dispersionGridSize[1];
        nz = dispersionGridSize[2];
    }
}

void CpuCalcNonbondedForceKernel::computeParameters(ContextImpl& context, bool offsetsOnly) {
    bool paramChanged = false;
    for (int i = 0; i < paramNames.size(); i++) {
        double value = context.getParameter(paramNames[i]);
        if (value != paramValues[i]) {
            paramValues[i] = value;;
            paramChanged = true;
        }
    }
    if (!paramChanged && offsetsOnly)
        return;

    // Compute particle parameters.

    if (hasParticleOffsets || !offsetsOnly) {
        double sumSquaredCharges = 0.0;
        for (int i = 0; i < numParticles; i++) {
            double charge = baseParticleParams[i][0];
            double sigma = baseParticleParams[i][1];
            double epsilon = baseParticleParams[i][2];
            for (auto& offset : particleParamOffsets[i]) {
                double value = paramValues[get<3>(offset)];
                charge += value*get<0>(offset);
                sigma += value*get<1>(offset);
                epsilon += value*get<2>(offset);
            }
            charges[i] = (float) charge;
            particleParams[i] = make_pair((float) (0.5*sigma), (float) (2.0*sqrt(epsilon)));
            C6params[i] = 8.0*pow(particleParams[i].first, 3.0) * particleParams[i].second;
            sumSquaredCharges += charge*charge;
        }
        if (nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME) {
            ewaldSelfEnergy = -ONE_4PI_EPS0*ewaldAlpha*sumSquaredCharges/sqrt(M_PI);
            if (nonbondedMethod == LJPME)
                for (int atom = 0; atom < numParticles; atom++)
                    ewaldSelfEnergy += pow(ewaldDispersionAlpha, 6.0) * C6params[atom]*C6params[atom] / 12.0;
        }
        else
            ewaldSelfEnergy = 0.0;
        chargePosqIndex = data.requestPosqIndex();
        ljPosqIndex = data.requestPosqIndex();
    }

    // Compute exception parameters.

    if (hasExceptionOffsets || !offsetsOnly) {
        for (int i = 0; i < num14; i++) {
            double chargeProd = baseExceptionParams[i][0];
            double sigma = baseExceptionParams[i][1];
            double epsilon = baseExceptionParams[i][2];
            for (auto& offset : exceptionParamOffsets[i]) {
                double value = paramValues[get<3>(offset)];
                chargeProd += value*get<0>(offset);
                sigma += value*get<1>(offset);
                epsilon += value*get<2>(offset);
            }
            bonded14ParamArray[i][0] = sigma;
            bonded14ParamArray[i][1] = 4.0*epsilon;
            bonded14ParamArray[i][2] = chargeProd;
        }
    }
}

CpuCalcCustomNonbondedForceKernel::CpuCalcCustomNonbondedForceKernel(string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcCustomNonbondedForceKernel(name, platform), data(data), forceCopy(NULL), nonbonded(NULL) {
}

CpuCalcCustomNonbondedForceKernel::~CpuCalcCustomNonbondedForceKernel() {
    if (nonbonded != NULL)
        delete nonbonded;
    if (forceCopy != NULL)
        delete forceCopy;
}

void CpuCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {

    // Record the exclusions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }

    // Build the arrays.

    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        force.getParticleParameters(i, particleParamArray[i]);
    nonbondedMethod = CalcCustomNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff) {
        data.requestNeighborList(0.0, 0.0, true, exclusions);
        useSwitchingFunction = false;
    }
    else {
        data.requestNeighborList(nonbondedCutoff, 0.25*nonbondedCutoff, true, exclusions);
        useSwitchingFunction = force.getUseSwitchingFunction();
        switchingDistance = force.getSwitchingDistance();
    }

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Record information for the long range correction.

    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic && force.getUseLongRangeCorrection()) {
        forceCopy = new CustomNonbondedForce(force);
        hasInitializedLongRangeCorrection = false;
    }
    else {
        longRangeCoefficient = 0.0;
        hasInitializedLongRangeCorrection = true;
    }

    // Record the interaction groups.

    for (int i = 0; i < force.getNumInteractionGroups(); i++) {
        set<int> set1, set2;
        force.getInteractionGroupParameters(i, set1, set2);
        interactionGroups.push_back(make_pair(set1, set2));
    }
    data.isPeriodic |= (nonbondedMethod == CutoffPeriodic);

    // Create the interaction.

    createInteraction(force);
}

void CpuCalcCustomNonbondedForceKernel::createInteraction(const CustomNonbondedForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the various expressions used to calculate the force.

    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r");
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParameterNames.push_back(force.getGlobalParameterName(i));
        globalParamValues[force.getGlobalParameterName(i)] = force.getGlobalParameterDefaultValue(i);
    }
    set<string> particleVariables, pairVariables;
    particleVariables.insert("r");
    pairVariables.insert("r");
    for (auto& name : parameterNames) {
        particleVariables.insert(name);
        pairVariables.insert(name+"1");
        pairVariables.insert(name+"2");
    }
    particleVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    pairVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    vector<Lepton::ParsedExpression> computedValueExpressions, energyParamDerivExpressions;
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, exp;
        force.getComputedValueParameters(i, name, exp);
        Lepton::ParsedExpression parsed = Lepton::Parser::parse(exp, functions);
        validateVariables(parsed.getRootNode(), particleVariables);
        computedValueNames.push_back(name);
        computedValueExpressions.push_back(parsed);
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(energyExpression.differentiate(param));
    }
    for (auto& name : computedValueNames) {
        pairVariables.insert(name+"1");
        pairVariables.insert(name+"2");
    }
    validateVariables(energyExpression.getRootNode(), pairVariables);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;

    // Create the object that computes the interaction.

    nonbonded = createCpuCustomNonbondedForce(data.threads, *data.neighborList);
    nonbonded->initialize(energyExpression, forceExpression, parameterNames, exclusions, energyParamDerivExpressions,
            computedValueNames, computedValueExpressions);
    if (interactionGroups.size() > 0)
        nonbonded->setInteractionGroups(interactionGroups);
}

double CpuCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    Vec3* boxVectors = extractBoxVectors(context);
    double energy = 0;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff)
        nonbonded->setUseCutoff(nonbondedCutoff);
    if (periodic) {
        double minAllowedSize = 2*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded->setPeriodic(boxVectors);
    }
    bool globalParamsChanged = false;
    for (auto& name : globalParameterNames) {
        double value = context.getParameter(name);
        if (globalParamValues[name] != value)
            globalParamsChanged = true;
        globalParamValues[name] = value;
    }
    if (useSwitchingFunction)
        nonbonded->setUseSwitchingFunction(switchingDistance);
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    nonbonded->calculatePairIxn(numParticles, &data.posq[0], posData, particleParamArray, globalParamValues, data.threadForce, includeForces, includeEnergy, energy, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    
    // Add in the long range correction.
    
    if (!hasInitializedLongRangeCorrection) {
        longRangeCorrectionData = CustomNonbondedForceImpl::prepareLongRangeCorrection(*forceCopy, data.threads.getNumThreads());
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, data.threads);
        hasInitializedLongRangeCorrection = true;
    }
    else if (globalParamsChanged && forceCopy != NULL)
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, data.threads);
    double volume = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
    energy += longRangeCoefficient/volume;
    for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += longRangeCoefficientDerivs[i]/volume;
    return energy;
}

void CpuCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> parameters;
    for (int i = firstParticle; i <= lastParticle; ++i) {
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = parameters[j];
    }
    
    // If necessary, recompute the long range correction.
    
    if (forceCopy != NULL) {
        longRangeCorrectionData = CustomNonbondedForceImpl::prepareLongRangeCorrection(force, data.threads.getNumThreads());
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, data.threads);
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete nonbonded;
        nonbonded = NULL;
        createInteraction(force);
    }
}

CpuCalcGBSAOBCForceKernel::~CpuCalcGBSAOBCForceKernel() {
}

void CpuCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    posqIndex = data.requestPosqIndex();
    int numParticles = system.getNumParticles();
    particleParams.resize(numParticles);
    charges.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        charges[i] = (float) charge;
        radius -= 0.009;
        particleParams[i] = make_pair((float) radius, (float) (scalingFactor*radius));
    }
    obc.setParticleParameters(particleParams);
    obc.setSolventDielectric((float) force.getSolventDielectric());
    obc.setSoluteDielectric((float) force.getSoluteDielectric());
    obc.setSurfaceAreaEnergy((float) force.getSurfaceAreaEnergy());
    if (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff)
        obc.setUseCutoff((float) force.getCutoffDistance());
    data.isPeriodic |= (force.getNonbondedMethod() == GBSAOBCForce::CutoffPeriodic);
}

double CpuCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    copyChargesToPosq(context, charges, posqIndex);
    if (data.isPeriodic) {
        Vec3& boxSize = extractBoxSize(context);
        float floatBoxSize[3] = {(float) boxSize[0], (float) boxSize[1], (float) boxSize[2]};
        obc.setPeriodic(floatBoxSize);
    }
    double energy = 0.0;
    obc.computeForce(data.posq, data.threadForce, includeEnergy ? &energy : NULL, data.threads);
    return energy;
}

void CpuCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
    int numParticles = force.getNumParticles();
    if (numParticles != obc.getParticleParameters().size())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    posqIndex = data.requestPosqIndex();
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        charges[i] = (float) charge;
        radius -= 0.009;
        particleParams[i] = make_pair((float) radius, (float) (scalingFactor*radius));
    }
    obc.setParticleParameters(particleParams);
}

CpuCalcCustomGBForceKernel::~CpuCalcCustomGBForceKernel() {
    if (ixn != NULL)
        delete ixn;
    if (neighborList != NULL)
        delete neighborList;
}

void CpuCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    if (force.getNumComputedValues() > 0) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(0, name, expression, type);
        if (type == CustomGBForce::SingleParticle)
            throw OpenMMException("CpuPlatform requires that the first computed value for a CustomGBForce be of type ParticlePair or ParticlePairNoExclusions.");
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            force.getComputedValueParameters(i, name, expression, type);
            if (type != CustomGBForce::SingleParticle)
                throw OpenMMException("CpuPlatform requires that a CustomGBForce only have one computed value of type ParticlePair or ParticlePairNoExclusions.");
        }
    }

    // Record the exclusions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }

    // Build the arrays.

    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        force.getParticleParameters(i, particleParamArray[i]);
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        particleParameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    nonbondedMethod = CalcCustomGBForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod != NoCutoff)
        neighborList = new CpuNeighborList(4);
    data.isPeriodic |= (force.getNonbondedMethod() == CustomGBForce::CutoffPeriodic);

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the interaction.

    createInteraction(force);
}

void CpuCalcCustomGBForceKernel::createInteraction(const CustomGBForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expressions for computed values.

    valueTypes.clear();
    valueNames.clear();
    energyParamDerivNames.clear();
    vector<vector<Lepton::CompiledExpression> > valueDerivExpressions(force.getNumComputedValues());
    vector<vector<Lepton::CompiledExpression> > valueGradientExpressions(force.getNumComputedValues());
    vector<vector<Lepton::CompiledExpression> > valueParamDerivExpressions(force.getNumComputedValues());
    vector<Lepton::CompiledExpression> valueExpressions;
    vector<Lepton::CompiledExpression> energyExpressions;
    set<string> particleVariables, pairVariables;
    pairVariables.insert("r");
    particleVariables.insert("x");
    particleVariables.insert("y");
    particleVariables.insert("z");
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        particleVariables.insert(particleParameterNames[i]);
        pairVariables.insert(particleParameterNames[i]+"1");
        pairVariables.insert(particleParameterNames[i]+"2");
    }
    particleVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    pairVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(i, name, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        valueExpressions.push_back(ex.createCompiledExpression());
        valueTypes.push_back(type);
        valueNames.push_back(name);
        if (i == 0) {
            valueDerivExpressions[i].push_back(ex.differentiate("r").createCompiledExpression());
            validateVariables(ex.getRootNode(), pairVariables);
        }
        else {
            valueGradientExpressions[i].push_back(ex.differentiate("x").createCompiledExpression());
            valueGradientExpressions[i].push_back(ex.differentiate("y").createCompiledExpression());
            valueGradientExpressions[i].push_back(ex.differentiate("z").createCompiledExpression());
            for (int j = 0; j < i; j++)
                valueDerivExpressions[i].push_back(ex.differentiate(valueNames[j]).createCompiledExpression());
            validateVariables(ex.getRootNode(), particleVariables);
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++) {
            string param = force.getEnergyParameterDerivativeName(j);
            energyParamDerivNames.push_back(param);
            valueParamDerivExpressions[i].push_back(ex.differentiate(param).createCompiledExpression());
        }
        particleVariables.insert(name);
        pairVariables.insert(name+"1");
        pairVariables.insert(name+"2");
    }

    // Parse the expressions for energy terms.

    energyTypes.clear();
    vector<vector<Lepton::CompiledExpression> > energyDerivExpressions(force.getNumEnergyTerms());
    vector<vector<Lepton::CompiledExpression> > energyGradientExpressions(force.getNumEnergyTerms());
    vector<vector<Lepton::CompiledExpression> > energyParamDerivExpressions(force.getNumEnergyTerms());
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        energyExpressions.push_back(ex.createCompiledExpression());
        energyTypes.push_back(type);
        if (type != CustomGBForce::SingleParticle)
            energyDerivExpressions[i].push_back(ex.differentiate("r").createCompiledExpression());
        for (int j = 0; j < force.getNumComputedValues(); j++) {
            if (type == CustomGBForce::SingleParticle) {
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]).createCompiledExpression());
                energyGradientExpressions[i].push_back(ex.differentiate("x").createCompiledExpression());
                energyGradientExpressions[i].push_back(ex.differentiate("y").createCompiledExpression());
                energyGradientExpressions[i].push_back(ex.differentiate("z").createCompiledExpression());
                validateVariables(ex.getRootNode(), particleVariables);
            }
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"1").createCompiledExpression());
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"2").createCompiledExpression());
                validateVariables(ex.getRootNode(), pairVariables);
            }
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
            energyParamDerivExpressions[i].push_back(ex.differentiate(force.getEnergyParameterDerivativeName(j)).createCompiledExpression());
    }

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
    ixn = new CpuCustomGBForce(numParticles, exclusions, valueExpressions, valueDerivExpressions, valueGradientExpressions, valueParamDerivExpressions,
        valueNames, valueTypes, energyExpressions, energyDerivExpressions, energyGradientExpressions, energyParamDerivExpressions, energyTypes,
        particleParameterNames, data.threads);
}

double CpuCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    Vec3* boxVectors = extractBoxVectors(context);
    if (data.isPeriodic)
        ixn->setPeriodic(extractBoxSize(context));
    if (nonbondedMethod != NoCutoff) {
        vector<set<int> > noExclusions(numParticles);
        neighborList->computeNeighborList(numParticles, data.posq, noExclusions, boxVectors, data.isPeriodic, nonbondedCutoff, data.threads);
        ixn->setUseCutoff(nonbondedCutoff, *neighborList);
    }
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    ixn->calculateIxn(numParticles, &data.posq[0], particleParamArray, globalParameters, data.threadForce, includeForces, includeEnergy, energy, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void CpuCalcCustomGBForceKernel::copyParametersToContext(ContextImpl& context, const CustomGBForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<double>(parameters[j]);
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete ixn;
        ixn = NULL;
        createInteraction(force);
    }
}

CpuCalcCustomManyParticleForceKernel::~CpuCalcCustomManyParticleForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void CpuCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {

    // Build the arrays.

    numParticles = system.getNumParticles();
    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        int type;
        force.getParticleParameters(i, particleParamArray[i], type);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the interaction.

    ixn = new CpuCustomManyParticleForce(force, data.threads);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    cutoffDistance = force.getCutoffDistance();
    data.isPeriodic |= (nonbondedMethod == CutoffPeriodic);
}

double CpuCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    if (nonbondedMethod == CutoffPeriodic) {
        Vec3* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 2*cutoffDistance;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        ixn->setPeriodic(boxVectors);
    }
    double energy = 0;
    ixn->calculateIxn(data.posq, particleParamArray, globalParameters, data.threadForce, includeForces, includeEnergy, energy);
    return energy;
}

void CpuCalcCustomManyParticleForceKernel::copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        int type;
        force.getParticleParameters(i, parameters, type);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<double>(parameters[j]);
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete ixn;
        ixn = NULL;
        ixn = new CpuCustomManyParticleForce(force, data.threads);
    }
}

CpuCalcGayBerneForceKernel::~CpuCalcGayBerneForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void CpuCalcGayBerneForceKernel::initialize(const System& system, const GayBerneForce& force) {
    ixn = new CpuGayBerneForce(force);
    data.isPeriodic |= (force.getNonbondedMethod() == GayBerneForce::CutoffPeriodic);
    if (force.getNonbondedMethod() != GayBerneForce::NoCutoff) {
        double cutoff = force.getCutoffDistance();
        data.requestNeighborList(cutoff, 0.1*cutoff, true, ixn->getExclusions());
    }
}

double CpuCalcGayBerneForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return ixn->calculateForce(extractPositions(context), extractForces(context), data.threadForce, extractBoxVectors(context), data);
}

void CpuCalcGayBerneForceKernel::copyParametersToContext(ContextImpl& context, const GayBerneForce& force) {
    delete ixn;
    ixn = NULL;
    ixn = new CpuGayBerneForce(force);
}

CpuIntegrateLangevinMiddleStepKernel::~CpuIntegrateLangevinMiddleStepKernel() {
    if (dynamics)
        delete dynamics;
}

void CpuIntegrateLangevinMiddleStepKernel::initialize(const System& system, const LangevinMiddleIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    data.random.initialize(integrator.getRandomNumberSeed(), data.threads.getNumThreads());
}

void CpuIntegrateLangevinMiddleStepKernel::execute(ContextImpl& context, const LangevinMiddleIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new CpuLangevinMiddleDynamics(context.getSystem().getNumParticles(), stepSize, friction, temperature, data.threads, data.random);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context, posData, velData, masses, integrator.getConstraintTolerance());
    ReferencePlatform::PlatformData* refData = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    refData->time += stepSize;
    refData->stepCount++;
}

double CpuIntegrateLangevinMiddleStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinMiddleIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.0);
}
