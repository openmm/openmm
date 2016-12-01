/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2016 Stanford University and the Authors.      *
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
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/internal/vectorize.h"
#include "RealVec.h"
#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"

using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(RealVec*) data->periodicBoxSize;
}

static RealVec* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealVec*) data->periodicBoxVectors;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
}

static map<string, double>& extractEnergyParameterDerivatives(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((map<string, double>*) data->energyParameterDerivatives);
}

/**
 * Make sure an expression doesn't use any undefined variables.
 */
static void validateVariables(const Lepton::ExpressionTreeNode& node, const set<string>& variables) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("Unknown variable in expression: "+op.getName());
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        validateVariables(node.getChildren()[i], variables);
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.
    
    vector<RealVec> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (masses[i] > 0)
            shiftedVel[i] = velData[i]+forceData[i]*(timeShift/masses[i]);
        else
            shiftedVel[i] = velData[i];
    }
    
    // Apply constraints to them.
    
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; i++)
        inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
    extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}

class CpuCalcForcesAndEnergyKernel::SumForceTask : public ThreadPool::Task {
public:
    SumForceTask(int numParticles, vector<RealVec>& forceData, CpuPlatform::PlatformData& data) : numParticles(numParticles), forceData(forceData), data(data) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        // Sum the contributions to forces that have been calculated by different threads.
        
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        for (int i = start; i < end; i++) {
            fvec4 f(0.0f);
            for (int j = 0; j < numThreads; j++)
                f += fvec4(&data.threadForce[j][4*i]);
            forceData[i][0] += f[0];
            forceData[i][1] += f[1];
            forceData[i][2] += f[2];
        }
    }
    int numParticles;
    vector<RealVec>& forceData;
    CpuPlatform::PlatformData& data;
};

class CpuCalcForcesAndEnergyKernel::InitForceTask : public ThreadPool::Task {
public:
    InitForceTask(int numParticles, ContextImpl& context, CpuPlatform::PlatformData& data) : numParticles(numParticles), positionsValid(true), context(context), data(data) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        // Convert the positions to single precision and apply periodic boundary conditions

        AlignedArray<float>& posq = data.posq;
        vector<RealVec>& posData = extractPositions(context);
        RealVec* boxVectors = extractBoxVectors(context);
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
                    RealVec pos = posData[i];
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
                        RealOpenMM x = posData[i][j];
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
    }
    int numParticles;
    bool positionsValid;
    ContextImpl& context;
    CpuPlatform::PlatformData& data;
};

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
    InitForceTask task(numParticles, context, data);
    data.threads.execute(task);
    data.threads.waitForThreads();
    if (!task.positionsValid)
        throw OpenMMException("Particle coordinate is nan");

    // Determine whether we need to recompute the neighbor list.
        
    if (data.neighborList != NULL) {
        double padding = data.paddedCutoff-data.cutoff;;
        bool needRecompute = false;
        double closeCutoff2 = 0.25*padding*padding;
        double farCutoff2 = 0.5*padding*padding;
        int maxNumMoved = numParticles/10;
        vector<int> moved;
        vector<RealVec>& posData = extractPositions(context);
        for (int i = 0; i < numParticles; i++) {
            RealVec delta = posData[i]-lastPositions[i];
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
                    RealVec delta = posData[moved[i]]-posData[moved[j]];
                    if (delta.dot(delta) < cutoff2) {
                        // These particles should interact.  See if they are in the neighbor list.
                        
                        RealVec oldDelta = lastPositions[moved[i]]-lastPositions[moved[j]];
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
    
    SumForceTask task(context.getSystem().getNumParticles(), extractForces(context), data);
    data.threads.execute(task);
    data.threads.waitForThreads();
    return referenceKernel.getAs<ReferenceCalcForcesAndEnergyKernel>().finishComputation(context, includeForce, includeEnergy, groups, valid);
}

CpuCalcHarmonicAngleForceKernel::~CpuCalcHarmonicAngleForceKernel() {
    if (angleIndexArray != NULL) {
        for (int i = 0; i < numAngles; i++) {
            delete[] angleIndexArray[i];
            delete[] angleParamArray[i];
        }
        delete[] angleIndexArray;
        delete[] angleParamArray;
    }
}

void CpuCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    angleIndexArray = new int*[numAngles];
    for (int i = 0; i < numAngles; i++)
        angleIndexArray[i] = new int[3];
    angleParamArray = new RealOpenMM*[numAngles];
    for (int i = 0; i < numAngles; i++)
        angleParamArray[i] = new RealOpenMM[2];
    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        angleParamArray[i][0] = (RealOpenMM) angle;
        angleParamArray[i][1] = (RealOpenMM) k;
    }
    bondForce.initialize(system.getNumParticles(), numAngles, 3, angleIndexArray, data.threads);
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double CpuCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceAngleBondIxn angleBond;
    if (usePeriodic)
        angleBond.setPeriodic(extractBoxVectors(context));
    bondForce.calculateForce(posData, angleParamArray, forceData, includeEnergy ? &energy : NULL, angleBond);
    return energy;
}

void CpuCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        if (particle1 != angleIndexArray[i][0] || particle2 != angleIndexArray[i][1] || particle3 != angleIndexArray[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        angleParamArray[i][0] = (RealOpenMM) angle;
        angleParamArray[i][1] = (RealOpenMM) k;
    }
}

CpuCalcPeriodicTorsionForceKernel::~CpuCalcPeriodicTorsionForceKernel() {
    if (torsionIndexArray != NULL) {
        for (int i = 0; i < numTorsions; i++) {
            delete[] torsionIndexArray[i];
            delete[] torsionParamArray[i];
        }
        delete[] torsionIndexArray;
        delete[] torsionParamArray;
    }
}

void CpuCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray = new int*[numTorsions];
    for (int i = 0; i < numTorsions; i++)
        torsionIndexArray[i] = new int[4];
    torsionParamArray = new RealOpenMM*[numTorsions];
    for (int i = 0; i < numTorsions; i++)
        torsionParamArray[i] = new RealOpenMM[3];
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = (RealOpenMM) k;
        torsionParamArray[i][1] = (RealOpenMM) phase;
        torsionParamArray[i][2] = (RealOpenMM) periodicity;
    }
    bondForce.initialize(system.getNumParticles(), numTorsions, 4, torsionIndexArray, data.threads);
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double CpuCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceProperDihedralBond periodicTorsionBond;
    if (usePeriodic)
        periodicTorsionBond.setPeriodic(extractBoxVectors(context));
    bondForce.calculateForce(posData, torsionParamArray, forceData, includeEnergy ? &energy : NULL, periodicTorsionBond);
    return energy;
}

void CpuCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        torsionParamArray[i][0] = (RealOpenMM) k;
        torsionParamArray[i][1] = (RealOpenMM) phase;
        torsionParamArray[i][2] = (RealOpenMM) periodicity;
    }
}

CpuCalcRBTorsionForceKernel::~CpuCalcRBTorsionForceKernel() {
    if (torsionIndexArray != NULL) {
        for (int i = 0; i < numTorsions; i++) {
            delete[] torsionIndexArray[i];
            delete[] torsionParamArray[i];
        }
        delete[] torsionIndexArray;
        delete[] torsionParamArray;
    }
}

void CpuCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray = new int*[numTorsions];
    for (int i = 0; i < numTorsions; i++)
        torsionIndexArray[i] = new int[4];
    torsionParamArray = new RealOpenMM*[numTorsions];
    for (int i = 0; i < numTorsions; i++)
        torsionParamArray[i] = new RealOpenMM[6];
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = (RealOpenMM) c0;
        torsionParamArray[i][1] = (RealOpenMM) c1;
        torsionParamArray[i][2] = (RealOpenMM) c2;
        torsionParamArray[i][3] = (RealOpenMM) c3;
        torsionParamArray[i][4] = (RealOpenMM) c4;
        torsionParamArray[i][5] = (RealOpenMM) c5;
    }
    bondForce.initialize(system.getNumParticles(), numTorsions, 4, torsionIndexArray, data.threads);
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double CpuCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
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
        torsionParamArray[i][0] = (RealOpenMM) c0;
        torsionParamArray[i][1] = (RealOpenMM) c1;
        torsionParamArray[i][2] = (RealOpenMM) c2;
        torsionParamArray[i][3] = (RealOpenMM) c3;
        torsionParamArray[i][4] = (RealOpenMM) c4;
        torsionParamArray[i][5] = (RealOpenMM) c5;
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

bool isVec8Supported();
CpuNonbondedForce* createCpuNonbondedForceVec4();
CpuNonbondedForce* createCpuNonbondedForceVec8();

CpuCalcNonbondedForceKernel::CpuCalcNonbondedForceKernel(string name, const Platform& platform, CpuPlatform::PlatformData& data) : CalcNonbondedForceKernel(name, platform),
        data(data), bonded14IndexArray(NULL), bonded14ParamArray(NULL), hasInitializedPme(false), nonbonded(NULL) {
    if (isVec8Supported())
        nonbonded = createCpuNonbondedForceVec8();
    else
        nonbonded = createCpuNonbondedForceVec4();
}

CpuCalcNonbondedForceKernel::~CpuCalcNonbondedForceKernel() {
    if (bonded14ParamArray != NULL) {
        for (int i = 0; i < num14; i++) {
            delete[] bonded14IndexArray[i];
            delete[] bonded14ParamArray[i];
        }
        delete[] bonded14IndexArray;
        delete[] bonded14ParamArray;
    }
    if (nonbonded != NULL)
        delete nonbonded;
}

void CpuCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {

    // Identify which exceptions are 1-4 interactions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
        if (chargeProd != 0.0 || epsilon != 0.0)
            nb14s.push_back(i);
    }

    // Record the particle parameters.

    num14 = nb14s.size();
    bonded14IndexArray = new int*[num14];
    for (int i = 0; i < num14; i++)
        bonded14IndexArray[i] = new int[2];
    bonded14ParamArray = new double*[num14];
    for (int i = 0; i < num14; i++)
        bonded14ParamArray[i] = new double[3];
    particleParams.resize(numParticles);
    double sumSquaredCharges = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, depth;
        force.getParticleParameters(i, charge, radius, depth);
        data.posq[4*i+3] = (float) charge;
        particleParams[i] = make_pair((float) (0.5*radius), (float) (2.0*sqrt(depth)));
        sumSquaredCharges += charge*charge;
    }
    
    // Recorded exception parameters.
    
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        double charge, radius, depth;
        force.getExceptionParameters(nb14s[i], particle1, particle2, charge, radius, depth);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
        bonded14ParamArray[i][0] = static_cast<RealOpenMM>(radius);
        bonded14ParamArray[i][1] = static_cast<RealOpenMM>(4.0*depth);
        bonded14ParamArray[i][2] = static_cast<RealOpenMM>(charge);
    }
    bondForce.initialize(system.getNumParticles(), num14, 2, bonded14IndexArray, data.threads);
    
    // Record other parameters.
    
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff)
        useSwitchingFunction = false;
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
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSize[0], gridSize[1], gridSize[2]);
        ewaldAlpha = alpha;
    }
    if (nonbondedMethod == Ewald || nonbondedMethod == PME)
        ewaldSelfEnergy = -ONE_4PI_EPS0*ewaldAlpha*sumSquaredCharges/sqrt(M_PI);
    else
        ewaldSelfEnergy = 0.0;
    rfDielectric = force.getReactionFieldDielectric();
    if (force.getUseDispersionCorrection())
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
    data.isPeriodic = (nonbondedMethod == CutoffPeriodic || nonbondedMethod == Ewald || nonbondedMethod == PME);
}

double CpuCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    if (!hasInitializedPme) {
        hasInitializedPme = true;
        useOptimizedPme = false;
        if (nonbondedMethod == PME) {
            // If available, use the optimized PME implementation.

            vector<string> kernelNames;
            kernelNames.push_back("CalcPmeReciprocalForce");
            useOptimizedPme = getPlatform().supportsKernels(kernelNames);
            if (useOptimizedPme) {
                optimizedPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), context);
                optimizedPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSize[0], gridSize[1], gridSize[2], numParticles, ewaldAlpha);
            }
        }
    }
    AlignedArray<float>& posq = data.posq;
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealVec* boxVectors = extractBoxVectors(context);
    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    if (nonbondedMethod != NoCutoff)
        nonbonded->setUseCutoff(nonbondedCutoff, *data.neighborList, rfDielectric);
    if (data.isPeriodic) {
        RealVec* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 1.999999*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded->setPeriodic(boxVectors);
    }
    if (ewald)
        nonbonded->setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        nonbonded->setUsePME(ewaldAlpha, gridSize);
    if (useSwitchingFunction)
        nonbonded->setUseSwitchingFunction(switchingDistance);
    double nonbondedEnergy = 0;
    if (includeDirect)
        nonbonded->calculateDirectIxn(numParticles, &posq[0], posData, particleParams, exclusions, data.threadForce, includeEnergy ? &nonbondedEnergy : NULL, data.threads);
    if (includeReciprocal) {
        if (useOptimizedPme) {
            PmeIO io(&posq[0], &data.threadForce[0][0], numParticles);
            Vec3 periodicBoxVectors[3] = {boxVectors[0], boxVectors[1], boxVectors[2]};
            optimizedPme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, periodicBoxVectors, includeEnergy);
            nonbondedEnergy += optimizedPme.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
        }
        else
            nonbonded->calculateReciprocalIxn(numParticles, &posq[0], posData, particleParams, exclusions, forceData, includeEnergy ? &nonbondedEnergy : NULL);
    }
    energy += nonbondedEnergy;
    if (includeDirect) {
        ReferenceLJCoulomb14 nonbonded14;
        bondForce.calculateForce(posData, bonded14ParamArray, forceData, includeEnergy ? &energy : NULL, nonbonded14);
        if (data.isPeriodic)
            energy += dispersionCoefficient/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]);
    }
    return energy;
}

void CpuCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
    if (force.getNumParticles() != numParticles)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (chargeProd != 0.0 || epsilon != 0.0)
            nb14s.push_back(i);
    }
    if (nb14s.size() != num14)
        throw OpenMMException("updateParametersInContext: The number of non-excluded exceptions has changed");

    // Record the values.

    double sumSquaredCharges = 0.0;
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, depth;
        force.getParticleParameters(i, charge, radius, depth);
        data.posq[4*i+3] = (float) charge;
        particleParams[i] = make_pair((float) (0.5*radius), (float) (2.0*sqrt(depth)));
        sumSquaredCharges += charge*charge;
    }
    if (nonbondedMethod == Ewald || nonbondedMethod == PME)
        ewaldSelfEnergy = -ONE_4PI_EPS0*ewaldAlpha*sumSquaredCharges/sqrt(M_PI);
    else
        ewaldSelfEnergy = 0.0;
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        double charge, radius, depth;
        force.getExceptionParameters(nb14s[i], particle1, particle2, charge, radius, depth);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
        bonded14ParamArray[i][0] = static_cast<RealOpenMM>(radius);
        bonded14ParamArray[i][1] = static_cast<RealOpenMM>(4.0*depth);
        bonded14ParamArray[i][2] = static_cast<RealOpenMM>(charge);
    }
    
    // Recompute the coefficient for the dispersion correction.

    NonbondedForce::NonbondedMethod method = force.getNonbondedMethod();
    if (force.getUseDispersionCorrection() && (method == NonbondedForce::CutoffPeriodic || method == NonbondedForce::Ewald || method == NonbondedForce::PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
}

void CpuCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME)
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

CpuCalcCustomNonbondedForceKernel::CpuCalcCustomNonbondedForceKernel(string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcCustomNonbondedForceKernel(name, platform), data(data), forceCopy(NULL), nonbonded(NULL) {
}

CpuCalcCustomNonbondedForceKernel::~CpuCalcCustomNonbondedForceKernel() {
    if (particleParamArray != NULL) {
        for (int i = 0; i < numParticles; i++)
            delete[] particleParamArray[i];
        delete[] particleParamArray;
    }
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

    int numParameters = force.getNumPerParticleParameters();
    particleParamArray = new double*[numParticles];
    for (int i = 0; i < numParticles; i++)
        particleParamArray[i] = new double[numParameters];
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = parameters[j];
    }
    nonbondedMethod = CalcCustomNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff)
        useSwitchingFunction = false;
    else {
        data.requestNeighborList(nonbondedCutoff, 0.25*nonbondedCutoff, true, exclusions);
        useSwitchingFunction = force.getUseSwitchingFunction();
        switchingDistance = force.getSwitchingDistance();
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the various expressions used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    Lepton::CompiledExpression energyExpression = expression.createCompiledExpression();
    Lepton::CompiledExpression forceExpression = expression.differentiate("r").createCompiledExpression();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParameterNames.push_back(force.getGlobalParameterName(i));
        globalParamValues[force.getGlobalParameterName(i)] = force.getGlobalParameterDefaultValue(i);
    }
    std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(expression.differentiate(param).createCompiledExpression());
    }
    set<string> variables;
    variables.insert("r");
    for (int i = 0; i < numParameters; i++) {
        variables.insert(parameterNames[i]+"1");
        variables.insert(parameterNames[i]+"2");
    }
    variables.insert(globalParameterNames.begin(), globalParameterNames.end());
    validateVariables(expression.getRootNode(), variables);

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
    
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
    data.isPeriodic = (nonbondedMethod == CutoffPeriodic);
    nonbonded = new CpuCustomNonbondedForce(energyExpression, forceExpression, parameterNames, exclusions, energyParamDerivExpressions, data.threads);
    if (interactionGroups.size() > 0)
        nonbonded->setInteractionGroups(interactionGroups);
}

double CpuCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealVec* boxVectors = extractBoxVectors(context);
    double energy = 0;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff)
        nonbonded->setUseCutoff(nonbondedCutoff, *data.neighborList);
    if (periodic) {
        double minAllowedSize = 2*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded->setPeriodic(boxVectors);
    }
    bool globalParamsChanged = false;
    for (int i = 0; i < (int) globalParameterNames.size(); i++) {
        double value = context.getParameter(globalParameterNames[i]);
        if (globalParamValues[globalParameterNames[i]] != value)
            globalParamsChanged = true;
        globalParamValues[globalParameterNames[i]] = value;
    }
    if (useSwitchingFunction)
        nonbonded->setUseSwitchingFunction(switchingDistance);
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    nonbonded->calculatePairIxn(numParticles, &data.posq[0], posData, particleParamArray, 0, globalParamValues, data.threadForce, includeForces, includeEnergy, energy, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    
    // Add in the long range correction.
    
    if (!hasInitializedLongRangeCorrection || (globalParamsChanged && forceCopy != NULL)) {
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
    }
    double volume = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
    energy += longRangeCoefficient/volume;
    for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += longRangeCoefficientDerivs[i]/volume;
    return energy;
}

void CpuCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = parameters[j];
    }
    
    // If necessary, recompute the long range correction.
    
    if (forceCopy != NULL) {
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }
}

CpuCalcGBSAOBCForceKernel::~CpuCalcGBSAOBCForceKernel() {
}

void CpuCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    int numParticles = system.getNumParticles();
    particleParams.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        data.posq[4*i+3] = (float) charge;
        radius -= 0.009;
        particleParams[i] = make_pair((float) radius, (float) (scalingFactor*radius));
    }
    obc.setParticleParameters(particleParams);
    obc.setSolventDielectric((float) force.getSolventDielectric());
    obc.setSoluteDielectric((float) force.getSoluteDielectric());
    obc.setSurfaceAreaEnergy((float) force.getSurfaceAreaEnergy());
    if (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff)
        obc.setUseCutoff((float) force.getCutoffDistance());
    data.isPeriodic = (force.getNonbondedMethod() == GBSAOBCForce::CutoffPeriodic);
}

double CpuCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (data.isPeriodic) {
        RealVec& boxSize = extractBoxSize(context);
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

    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        data.posq[4*i+3] = (float) charge;
        radius -= 0.009;
        particleParams[i] = make_pair((float) radius, (float) (scalingFactor*radius));
    }
    obc.setParticleParameters(particleParams);
}

CpuCalcCustomGBForceKernel::~CpuCalcCustomGBForceKernel() {
    if (particleParamArray != NULL) {
        for (int i = 0; i < numParticles; i++)
            delete[] particleParamArray[i];
        delete[] particleParamArray;
    }
    if (ixn != NULL)
        delete ixn;
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

    int numPerParticleParameters = force.getNumPerParticleParameters();
    particleParamArray = new double*[numParticles];
    for (int i = 0; i < numParticles; i++)
        particleParamArray[i] = new double[numPerParticleParameters];
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numPerParticleParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    for (int i = 0; i < numPerParticleParameters; i++)
        particleParameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    nonbondedMethod = CalcCustomGBForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    if (nonbondedMethod != NoCutoff)
        data.requestNeighborList(nonbondedCutoff, 0.25*nonbondedCutoff, force.getNumExclusions() > 0, exclusions);

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expressions for computed values.

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
    for (int i = 0; i < numPerParticleParameters; i++) {
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

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
    ixn = new CpuCustomGBForce(numParticles, exclusions, valueExpressions, valueDerivExpressions, valueGradientExpressions, valueParamDerivExpressions,
        valueNames, valueTypes, energyExpressions, energyDerivExpressions, energyGradientExpressions, energyParamDerivExpressions, energyTypes,
        particleParameterNames, data.threads);
    data.isPeriodic = (force.getNonbondedMethod() == CustomGBForce::CutoffPeriodic);
}

double CpuCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    RealVec* boxVectors = extractBoxVectors(context);
    if (data.isPeriodic)
        ixn->setPeriodic(extractBoxSize(context));
    if (nonbondedMethod != NoCutoff) {
        vector<set<int> > noExclusions(numParticles);
        ixn->setUseCutoff(nonbondedCutoff, *data.neighborList);
    }
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
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
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
}

CpuCalcCustomManyParticleForceKernel::~CpuCalcCustomManyParticleForceKernel() {
    if (particleParamArray != NULL) {
        for (int i = 0; i < numParticles; i++)
            delete[] particleParamArray[i];
        delete[] particleParamArray;
    }
    if (ixn != NULL)
        delete ixn;
}

void CpuCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {

    // Build the arrays.

    numParticles = system.getNumParticles();
    int numParticleParameters = force.getNumPerParticleParameters();
    particleParamArray = new double*[numParticles];
    for (int i = 0; i < numParticles; i++)
        particleParamArray[i] = new double[numParticleParameters];
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        int type;
        force.getParticleParameters(i, parameters, type);
        for (int j = 0; j < numParticleParameters; j++)
            particleParamArray[i][j] = parameters[j];
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    ixn = new CpuCustomManyParticleForce(force, data.threads);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    cutoffDistance = force.getCutoffDistance();
    data.isPeriodic = (nonbondedMethod == CutoffPeriodic);
}

double CpuCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    if (nonbondedMethod == CutoffPeriodic) {
        RealVec* boxVectors = extractBoxVectors(context);
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
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
}

CpuCalcGayBerneForceKernel::~CpuCalcGayBerneForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void CpuCalcGayBerneForceKernel::initialize(const System& system, const GayBerneForce& force) {
    ixn = new CpuGayBerneForce(force);
    data.isPeriodic = (force.getNonbondedMethod() == GayBerneForce::CutoffPeriodic);
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

CpuIntegrateLangevinStepKernel::~CpuIntegrateLangevinStepKernel() {
    if (dynamics)
        delete dynamics;
}

void CpuIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    data.random.initialize(integrator.getRandomNumberSeed(), data.threads.getNumThreads());
}

void CpuIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new CpuLangevinDynamics(context.getSystem().getNumParticles(), stepSize, friction, temperature, data.threads, data.random);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    ReferencePlatform::PlatformData* refData = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    refData->time += stepSize;
    refData->stepCount++;
}

double CpuIntegrateLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}
