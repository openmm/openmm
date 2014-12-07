/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2014 Stanford University and the Authors.      *
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

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
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
    InitForceTask(int numParticles, ContextImpl& context, CpuPlatform::PlatformData& data) : numParticles(numParticles), context(context), data(data) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        // Convert the positions to single precision and apply periodic boundary conditions

        AlignedArray<float>& posq = data.posq;
        vector<RealVec>& posData = extractPositions(context);
        RealVec boxSize = extractBoxSize(context);
        double invBoxSize[3] = {1/boxSize[0], 1/boxSize[1], 1/boxSize[2]};
        int numParticles = context.getSystem().getNumParticles();
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (data.isPeriodic)
            for (int i = start; i < end; i++)
                for (int j = 0; j < 3; j++) {
                    RealOpenMM x = posData[i][j];
                    double base = floor(x*invBoxSize[j])*boxSize[j];
                    posq[4*i+j] = (float) (x-base);
                }
        else
            for (int i = start; i < end; i++) {
                posq[4*i] = (float) posData[i][0];
                posq[4*i+1] = (float) posData[i][1];
                posq[4*i+2] = (float) posData[i][2];
            }

        // Clear the forces.

        fvec4 zero(0.0f);
        for (int j = 0; j < numParticles; j++)
            zero.store(&data.threadForce[threadIndex][j*4]);
    }
    int numParticles;
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
}

void CpuCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    referenceKernel.getAs<ReferenceCalcForcesAndEnergyKernel>().beginComputation(context, includeForce, includeEnergy, groups);
    
    // Convert positions to single precision and clear the forces.
    
    InitForceTask task(context.getSystem().getNumParticles(), context, data);
    data.threads.execute(task);
    data.threads.waitForThreads();
}

double CpuCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForce, bool includeEnergy, int groups) {
    // Sum the forces from all the threads.
    
    SumForceTask task(context.getSystem().getNumParticles(), extractForces(context), data);
    data.threads.execute(task);
    data.threads.waitForThreads();
    return referenceKernel.getAs<ReferenceCalcForcesAndEnergyKernel>().finishComputation(context, includeForce, includeEnergy, groups);
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
}

double CpuCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceProperDihedralBond periodicTorsionBond;
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
}

double CpuCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceRbDihedralBond rbTorsionBond;
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
        data(data), bonded14IndexArray(NULL), bonded14ParamArray(NULL), hasInitializedPme(false), neighborList(NULL), nonbonded(NULL) {
    if (isVec8Supported()) {
        neighborList = new CpuNeighborList(8);
        nonbonded = createCpuNonbondedForceVec8();
    }
    else {
        neighborList = new CpuNeighborList(4);
        nonbonded = createCpuNonbondedForceVec4();
    }
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
    if (neighborList != NULL)
        delete neighborList;
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
    
    // Record other parameters.
    
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff)
        useSwitchingFunction = false;
    else {
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
    lastPositions.resize(numParticles, Vec3(1e10, 1e10, 1e10));
    data.isPeriodic = (nonbondedMethod == CutoffPeriodic || nonbondedMethod == Ewald || nonbondedMethod == PME);
}

double CpuCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    if (!hasInitializedPme) {
        hasInitializedPme = true;
        useOptimizedPme = false;
        if (nonbondedMethod == PME) {
            // If available, use the optimized PME implementation.

            try {
                optimizedPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), context);
                optimizedPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSize[0], gridSize[1], gridSize[2], numParticles, ewaldAlpha);
                useOptimizedPme = true;
            }
            catch (OpenMMException& ex) {
                // The CPU PME plugin isn't available.
            }
        }
    }
    AlignedArray<float>& posq = data.posq;
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealVec boxSize = extractBoxSize(context);
    float floatBoxSize[3] = {(float) boxSize[0], (float) boxSize[1], (float) boxSize[2]};
    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    if (nonbondedMethod != NoCutoff) {
        // Determine whether we need to recompute the neighbor list.
        
        double padding = 0.15*nonbondedCutoff;
        bool needRecompute = false;
        double closeCutoff2 = 0.25*padding*padding;
        double farCutoff2 = 0.5*padding*padding;
        int maxNumMoved = numParticles/10;
        vector<int> moved;
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
            double cutoff2 = nonbondedCutoff*nonbondedCutoff;
            double paddedCutoff2 = (nonbondedCutoff+padding)*(nonbondedCutoff+padding);
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
            neighborList->computeNeighborList(numParticles, posq, exclusions, floatBoxSize, data.isPeriodic, nonbondedCutoff+padding, data.threads);
            lastPositions = posData;
        }
        nonbonded->setUseCutoff(nonbondedCutoff, *neighborList, rfDielectric);
    }
    if (data.isPeriodic) {
        double minAllowedSize = 1.999999*nonbondedCutoff;
        if (boxSize[0] < minAllowedSize || boxSize[1] < minAllowedSize || boxSize[2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded->setPeriodic(floatBoxSize);
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
            Vec3 periodicBoxSize(boxSize[0], boxSize[1], boxSize[2]);
            optimizedPme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, periodicBoxSize, includeEnergy);
            nonbondedEnergy += optimizedPme.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
        }
        else
            nonbonded->calculateReciprocalIxn(numParticles, &posq[0], posData, particleParams, exclusions, forceData, includeEnergy ? &nonbondedEnergy : NULL);
    }
    energy += nonbondedEnergy;
    if (includeDirect) {
        ReferenceBondForce refBondForce;
        ReferenceLJCoulomb14 nonbonded14;
        refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, includeEnergy ? &energy : NULL, nonbonded14);
        if (data.isPeriodic)
            energy += dispersionCoefficient/(boxSize[0]*boxSize[1]*boxSize[2]);
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

CpuCalcCustomNonbondedForceKernel::CpuCalcCustomNonbondedForceKernel(string name, const Platform& platform, CpuPlatform::PlatformData& data) :
            CalcCustomNonbondedForceKernel(name, platform), data(data), forceCopy(NULL), neighborList(NULL), nonbonded(NULL) {
}

CpuCalcCustomNonbondedForceKernel::~CpuCalcCustomNonbondedForceKernel() {
    if (particleParamArray != NULL) {
        for (int i = 0; i < numParticles; i++)
            delete[] particleParamArray[i];
        delete[] particleParamArray;
    }
    if (neighborList != NULL)
        delete neighborList;
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
        neighborList = new CpuNeighborList(4);
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
    nonbonded = new CpuCustomNonbondedForce(energyExpression, forceExpression, parameterNames, exclusions, data.threads);
    if (interactionGroups.size() > 0)
        nonbonded->setInteractionGroups(interactionGroups);
}

double CpuCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealVec& box = extractBoxSize(context);
    float floatBoxSize[3] = {(float) box[0], (float) box[1], (float) box[2]};
    double energy = 0;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        neighborList->computeNeighborList(numParticles, data.posq, exclusions, floatBoxSize, data.isPeriodic, nonbondedCutoff, data.threads);
        nonbonded->setUseCutoff(nonbondedCutoff, *neighborList);
    }
    if (periodic) {
        double minAllowedSize = 2*nonbondedCutoff;
        if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded->setPeriodic(box);
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
    nonbonded->calculatePairIxn(numParticles, &data.posq[0], posData, particleParamArray, 0, globalParamValues, data.threadForce, includeForces, includeEnergy, energy);
    
    // Add in the long range correction.
    
    if (!hasInitializedLongRangeCorrection || (globalParamsChanged && forceCopy != NULL)) {
        longRangeCoefficient = CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner());
        hasInitializedLongRangeCorrection = true;
    }
    energy += longRangeCoefficient/(box[0]*box[1]*box[2]);
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
        longRangeCoefficient = CustomNonbondedForceImpl::calcLongRangeCorrection(force, context.getOwner());
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
    if (neighborList != NULL)
        delete neighborList;
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
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new CpuNeighborList(4);

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expressions for computed values.

    vector<vector<Lepton::CompiledExpression> > valueDerivExpressions(force.getNumComputedValues());
    vector<vector<Lepton::CompiledExpression> > valueGradientExpressions(force.getNumComputedValues());
    vector<Lepton::CompiledExpression> valueExpressions;
    vector<Lepton::CompiledExpression> energyExpressions;
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(i, name, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        valueExpressions.push_back(ex.createCompiledExpression());
        valueTypes.push_back(type);
        valueNames.push_back(name);
        if (i == 0)
            valueDerivExpressions[i].push_back(ex.differentiate("r").createCompiledExpression());
        else {
            valueGradientExpressions[i].push_back(ex.differentiate("x").createCompiledExpression());
            valueGradientExpressions[i].push_back(ex.differentiate("y").createCompiledExpression());
            valueGradientExpressions[i].push_back(ex.differentiate("z").createCompiledExpression());
            for (int j = 0; j < i; j++)
                valueDerivExpressions[i].push_back(ex.differentiate(valueNames[j]).createCompiledExpression());
        }
    }

    // Parse the expressions for energy terms.

    vector<vector<Lepton::CompiledExpression> > energyDerivExpressions(force.getNumEnergyTerms());
    vector<vector<Lepton::CompiledExpression> > energyGradientExpressions(force.getNumEnergyTerms());
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
            }
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"1").createCompiledExpression());
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"2").createCompiledExpression());
            }
        }
    }

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
    ixn = new CpuCustomGBForce(numParticles, exclusions, valueExpressions, valueDerivExpressions, valueGradientExpressions, valueNames, valueTypes, energyExpressions,
        energyDerivExpressions, energyGradientExpressions, energyTypes, particleParameterNames, data.threads);
    data.isPeriodic = (force.getNonbondedMethod() == CustomGBForce::CutoffPeriodic);
}

double CpuCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    RealVec& box = extractBoxSize(context);
    float floatBoxSize[3] = {(float) box[0], (float) box[1], (float) box[2]};
    if (data.isPeriodic)
        ixn->setPeriodic(extractBoxSize(context));
    if (nonbondedMethod != NoCutoff) {
        vector<set<int> > noExclusions(numParticles);
        neighborList->computeNeighborList(numParticles, data.posq, exclusions, floatBoxSize, data.isPeriodic, nonbondedCutoff, data.threads);
        ixn->setUseCutoff(nonbondedCutoff, *neighborList);
    }
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ixn->calculateIxn(numParticles, &data.posq[0], particleParamArray, globalParameters, data.threadForce, includeForces, includeEnergy, energy);
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
}

double CpuCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    if (nonbondedMethod == CutoffPeriodic) {
        RealVec& box = extractBoxSize(context);
        double minAllowedSize = 2*cutoffDistance;
        if (box[0] < minAllowedSize || box[1] < minAllowedSize || box[2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        ixn->setPeriodic(box);
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
        RealOpenMM tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        dynamics = new CpuLangevinDynamics(context.getSystem().getNumParticles(), stepSize, tau, temperature, data.threads, data.random);
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
