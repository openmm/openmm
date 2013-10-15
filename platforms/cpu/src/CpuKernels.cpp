/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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
#include "ReferenceLJCoulomb14.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "RealVec.h"

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

CpuCalcNonbondedForceKernel::~CpuCalcNonbondedForceKernel() {
    if (bonded14ParamArray != NULL) {
        for (int i = 0; i < num14; i++) {
            delete[] bonded14IndexArray[i];
            delete[] bonded14ParamArray[i];
        }
        delete bonded14IndexArray;
        delete bonded14ParamArray;
    }
}

void CpuCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {

    // Identify which exceptions are 1-4 interactions.

    numParticles = force.getNumParticles();
    posq.resize(4*numParticles, 0);
    forces.resize(4*numParticles, 0);
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
        posq[4*i+3] = (float) charge;
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
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealVec boxSize = extractBoxSize(context);
    float floatBoxSize[3] = {(float) boxSize[0], (float) boxSize[1], (float) boxSize[2]};
    double energy = ewaldSelfEnergy;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    
    // Convert the positions to single precision.
    
    if (periodic)
        for (int i = 0; i < numParticles; i++)
            for (int j = 0; j < 3; j++) {
                RealOpenMM x = posData[i][j];
                double base = floor(x/boxSize[j]+0.5)*boxSize[j];
                posq[4*i+j] = (float) (x-base);
            }
    else
        for (int i = 0; i < numParticles; i++) {
            posq[4*i] = (float) posData[i][0];
            posq[4*i+1] = (float) posData[i][1];
            posq[4*i+2] = (float) posData[i][2];
        }
    for (int i = 0; i < 4*numParticles; i++)
        forces[i] = 0.0f;
    if (nonbondedMethod != NoCutoff) {
        // Determine whether we need to recompute the neighbor list.
        
        double padding = 0.1*nonbondedCutoff;
        bool needRecompute = false;
        for (int i = 0; i < numParticles; i++) {
            RealVec delta = posData[i]-lastPositions[i];
            if (delta.dot(delta) > 0.25*padding*padding) {
                needRecompute = true;
                break;
            }
        }
        if (needRecompute) {
            neighborList.computeNeighborList(numParticles, posq, exclusions, floatBoxSize, periodic || ewald || pme, nonbondedCutoff+padding);
            lastPositions = posData;
        }
        nonbonded.setUseCutoff(nonbondedCutoff, neighborList.getNeighbors(), rfDielectric);
    }
    if (periodic || ewald || pme) {
        double minAllowedSize = 1.999999*nonbondedCutoff;
        if (boxSize[0] < minAllowedSize || boxSize[1] < minAllowedSize || boxSize[2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        nonbonded.setPeriodic(floatBoxSize);
    }
    if (ewald)
        nonbonded.setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        nonbonded.setUsePME(ewaldAlpha, gridSize);
    if (useSwitchingFunction)
        nonbonded.setUseSwitchingFunction(switchingDistance);
    float nonbondedEnergy = 0;
    if (includeDirect)
        nonbonded.calculateDirectIxn(numParticles, &posq[0], particleParams, exclusions, &forces[0], includeEnergy ? &nonbondedEnergy : NULL);
    if (includeReciprocal) {
        if (useOptimizedPme) {
            PmeIO io(&posq[0], &forces[0], numParticles);
            Vec3 periodicBoxSize(boxSize[0], boxSize[1], boxSize[2]);
            optimizedPme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, periodicBoxSize, includeEnergy);
            optimizedPme.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
        }
        else
            nonbonded.calculateReciprocalIxn(numParticles, &posq[0], posData, particleParams, exclusions, forceData, includeEnergy ? &nonbondedEnergy : NULL);
    }
    energy += nonbondedEnergy;
    for (int i = 0; i < numParticles; i++) {
        forceData[i][0] += forces[4*i];
        forceData[i][1] += forces[4*i+1];
        forceData[i][2] += forces[4*i+2];
    }
    if (includeDirect) {
        ReferenceBondForce refBondForce;
        ReferenceLJCoulomb14 nonbonded14;
        refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, includeEnergy ? &energy : NULL, nonbonded14);
        if (periodic || ewald || pme)
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
        posq[4*i+3] = (float) charge;
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
