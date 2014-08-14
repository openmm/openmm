
/* Portions copyright (c) 2009-2014 Stanford University and Simbios.
 * Contributors: Peter Eastman
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <string.h>
#include <sstream>
#include <utility>

#include "SimTKOpenMMCommon.h"
#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "CpuCustomManyParticleForce.h"
#include "ReferenceTabulatedFunction.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "lepton/CustomFunction.h"
#include "gmx_atomic.h"

using namespace OpenMM;
using namespace std;

class CpuCustomManyParticleForce::ComputeForceTask : public ThreadPool::Task {
public:
    ComputeForceTask(CpuCustomManyParticleForce& owner) : owner(owner) {
    }
    void execute(ThreadPool& threads, int threadIndex) {
        owner.threadComputeForce(threads, threadIndex);
    }
    CpuCustomManyParticleForce& owner;
};

CpuCustomManyParticleForce::CpuCustomManyParticleForce(const CustomManyParticleForce& force, ThreadPool& threads) :
            threads(threads), useCutoff(false), usePeriodic(false), neighborList(NULL) {
    numParticlesPerSet = force.getNumParticlesPerSet();
    numPerParticleParameters = force.getNumPerParticleParameters();
    
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the objects used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpr = CustomManyParticleForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    for (int i = 0; i < threads.getNumThreads(); i++)
        threadData.push_back(new ThreadData(force, energyExpr, distances, angles, dihedrals));
    if (force.getNonbondedMethod() != CustomManyParticleForce::NoCutoff)
        setUseCutoff(force.getCutoffDistance());

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
    
    // Record exclusions.
    
    exclusions.resize(force.getNumParticles());
    for (int i = 0; i < (int) force.getNumExclusions(); i++) {
        int p1, p2;
        force.getExclusionParticles(i, p1, p2);
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
    }
    
    // Record information about type filters.
    
    CustomManyParticleForceImpl::buildFilterArrays(force, numTypes, particleTypes, orderIndex, particleOrder);
}

CpuCustomManyParticleForce::~CpuCustomManyParticleForce() {
    if (neighborList != NULL)
        delete neighborList;
    for (int i = 0; i < (int) threadData.size(); i++)
        delete threadData[i];
}

void CpuCustomManyParticleForce::calculateIxn(AlignedArray<float>& posq, vector<RealVec>& atomCoordinates, RealOpenMM** particleParameters,
                                                  const map<string, double>& globalParameters, vector<AlignedArray<float> >& threadForce,
                                                  bool includeForces, bool includeEnergy, double& energy) {
    // Record the parameters for the threads.
    
    this->numParticles = atomCoordinates.size();
    this->posq = &posq[0];
    this->atomCoordinates = &atomCoordinates[0];
    this->particleParameters = particleParameters;
    this->globalParameters = &globalParameters;
    this->threadForce = &threadForce;
    this->includeForces = includeForces;
    this->includeEnergy = includeEnergy;
    gmx_atomic_t counter;
    gmx_atomic_set(&counter, 0);
    this->atomicCounter = &counter;
    if (useCutoff) {
        // Construct a neighbor list.
        
        float boxSizeFloat[] = {(float) periodicBoxSize[0], (float) periodicBoxSize[1], (float) periodicBoxSize[2]};
        neighborList->computeNeighborList(numParticles, posq, exclusions, boxSizeFloat, usePeriodic, cutoffDistance, threads);
    }
    
    // Signal the threads to start running and wait for them to finish.
    
    ComputeForceTask task(*this);
    threads.execute(task);
    threads.waitForThreads();
    
    // Combine the energies from all the threads.
    
    if (includeEnergy) {
        int numThreads = threads.getNumThreads();
        for (int i = 0; i < numThreads; i++)
            energy += threadData[i]->energy;
    }
}

void CpuCustomManyParticleForce::threadComputeForce(ThreadPool& threads, int threadIndex) {
    vector<int> particleIndices(numParticlesPerSet);
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    float* forces = &(*threadForce)[threadIndex][0];
    ThreadData& data = *threadData[threadIndex];
    data.energy = 0;
    for (map<string, double>::const_iterator iter = globalParameters->begin(); iter != globalParameters->end(); ++iter)
        data.expressionSet.setVariable(data.expressionSet.getVariableIndex(iter->first), iter->second);
    if (useCutoff) {
        // Loop over interactions from the neighbor list.
        
        while (true) {
            int blockIndex = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (blockIndex >= neighborList->getNumBlocks())
                break;
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
            int numNeighbors = neighbors.size();
            for (int i = 0; i < 4; i++) {
                particleIndices[0] = neighborList->getSortedAtoms()[4*blockIndex+i];
                
                // Build a filtered list of neighbors after removing exclusions.  We'll check for actual exclusions
                // again later, but the neighbor list also includes padding atoms that it marks as exclusions, so
                // we need to remove those now.
                
                vector<int> particles;
                for (int j = 0; j < numNeighbors; j++)
                    if ((exclusions[j] & (1<<i)) == 0)
                        particles.push_back(neighbors[j]);
                loopOverInteractions(particles, particleIndices, 1, 0, particleParameters, forces, data, boxSize, invBoxSize);
            }
        }
    }
    else {
        // Loop over all possible sets of particles.
        
        vector<int> particles(numParticles);
        for (int i = 0; i < numParticles; i++)
            particles[i] = i;
        while (true) {
            int i = gmx_atomic_fetch_add(reinterpret_cast<gmx_atomic_t*>(atomicCounter), 1);
            if (i >= numParticles)
                break;
            particleIndices[0] = i;
            loopOverInteractions(particles, particleIndices, 1, i+1, particleParameters, forces, data, boxSize, invBoxSize);
        }
    }
}

void CpuCustomManyParticleForce::setUseCutoff(RealOpenMM distance) {
    useCutoff = true;
    cutoffDistance = distance;
    if (neighborList == NULL)
        neighborList = new CpuNeighborList(4);
}

void CpuCustomManyParticleForce::setPeriodic(RealVec& boxSize) {
    assert(useCutoff);
    assert(boxSize[0] >= 2.0*cutoffDistance);
    assert(boxSize[1] >= 2.0*cutoffDistance);
    assert(boxSize[2] >= 2.0*cutoffDistance);
    usePeriodic = true;
    periodicBoxSize[0] = boxSize[0];
    periodicBoxSize[1] = boxSize[1];
    periodicBoxSize[2] = boxSize[2];
}

void CpuCustomManyParticleForce::loopOverInteractions(vector<int>& availableParticles, vector<int>& particleSet, int loopIndex, int startIndex,
                                                          RealOpenMM** particleParameters, float* forces, ThreadData& data, const fvec4& boxSize, const fvec4& invBoxSize) {
    int numParticles = availableParticles.size();
    double cutoff2 = cutoffDistance*cutoffDistance;
    for (int i = startIndex; i < numParticles; i++) {
        int particle = availableParticles[i];
        
        // Check whether this particle can actually participate in interactions with the others found so far.
        
        bool include = true;
        if (useCutoff) {
            fvec4 deltaR;
            fvec4 pos1(posq+4*particle);
            float r2;
            for (int j = 0; j < loopIndex && include; j++) {
                fvec4 pos2(posq+4*particleSet[j]);
                getDeltaR(pos1, pos2, deltaR, r2, boxSize, invBoxSize);
                include &= (r2 < cutoff2);
            }
        }
        for (int j = 0; j < loopIndex && include; j++)
            include &= (exclusions[particle].find(particleSet[j]) == exclusions[particle].end());
        if (include) {
            particleSet[loopIndex] = availableParticles[i];
            if (loopIndex == numParticlesPerSet-1)
                calculateOneIxn(particleSet, particleParameters, forces, data, boxSize, invBoxSize);
            else
                loopOverInteractions(availableParticles, particleSet, loopIndex+1, i+1, particleParameters, forces, data, boxSize, invBoxSize);
        }
    }
}

void CpuCustomManyParticleForce::calculateOneIxn(vector<int>& particleSet, RealOpenMM** particleParameters, float* forces, ThreadData& data, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Select the ordering to use for the particles.
    
    vector<int> permutedParticles(numParticlesPerSet);
    if (particleOrder.size() == 1) {
        // There are no filters, so we don't need to worry about ordering.
        
        permutedParticles = particleSet;
    }
    else {
        int index = 0;
        for (int i = numParticlesPerSet-1; i >= 0; i--)
            index = particleTypes[particleSet[i]]+numTypes*index;
        int order = orderIndex[index];
        if (order == -1)
            return;
        for (int i = 0; i < numParticlesPerSet; i++)
            permutedParticles[i] = particleSet[particleOrder[order][i]];
    }

    // Record per-particle parameters.
    
    CompiledExpressionSet& expressionSet = data.expressionSet;
    for (int i = 0; i < numParticlesPerSet; i++)
        for (int j = 0; j < numPerParticleParameters; j++)
            expressionSet.setVariable(data.particleParamIndices[i][j], particleParameters[permutedParticles[i]][j]);
    
    // Compute inter-particle deltas.
    
    int numDeltas = data.deltaPairs.size();
    RealOpenMM delta[numDeltas][ReferenceForce::LastDeltaRIndex];
    for (int i = 0; i < numDeltas; i++)
        computeDelta(permutedParticles[data.deltaPairs[i].first], permutedParticles[data.deltaPairs[i].second], delta[i], atomCoordinates);
    
    // Compute all of the variables the energy can depend on.

    for (int i = 0; i < (int) data.particleTerms.size(); i++) {
        const ParticleTermInfo& term = data.particleTerms[i];
        expressionSet.setVariable(term.variableIndex, atomCoordinates[permutedParticles[term.atom]][term.component]);
    }
    for (int i = 0; i < (int) data.distanceTerms.size(); i++) {
        const DistanceTermInfo& term = data.distanceTerms[i];
        expressionSet.setVariable(term.variableIndex, delta[term.delta][ReferenceForce::RIndex]);
    }
    for (int i = 0; i < (int) data.angleTerms.size(); i++) {
        const AngleTermInfo& term = data.angleTerms[i];
        expressionSet.setVariable(term.variableIndex, computeAngle(delta[term.delta1], delta[term.delta2], term.delta1Sign*term.delta2Sign));
    }
    for (int i = 0; i < (int) data.dihedralTerms.size(); i++) {
        const DihedralTermInfo& term = data.dihedralTerms[i];
        RealOpenMM dotDihedral, signOfDihedral;
        RealOpenMM* crossProduct[] = {term.cross1, term.cross2};
        expressionSet.setVariable(term.variableIndex, ReferenceBondIxn::getDihedralAngleBetweenThreeVectors(delta[term.delta1], delta[term.delta2], delta[term.delta3], crossProduct, &dotDihedral, delta[term.delta1], &signOfDihedral, 1));
    }
    
    if (includeForces) {
        // Apply forces based on individual particle coordinates.

        AlignedArray<fvec4> f(numParticlesPerSet);
        for (int i = 0; i < numParticlesPerSet; i++)
            f[i] = fvec4(0.0f);
        for (int i = 0; i < (int) data.particleTerms.size(); i++) {
            const ParticleTermInfo& term = data.particleTerms[i];
            float temp[4];
            f[term.atom].store(temp);
            temp[term.component] -= term.forceExpression.evaluate();
            f[term.atom] = fvec4(temp);
        }

        // Apply forces based on distances.

        for (int i = 0; i < (int) data.distanceTerms.size(); i++) {
            const DistanceTermInfo& term = data.distanceTerms[i];
            RealOpenMM dEdR = (RealOpenMM) (term.forceExpression.evaluate()*term.deltaSign/(delta[term.delta][ReferenceForce::RIndex]));
            fvec4 force = -dEdR*fvec4((float) delta[term.delta][0], (float) delta[term.delta][1], (float) delta[term.delta][2], 0.0f);
            f[term.p1] -= force;
            f[term.p2] += force;
        }

        // Apply forces based on angles.

        for (int i = 0; i < (int) data.angleTerms.size(); i++) {
            const AngleTermInfo& term = data.angleTerms[i];
            RealOpenMM dEdTheta = (RealOpenMM) term.forceExpression.evaluate();
            fvec4 delta1((float) delta[term.delta1][0], (float) delta[term.delta1][1], (float) delta[term.delta1][2], 0.0f);
            fvec4 delta2((float) delta[term.delta2][0], (float) delta[term.delta2][1], (float) delta[term.delta2][2], 0.0f);
            fvec4 thetaCross = cross(delta1, delta2);
            float lengthThetaCross = sqrtf(dot3(thetaCross, thetaCross));
            if (lengthThetaCross < 1.0e-6f)
                lengthThetaCross = 1.0e-6f;
            RealOpenMM termA = dEdTheta*term.delta2Sign/(delta[term.delta1][ReferenceForce::R2Index]*lengthThetaCross);
            RealOpenMM termC = -dEdTheta*term.delta1Sign/(delta[term.delta2][ReferenceForce::R2Index]*lengthThetaCross);
            fvec4 deltaCross1 = cross(delta1, thetaCross);
            fvec4 deltaCross2 = cross(delta2, thetaCross);
            fvec4 force1 = termA*deltaCross1;
            fvec4 force3 = termC*deltaCross2;
            fvec4 force2 = -(force1+force3);
            f[term.p1] += force1;
            f[term.p2] += force2;
            f[term.p3] += force3;
        }

        // Apply forces based on dihedrals.

        for (int i = 0; i < (int) data.dihedralTerms.size(); i++) {
            const DihedralTermInfo& term = data.dihedralTerms[i];
            RealOpenMM dEdTheta = (RealOpenMM) term.forceExpression.evaluate();
            RealOpenMM internalF[4][3];
            RealOpenMM forceFactors[4];
            RealOpenMM normCross1 = DOT3(term.cross1, term.cross1);
            RealOpenMM normBC = delta[term.delta2][ReferenceForce::RIndex];
            forceFactors[0] = (-dEdTheta*normBC)/normCross1;
            RealOpenMM normCross2 = DOT3(term.cross2, term.cross2);
            forceFactors[3] = (dEdTheta*normBC)/normCross2;
            forceFactors[1] = DOT3(delta[term.delta1], delta[term.delta2]);
            forceFactors[1] /= delta[term.delta2][ReferenceForce::R2Index];
            forceFactors[2] = DOT3(delta[term.delta3], delta[term.delta2]);
            forceFactors[2] /= delta[term.delta2][ReferenceForce::R2Index];
            fvec4 force1 = forceFactors[0]*fvec4((float) term.cross1[0], (float) term.cross1[1], (float) term.cross1[2], 0.0f);
            fvec4 force4 = forceFactors[3]*fvec4((float) term.cross2[0], (float) term.cross2[1], (float) term.cross2[2], 0.0f);
            fvec4 s = forceFactors[1]*force1 - forceFactors[2]*force4;
            f[term.p1] += force1;
            f[term.p2] -= force1-s;
            f[term.p3] -= force4+s;
            f[term.p4] += force4;
        }

        // Store the forces.

        for (int i = 0; i < numParticlesPerSet; i++) {
            int index = permutedParticles[i];
            (fvec4(forces+4*index)+f[i]).store(forces+4*index);
        }
    }

    // Add the energy

    if (includeEnergy)
        data.energy += data.energyExpression.evaluate();
}

void CpuCustomManyParticleForce::computeDelta(int atom1, int atom2, RealOpenMM* delta, const OpenMM::RealVec* atomCoordinates) const {
    if (usePeriodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom1], atomCoordinates[atom2], periodicBoxSize, delta);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom1], atomCoordinates[atom2], delta);
}

RealOpenMM CpuCustomManyParticleForce::computeAngle(RealOpenMM* vec1, RealOpenMM* vec2, float sign) {
    RealOpenMM dot = DOT3(vec1, vec2)*sign;
    RealOpenMM cosine = dot/SQRT((vec1[ReferenceForce::R2Index]*vec2[ReferenceForce::R2Index]));
    RealOpenMM angle;
    if (cosine >= 1)
        angle = 0;
    else if (cosine <= -1)
        angle = PI_M;
    else
        angle = ACOS(cosine);
    return angle;
}

void CpuCustomManyParticleForce::getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (usePeriodic) {
        fvec4 base = round(deltaR*invBoxSize)*boxSize;
        deltaR = deltaR-base;
    }
    r2 = dot3(deltaR, deltaR);
}

CpuCustomManyParticleForce::ParticleTermInfo::ParticleTermInfo(const string& name, int atom, int component, const Lepton::CompiledExpression& forceExpression, ThreadData& data) :
        name(name), atom(atom), component(component), forceExpression(forceExpression) {
    variableIndex = data.expressionSet.getVariableIndex(name);
}

CpuCustomManyParticleForce::DistanceTermInfo::DistanceTermInfo(const string& name, const vector<int>& atoms, const Lepton::CompiledExpression& forceExpression, ThreadData& data) :
        name(name), p1(atoms[0]), p2(atoms[1]), forceExpression(forceExpression) {
    variableIndex = data.expressionSet.getVariableIndex(name);
    data.requestDeltaPair(p1, p2, delta, deltaSign, true);
}

CpuCustomManyParticleForce::AngleTermInfo::AngleTermInfo(const string& name, const vector<int>& atoms, const Lepton::CompiledExpression& forceExpression, ThreadData& data) :
        name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), forceExpression(forceExpression) {
    variableIndex = data.expressionSet.getVariableIndex(name);
    data.requestDeltaPair(p1, p2,delta1, delta1Sign, true);
    data.requestDeltaPair(p3, p2, delta2, delta2Sign, true);
}

CpuCustomManyParticleForce::DihedralTermInfo::DihedralTermInfo(const string& name, const vector<int>& atoms, const Lepton::CompiledExpression& forceExpression, ThreadData& data) :
        name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), p4(atoms[3]), forceExpression(forceExpression) {
    variableIndex = data.expressionSet.getVariableIndex(name);
    float sign;
    data.requestDeltaPair(p2, p1, delta1, sign, false);
    data.requestDeltaPair(p2, p3, delta2, sign, false);
    data.requestDeltaPair(p4, p3, delta3, sign, false);
}

CpuCustomManyParticleForce::ThreadData::ThreadData(const CustomManyParticleForce& force, Lepton::ParsedExpression& energyExpr,
            map<string, vector<int> >& distances, map<string, vector<int> >& angles, map<string, vector<int> >& dihedrals) {
    int numParticlesPerSet = force.getNumParticlesPerSet();
    int numPerParticleParameters = force.getNumPerParticleParameters();
    particleParamIndices.resize(numParticlesPerSet);
    energyExpression = energyExpr.createCompiledExpression();
    expressionSet.registerExpression(energyExpression);

    // Differentiate the energy to get expressions for the force.

    for (int i = 0; i < numParticlesPerSet; i++) {
        stringstream xname, yname, zname;
        xname << 'x' << (i+1);
        yname << 'y' << (i+1);
        zname << 'z' << (i+1);
        particleTerms.push_back(CpuCustomManyParticleForce::ParticleTermInfo(xname.str(), i, 0, energyExpr.differentiate(xname.str()).optimize().createCompiledExpression(), *this));
        particleTerms.push_back(CpuCustomManyParticleForce::ParticleTermInfo(yname.str(), i, 1, energyExpr.differentiate(yname.str()).optimize().createCompiledExpression(), *this));
        particleTerms.push_back(CpuCustomManyParticleForce::ParticleTermInfo(zname.str(), i, 2, energyExpr.differentiate(zname.str()).optimize().createCompiledExpression(), *this));
        for (int j = 0; j < numPerParticleParameters; j++) {
            stringstream paramname;
            paramname << force.getPerParticleParameterName(j) << (i+1);
            particleParamIndices[i].push_back(expressionSet.getVariableIndex(paramname.str()));
        }
    }
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter)
        dihedralTerms.push_back(CpuCustomManyParticleForce::DihedralTermInfo(iter->first, iter->second, energyExpr.differentiate(iter->first).optimize().createCompiledExpression(), *this));
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter)
        distanceTerms.push_back(CpuCustomManyParticleForce::DistanceTermInfo(iter->first, iter->second, energyExpr.differentiate(iter->first).optimize().createCompiledExpression(), *this));
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter)
        angleTerms.push_back(CpuCustomManyParticleForce::AngleTermInfo(iter->first, iter->second, energyExpr.differentiate(iter->first).optimize().createCompiledExpression(), *this));
    for (int i = 0; i < particleTerms.size(); i++)
        expressionSet.registerExpression(particleTerms[i].forceExpression);
    for (int i = 0; i < distanceTerms.size(); i++)
        expressionSet.registerExpression(distanceTerms[i].forceExpression);
    for (int i = 0; i < angleTerms.size(); i++)
        expressionSet.registerExpression(angleTerms[i].forceExpression);
    for (int i = 0; i < dihedralTerms.size(); i++)
        expressionSet.registerExpression(dihedralTerms[i].forceExpression);
}

void CpuCustomManyParticleForce::ThreadData::requestDeltaPair(int p1, int p2, int& pairIndex, float& pairSign, bool allowReversed) {
    for (int i = 0; i < (int) deltaPairs.size(); i++) {
        if (deltaPairs[i].first == p1 && deltaPairs[i].second == p2) {
            pairIndex = i;
            pairSign = 1;
            return;
        }
        if (deltaPairs[i].first == p2 && deltaPairs[i].second == p1 && allowReversed) {
            pairIndex = i;
            pairSign = -1;
            return;
        }
    }
    pairIndex = deltaPairs.size();
    pairSign = 1;
    deltaPairs.push_back(make_pair(p1, p2));
}
