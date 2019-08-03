
/* Portions copyright (c) 2009-2018 Stanford University and Simbios.
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

#include "SimTKOpenMMUtilities.h"
#include "ReferenceForce.h"
#include "CpuCustomManyParticleForce.h"
#include "ReferenceTabulatedFunction.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "lepton/CustomFunction.h"

using namespace OpenMM;
using namespace std;

CpuCustomManyParticleForce::CpuCustomManyParticleForce(const CustomManyParticleForce& force, ThreadPool& threads) :
            threads(threads), useCutoff(false), usePeriodic(false), neighborList(NULL) {
    numParticles = force.getNumParticles();
    numParticlesPerSet = force.getNumParticlesPerSet();
    numPerParticleParameters = force.getNumPerParticleParameters();
    centralParticleMode = (force.getPermutationMode() == CustomManyParticleForce::UniqueCentralParticle);
    
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

    for (auto& function : functions)
        delete function.second;
    
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
    for (auto data : threadData)
        delete data;
}

void CpuCustomManyParticleForce::calculateIxn(AlignedArray<float>& posq, vector<vector<double> >& particleParameters,
                                                  const map<string, double>& globalParameters, vector<AlignedArray<float> >& threadForce,
                                                  bool includeForces, bool includeEnergy, double& energy) {
    // Record the parameters for the threads.
    
    this->posq = &posq[0];
    this->particleParameters = &particleParameters[0];
    this->globalParameters = &globalParameters;
    this->threadForce = &threadForce;
    this->includeForces = includeForces;
    this->includeEnergy = includeEnergy;
    atomicCounter = 0;
    if (useCutoff) {
        // Construct a neighbor list.  We use CpuNeighborList to do this, but then copy the result
        // into a new data structure.  This is needed because in UniqueCentralParticle mode, the
        // the neighbor list needs to include symmetric pairs.
        
        particleNeighbors.resize(numParticles);
        for (int i = 0; i < numParticles; i++)
            particleNeighbors[i].clear();
        neighborList->computeNeighborList(numParticles, posq, exclusions, periodicBoxVectors, usePeriodic, cutoffDistance, threads);
        for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
            const vector<int>& neighbors = neighborList->getBlockNeighbors(blockIndex);
            const vector<char>& exclusions = neighborList->getBlockExclusions(blockIndex);
            int numNeighbors = neighbors.size();
            for (int i = 0; i < 4; i++) {
                int p1 = neighborList->getSortedAtoms()[4*blockIndex+i];
                for (int j = 0; j < numNeighbors; j++) {
                    if ((exclusions[j] & (1<<i)) == 0) {
                        int p2 = neighbors[j];
                        particleNeighbors[p1].push_back(p2);
                        if (centralParticleMode)
                            particleNeighbors[p2].push_back(p1);
                    }
                }
            }
        }
    }
    
    // Signal the threads to start running and wait for them to finish.
    
    threads.execute([&] (ThreadPool& threads, int threadIndex) { threadComputeForce(threads, threadIndex); });
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
    fvec4 boxSize(periodicBoxVectors[0][0], periodicBoxVectors[1][1], periodicBoxVectors[2][2], 0);
    fvec4 invBoxSize(recipBoxSize[0], recipBoxSize[1], recipBoxSize[2], 0);
    float* forces = &(*threadForce)[threadIndex][0];
    ThreadData& data = *threadData[threadIndex];
    data.energy = 0;
    for (auto& param : *globalParameters)
        data.expressionSet.setVariable(data.expressionSet.getVariableIndex(param.first), param.second);
    if (useCutoff) {
        // Loop over interactions from the neighbor list.
        
        while (true) {
            int i = atomicCounter++;
            if (i >= numParticles)
                break;
            particleIndices[0] = i;
            loopOverInteractions(particleNeighbors[i], particleIndices, 1, 0, particleParameters, forces, data, boxSize, invBoxSize);
        }
    }
    else {
        // Loop over all possible sets of particles.
        
        vector<int> particles(numParticles);
        for (int i = 0; i < numParticles; i++)
            particles[i] = i;
        while (true) {
            int i = atomicCounter++;
            if (i >= numParticles)
                break;
            particleIndices[0] = i;
            int startIndex = (centralParticleMode ? 0 : i+1);
            loopOverInteractions(particles, particleIndices, 1, startIndex, particleParameters, forces, data, boxSize, invBoxSize);
        }
    }
}

void CpuCustomManyParticleForce::setUseCutoff(double distance) {
    useCutoff = true;
    cutoffDistance = distance;
    if (neighborList == NULL)
        neighborList = new CpuNeighborList(4);
}

void CpuCustomManyParticleForce::setPeriodic(Vec3* periodicBoxVectors) {
    assert(useCutoff);
    assert(periodicBoxVectors[0][0] >= 2.0*cutoffDistance);
    assert(periodicBoxVectors[1][1] >= 2.0*cutoffDistance);
    assert(periodicBoxVectors[2][2] >= 2.0*cutoffDistance);
    usePeriodic = true;
    this->periodicBoxVectors[0] = periodicBoxVectors[0];
    this->periodicBoxVectors[1] = periodicBoxVectors[1];
    this->periodicBoxVectors[2] = periodicBoxVectors[2];
    recipBoxSize[0] = (float) (1.0/periodicBoxVectors[0][0]);
    recipBoxSize[1] = (float) (1.0/periodicBoxVectors[1][1]);
    recipBoxSize[2] = (float) (1.0/periodicBoxVectors[2][2]);
    periodicBoxVec4.resize(3);
    periodicBoxVec4[0] = fvec4(periodicBoxVectors[0][0], periodicBoxVectors[0][1], periodicBoxVectors[0][2], 0);
    periodicBoxVec4[1] = fvec4(periodicBoxVectors[1][0], periodicBoxVectors[1][1], periodicBoxVectors[1][2], 0);
    periodicBoxVec4[2] = fvec4(periodicBoxVectors[2][0], periodicBoxVectors[2][1], periodicBoxVectors[2][2], 0);
    triclinic = (periodicBoxVectors[0][1] != 0.0 || periodicBoxVectors[0][2] != 0.0 ||
                 periodicBoxVectors[1][0] != 0.0 || periodicBoxVectors[1][2] != 0.0 ||
                 periodicBoxVectors[2][0] != 0.0 || periodicBoxVectors[2][1] != 0.0);
}

void CpuCustomManyParticleForce::loopOverInteractions(vector<int>& availableParticles, vector<int>& particleSet, int loopIndex, int startIndex,
                                                      vector<double>* particleParameters, float* forces, ThreadData& data, const fvec4& boxSize, const fvec4& invBoxSize) {
    int numParticles = availableParticles.size();
    double cutoff2 = cutoffDistance*cutoffDistance;
    int checkRange = (centralParticleMode ? 1 : loopIndex);
    for (int i = startIndex; i < numParticles; i++) {
        int particle = availableParticles[i];
        
        // Check whether this particle can actually participate in interactions with the others found so far.
        
        bool include = true;
        if (useCutoff) {
            fvec4 deltaR;
            fvec4 pos1(posq+4*particle);
            float r2;
            for (int j = 0; j < checkRange && include; j++) {
                fvec4 pos2(posq+4*particleSet[j]);
                computeDelta(pos1, pos2, deltaR, r2, boxSize, invBoxSize);
                include &= (r2 < cutoff2);
            }
        }
        for (int j = 0; j < loopIndex && include; j++)
            include &= (exclusions[particle].find(particleSet[j]) == exclusions[particle].end());
        if (include) {
            if (loopIndex > 0 && availableParticles[i] == particleSet[0])
                continue;
            particleSet[loopIndex] = availableParticles[i];
            if (loopIndex == numParticlesPerSet-1)
                calculateOneIxn(particleSet, particleParameters, forces, data, boxSize, invBoxSize);
            else
                loopOverInteractions(availableParticles, particleSet, loopIndex+1, i+1, particleParameters, forces, data, boxSize, invBoxSize);
        }
    }
}

void CpuCustomManyParticleForce::calculateOneIxn(vector<int>& particleSet, vector<double>* particleParameters, float* forces, ThreadData& data, const fvec4& boxSize, const fvec4& invBoxSize) {
    // Select the ordering to use for the particles.
    
    vector<int>& permutedParticles = data.permutedParticles;
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
    AlignedArray<fvec4>& delta = data.delta;
    AlignedArray<fvec4>& cross1 = data.cross1;
    AlignedArray<fvec4>& cross2 = data.cross2;
    vector<float>& normDelta = data.normDelta;
    vector<float>& norm2Delta = data.norm2Delta;
    for (int i = 0; i < numDeltas; i++) {
        int p1 = permutedParticles[data.deltaPairs[i].first];
        int p2 = permutedParticles[data.deltaPairs[i].second];
        computeDelta(fvec4(posq+4*p1), fvec4(posq+4*p2), delta[i], norm2Delta[i], boxSize, invBoxSize);
        normDelta[i] = sqrtf(norm2Delta[i]);
    }
    
    // Compute all of the variables the energy can depend on.

    for (auto& term : data.particleTerms)
        expressionSet.setVariable(term.variableIndex, posq[4*permutedParticles[term.atom]+term.component]);
    for (auto& term : data.distanceTerms)
        expressionSet.setVariable(term.variableIndex, normDelta[term.delta]);
    for (auto& term : data.angleTerms)
        expressionSet.setVariable(term.variableIndex, computeAngle(delta[term.delta1], delta[term.delta2], norm2Delta[term.delta1], norm2Delta[term.delta2], term.delta1Sign*term.delta2Sign));
    for (int i = 0; i < (int) data.dihedralTerms.size(); i++) {
        const DihedralTermInfo& term = data.dihedralTerms[i];
        expressionSet.setVariable(term.variableIndex, getDihedralAngleBetweenThreeVectors(delta[term.delta1], delta[term.delta2], delta[term.delta3], cross1[i], cross2[i], delta[term.delta1]));
    }
    
    if (includeForces) {
        // Apply forces based on individual particle coordinates.

        AlignedArray<fvec4>& f = data.f;
        for (int i = 0; i < numParticlesPerSet; i++)
            f[i] = fvec4(0.0f);
        for (auto& term : data.particleTerms) {
            float temp[4];
            f[term.atom].store(temp);
            temp[term.component] -= term.forceExpression.evaluate();
            f[term.atom] = fvec4(temp);
        }

        // Apply forces based on distances.

        for (auto& term : data.distanceTerms) {
            float dEdR = (float) (term.forceExpression.evaluate()*term.deltaSign/(normDelta[term.delta]));
            fvec4 force = -dEdR*delta[term.delta];
            f[term.p1] -= force;
            f[term.p2] += force;
        }

        // Apply forces based on angles.

        for (auto& term : data.angleTerms) {
            float dEdTheta = (float) term.forceExpression.evaluate();
            fvec4 thetaCross = cross(delta[term.delta1], delta[term.delta2]);
            float lengthThetaCross = sqrtf(dot3(thetaCross, thetaCross));
            if (lengthThetaCross < 1.0e-6f)
                lengthThetaCross = 1.0e-6f;
            float termA = dEdTheta*term.delta2Sign/(norm2Delta[term.delta1]*lengthThetaCross);
            float termC = -dEdTheta*term.delta1Sign/(norm2Delta[term.delta2]*lengthThetaCross);
            fvec4 deltaCross1 = cross(delta[term.delta1], thetaCross);
            fvec4 deltaCross2 = cross(delta[term.delta2], thetaCross);
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
            float dEdTheta = (float) term.forceExpression.evaluate();
            float normCross1 = dot3(cross1[i], cross1[i]);
            float normBC = normDelta[term.delta2];
            float forceFactors[4];
            forceFactors[0] = (-dEdTheta*normBC)/normCross1;
            float normCross2 = dot3(cross2[i], cross2[i]);
            forceFactors[3] = (dEdTheta*normBC)/normCross2;
            forceFactors[1] = dot3(delta[term.delta1], delta[term.delta2]);
            forceFactors[1] /= norm2Delta[term.delta2];
            forceFactors[2] = dot3(delta[term.delta3], delta[term.delta2]);
            forceFactors[2] /= norm2Delta[term.delta2];
            fvec4 force1 = forceFactors[0]*cross1[i];
            fvec4 force4 = forceFactors[3]*cross2[i];
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

void CpuCustomManyParticleForce::computeDelta(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const {
    deltaR = posJ-posI;
    if (usePeriodic) {
        if (triclinic) {
            deltaR -= periodicBoxVec4[2]*floorf(deltaR[2]*recipBoxSize[2]+0.5f);
            deltaR -= periodicBoxVec4[1]*floorf(deltaR[1]*recipBoxSize[1]+0.5f);
            deltaR -= periodicBoxVec4[0]*floorf(deltaR[0]*recipBoxSize[0]+0.5f);
        }
        else {
            fvec4 base = round(deltaR*invBoxSize)*boxSize;
            deltaR = deltaR-base;
        }
    }
    r2 = dot3(deltaR, deltaR);
}

float CpuCustomManyParticleForce::computeAngle(const fvec4& vi, const fvec4& vj, float v2i, float v2j, float sign) {
    float dot = dot3(vi, vj)*sign;
    float cosine = dot/sqrtf(v2i*v2j);
    if (cosine > 0.99f || cosine < -0.99f) {
        // We're close to the singularity in acos(), so take the cross product and use asin() instead.

        fvec4 cross12 = cross(vi, vj);
        float scale = v2i*v2j;
        float angle = asinf(sqrtf(dot3(cross12, cross12)/scale));
        if (cosine < 0.0f)
            angle = (float) (M_PI-angle);
        return angle;
    }
    return acosf(cosine);
}

float CpuCustomManyParticleForce::getDihedralAngleBetweenThreeVectors(const fvec4& v1, const fvec4& v2, const fvec4& v3, fvec4& cross1, fvec4& cross2, const fvec4& signVector) {
    cross1 = cross(v1, v2);
    cross2 = cross(v2, v3);
    float angle = computeAngle(cross1, cross2, dot3(cross1, cross1), dot3(cross2, cross2), 1.0f);
    float dotProduct = dot3(signVector, cross2);
    if (dotProduct < 0) 
        angle = -angle;
    return angle;
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
    permutedParticles.resize(numParticlesPerSet);
    f.resize(numParticlesPerSet);
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
    for (auto& term : dihedrals)
        dihedralTerms.push_back(CpuCustomManyParticleForce::DihedralTermInfo(term.first, term.second, energyExpr.differentiate(term.first).optimize().createCompiledExpression(), *this));
    for (auto& term : distances)
        distanceTerms.push_back(CpuCustomManyParticleForce::DistanceTermInfo(term.first, term.second, energyExpr.differentiate(term.first).optimize().createCompiledExpression(), *this));
    for (auto& term : angles)
        angleTerms.push_back(CpuCustomManyParticleForce::AngleTermInfo(term.first, term.second, energyExpr.differentiate(term.first).optimize().createCompiledExpression(), *this));
    for (auto& term : particleTerms)
        expressionSet.registerExpression(term.forceExpression);
    for (auto& term : distanceTerms)
        expressionSet.registerExpression(term.forceExpression);
    for (auto& term : angleTerms)
        expressionSet.registerExpression(term.forceExpression);
    for (auto& term : dihedralTerms)
        expressionSet.registerExpression(term.forceExpression);
    int numDeltas = deltaPairs.size();
    delta.resize(numDeltas);
    normDelta.resize(numDeltas);
    norm2Delta.resize(numDeltas);
    cross1.resize(numDeltas);
    cross2.resize(numDeltas);
    
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
