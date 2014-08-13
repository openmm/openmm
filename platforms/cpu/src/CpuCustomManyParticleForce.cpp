
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

using namespace OpenMM;
using namespace std;

CpuCustomManyParticleForce::CpuCustomManyParticleForce(const CustomManyParticleForce& force, ThreadPool& threads) :
            threads(threads), useCutoff(false), usePeriodic(false), neighborList(NULL) {
    numParticlesPerSet = force.getNumParticlesPerSet();
    numPerParticleParameters = force.getNumPerParticleParameters();
    
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the object used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpr = CustomManyParticleForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    energyExpression = energyExpr.createCompiledExpression();
    expressionSet.registerExpression(energyExpression);
    vector<string> particleParameterNames;
    if (force.getNonbondedMethod() != CustomManyParticleForce::NoCutoff)
        setUseCutoff(force.getCutoffDistance());

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
    
    // Differentiate the energy to get expressions for the force.

    particleParamIndices.resize(numParticlesPerSet);
    for (int i = 0; i < numParticlesPerSet; i++) {
        stringstream xname, yname, zname;
        xname << 'x' << (i+1);
        yname << 'y' << (i+1);
        zname << 'z' << (i+1);
        particleTerms.push_back(CpuCustomManyParticleForce::ParticleTermInfo(xname.str(), i, 0, energyExpr.differentiate(xname.str()).optimize().createCompiledExpression(), expressionSet));
        particleTerms.push_back(CpuCustomManyParticleForce::ParticleTermInfo(yname.str(), i, 1, energyExpr.differentiate(yname.str()).optimize().createCompiledExpression(), expressionSet));
        particleTerms.push_back(CpuCustomManyParticleForce::ParticleTermInfo(zname.str(), i, 2, energyExpr.differentiate(zname.str()).optimize().createCompiledExpression(), expressionSet));
        for (int j = 0; j < numPerParticleParameters; j++) {
            stringstream paramname;
            paramname << force.getPerParticleParameterName(j) << (i+1);
            particleParamIndices[i].push_back(expressionSet.getVariableIndex(paramname.str()));
        }
    }
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter)
        distanceTerms.push_back(CpuCustomManyParticleForce::DistanceTermInfo(iter->first, iter->second, energyExpr.differentiate(iter->first).optimize().createCompiledExpression(), expressionSet));
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter)
        angleTerms.push_back(CpuCustomManyParticleForce::AngleTermInfo(iter->first, iter->second, energyExpr.differentiate(iter->first).optimize().createCompiledExpression(), expressionSet));
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter)
        dihedralTerms.push_back(CpuCustomManyParticleForce::DihedralTermInfo(iter->first, iter->second, energyExpr.differentiate(iter->first).optimize().createCompiledExpression(), expressionSet));
    for (int i = 0; i < particleTerms.size(); i++)
        expressionSet.registerExpression(particleTerms[i].forceExpression);
    for (int i = 0; i < distanceTerms.size(); i++)
        expressionSet.registerExpression(distanceTerms[i].forceExpression);
    for (int i = 0; i < angleTerms.size(); i++)
        expressionSet.registerExpression(angleTerms[i].forceExpression);
    for (int i = 0; i < dihedralTerms.size(); i++)
        expressionSet.registerExpression(dihedralTerms[i].forceExpression);
    
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
}

void CpuCustomManyParticleForce::calculateIxn(AlignedArray<float>& posq, vector<RealVec>& atomCoordinates, RealOpenMM** particleParameters,
                                                  const map<string, double>& globalParameters, vector<RealVec>& forces,
                                                  RealOpenMM* totalEnergy) {
    for (map<string, double>::const_iterator iter = globalParameters.begin(); iter != globalParameters.end(); ++iter)
        expressionSet.setVariable(expressionSet.getVariableIndex(iter->first), iter->second);
    int numParticles = atomCoordinates.size();
    vector<int> particleIndices(numParticlesPerSet);
    fvec4 boxSize(periodicBoxSize[0], periodicBoxSize[1], periodicBoxSize[2], 0);
    fvec4 invBoxSize((1/periodicBoxSize[0]), (1/periodicBoxSize[1]), (1/periodicBoxSize[2]), 0);
    if (useCutoff) {
        // Construct a neighbor list.
        
        float boxSizeFloat[] = {(float) periodicBoxSize[0], (float) periodicBoxSize[1], (float) periodicBoxSize[2]};
        neighborList->computeNeighborList(numParticles, posq, exclusions, boxSizeFloat, usePeriodic, cutoffDistance, threads);
        for (int blockIndex = 0; blockIndex < neighborList->getNumBlocks(); blockIndex++) {
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
                loopOverInteractions(particles, particleIndices, 1, 0, &posq[0], atomCoordinates, particleParameters, forces, totalEnergy, boxSize, invBoxSize);
            }
        }
    }
    else {
        // Loop over all possible sets of particles.
        
        vector<int> particles(numParticles);
        for (int i = 0; i < numParticles; i++)
            particles[i] = i;
        for (int i = 0; i < numParticles; i++) {
            particleIndices[0] = i;
            loopOverInteractions(particles, particleIndices, 1, i+1, &posq[0], atomCoordinates, particleParameters, forces, totalEnergy, boxSize, invBoxSize);
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

void CpuCustomManyParticleForce::loopOverInteractions(vector<int>& availableParticles, vector<int>& particleSet, int loopIndex, int startIndex, float* posq, vector<OpenMM::RealVec>& atomCoordinates,
                                                          RealOpenMM** particleParameters, vector<OpenMM::RealVec>& forces, RealOpenMM* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
    int numParticles = availableParticles.size();
//    double cutoff2 = cutoffDistance*cutoffDistance;
    for (int i = startIndex; i < numParticles; i++) {
        int particle = availableParticles[i];
        
        // Check whether this particle can actually participate in interactions with the others found so far.
        
        bool include = true;
        if (useCutoff) {
//            fvec4 deltaR;
//            fvec4 pos1(posq+4*particle);
//            float r2;
//            for (int j = 0; j < loopIndex && include; j++) {
//                fvec4 pos2(posq+4*particleSet[j]);
//                getDeltaR(pos1, pos2, deltaR, r2, boxSize, invBoxSize);
//                include &= (r2 < cutoff2);
//            }
            RealOpenMM delta[ReferenceForce::LastDeltaRIndex];
            for (int j = 0; j < loopIndex && include; j++) {
                computeDelta(particle, particleSet[j], delta, atomCoordinates);
                include &= (delta[ReferenceForce::RIndex] < cutoffDistance);
            }
        }
        for (int j = 0; j < loopIndex && include; j++)
            include &= (exclusions[particle].find(particleSet[j]) == exclusions[particle].end());
        if (include) {
            particleSet[loopIndex] = availableParticles[i];
            if (loopIndex == numParticlesPerSet-1)
                calculateOneIxn(particleSet, posq, atomCoordinates, particleParameters, forces, totalEnergy, boxSize, invBoxSize);
            else
                loopOverInteractions(availableParticles, particleSet, loopIndex+1, i+1, posq, atomCoordinates, particleParameters, forces, totalEnergy, boxSize, invBoxSize);
        }
    }
}

void CpuCustomManyParticleForce::calculateOneIxn(vector<int>& particleSet, float* posq, vector<RealVec>& atomCoordinates, RealOpenMM** particleParameters, vector<RealVec>& forces, RealOpenMM* totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) {
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
    
    for (int i = 0; i < numParticlesPerSet; i++)
        for (int j = 0; j < numPerParticleParameters; j++)
            expressionSet.setVariable(particleParamIndices[i][j], particleParameters[permutedParticles[i]][j]);
    
    // Compute all of the variables the energy can depend on.

    for (int i = 0; i < (int) particleTerms.size(); i++) {
        const ParticleTermInfo& term = particleTerms[i];
        expressionSet.setVariable(term.variableIndex, atomCoordinates[term.atom][term.component]);
    }
    for (int i = 0; i < (int) distanceTerms.size(); i++) {
        const DistanceTermInfo& term = distanceTerms[i];
        computeDelta(permutedParticles[term.p1], permutedParticles[term.p2], term.delta, atomCoordinates);
        expressionSet.setVariable(term.variableIndex, term.delta[ReferenceForce::RIndex]);
    }
    for (int i = 0; i < (int) angleTerms.size(); i++) {
        const AngleTermInfo& term = angleTerms[i];
        computeDelta(permutedParticles[term.p1], permutedParticles[term.p2], term.delta1, atomCoordinates);
        computeDelta(permutedParticles[term.p3], permutedParticles[term.p2], term.delta2, atomCoordinates);
        expressionSet.setVariable(term.variableIndex, computeAngle(term.delta1, term.delta2));
    }
    for (int i = 0; i < (int) dihedralTerms.size(); i++) {
        const DihedralTermInfo& term = dihedralTerms[i];
        computeDelta(permutedParticles[term.p2], permutedParticles[term.p1], term.delta1, atomCoordinates);
        computeDelta(permutedParticles[term.p2], permutedParticles[term.p3], term.delta2, atomCoordinates);
        computeDelta(permutedParticles[term.p4], permutedParticles[term.p3], term.delta3, atomCoordinates);
        RealOpenMM dotDihedral, signOfDihedral;
        RealOpenMM* crossProduct[] = {term.cross1, term.cross2};
        expressionSet.setVariable(term.variableIndex, ReferenceBondIxn::getDihedralAngleBetweenThreeVectors(term.delta1, term.delta2, term.delta3, crossProduct, &dotDihedral, term.delta1, &signOfDihedral, 1));
    }
    
    // Apply forces based on individual particle coordinates.
    
    for (int i = 0; i < (int) particleTerms.size(); i++) {
        const ParticleTermInfo& term = particleTerms[i];
        forces[permutedParticles[term.atom]][term.component] -= term.forceExpression.evaluate();
    }

    // Apply forces based on distances.

    for (int i = 0; i < (int) distanceTerms.size(); i++) {
        const DistanceTermInfo& term = distanceTerms[i];
        RealOpenMM dEdR = (RealOpenMM) (term.forceExpression.evaluate()/(term.delta[ReferenceForce::RIndex]));
        for (int i = 0; i < 3; i++) {
           RealOpenMM force  = -dEdR*term.delta[i];
           forces[permutedParticles[term.p1]][i] -= force;
           forces[permutedParticles[term.p2]][i] += force;
        }
    }

    // Apply forces based on angles.

    for (int i = 0; i < (int) angleTerms.size(); i++) {
        const AngleTermInfo& term = angleTerms[i];
        RealOpenMM dEdTheta = (RealOpenMM) term.forceExpression.evaluate();
        RealOpenMM thetaCross[ReferenceForce::LastDeltaRIndex];
        SimTKOpenMMUtilities::crossProductVector3(term.delta1, term.delta2, thetaCross);
        RealOpenMM lengthThetaCross = SQRT(DOT3(thetaCross, thetaCross));
        if (lengthThetaCross < 1.0e-06)
            lengthThetaCross = (RealOpenMM) 1.0e-06;
        RealOpenMM termA = dEdTheta/(term.delta1[ReferenceForce::R2Index]*lengthThetaCross);
        RealOpenMM termC = -dEdTheta/(term.delta2[ReferenceForce::R2Index]*lengthThetaCross);
        RealOpenMM deltaCrossP[3][3];
        SimTKOpenMMUtilities::crossProductVector3(term.delta1, thetaCross, deltaCrossP[0]);
        SimTKOpenMMUtilities::crossProductVector3(term.delta2, thetaCross, deltaCrossP[2]);
        for (int i = 0; i < 3; i++) {
            deltaCrossP[0][i] *= termA;
            deltaCrossP[2][i] *= termC;
            deltaCrossP[1][i] = -(deltaCrossP[0][i]+deltaCrossP[2][i]);
        }
        for (int i = 0; i < 3; i++) {
            forces[permutedParticles[term.p1]][i] += deltaCrossP[0][i];
            forces[permutedParticles[term.p2]][i] += deltaCrossP[1][i];
            forces[permutedParticles[term.p3]][i] += deltaCrossP[2][i];
        }
    }

    // Apply forces based on dihedrals.

    for (int i = 0; i < (int) dihedralTerms.size(); i++) {
        const DihedralTermInfo& term = dihedralTerms[i];
        RealOpenMM dEdTheta = (RealOpenMM) term.forceExpression.evaluate();
        RealOpenMM internalF[4][3];
        RealOpenMM forceFactors[4];
        RealOpenMM normCross1 = DOT3(term.cross1, term.cross1);
        RealOpenMM normBC = term.delta2[ReferenceForce::RIndex];
        forceFactors[0] = (-dEdTheta*normBC)/normCross1;
        RealOpenMM normCross2 = DOT3(term.cross2, term.cross2);
                   forceFactors[3] = (dEdTheta*normBC)/normCross2;
                   forceFactors[1] = DOT3(term.delta1, term.delta2);
                   forceFactors[1] /= term.delta2[ReferenceForce::R2Index];
                   forceFactors[2] = DOT3(term.delta3, term.delta2);
                   forceFactors[2] /= term.delta2[ReferenceForce::R2Index];
        for (int i = 0; i < 3; i++) {
            internalF[0][i] = forceFactors[0]*term.cross1[i];
            internalF[3][i] = forceFactors[3]*term.cross2[i];
            RealOpenMM s = forceFactors[1]*internalF[0][i] - forceFactors[2]*internalF[3][i];
            internalF[1][i] = internalF[0][i] - s;
            internalF[2][i] = internalF[3][i] + s;
        }
        for (int i = 0; i < 3; i++) {
            forces[permutedParticles[term.p1]][i] += internalF[0][i];
            forces[permutedParticles[term.p2]][i] -= internalF[1][i];
            forces[permutedParticles[term.p3]][i] -= internalF[2][i];
            forces[permutedParticles[term.p4]][i] += internalF[3][i];
        }
    }

    // Add the energy

    if (totalEnergy)
        *totalEnergy += (RealOpenMM) energyExpression.evaluate();
}

void CpuCustomManyParticleForce::computeDelta(int atom1, int atom2, RealOpenMM* delta, vector<RealVec>& atomCoordinates) const {
    if (usePeriodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom1], atomCoordinates[atom2], periodicBoxSize, delta);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom1], atomCoordinates[atom2], delta);
}

RealOpenMM CpuCustomManyParticleForce::computeAngle(RealOpenMM* vec1, RealOpenMM* vec2) {
    RealOpenMM dot = DOT3(vec1, vec2);
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
