
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
#include "ReferenceCustomManyParticleIxn.h"
#include "ReferenceTabulatedFunction.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "lepton/CustomFunction.h"

using namespace OpenMM;
using namespace std;

ReferenceCustomManyParticleIxn::ReferenceCustomManyParticleIxn(const CustomManyParticleForce& force) : useCutoff(false), usePeriodic(false) {
    numParticlesPerSet = force.getNumParticlesPerSet();
    numPerParticleParameters = force.getNumPerParticleParameters();
    centralParticleMode = (force.getPermutationMode() == CustomManyParticleForce::UniqueCentralParticle);
    
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the object used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpr = CustomManyParticleForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    energyExpression = energyExpr.createProgram();
    vector<string> particleParameterNames;
    if (force.getNonbondedMethod() != CustomManyParticleForce::NoCutoff)
        setUseCutoff(force.getCutoffDistance());

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;

    // Differentiate the energy to get expressions for the force.

    particleParamNames.resize(numParticlesPerSet);
    for (int i = 0; i < numParticlesPerSet; i++) {
        stringstream xname, yname, zname;
        xname << 'x' << (i+1);
        yname << 'y' << (i+1);
        zname << 'z' << (i+1);
        particleTerms.push_back(ReferenceCustomManyParticleIxn::ParticleTermInfo(xname.str(), i, 0, energyExpr.differentiate(xname.str()).optimize().createProgram()));
        particleTerms.push_back(ReferenceCustomManyParticleIxn::ParticleTermInfo(yname.str(), i, 1, energyExpr.differentiate(yname.str()).optimize().createProgram()));
        particleTerms.push_back(ReferenceCustomManyParticleIxn::ParticleTermInfo(zname.str(), i, 2, energyExpr.differentiate(zname.str()).optimize().createProgram()));
        for (int j = 0; j < numPerParticleParameters; j++) {
            stringstream paramname;
            paramname << force.getPerParticleParameterName(j) << (i+1);
            particleParamNames[i].push_back(paramname.str());
        }
    }
    for (auto& term : distances)
        distanceTerms.push_back(ReferenceCustomManyParticleIxn::DistanceTermInfo(term.first, term.second, energyExpr.differentiate(term.first).optimize().createProgram()));
    for (auto& term : angles)
        angleTerms.push_back(ReferenceCustomManyParticleIxn::AngleTermInfo(term.first, term.second, energyExpr.differentiate(term.first).optimize().createProgram()));
    for (auto& term : dihedrals)
        dihedralTerms.push_back(ReferenceCustomManyParticleIxn::DihedralTermInfo(term.first, term.second, energyExpr.differentiate(term.first).optimize().createProgram()));
    
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

ReferenceCustomManyParticleIxn::~ReferenceCustomManyParticleIxn() {
}

void ReferenceCustomManyParticleIxn::calculateIxn(vector<Vec3>& atomCoordinates, vector<vector<double> >& particleParameters,
                                                  const map<string, double>& globalParameters, vector<Vec3>& forces,
                                                  double* totalEnergy) const {
    map<string, double> variables = globalParameters;
    vector<int> particles(numParticlesPerSet);
    loopOverInteractions(particles, 0, atomCoordinates, particleParameters, variables, forces, totalEnergy);
}

void ReferenceCustomManyParticleIxn::setUseCutoff(double distance) {
    useCutoff = true;
    cutoffDistance = distance;
}

void ReferenceCustomManyParticleIxn::setPeriodic(Vec3* vectors) {
    assert(useCutoff);
    assert(vectors[0][0] >= 2.0*cutoffDistance);
    assert(vectors[1][1] >= 2.0*cutoffDistance);
    assert(vectors[2][2] >= 2.0*cutoffDistance);
    usePeriodic = true;
    periodicBoxVectors[0] = vectors[0];
    periodicBoxVectors[1] = vectors[1];
    periodicBoxVectors[2] = vectors[2];
}

void ReferenceCustomManyParticleIxn::loopOverInteractions(vector<int>& particles, int loopIndex, vector<OpenMM::Vec3>& atomCoordinates,
                                                          vector<vector<double> >& particleParameters, map<string, double>& variables, vector<OpenMM::Vec3>& forces,
                                                          double* totalEnergy) const {
    int numParticles = atomCoordinates.size();
    int firstPartialLoop = (centralParticleMode ? 2 : 1);
    int start = (loopIndex < firstPartialLoop ? 0 : particles[loopIndex-1]+1);
    for (int i = start; i < numParticles; i++) {
        if (loopIndex > 0 && i == particles[0])
            continue;
        particles[loopIndex] = i;
        if (loopIndex == numParticlesPerSet-1)
            calculateOneIxn(particles, atomCoordinates, particleParameters, variables, forces, totalEnergy);
        else
            loopOverInteractions(particles, loopIndex+1, atomCoordinates, particleParameters, variables, forces, totalEnergy);
    }
}

void ReferenceCustomManyParticleIxn::calculateOneIxn(const vector<int>& particles, vector<Vec3>& atomCoordinates,
                        vector<vector<double> >& particleParameters, map<string, double>& variables, vector<Vec3>& forces, double* totalEnergy) const {
    // Select the ordering to use for the particles.
    
    vector<int> permutedParticles(numParticlesPerSet);
    if (particleOrder.size() == 1) {
        // There are no filters, so we don't need to worry about ordering.
        
        permutedParticles = particles;
    }
    else {
        int index = 0;
        for (int i = numParticlesPerSet-1; i >= 0; i--)
            index = particleTypes[particles[i]]+numTypes*index;
        int order = orderIndex[index];
        if (order == -1)
            return;
        for (int i = 0; i < numParticlesPerSet; i++)
            permutedParticles[i] = particles[particleOrder[order][i]];
    }
    
    // Decide whether to include this interaction.
    
    for (int i = 0; i < numParticlesPerSet; i++) {
        int p1 = permutedParticles[i];
        for (int j = i+1; j < numParticlesPerSet; j++) {
            int p2 = permutedParticles[j];
            if (exclusions[p1].find(p2) != exclusions[p1].end())
                return;
            if (useCutoff && (i == 0 || !centralParticleMode)) {
                double delta[ReferenceForce::LastDeltaRIndex];
                computeDelta(p1, p2, delta, atomCoordinates);
                if (delta[ReferenceForce::RIndex] >= cutoffDistance)
                    return;
            }
        }
    }

    // Record per-particle parameters.
    
    for (int i = 0; i < numParticlesPerSet; i++)
        for (int j = 0; j < numPerParticleParameters; j++)
            variables[particleParamNames[i][j]] = particleParameters[permutedParticles[i]][j];
    
    // Compute all of the variables the energy can depend on.

    for (auto& term : particleTerms)
        variables[term.name] = atomCoordinates[permutedParticles[term.atom]][term.component];
    for (auto& term : distanceTerms) {
        computeDelta(permutedParticles[term.p1], permutedParticles[term.p2], term.delta, atomCoordinates);
        variables[term.name] = term.delta[ReferenceForce::RIndex];
    }
    for (auto& term : angleTerms) {
        computeDelta(permutedParticles[term.p1], permutedParticles[term.p2], term.delta1, atomCoordinates);
        computeDelta(permutedParticles[term.p3], permutedParticles[term.p2], term.delta2, atomCoordinates);
        variables[term.name] = computeAngle(term.delta1, term.delta2);
    }
    for (auto& term : dihedralTerms) {
        computeDelta(permutedParticles[term.p2], permutedParticles[term.p1], term.delta1, atomCoordinates);
        computeDelta(permutedParticles[term.p2], permutedParticles[term.p3], term.delta2, atomCoordinates);
        computeDelta(permutedParticles[term.p4], permutedParticles[term.p3], term.delta3, atomCoordinates);
        double dotDihedral, signOfDihedral;
        double* crossProduct[] = {term.cross1, term.cross2};
        variables[term.name] = ReferenceBondIxn::getDihedralAngleBetweenThreeVectors(term.delta1, term.delta2, term.delta3, crossProduct, &dotDihedral, term.delta1, &signOfDihedral, 1);
    }
    
    // Apply forces based on individual particle coordinates.
    
    for (auto& term : particleTerms)
        forces[permutedParticles[term.atom]][term.component] -= term.forceExpression.evaluate(variables);

    // Apply forces based on distances.

    for (auto& term : distanceTerms) {
        double dEdR = term.forceExpression.evaluate(variables)/(term.delta[ReferenceForce::RIndex]);
        for (int i = 0; i < 3; i++) {
           double force  = -dEdR*term.delta[i];
           forces[permutedParticles[term.p1]][i] -= force;
           forces[permutedParticles[term.p2]][i] += force;
        }
    }

    // Apply forces based on angles.

    for (auto& term : angleTerms) {
        double dEdTheta = term.forceExpression.evaluate(variables);
        double thetaCross[ReferenceForce::LastDeltaRIndex];
        SimTKOpenMMUtilities::crossProductVector3(term.delta1, term.delta2, thetaCross);
        double lengthThetaCross = sqrt(DOT3(thetaCross, thetaCross));
        if (lengthThetaCross < 1.0e-06)
            lengthThetaCross = 1.0e-06;
        double termA = dEdTheta/(term.delta1[ReferenceForce::R2Index]*lengthThetaCross);
        double termC = -dEdTheta/(term.delta2[ReferenceForce::R2Index]*lengthThetaCross);
        double deltaCrossP[3][3];
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

    for (auto& term : dihedralTerms) {
        double dEdTheta = term.forceExpression.evaluate(variables);
        double internalF[4][3];
        double forceFactors[4];
        double normCross1 = DOT3(term.cross1, term.cross1);
        double normBC = term.delta2[ReferenceForce::RIndex];
        forceFactors[0] = (-dEdTheta*normBC)/normCross1;
        double normCross2 = DOT3(term.cross2, term.cross2);
        forceFactors[3] = (dEdTheta*normBC)/normCross2;
        forceFactors[1] = DOT3(term.delta1, term.delta2);
        forceFactors[1] /= term.delta2[ReferenceForce::R2Index];
        forceFactors[2] = DOT3(term.delta3, term.delta2);
        forceFactors[2] /= term.delta2[ReferenceForce::R2Index];
        for (int i = 0; i < 3; i++) {
            internalF[0][i] = forceFactors[0]*term.cross1[i];
            internalF[3][i] = forceFactors[3]*term.cross2[i];
            double s = forceFactors[1]*internalF[0][i] - forceFactors[2]*internalF[3][i];
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
        *totalEnergy += energyExpression.evaluate(variables);
}

void ReferenceCustomManyParticleIxn::computeDelta(int atom1, int atom2, double* delta, vector<Vec3>& atomCoordinates) const {
    if (usePeriodic)
        ReferenceForce::getDeltaRPeriodic(atomCoordinates[atom1], atomCoordinates[atom2], periodicBoxVectors, delta);
    else
        ReferenceForce::getDeltaR(atomCoordinates[atom1], atomCoordinates[atom2], delta);
}

double ReferenceCustomManyParticleIxn::computeAngle(double* vec1, double* vec2) {
    double dot = DOT3(vec1, vec2);
    double cosine = dot/sqrt((vec1[ReferenceForce::R2Index]*vec2[ReferenceForce::R2Index]));
    double angle;
    if (cosine >= 1)
        angle = 0;
    else if (cosine <= -1)
        angle = PI_M;
    else
        angle = acos(cosine);
    return angle;
}
