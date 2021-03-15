/* Portions copyright (c) 2009-2021 Stanford University and Simbios.
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
#include "ReferencePointFunctions.h"
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

    // Create implementations of point functions.

    functions["pointdistance"] = new ReferencePointDistanceFunction(force.usesPeriodicBoundaryConditions(), &periodicBoxVectors);
    functions["pointangle"] = new ReferencePointAngleFunction(force.usesPeriodicBoundaryConditions(), &periodicBoxVectors);
    functions["pointdihedral"] = new ReferencePointDihedralFunction(force.usesPeriodicBoundaryConditions(), &periodicBoxVectors);

    // Parse the expression and create the object used to calculate the interaction.

    Lepton::ParsedExpression energyExpr = CustomManyParticleForceImpl::prepareExpression(force, functions);
    energyExpression = energyExpr.createProgram();
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
    periodicBoxVectors = vectors;
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
    
    // Record particle coordinates.

    for (auto& term : particleTerms)
        variables[term.name] = atomCoordinates[permutedParticles[term.atom]][term.component];

    // Apply forces based on particle coordinates.

    for (auto& term : particleTerms)
        forces[permutedParticles[term.atom]][term.component] -= term.forceExpression.evaluate(variables);

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
