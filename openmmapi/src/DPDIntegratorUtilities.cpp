/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
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

#include "openmm/internal/DPDIntegratorUtilities.h"
#include "openmm/OpenMMException.h"
#include <map>
#include <cstdio>

using namespace OpenMM;
using namespace std;

void DPDIntegratorUtilities::createTypeTables(const DPDIntegrator& integrator, int numParticles, int& numTypes, vector<int>& particleType,
        vector<vector<double> >& friction, vector<vector<double> >& cutoff, double& maxCutoff) {
    // Identify all unique types and record the types of particles.

    map<int, int> typeMap;
    particleType.resize(numParticles);
    for (int i = 0; i < numParticles; i++) {
        int type = integrator.getParticleType(i);
        auto mapping = typeMap.find(type);
        if (mapping == typeMap.end()) {
            particleType[i] = typeMap.size();
            typeMap[type] = particleType[i];
        }
        else
            particleType[i] = mapping->second;
    }
    numTypes = typeMap.size();

    // Build tables of friction and cutoff.

    friction.clear();
    cutoff.clear();
    friction.resize(numTypes, vector<double>(numTypes, integrator.getDefaultFriction()));
    cutoff.resize(numTypes, vector<double>(numTypes, integrator.getDefaultCutoff()));
    map<int, int> inverseMap;
    for (auto mapping : typeMap)
        inverseMap[mapping.second] = mapping.first;
    maxCutoff = integrator.getDefaultCutoff();
    for (int i = 0; i < integrator.getNumTypePairs(); i++) {
        int t1, t2;
        double f, c;
        integrator.getTypePairParameters(i, t1, t2, f, c);
        if (f < 0.0)
            throw OpenMMException("DPDIntegrator: friction cannot be negative");
        if (c <= 0.0)
            throw OpenMMException("DPDIntegrator: cutoff must be positive");
        int type1 = inverseMap[t1];
        int type2 = inverseMap[t2];
        friction[type1][type2] = f;
        friction[type2][type1] = f;
        cutoff[type1][type2] = c;
        cutoff[type2][type1] = c;
        if (c > maxCutoff)
            maxCutoff = c;
    }
}
