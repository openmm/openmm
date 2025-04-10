#ifndef OPENMM_REFERENCECONSTANTPOTENTIAL_H_
#define OPENMM_REFERENCECONSTANTPOTENTIAL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
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

#include <array>

#include "ReferenceNeighborList.h"
#include "ReferencePME.h"
#include "tnt_array2d.h"
#include "jama_cholesky.h"

namespace OpenMM {

class ReferenceConstantPotentialMatrix {
    friend class ReferenceConstantPotential;

private:
    bool valid;
    Vec3 boxVectors[3];
    JAMA::Cholesky<double> capacitance;
    std::vector<double> constraintVector;

public:
    ReferenceConstantPotentialMatrix(int numElectrodeParticles);
    void invalidate();
};

class ReferenceConstantPotentialCG {
    friend class ReferenceConstantPotential;

private:
    bool valid;
    Vec3 boxVectors[3];
    std::vector<double> precondVector;
    double precondScale;

public:
    ReferenceConstantPotentialCG(int numElectrodeParticles);
    void invalidate();
};

class ReferenceConstantPotential {
private:
    double nonbondedCutoff, ewaldAlpha, cgErrorTol, chargeTarget;
    const NeighborList* neighborList;
    Vec3 boxVectors[3], externalField;
    bool exceptionsArePeriodic, useChargeConstraint;
    int gridSize[3];

    static const int PotentialIndex = 0;
    static const int GaussianWidthIndex = 1;
    static const int ThomasFermiScaleIndex = 2;

public:
    ReferenceConstantPotential(double nonbondedCutoff,
        const NeighborList* neighborList, const Vec3* boxVectors,
        bool exceptionsArePeriodic, double ewaldAlpha, const int* gridSize,
        double cgErrorTol, bool useChargeConstraint, double chargeTarget,
        Vec3 externalField);
    void updateMatrix(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, const std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotentialMatrix* matrix);
    void updateCG(int numParticles, int numElectrodeParticles,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotentialCG* cg);
    void execute(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<Vec3>& forceData,
        std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        double* energy, ReferenceConstantPotentialMatrix* matrix,
        ReferenceConstantPotentialCG* cg);
    void getCharges(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotentialMatrix* matrix,
        ReferenceConstantPotentialCG* cg);

private:
    void solve(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotentialMatrix* matrix,
        ReferenceConstantPotentialCG* cg, pme_t pmeData);
    void getEnergyForces(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<Vec3>& forceData,
        std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        double* energy, pme_t pmeData);
    void getDerivatives(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& electrodeIndexMap,
        const std::vector<int>& electrodeIndices,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        std::vector<double>& chargeDerivatives, pme_t pmeData);
};

} // namespace OpenMM

#endif // OPENMM_REFERENCECONSTANTPOTENTIAL_H_
