#ifndef OPENMM_REFERENCECONSTANTPOTENTIAL_H_
#define OPENMM_REFERENCECONSTANTPOTENTIAL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
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

class ReferenceConstantPotential;

/**
 * A generic charge solver for the constant potential method.
 */
class ReferenceConstantPotentialSolver {
protected:
    bool valid;

public:
    /**
     * Creates a ReferenceConstantPotentialSolver.
     */
    ReferenceConstantPotentialSolver();
    virtual ~ReferenceConstantPotentialSolver();
    /**
     * Mark precomputed data stored by the solver as invalid due to a change in
     * electrode parameters.
     */
    void invalidate();
    /**
     * Updates precomputed data stored by the solver.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param conp                   constant potential derivative evaluation class
     */
    virtual void update(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, const std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotential& conp) = 0;
    /**
     * Solves for charges.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                output particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param conp                   constant potential derivative evaluation class
     * @param pmeData                reference PME solver
     */
    virtual void solve(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotential& conp, pme_t pmeData) = 0;
};

/**
 * A constant potential solver using direct inversion of the Coulomb matrix.
 * Suitable only when electrode particle positions are fixed.
 */
class ReferenceConstantPotentialMatrixSolver : public ReferenceConstantPotentialSolver {
private:
    Vec3 boxVectors[3];
    std::vector<Vec3> electrodePosData;
    JAMA::Cholesky<double> capacitance;
    std::vector<double> constraintVector;

public:
    /**
     * Creates a ReferenceConstantPotentialMatrixSolver.
     * 
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     */
    ReferenceConstantPotentialMatrixSolver(int numElectrodeParticles);
    /**
     * Updates precomputed data stored by the solver.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param conp                   constant potential derivative evaluation class
     */
    void update(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, const std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotential& conp);
    /**
     * Solves for charges.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                output particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param conp                   constant potential derivative evaluation class
     * @param pmeData                reference PME solver
     */
    void solve(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotential& conp, pme_t pmeData);
};

/**
 * A constant potential solver using the conjugate gradient method.  Suitable
 * for both fixed and variable electrode particle positions.
 */
class ReferenceConstantPotentialCGSolver : public ReferenceConstantPotentialSolver {
private:
    Vec3 boxVectors[3];
    bool precond;
    std::vector<double> precondVector;

public:
    /**
     * Creates a ReferenceConstantPotentialCGSolver.
     * 
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param precond                whether or not to use a preconditioner
     */
    ReferenceConstantPotentialCGSolver(int numElectrodeParticles, bool precond);
    /**
     * Updates precomputed data stored by the solver.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param conp                   constant potential derivative evaluation class
     */
    void update(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, const std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotential& conp);
    /**
     * Solves for charges.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                output particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param conp                   constant potential derivative evaluation class
     * @param pmeData                reference PME solver
     */
    void solve(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotential& conp, pme_t pmeData);
};

/**
 * Performs energy, force, and charge derivative calculations for the reference
 * kernel for ConstantPotentialForce.
 */
class ReferenceConstantPotential {
    friend class ReferenceConstantPotentialMatrixSolver;
    friend class ReferenceConstantPotentialCGSolver;

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
    /**
     * Creates a ReferenceConstantPotential.
     * 
     * @param nonbondedCutoff        direct space cutoff
     * @param neighborList           neighbor list for direct space calculation
     * @param boxVectors             periodic box vectors
     * @param exceptionsArePeriodic  whether or not exceptions use periodic boundary conditions
     * @param ewaldAlpha             Ewald reciprocal Gaussian width parameter
     * @param gridSize               Ewald mesh dimensions
     * @param cgErrorTol             constant potential conjugate gradient error tolerance
     * @param useChargeConstraint    whether or not to constrain total charge
     * @param chargeTarget           target sum of charges on electrode particles only
     * @param externalField          electric field vector
     */
    ReferenceConstantPotential(double nonbondedCutoff,
        const NeighborList* neighborList, const Vec3* boxVectors,
        bool exceptionsArePeriodic, double ewaldAlpha, const int* gridSize,
        double cgErrorTol, bool useChargeConstraint, double chargeTarget,
        Vec3 externalField);
    /**
     * Solves for charges and computes energies and forces.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param forceData              output forces on particles
     * @param charges                output particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param energy                 output system energy
     * @param solver                 charge solver implementation
     */
    void execute(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<Vec3>& forceData,
        std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        double* energy, ReferenceConstantPotentialSolver* solver);
    /**
     * Solves for charges without computing energies and forces.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                output particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param solver                 charge solver implementation
     */
    void getCharges(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        ReferenceConstantPotentialSolver* solver);

private:
    /**
     * Computes energies and forces for fixed (solved) charges.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param forceData              output forces on particles
     * @param charges                particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param energy                 output system energy
     * @param pmeData                reference PME solver
     */
    void getEnergyForces(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<Vec3>& forceData,
        std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        double* energy, pme_t pmeData);
    /**
     * Computes energy derivatives with respect to charges.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posData                particle positions
     * @param charges                particle charges
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param electrodeParamArray    electrode particle parameters
     * @param chargeDerivatives      output charge derivatives
     * @param pmeData                reference PME solver
     */
    void getDerivatives(int numParticles, int numElectrodeParticles,
        const std::vector<Vec3>& posData, std::vector<double>& charges,
        const std::vector<std::set<int>>& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<std::array<double, 3> >& electrodeParamArray,
        std::vector<double>& chargeDerivatives, pme_t pmeData);
};

} // namespace OpenMM

#endif // OPENMM_REFERENCECONSTANTPOTENTIAL_H_
