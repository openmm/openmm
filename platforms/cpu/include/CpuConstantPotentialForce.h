#ifndef OPENMM_CPUCONSTANTPOTENTIALFORCE_H_
#define OPENMM_CPUCONSTANTPOTENTIALFORCE_H_

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

#include "CpuNeighborList.h"
#include "openmm/Kernel.h"
#include "openmm/kernels.h"
#include "openmm/internal/vectorize.h"
#include "tnt_array2d.h"
#include "jama_cholesky.h"

namespace OpenMM {

class CpuConstantPotentialForce;

/**
 * A generic charge solver for the constant potential method.
 */
class CpuConstantPotentialSolver {
public:
    /**
     * Creates a CpuConstantPotentialSolver.
     * 
     * @param numParticles           the number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     */
    CpuConstantPotentialSolver(int numParticles, int numElectrodeParticles);
    virtual ~CpuConstantPotentialSolver();
    /**
     * Mark precomputed data stored by the solver as invalid due to a change in
     * electrode parameters.
     */
    void invalidate();
    /**
     * Mark any existing solution stored by the solver as invalid due to a
     * change in system parameters.
     */
    void discardSavedSolution();
    /**
     * Solves for charges if needed.
     * 
     * @param conp       CPU constant potential force implementation
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    void solve(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel);
    /**
     * Solves for charges.
     * 
     * @param conp       CPU constant potential force implementation
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    virtual void solveImpl(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel) = 0;

protected:
    int numParticles, numElectrodeParticles;
    bool valid, hasSavedSolution;
    Vec3 savedBoxVectors[3];
    std::vector<Vec3> savedPositions;
    std::vector<float> savedCharges;
};

/**
 * A constant potential solver using direct inversion of the Coulomb matrix.
 * Suitable only when electrode particle positions are fixed.
 */
class CpuConstantPotentialMatrixSolver : public CpuConstantPotentialSolver {
private:
    Vec3 boxVectors[3];
    std::vector<Vec3> electrodePosData;
    JAMA::Cholesky<float> capacitance;
    std::vector<float> constraintVector;
    TNT::Array1D<float> b;

public:
    /**
     * Creates a CpuConstantPotentialMatrixSolver.
     * 
     * @param numParticles           the number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     */
    CpuConstantPotentialMatrixSolver(int numParticles, int numElectrodeParticles);
    /**
     * Solves for charges.
     * 
     * @param conp       CPU constant potential force implementation
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    void solveImpl(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel);

private:
    /**
     * Ensures that precomputed data stored by the solver is valid.
     * 
     * @param conp       CPU constant potential force implementation
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    void ensureValid(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel);
};

/**
 * A constant potential solver using the conjugate gradient method.  Suitable
 * for both fixed and variable electrode particle positions.
 */
class CpuConstantPotentialCGSolver : public CpuConstantPotentialSolver {
private:
    Vec3 boxVectors[3];
    bool precondRequested, precondActivated;
    std::vector<double> precondVector;
    std::vector<float> q;
    std::vector<float> grad;
    std::vector<float> projGrad;
    std::vector<float> precGrad;
    std::vector<float> qStep;
    std::vector<float> gradStep;
    std::vector<float> grad0;
    std::vector<float> qLast;

public:
    /**
     * Creates a CpuConstantPotentialCGSolver.
     * 
     * @param numParticles           the number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param precond                whether or not to use a preconditioner
     */
    CpuConstantPotentialCGSolver(int numParticles, int numElectrodeParticles, bool precond);
    /**
     * Solves for charges.
     * 
     * @param conp       CPU constant potential force implementation
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    void solveImpl(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel);

private:
    /**
     * Ensures that precomputed data stored by the solver is valid.
     * 
     * @param conp       CPU constant potential force implementation
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    void ensureValid(CpuConstantPotentialForce& conp, ThreadPool& threads, Kernel& pmeKernel);
};

/**
 * Performs energy, force, and charge derivative calculations for the CPU kernel
 * for ConstantPotentialForce.
 */
class CpuConstantPotentialForce {
    friend class CpuConstantPotentialSolver;
    friend class CpuConstantPotentialMatrixSolver;
    friend class CpuConstantPotentialCGSolver;

public:
    CpuConstantPotentialForce();
    virtual ~CpuConstantPotentialForce();

    /**
     * Initializes CpuConstantPotentialForce with data that cannot vary after
     * context creation, or pointers to data that can vary.
     * 
     * @param numParticles           the total number of particles
     * @param numElectrodeParticles  the number of electrode (fluctuating-charge) particles
     * @param posqIndex              index to identify charges loaded into the posq array
     * @param nonbondedCutoff        direct space cutoff
     * @param ewaldAlpha             Ewald reciprocal Gaussian width parameter
     * @param cgErrorTol             constant potential conjugate gradient error tolerance
     * @param gridSize               Ewald mesh dimensions
     * @param exceptionsArePeriodic  whether or not exceptions use periodic boundary conditions
     * @param useChargeConstraint    whether or not to constrain total charge
     * @param neighborList           neighbor list for direct space calculation
     * @param solver                 charge solver implementation
     * @param exclusions             particle exclusions
     * @param sysToElec              mapping from system particle indices to electrode particle indices
     * @param elecToSys              mapping from electrode particle indices to system particle indices
     * @param sysElec                mapping from system particle indices to electrode indices
     * @param elecElec               mapping from electrode particle indices to electrode indices
     * @param electrodeParams        electrode parameters
     * @param chargeTarget           target sum of charges on electrode particles only
     * @param externalField          electric field vector
     */
    void initialize(int numParticles, int numElectrodeParticles, int posqIndex,
        float nonbondedCutoff, float ewaldAlpha, float cgErrorTol,
        const int* gridSize, bool exceptionsArePeriodic,
        bool useChargeConstraint, const CpuNeighborList& neighborList,
        CpuConstantPotentialSolver* solver,
        const std::vector<std::set<int> >& exclusions,
        const std::vector<int>& sysToElec, const std::vector<int>& elecToSys,
        const std::vector<int>& sysElec, const std::vector<int>& elecElec,
        const std::vector<std::array<double, 3> >& electrodeParams,
        float chargeTarget, const float* externalField);
    /**
     * Updates CpuConstantPotentialForce after data has been modified.
     * 
     * @param chargeTarget    target sum of charges on electrode particles only
     * @param externalField   electric field vector
     * @param firstElectrode  index of the first electrode whose parameters may have been modified
     * @param lastElectrode   index of the last electrode whose parameters may have been modified
     */
    void update(float chargeTarget, const float* externalField,
        int firstElectrode, int lastElectrode);
    /**
     * Solves for charges and computes energies and forces.
     * 
     * @param boxVectors   periodic box vectors
     * @param posData      particle positions
     * @param charges      output particle charges
     * @param posq         wrapped single-precision position and charge array
     * @param threadForce  per-thread force accumulation arrays
     * @param energy       output system energy
     * @param threads      thread pool for parallel evaluation
     * @param pmeKernel    CPU PME solver kernel
     */
    void execute(const Vec3* boxVectors, const std::vector<Vec3>& posData,
        std::vector<float>& charges, float* posq,
        std::vector<AlignedArray<float> >& threadForce, double* energy,
        ThreadPool& threads, Kernel& pmeKernel);
    /**
     * Solves for charges without computing energies and forces.
     * 
     * @param boxVectors   periodic box vectors
     * @param posData      particle positions
     * @param charges      output particle charges
     * @param posq         wrapped single-precision position and charge array
     * @param threadForce  per-thread force accumulation arrays
     * @param threads      thread pool for parallel evaluation
     * @param pmeKernel    CPU PME solver kernel
     */
    void getCharges(const Vec3* boxVectors, const std::vector<Vec3>& posData,
        std::vector<float>& charges, float* posq,
        std::vector<AlignedArray<float> >& threadForce, ThreadPool& threads,
        Kernel& pmeKernel);

protected:
    /**
     * Precomputes values and updates pointers to data accessed by threads.
     * 
     * @param boxVectors   periodic box vectors
     * @param posData      particle positions
     * @param posq         wrapped single-precision position and charge array
     * @param threadForce  per-thread force accumulation arrays
     */
    void setThreadData(const Vec3* boxVectors, const std::vector<Vec3>& posData, float* posq, std::vector<AlignedArray<float> >& threadForce);
    /**
     * Extracts electrode particle charges from posq and places them in the
     * output charge array.
     * 
     * @param charges  output particle charges
     */
    void saveCharges(std::vector<float>& charges);
    /**
     * Computes energies and forces for fixed (solved) charges.
     * 
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     * @param energy     output system energy
     */
    void getEnergyForces(ThreadPool& threads, Kernel& pmeKernel, double* energy);
    /**
     * Thread worker for getEnergyForcesThread().
     * 
     * @param threads      thread pool for parallel evaluation
     * @param threadIndex  index of the thread in the thread pool
     * @param energy       whether or not to compute system energy
     */
    void getEnergyForcesThread(ThreadPool& threads, int threadIndex, bool includeEnergy);
    /**
     * Computes energy derivatives with respect to charges.
     * 
     * @param threads    thread pool for parallel evaluation
     * @param pmeKernel  CPU PME solver kernel
     */
    void getDerivatives(ThreadPool& threads, Kernel& pmeKernel);
    /**
     * Thread worker for getDerivatives().
     * 
     * @param threads      thread pool for parallel evaluation
     * @param threadIndex  index of the thread in the thread pool
     * @param plasmaTerm   precomputed neutralizing plasma term for derivatives
     */
    void getDerivativesThread(ThreadPool& threads, int threadIndex, float plasmaTerm);
    /**
     * Calculates a displacement between a pair of exception particles.
     * 
     * @param posI    the position of the first particle
     * @param posJ    the position of the second particle
     * @param deltaR  the displacement vector from particle I to particle J
     * @param r2      the square magnitude of the displacement vector
     */
    void getExceptionDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2) const;

    /**
     * Computes the direct space contribution to energies and forces for one
     * block of particles.
     * 
     * @param blockIndex  the index of the block
     * @param forces      output forces
     * @param energy      output energy
     */
    virtual void getEnergyForcesBlock(int blockIndex, float* forces, double* energy) = 0;
    /**
     * Computes the direct space contribution to charge derivatives for one
     * block of particles.
     * 
     * @param blockIndex   the index of the block
     * @param derivatives  output derivatives
     */
    virtual void getDerivativesBlock(int blockIndex, float* derivatives) = 0;
private:
    /**
     * Initializes energy and force lookup tables for a pair of electrodes.
     * 
     * @param ie  the index of the first electrode
     * @param je  the index of the second electrode
     */
    void initializeLookupTables(int ie, int je);

protected:
    static const int PotentialIndex;
    static const int GaussianWidthIndex;
    static const int ThomasFermiScaleIndex;

    static const int NUM_TABLE_POINTS;

    static const double TWO_OVER_SQRT_PI;
    static const double SELF_ALPHA_SCALE;
    static const double SELF_ETA_SCALE;
    static const double SELF_TF_SCALE;

    int numParticles, numElectrodeParticles, numElectrodes, posqIndex;
    float nonbondedCutoff, ewaldAlpha, cgErrorTol;
    int gridSize[3];
    bool exceptionsArePeriodic, useChargeConstraint;
    const CpuNeighborList* neighborList;
    CpuConstantPotentialSolver* solver;

    const std::set<int>* exclusions;
    const int* sysToElec;
    const int* elecToSys;
    const int* sysElec;
    const int* elecElec;
    const std::array<double, 3>* electrodeParams;

    std::vector<float> electrodePotentials;
    std::vector<float> electrodeSelfScales;

    float chargeTarget;
    fvec4 externalField;

    Vec3 boxVectors[3];
    AlignedArray<fvec4> boxVectorsVec4;
    fvec4 boxSize;
    fvec4 recipBoxSize;
    bool triclinic;
    float plasmaScale;

    const Vec3* posData;
    float* posq;
    std::vector<double> threadEnergy;
    std::vector<AlignedArray<float> >* threadForce;
    std::vector<float> chargeDerivatives;
    std::vector<float> energyLookupTable;
    std::vector<float> forceLookupTable;
    float tableScale;

    std::atomic<int> atomicBlockCounter, atomicParticleCounter;
};

class CpuConstantPotentialPmeIO : public CalcPmeReciprocalForceKernel::IO {
public:
    CpuConstantPotentialPmeIO(float* posq, float* force, float* chargeDerivatives, int numParticles, int numElectrodeParticles);
    float* getPosq();
    void setForce(float* f);
    void setChargeDerivatives(float* derivatives);
private:
    float* posq;
    float* force;
    float* chargeDerivatives;
    int numParticles;
    int numElectrodeParticles;
};

} // namespace OpenMM

#endif // OPENMM_CPUCONSTANTPOTENTIALFORCE_H_
