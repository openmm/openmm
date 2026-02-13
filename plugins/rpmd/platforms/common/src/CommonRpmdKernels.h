#ifndef COMMON_RPMD_KERNELS_H_
#define COMMON_RPMD_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2011-2021 Stanford University and the Authors.      *
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

#include "openmm/RpmdKernels.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeArray.h"
#include <map>

namespace OpenMM {

/**
 * This kernel is invoked by RPMDIntegrator to take one time step, and to get and
 * set the state of system copies.
 */
class CommonIntegrateRPMDStepKernel : public IntegrateRPMDStepKernel {
public:
    CommonIntegrateRPMDStepKernel(const std::string& name, const Platform& platform, ComputeContext& cc) :
            IntegrateRPMDStepKernel(name, platform), cc(cc), hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the RPMDIntegrator this kernel will be used for
     */
    void initialize(const System& system, const RPMDIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     * @param forcesAreValid if the context has been modified since the last time step, this will be
     *                       false to show that cached forces are invalid and must be recalculated
     */
    void execute(ContextImpl& context, const RPMDIntegrator& integrator, bool forcesAreValid);
    /**
     * Compute the kinetic energy.
     * 
     * @param context        the context in which to execute this kernel
     * @param integrator     the RPMDIntegrator this kernel is being used for
     */
    double computeKineticEnergy(ContextImpl& context, const RPMDIntegrator& integrator);
    /**
     * Set the positions of all particles in one copy of the system.
     */
    void setPositions(int copy, const std::vector<Vec3>& positions);
    /**
     * Set the velocities of all particles in one copy of the system.
     */
    void setVelocities(int copy, const std::vector<Vec3>& velocities);
    /**
     * Copy positions and velocities for one copy into the context.
     */
    void copyToContext(int copy, ContextImpl& context);
private:
    void initializeKernels(ContextImpl& context);
    void computeForces(ContextImpl& context);
    /**
     * Download bead positions from GPU for batched force evaluation.
     * 
     * @param beadIndex    which bead to download (0 to numCopies-1)
     * @param outPositions output vector to fill with positions
     */
    void downloadPositionsFromGPU(int beadIndex, std::vector<Vec3>& outPositions);
    /**
     * Upload bead forces to GPU after batched force evaluation.
     * 
     * @param beadIndex which bead to upload (0 to numCopies-1)
     * @param inForces  input vector with forces
     */
    void uploadForcesToGPU(int beadIndex, const std::vector<Vec3>& inForces);
    /**
     * Upload all bead forces to GPU at once after batched force evaluation.
     * This is more efficient than calling uploadForcesToGPU multiple times.
     * 
     * @param allBeadForces vector of force vectors for all beads [numCopies][numParticles]
     */
    void uploadAllForcesToGPU(const std::vector<std::vector<Vec3>>& allBeadForces);
    /**
     * Apply the Bussi stochastic velocity rescaling thermostat to the centroid mode.
     * This is used for PILE_G mode where Bussi thermostat is applied to centroid only.
     */
    void applyBussiCentroidThermostat(const System& system, const RPMDIntegrator& integrator, double halfdt);
    /**
     * Apply the Bussi stochastic velocity rescaling thermostat to classical particles.
     * This is used in hybrid mode for classical particle thermostating.
     */
    void applyBussiClassicalThermostat(const System& system, const RPMDIntegrator& integrator, double halfdt);
    std::string createFFT(int size, const std::string& variable, bool forward);
    ComputeContext& cc;
    bool hasInitializedKernels;
    int numCopies, numParticles, workgroupSize;
    std::map<int, int> groupsByCopies;
    int groupsNotContracted;
    ComputeArray forces;
    ComputeArray positions;
    ComputeArray velocities;
    ComputeArray contractedForces;
    ComputeArray contractedPositions;
    ComputeArray centroidKE;
    ComputeKernel pileKernel, stepKernel, velocitiesKernel, copyToContextKernel, copyFromContextKernel, translateKernel;
    ComputeKernel computeCentroidKEKernel, applyBussiScalingKernel;
    std::map<int, ComputeKernel> positionContractionKernels;
    std::map<int, ComputeKernel> forceContractionKernels;
    // Hybrid mode support: quantum vs classical particles (uniform memory layout)
    // All particles stored with all beads; classical beads synced to stay identical
    int numQuantumParticles, numClassicalParticles;
    bool hybridMode;  // True if any classical particles exist
    ComputeArray isQuantum;         // Per-particle: 1=quantum, 0=classical
    std::vector<int> quantumParticleIndices;   // Original indices of quantum particles
    std::vector<int> classicalParticleIndices; // Original indices of classical particles
    
    // Kernels for hybrid mode
    ComputeArray classicalKE;  // For Bussi thermostat on classical particles
    ComputeKernel pileKernelHybrid;  // PILE thermostat for quantum particles only
    ComputeKernel stepKernelHybrid;  // Integration step (quantum: FFT, classical: Verlet)
    ComputeKernel velocitiesKernelHybrid;  // Velocity update for all particles
    ComputeKernel syncClassicalBeadsKernel;  // Sync bead 0 to beads 1..N-1 for classical particles
    ComputeKernel applyClassicalThermostatKernel;  // Langevin for classical particles
    ComputeKernel computeClassicalKEKernel;  // KE reduction for Bussi thermostat
    ComputeKernel applyBussiClassicalScalingKernel;  // Bussi velocity scaling
};

} // namespace OpenMM

#endif /*COMMON_RPMD_KERNELS_H_*/
