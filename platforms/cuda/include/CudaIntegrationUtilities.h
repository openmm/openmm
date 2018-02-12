#ifndef OPENMM_CUDAINTEGRATIONUTILITIES_H_
#define OPENMM_CUDAINTEGRATIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2018 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/System.h"
#include "CudaContext.h"
#include "windowsExportCuda.h"
#include <iosfwd>

namespace OpenMM {

/**
 * This class implements features that are used by many different integrators, including
 * common workspace arrays, random number generation, and enforcing constraints.
 */

class OPENMM_EXPORT_CUDA CudaIntegrationUtilities {
public:
    CudaIntegrationUtilities(CudaContext& context, const System& system);
    ~CudaIntegrationUtilities();
    /**
     * Get the array which contains position deltas.
     */
    CudaArray& getPosDelta() {
        return posDelta;
    }
    /**
     * Get the array which contains random values.  Each element is a float4, whose components
     * are independent, normally distributed random numbers with mean 0 and variance 1.
     */
    CudaArray& getRandom() {
        return random;
    }
    /**
     * Get the array which contains the current step size.
     */
    CudaArray& getStepSize() {
        return stepSize;
    }
    /**
     * Set the size to use for the next step.
     */
    void setNextStepSize(double size);
    /**
     * Get the size that was used for the last step.
     */
    double getLastStepSize();
    /**
     * Apply constraints to the atom positions.
     *
     * @param tol             the constraint tolerance
     */
    void applyConstraints(double tol);
    /**
     * Apply constraints to the atom velocities.
     *
     * @param tol             the constraint tolerance
     */
    void applyVelocityConstraints(double tol);
    /**
     * Initialize the random number generator.
     */
    void initRandomNumberGenerator(unsigned int randomNumberSeed);
    /**
     * Ensure that sufficient random numbers are available in the array, and generate new ones if not.
     *
     * @param numValues     the number of random float4's that will be required
     * @return the index in the array at which to start reading
     */
    int prepareRandomNumbers(int numValues);
    /**
     * Compute the positions of virtual sites.
     */
    void computeVirtualSites();
    /**
     * Distribute forces from virtual sites to the atoms they are based on.
     */
    void distributeForcesFromVirtualSites();
    /**
     * Create a checkpoint recording the current state of the random number generator.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(std::istream& stream);
    /**
     * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
     * for a leapfrog integrator.
     * 
     * @param timeShift   the amount by which to shift the velocities in time
     */
    double computeKineticEnergy(double timeShift);
private:
    void applyConstraints(bool constrainVelocities, double tol);
    CudaContext& context;
    CUfunction settlePosKernel, settleVelKernel;
    CUfunction shakePosKernel, shakeVelKernel;
    CUfunction ccmaDirectionsKernel;
    CUfunction ccmaPosForceKernel, ccmaVelForceKernel;
    CUfunction ccmaMultiplyKernel;
    CUfunction ccmaUpdateKernel;
    CUfunction vsitePositionKernel, vsiteForceKernel;
    CUfunction randomKernel, timeShiftKernel;
    CudaArray posDelta;
    CudaArray settleAtoms;
    CudaArray settleParams;
    CudaArray shakeAtoms;
    CudaArray shakeParams;
    CudaArray random;
    CudaArray randomSeed;
    CudaArray stepSize;
    CudaArray ccmaAtoms;
    CudaArray ccmaDistance;
    CudaArray ccmaReducedMass;
    CudaArray ccmaAtomConstraints;
    CudaArray ccmaNumAtomConstraints;
    CudaArray ccmaConstraintMatrixColumn;
    CudaArray ccmaConstraintMatrixValue;
    CudaArray ccmaDelta1;
    CudaArray ccmaDelta2;
    CudaArray ccmaConverged;
    int* ccmaConvergedMemory;
    CUdeviceptr ccmaConvergedDeviceMemory;
    CUevent ccmaEvent;
    CudaArray vsite2AvgAtoms;
    CudaArray vsite2AvgWeights;
    CudaArray vsite3AvgAtoms;
    CudaArray vsite3AvgWeights;
    CudaArray vsiteOutOfPlaneAtoms;
    CudaArray vsiteOutOfPlaneWeights;
    CudaArray vsiteLocalCoordsIndex;
    CudaArray vsiteLocalCoordsAtoms;
    CudaArray vsiteLocalCoordsWeights;
    CudaArray vsiteLocalCoordsPos;
    CudaArray vsiteLocalCoordsStartIndex;
    int randomPos;
    int lastSeed, numVsites;
    double2 lastStepSize;
    struct ShakeCluster;
    struct ConstraintOrderer;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAINTEGRATIONUTILITIES_H_*/
