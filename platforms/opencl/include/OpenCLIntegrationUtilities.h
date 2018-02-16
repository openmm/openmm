#ifndef OPENMM_OPENCLINTEGRATIONUTILITIES_H_
#define OPENMM_OPENCLINTEGRATIONUTILITIES_H_

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
#include "OpenCLContext.h"
#include "windowsExportOpenCL.h"
#include <iosfwd>

namespace OpenMM {

/**
 * This class implements features that are used by many different integrators, including
 * common workspace arrays, random number generation, and enforcing constraints.
 */

class OPENMM_EXPORT_OPENCL OpenCLIntegrationUtilities {
public:
    OpenCLIntegrationUtilities(OpenCLContext& context, const System& system);
    /**
     * Get the array which contains position deltas.
     */
    OpenCLArray& getPosDelta() {
        return posDelta;
    }
    /**
     * Get the array which contains random values.  Each element is a float4, whose components
     * are independent, normally distributed random numbers with mean 0 and variance 1.
     */
    OpenCLArray& getRandom() {
        return random;
    }
    /**
     * Get the array which contains the current step size.
     */
    OpenCLArray& getStepSize() {
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
    OpenCLContext& context;
    cl::Kernel settlePosKernel, settleVelKernel;
    cl::Kernel shakePosKernel, shakeVelKernel;
    cl::Kernel ccmaDirectionsKernel;
    cl::Kernel ccmaPosForceKernel, ccmaVelForceKernel;
    cl::Kernel ccmaMultiplyKernel;
    cl::Kernel ccmaPosUpdateKernel, ccmaVelUpdateKernel;
    cl::Kernel vsitePositionKernel, vsiteForceKernel, vsiteAddForcesKernel;
    cl::Kernel randomKernel, timeShiftKernel;
    OpenCLArray posDelta;
    OpenCLArray settleAtoms;
    OpenCLArray settleParams;
    OpenCLArray shakeAtoms;
    OpenCLArray shakeParams;
    OpenCLArray random;
    OpenCLArray randomSeed;
    OpenCLArray stepSize;
    OpenCLArray ccmaAtoms;
    OpenCLArray ccmaDistance;
    OpenCLArray ccmaReducedMass;
    OpenCLArray ccmaAtomConstraints;
    OpenCLArray ccmaNumAtomConstraints;
    OpenCLArray ccmaConstraintMatrixColumn;
    OpenCLArray ccmaConstraintMatrixValue;
    OpenCLArray ccmaDelta1;
    OpenCLArray ccmaDelta2;
    OpenCLArray ccmaConverged;
    OpenCLArray ccmaConvergedHostBuffer;
    OpenCLArray vsite2AvgAtoms;
    OpenCLArray vsite2AvgWeights;
    OpenCLArray vsite3AvgAtoms;
    OpenCLArray vsite3AvgWeights;
    OpenCLArray vsiteOutOfPlaneAtoms;
    OpenCLArray vsiteOutOfPlaneWeights;
    OpenCLArray vsiteLocalCoordsIndex;
    OpenCLArray vsiteLocalCoordsAtoms;
    OpenCLArray vsiteLocalCoordsWeights;
    OpenCLArray vsiteLocalCoordsPos;
    OpenCLArray vsiteLocalCoordsStartIndex;
    int randomPos;
    int lastSeed, numVsites;
    bool hasInitializedPosConstraintKernels, hasInitializedVelConstraintKernels, ccmaUseDirectBuffer, hasOverlappingVsites;
    mm_double2 lastStepSize;
    struct ShakeCluster;
    struct ConstraintOrderer;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLINTEGRATIONUTILITIES_H_*/
