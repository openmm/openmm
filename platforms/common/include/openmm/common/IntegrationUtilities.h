#ifndef OPENMM_INTEGRATIONUTILITIES_H_
#define OPENMM_INTEGRATIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2022 Stanford University and the Authors.      *
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

#include "openmm/common/ComputeArray.h"
#include "openmm/common/ComputeKernel.h"
#include "openmm/common/ComputeVectorTypes.h"
#include "openmm/System.h"
#include <iosfwd>
#include <map>

namespace OpenMM {

class ComputeContext;

/**
 * This class implements features that are used by many different integrators, including
 * common workspace arrays, random number generation, and enforcing constraints.
 */

class OPENMM_EXPORT_COMMON IntegrationUtilities {
public:
    IntegrationUtilities(ComputeContext& context, const System& system);
    virtual ~IntegrationUtilities() {
    }
    /**
     * Get the array which contains position deltas.  These are the amounts by
     * which the position of each atom will change in the current step.  The actual
     * positions should not be modified until after constraints have been applied.
     */
    virtual ArrayInterface& getPosDelta() = 0;
    /**
     * Get the array which contains random values.  Each element is a float4 whose components
     * are independent, normally distributed random numbers with mean 0 and variance 1.
     * Be sure to call initRandomNumberGenerator() and prepareRandomNumbers() before
     * accessing this array.
     */
    virtual ArrayInterface& getRandom() = 0;
    /**
     * Get the array which contains the current step size.
     */
    virtual ArrayInterface& getStepSize() = 0;
    /**
     * Set the size to use for the next step.
     */
    void setNextStepSize(double size);
    /**
     * Get the size that was used for the last step.
     */
    double getLastStepSize();
    /**
     * Apply constraints to the atom positions.  When calling this method, the
     * context's array of positions should contain the positions at the start of the
     * step, and the array returned by getPosDelta() should contain the intended
     * change to each position.  This method modifies the position deltas so that,
     * once they are added to the positions, constraints will be satisfied.
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
     * Initialize the random number generator.  This should be called once when the
     * context is first created.  Subsequent calls will be ignored if the random
     * seed is the same as on the first call, or throw an exception if the random
     * seed is different.
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
    virtual void distributeForcesFromVirtualSites() = 0;
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
    /**
     * Compute the current velocities, shifting them in time to account for a leapfrog integrator.
     * 
     * @param timeShift   the amount by which to shift the velocities in time
     * @param velocities  the shifted velocities are returned in this
     */
    void computeShiftedVelocities(double timeShift, std::vector<Vec3>& velocities);
protected:
    virtual void applyConstraintsImpl(bool constrainVelocities, double tol) = 0;
    ComputeContext& context;
    ComputeKernel settlePosKernel, settleVelKernel;
    ComputeKernel shakePosKernel, shakeVelKernel;
    ComputeKernel ccmaDirectionsKernel, ccmaPosForceKernel, ccmaVelForceKernel;
    ComputeKernel ccmaMultiplyKernel, ccmaUpdateKernel, ccmaFullKernel;
    ComputeKernel vsitePositionKernel, vsiteForceKernel, vsiteSaveForcesKernel;
    ComputeKernel randomKernel, timeShiftKernel;
    ComputeArray posDelta;
    ComputeArray settleAtoms;
    ComputeArray settleParams;
    ComputeArray shakeAtoms;
    ComputeArray shakeParams;
    ComputeArray random;
    ComputeArray randomSeed;
    ComputeArray stepSize;
    ComputeArray ccmaAtoms;
    ComputeArray ccmaConstraintAtoms;
    ComputeArray ccmaDistance;
    ComputeArray ccmaReducedMass;
    ComputeArray ccmaAtomConstraints;
    ComputeArray ccmaNumAtomConstraints;
    ComputeArray ccmaConstraintMatrixColumn;
    ComputeArray ccmaConstraintMatrixValue;
    ComputeArray ccmaDelta1;
    ComputeArray ccmaDelta2;
    ComputeArray ccmaConverged;
    ComputeArray vsite2AvgAtoms;
    ComputeArray vsite2AvgWeights;
    ComputeArray vsite3AvgAtoms;
    ComputeArray vsite3AvgWeights;
    ComputeArray vsiteOutOfPlaneAtoms;
    ComputeArray vsiteOutOfPlaneWeights;
    ComputeArray vsiteLocalCoordsIndex;
    ComputeArray vsiteLocalCoordsAtoms;
    ComputeArray vsiteLocalCoordsWeights;
    ComputeArray vsiteLocalCoordsPos;
    ComputeArray vsiteLocalCoordsStartIndex;
    ComputeArray vsiteStage;
    int randomPos, lastSeed, numVsites, numVsiteStages;
    bool hasOverlappingVsites;
    mm_double2 lastStepSize;
    struct ShakeCluster;
    struct ConstraintOrderer;
};

} // namespace OpenMM

#endif /*OPENMM_INTEGRATIONUTILITIES_H_*/
