#ifndef OPENMM_NONBONDEDUTILITIES_H_
#define OPENMM_NONBONDEDUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
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

#include "openmm/common/ArrayInterface.h"
#include "openmm/common/ComputeParameterInfo.h"
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This class provides a generic interface for calculating nonbonded interactions.  Clients only need
 * to provide the code for evaluating a single interaction and the list of parameters it depends on.
 * A complete kernel is then synthesized using an appropriate algorithm to evaluate all interactions on
 * all atoms.  Call addInteraction() to define a nonbonded interaction, and addParameter() to define
 * per-particle parameters that the interaction depends on.
 *
 * During each force or energy evaluation, the following sequence of steps takes place:
 *
 * 1. Data structures (e.g. neighbor lists) are calculated to allow nonbonded interactions to be evaluated
 * quickly.
 *
 * 2. calcForcesAndEnergy() is called on each ForceImpl in the System.
 *
 * 3. Finally, the default interaction kernel is invoked to calculate all interactions that were added
 * to it.
 *
 * This sequence means that the default interaction kernel may depend on quantities that were calculated
 * by ForceImpls during calcForcesAndEnergy().
 */

class OPENMM_EXPORT_COMMON NonbondedUtilities {
public:
    virtual ~NonbondedUtilities() {
    }
    /**
     * Add a nonbonded interaction to be evaluated by the default interaction kernel.
     *
     * @param usesCutoff     specifies whether a cutoff should be applied to this interaction
     * @param usesPeriodic   specifies whether periodic boundary conditions should be applied to this interaction
     * @param usesExclusions specifies whether this interaction uses exclusions.  If this is true, it must have identical exclusions to every other interaction.
     * @param cutoffDistance the cutoff distance for this interaction (ignored if usesCutoff is false)
     * @param exclusionList  for each atom, specifies the list of other atoms whose interactions should be excluded
     * @param kernel         the code to evaluate the interaction
     * @param forceGroup     the force group in which the interaction should be calculated
     */
    virtual void addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance, const std::vector<std::vector<int> >& exclusionList, const std::string& kernel, int forceGroup) = 0;
    /**
     * Add a per-atom parameter that the default interaction kernel may depend on.
     */
    virtual void addParameter(ComputeParameterInfo parameter) = 0;
    /**
     * Add an array (other than a per-atom parameter) that should be passed as an argument to the default interaction kernel.
     */
    virtual void addArgument(ComputeParameterInfo parameter) = 0;
    /**
     * Register that the interaction kernel will be computing the derivative of the potential energy
     * with respect to a parameter.
     * 
     * @param param   the name of the parameter
     * @return the variable that will be used to accumulate the derivative.  Any code you pass to addInteraction() should
     * add its contributions to this variable.
     */
    virtual std::string addEnergyParameterDerivative(const std::string& param) = 0;
    /**
     * Get the number of force buffers required for nonbonded forces.
     */
    virtual int getNumForceBuffers() const = 0;
    /**
     * Get whether a cutoff is being used.
     */
    virtual bool getUseCutoff() = 0;
    /**
     * Get whether periodic boundary conditions are being used.
     */
    virtual bool getUsePeriodic() = 0;
    /**
     * Get the number of thread blocks used for computing nonbonded forces.
     */
    virtual int getNumForceThreadBlocks() = 0;
    /**
     * Get the size of each thread block used for computing nonbonded forces.
     */
    virtual int getForceThreadBlockSize() = 0;
    /**
     * Get the maximum cutoff distance used by any interaction.
     */
    virtual double getMaxCutoffDistance() = 0;
    /**
     * Given a nonbonded cutoff, get the padded cutoff distance used in computing
     * the neighbor list.
     */
    virtual double padCutoff(double cutoff) = 0;
    /**
     * Get the array containing the center of each atom block.
     */
    virtual ArrayInterface& getBlockCenters() = 0;
    /**
     * Get the array containing the dimensions of each atom block.
     */
    virtual ArrayInterface& getBlockBoundingBoxes() = 0;
    /**
     * Get the array whose first element contains the number of tiles with interactions.
     */
    virtual ArrayInterface& getInteractionCount() = 0;
    /**
     * Get the array containing tiles with interactions.
     */
    virtual ArrayInterface& getInteractingTiles() = 0;
    /**
     * Get the array containing the atoms in each tile with interactions.
     */
    virtual ArrayInterface& getInteractingAtoms() = 0;
    /**
     * Get the array containing exclusion flags.
     */
    virtual ArrayInterface& getExclusions() = 0;
    /**
     * Get the array containing tiles with exclusions.
     */
    virtual ArrayInterface& getExclusionTiles() = 0;
    /**
     * Get the array containing the index into the exclusion array for each tile.
     */
    virtual ArrayInterface& getExclusionIndices() = 0;
    /**
     * Get the array listing where the exclusion data starts for each row.
     */
    virtual ArrayInterface& getExclusionRowIndices() = 0;
    /**
     * Get the array containing a flag for whether the neighbor list was rebuilt
     * on the most recent call to prepareInteractions().
     */
    virtual ArrayInterface& getRebuildNeighborList() = 0;
    /**
     * Get the index of the first tile this context is responsible for processing.
     */
    virtual int getStartTileIndex() const = 0;
    /**
     * Get the total number of tiles this context is responsible for processing.
     */
    virtual int getNumTiles() const = 0;
    /**
     * Set whether to add padding to the cutoff distance when building the neighbor list.
     * This increases the size of the neighbor list (and thus the cost of computing interactions),
     * but also means we don't need to rebuild it every time step.  The default value is true,
     * since usually this improves performance.  For very expensive interactions, however,
     * it may be better to set this to false.
     */
    virtual void setUsePadding(bool padding) = 0;
    /**
     * Initialize this object in preparation for a simulation.
     */
    virtual void initialize(const System& system) = 0;
    /**
     * Prepare to compute interactions.  This updates the neighbor list.
     */
    virtual void prepareInteractions(int forceGroups) = 0;
    /**
     * Compute the nonbonded interactions.
     * 
     * @param forceGroups    the flags specifying which force groups to include
     * @param includeForces  whether to compute forces
     * @param includeEnergy  whether to compute the potential energy
     */
    virtual void computeInteractions(int forceGroups, bool includeForces, bool includeEnergy) = 0;
    /**
     * Set the source code for the main kernel.  It only needs to be changed in very unusual circumstances.
     */
    virtual void setKernelSource(const std::string& source) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_NONBONDEDUTILITIES_H_*/
