#ifndef OPENMM_OPENCLNONBONDEDUTILITIES_H_
#define OPENMM_OPENCLNONBONDEDUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
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

#include "OpenCLContext.h"
#include "openmm/System.h"
#include <string>
#include <vector>

namespace OpenMM {

class OpenCLCompact;

/**
 * This class implements features that are used by several different force.  It provides
 * a generic interface for calculating nonbonded interactions.
 */

class OpenCLNonbondedUtilities {
public:
    OpenCLNonbondedUtilities(OpenCLContext& context);
    ~OpenCLNonbondedUtilities();
    /**
     * Add a nonbonded interaction.
     *
     * @param usesCutoff     specifies whether a cutoff should be applied to this interaction
     * @param usesPeriodic   specifies whether periodic boundary conditions should be applied to this interaction
     * @param cutoffDistance the cutoff distance for this interaction (ignored if usesCutoff is false)
     * @param exclusionList  for each atom, specifies the list of other atoms whose interactions should be excluded
     */
    void addInteraction(bool usesCutoff, bool usesPeriodic, double cutoffDistance, const std::vector<std::vector<int> >& exclusionList);
    /**
     * Add a per-atom parameter that interactions may depend on.
     *
     * @param name      the name of the parameter
     * @param type      the data type of the parameter
     * @param size      the size of the parameter in bytes
     * @param buffer    the buffer containing the parameter values
     */
    void addParameter(const std::string& name, const std::string& type, int size, cl::Buffer& buffer);
    /**
     * Initialize this object in preparation for a simulation.
     */
    void initialize(const System& system);
    /**
     * Get the number of force buffers required for nonbonded forces.
     */
    int getNumForceBuffers() {
        return numForceBuffers;
    }
    /**
     * Get the periodic box size.
     */
    mm_float4 getPeriodicBoxSize() {
        return periodicBoxSize;
    }
    /**
     * Prepare to compute interactions.  This updates the neighbor list.
     */
    void prepareInteractions();
    /**
     * Compute the nonbonded interactions.  This will only be executed once after each call to
     * prepareInteractions().  Additional calls return immediately without doing anything.
     */
    void computeInteractions();
    /**
     * Get the array containing the center of each atom block.
     */
    OpenCLArray<mm_float4>& getBlockCenters() {
        return *blockCenter;
    }
    /**
     * Get the array containing the dimensions of each atom block.
     */
    OpenCLArray<mm_float4>& getBlockBoundingBoxes() {
        return *blockBoundingBox;
    }
    /**
     * Get the array containing the full set of tiles.
     */
    OpenCLArray<cl_uint>& getTiles() {
        return *tiles;
    }
    /**
     * Get the array whose first element contains the number of tiles with interactions.
     */
    OpenCLArray<cl_uint>& getInteractionCount() {
        return *interactionCount;
    }
    /**
     * Get the array containing tiles with interactions.
     */
    OpenCLArray<cl_uint>& getInteractingTiles() {
        return *interactingTiles;
    }
    /**
     * Get the array containing flags for tiles with interactions.
     */
    OpenCLArray<cl_uint>& getInteractionFlags() {
        return *interactionFlags;
    }
private:
    class ParameterInfo;
    OpenCLContext& context;
    cl::Kernel forceKernel;
    cl::Kernel findBlockBoundsKernel;
    cl::Kernel findInteractingBlocksKernel;
    cl::Kernel findInteractionsWithinBlocksKernel;
    OpenCLArray<cl_uint>* tiles;
    OpenCLArray<cl_uint>* exclusionIndex;
    OpenCLArray<cl_uint>* exclusions;
    OpenCLArray<cl_uint>* interactingTiles;
    OpenCLArray<cl_uint>* interactionFlags;
    OpenCLArray<cl_uint>* interactionCount;
    OpenCLArray<mm_float4>* blockCenter;
    OpenCLArray<mm_float4>* blockBoundingBox;
    std::vector<std::vector<int> > atomExclusions;
    std::vector<ParameterInfo> parameters;
    OpenCLCompact* compact;
    double cutoff;
    bool useCutoff, usePeriodic, forceBufferPerAtomBlock, hasComputedInteractions;
    int numForceBuffers;
    mm_float4 periodicBoxSize;
};

class OpenCLNonbondedUtilities::ParameterInfo {
public:
    ParameterInfo(const std::string& name, const std::string& type, int size, cl::Buffer& buffer) :
            name(name), type(type), size(size), buffer(&buffer) {
    }
    std::string name;
    std::string type;
    int size;
    cl::Buffer* buffer;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLNONBONDEDUTILITIES_H_*/
