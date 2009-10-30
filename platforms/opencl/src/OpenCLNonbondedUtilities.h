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
 * This class provides a generic interface for calculating nonbonded interactions.  It does this in two
 * ways.  First, it can be used to create Kernels that evaluate nonbonded interactions.  Clients
 * only need to provide the code for evaluating a single interaction and the list of parameters it depends on.
 * A complete kernel is then synthesized using an appropriate algorithm to evaluate all interactions on all
 * atoms.
 *
 * Second, this class itself creates and invokes a single "default" interaction kernel, allowing several
 * different forces to be evaluated at once for greater efficiency.  Call addInteraction() and addParameter()
 * to add interactions to this default kernel.
 *
 * During each force or energy evaluation, the following sequence of steps takes place:
 *
 * 1. Data structures (e.g. neighbor lists) are calculated to allow nonbonded interactions to be evaluated
 * quickly.
 *
 * 2. calcForces() or calcEnergy() is called on each ForceImpl in the System.
 *
 * 3. Finally, the default interaction kernel is invoked to calculate all interactions that were added
 * to it.
 *
 * This sequence means that the default interaction kernel may depend on quantities that were calculated
 * by ForceImpls during calcForces() or calcEnergy().
 */

class OpenCLNonbondedUtilities {
public:
    class ParameterInfo;
    OpenCLNonbondedUtilities(OpenCLContext& context);
    ~OpenCLNonbondedUtilities();
    /**
     * Add a nonbonded interaction to be evaluated by the default interaction kernel.
     *
     * @param usesCutoff     specifies whether a cutoff should be applied to this interaction
     * @param usesPeriodic   specifies whether periodic boundary conditions should be applied to this interaction
     * @param usesExclusions specifies whether this interaction uses exclusions.  If this is true, it must have identical exclusions to every other interaction.
     * @param cutoffDistance the cutoff distance for this interaction (ignored if usesCutoff is false)
     * @param exclusionList  for each atom, specifies the list of other atoms whose interactions should be excluded
     * @param kernel         the code to evaluate the interaction
     */
    void addInteraction(bool usesCutoff, bool usesPeriodic, bool usesExclusions, double cutoffDistance, const std::vector<std::vector<int> >& exclusionList, const std::string& kernel);
    /**
     * Add a per-atom parameter that the default interaction kernel may depend on.
     */
    void addParameter(const ParameterInfo& parameter);
    /**
     * Add an array (other than a per-atom parameter) that should be passed as an argument to the default interaction kernel.
     */
    void addArgument(const ParameterInfo& parameter);
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
     * Get whether a cutoff is being used.
     */
    bool getUseCutoff() {
        return useCutoff;
    }
    /**
     * Get whether periodic boundary conditions are being used.
     */
    bool getUsePeriodic() {
        return usePeriodic;
    }
    /**
     * Get whether there is one force buffer per atom block.
     */
    bool getForceBufferPerAtomBlock() {
        return forceBufferPerAtomBlock;
    }
    /**
     * Get the cutoff distance.
     */
    double getCutoffDistance() {
        return cutoff;
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
    /**
     * Create a Kernel for evaluating a nonbonded interaction.  Cutoffs and periodic boundary conditions
     * are assumed to be the same as those for the default interaction Kernel, since this kernel will use
     * the same neighbor list.
     * 
     * @param source        the source code for evaluating the force and energy
     * @param params        the per-atom parameters this kernel may depend on
     * @param arguments     arrays (other than per-atom parameters) that should be passed as arguments to the kernel
     * @param useExclusions specifies whether exclusions are applied to this interaction
     */
    cl::Kernel createInteractionKernel(const std::string& source, const std::vector<ParameterInfo>& params, const std::vector<ParameterInfo>& arguments, bool useExclusions) const;
private:
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
    std::vector<ParameterInfo> arguments;
    OpenCLCompact* compact;
    std::string kernelSource;
    std::map<std::string, std::string> kernelDefines;
    double cutoff;
    bool useCutoff, usePeriodic, forceBufferPerAtomBlock;
    int numForceBuffers;
    mm_float4 periodicBoxSize;
};

/**
 * This class stores information about a per-atom parameter that may be used in a nonbonded kernel.
 */

class OpenCLNonbondedUtilities::ParameterInfo {
public:
    /**
     * Create a ParameterInfo object.
     *
     * @param name      the name of the parameter
     * @param type      the data type of the parameter
     * @param size      the size of the parameter in bytes
     * @param buffer    the buffer containing the parameter values
     */
    ParameterInfo(const std::string& name, const std::string& type, int size, cl::Buffer& buffer) :
            name(name), type(type), size(size), buffer(&buffer) {
    }
    const std::string& getName() const {
        return name;
    }
    const std::string& getType() const {
        return type;
    }
    int getSize() const {
        return size;
    }
    cl::Buffer& getBuffer() const {
        return *buffer;
    }
private:
    std::string name;
    std::string type;
    int size;
    cl::Buffer* buffer;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLNONBONDEDUTILITIES_H_*/
