#ifndef OPENMM_COMPUTECONTEXT_H_
#define OPENMM_COMPUTECONTEXT_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019 Stanford University and the Authors.           *
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

#include "openmm/common/BondedUtilities.h"
#include "openmm/common/ComputeForceInfo.h"
#include "openmm/common/ComputeProgram.h"
#include "openmm/common/ComputeVectorTypes.h"
#include <map>
#include <string>

namespace OpenMM {

class ArrayInterface;

/**
 * This abstract class defines the interface by which platforms compile and execute
 * kernels.  It also manages the arrays use for storing standard information, like
 * positions and forces.
 */

class ComputeContext {
public:
    virtual ~ComputeContext() {
    }
    /**
     * Add a ComputeForceInfo to this context.  Force kernels call this during initialization
     * to provide information about particular forces.
     */
    virtual void addForce(ComputeForceInfo* force) = 0;
    /**
     * Set this as the current context for the calling thread.  This should be called before
     * doing any computation when you do not know what other code has just been executing on
     * the thread.  Platforms that rely on binding contexts to threads (such as CUDA) need to
     * implement this.
     */
    virtual void setAsCurrent() {
    }
    /**
     * Get the number of contexts being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one ComputeContext is created for each device.
     */
    virtual int getNumContexts() const = 0;
    /**
     * Get the index of this context in the list of ones being used for the current simulation.
     * This is relevant when a simulation is parallelized across multiple devices.  In that case,
     * one ComputeContext is created for each device.
     */
    virtual int getContextIndex() const = 0;
    /**
     * Construct an uninitialized array of the appropriate class for this platform.  The returned
     * value should be created on the heap with the "new" operator.
     */
    virtual ArrayInterface* createArray() = 0;
    /**
     * Compile source code to create a ComputeProgram.
     *
     * @param source             the source code of the program
     * @param defines            a set of preprocessor definitions (name, value) to define when compiling the program
     */
    virtual ComputeProgram compileProgram(const std::string source, const std::map<std::string, std::string>& defines=std::map<std::string, std::string>()) = 0;
    /**
     * Get the SIMD width of the device being used.
     */
    virtual int getSIMDWidth() const = 0;
    /**
     * Get whether the device being used supports 64 bit atomic operations on global memory.
     */
    virtual bool getSupports64BitGlobalAtomics() const = 0;
    /**
     * Get whether the device being used supports double precision math.
     */
    virtual bool getSupportsDoublePrecision() const = 0;
    /**
     * Get whether double precision is being used.
     */
    virtual bool getUseDoublePrecision() const = 0;
    /**
     * Get whether mixed precision is being used.
     */
    virtual bool getUseMixedPrecision() const = 0;
    /**
     * Get the number of atoms in the system.
     */
    virtual int getNumAtoms() const = 0;
    /**
     * Get the number of atoms, rounded up to a multiple of 32.  This is the actual size of
     * most arrays with one element per atom.
     */
    virtual int getPaddedNumAtoms() const = 0;
    /**
     * Get the array which contains the position (the xyz components) and charge (the w component) of each atom.
     */
    virtual ArrayInterface& getPosq() = 0;
    /**
     * Get the array which contains a correction to the position of each atom.  This only exists if getUseMixedPrecision() returns true.
     */
    virtual ArrayInterface& getPosqCorrection() = 0;
    /**
     * Get the array which contains the velocity (the xyz components) and inverse mass (the w component) of each atom.
     */
    virtual ArrayInterface& getVelm() = 0;
    /**
     * Get the array which contains a contribution to each force represented as 64 bit fixed point.
     */
    virtual ArrayInterface& getLongForceBuffer() = 0;
    /**
     * Get the array which contains the buffer in which energy is computed.
     */
    virtual ArrayInterface& getEnergyBuffer() = 0;
    /**
     * Get the array which contains the buffer in which derivatives of the energy with respect to parameters are computed.
     */
    virtual ArrayInterface& getEnergyParamDerivBuffer() = 0;
    /**
     * Replace all occurrences of a list of substrings.
     *
     * @param input   a string to process
     * @param replacements a set of strings that should be replaced with new strings wherever they appear in the input string
     * @return a new string produced by performing the replacements
     */
    std::string replaceStrings(const std::string& input, const std::map<std::string, std::string>& replacements) const;
    /**
     * Convert a number to a string in a format suitable for including in a kernel.
     * This takes into account whether the context uses single or double precision.
     */
    std::string doubleToString(double value) const;
    /**
     * Convert a number to a string in a format suitable for including in a kernel.
     */
    std::string intToString(int value) const;
    /**
     * Get whether the periodic box is triclinic.
     */
    virtual bool getBoxIsTriclinic() const = 0;
    /**
     * Get the BondedUtilities for this context.
     */
    virtual BondedUtilities& getBondedUtilities() = 0;
    /**
     * Mark that the current molecule definitions (and hence the atom order) may be invalid.
     * This should be called whenever force field parameters change.  It will cause the definitions
     * and order to be revalidated.
     * 
     * If you know which force has changed, calling the alternate form that takes a ComputeForceInfo
     * is more efficient.
     */
    virtual void invalidateMolecules() = 0;
    /**
     * Mark that the current molecule definitions from one particular force (and hence the atom order)
     * may be invalid.  This should be called whenever force field parameters change.  It will cause the
     * definitions and order to be revalidated.
     */
    virtual bool invalidateMolecules(ComputeForceInfo* force) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_COMPUTECONTEXT_H_*/
