#ifndef OPENMM_BONDEDUTILITIES_H_
#define OPENMM_BONDEDUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011-2019 Stanford University and the Authors.      *
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
#include <string>
#include <vector>

namespace OpenMM {

/**
 * This abstract class defines an interface for computing bonded interactions.  Call
 * getBondedUtilities() on a ComputeContext to get the BondedUtilities object for that
 * context.
 * 
 * This class provides a generic mechanism for evaluating bonded interactions.  You write only
 * the source code needed to compute one interaction, and this object takes care of creating
 * and executing a complete kernel that loops over bonds, evaluates each one, and accumulates
 * the resulting forces and energies.  This offers two advantages.  First, it simplifies the
 * task of writing a new Force.  Second, it allows multiple forces to be evaluated by a single
 * kernel, which reduces overhead and improves performance.
 * 
 * A "bonded interaction" means an interaction that affects a small, fixed set of particles.
 * The interaction energy may depend on the positions of only those particles, and the list of
 * particles forming a "bond" may not change with time.  Examples of bonded interactions
 * include HarmonicBondForce, HarmonicAngleForce, and PeriodicTorsionForce.
 * 
 * To create a bonded interaction, call addInteraction().  You pass to it a block of source
 * code for evaluating the interaction.  The inputs and outputs for that source code are as
 * follows:
 * 
 * <ol>
 * <li>The index of the bond being evaluated will have been stored in the unsigned int variable "index".</li>
 * <li>The indices of the atoms forming that bond will have been stored in the unsigned int variables "atom1",
 * "atom2", ....</li>
 * <li>The positions of those atoms will have been stored in the real4 variables "pos1", "pos2", ....</li>
 * <li>A real variable called "energy" will exist.  Your code should add the potential energy of the
 * bond to that variable.</li>
 * <li>Your code should define real3 variables called "force1", "force2", ... that contain the force to
 * apply to each atom.</li>
 * </ol>
 * 
 * As a simple example, the following source code would be used to implement a pairwise interaction of
 * the form E=r^2:
 * 
 * <tt><pre>
 * real4 delta = pos2-pos1;
 * energy += delta.x*delta.x + delta.y*delta.y + delta.z*delta.z;
 * real3 force1 = 2.0f*delta;
 * real3 force2 = -2.0f*delta;
 * </pre></tt>
 * 
 * Interactions will often depend on parameters or other data.  Call addArgument() to provide the data
 * to this class.  It will be passed to the interaction kernel as an argument, and you can refer to it
 * from your interaction code.
 */

class OPENMM_EXPORT_COMMON BondedUtilities {
public:
    virtual ~BondedUtilities() {
    }
    /**
     * Add a bonded interaction.
     *
     * @param atoms    this should have one entry for each bond, and that entry should contain the list
     *                 of atoms involved in the bond.  Every entry must have the same number of atoms.
     * @param source   the code to evaluate the interaction
     * @param group    the force group in which the interaction should be calculated
     */
    virtual void addInteraction(const std::vector<std::vector<int> >& atoms, const std::string& source, int group) = 0;
    /**
     * Add an argument that should be passed to the interaction kernel.
     * 
     * @param data    the array containing the data to pass
     * @param type    the data type contained in the memory (e.g. "float4")
     * @return the name that will be used for the argument.  Any code you pass to addInteraction() should
     * refer to it by this name.
     */
    virtual std::string addArgument(ArrayInterface& data, const std::string& type) = 0;
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
     * Add some code that should be included in the program, before the start of the kernel.
     * This can be used, for example, to define functions that will be called by the kernel.
     * 
     * @param source   the code to include
     */
    virtual void addPrefixCode(const std::string& source) = 0;
};

} // namespace OpenMM

#endif /*OPENMM_BONDEDUTILITIES_H_*/
