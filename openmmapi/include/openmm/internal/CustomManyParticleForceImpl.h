#ifndef OPENMM_CUSTOMMANYPARTICLEFORCEIMPL_H_
#define OPENMM_CUSTOMMANYPARTICLEFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2014 Stanford University and the Authors.      *
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

#include "ForceImpl.h"
#include "openmm/CustomManyParticleForce.h"
#include "openmm/Kernel.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/ParsedExpression.h"
#include <utility>
#include <map>
#include <set>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of CustomManyParticleForce.
 */

class OPENMM_EXPORT CustomManyParticleForceImpl : public ForceImpl {
public:
    CustomManyParticleForceImpl(const CustomManyParticleForce& owner);
    ~CustomManyParticleForceImpl();
    void initialize(ContextImpl& context);
    const CustomManyParticleForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    void updateParametersInContext(ContextImpl& context);
    /**
     * This is a utility routine that parses the energy expression, identifies the angles and dihedrals
     * in it, and replaces them with variables.
     *
     * @param force     the CustomManyParticleForce to process
     * @param functions definitions of custom function that may appear in the expression
     * @param distances on exit, this will contain an entry for each distance used in the expression.  The key is the name
     *                  of the corresponding variable, and the value is the list of particle indices.
     * @param angles    on exit, this will contain an entry for each angle used in the expression.  The key is the name
     *                  of the corresponding variable, and the value is the list of particle indices.
     * @param dihedrals on exit, this will contain an entry for each dihedral used in the expression.  The key is the name
     *                  of the corresponding variable, and the value is the list of particle indices.
     * @return a Parsed expression for the energy
     */
    static Lepton::ParsedExpression prepareExpression(const CustomManyParticleForce& force, const std::map<std::string, Lepton::CustomFunction*>& functions, std::map<std::string, std::vector<int> >& distances,
            std::map<std::string, std::vector<int> >& angles, std::map<std::string, std::vector<int> >& dihedrals);
    /**
     * Analyze the type filters for a force and build a set of arrays that can be used for reordering the
     * particles in an interaction.
     * 
     * @param force          the CustomManyParticleForce to process
     * @param numTypes       on exit, the number of unique particle types
     * @param particleTypes  on exit, this contains a type code for each particle.  These codes are <i>not</i> necessarily the
     *                       same as the types assigned by the force.  They are guaranteed to be successive integers starting from 0,
     *                       whereas the force may have used arbitrary integers.
     * @param orderIndex     on exit, this contains a lookup table for selecting the particle order for an interaction.
     *                       orderIndex[t1+numTypes*t2+numTypes*numTypes*t3+...] is the index of the order to use, where t1, t2, etc. are the type codes
     *                       of the particles involved in the interaction.  If this equals -1, the interaction should be omitted.
     * @param particleOrder  on exit, particleOrder[i][j] tells which particle to use as the j'th particle, where i is the value found in orderIndex.
     */
    static void buildFilterArrays(const CustomManyParticleForce& force, int& numTypes, std::vector<int>& particleTypes, std::vector<int>& orderIndex, std::vector<std::vector<int> >& particleOrder);
private:
    class FunctionPlaceholder;
    static Lepton::ExpressionTreeNode replaceFunctions(const Lepton::ExpressionTreeNode& node, std::map<std::string, int> atoms,
            std::map<std::string, std::vector<int> >& distances, std::map<std::string, std::vector<int> >& angles,
            std::map<std::string, std::vector<int> >& dihedrals, std::set<std::string>& variables);
    static void generatePermutations(std::vector<int>& values, int numFixed, std::vector<std::vector<int> >& result);
    const CustomManyParticleForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMMANYPARTICLEFORCEIMPL_H_*/
