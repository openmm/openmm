#ifndef OPENMM_CUSTOMCENTROIDBONDFORCEIMPL_H_
#define OPENMM_CUSTOMCENTROIDBONDFORCEIMPL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2015 Stanford University and the Authors.      *
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
#include "openmm/CustomCentroidBondForce.h"
#include "openmm/Kernel.h"
#include "openmm/System.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/ParsedExpression.h"
#include <utility>
#include <map>
#include <set>
#include <string>

namespace OpenMM {

/**
 * This is the internal implementation of CustomCentroidBondForce.
 */

class OPENMM_EXPORT CustomCentroidBondForceImpl : public ForceImpl {
public:
    CustomCentroidBondForceImpl(const CustomCentroidBondForce& owner);
    ~CustomCentroidBondForceImpl();
    void initialize(ContextImpl& context);
    const CustomCentroidBondForce& getOwner() const {
        return owner;
    }
    void updateContextState(ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters();
    std::vector<std::string> getKernelNames();
    std::vector<std::pair<int, int> > getBondedParticles() const;
    void updateParametersInContext(ContextImpl& context);
    /**
     * This is a utility routine that parses the energy expression, identifies the angles and dihedrals
     * in it, and replaces them with variables.
     *
     * @param force     the CustomCentroidBondForce to process
     * @param functions definitions of custom function that may appear in the expression
     * @param distances on exit, this will contain an entry for each distance used in the expression.  The key is the name
     *                  of the corresponding variable, and the value is the list of particle indices.
     * @param angles    on exit, this will contain an entry for each angle used in the expression.  The key is the name
     *                  of the corresponding variable, and the value is the list of particle indices.
     * @param dihedrals on exit, this will contain an entry for each dihedral used in the expression.  The key is the name
     *                  of the corresponding variable, and the value is the list of particle indices.
     * @return a Parsed expression for the energy
     */
    static Lepton::ParsedExpression prepareExpression(const CustomCentroidBondForce& force, const std::map<std::string, Lepton::CustomFunction*>& functions, std::map<std::string, std::vector<int> >& distances,
            std::map<std::string, std::vector<int> >& angles, std::map<std::string, std::vector<int> >& dihedrals);
    /**
     * Compute the normalized weights to use for each particle in each group.
     *
     * @param force     the CustomCentroidBondForce to process
     * @param system    the System it is part of
     * @param weights   on exit, weights[i][j] contains the normalized weight for particle j in group i.
     */
    static void computeNormalizedWeights(const CustomCentroidBondForce& force, const System& system, std::vector<std::vector<double> >& weights);
private:
    class FunctionPlaceholder;
    static Lepton::ExpressionTreeNode replaceFunctions(const Lepton::ExpressionTreeNode& node, std::map<std::string, int> atoms,
            std::map<std::string, std::vector<int> >& distances, std::map<std::string, std::vector<int> >& angles,
            std::map<std::string, std::vector<int> >& dihedrals, std::set<std::string>& variables);
    void addBondsBetweenGroups(int group1, int group2, std::vector<std::pair<int, int> >& bonds) const;
    const CustomCentroidBondForce& owner;
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_CUSTOMCENTROIDBONDFORCEIMPL_H_*/
