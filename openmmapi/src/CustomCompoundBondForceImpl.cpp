/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2021 Stanford University and the Authors.      *
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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/kernels.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include <sstream>

using namespace OpenMM;
using Lepton::CustomFunction;
using Lepton::ExpressionTreeNode;
using Lepton::Operation;
using Lepton::ParsedExpression;
using std::map;
using std::pair;
using std::vector;
using std::set;
using std::string;
using std::stringstream;

/**
 * This class serves as a placeholder for angles and dihedrals in expressions.
 */
class CustomCompoundBondForceImpl::FunctionPlaceholder : public CustomFunction {
public:
    int numArguments;
    FunctionPlaceholder(int numArguments) : numArguments(numArguments) {
    }
    int getNumArguments() const {
        return numArguments;
    }
    double evaluate(const double* arguments) const {
        return 0.0;
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        return 0.0;
    }
    CustomFunction* clone() const {
        return new FunctionPlaceholder(numArguments);
    }
};

CustomCompoundBondForceImpl::CustomCompoundBondForceImpl(const CustomCompoundBondForce& owner) : owner(owner) {
}

CustomCompoundBondForceImpl::~CustomCompoundBondForceImpl() {
}

void CustomCompoundBondForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomCompoundBondForceKernel::Name(), context);

    // Check for errors in the specification of parameters and exclusions.

    const System& system = context.getSystem();
    vector<int> particles;
    vector<double> parameters;
    int numBondParameters = owner.getNumPerBondParameters();
    for (int i = 0; i < owner.getNumBonds(); i++) {
        owner.getBondParameters(i, particles, parameters);
        for (int particle : particles)
            if (particle < 0 || particle >= system.getNumParticles()) {
                stringstream msg;
                msg << "CustomCompoundBondForce: Illegal particle index for a bond: ";
                msg << particle;
                throw OpenMMException(msg.str());
            }
        if (parameters.size() != numBondParameters) {
            stringstream msg;
            msg << "CustomCompoundBondForce: Wrong number of parameters for bond ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    kernel.getAs<CalcCustomCompoundBondForceKernel>().initialize(context.getSystem(), owner);
}

double CustomCompoundBondForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCustomCompoundBondForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomCompoundBondForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomCompoundBondForceKernel::Name());
    return names;
}

map<string, double> CustomCompoundBondForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

ParsedExpression CustomCompoundBondForceImpl::prepareExpression(const CustomCompoundBondForce& force, const map<string, CustomFunction*>& customFunctions) {
    CustomCompoundBondForceImpl::FunctionPlaceholder distance(2);
    CustomCompoundBondForceImpl::FunctionPlaceholder angle(3);
    CustomCompoundBondForceImpl::FunctionPlaceholder dihedral(4);
    CustomCompoundBondForceImpl::FunctionPlaceholder pointdistance(6);
    CustomCompoundBondForceImpl::FunctionPlaceholder pointangle(9);
    CustomCompoundBondForceImpl::FunctionPlaceholder pointdihedral(12);
    map<string, CustomFunction*> functions = customFunctions;
    functions["distance"] = &distance;
    functions["angle"] = &angle;
    functions["dihedral"] = &dihedral;
    if (functions.find("pointdistance") == functions.end())
        functions["pointdistance"] = &pointdistance;
    if (functions.find("pointangle") == functions.end())
        functions["pointangle"] = &pointangle;
    if (functions.find("pointdihedral") == functions.end())
        functions["pointdihedral"] = &pointdihedral;
    ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions);
    map<string, int> atoms;
    set<string> variables;
    for (int i = 0; i < force.getNumParticlesPerBond(); i++) {
        stringstream name, x, y, z;
        name << 'p' << (i+1);
        x << 'x' << (i+1);
        y << 'y' << (i+1);
        z << 'z' << (i+1);
        atoms[name.str()] = i;
        variables.insert(x.str());
        variables.insert(y.str());
        variables.insert(z.str());
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        variables.insert(force.getGlobalParameterName(i));
    for (int i = 0; i < force.getNumPerBondParameters(); i++)
        variables.insert(force.getPerBondParameterName(i));
    return ParsedExpression(replaceFunctions(expression.getRootNode(), atoms, functions, variables)).optimize();
}

ExpressionTreeNode CustomCompoundBondForceImpl::replaceFunctions(const ExpressionTreeNode& node, map<string, int> atoms,
        const map<string, CustomFunction*>& functions, set<string>& variables) {
    const Operation& op = node.getOperation();
    if (op.getId() == Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("CustomCompoundBondForce: Unknown variable '"+op.getName()+"'");
    vector<ExpressionTreeNode> children;
    if (op.getId() != Operation::CUSTOM || (op.getName() != "distance" && op.getName() != "angle" && op.getName() != "dihedral")) {
        // The arguments are not particle identifiers, so process its children.

        for (auto& child : node.getChildren())
            children.push_back(replaceFunctions(child, atoms, functions, variables));
        return ExpressionTreeNode(op.clone(), children);
    }
    const Operation::Custom& custom = static_cast<const Operation::Custom&>(op);

    // Identify the atoms this term is based on.

    int numArgs = custom.getNumArguments();
    vector<int> indices(numArgs);
    for (int i = 0; i < numArgs; i++) {
        map<string, int>::const_iterator iter = atoms.find(node.getChildren()[i].getOperation().getName());
        if (iter == atoms.end())
            throw OpenMMException("CustomCompoundBondForce: Unknown particle '"+node.getChildren()[i].getOperation().getName()+"'");
        indices[i] = iter->second;
    }

    // Replace it by the corresponding point based function.

    for (int i = 0; i < numArgs; i++) {
        stringstream x, y, z;
        x << 'x' << (indices[i]+1);
        y << 'y' << (indices[i]+1);
        z << 'z' << (indices[i]+1);
        children.push_back(ExpressionTreeNode(new Operation::Variable(x.str())));
        children.push_back(ExpressionTreeNode(new Operation::Variable(y.str())));
        children.push_back(ExpressionTreeNode(new Operation::Variable(z.str())));
    }
    if (op.getName() == "distance")
        return ExpressionTreeNode(new Operation::Custom("pointdistance", functions.at("pointdistance")->clone()), children);
    if (op.getName() == "angle")
        return ExpressionTreeNode(new Operation::Custom("pointangle", functions.at("pointangle")->clone()), children);
    if (op.getName() == "dihedral")
        return ExpressionTreeNode(new Operation::Custom("pointdihedral", functions.at("pointdihedral")->clone()), children);
    throw OpenMMException("Internal error.  Unexpected function '"+op.getName()+"'");
}

void CustomCompoundBondForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCustomCompoundBondForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}
