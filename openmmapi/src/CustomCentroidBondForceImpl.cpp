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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCentroidBondForceImpl.h"
#include "openmm/kernels.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include <sstream>
#include <utility>

using namespace OpenMM;
using namespace std;
using Lepton::CustomFunction;
using Lepton::ExpressionTreeNode;
using Lepton::Operation;
using Lepton::ParsedExpression;

/**
 * This class serves as a placeholder for angles and dihedrals in expressions.
 */
class CustomCentroidBondForceImpl::FunctionPlaceholder : public CustomFunction {
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

CustomCentroidBondForceImpl::CustomCentroidBondForceImpl(const CustomCentroidBondForce& owner) : owner(owner) {
}

CustomCentroidBondForceImpl::~CustomCentroidBondForceImpl() {
}

void CustomCentroidBondForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomCentroidBondForceKernel::Name(), context);

    // Check for errors in the specification of parameters and exclusions.

    const System& system = context.getSystem();
    vector<int> particles;
    vector<double> weights;
    for (int i = 0; i < owner.getNumGroups(); i++) {
        owner.getGroupParameters(i, particles, weights);
        for (int particle : particles)
            if (particle < 0 || particle >= system.getNumParticles()) {
                stringstream msg;
                msg << "CustomCentroidBondForce: Illegal particle index for a group: ";
                msg << particle;
                throw OpenMMException(msg.str());
            }
        if (weights.size() != particles.size() && weights.size() > 0) {
            stringstream msg;
            msg << "CustomCentroidBondForce: Wrong number of weights for group ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    vector<int> groups;
    vector<double> parameters;
    int numBondParameters = owner.getNumPerBondParameters();
    for (int i = 0; i < owner.getNumBonds(); i++) {
        owner.getBondParameters(i, groups, parameters);
        for (int group : groups)
            if (group < 0 || group >= owner.getNumGroups()) {
                stringstream msg;
                msg << "CustomCentroidBondForce: Illegal group index for a bond: ";
                msg << group;
                throw OpenMMException(msg.str());
            }
        if (parameters.size() != numBondParameters) {
            stringstream msg;
            msg << "CustomCentroidBondForce: Wrong number of parameters for bond ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    kernel.getAs<CalcCustomCentroidBondForceKernel>().initialize(context.getSystem(), owner);
}

double CustomCentroidBondForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCustomCentroidBondForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomCentroidBondForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomCentroidBondForceKernel::Name());
    return names;
}

map<string, double> CustomCentroidBondForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

ParsedExpression CustomCentroidBondForceImpl::prepareExpression(const CustomCentroidBondForce& force, const map<string, CustomFunction*>& customFunctions, map<string, vector<int> >& distances,
        map<string, vector<int> >& angles, map<string, vector<int> >& dihedrals) {
    CustomCentroidBondForceImpl::FunctionPlaceholder custom(1);
    CustomCentroidBondForceImpl::FunctionPlaceholder distance(2);
    CustomCentroidBondForceImpl::FunctionPlaceholder angle(3);
    CustomCentroidBondForceImpl::FunctionPlaceholder dihedral(4);
    map<string, CustomFunction*> functions = customFunctions;
    functions["distance"] = &distance;
    functions["angle"] = &angle;
    functions["dihedral"] = &dihedral;
    ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions);
    map<string, int> groups;
    set<string> variables;
    for (int i = 0; i < force.getNumGroupsPerBond(); i++) {
        stringstream name, x, y, z;
        name << 'g' << (i+1);
        x << 'x' << (i+1);
        y << 'y' << (i+1);
        z << 'z' << (i+1);
        groups[name.str()] = i;
        variables.insert(x.str());
        variables.insert(y.str());
        variables.insert(z.str());
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        variables.insert(force.getGlobalParameterName(i));
    for (int i = 0; i < force.getNumPerBondParameters(); i++)
        variables.insert(force.getPerBondParameterName(i));
    return ParsedExpression(replaceFunctions(expression.getRootNode(), groups, distances, angles, dihedrals, variables)).optimize();
}

ExpressionTreeNode CustomCentroidBondForceImpl::replaceFunctions(const ExpressionTreeNode& node, map<string, int> groups,
        map<string, vector<int> >& distances, map<string, vector<int> >& angles, map<string, vector<int> >& dihedrals, set<string>& variables) {
    const Operation& op = node.getOperation();
    if (op.getId() == Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("CustomCentroidBondForce: Unknown variable '"+op.getName()+"'");
    if (op.getId() != Operation::CUSTOM || (op.getName() != "distance" && op.getName() != "angle" && op.getName() != "dihedral"))
    {
        // This is not an angle or dihedral, so process its children.

        vector<ExpressionTreeNode> children;
        for (auto& child : node.getChildren())
            children.push_back(replaceFunctions(child, groups, distances, angles, dihedrals, variables));
        return ExpressionTreeNode(op.clone(), children);
    }
    const Operation::Custom& custom = static_cast<const Operation::Custom&>(op);

    // Identify the groups this term is based on.

    int numArgs = custom.getNumArguments();
    vector<int> indices(numArgs);
    for (int i = 0; i < numArgs; i++) {
        map<string, int>::const_iterator iter = groups.find(node.getChildren()[i].getOperation().getName());
        if (iter == groups.end())
            throw OpenMMException("CustomCentroidBondForce: Unknown group '"+node.getChildren()[i].getOperation().getName()+"'");
        indices[i] = iter->second;
    }

    // Select a name for the variable and add it to the appropriate map.

    stringstream variable;
    if (numArgs == 2)
        variable << "distance";
    else if (numArgs == 3)
        variable << "angle";
    else
        variable << "dihedral";
    for (int i = 0; i < numArgs; i++)
        variable << indices[i];
    string name = variable.str();
    if (numArgs == 2)
        distances[name] = indices;
    else if (numArgs == 3)
        angles[name] = indices;
    else
        dihedrals[name] = indices;

    // Return a new node that represents it as a simple variable.

    return ExpressionTreeNode(new Operation::Variable(name));
}

vector<pair<int, int> > CustomCentroidBondForceImpl::getBondedParticles() const {
    vector<pair<int, int> > bonds;
    for (int i = 0; i < owner.getNumBonds(); i++) {
        vector<int> groups;
        vector<double> parameters;
        owner.getBondParameters(i, groups, parameters);
        for (int j = 1; j < groups.size(); j++)
            for (int k = 0; k < j; k++)
                addBondsBetweenGroups(j, k, bonds);
    }
    return bonds;
}

void CustomCentroidBondForceImpl::addBondsBetweenGroups(int group1, int group2, vector<pair<int, int> >& bonds) const {
    vector<int> atoms1;
    vector<int> atoms2;
    vector<double> weights;
    owner.getGroupParameters(group1, atoms1, weights);
    owner.getGroupParameters(group2, atoms2, weights);
    for (int i = 0; i < atoms1.size(); i++)
        for (int j = 0; j < atoms2.size(); j++)
            bonds.push_back(make_pair(atoms1[i], atoms2[j]));
}

void CustomCentroidBondForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCustomCentroidBondForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

void CustomCentroidBondForceImpl::computeNormalizedWeights(const CustomCentroidBondForce& force, const System& system, vector<vector<double> >& weights) {
    int numGroups = force.getNumGroups();
    weights.resize(numGroups);
    for (int i = 0; i < numGroups; i++) {
        vector<int> particles;
        vector<double> groupWeights;
        force.getGroupParameters(i, particles, groupWeights);
        int numParticles = particles.size();

        // If weights were not specified, use particle masses.

        if (groupWeights.size() == 0) {
            groupWeights.resize(numParticles);
            for (int j = 0; j < numParticles; j++)
                groupWeights[j] = system.getParticleMass(particles[j]);
        }

        // Normalize the weights.

        double total = 0;
        for (int j = 0; j < numParticles; j++)
            total += groupWeights[j];
        if (total == 0.0) {
            stringstream msg;
            msg << "CustomCentroidBondForce: Weights for group ";
            msg << i;
            msg << " add to 0";
            throw OpenMMException(msg.str());
        }
        weights[i].resize(numParticles);
        for (int j = 0; j < numParticles; j++)
            weights[i][j] = groupWeights[j]/total;
    }
}
