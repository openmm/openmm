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

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
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
class CustomManyParticleForceImpl::FunctionPlaceholder : public CustomFunction {
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

CustomManyParticleForceImpl::CustomManyParticleForceImpl(const CustomManyParticleForce& owner) : owner(owner) {
}

CustomManyParticleForceImpl::~CustomManyParticleForceImpl() {
}

void CustomManyParticleForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomManyParticleForceKernel::Name(), context);

    // Check for errors in the specification of parameters and exclusions.

    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("CustomManyParticleForce must have exactly as many particles as the System it belongs to.");
    vector<set<int> > exclusions(owner.getNumParticles());
    vector<double> parameters;
    int type;
    int numParameters = owner.getNumPerParticleParameters();
    for (int i = 0; i < owner.getNumParticles(); i++) {
        owner.getParticleParameters(i, parameters, type);
        if (parameters.size() != numParameters) {
            stringstream msg;
            msg << "CustomManyParticleForce: Wrong number of parameters for particle ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    for (int i = 0; i < owner.getNumExclusions(); i++) {
        int particle1, particle2;
        owner.getExclusionParticles(i, particle1, particle2);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "CustomManyParticleForce: Illegal particle index for an exclusion: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "CustomManyParticleForce: Illegal particle index for an exclusion: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (exclusions[particle1].count(particle2) > 0 || exclusions[particle2].count(particle1) > 0) {
            stringstream msg;
            msg << "CustomManyParticleForce: Multiple exclusions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }
    if (owner.getNonbondedMethod() == CustomManyParticleForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("CustomManyParticleForce: The cutoff distance cannot be greater than half the periodic box size.");
    }
    kernel.getAs<CalcCustomManyParticleForceKernel>().initialize(context.getSystem(), owner);
}

double CustomManyParticleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCustomManyParticleForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomManyParticleForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomManyParticleForceKernel::Name());
    return names;
}

map<string, double> CustomManyParticleForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

ParsedExpression CustomManyParticleForceImpl::prepareExpression(const CustomManyParticleForce& force, const map<string, CustomFunction*>& customFunctions, map<string, vector<int> >& distances,
        map<string, vector<int> >& angles, map<string, vector<int> >& dihedrals) {
    CustomManyParticleForceImpl::FunctionPlaceholder custom(1);
    CustomManyParticleForceImpl::FunctionPlaceholder distance(2);
    CustomManyParticleForceImpl::FunctionPlaceholder angle(3);
    CustomManyParticleForceImpl::FunctionPlaceholder dihedral(4);
    map<string, CustomFunction*> functions = customFunctions;
    functions["distance"] = &distance;
    functions["angle"] = &angle;
    functions["dihedral"] = &dihedral;
    ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions);
    map<string, int> atoms;
    set<string> variables;
    for (int i = 0; i < force.getNumParticlesPerSet(); i++) {
        stringstream name, x, y, z;
        name << 'p' << (i+1);
        x << 'x' << (i+1);
        y << 'y' << (i+1);
        z << 'z' << (i+1);
        atoms[name.str()] = i;
        variables.insert(x.str());
        variables.insert(y.str());
        variables.insert(z.str());
        for (int j = 0; j < force.getNumPerParticleParameters(); j++) {
            stringstream param;
            param << force.getPerParticleParameterName(j) << (i+1);
            variables.insert(param.str());
        }
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        variables.insert(force.getGlobalParameterName(i));
    return ParsedExpression(replaceFunctions(expression.getRootNode(), atoms, distances, angles, dihedrals, variables)).optimize();
}

ExpressionTreeNode CustomManyParticleForceImpl::replaceFunctions(const ExpressionTreeNode& node, map<string, int> atoms,
        map<string, vector<int> >& distances, map<string, vector<int> >& angles, map<string, vector<int> >& dihedrals, set<string>& variables) {
    const Operation& op = node.getOperation();
    if (op.getId() == Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("CustomManyParticleForce: Unknown variable '"+op.getName()+"'");
    if (op.getId() != Operation::CUSTOM || (op.getName() != "distance" && op.getName() != "angle" && op.getName() != "dihedral"))
    {
        // This is not an angle or dihedral, so process its children.

        vector<ExpressionTreeNode> children;
        for (auto& child : node.getChildren())
            children.push_back(replaceFunctions(child, atoms, distances, angles, dihedrals, variables));
        return ExpressionTreeNode(op.clone(), children);
    }
    const Operation::Custom& custom = static_cast<const Operation::Custom&>(op);

    // Identify the atoms this term is based on.

    int numArgs = custom.getNumArguments();
    vector<int> indices(numArgs);
    for (int i = 0; i < numArgs; i++) {
        map<string, int>::const_iterator iter = atoms.find(node.getChildren()[i].getOperation().getName());
        if (iter == atoms.end())
            throw OpenMMException("CustomManyParticleForce: Unknown particle '"+node.getChildren()[i].getOperation().getName()+"'");
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

void CustomManyParticleForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCustomManyParticleForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}

void CustomManyParticleForceImpl::buildFilterArrays(const CustomManyParticleForce& force, int& numTypes, vector<int>& particleTypes, vector<int>& orderIndex, vector<vector<int> >& particleOrder) {
    // Build a canonical list of type codes.
    
    int numParticles = force.getNumParticles();
    int numParticlesPerSet = force.getNumParticlesPerSet();
    particleTypes.resize(numParticles);
    map<int, int> typeMap;
    for (int i = 0; i < numParticles; i++) {
        vector<double> params;
        int type;
        force.getParticleParameters(i, params, type);
        map<int, int>::const_iterator element = typeMap.find(type);
        if (element == typeMap.end()) {
            int newType = typeMap.size();
            typeMap[type] = newType;
            particleTypes[i] = newType;
        }
        else
            particleTypes[i] = element->second;
    }
    numTypes = typeMap.size();
    int numIndices = 1;
    for (int i = 0; i < numParticlesPerSet; i++)
        numIndices *= numTypes;
    orderIndex.resize(numIndices, 0);
    
    // Find the allowed type codes for each particle in an interaction.
    
    vector<set<int> > allowedTypes(numParticlesPerSet);
    bool anyFilters = false;
    for (int i = 0; i < numParticlesPerSet; i++) {
        set<int> types;
        force.getTypeFilter(i, types);
        if (types.size() == 0)
            for (int j = 0; j < numTypes; j++)
                allowedTypes[i].insert(j);
        else {
            for (int type : types)
                if (typeMap.find(type) != typeMap.end())
                    allowedTypes[i].insert(typeMap[type]);
            if (allowedTypes[i].size() < numTypes)
                anyFilters = true;
        }
    }
    
    // If there are no filters, reordering is unnecessary.
    
    if (!anyFilters) {
        particleOrder.resize(1);
        particleOrder[0].resize(numParticlesPerSet);
        for (int i = 0; i < numParticlesPerSet; i++)
            particleOrder[0][i] = i;
        return;
    }
    
    // Build a list of every possible permutation of the particles.
    
    particleOrder.clear();
    vector<int> values;
    for (int i = 0; i < numParticlesPerSet; i++)
        values.push_back(i);
    generatePermutations(values, force.getPermutationMode() == CustomManyParticleForce::SinglePermutation ? 0 : 1, particleOrder);
    int numOrders = particleOrder.size();
    
    // Now we need to loop over every possible sequence of type codes, and for each one figure out which order to use.
    
    for (int i = 0; i < numIndices; i++) {
        vector<int> types(numParticlesPerSet);
        int temp = i;
        for (int j = 0; j < numParticlesPerSet; j++) {
            types[j] = temp%numTypes;
            temp /= numTypes;
        }
        
        // Loop over possible orders until we find one that matches the filters.
        
        int order = -1;
        for (int j = 0; j < numOrders && order == -1; j++) {
            bool matches = true;
            for (int k = 0; k < numParticlesPerSet && matches; k++)
                if (allowedTypes[k].find(types[particleOrder[j][k]]) == allowedTypes[k].end())
                    matches = false;
            if (matches)
                order = j;
        }
        orderIndex[i] = order;
    }
}

void CustomManyParticleForceImpl::generatePermutations(vector<int>& values, int numFixed, vector<vector<int> >& result) {
    int numValues = values.size();
    if (numFixed == numValues) {
        result.push_back(values);
        return;
    }
    for (int i = numFixed; i < numValues; i++) {
        int v1 = values[numFixed];
        int v2 = values[i];
        values[numFixed] = v2;
        values[i] = v1;
        generatePermutations(values, numFixed+1, result);
        values[numFixed] = v1;
        values[i] = v2;
    }
}