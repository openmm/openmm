/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include "openmm/internal/SplineFitter.h"
#include "openmm/kernels.h"
#include "lepton/CustomFunction.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include <cmath>
#include <sstream>

using namespace OpenMM;
using std::map;
using std::pair;
using std::vector;
using std::set;
using std::string;
using std::stringstream;

CustomNonbondedForceImpl::CustomNonbondedForceImpl(const CustomNonbondedForce& owner) : owner(owner) {
}

CustomNonbondedForceImpl::~CustomNonbondedForceImpl() {
}

void CustomNonbondedForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomNonbondedForceKernel::Name(), context);

    // Check for errors in the specification of parameters and exclusions.

    const System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("CustomNonbondedForce must have exactly as many particles as the System it belongs to.");
    if (owner.getUseSwitchingFunction()) {
        if (owner.getSwitchingDistance() < 0 || owner.getSwitchingDistance() >= owner.getCutoffDistance())
            throw OpenMMException("CustomNonbondedForce: Switching distance must satisfy 0 <= r_switch < r_cutoff");
    }
    vector<set<int> > exclusions(owner.getNumParticles());
    vector<double> parameters;
    int numParameters = owner.getNumPerParticleParameters();
    for (int i = 0; i < owner.getNumParticles(); i++) {
        owner.getParticleParameters(i, parameters);
        if (parameters.size() != numParameters) {
            stringstream msg;
            msg << "CustomNonbondedForce: Wrong number of parameters for particle ";
            msg << i;
            throw OpenMMException(msg.str());
        }
    }
    for (int i = 0; i < owner.getNumExclusions(); i++) {
        int particle1, particle2;
        owner.getExclusionParticles(i, particle1, particle2);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "CustomNonbondedForce: Illegal particle index for an exclusion: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "CustomNonbondedForce: Illegal particle index for an exclusion: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (exclusions[particle1].count(particle2) > 0 || exclusions[particle2].count(particle1) > 0) {
            stringstream msg;
            msg << "CustomNonbondedForce: Multiple exclusions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }
    if (owner.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("CustomNonbondedForce: The cutoff distance cannot be greater than half the periodic box size.");
    }
    kernel.getAs<CalcCustomNonbondedForceKernel>().initialize(context.getSystem(), owner);
}

double CustomNonbondedForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCustomNonbondedForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomNonbondedForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomNonbondedForceKernel::Name());
    return names;
}

map<string, double> CustomNonbondedForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

void CustomNonbondedForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCustomNonbondedForceKernel>().copyParametersToContext(context, owner);
}

class CustomNonbondedForceImpl::TabulatedFunction : public Lepton::CustomFunction {
public:
    TabulatedFunction(double min, double max, const vector<double>& values) :
            min(min), max(max), values(values) {
        int numValues = values.size();
        x.resize(numValues);
        for (int i = 0; i < numValues; i++)
            x[i] = min+i*(max-min)/(numValues-1);
        SplineFitter::createNaturalSpline(x, values, derivs);
    }
    int getNumArguments() const {
        return 1;
    }
    double evaluate(const double* arguments) const {
        double t = arguments[0];
        if (t < min || t > max)
            return 0.0;
        return SplineFitter::evaluateSpline(x, values, derivs, t);
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        double t = arguments[0];
        if (t < min || t > max)
            return 0.0;
        return SplineFitter::evaluateSplineDerivative(x, values, derivs, t);
    }
    CustomFunction* clone() const {
        return new TabulatedFunction(min, max, values);
    }
    double min, max;
    vector<double> x, values, derivs;
};

double CustomNonbondedForceImpl::calcLongRangeCorrection(const CustomNonbondedForce& force, const Context& context) {
    if (force.getNonbondedMethod() == CustomNonbondedForce::NoCutoff || force.getNonbondedMethod() == CustomNonbondedForce::CutoffNonPeriodic)
        return 0.0;
    
    // Identify all particle classes (defined by parameters), and count the number of
    // particles in each class.

    map<vector<double>, int> classCounts;
    for (int i = 0; i < force.getNumParticles(); i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        map<vector<double>, int>::iterator entry = classCounts.find(parameters);
        if (entry == classCounts.end())
            classCounts[parameters] = 1;
        else
            entry->second++;
    }
    
    // Parse the energy expression.
    
    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        force.getFunctionParameters(i, name, values, min, max);
        functions[name] = new TabulatedFunction(min, max, values);
    }
    Lepton::ExpressionProgram expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize().createProgram();

    // Loop over all pairs of classes to compute the coefficient.

    double sum = 0;
    for (map<vector<double>, int>::const_iterator entry = classCounts.begin(); entry != classCounts.end(); ++entry) {
        int count = (entry->second*(entry->second+1))/2;
        sum += count*integrateInteraction(expression, entry->first, entry->first, force, context);
    }
    for (map<vector<double>, int>::const_iterator class1 = classCounts.begin(); class1 != classCounts.end(); ++class1)
        for (map<vector<double>, int>::const_iterator class2 = classCounts.begin(); class2 != class1; ++class2) {
            int count = class1->second*class2->second;
            sum += count*integrateInteraction(expression, class1->first, class2->first, force, context);
        }
    int numParticles = force.getNumParticles();
    int numInteractions = (numParticles*(numParticles+1))/2;
    sum /= numInteractions;
    return 2*M_PI*numParticles*numParticles*sum;
}

double CustomNonbondedForceImpl::integrateInteraction(const Lepton::ExpressionProgram& expression, const vector<double>& params1, const vector<double>& params2,
        const CustomNonbondedForce& force, const Context& context) {
    map<std::string, double> variables;
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        stringstream name1, name2;
        name1 << force.getPerParticleParameterName(i) << 1;
        name2 << force.getPerParticleParameterName(i) << 2;
        variables[name1.str()] = params1[i];
        variables[name2.str()] = params2[i];
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        variables[name] = context.getParameter(name);
    }
    
    // To integrate from r_cutoff to infinity, make the change of variables x=r_cutoff/r and integrate from 0 to 1.
    // This introduces another r^2 into the integral, which along with the r^2 in the formula for the correction
    // means we multiply the function by r^4.  Use the midpoint method.

    double cutoff = force.getCutoffDistance();
    variables["r"] = 2*cutoff;
    double sum = expression.evaluate(variables);
    int numPoints = 1;
    for (int iteration = 0; iteration < 10; iteration++) {
        double oldSum = sum;
        double newSum = 0;
        numPoints *= 3;
        for (int i = 0; i < numPoints; i++) {
            if (i%3 == 1)
                continue;
            double x = (i+0.5)/numPoints;
            double r = cutoff/x;
            variables["r"] = r;
            double r2 = r*r;
            newSum += expression.evaluate(variables)*r2*r2;
        }
        sum = newSum/numPoints + oldSum/3;
        if (iteration > 2 && (fabs((sum-oldSum)/sum) < 1e-5 || sum == 0))
            break;
    }
    return sum/cutoff;
}
