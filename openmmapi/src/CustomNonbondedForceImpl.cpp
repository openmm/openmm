/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include "openmm/internal/SplineFitter.h"
#include "openmm/kernels.h"
#include "ReferenceTabulatedFunction.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include <cmath>
#include <sstream>
#include <utility>
#include <algorithm>

using namespace OpenMM;
using namespace std;

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
    // Check that all interaction groups only specify particles that have been defined.
    for (int group = 0; group < owner.getNumInteractionGroups(); group++) {
        set<int> set1, set2;
        owner.getInteractionGroupParameters(group, set1, set2);
        for (set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
            if ((*it < 0) || (*it >= owner.getNumParticles())) {
                stringstream msg;
                msg << "CustomNonbondedForce: Interaction group " << group << " set1 contains a particle index (" << *it << ") "
                    << "not present in system (" << owner.getNumParticles() << " particles).";
                throw OpenMMException(msg.str());
            }
        for (set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
            if ((*it < 0) || (*it >= owner.getNumParticles())) {
                stringstream msg;
                msg << "CustomNonbondedForce: Interaction group " << group << " set2 contains a particle index (" << *it << ") "
                    << "not present in system (" << owner.getNumParticles() << " particles).";
                throw OpenMMException(msg.str());
            }
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
    context.systemChanged();
}

void CustomNonbondedForceImpl::calcLongRangeCorrection(const CustomNonbondedForce& force, const Context& context, double& coefficient, vector<double>& derivatives) {
    if (force.getNonbondedMethod() == CustomNonbondedForce::NoCutoff || force.getNonbondedMethod() == CustomNonbondedForce::CutoffNonPeriodic) {
        coefficient = 0.0;
        return;
    }
    
    // Identify all particle classes (defined by parameters), and record the class of each particle.
    
    int numParticles = force.getNumParticles();
    vector<vector<double> > classes;
    map<vector<double>, int> classIndex;
    vector<int> atomClass(numParticles);
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        map<vector<double>, int>::iterator entry = classIndex.find(parameters);
        if (entry == classIndex.end()) {
            classIndex[parameters] = classes.size();
            atomClass[i] = classes.size();
            classes.push_back(parameters);
        }
        else
            atomClass[i] = entry->second;
    }
    int numClasses = classes.size();
    
    // Count the total number of particle pairs for each pair of classes.
    
    map<pair<int, int>, long long int> interactionCount;
    if (force.getNumInteractionGroups() == 0) {
        // Count the particles of each class.
        
        vector<long long int> classCounts(numClasses, 0);
        for (int i = 0; i < numParticles; i++)
            classCounts[atomClass[i]]++;
        for (int i = 0; i < numClasses; i++) {
            interactionCount[make_pair(i, i)] = (classCounts[i]*(classCounts[i]+1))/2;
            for (int j = i+1; j < numClasses; j++)
                interactionCount[make_pair(i, j)] = classCounts[i]*classCounts[j];
        }
    }
    else {
        // Initialize the counts to 0.
        
        for (int i = 0; i < numClasses; i++) {
            for (int j = i; j < numClasses; j++)
                interactionCount[make_pair(i, j)] = 0;
        }
        
        // Loop over interaction groups and count the interactions in each one.
        
        for (int group = 0; group < force.getNumInteractionGroups(); group++) {
            set<int> set1, set2;
            force.getInteractionGroupParameters(group, set1, set2);
            for (set<int>::const_iterator a1 = set1.begin(); a1 != set1.end(); ++a1)
                for (set<int>::const_iterator a2 = set2.begin(); a2 != set2.end(); ++a2) {
                    if (*a1 >= *a2 && set1.find(*a2) != set1.end() && set2.find(*a1) != set2.end())
                        continue;
                    int class1 = atomClass[*a1];
                    int class2 = atomClass[*a2];
                    interactionCount[make_pair(min(class1, class2), max(class1, class2))]++;
                }
        }
    }
    
    // Compute the coefficient.
    
    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
    double nPart = (double) numParticles;
    double numInteractions = (nPart*(nPart+1))/2;
    Lepton::CompiledExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).createCompiledExpression();
    vector<string> paramNames;
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        stringstream name1, name2;
        name1 << force.getPerParticleParameterName(i) << 1;
        name2 << force.getPerParticleParameterName(i) << 2;
        paramNames.push_back(name1.str());
        paramNames.push_back(name2.str());
    }
    double sum = 0;
    for (int i = 0; i < numClasses; i++)
        for (int j = i; j < numClasses; j++)
            sum += interactionCount[make_pair(i, j)]*integrateInteraction(expression, classes[i], classes[j], force, context, paramNames);
    sum /= numInteractions;
    coefficient = 2*M_PI*nPart*nPart*sum;
    
    // Now do the same for parameter derivatives.
    
    int numDerivs = force.getNumEnergyParameterDerivatives();
    derivatives.resize(numDerivs);
    for (int k = 0; k < numDerivs; k++) {
        expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).differentiate(force.getEnergyParameterDerivativeName(k)).createCompiledExpression();
        sum = 0;
        for (int i = 0; i < numClasses; i++)
            for (int j = i; j < numClasses; j++)
                sum += interactionCount[make_pair(i, j)]*integrateInteraction(expression, classes[i], classes[j], force, context, paramNames);
        sum /= numInteractions;
        derivatives[k] = 2*M_PI*nPart*nPart*sum;
    }
}

double CustomNonbondedForceImpl::integrateInteraction(Lepton::CompiledExpression& expression, const vector<double>& params1, const vector<double>& params2,
        const CustomNonbondedForce& force, const Context& context, const vector<string>& paramNames) {
    const set<string>& variables = expression.getVariables();
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        if (variables.find(paramNames[2*i]) != variables.end())
            expression.getVariableReference(paramNames[2*i]) = params1[i];
        if (variables.find(paramNames[2*i+1]) != variables.end())
            expression.getVariableReference(paramNames[2*i+1]) = params2[i];
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        if (variables.find(name) != variables.end())
            expression.getVariableReference(name) = context.getParameter(name);
    }
    
    // To integrate from r_cutoff to infinity, make the change of variables x=r_cutoff/r and integrate from 0 to 1.
    // This introduces another r^2 into the integral, which along with the r^2 in the formula for the correction
    // means we multiply the function by r^4.  Use the midpoint method.

    double* rPointer;
    try {
        rPointer = &expression.getVariableReference("r");
    }
    catch (exception& ex) {
        throw OpenMMException("CustomNonbondedForce: Cannot use long range correction with a force that does not depend on r.");
    }
    double cutoff = force.getCutoffDistance();
    double sum = 0;
    int numPoints = 1;
    for (int iteration = 0; ; iteration++) {
        double oldSum = sum;
        double newSum = 0;
        for (int i = 0; i < numPoints; i++) {
            if (i%3 == 1)
                continue;
            double x = (i+0.5)/numPoints;
            double r = cutoff/x;
            *rPointer = r;
            double r2 = r*r;
            newSum += expression.evaluate()*r2*r2;
        }
        sum = newSum/numPoints + oldSum/3;
        if (iteration > 2 && (fabs((sum-oldSum)/sum) < 1e-5 || sum == 0))
            break;
        if (iteration == 8)
            throw OpenMMException("CustomNonbondedForce: Long range correction did not converge.  Does the energy go to 0 faster than 1/r^2?");
        numPoints *= 3;
    }
    
    // If a switching function is used, integrate over the switching interval.
    
    double sum2 = 0;
    if (force.getUseSwitchingFunction()) {
        double rswitch = force.getSwitchingDistance();
        sum2 = 0;
        numPoints = 1;
        for (int iteration = 0; ; iteration++) {
            double oldSum = sum2;
            double newSum = 0;
            for (int i = 0; i < numPoints; i++) {
                if (i%3 == 1)
                    continue;
                double x = (i+0.5)/numPoints;
                double r = rswitch+x*(cutoff-rswitch);
                double switchValue = x*x*x*(10+x*(-15+x*6));
                *rPointer = r;
                newSum += switchValue*expression.evaluate()*r*r;
            }
            sum2 = newSum/numPoints + oldSum/3;
            if (iteration > 2 && (fabs((sum2-oldSum)/sum2) < 1e-5 || sum2 == 0))
                break;
            if (iteration == 8)
                throw OpenMMException("CustomNonbondedForce: Long range correction did not converge.  Is the energy finite everywhere in the switching interval?");
            numPoints *= 3;
        }
        sum2 *= cutoff-rswitch;
    }
    return sum/cutoff+sum2;
}
