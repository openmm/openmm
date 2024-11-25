/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "openmm/internal/Messages.h"
#include "openmm/internal/SplineFitter.h"
#include "openmm/kernels.h"
#include "ReferenceTabulatedFunction.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include <atomic>
#include <cmath>
#include <sstream>
#include <utility>
#include <algorithm>

using namespace OpenMM;
using namespace std;

CustomNonbondedForceImpl::CustomNonbondedForceImpl(const CustomNonbondedForce& owner) : owner(owner) {
    forceGroup = owner.getForceGroup();
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
        int minp = min(particle1, particle2);
        int maxp = max(particle1, particle2);
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
        if (exclusions[minp].count(maxp) > 0) {
            stringstream msg;
            msg << "CustomNonbondedForce: Multiple exclusions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exclusions[minp].insert(maxp);
    }
    if (owner.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("CustomNonbondedForce: "+Messages::cutoffTooLarge);
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
    if (owner.getNumEnergyParameterDerivatives() > 0 && owner.getNumComputedValues() > 0)
        throw OpenMMException("CustomNonbondedForce: Cannot compute parameter derivatives for a force that uses computed values.");

    kernel.getAs<CalcCustomNonbondedForceKernel>().initialize(context.getSystem(), owner);
}

double CustomNonbondedForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<forceGroup)) != 0)
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

void CustomNonbondedForceImpl::updateParametersInContext(ContextImpl& context, int firstParticle, int lastParticle) {
    kernel.getAs<CalcCustomNonbondedForceKernel>().copyParametersToContext(context, owner, firstParticle, lastParticle);
    context.systemChanged();
}

CustomNonbondedForceImpl::LongRangeCorrectionData CustomNonbondedForceImpl::prepareLongRangeCorrection(const CustomNonbondedForce& force, int numThreads) {
    LongRangeCorrectionData data;
    data.method = force.getNonbondedMethod();
    if (data.method == CustomNonbondedForce::NoCutoff || data.method == CustomNonbondedForce::CutoffNonPeriodic)
        return data;
    
    // Identify all particle classes (defined by parameters), and record the class of each particle.
    
    int numParticles = force.getNumParticles();
    map<vector<double>, int> classIndex;
    vector<int> atomClass(numParticles);
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        map<vector<double>, int>::iterator entry = classIndex.find(parameters);
        if (entry == classIndex.end()) {
            classIndex[parameters] = data.classes.size();
            atomClass[i] = data.classes.size();
            data.classes.push_back(parameters);
        }
        else
            atomClass[i] = entry->second;
    }
    int numClasses = data.classes.size();
    
    // Count the total number of particle pairs for each pair of classes.
    
    if (force.getNumInteractionGroups() == 0) {
        // Count the particles of each class.
        
        vector<long long int> classCounts(numClasses, 0);
        for (int i = 0; i < numParticles; i++)
            classCounts[atomClass[i]]++;
        for (int i = 0; i < numClasses; i++) {
            data.interactionCount[make_pair(i, i)] = (classCounts[i]*(classCounts[i]+1))/2;
            for (int j = i+1; j < numClasses; j++)
                data.interactionCount[make_pair(i, j)] = classCounts[i]*classCounts[j];
        }
    }
    else {
        // Initialize the counts to 0.
        
        for (int i = 0; i < numClasses; i++) {
            for (int j = i; j < numClasses; j++)
                data.interactionCount[make_pair(i, j)] = 0;
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
                    data.interactionCount[make_pair(min(class1, class2), max(class1, class2))]++;
                }
        }
    }
    
    // Prepare for evaluating the expressions.
    
    int width = Lepton::CompiledVectorExpression::getAllowedWidths().back();
    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
    Lepton::CompiledVectorExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).createCompiledVectorExpression(width);
    for (int i = 0; i < numThreads; i++)
        data.energyExpression.push_back(energyExpression);
    data.derivExpressions.resize(numThreads);
    for (int k = 0; k < force.getNumEnergyParameterDerivatives(); k++) {
        Lepton::CompiledVectorExpression derivExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).differentiate(force.getEnergyParameterDerivativeName(k)).createCompiledVectorExpression(width);
        for (int i = 0; i < numThreads; i++)
            data.derivExpressions[i].push_back(derivExpression);
    }
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        string name = force.getPerParticleParameterName(i);
        data.paramNames.push_back(name+"1");
        data.paramNames.push_back(name+"2");
    }
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, exp;
        force.getComputedValueParameters(i, name, exp);
        data.computedValueNames.push_back(name+"1");
        data.computedValueNames.push_back(name+"2");
        data.computedValueExpressions.push_back(Lepton::Parser::parse(exp, functions).createCompiledExpression());
    }
    return data;
}

void CustomNonbondedForceImpl::calcLongRangeCorrection(const CustomNonbondedForce& force, LongRangeCorrectionData& data, const Context& context, double& coefficient, vector<double>& derivatives, ThreadPool& threads) {
    if (data.method == CustomNonbondedForce::NoCutoff || data.method == CustomNonbondedForce::CutoffNonPeriodic) {
        coefficient = 0.0;
        return;
    }
    
    // Calculate the computed values for all atom classes.
    
    int numClasses = data.classes.size();
    vector<vector<double> > computedValues(numClasses, vector<double>(force.getNumComputedValues()));
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        Lepton::CompiledExpression& expression = data.computedValueExpressions[i];
        const set<string>& variables = expression.getVariables();
        for (int j = 0; j < force.getNumGlobalParameters(); j++) {
            const string& name = force.getGlobalParameterName(j);
            if (variables.find(name) != variables.end())
                expression.getVariableReference(name) = context.getParameter(name);
        }
        for (int j = 0; j < numClasses; j++) {
            for (int k = 0; k < force.getNumPerParticleParameters(); k++) {
                const string& name = force.getPerParticleParameterName(k);
                if (variables.find(name) != variables.end())
                    expression.getVariableReference(name) = data.classes[j][k];
            }
            computedValues[j][i] = expression.evaluate();
        }
    }

    // Compute the coefficient.  Use multiple threads to compute the integrals in parallel.

    double nPart = (double) context.getSystem().getNumParticles();
    double numInteractions = (nPart*(nPart+1))/2;
    vector<double> threadSum(threads.getNumThreads(), 0.0);
    atomic<int> atomicCounter(0);
    threads.execute([&] (ThreadPool& threads, int threadIndex) {
        Lepton::CompiledVectorExpression& expression = data.energyExpression[threadIndex];
        while (true) {
            int i = atomicCounter++;
            if (i >= numClasses)
                break;
            for (int j = i; j < numClasses; j++)
                threadSum[threadIndex] += data.interactionCount.at(make_pair(i, j))*integrateInteraction(expression, data.classes[i], data.classes[j],
                        computedValues[i], computedValues[j], force, context, data.paramNames, data.computedValueNames);
        }
    });
    threads.waitForThreads();
    double sum = 0;
    for (int i = 0; i < threadSum.size(); i++)
        sum += threadSum[i];
    sum /= numInteractions;
    coefficient = 2*M_PI*nPart*nPart*sum;
    
    // Now do the same for parameter derivatives.
    
    int numDerivs = data.derivExpressions[0].size();
    derivatives.resize(numDerivs);
    for (int k = 0; k < numDerivs; k++) {
        atomicCounter = 0;
        threads.execute([&] (ThreadPool& threads, int threadIndex) {
            threadSum[threadIndex] = 0;
            Lepton::CompiledVectorExpression& expression = data.derivExpressions[threadIndex][k];
            while (true) {
                int i = atomicCounter++;
                if (i >= numClasses)
                    break;
                for (int j = i; j < numClasses; j++)
                    threadSum[threadIndex] += data.interactionCount.at(make_pair(i, j))*integrateInteraction(expression, data.classes[i], data.classes[j],
                            computedValues[i], computedValues[j], force, context, data.paramNames, data.computedValueNames);
            }
        });
        threads.waitForThreads();
        sum = 0;
        for (int i = 0; i < threadSum.size(); i++)
            sum += threadSum[i];
        sum /= numInteractions;
        derivatives[k] = 2*M_PI*nPart*nPart*sum;
    }
}

double CustomNonbondedForceImpl::integrateInteraction(Lepton::CompiledVectorExpression& expression, const vector<double>& params1, const vector<double>& params2,
        const vector<double>& computedValues1, const vector<double>& computedValues2, const CustomNonbondedForce& force, const Context& context,
        const vector<string>& paramNames, const vector<string>& computedValueNames) {
    int width = expression.getWidth();
    const set<string>& variables = expression.getVariables();
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        if (variables.find(paramNames[2*i]) != variables.end()) {
            float* pointer = expression.getVariablePointer(paramNames[2*i]);
            for (int j = 0; j < width; j++)
                pointer[j] = params1[i];
        }
        if (variables.find(paramNames[2*i+1]) != variables.end()) {
            float* pointer = expression.getVariablePointer(paramNames[2*i+1]);
            for (int j = 0; j < width; j++)
                pointer[j] = params2[i];
        }
    }
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        if (variables.find(computedValueNames[2*i]) != variables.end()) {
            float* pointer = expression.getVariablePointer(computedValueNames[2*i]);
            for (int j = 0; j < width; j++)
                pointer[j] = computedValues1[i];
        }
        if (variables.find(computedValueNames[2*i+1]) != variables.end()) {
            float* pointer = expression.getVariablePointer(computedValueNames[2*i+1]);
            for (int j = 0; j < width; j++)
                pointer[j] = computedValues2[i];
        }
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        if (variables.find(name) != variables.end()) {
            float* pointer = expression.getVariablePointer(name);
            for (int j = 0; j < width; j++)
                pointer[j] = context.getParameter(name);
        }
    }

    // To integrate from r_cutoff to infinity, make the change of variables x=r_cutoff/r and integrate from 0 to 1.
    // This introduces another r^2 into the integral, which along with the r^2 in the formula for the correction
    // means we multiply the function by r^4.  Use the midpoint method.

    float* r;
    try {
        r = expression.getVariablePointer("r");
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
        int element = 0;
        for (int i = 0; i < numPoints; i++) {
            if (i%3 != 1) {
                double x = (i+0.5)/numPoints;
                r[element++] = cutoff/x;
                if (element == width || i == numPoints-1) {
                    const float* result = expression.evaluate();
                    for (int j = 0; j < element; j++) {
                        float r2 = r[j]*r[j];
                        newSum += result[j]*r2*r2;
                    }
                    element = 0;
                }
            }
        }
        sum = newSum/numPoints + oldSum/3;
        double relativeChange = fabs((sum-oldSum)/sum);
        if (iteration > 2 && (relativeChange < 1e-5 || sum == 0))
            break;
        if (iteration == 10 || (iteration > 7 && relativeChange > 1e-3))
            throw OpenMMException("CustomNonbondedForce: Long range correction did not converge.  Does the energy go to 0 faster than 1/r^2?");
        numPoints *= 3;
    }

    // If a switching function is used, integrate over the switching interval.

    double sum2 = 0;
    if (force.getUseSwitchingFunction()) {
        double rswitch = force.getSwitchingDistance();
        sum2 = 0;
        numPoints = 1;
        vector<double> switchValue(width);
        for (int iteration = 0; ; iteration++) {
            double oldSum = sum2;
            double newSum = 0;
            int element = 0;
            for (int i = 0; i < numPoints; i++) {
                if (i%3 != 1) {
                    double x = (i+0.5)/numPoints;
                    switchValue[element] = x*x*x*(10+x*(-15+x*6));
                    r[element++] = rswitch+x*(cutoff-rswitch);
                    if (element == width || i == numPoints-1) {
                        const float* result = expression.evaluate();
                        for (int j = 0; j < element; j++)
                            newSum += switchValue[j]*result[j]*r[j]*r[j];
                        element = 0;
                    }
                }
            }
            sum2 = newSum/numPoints + oldSum/3;
            double relativeChange = fabs((sum2-oldSum)/sum2);
            if (iteration > 2 && (relativeChange < 1e-5 || sum2 == 0))
                break;
            if (iteration == 10 || (iteration > 7 && relativeChange > 1e-3))
                throw OpenMMException("CustomNonbondedForce: Long range correction did not converge.  Is the energy finite everywhere in the switching interval?");
            numPoints *= 3;
        }
        sum2 *= cutoff-rswitch;
    }
    return sum/cutoff+sum2;
}
