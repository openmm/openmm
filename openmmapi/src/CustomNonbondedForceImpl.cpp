/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
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
    return {CalcCustomNonbondedForceKernel::Name()};
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

/**
 * Index into a packed table holding the upper triangle of a numClasses x numClasses
 * matrix, for i <= j.
 */
static size_t triangularIndex(int i, int j, int numClasses) {
    return (size_t) i*numClasses - ((size_t) i*(i-1))/2 + (j-i);
}

static size_t triangularSize(int numClasses) {
    return ((size_t) numClasses*(numClasses+1))/2;
}

/**
 * Identify the particle classes and count the interactions between each pair of them.
 * This is the only part of the long range correction data that depends on the per-particle
 * parameters, so it is all that needs to be recomputed when they change.
 */
static void computeParticleClasses(const CustomNonbondedForce& force, CustomNonbondedForceImpl::LongRangeCorrectionData& data) {
    // Identify all particle classes (defined by parameters), and record the class of each particle.

    data.classes.clear();
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

    data.interactionCount.assign(triangularSize(numClasses), 0);
    
    if (force.getNumInteractionGroups() == 0) {
        // Count the particles of each class.
        
        vector<long long int> classCounts(numClasses, 0);
        for (int i = 0; i < numParticles; i++)
            classCounts[atomClass[i]]++;
        for (int i = 0; i < numClasses; i++) {
            data.interactionCount[triangularIndex(i, i, numClasses)] = (classCounts[i]*(classCounts[i]+1))/2;
            for (int j = i+1; j < numClasses; j++)
                data.interactionCount[triangularIndex(i, j, numClasses)] = classCounts[i]*classCounts[j];
        }
    }
    else {
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
                    data.interactionCount[triangularIndex(min(class1, class2), max(class1, class2), numClasses)]++;
                }
        }
    }
}

CustomNonbondedForceImpl::LongRangeCorrectionData CustomNonbondedForceImpl::prepareLongRangeCorrection(const CustomNonbondedForce& force, int numThreads) {
    LongRangeCorrectionData data;
    data.method = force.getNonbondedMethod();
    data.cutoffDistance = force.getCutoffDistance();
    data.switchingDistance = force.getSwitchingDistance();
    data.useSwitchingFunction = force.getUseSwitchingFunction();
    if (data.method == CustomNonbondedForce::NoCutoff || data.method == CustomNonbondedForce::CutoffNonPeriodic)
        return data;
        computeParticleClasses(force, data);

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
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        data.globalParameterNames.push_back(force.getGlobalParameterName(i));
    }
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        string name = force.getPerParticleParameterName(i);
        data.perParticleParameterNames.push_back(name);
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

    // Record the state of the tabulated functions, so we can tell whether the expressions
    // need to be recompiled later.

    for (int i = 0; i < force.getNumFunctions(); i++)
        data.tabulatedFunctionUpdateCount.push_back(force.getTabulatedFunction(i).getUpdateCount());
    return data;
}

void CustomNonbondedForceImpl::updateLongRangeCorrection(const CustomNonbondedForce& force, LongRangeCorrectionData& data, int numThreads) {
    // The compiled expressions depend only on the energy function and the tabulated functions,
    // neither of which can be changed by updateParametersInContext() (aside from the contents
    // of a tabulated function).  If they are still valid, we only need to recompute the classes.

    bool canReuse = (!data.energyExpression.empty() && data.method == force.getNonbondedMethod() &&
                     data.tabulatedFunctionUpdateCount.size() == force.getNumFunctions());
    for (int i = 0; canReuse && i < force.getNumFunctions(); i++)
        canReuse = (data.tabulatedFunctionUpdateCount[i] == force.getTabulatedFunction(i).getUpdateCount());
    if (!canReuse) {
        data = prepareLongRangeCorrection(force, numThreads);
        return;
    }
    computeParticleClasses(force, data);
}

/**
 * The most memory to spend on remembering integrals.  There is one table for the energy
 * and one for each energy parameter derivative, so the number of classes this allows for
 * depends on how many of those there are.  Beyond it the integrals are simply recomputed,
 * exactly as they were before the tables existed.
 */
static const size_t MAX_CACHE_BYTES = 64*1024*1024;

/**
 * Sum interactionCount*integral over all pairs of classes, reusing integrals computed on
 * the previous call.  An integral depends only on the parameters of the two classes and
 * the global parameters, so when per-particle parameters change, only the pairs involving
 * a class that was not present before need to be integrated again.
 *
 * oldIndex maps each current class to its index on the previous call, or -1 if it is new,
 * so a reusable integral costs one array lookup.
 */
double CustomNonbondedForceImpl::sumIntegrals(const function<Lepton::CompiledVectorExpression&(int)>& getExpression,
        const vector<int>& oldIndex, int numOldClasses, const vector<double>& oldIntegrals, vector<double>& newIntegrals,
        bool remember, LongRangeCorrectionData& data, const vector<vector<double> >& computedValues,
        const Context& context, ThreadPool& threads) {
    int numClasses = data.classes.size();
    bool reuse = (!oldIntegrals.empty() && numOldClasses > 0);

    // Collect the integrals we can reuse, and record which ones are missing.

    vector<pair<int, int> > pairs;
    vector<double> values;
    vector<int> missing;
    for (int i = 0; i < numClasses; i++)
        for (int j = i; j < numClasses; j++) {
            if (data.interactionCount[triangularIndex(i, j, numClasses)] == 0)
                continue;
            pairs.push_back(make_pair(i, j));
            values.push_back(0.0);
            if (reuse && oldIndex[i] >= 0 && oldIndex[j] >= 0) {
                int oldI = min(oldIndex[i], oldIndex[j]);
                int oldJ = max(oldIndex[i], oldIndex[j]);
                values.back() = oldIntegrals[triangularIndex(oldI, oldJ, numOldClasses)];
                continue;
            }
            missing.push_back(pairs.size()-1);
        }

    // Compute the missing integrals in parallel.  Each thread writes to its own elements
    // of values, so no synchronization is needed.

    atomic<int> atomicCounter(0);
    threads.execute([&] (ThreadPool& threads, int threadIndex) {
        Lepton::CompiledVectorExpression& expression = getExpression(threadIndex);
        while (true) {
            int k = atomicCounter++;
            if (k >= (int) missing.size())
                break;
            int i = pairs[missing[k]].first;
            int j = pairs[missing[k]].second;
            values[missing[k]] = integrateInteraction(expression, data.classes[i], data.classes[j],
                    computedValues[i], computedValues[j], data, context);
        }
    });
    threads.waitForThreads();

    // Record the integrals for the next call and add everything up.

    if (remember)
        newIntegrals.assign(triangularSize(numClasses), 0.0);
    else
        newIntegrals.clear();
    double sum = 0;
    for (int k = 0; k < (int) pairs.size(); k++) {
        size_t index = triangularIndex(pairs[k].first, pairs[k].second, numClasses);
        sum += data.interactionCount[index]*values[k];
        if (remember)
            newIntegrals[index] = values[k];
    }
    return sum;
}

void CustomNonbondedForceImpl::calcLongRangeCorrection(LongRangeCorrectionData& data, const Context& context, double& coefficient, vector<double>& derivatives, ThreadPool& threads) {
    if (data.method == CustomNonbondedForce::NoCutoff || data.method == CustomNonbondedForce::CutoffNonPeriodic) {
        coefficient = 0.0;
        return;
    }
    
    // Calculate the computed values for all atom classes.
    
    int numClasses = data.classes.size();
    vector<vector<double> > computedValues(numClasses, vector<double>(data.computedValueExpressions.size()));
    for (int i = 0; i < data.computedValueExpressions.size(); i++) {
        Lepton::CompiledExpression& expression = data.computedValueExpressions[i];
        const set<string>& variables = expression.getVariables();
        for (int j = 0; j < data.globalParameterNames.size(); j++) {
            const string& name = data.globalParameterNames[j];
            if (variables.find(name) != variables.end())
                expression.getVariableReference(name) = context.getParameter(name);
        }
        for (int j = 0; j < numClasses; j++) {
            for (int k = 0; k < data.perParticleParameterNames.size(); k++) {
                const string& name = data.perParticleParameterNames[k];
                if (variables.find(name) != variables.end())
                    expression.getVariableReference(name) = data.classes[j][k];
            }
            computedValues[j][i] = expression.evaluate();
        }
    }

    // Work out which classes were also present on the previous call.  Their integrals are
    // unchanged and can be reused.  This costs one lookup per class, not per pair.

    vector<double> globalValues;
    for (int i = 0; i < (int) data.globalParameterNames.size(); i++)
        globalValues.push_back(context.getParameter(data.globalParameterNames[i]));
    int numOldClasses = data.cachedClasses.size();
    vector<int> oldIndex(numClasses, -1);
    if (numOldClasses > 0 && data.cachedGlobalValues == globalValues) {
        map<vector<double>, int> oldClassIndex;
        for (int i = 0; i < numOldClasses; i++)
            oldClassIndex[data.cachedClasses[i]] = i;
        for (int i = 0; i < numClasses; i++) {
            map<vector<double>, int>::const_iterator entry = oldClassIndex.find(data.classes[i]);
            if (entry != oldClassIndex.end())
                oldIndex[i] = entry->second;
        }
    }

    // Decide whether to remember the integrals for the next call.  There is one table for
    // the energy and one for each parameter derivative, and they are held to a fixed
    // budget between them.

    int numDerivs = data.derivExpressions[0].size();
    size_t tableBytes = triangularSize(numClasses)*sizeof(double)*(numDerivs+1);
    bool remember = (tableBytes <= MAX_CACHE_BYTES);

    // Compute the coefficient.  The integrals are computed in parallel.

    double nPart = (double) context.getSystem().getNumParticles();
    double numInteractions = (nPart*(nPart+1))/2;
    vector<double> newIntegrals;
    double sum = sumIntegrals([&] (int threadIndex) -> Lepton::CompiledVectorExpression& { return data.energyExpression[threadIndex]; },
            oldIndex, numOldClasses, data.cachedIntegrals, newIntegrals, remember, data, computedValues, context, threads);
    sum /= numInteractions;
    coefficient = 2*M_PI*nPart*nPart*sum;

    // Now do the same for parameter derivatives.

    derivatives.resize(numDerivs);
    data.cachedDerivIntegrals.resize(numDerivs);
    vector<vector<double> > newDerivIntegrals(numDerivs);
    for (int k = 0; k < numDerivs; k++) {
        sum = sumIntegrals([&] (int threadIndex) -> Lepton::CompiledVectorExpression& { return data.derivExpressions[threadIndex][k]; },
                oldIndex, numOldClasses, data.cachedDerivIntegrals[k], newDerivIntegrals[k], remember, data, computedValues, context, threads);
        sum /= numInteractions;
        derivatives[k] = 2*M_PI*nPart*nPart*sum;
    }

    // Remember what we computed, for the next call.

    data.cachedClasses = data.classes;
    data.cachedGlobalValues = globalValues;
    data.cachedIntegrals.swap(newIntegrals);
    for (int k = 0; k < numDerivs; k++)
        data.cachedDerivIntegrals[k].swap(newDerivIntegrals[k]);
}


double CustomNonbondedForceImpl::integrateInteraction(Lepton::CompiledVectorExpression& expression, const vector<double>& params1, const vector<double>& params2,
        const vector<double>& computedValues1, const vector<double>& computedValues2, const LongRangeCorrectionData& data, const Context& context) {
    int width = expression.getWidth();
    const set<string>& variables = expression.getVariables();
    const vector<string>& paramNames = data.paramNames;
    for (int i = 0; i < data.perParticleParameterNames.size(); i++) {
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
    const vector<string>& computedValueNames = data.computedValueNames;
    for (int i = 0; i < data.computedValueExpressions.size(); i++) {
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
    for (int i = 0; i < data.globalParameterNames.size(); i++) {
        const string& name = data.globalParameterNames[i];
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
    double cutoff = data.cutoffDistance;
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
    if (data.useSwitchingFunction) {
        double rswitch = data.switchingDistance;
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
