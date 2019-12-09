/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
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

#include "openmm/common/CommonKernels.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Context.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCentroidBondForceImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "CommonKernelSources.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "SimTKOpenMMRealType.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>

using namespace OpenMM;
using namespace std;
using namespace Lepton;

static void setPeriodicBoxArgs(ComputeContext& cc, ComputeKernel kernel, int index) {
    Vec3 a, b, c;
    cc.getPeriodicBoxVectors(a, b, c);
    if (cc.getUseDoublePrecision()) {
        kernel->setArg(index++, mm_double4(a[0], b[1], c[2], 0.0));
        kernel->setArg(index++, mm_double4(1.0/a[0], 1.0/b[1], 1.0/c[2], 0.0));
        kernel->setArg(index++, mm_double4(a[0], a[1], a[2], 0.0));
        kernel->setArg(index++, mm_double4(b[0], b[1], b[2], 0.0));
        kernel->setArg(index, mm_double4(c[0], c[1], c[2], 0.0));
    }
    else {
        kernel->setArg(index++, mm_float4((float) a[0], (float) b[1], (float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4(1.0f/(float) a[0], 1.0f/(float) b[1], 1.0f/(float) c[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) a[0], (float) a[1], (float) a[2], 0.0f));
        kernel->setArg(index++, mm_float4((float) b[0], (float) b[1], (float) b[2], 0.0f));
        kernel->setArg(index, mm_float4((float) c[0], (float) c[1], (float) c[2], 0.0f));
    }
}

static bool isZeroExpression(const Lepton::ParsedExpression& expression) {
    const Lepton::Operation& op = expression.getRootNode().getOperation();
    if (op.getId() != Lepton::Operation::CONSTANT)
        return false;
    return (dynamic_cast<const Lepton::Operation::Constant&>(op).getValue() == 0.0);
}

static pair<ExpressionTreeNode, string> makeVariable(const string& name, const string& value) {
    return make_pair(ExpressionTreeNode(new Operation::Variable(name)), value);
}

class CommonCalcHarmonicBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const HarmonicBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(index, particle1, particle2, length, k);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double length1, length2, k1, k2;
        force.getBondParameters(group1, particle1, particle2, length1, k1);
        force.getBondParameters(group2, particle1, particle2, length2, k2);
        return (length1 == length2 && k1 == k2);
    }
private:
    const HarmonicBondForce& force;
};

void CommonCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params.initialize<mm_float2>(cc, numBonds, "bondParams");
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = mm_float2((float) length, (float) k);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::harmonicBondForce;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::bondForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(startIndex+i, atom1, atom2, length, k);
        paramVector[i] = mm_float2((float) length, (float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}
class CommonCalcCustomBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        vector<double> parameters;
        force.getBondParameters(index, particle1, particle2, parameters);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        vector<double> parameters1, parameters2;
        force.getBondParameters(group1, particle1, particle2, parameters1);
        force.getBondParameters(group2, particle1, particle2, parameters2);
        for (int i = 0; i < (int) parameters1.size(); i++)
            if (parameters1[i] != parameters2[i])
                return false;
        return true;
    }
private:
    const CustomBondForce& force;
};

CommonCalcCustomBondForceKernel::~CommonCalcCustomBondForceKernel() {
    cc.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CommonCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = new ComputeParameterSet(cc, force.getNumPerBondParameters(), numBonds, "customBondParams");
    vector<vector<float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["real dEdR = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["r"] = "r";
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customBondGlobals");
        globals.upload(globalParamValues);
        string argName = cc.getBondedUtilities().addArgument(globals, "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cc.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cc.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        string argName = cc.getBondedUtilities().addArgument(parameter.getArray(), parameter.getType());
        compute<<parameter.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cc.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::bondForce, replacements), force.getForceGroup());
}

double CommonCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    return 0.0;
}

void CommonCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<vector<float> > paramVector(numBonds);
    vector<double> parameters;
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        force.getBondParameters(startIndex+i, atom1, atom2, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcHarmonicAngleForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const HarmonicAngleForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(index, particle1, particle2, particle3, angle, k);
        particles.resize(3);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3;
        double angle1, angle2, k1, k2;
        force.getAngleParameters(group1, particle1, particle2, particle3, angle1, k1);
        force.getAngleParameters(group2, particle1, particle2, particle3, angle2, k2);
        return (angle1 == angle2 && k1 == k2);
    }
private:
    const HarmonicAngleForce& force;
};

void CommonCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params.initialize<mm_float2>(cc, numAngles, "angleParams");
    vector<mm_float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], angle, k);
        paramVector[i] = mm_float2((float) angle, (float) k);

    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::harmonicAngleForce;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::angleForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0)
        return;
    
    // Record the per-angle parameters.
    
    vector<mm_float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        int atom1, atom2, atom3;
        double angle, k;
        force.getAngleParameters(startIndex+i, atom1, atom2, atom3, angle, k);
        paramVector[i] = mm_float2((float) angle, (float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

class CommonCalcCustomAngleForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomAngleForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3;
        vector<double> parameters;
        force.getAngleParameters(index, particle1, particle2, particle3, parameters);
        particles.resize(3);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3;
        vector<double> parameters1, parameters2;
        force.getAngleParameters(group1, particle1, particle2, particle3, parameters1);
        force.getAngleParameters(group2, particle1, particle2, particle3, parameters2);
        for (int i = 0; i < (int) parameters1.size(); i++)
            if (parameters1[i] != parameters2[i])
                return false;
        return true;
    }
private:
    const CustomAngleForce& force;
};

CommonCalcCustomAngleForceKernel::~CommonCalcCustomAngleForceKernel() {
    cc.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CommonCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = new ComputeParameterSet(cc, force.getNumPerAngleParameters(), numAngles, "customAngleParams");
    vector<vector<float> > paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        vector<double> parameters;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("theta").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["real dEdAngle = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["theta"] = "theta";
    for (int i = 0; i < force.getNumPerAngleParameters(); i++) {
        const string& name = force.getPerAngleParameterName(i);
        variables[name] = "angleParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customAngleGlobals");
        globals.upload(globalParamValues);
        string argName = cc.getBondedUtilities().addArgument(globals, "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cc.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cc.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        string argName = cc.getBondedUtilities().addArgument(parameter.getArray(), parameter.getType());
        compute<<parameter.getType()<<" angleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cc.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::angleForce, replacements), force.getForceGroup());
}

double CommonCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    return 0.0;
}

void CommonCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0)
        return;
    
    // Record the per-angle parameters.
    
    vector<vector<float> > paramVector(numAngles);
    vector<double> parameters;
    for (int i = 0; i < numAngles; i++) {
        int atom1, atom2, atom3;
        force.getAngleParameters(startIndex+i, atom1, atom2, atom3, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcPeriodicTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const PeriodicTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, periodicity, phase, k);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4, periodicity1, periodicity2;
        double phase1, phase2, k1, k2;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, periodicity1, phase1, k1);
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, periodicity2, phase2, k2);
        return (periodicity1 == periodicity2 && phase1 == phase2 && k1 == k2);
    }
private:
    const PeriodicTorsionForce& force;
};

void CommonCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params.initialize<mm_float4>(cc, numTorsions, "periodicTorsionParams");
    vector<mm_float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], periodicity, phase, k);
        paramVector[i] = mm_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::periodicTorsionForce;
    replacements["PARAMS"] = cc.getBondedUtilities().addArgument(params, "float4");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<mm_float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4, periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, periodicity, phase, k);
        paramVector[i] = mm_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

class CommonCalcRBTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const RBTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        double c0a, c0b, c1a, c1b, c2a, c2b, c3a, c3b, c4a, c4b, c5a, c5b;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, c0a, c1a, c2a, c3a, c4a, c5a);
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, c0b, c1b, c2b, c3b, c4b, c5b);
        return (c0a == c0b && c1a == c1b && c2a == c2b && c3a == c3b && c4a == c4b && c5a == c5b);
    }
private:
    const RBTorsionForce& force;
};

void CommonCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params1.initialize<mm_float4>(cc, numTorsions, "rbTorsionParams1");
    params2.initialize<mm_float2>(cc, numTorsions, "rbTorsionParams2");
    vector<mm_float4> paramVector1(numTorsions);
    vector<mm_float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], c0, c1, c2, c3, c4, c5);
        paramVector1[i] = mm_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = mm_float2((float) c4, (float) c5);

    }
    params1.upload(paramVector1);
    params2.upload(paramVector2);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CommonKernelSources::rbTorsionForce;
    replacements["PARAMS1"] = cc.getBondedUtilities().addArgument(params1, "float4");
    replacements["PARAMS2"] = cc.getBondedUtilities().addArgument(params2, "float2");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<mm_float4> paramVector1(numTorsions);
    vector<mm_float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        paramVector1[i] = mm_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = mm_float2((float) c4, (float) c5);
    }
    params1.upload(paramVector1);
    params2.upload(paramVector2);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules();
}

class CommonCalcCustomTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2, particle3, particle4;
        vector<double> parameters;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, parameters);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        vector<double> parameters1, parameters2;
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, parameters1);
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, parameters2);
        for (int i = 0; i < (int) parameters1.size(); i++)
            if (parameters1[i] != parameters2[i])
                return false;
        return true;
    }
private:
    const CustomTorsionForce& force;
};

CommonCalcCustomTorsionForceKernel::~CommonCalcCustomTorsionForceKernel() {
    if (params != NULL)
        delete params;
}

void CommonCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params = new ComputeParameterSet(cc, force.getNumPerTorsionParameters(), numTorsions, "customTorsionParams");
    vector<vector<float> > paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        vector<double> parameters;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("theta").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["real dEdAngle = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["theta"] = "theta";
    for (int i = 0; i < force.getNumPerTorsionParameters(); i++) {
        const string& name = force.getPerTorsionParameterName(i);
        variables[name] = "torsionParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customTorsionGlobals");
        globals.upload(globalParamValues);
        string argName = cc.getBondedUtilities().addArgument(globals, "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cc.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cc.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        string argName = cc.getBondedUtilities().addArgument(parameter.getArray(), parameter.getType());
        compute<<parameter.getType()<<" torsionParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cc.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::torsionForce, replacements), force.getForceGroup());
}

double CommonCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    return 0.0;
}

void CommonCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<vector<float> > paramVector(numTorsions);
    vector<double> parameters;
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcCMAPTorsionForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CMAPTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int map, a1, a2, a3, a4, b1, b2, b3, b4;
        force.getTorsionParameters(index, map, a1, a2, a3, a4, b1, b2, b3, b4);
        particles.resize(8);
        particles[0] = a1;
        particles[1] = a2;
        particles[2] = a3;
        particles[3] = a4;
        particles[4] = b1;
        particles[5] = b2;
        particles[6] = b3;
        particles[7] = b4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int map1, map2, a1, a2, a3, a4, b1, b2, b3, b4;
        force.getTorsionParameters(group1, map1, a1, a2, a3, a4, b1, b2, b3, b4);
        force.getTorsionParameters(group2, map2, a1, a2, a3, a4, b1, b2, b3, b4);
        return (map1 == map2);
    }
private:
    const CMAPTorsionForce& force;
};

void CommonCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    int numMaps = force.getNumMaps();
    vector<mm_float4> coeffVec;
    mapPositionsVec.resize(numMaps);
    vector<double> energy;
    vector<vector<double> > c;
    int currentPosition = 0;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        mapPositionsVec[i] = mm_int2(currentPosition, size);
        currentPosition += 4*size*size;
        for (int j = 0; j < size*size; j++) {
            coeffVec.push_back(mm_float4((float) c[j][0], (float) c[j][1], (float) c[j][2], (float) c[j][3]));
            coeffVec.push_back(mm_float4((float) c[j][4], (float) c[j][5], (float) c[j][6], (float) c[j][7]));
            coeffVec.push_back(mm_float4((float) c[j][8], (float) c[j][9], (float) c[j][10], (float) c[j][11]));
            coeffVec.push_back(mm_float4((float) c[j][12], (float) c[j][13], (float) c[j][14], (float) c[j][15]));
        }
    }
    vector<vector<int> > atoms(numTorsions, vector<int>(8));
    vector<int> torsionMapsVec(numTorsions);
    for (int i = 0; i < numTorsions; i++)
        force.getTorsionParameters(startIndex+i, torsionMapsVec[i], atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], atoms[i][6], atoms[i][7]);
    coefficients.initialize<mm_float4>(cc, coeffVec.size(), "cmapTorsionCoefficients");
    mapPositions.initialize<mm_int2>(cc, numMaps, "cmapTorsionMapPositions");
    torsionMaps.initialize<int>(cc, numTorsions, "cmapTorsionMaps");
    coefficients.upload(coeffVec);
    mapPositions.upload(mapPositionsVec);
    torsionMaps.upload(torsionMapsVec);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COEFF"] = cc.getBondedUtilities().addArgument(coefficients, "float4");
    replacements["MAP_POS"] = cc.getBondedUtilities().addArgument(mapPositions, "int2");
    replacements["MAPS"] = cc.getBondedUtilities().addArgument(torsionMaps, "int");
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::cmapTorsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CommonCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    int numMaps = force.getNumMaps();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (mapPositions.getSize() != numMaps)
        throw OpenMMException("updateParametersInContext: The number of maps has changed");
    if (torsionMaps.getSize() != numTorsions)
        throw OpenMMException("updateParametersInContext: The number of CMAP torsions has changed");

    // Update the maps.

    vector<mm_float4> coeffVec;
    vector<double> energy;
    vector<vector<double> > c;
    int currentPosition = 0;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        if (size != mapPositionsVec[i].y)
            throw OpenMMException("updateParametersInContext: The size of a map has changed");
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        currentPosition += 4*size*size;
        for (int j = 0; j < size*size; j++) {
            coeffVec.push_back(mm_float4((float) c[j][0], (float) c[j][1], (float) c[j][2], (float) c[j][3]));
            coeffVec.push_back(mm_float4((float) c[j][4], (float) c[j][5], (float) c[j][6], (float) c[j][7]));
            coeffVec.push_back(mm_float4((float) c[j][8], (float) c[j][9], (float) c[j][10], (float) c[j][11]));
            coeffVec.push_back(mm_float4((float) c[j][12], (float) c[j][13], (float) c[j][14], (float) c[j][15]));
        }
    }
    coefficients.upload(coeffVec);

    // Update the indices.

    vector<int> torsionMapsVec(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int index[8];
        force.getTorsionParameters(i, torsionMapsVec[i], index[0], index[1], index[2], index[3], index[4], index[5], index[6], index[7]);
    }
    torsionMaps.upload(torsionMapsVec);
}
class CommonCalcCustomExternalForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomExternalForce& force, int numParticles) : force(force), indices(numParticles, -1) {
        vector<double> params;
        for (int i = 0; i < force.getNumParticles(); i++) {
            int particle;
            force.getParticleParameters(i, particle, params);
            indices[particle] = i;
        }
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        particle1 = indices[particle1];
        particle2 = indices[particle2];
        if (particle1 == -1 && particle2 == -1)
            return true;
        if (particle1 == -1 || particle2 == -1)
            return false;
        int temp;
        vector<double> params1;
        vector<double> params2;
        force.getParticleParameters(particle1, temp, params1);
        force.getParticleParameters(particle2, temp, params2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        return true;
    }
private:
    const CustomExternalForce& force;
    vector<int> indices;
};

CommonCalcCustomExternalForceKernel::~CommonCalcCustomExternalForceKernel() {
    cc.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CommonCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumParticles()/numContexts;
    numParticles = endIndex-startIndex;
    if (numParticles == 0)
        return;
    vector<vector<int> > atoms(numParticles, vector<int>(1));
    params = new ComputeParameterSet(cc, force.getNumPerParticleParameters(), numParticles, "customExternalParams");
    vector<vector<float> > paramVector(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(startIndex+i, atoms[i][0], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force, system.getNumParticles());
    cc.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    map<string, Lepton::CustomFunction*> customFunctions;
    customFunctions["periodicdistance"] = cc.getExpressionUtilities().getPeriodicDistancePlaceholder();
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), customFunctions).optimize();
    Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x").optimize();
    Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y").optimize();
    Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["real dEdX = "] = forceExpressionX;
    expressions["real dEdY = "] = forceExpressionY;
    expressions["real dEdZ = "] = forceExpressionZ;

    // Create the kernels.

    map<string, string> variables;
    variables["x"] = "pos1.x";
    variables["y"] = "pos1.y";
    variables["z"] = "pos1.z";
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        variables[name] = "particleParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customExternalGlobals");
        globals.upload(globalParamValues);
        string argName = cc.getBondedUtilities().addArgument(globals, "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cc.intToString(i)+"]";
            variables[name] = value;
        }
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        string argName = cc.getBondedUtilities().addArgument(parameter.getArray(), parameter.getType());
        compute<<parameter.getType()<<" particleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cc.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    cc.getBondedUtilities().addInteraction(atoms, cc.replaceStrings(CommonKernelSources::customExternalForce, replacements), force.getForceGroup());
}

double CommonCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    return 0.0;
}

void CommonCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumParticles()/numContexts;
    if (numParticles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    if (numParticles == 0)
        return;
    
    // Record the per-particle parameters.
    
    vector<vector<float> > paramVector(numParticles);
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        int particle;
        force.getParticleParameters(startIndex+i, particle, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcCustomCompoundBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomCompoundBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        vector<double> parameters;
        force.getBondParameters(index, particles, parameters);
    }
    bool areGroupsIdentical(int group1, int group2) {
        vector<int> particles;
        vector<double> parameters1, parameters2;
        force.getBondParameters(group1, particles, parameters1);
        force.getBondParameters(group2, particles, parameters2);
        for (int i = 0; i < (int) parameters1.size(); i++)
            if (parameters1[i] != parameters2[i])
                return false;
        return true;
    }
private:
    const CustomCompoundBondForce& force;
};

CommonCalcCustomCompoundBondForceKernel::~CommonCalcCustomCompoundBondForceKernel() {
    cc.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CommonCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    int particlesPerBond = force.getNumParticlesPerBond();
    vector<vector<int> > atoms(numBonds, vector<int>(particlesPerBond));
    params = new ComputeParameterSet(cc, force.getNumPerBondParameters(), numBonds, "customCompoundBondParams");
    vector<vector<float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        force.getBondParameters(startIndex+i, atoms[i], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        string arrayName = cc.getBondedUtilities().addArgument(tabulatedFunctions[i], width == 1 ? "float" : "float"+cc.intToString(width));
        functionDefinitions.push_back(make_pair(name, arrayName));
    }
    
    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    map<string, string> variables;
    for (int i = 0; i < particlesPerBond; i++) {
        string index = cc.intToString(i+1);
        variables["x"+index] = "pos"+index+".x";
        variables["y"+index] = "pos"+index+".y";
        variables["z"+index] = "pos"+index+".z";
    }
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customCompoundBondGlobals");
        globals.upload(globalParamValues);
        string argName = cc.getBondedUtilities().addArgument(globals, "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cc.intToString(i)+"]";
            variables[name] = value;
        }
    }

    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
    // and dihedrals the expression depends on.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomCompoundBondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    map<string, Lepton::ParsedExpression> forceExpressions;
    set<string> computedDeltas;
    vector<string> atomNames, posNames;
    for (int i = 0; i < particlesPerBond; i++) {
        string index = cc.intToString(i+1);
        atomNames.push_back("P"+index);
        posNames.push_back("pos"+index);
    }
    stringstream compute;
    int index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = SQRT(delta"<<deltaName<<".w);\n";
        variables[distance.first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cc.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = ccb_computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        variables[angle.first] = angleName;
        forceExpressions["real dEdAngle"+cc.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& atoms = dihedral.second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        string dihedralName = "dihedral_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]]+atomNames[atoms[3]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            compute<<"real4 delta"<<deltaName3<<" = ccb_delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[3]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName3);
        }
        compute<<"real4 "<<crossName1<<" = ccb_computeCross(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        compute<<"real4 "<<crossName2<<" = ccb_computeCross(delta"<<deltaName2<<", delta"<<deltaName3<<");\n";
        compute<<"real "<<dihedralName<<" = ccb_computeAngle("<<crossName1<<", "<<crossName2<<");\n";
        compute<<dihedralName<<" *= (delta"<<deltaName1<<".x*"<<crossName2<<".x + delta"<<deltaName1<<".y*"<<crossName2<<".y + delta"<<deltaName1<<".z*"<<crossName2<<".z < 0 ? -1 : 1);\n";
        variables[dihedral.first] = dihedralName;
        forceExpressions["real dEdDihedral"+cc.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        string argName = cc.getBondedUtilities().addArgument(parameter.getArray(), parameter.getType());
        compute<<parameter.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cc.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        forceExpressions[derivVariable+" += "] = derivExpression;
    }
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    vector<string> forceNames;
    for (int i = 0; i < particlesPerBond; i++) {
        string istr = cc.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = make_real3(0);\n";
        compute<<"{\n";
        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
        map<string, Lepton::ParsedExpression> expressions;
        if (!isZeroExpression(forceExpressionX))
            expressions[forceName+".x -= "] = forceExpressionX;
        if (!isZeroExpression(forceExpressionY))
            expressions[forceName+".y -= "] = forceExpressionY;
        if (!isZeroExpression(forceExpressionZ))
            expressions[forceName+".z -= "] = forceExpressionZ;
        if (expressions.size() > 0)
            compute<<cc.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cc.intToString(index)+"/r_"+deltaName+")*trimTo3(delta"+deltaName+")";
        compute<<forceNames[atoms[0]]<<" += "<<"-"<<value<<";\n";
        compute<<forceNames[atoms[1]]<<" += "<<value<<";\n";
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute<<"{\n";
        compute<<"real3 crossProd = trimTo3(cross(delta"<<deltaName2<<", delta"<<deltaName1<<"));\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), (real) 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(trimTo3(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cc.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(trimTo3(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cc.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real3 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[atoms[0]]<<" += deltaCross0;\n";
        compute<<forceNames[atoms[1]]<<" += deltaCross1;\n";
        compute<<forceNames[atoms[2]]<<" += deltaCross2;\n";
        compute<<"}\n";
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& atoms = dihedral.second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        compute<<"{\n";
        compute<<"real r = SQRT(delta"<<deltaName2<<".w);\n";
        compute<<"real4 ff;\n";
        compute<<"ff.x = (-dEdDihedral"<<cc.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cc.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real3 internalF0 = ff.x*trimTo3("<<crossName1<<");\n";
        compute<<"real3 internalF3 = ff.w*trimTo3("<<crossName2<<");\n";
        compute<<"real3 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[atoms[0]]<<" += internalF0;\n";
        compute<<forceNames[atoms[1]]<<" += s-internalF0;\n";
        compute<<forceNames[atoms[2]]<<" += -s-internalF3;\n";
        compute<<forceNames[atoms[3]]<<" += internalF3;\n";
        compute<<"}\n";
        index++;
    }
    cc.getBondedUtilities().addInteraction(atoms, compute.str(), force.getForceGroup());
    map<string, string> replacements;
    replacements["M_PI"] = cc.doubleToString(M_PI);
    cc.getBondedUtilities().addPrefixCode(cc.replaceStrings(CommonKernelSources::customCompoundBond, replacements));
}

double CommonCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    return 0.0;
}

void CommonCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    cc.setAsCurrent();
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<vector<float> > paramVector(numBonds);
    vector<int> particles;
    vector<double> parameters;
    for (int i = 0; i < numBonds; i++) {
        force.getBondParameters(startIndex+i, particles, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcCustomCentroidBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomCentroidBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        vector<double> parameters;
        vector<int> groups;
        force.getBondParameters(index, groups, parameters);
        for (int group : groups) {
            vector<int> groupParticles;
            vector<double> weights;
            force.getGroupParameters(group, groupParticles, weights);
            particles.insert(particles.end(), groupParticles.begin(), groupParticles.end());
        }
    }
    bool areGroupsIdentical(int group1, int group2) {
        vector<int> groups1, groups2;
        vector<double> parameters1, parameters2;
        force.getBondParameters(group1, groups1, parameters1);
        force.getBondParameters(group2, groups2, parameters2);
        for (int i = 0; i < (int) parameters1.size(); i++)
            if (parameters1[i] != parameters2[i])
                return false;
        for (int i = 0; i < groups1.size(); i++) {
            vector<int> groupParticles;
            vector<double> weights1, weights2;
            force.getGroupParameters(groups1[i], groupParticles, weights1);
            force.getGroupParameters(groups2[i], groupParticles, weights2);
            if (weights1.size() != weights2.size())
                return false;
            for (int j = 0; j < weights1.size(); j++)
                if (weights1[j] != weights2[j])
                    return false;
        }
        return true;
    }
private:
    const CustomCentroidBondForce& force;
};

CommonCalcCustomCentroidBondForceKernel::~CommonCalcCustomCentroidBondForceKernel() {
    cc.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CommonCalcCustomCentroidBondForceKernel::initialize(const System& system, const CustomCentroidBondForce& force) {
    cc.setAsCurrent();
    numBonds = force.getNumBonds();
    if (numBonds == 0)
        return;
    if (!cc.getSupports64BitGlobalAtomics())
        throw OpenMMException("CustomCentroidBondForce requires a device that supports 64 bit atomic operations");
    info = new ForceInfo(force);
    cc.addForce(info);
    
    // Record the groups.
    
    numGroups = force.getNumGroups();
    vector<int> groupParticleVec;
    vector<double> groupWeightVec;
    vector<int> groupOffsetVec;
    groupOffsetVec.push_back(0);
    for (int i = 0; i < numGroups; i++) {
        vector<int> particles;
        vector<double> weights;
        force.getGroupParameters(i, particles, weights);
        groupParticleVec.insert(groupParticleVec.end(), particles.begin(), particles.end());
        groupOffsetVec.push_back(groupParticleVec.size());
    }
    vector<vector<double> > normalizedWeights;
    CustomCentroidBondForceImpl::computeNormalizedWeights(force, system, normalizedWeights);
    for (int i = 0; i < numGroups; i++)
        groupWeightVec.insert(groupWeightVec.end(), normalizedWeights[i].begin(), normalizedWeights[i].end());
    groupParticles.initialize<int>(cc, groupParticleVec.size(), "groupParticles");
    groupParticles.upload(groupParticleVec);
    if (cc.getUseDoublePrecision()) {
        groupWeights.initialize<double>(cc, groupParticleVec.size(), "groupWeights");
        centerPositions.initialize<mm_double4>(cc, numGroups, "centerPositions");
    }
    else {
        groupWeights.initialize<float>(cc, groupParticleVec.size(), "groupWeights");
        centerPositions.initialize<mm_float4>(cc, numGroups, "centerPositions");
    }
    groupWeights.upload(groupWeightVec, true);
    groupOffsets.initialize<int>(cc, groupOffsetVec.size(), "groupOffsets");
    groupOffsets.upload(groupOffsetVec);
    groupForces.initialize<long long>(cc, numGroups*3, "groupForces");
    cc.addAutoclearBuffer(groupForces);
    
    // Record the bonds.
    
    int groupsPerBond = force.getNumGroupsPerBond();
    vector<int> bondGroupVec(numBonds*groupsPerBond);
    params = new ComputeParameterSet(cc, force.getNumPerBondParameters(), numBonds, "customCentroidBondParams");
    vector<vector<float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<int> groups;
        vector<double> parameters;
        force.getBondParameters(i, groups, parameters);
        for (int j = 0; j < groups.size(); j++)
            bondGroupVec[i+j*numBonds] = groups[j];
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    bondGroups.initialize<int>(cc, bondGroupVec.size(), "bondGroups");
    bondGroups.upload(bondGroupVec);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream extraArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cc.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        extraArgs << ", GLOBAL const float";
        if (width > 1)
            extraArgs << width;
        extraArgs << "* RESTRICT " << arrayName;
    }
    
    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    map<string, string> variables;
    for (int i = 0; i < groupsPerBond; i++) {
        string index = cc.intToString(i+1);
        variables["x"+index] = "pos"+index+".x";
        variables["y"+index] = "pos"+index+".y";
        variables["z"+index] = "pos"+index+".z";
    }
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    needEnergyParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    if (needEnergyParamDerivs)
        extraArgs << ", GLOBAL mixed* RESTRICT energyParamDerivs";
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customCentroidBondGlobals");
        globals.upload(globalParamValues);
        extraArgs << ", GLOBAL const float* RESTRICT globals";
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cc.intToString(i)+"]";
            variables[name] = value;
        }
    }

    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
    // and dihedrals the expression depends on.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomCentroidBondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    map<string, Lepton::ParsedExpression> forceExpressions;
    set<string> computedDeltas;
    vector<string> atomNames, posNames;
    for (int i = 0; i < groupsPerBond; i++) {
        string index = cc.intToString(i+1);
        atomNames.push_back("P"+index);
        posNames.push_back("pos"+index);
    }
    stringstream compute, initParamDerivs, saveParamDerivs;
    for (int i = 0; i < groupsPerBond; i++) {
        compute<<"int group"<<(i+1)<<" = bondGroups[index+"<<(i*numBonds)<<"];\n";
        compute<<"real4 pos"<<(i+1)<<" = centerPositions[group"<<(i+1)<<"];\n";
    }
    int index = 0;
    for (auto& distance : distances) {
        const vector<int>& groups = distance.second;
        string deltaName = atomNames[groups[0]]+atomNames[groups[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = delta("<<posNames[groups[0]]<<", "<<posNames[groups[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
        variables[distance.first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cc.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& groups = angle.second;
        string deltaName1 = atomNames[groups[1]]+atomNames[groups[0]];
        string deltaName2 = atomNames[groups[1]]+atomNames[groups[2]];
        string angleName = "angle_"+atomNames[groups[0]]+atomNames[groups[1]]+atomNames[groups[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[groups[1]]<<", "<<posNames[groups[0]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[groups[1]]<<", "<<posNames[groups[2]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        variables[angle.first] = angleName;
        forceExpressions["real dEdAngle"+cc.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& groups = dihedral.second;
        string deltaName1 = atomNames[groups[0]]+atomNames[groups[1]];
        string deltaName2 = atomNames[groups[2]]+atomNames[groups[1]];
        string deltaName3 = atomNames[groups[2]]+atomNames[groups[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        string dihedralName = "dihedral_"+atomNames[groups[0]]+atomNames[groups[1]]+atomNames[groups[2]]+atomNames[groups[3]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[groups[0]]<<", "<<posNames[groups[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[groups[2]]<<", "<<posNames[groups[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            compute<<"real4 delta"<<deltaName3<<" = delta("<<posNames[groups[2]]<<", "<<posNames[groups[3]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName3);
        }
        compute<<"real4 "<<crossName1<<" = computeCross(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        compute<<"real4 "<<crossName2<<" = computeCross(delta"<<deltaName2<<", delta"<<deltaName3<<");\n";
        compute<<"real "<<dihedralName<<" = computeAngle("<<crossName1<<", "<<crossName2<<");\n";
        compute<<dihedralName<<" *= (delta"<<deltaName1<<".x*"<<crossName2<<".x + delta"<<deltaName1<<".y*"<<crossName2<<".y + delta"<<deltaName1<<".z*"<<crossName2<<".z < 0 ? -1 : 1);\n";
        variables[dihedral.first] = dihedralName;
        forceExpressions["real dEdDihedral"+cc.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        extraArgs<<", GLOBAL const "<<parameter.getType()<<"* RESTRICT globalParams"<<i;
        compute<<parameter.getType()<<" bondParams"<<(i+1)<<" = globalParams"<<i<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    if (needEnergyParamDerivs) {
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string paramName = force.getEnergyParameterDerivativeName(i);
            cc.addEnergyParameterDerivative(paramName);
            Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
            forceExpressions[string("energyParamDeriv")+cc.intToString(i)+" += "] = derivExpression;
            initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
        }
        const vector<string>& allParamDerivNames = cc.getEnergyParamDerivNames();
        int numDerivs = allParamDerivNames.size();
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
            for (int index = 0; index < numDerivs; index++)
                if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                    saveParamDerivs << "energyParamDerivs[GLOBAL_ID*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
    }
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to groups.

    vector<string> forceNames;
    for (int i = 0; i < groupsPerBond; i++) {
        string istr = cc.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = make_real3(0);\n";
        compute<<"{\n";
        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
        map<string, Lepton::ParsedExpression> expressions;
        if (!isZeroExpression(forceExpressionX))
            expressions[forceName+".x -= "] = forceExpressionX;
        if (!isZeroExpression(forceExpressionY))
            expressions[forceName+".y -= "] = forceExpressionY;
        if (!isZeroExpression(forceExpressionZ))
            expressions[forceName+".z -= "] = forceExpressionZ;
        if (expressions.size() > 0)
            compute<<cc.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (auto& distance : distances) {
        const vector<int>& groups = distance.second;
        string deltaName = atomNames[groups[0]]+atomNames[groups[1]];
        string value = "(dEdDistance"+cc.intToString(index)+"/r_"+deltaName+")*trimTo3(delta"+deltaName+")";
        compute<<forceNames[groups[0]]<<" += "<<"-"<<value<<";\n";
        compute<<forceNames[groups[1]]<<" += "<<value<<";\n";
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& groups = angle.second;
        string deltaName1 = atomNames[groups[1]]+atomNames[groups[0]];
        string deltaName2 = atomNames[groups[1]]+atomNames[groups[2]];
        compute<<"{\n";
        compute<<"real3 crossProd = trimTo3(cross(delta"<<deltaName2<<", delta"<<deltaName1<<"));\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), (real) 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(trimTo3(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cc.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(trimTo3(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cc.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real3 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[groups[0]]<<" += deltaCross0;\n";
        compute<<forceNames[groups[1]]<<" += deltaCross1;\n";
        compute<<forceNames[groups[2]]<<" += deltaCross2;\n";
        compute<<"}\n";
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& groups = dihedral.second;
        string deltaName1 = atomNames[groups[0]]+atomNames[groups[1]];
        string deltaName2 = atomNames[groups[2]]+atomNames[groups[1]];
        string deltaName3 = atomNames[groups[2]]+atomNames[groups[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        compute<<"{\n";
        compute<<"real r = sqrt(delta"<<deltaName2<<".w);\n";
        compute<<"real4 ff;\n";
        compute<<"ff.x = (-dEdDihedral"<<cc.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cc.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real3 internalF0 = ff.x*trimTo3("<<crossName1<<");\n";
        compute<<"real3 internalF3 = ff.w*trimTo3("<<crossName2<<");\n";
        compute<<"real3 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[groups[0]]<<" += internalF0;\n";
        compute<<forceNames[groups[1]]<<" += s-internalF0;\n";
        compute<<forceNames[groups[2]]<<" += -s-internalF3;\n";
        compute<<forceNames[groups[3]]<<" += internalF3;\n";
        compute<<"}\n";
        index++;
    }
    
    // Save the forces to global memory.
    
    for (int i = 0; i < groupsPerBond; i++) {
        compute<<"ATOMIC_ADD(&groupForce[group"<<(i+1)<<"], (mm_ulong) ((mm_long) (force"<<(i+1)<<".x*0x100000000)));\n";
        compute<<"ATOMIC_ADD(&groupForce[group"<<(i+1)<<"+numParticleGroups], (mm_ulong) ((mm_long) (force"<<(i+1)<<".y*0x100000000)));\n";
        compute<<"ATOMIC_ADD(&groupForce[group"<<(i+1)<<"+numParticleGroups*2], (mm_ulong) ((mm_long) (force"<<(i+1)<<".z*0x100000000)));\n";
        compute<<"MEM_FENCE;\n";
    }
    map<string, string> replacements;
    replacements["M_PI"] = cc.doubleToString(M_PI);
    replacements["NUM_BONDS"] = cc.intToString(numBonds);
    replacements["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    replacements["EXTRA_ARGS"] = extraArgs.str();
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
    replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customCentroidBond, replacements));
    computeCentersKernel = program->createKernel("computeGroupCenters");
    computeCentersKernel->addArg(numGroups);
    computeCentersKernel->addArg(cc.getPosq());
    computeCentersKernel->addArg(groupParticles);
    computeCentersKernel->addArg(groupWeights);
    computeCentersKernel->addArg(groupOffsets);
    computeCentersKernel->addArg(centerPositions);
    groupForcesKernel = program->createKernel("computeGroupForces");
    groupForcesKernel->addArg(numGroups);
    groupForcesKernel->addArg(groupForces);
    groupForcesKernel->addArg(); // Energy buffer hasn't been created yet
    groupForcesKernel->addArg(centerPositions);
    groupForcesKernel->addArg(bondGroups);
    for (int i = 0; i < 5; i++)
        groupForcesKernel->addArg(); // Periodic box information will be set just before it is executed.
    if (needEnergyParamDerivs)
        groupForcesKernel->addArg(); // Deriv buffer hasn't been created yet.
    for (auto& function : tabulatedFunctions)
        groupForcesKernel->addArg(function);
    if (globals.isInitialized())
        groupForcesKernel->addArg(globals);
    for (auto& parameter : params->getParameterInfos())
        groupForcesKernel->addArg(parameter.getArray());
    applyForcesKernel = program->createKernel("applyForcesToAtoms");
    applyForcesKernel->addArg(numGroups);
    applyForcesKernel->addArg(groupParticles);
    applyForcesKernel->addArg(groupWeights);
    applyForcesKernel->addArg(groupOffsets);
    applyForcesKernel->addArg(groupForces);
    applyForcesKernel->addArg();
}

double CommonCalcCustomCentroidBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (numBonds == 0)
        return 0.0;
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    computeCentersKernel->execute(32*numGroups);
    groupForcesKernel->setArg(2, cc.getEnergyBuffer());
    setPeriodicBoxArgs(cc, groupForcesKernel, 5);
    if (needEnergyParamDerivs)
        groupForcesKernel->setArg(10, cc.getEnergyParamDerivBuffer());
    groupForcesKernel->execute(numBonds);
    applyForcesKernel->setArg(5, cc.getLongForceBuffer());
    applyForcesKernel->execute(32*numGroups);
    return 0.0;
}

void CommonCalcCustomCentroidBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force) {
    cc.setAsCurrent();
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<vector<float> > paramVector(numBonds);
    vector<int> particles;
    vector<double> parameters;
    for (int i = 0; i < numBonds; i++) {
        force.getBondParameters(i, particles, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcCustomManyParticleForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomManyParticleForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        vector<double> params1, params2;
        int type1, type2;
        force.getParticleParameters(particle1, params1, type1);
        force.getParticleParameters(particle2, params2, type2);
        if (type1 != type2)
            return false;
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExclusions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        force.getExclusionParticles(index, particle1, particle2);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        return true;
    }
private:
    const CustomManyParticleForce& force;
};

CommonCalcCustomManyParticleForceKernel::~CommonCalcCustomManyParticleForceKernel() {
    cc.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CommonCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {
    cc.setAsCurrent();
    if (!cc.getSupports64BitGlobalAtomics())
        throw OpenMMException("CustomManyParticleForce requires a device that supports 64 bit atomic operations");
    int numParticles = force.getNumParticles();
    int particlesPerSet = force.getNumParticlesPerSet();
    bool centralParticleMode = (force.getPermutationMode() == CustomManyParticleForce::UniqueCentralParticle);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    forceWorkgroupSize = 128;
    findNeighborsWorkgroupSize = (cc.getSIMDWidth() >= 32 ? 128 : 32);
    
    // Record parameter values.
    
    params = new ComputeParameterSet(cc, force.getNumPerParticleParameters(), numParticles, "customManyParticleParameters");
    vector<vector<float> > paramVector(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        int type;
        force.getParticleParameters(i, parameters, type);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cc.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        tableArgs << ", GLOBAL const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* RESTRICT " << arrayName;
    }
    
    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    vector<pair<ExpressionTreeNode, string> > variables;
    for (int i = 0; i < particlesPerSet; i++) {
        string index = cc.intToString(i+1);
        variables.push_back(makeVariable("x"+index, "pos"+index+".x"));
        variables.push_back(makeVariable("y"+index, "pos"+index+".y"));
        variables.push_back(makeVariable("z"+index, "pos"+index+".z"));
    }
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        for (int j = 0; j < particlesPerSet; j++) {
            string index = cc.intToString(j+1);
            variables.push_back(makeVariable(name+index, "params"+params->getParameterSuffix(i, index)));
        }
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cc, force.getNumGlobalParameters(), "customManyParticleGlobals");
        globals.upload(globalParamValues);
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cc.intToString(i)+"]";
            variables.push_back(makeVariable(name, value));
        }
    }
    
    // Build data structures for type filters.
    
    vector<int> particleTypesVec;
    vector<int> orderIndexVec;
    vector<std::vector<int> > particleOrderVec;
    int numTypes;
    CustomManyParticleForceImpl::buildFilterArrays(force, numTypes, particleTypesVec, orderIndexVec, particleOrderVec);
    bool hasTypeFilters = (particleOrderVec.size() > 1);
    if (hasTypeFilters) {
        particleTypes.initialize<int>(cc, particleTypesVec.size(), "customManyParticleTypes");
        orderIndex.initialize<int>(cc, orderIndexVec.size(), "customManyParticleOrderIndex");
        particleOrder.initialize<int>(cc, particleOrderVec.size()*particlesPerSet, "customManyParticleOrder");
        particleTypes.upload(particleTypesVec);
        orderIndex.upload(orderIndexVec);
        vector<int> flattenedOrder(particleOrder.getSize());
        for (int i = 0; i < (int) particleOrderVec.size(); i++)
            for (int j = 0; j < particlesPerSet; j++)
                flattenedOrder[i*particlesPerSet+j] = particleOrderVec[i][j];
        particleOrder.upload(flattenedOrder);
    }
    
    // Build data structures for exclusions.
    
    if (force.getNumExclusions() > 0) {
        vector<vector<int> > particleExclusions(numParticles);
        for (int i = 0; i < force.getNumExclusions(); i++) {
            int p1, p2;
            force.getExclusionParticles(i, p1, p2);
            particleExclusions[p1].push_back(p2);
            particleExclusions[p2].push_back(p1);
        }
        vector<int> exclusionsVec;
        vector<int> exclusionStartIndexVec(numParticles+1);
        exclusionStartIndexVec[0] = 0;
        for (int i = 0; i < numParticles; i++) {
            sort(particleExclusions[i].begin(), particleExclusions[i].end());
            exclusionsVec.insert(exclusionsVec.end(), particleExclusions[i].begin(), particleExclusions[i].end());
            exclusionStartIndexVec[i+1] = exclusionsVec.size();
        }
        exclusions.initialize<int>(cc, exclusionsVec.size(), "customManyParticleExclusions");
        exclusionStartIndex.initialize<int>(cc, exclusionStartIndexVec.size(), "customManyParticleExclusionStart");
        exclusions.upload(exclusionsVec);
        exclusionStartIndex.upload(exclusionStartIndexVec);
    }
    
    // Build data structures for the neighbor list.
    
    int numAtomBlocks = cc.getPaddedNumAtoms()/32;
    if (nonbondedMethod != NoCutoff) {
        int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter.initialize(cc, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(cc, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        numNeighborPairs.initialize<int>(cc, 1, "customManyParticleNumNeighborPairs");
        neighborStartIndex.initialize<int>(cc, numParticles+1, "customManyParticleNeighborStartIndex");
        numNeighborsForAtom.initialize<int>(cc, numParticles, "customManyParticleNumNeighborsForAtom");

        // Select a size for the array that holds the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxNeighborPairs = 150*numParticles;
        neighborPairs.initialize<mm_int2>(cc, maxNeighborPairs, "customManyParticleNeighborPairs");
        neighbors.initialize<int>(cc, maxNeighborPairs, "customManyParticleNeighbors");
    }

    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
    // and dihedrals the expression depends on.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomManyParticleForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    map<string, Lepton::ParsedExpression> forceExpressions;
    set<string> computedDeltas;
    vector<string> atomNames, posNames;
    for (int i = 0; i < particlesPerSet; i++) {
        string index = cc.intToString(i+1);
        atomNames.push_back("P"+index);
        posNames.push_back("pos"+index);
    }
    stringstream compute;
    int index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
        variables.push_back(makeVariable(distance.first, "r_"+deltaName));
        forceExpressions["real dEdDistance"+cc.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        variables.push_back(makeVariable(angle.first, angleName));
        forceExpressions["real dEdAngle"+cc.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& atoms = dihedral.second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        string dihedralName = "dihedral_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]]+atomNames[atoms[3]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            compute<<"real4 delta"<<deltaName3<<" = delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[3]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName3);
        }
        compute<<"real4 "<<crossName1<<" = computeCross(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        compute<<"real4 "<<crossName2<<" = computeCross(delta"<<deltaName2<<", delta"<<deltaName3<<");\n";
        compute<<"real "<<dihedralName<<" = computeAngle("<<crossName1<<", "<<crossName2<<");\n";
        compute<<dihedralName<<" *= (delta"<<deltaName1<<".x*"<<crossName2<<".x + delta"<<deltaName1<<".y*"<<crossName2<<".y + delta"<<deltaName1<<".z*"<<crossName2<<".z < 0 ? -1 : 1);\n";
        variables.push_back(makeVariable(dihedral.first, dihedralName));
        forceExpressions["real dEdDihedral"+cc.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        compute<<parameter.getType()<<" params"<<(i+1)<<" = global_params"<<(i+1)<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Apply forces to atoms.

    vector<string> forceNames;
    for (int i = 0; i < particlesPerSet; i++) {
        string istr = cc.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = make_real3(0);\n";
        compute<<"{\n";
        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
        map<string, Lepton::ParsedExpression> expressions;
        if (!isZeroExpression(forceExpressionX))
            expressions[forceName+".x -= "] = forceExpressionX;
        if (!isZeroExpression(forceExpressionY))
            expressions[forceName+".y -= "] = forceExpressionY;
        if (!isZeroExpression(forceExpressionZ))
            expressions[forceName+".z -= "] = forceExpressionZ;
        if (expressions.size() > 0)
            compute<<cc.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cc.intToString(index)+"/r_"+deltaName+")*trimTo3(delta"+deltaName+")";
        compute<<forceNames[atoms[0]]<<" += "<<"-"<<value<<";\n";
        compute<<forceNames[atoms[1]]<<" += "<<value<<";\n";
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute<<"{\n";
        compute<<"real3 crossProd = trimTo3(cross(delta"<<deltaName2<<", delta"<<deltaName1<<"));\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), (real) 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(trimTo3(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cc.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(trimTo3(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cc.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real3 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[atoms[0]]<<" += deltaCross0;\n";
        compute<<forceNames[atoms[1]]<<" += deltaCross1;\n";
        compute<<forceNames[atoms[2]]<<" += deltaCross2;\n";
        compute<<"}\n";
        index++;
    }
    index = 0;
    for (auto& dihedral : dihedrals) {
        const vector<int>& atoms = dihedral.second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        compute<<"{\n";
        compute<<"real r = sqrt(delta"<<deltaName2<<".w);\n";
        compute<<"real4 ff;\n";
        compute<<"ff.x = (-dEdDihedral"<<cc.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cc.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real3 internalF0 = ff.x*trimTo3("<<crossName1<<");\n";
        compute<<"real3 internalF3 = ff.w*trimTo3("<<crossName2<<");\n";
        compute<<"real3 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[atoms[0]]<<" += internalF0;\n";
        compute<<forceNames[atoms[1]]<<" += s-internalF0;\n";
        compute<<forceNames[atoms[2]]<<" += -s-internalF3;\n";
        compute<<forceNames[atoms[3]]<<" += internalF3;\n";
        compute<<"}\n";
        index++;
    }
    
    // Store forces to global memory.
    
    for (int i = 0; i < particlesPerSet; i++)
        compute<<"storeForce(atom"<<(i+1)<<", "<<forceNames[i]<<", forceBuffers);\n";
    
    // Create other replacements that depend on the number of particles per set.
    
    stringstream numCombinations, atomsForCombination, isValidCombination, permute, loadData, verifyCutoff, verifyExclusions;
    if (hasTypeFilters) {
        permute<<"int particleSet[] = {";
        for (int i = 0; i < particlesPerSet; i++) {
            permute<<"p"<<(i+1);
            if (i < particlesPerSet-1)
                permute<<", ";
        }
        permute<<"};\n";
    }
    for (int i = 0; i < particlesPerSet; i++) {
        if (hasTypeFilters)
            permute<<"int atom"<<(i+1)<<" = particleSet[particleOrder["<<particlesPerSet<<"*order+"<<i<<"]];\n";
        else
            permute<<"int atom"<<(i+1)<<" = p"<<(i+1)<<";\n";
        loadData<<"real3 pos"<<(i+1)<<" = trimTo3(posq[atom"<<(i+1)<<"]);\n";
        for (int j = 0; j < (int) params->getParameterInfos().size(); j++)
            loadData<<params->getParameterInfos()[j].getType()<<" params"<<(j+1)<<(i+1)<<" = global_params"<<(j+1)<<"[atom"<<(i+1)<<"];\n";
    }
    if (centralParticleMode) {
        for (int i = 1; i < particlesPerSet; i++) {
            if (i > 1)
                isValidCombination<<" && p"<<(i+1)<<">p"<<i<<" && ";
            isValidCombination<<"p"<<(i+1)<<"!=p1";
        }
    }
    else {
        for (int i = 2; i < particlesPerSet; i++) {
            if (i > 2)
                isValidCombination<<" && ";
            isValidCombination<<"a"<<(i+1)<<">a"<<i;
        }
    }
    atomsForCombination<<"int tempIndex = index;\n";
    for (int i = 1; i < particlesPerSet; i++) {
        if (i > 1)
            numCombinations<<"*";
        numCombinations<<"numNeighbors";
        if (centralParticleMode)
            atomsForCombination<<"int a"<<(i+1)<<" = tempIndex%numNeighbors;\n";
        else
            atomsForCombination<<"int a"<<(i+1)<<" = 1+tempIndex%numNeighbors;\n";
        if (i < particlesPerSet-1)
            atomsForCombination<<"tempIndex /= numNeighbors;\n";
    }
    if (particlesPerSet > 2) {
        if (centralParticleMode)
            atomsForCombination<<"a2 = (a3%2 == 0 ? a2 : numNeighbors-a2-1);\n";
        else
            atomsForCombination<<"a2 = (a3%2 == 0 ? a2 : numNeighbors-a2+1);\n";
    }
    for (int i = 1; i < particlesPerSet; i++) {
        if (nonbondedMethod == NoCutoff) {
            if (centralParticleMode)
                atomsForCombination<<"int p"<<(i+1)<<" = a"<<(i+1)<<";\n";
            else
                atomsForCombination<<"int p"<<(i+1)<<" = p1+a"<<(i+1)<<";\n";
        }
        else {
            if (centralParticleMode)
                atomsForCombination<<"int p"<<(i+1)<<" = neighbors[firstNeighbor+a"<<(i+1)<<"];\n";
            else
                atomsForCombination<<"int p"<<(i+1)<<" = neighbors[firstNeighbor-1+a"<<(i+1)<<"];\n";
        }
    }
    if (nonbondedMethod != NoCutoff) {
        for (int i = 1; i < particlesPerSet; i++)
            verifyCutoff<<"real3 pos"<<(i+1)<<" = trimTo3(posq[p"<<(i+1)<<"]);\n";
        if (!centralParticleMode) {
            for (int i = 1; i < particlesPerSet; i++) {
                for (int j = i+1; j < particlesPerSet; j++)
                    verifyCutoff<<"includeInteraction &= (delta(pos"<<(i+1)<<", pos"<<(j+1)<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ).w < CUTOFF_SQUARED);\n";
            }
        }
    }
    if (force.getNumExclusions() > 0) {
        int startCheckFrom = (nonbondedMethod == NoCutoff ? 0 : 1);
        for (int i = startCheckFrom; i < particlesPerSet; i++)
            for (int j = i+1; j < particlesPerSet; j++)
                verifyExclusions<<"includeInteraction &= !isInteractionExcluded(p"<<(i+1)<<", p"<<(j+1)<<", exclusions, exclusionStartIndex);\n";
    }
    string computeTypeIndex = "particleTypes[p"+cc.intToString(particlesPerSet)+"]";
    for (int i = particlesPerSet-2; i >= 0; i--)
        computeTypeIndex = "particleTypes[p"+cc.intToString(i+1)+"]+"+cc.intToString(numTypes)+"*("+computeTypeIndex+")";
    
    // Create replacements for extra arguments.
    
    stringstream extraArgs;
    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", GLOBAL const float* globals";
    for (int i = 0; i < (int) params->getParameterInfos().size(); i++) {
        ComputeParameterInfo& parameter = params->getParameterInfos()[i];
        extraArgs<<", GLOBAL const "<<parameter.getType()<<"* RESTRICT global_params"<<(i+1);
    }

    // Create the kernels.

    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = compute.str();
    replacements["NUM_CANDIDATE_COMBINATIONS"] = numCombinations.str();
    replacements["FIND_ATOMS_FOR_COMBINATION_INDEX"] = atomsForCombination.str();
    replacements["IS_VALID_COMBINATION"] = isValidCombination.str();
    replacements["VERIFY_CUTOFF"] = verifyCutoff.str();
    replacements["VERIFY_EXCLUSIONS"] = verifyExclusions.str();
    replacements["PERMUTE_ATOMS"] = permute.str();
    replacements["LOAD_PARTICLE_DATA"] = loadData.str();
    replacements["COMPUTE_TYPE_INDEX"] = computeTypeIndex;
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string, string> defines;
    if (nonbondedMethod != NoCutoff)
        defines["USE_CUTOFF"] = "1";
    if (nonbondedMethod == CutoffPeriodic)
        defines["USE_PERIODIC"] = "1";
    if (centralParticleMode)
        defines["USE_CENTRAL_PARTICLE"] = "1";
    if (hasTypeFilters)
        defines["USE_FILTERS"] = "1";
    if (force.getNumExclusions() > 0)
        defines["USE_EXCLUSIONS"] = "1";
    defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    defines["M_PI"] = cc.doubleToString(M_PI);
    defines["CUTOFF_SQUARED"] = cc.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    defines["TILE_SIZE"] = cc.intToString(32);
    defines["NUM_BLOCKS"] = cc.intToString(numAtomBlocks);
    defines["FIND_NEIGHBORS_WORKGROUP_SIZE"] = cc.intToString(findNeighborsWorkgroupSize);
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customManyParticle, replacements), defines);
    forceKernel = program->createKernel("computeInteraction");
    blockBoundsKernel = program->createKernel("findBlockBounds");
    neighborsKernel = program->createKernel("findNeighbors");
    startIndicesKernel = program->createKernel("computeNeighborStartIndices");
    copyPairsKernel = program->createKernel("copyPairsToNeighborList");
    event = cc.createEvent();
}

double CommonCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        
        // Set arguments for the force kernel.
        
        forceKernel->addArg(cc.getLongForceBuffer());
        forceKernel->addArg(cc.getEnergyBuffer());
        forceKernel->addArg(cc.getPosq());
        for (int i = 0; i < 5; i++)
            forceKernel->addArg();
        setPeriodicBoxArgs(cc, forceKernel, 3);
        if (nonbondedMethod != NoCutoff) {
            forceKernel->addArg(neighbors);
            forceKernel->addArg(neighborStartIndex);
        }
        if (particleTypes.isInitialized()) {
            forceKernel->addArg(particleTypes);
            forceKernel->addArg(orderIndex);
            forceKernel->addArg(particleOrder);
        }
        if (exclusions.isInitialized()) {
            forceKernel->addArg(exclusions);
            forceKernel->addArg(exclusionStartIndex);
        }
        if (globals.isInitialized())
            forceKernel->addArg(globals);
        for (auto& parameter : params->getParameterInfos())
            forceKernel->addArg(parameter.getArray());
        for (auto& function : tabulatedFunctions)
            forceKernel->addArg(function);
        
        if (nonbondedMethod != NoCutoff) {
            // Set arguments for the block bounds kernel.

            for (int i = 0; i < 5; i++)
                blockBoundsKernel->addArg(); // Periodic box information will be set just before it is executed.
            blockBoundsKernel->addArg(cc.getPosq());
            blockBoundsKernel->addArg(blockCenter);
            blockBoundsKernel->addArg(blockBoundingBox);
            blockBoundsKernel->addArg(numNeighborPairs);

            // Set arguments for the neighbor list kernel.

            for (int i = 0; i < 5; i++)
                neighborsKernel->addArg(); // Periodic box information will be set just before it is executed.
            neighborsKernel->addArg(cc.getPosq());
            neighborsKernel->addArg(blockCenter);
            neighborsKernel->addArg(blockBoundingBox);
            neighborsKernel->addArg(neighborPairs);
            neighborsKernel->addArg(numNeighborPairs);
            neighborsKernel->addArg(numNeighborsForAtom);
            neighborsKernel->addArg(maxNeighborPairs);
            if (exclusions.isInitialized()) {
                neighborsKernel->addArg(exclusions);
                neighborsKernel->addArg(exclusionStartIndex);
            }
            
            // Set arguments for the kernel to find neighbor list start indices.
            
            startIndicesKernel->addArg(numNeighborsForAtom);
            startIndicesKernel->addArg(neighborStartIndex);
            startIndicesKernel->addArg(numNeighborPairs);
            startIndicesKernel->addArg(maxNeighborPairs);

            // Set arguments for the kernel to assemble the final neighbor list.
            
            copyPairsKernel->addArg(neighborPairs);
            copyPairsKernel->addArg(neighbors);
            copyPairsKernel->addArg(numNeighborPairs);
            copyPairsKernel->addArg(maxNeighborPairs);
            copyPairsKernel->addArg(numNeighborsForAtom);
            copyPairsKernel->addArg(neighborStartIndex);
       }
    }
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    while (true) {
        int* numPairs = (int*) cc.getPinnedBuffer();
        if (nonbondedMethod != NoCutoff) {
            setPeriodicBoxArgs(cc, forceKernel, 3);
            setPeriodicBoxArgs(cc, blockBoundsKernel, 0);
            setPeriodicBoxArgs(cc, neighborsKernel, 0);
            blockBoundsKernel->execute(cc.getPaddedNumAtoms()/32);
            neighborsKernel->execute(cc.getNumAtoms(), findNeighborsWorkgroupSize);

            // We need to make sure there was enough memory for the neighbor list.  Download the
            // information asynchronously so kernels can be running at the same time.

            numNeighborPairs.download(numPairs, false);
            event->enqueue();
            startIndicesKernel->execute(256, 256);
            copyPairsKernel->execute(maxNeighborPairs);
        }
        int maxThreads = min(cc.getNumAtoms()*forceWorkgroupSize, cc.getEnergyBuffer().getSize());
        forceKernel->execute(maxThreads, forceWorkgroupSize);
        if (nonbondedMethod != NoCutoff) {
            // Make sure there was enough memory for the neighbor list.

            event->wait();
            if (*numPairs > maxNeighborPairs) {
                // Resize the arrays and run the calculation again.

                maxNeighborPairs = (int) (1.1*(*numPairs));
                neighborPairs.resize(maxNeighborPairs);
                neighbors.resize(maxNeighborPairs);
                neighborsKernel->setArg(11, maxNeighborPairs);
                startIndicesKernel->setArg(3, maxNeighborPairs);
                copyPairsKernel->setArg(3, maxNeighborPairs);
                continue;
            }
        }
        break;
    }
    return 0.0;
}

void CommonCalcCustomManyParticleForceKernel::copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) {
    cc.setAsCurrent();
    int numParticles = force.getNumParticles();
    if (numParticles != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<vector<float> > paramVector(numParticles);
    vector<double> parameters;
    int type;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters, type);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info);
}

class CommonCalcGayBerneForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const GayBerneForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        int xparticle1, yparticle1;
        double sigma1, epsilon1, sx1, sy1, sz1, ex1, ey1, ez1;
        int xparticle2, yparticle2;
        double sigma2, epsilon2, sx2, sy2, sz2, ex2, ey2, ez2;
        force.getParticleParameters(particle1, sigma1, epsilon1, xparticle1, yparticle1, sx1, sy1, sz1, ex1, ey1, ez1);
        force.getParticleParameters(particle2, sigma2, epsilon2, xparticle2, yparticle2, sx2, sy2, sz2, ex2, ey2, ez2);
        return (sigma1 == sigma2 && epsilon1 == epsilon2 && sx1 == sx2 && sy1 == sy2 && sz1 == sz2 && ex1 == ex2 && ey1 == ey2 && ez1 == ez2);
    }
    int getNumParticleGroups() {
        return force.getNumExceptions()+force.getNumParticles();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        if (index < force.getNumExceptions()) {
            int particle1, particle2;
            double sigma, epsilon;
            force.getExceptionParameters(index, particle1, particle2, sigma, epsilon);
            particles.resize(2);
            particles[0] = particle1;
            particles[1] = particle2;
        }
        else {
            int particle = index-force.getNumExceptions();
            int xparticle, yparticle;
            double sigma, epsilon, sx, sy, sz, ex, ey, ez;
            force.getParticleParameters(particle, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
            particles.clear();
            particles.push_back(particle);
            if (xparticle > -1)
                particles.push_back(xparticle);
            if (yparticle > -1)
                particles.push_back(yparticle);
        }
    }
    bool areGroupsIdentical(int group1, int group2) {
        if (group1 < force.getNumExceptions() && group2 < force.getNumExceptions()) {
            int particle1, particle2;
            double sigma1, sigma2, epsilon1, epsilon2;
            force.getExceptionParameters(group1, particle1, particle2, sigma1, epsilon1);
            force.getExceptionParameters(group2, particle1, particle2, sigma2, epsilon2);
            return (sigma1 == sigma2 && epsilon1 == epsilon2);
        }
        return true;
    }
private:
    const GayBerneForce& force;
};

class CommonCalcGayBerneForceKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(CommonCalcGayBerneForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.sortAtoms();
    }
private:
    CommonCalcGayBerneForceKernel& owner;
};

void CommonCalcGayBerneForceKernel::initialize(const System& system, const GayBerneForce& force) {
    if (!cc.getSupports64BitGlobalAtomics())
        throw OpenMMException("GayBerneForce requires a device that supports 64 bit atomic operations");

    // Initialize interactions.

    int numParticles = force.getNumParticles();
    sigParams.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "sigParams");
    epsParams.initialize<mm_float2>(cc, cc.getPaddedNumAtoms(), "epsParams");
    scale.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "scale");
    axisParticleIndices.initialize<mm_int2>(cc, cc.getPaddedNumAtoms(), "axisParticleIndices");
    sortedParticles.initialize<int>(cc, cc.getPaddedNumAtoms(), "sortedParticles");
    aMatrix.initialize<float>(cc, 9*cc.getPaddedNumAtoms(), "aMatrix");
    bMatrix.initialize<float>(cc, 9*cc.getPaddedNumAtoms(), "bMatrix");
    gMatrix.initialize<float>(cc, 9*cc.getPaddedNumAtoms(), "gMatrix");
    vector<mm_float4> sigParamsVector(cc.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_float2> epsParamsVector(cc.getPaddedNumAtoms(), mm_float2(0, 0));
    vector<mm_float4> scaleVector(cc.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_int2> axisParticleVector(cc.getPaddedNumAtoms(), mm_int2(0, 0));
    isRealParticle.resize(cc.getPaddedNumAtoms());
    for (int i = 0; i < numParticles; i++) {
        int xparticle, yparticle;
        double sigma, epsilon, sx, sy, sz, ex, ey, ez;
        force.getParticleParameters(i, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
        axisParticleVector[i] = mm_int2(xparticle, yparticle);
        sigParamsVector[i] = mm_float4((float) (0.5*sigma), (float) (0.25*sx*sx), (float) (0.25*sy*sy), (float) (0.25*sz*sz));
        epsParamsVector[i] = mm_float2((float) sqrt(epsilon), (float) (0.125*(sx*sy + sz*sz)*sqrt(sx*sy)));
        scaleVector[i] = mm_float4((float) (1/sqrt(ex)), (float) (1/sqrt(ey)), (float) (1/sqrt(ez)), 0);
        isRealParticle[i] = (epsilon != 0.0);
    }
    sigParams.upload(sigParamsVector);
    epsParams.upload(epsParamsVector);
    scale.upload(scaleVector);
    axisParticleIndices.upload(axisParticleVector);
    
    // Record exceptions and exclusions.

    vector<mm_float2> exceptionParamsVec;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, sigma, epsilon);
        if (epsilon != 0.0) {
            exceptionParamsVec.push_back(mm_float2((float) sigma, (float) epsilon));
            exceptionAtoms.push_back(make_pair(particle1, particle2));
            isRealParticle[particle1] = true;
            isRealParticle[particle2] = true;
        }
        if (isRealParticle[particle1] && isRealParticle[particle2])
            excludedPairs.push_back(pair<int, int>(particle1, particle2));
    }
    numRealParticles = 0;
    for (int i = 0; i < isRealParticle.size(); i++)
        if (isRealParticle[i])
            numRealParticles++;
    int numExceptions = exceptionParamsVec.size();
    exclusions.initialize<int>(cc, max(1, (int) excludedPairs.size()), "exclusions");
    exclusionStartIndex.initialize<int>(cc, numRealParticles+1, "exclusionStartIndex");
    exceptionParticles.initialize<mm_int4>(cc, max(1, numExceptions), "exceptionParticles");
    exceptionParams.initialize<mm_float2>(cc, max(1, numExceptions), "exceptionParams");
    if (numExceptions > 0)
        exceptionParams.upload(exceptionParamsVec);
    
    // Create data structures used for the neighbor list.

    int numAtomBlocks = (numRealParticles+31)/32;
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    blockCenter.initialize(cc, numAtomBlocks, 4*elementSize, "blockCenter");
    blockBoundingBox.initialize(cc, numAtomBlocks, 4*elementSize, "blockBoundingBox");
    sortedPos.initialize(cc, numRealParticles, 4*elementSize, "sortedPos");
    maxNeighborBlocks = numRealParticles*2;
    neighbors.initialize<int>(cc, maxNeighborBlocks*32, "neighbors");
    neighborIndex.initialize<int>(cc, maxNeighborBlocks, "neighborIndex");
    neighborBlockCount.initialize<int>(cc, 1, "neighborBlockCount");
    event = cc.createEvent();

    // Create array for accumulating torques.
    
    torque.initialize<long long>(cc, 3*cc.getPaddedNumAtoms(), "torque");
    cc.addAutoclearBuffer(torque);

    // Create the kernels.
    
    nonbondedMethod = force.getNonbondedMethod();
    bool useCutoff = (nonbondedMethod != GayBerneForce::NoCutoff);
    bool usePeriodic = (nonbondedMethod == GayBerneForce::CutoffPeriodic);
    map<string, string> defines;
    defines["USE_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff*cutoff);
    if (useCutoff) {
        defines["USE_CUTOFF"] = 1;
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        
        // Compute the switching coefficients.
        
        if (force.getUseSwitchingFunction()) {
            defines["SWITCH_CUTOFF"] = cc.doubleToString(force.getSwitchingDistance());
            defines["SWITCH_C3"] = cc.doubleToString(10/pow(force.getSwitchingDistance()-cutoff, 3.0));
            defines["SWITCH_C4"] = cc.doubleToString(15/pow(force.getSwitchingDistance()-cutoff, 4.0));
            defines["SWITCH_C5"] = cc.doubleToString(6/pow(force.getSwitchingDistance()-cutoff, 5.0));
        }
    }
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    ComputeProgram program = cc.compileProgram(CommonKernelSources::gayBerne, defines);
    framesKernel = program->createKernel("computeEllipsoidFrames");
    blockBoundsKernel = program->createKernel("findBlockBounds");
    neighborsKernel = program->createKernel("findNeighbors");
    forceKernel = program->createKernel("computeForce");
    torqueKernel = program->createKernel("applyTorques");
    info = new ForceInfo(force);
    cc.addForce(info);
    cc.addReorderListener(new ReorderListener(*this));
}

double CommonCalcGayBerneForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        sortAtoms();
        framesKernel->addArg(numRealParticles);
        framesKernel->addArg(cc.getPosq());
        framesKernel->addArg(axisParticleIndices);
        framesKernel->addArg(sigParams);
        framesKernel->addArg(scale);
        framesKernel->addArg(aMatrix);
        framesKernel->addArg(bMatrix);
        framesKernel->addArg(gMatrix);
        framesKernel->addArg(sortedParticles);
        blockBoundsKernel->addArg(numRealParticles);
        for (int i = 0; i < 5; i++)
            blockBoundsKernel->addArg(); // Periodic box information will be set just before it is executed.
        blockBoundsKernel->addArg(sortedParticles);
        blockBoundsKernel->addArg(cc.getPosq());
        blockBoundsKernel->addArg(sortedPos);
        blockBoundsKernel->addArg(blockCenter);
        blockBoundsKernel->addArg(blockBoundingBox);
        blockBoundsKernel->addArg(neighborBlockCount);
        neighborsKernel->addArg(numRealParticles);
        neighborsKernel->addArg(maxNeighborBlocks);
        for (int i = 0; i < 5; i++)
            neighborsKernel->addArg(); // Periodic box information will be set just before it is executed.
        neighborsKernel->addArg(sortedPos);
        neighborsKernel->addArg(blockCenter);
        neighborsKernel->addArg(blockBoundingBox);
        neighborsKernel->addArg(neighbors);
        neighborsKernel->addArg(neighborIndex);
        neighborsKernel->addArg(neighborBlockCount);
        neighborsKernel->addArg(exclusions);
        neighborsKernel->addArg(exclusionStartIndex);
        forceKernel->addArg(cc.getLongForceBuffer());
        forceKernel->addArg(torque);
        forceKernel->addArg(numRealParticles);
        forceKernel->addArg((int) exceptionAtoms.size());
        forceKernel->addArg(cc.getEnergyBuffer());
        forceKernel->addArg(sortedPos);
        forceKernel->addArg(sigParams);
        forceKernel->addArg(epsParams);
        forceKernel->addArg(sortedParticles);
        forceKernel->addArg(aMatrix);
        forceKernel->addArg(bMatrix);
        forceKernel->addArg(gMatrix);
        forceKernel->addArg(exclusions);
        forceKernel->addArg(exclusionStartIndex);
        forceKernel->addArg(exceptionParticles);
        forceKernel->addArg(exceptionParams);
        if (nonbondedMethod != GayBerneForce::NoCutoff) {
            forceKernel->addArg(maxNeighborBlocks);
            forceKernel->addArg(neighbors);
            forceKernel->addArg(neighborIndex);
            forceKernel->addArg(neighborBlockCount);
            for (int i = 0; i < 5; i++)
                forceKernel->addArg(); // Periodic box information will be set just before it is executed.
        }
        torqueKernel->addArg(cc.getLongForceBuffer());
        torqueKernel->addArg(torque);
        torqueKernel->addArg(numRealParticles);
        torqueKernel->addArg(cc.getPosq());
        torqueKernel->addArg(axisParticleIndices);
        torqueKernel->addArg(sortedParticles);
    }
    framesKernel->execute(numRealParticles);
    setPeriodicBoxArgs(cc, blockBoundsKernel, 1);
    blockBoundsKernel->execute((numRealParticles+31)/32);
    if (nonbondedMethod == GayBerneForce::NoCutoff)
        forceKernel->execute(cc.getNonbondedUtilities().getNumForceThreadBlocks()*cc.getNonbondedUtilities().getForceThreadBlockSize());
    else {
        while (true) {
            setPeriodicBoxArgs(cc, neighborsKernel, 2);
            neighborsKernel->execute(numRealParticles);
            int* count = (int*) cc.getPinnedBuffer();
            neighborBlockCount.download(count, false);
            event->enqueue();
            setPeriodicBoxArgs(cc, forceKernel, 20);
            forceKernel->execute(cc.getNonbondedUtilities().getNumForceThreadBlocks()*cc.getNonbondedUtilities().getForceThreadBlockSize());
            event->wait();
            if (*count <= maxNeighborBlocks)
                break;
            
            // There wasn't enough room for the neighbor list, so we need to recreate it.

            maxNeighborBlocks = (int) ceil((*count)*1.1);
            neighbors.resize(maxNeighborBlocks*32);
            neighborIndex.resize(maxNeighborBlocks);
            neighborsKernel->setArg(10, neighbors);
            neighborsKernel->setArg(11, neighborIndex);
            forceKernel->setArg(17, neighbors);
            forceKernel->setArg(18, neighborIndex);
        }
    }
    torqueKernel->execute(numRealParticles);
    return 0.0;
}

void CommonCalcGayBerneForceKernel::copyParametersToContext(ContextImpl& context, const GayBerneForce& force) {
    // Make sure the new parameters are acceptable.
    
    if (force.getNumParticles() != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    vector<int> exceptions;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, sigma, epsilon);
        if (exceptionAtoms.size() > exceptions.size() && make_pair(particle1, particle2) == exceptionAtoms[exceptions.size()])
            exceptions.push_back(i);
        else if (epsilon != 0.0)
            throw OpenMMException("updateParametersInContext: The set of non-excluded exceptions has changed");
    }
    int numExceptions = exceptionAtoms.size();
    
    // Record the per-particle parameters.
    
    vector<mm_float4> sigParamsVector(cc.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_float2> epsParamsVector(cc.getPaddedNumAtoms(), mm_float2(0, 0));
    vector<mm_float4> scaleVector(cc.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    for (int i = 0; i < force.getNumParticles(); i++) {
        int xparticle, yparticle;
        double sigma, epsilon, sx, sy, sz, ex, ey, ez;
        force.getParticleParameters(i, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
        sigParamsVector[i] = mm_float4((float) (0.5*sigma), (float) (0.25*sx*sx), (float) (0.25*sy*sy), (float) (0.25*sz*sz));
        epsParamsVector[i] = mm_float2((float) sqrt(epsilon), (float) (0.125*(sx*sy + sz*sz)*sqrt(sx*sy)));
        scaleVector[i] = mm_float4((float) (1/sqrt(ex)), (float) (1/sqrt(ey)), (float) (1/sqrt(ez)), 0);
        if (epsilon != 0.0 && !isRealParticle[i])
            throw OpenMMException("updateParametersInContext: The set of ignored particles (ones with epsilon=0) has changed");
    }
    sigParams.upload(sigParamsVector);
    epsParams.upload(epsParamsVector);
    scale.upload(scaleVector);
    
    // Record the exceptions.
    
    if (numExceptions > 0) {
        vector<mm_float2> exceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            int atom1, atom2;
            double sigma, epsilon;
            force.getExceptionParameters(exceptions[i], atom1, atom2, sigma, epsilon);
            exceptionParamsVec[i] = mm_float2((float) sigma, (float) epsilon);
        }
        exceptionParams.upload(exceptionParamsVec);
    }
    cc.invalidateMolecules(info);
    sortAtoms();
}

void CommonCalcGayBerneForceKernel::sortAtoms() {
    // Sort the list of atoms by type to avoid thread divergence.  This is executed every time
    // the atoms are reordered.
    
    int nextIndex = 0;
    vector<int> particles(cc.getPaddedNumAtoms(), 0);
    const vector<int>& order = cc.getAtomIndex();
    vector<int> inverseOrder(order.size(), -1);
    for (int i = 0; i < cc.getNumAtoms(); i++) {
        int atom = order[i];
        if (isRealParticle[atom]) {
            inverseOrder[atom] = nextIndex;
            particles[nextIndex++] = atom;
        }
    }
    sortedParticles.upload(particles);
    
    // Update the list of exception particles.
    
    int numExceptions = exceptionAtoms.size();
    if (numExceptions > 0) {
        vector<mm_int4> exceptionParticlesVec(numExceptions);
        for (int i = 0; i < numExceptions; i++)
            exceptionParticlesVec[i] = mm_int4(exceptionAtoms[i].first, exceptionAtoms[i].second, inverseOrder[exceptionAtoms[i].first], inverseOrder[exceptionAtoms[i].second]);
        exceptionParticles.upload(exceptionParticlesVec);
    }
    
    // Rebuild the list of exclusions.
    
    vector<vector<int> > excludedAtoms(numRealParticles);
    for (int i = 0; i < excludedPairs.size(); i++) {
        int first = inverseOrder[min(excludedPairs[i].first, excludedPairs[i].second)];
        int second = inverseOrder[max(excludedPairs[i].first, excludedPairs[i].second)];
        excludedAtoms[first].push_back(second);
    }
    int index = 0;
    vector<int> exclusionVec(exclusions.getSize());
    vector<int> startIndexVec(exclusionStartIndex.getSize());
    for (int i = 0; i < numRealParticles; i++) {
        startIndexVec[i] = index;
        for (int j = 0; j < excludedAtoms[i].size(); j++)
            exclusionVec[index++] = excludedAtoms[i][j];
    }
    startIndexVec[numRealParticles] = index;
    exclusions.upload(exclusionVec);
    exclusionStartIndex.upload(startIndexVec);
}

void CommonIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cc.initializeContexts();
    cc.setAsCurrent();
    ComputeProgram program = cc.compileProgram(CommonKernelSources::verlet);
    kernel1 = program->createKernel("integrateVerletPart1");
    kernel2 = program->createKernel("integrateVerletPart2");
}

void CommonIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    cc.setAsCurrent();
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(integration.getStepSize());
        kernel1->addArg(cc.getPosq());
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel1->addArg(cc.getPosqCorrection());
        kernel2->addArg(numAtoms);
        kernel2->addArg(integration.getStepSize());
        kernel2->addArg(cc.getPosq());
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel2->addArg(cc.getPosqCorrection());
    }
    integration.setNextStepSize(dt);

    // Call the first integration kernel.

    kernel1->execute(numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2->execute(numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cc.setTime(cc.getTime()+dt);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cc.flushQueue();
#endif
}

double CommonIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void CommonIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    cc.initializeContexts();
    cc.setAsCurrent();
    ComputeProgram program = cc.compileProgram(CommonKernelSources::verlet);
    kernel1 = program->createKernel("integrateVerletPart1");
    kernel2 = program->createKernel("integrateVerletPart2");
    selectSizeKernel = program->createKernel("selectVerletStepSize");
    blockSize = min(256, system.getNumParticles());
}

double CommonIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    cc.setAsCurrent();
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(integration.getStepSize());
        kernel1->addArg(cc.getPosq());
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel1->addArg(cc.getPosqCorrection());
        kernel2->addArg(numAtoms);
        kernel2->addArg(integration.getStepSize());
        kernel2->addArg(cc.getPosq());
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel2->addArg(cc.getPosqCorrection());
        selectSizeKernel->addArg(numAtoms);
        selectSizeKernel->addArg(paddedNumAtoms);
        selectSizeKernel->addArg();
        selectSizeKernel->addArg();
        selectSizeKernel->addArg(cc.getIntegrationUtilities().getStepSize());
        selectSizeKernel->addArg(cc.getVelm());
        selectSizeKernel->addArg(cc.getLongForceBuffer());
    }

    // Select the step size to use.

    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    double maxStepSize = maxTime-cc.getTime();
    if (integrator.getMaximumStepSize() > 0)
        maxStepSize = min(integrator.getMaximumStepSize(), maxStepSize);
    float maxStepSizeFloat = (float) maxStepSize;
    if (useDouble) {
        selectSizeKernel->setArg(2, maxStepSize);
        selectSizeKernel->setArg(3, integrator.getErrorTolerance());
    }
    else {
        selectSizeKernel->setArg(2, maxStepSizeFloat);
        selectSizeKernel->setArg(3, (float) integrator.getErrorTolerance());
    }
    selectSizeKernel->execute(blockSize, blockSize);

    // Call the first integration kernel.

    kernel1->execute(numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2->execute(numAtoms);
    integration.computeVirtualSites();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cc.flushQueue();
#endif

    // Update the time and step count.

    double dt = cc.getIntegrationUtilities().getLastStepSize();
    double time = cc.getTime()+dt;
    if (useDouble) {
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
        if (dt == maxStepSizeFloat)
            time = maxTime; // Avoid round-off error
    }
    cc.setTime(time);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    return dt;
}

double CommonIntegrateVariableVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void CommonRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    cc.setAsCurrent();
    frequency = force.getFrequency();
    int numAtoms = cc.getNumAtoms();
    cmMomentum.initialize<mm_float3>(cc, cc.getPaddedNumAtoms(), "cmMomentum");
    double totalMass = 0.0;
    for (int i = 0; i < numAtoms; i++)
        totalMass += system.getParticleMass(i);
    map<string, string> defines;
    defines["INVERSE_TOTAL_MASS"] = cc.doubleToString(totalMass == 0 ? 0.0 : 1.0/totalMass);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::removeCM, defines);
    kernel1 = program->createKernel("calcCenterOfMassMomentum");
    kernel1->addArg(numAtoms);
    kernel1->addArg(cc.getVelm());
    kernel1->addArg(cmMomentum);
    kernel2 = program->createKernel("removeCenterOfMassMomentum");
    kernel2->addArg(numAtoms);
    kernel2->addArg(cc.getVelm());
    kernel2->addArg(cmMomentum);
}

void CommonRemoveCMMotionKernel::execute(ContextImpl& context) {
    cc.setAsCurrent();
    kernel1->execute(cc.getNumAtoms(), 64);
    kernel2->execute(cc.getNumAtoms(), 64);
}
