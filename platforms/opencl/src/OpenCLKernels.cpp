/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2010 Stanford University and the Authors.      *
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

#include "OpenCLKernels.h"
#include "OpenCLForceInfo.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomHbondForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "OpenCLExpressionUtilities.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLNonbondedUtilities.h"
#include "OpenCLKernelSources.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "openmm/internal/MSVC_erfc.h"
#include <cmath>
#include <set>

using namespace OpenMM;
using namespace std;

static string doubleToString(double value) {
    stringstream s;
    s.precision(8);
    s << scientific << value << "f";
    return s.str();
}

static string intToString(int value) {
    stringstream s;
    s << value;
    return s.str();
}

static bool isZeroExpression(const Lepton::ParsedExpression& expression) {
    const Lepton::Operation& op = expression.getRootNode().getOperation();
    if (op.getId() != Lepton::Operation::CONSTANT)
        return false;
    return (dynamic_cast<const Lepton::Operation::Constant&>(op).getValue() == 0.0);
}

void OpenCLCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void OpenCLCalcForcesAndEnergyKernel::beginForceComputation(ContextImpl& context) {
    if (cl.getNonbondedUtilities().getUseCutoff() && cl.getComputeForceCount()%100 == 0)
        cl.reorderAtoms();
    cl.setComputeForceCount(cl.getComputeForceCount()+1);
    cl.clearBuffer(cl.getForceBuffers());
    cl.getNonbondedUtilities().prepareInteractions();
}

void OpenCLCalcForcesAndEnergyKernel::finishForceComputation(ContextImpl& context) {
    cl.getNonbondedUtilities().computeInteractions();
    cl.reduceBuffer(cl.getForceBuffers(), cl.getNumForceBuffers());
}

void OpenCLCalcForcesAndEnergyKernel::beginEnergyComputation(ContextImpl& context) {
    if (cl.getNonbondedUtilities().getUseCutoff() && cl.getComputeForceCount()%100 == 0)
        cl.reorderAtoms();
    cl.setComputeForceCount(cl.getComputeForceCount()+1);
    cl.clearBuffer(cl.getEnergyBuffer());
    cl.getNonbondedUtilities().prepareInteractions();
}

double OpenCLCalcForcesAndEnergyKernel::finishEnergyComputation(ContextImpl& context) {
    cl.getNonbondedUtilities().computeInteractions();
    OpenCLArray<cl_float>& energy = cl.getEnergyBuffer();
    energy.download();
    double sum = 0.0f;
    for (int i = 0; i < energy.getSize(); i++)
        sum += energy[i];
    return sum;
}

void OpenCLUpdateStateDataKernel::initialize(const System& system) {
}

double OpenCLUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return cl.getTime();
}

void OpenCLUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    cl.setTime(time);
}

void OpenCLUpdateStateDataKernel::getPositions(ContextImpl& context, std::vector<Vec3>& positions) {
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    posq.download();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    mm_float4 periodicBoxSize = cl.getNonbondedUtilities().getPeriodicBoxSize();
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 pos = posq[i];
        mm_int4 offset = cl.getPosCellOffsets()[i];
        positions[order[i]] = Vec3(pos.x-offset.x*periodicBoxSize.x, pos.y-offset.y*periodicBoxSize.y, pos.z-offset.z*periodicBoxSize.z);
    }
}

void OpenCLUpdateStateDataKernel::setPositions(ContextImpl& context, const std::vector<Vec3>& positions) {
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    for (int i = 0; i < numParticles; ++i) {
        mm_float4& pos = posq[i];
        const Vec3& p = positions[order[i]];
        pos.x = (cl_float) p[0];
        pos.y = (cl_float) p[1];
        pos.z = (cl_float) p[2];
    }
    posq.upload();
    for (int i = 0; i < (int) cl.getPosCellOffsets().size(); i++)
        cl.getPosCellOffsets()[i] = mm_int4(0, 0, 0, 0);
}

void OpenCLUpdateStateDataKernel::getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) {
    OpenCLArray<mm_float4>& velm = cl.getVelm();
    velm.download();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 vel = velm[i];
        velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
    }
}

void OpenCLUpdateStateDataKernel::setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) {
    OpenCLArray<mm_float4>& velm = cl.getVelm();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    for (int i = 0; i < numParticles; ++i) {
        mm_float4& vel = velm[i];
        const Vec3& p = velocities[order[i]];
        vel.x = (cl_float) p[0];
        vel.y = (cl_float) p[1];
        vel.z = (cl_float) p[2];
    }
    velm.upload();
}

void OpenCLUpdateStateDataKernel::getForces(ContextImpl& context, std::vector<Vec3>& forces) {
    OpenCLArray<mm_float4>& force = cl.getForce();
    force.download();
    OpenCLArray<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    forces.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        mm_float4 f = force[i];
        forces[order[i]] = Vec3(f.x, f.y, f.z);
    }
}

class OpenCLBondForceInfo : public OpenCLForceInfo {
public:
    OpenCLBondForceInfo(int requiredBuffers, const HarmonicBondForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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

OpenCLCalcHarmonicBondForceKernel::~OpenCLCalcHarmonicBondForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    numBonds = force.getNumBonds();
    if (numBonds == 0)
        return;
    params = new OpenCLArray<mm_float2>(cl, numBonds, "bondParams");
    indices = new OpenCLArray<mm_int4>(cl, numBonds, "bondIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float2> paramVector(numBonds);
    vector<mm_int4> indicesVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
        indicesVector[i] = mm_int4(particle1, particle2, forceBufferCounter[particle1]++, forceBufferCounter[particle2]++);
    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLBondForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(OpenCLKernelSources::harmonicBondForce);
    kernel = cl::Kernel(program, "calcHarmonicBondForce");
}

void OpenCLCalcHarmonicBondForceKernel::executeForces(ContextImpl& context) {
    if (numBonds == 0)
        return;
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numBonds);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
        kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    }
    cl.executeKernel(kernel, numBonds);
}

double OpenCLCalcHarmonicBondForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLCustomBondForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomBondForceInfo(int requiredBuffers, const CustomBondForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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

OpenCLCalcCustomBondForceKernel::~OpenCLCalcCustomBondForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
    if (globals != NULL)
        delete globals;
}

void OpenCLCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    numBonds = force.getNumBonds();
    if (numBonds == 0)
        return;
    params = new OpenCLParameterSet(cl, force.getNumPerBondParameters(), numBonds, "customBondParams");
    indices = new OpenCLArray<mm_int4>(cl, numBonds, "customBondIndices");
    string extraArguments;
    if (force.getNumGlobalParameters() > 0) {
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customBondGlobals", false, CL_MEM_READ_ONLY);
        extraArguments += ", __constant float* globals";
    }
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<vector<cl_float> > paramVector(numBonds);
    vector<mm_int4> indicesVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int particle1, particle2;
        vector<double> parameters;
        force.getBondParameters(i, particle1, particle2, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
        indicesVector[i] = mm_int4(particle1, particle2, forceBufferCounter[particle1]++, forceBufferCounter[particle2]++);
    }
    params->setParameterValues(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLCustomBondForceInfo(maxBuffers, force));

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["float dEdR = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["r"] = "r";
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+intToString(i)+"]";
        variables[name] = value;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArguments += ", __global "+buffer.getType()+"* "+buffer.getName();
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<buffer.getName()<<"[index];\n";
    }
    vector<pair<string, string> > functions;
    compute << OpenCLExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["EXTRA_ARGUMENTS"] = extraArguments;
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customBondForce, replacements));
    kernel = cl::Kernel(program, "computeCustomBondForces");
}

void OpenCLCalcCustomBondForceKernel::executeForces(ContextImpl& context) {
    if (numBonds == 0)
        return;
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numBonds);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, indices->getDeviceBuffer());
        int nextIndex = 6;
        if (globals != NULL)
            kernel.setArg<cl::Buffer>(nextIndex++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            kernel.setArg<cl::Memory>(nextIndex++, buffer.getMemory());
        }
    }
    cl.executeKernel(kernel, numBonds);
}

double OpenCLCalcCustomBondForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLAngleForceInfo : public OpenCLForceInfo {
public:
    OpenCLAngleForceInfo(int requiredBuffers, const HarmonicAngleForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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

OpenCLCalcHarmonicAngleForceKernel::~OpenCLCalcHarmonicAngleForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    if (numAngles == 0)
        return;
    params = new OpenCLArray<mm_float2>(cl, numAngles, "angleParams");
    indices = new OpenCLArray<mm_int8>(cl, numAngles, "angleIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float2> paramVector(numAngles);
    vector<mm_int8> indicesVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        paramVector[i] = mm_float2((cl_float) angle, (cl_float) k);
        indicesVector[i] = mm_int8(particle1, particle2, particle3,
                forceBufferCounter[particle1]++, forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, 0, 0);

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLAngleForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(OpenCLKernelSources::harmonicAngleForce);
    kernel = cl::Kernel(program, "calcHarmonicAngleForce");
}

void OpenCLCalcHarmonicAngleForceKernel::executeForces(ContextImpl& context) {
    if (numAngles == 0)
        return;
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numAngles);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
        kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    }
    cl.executeKernel(kernel, numAngles);
}

double OpenCLCalcHarmonicAngleForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLCustomAngleForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomAngleForceInfo(int requiredBuffers, const CustomAngleForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumAngles();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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

OpenCLCalcCustomAngleForceKernel::~OpenCLCalcCustomAngleForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
    if (globals != NULL)
        delete globals;
}

void OpenCLCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    numAngles = force.getNumAngles();
    if (numAngles == 0)
        return;
    params = new OpenCLParameterSet(cl, force.getNumPerAngleParameters(), numAngles, "customAngleParams");
    indices = new OpenCLArray<mm_int8>(cl, numAngles, "customAngleIndices");
    string extraArguments;
    if (force.getNumGlobalParameters() > 0) {
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customAngleGlobals", false, CL_MEM_READ_ONLY);
        extraArguments += ", __constant float* globals";
    }
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<vector<cl_float> > paramVector(numAngles);
    vector<mm_int8> indicesVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        int particle1, particle2, particle3;
        vector<double> parameters;
        force.getAngleParameters(i, particle1, particle2, particle3, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
        indicesVector[i] = mm_int8(particle1, particle2, particle3, forceBufferCounter[particle1]++,
                forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, 0, 0);
    }
    params->setParameterValues(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLCustomAngleForceInfo(maxBuffers, force));

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("theta").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["float dEdAngle = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["theta"] = "theta";
    for (int i = 0; i < force.getNumPerAngleParameters(); i++) {
        const string& name = force.getPerAngleParameterName(i);
        variables[name] = "angleParams"+params->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+intToString(i)+"]";
        variables[name] = value;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArguments += ", __global "+buffer.getType()+"* "+buffer.getName();
        compute<<buffer.getType()<<" angleParams"<<(i+1)<<" = "<<buffer.getName()<<"[index];\n";
    }
    vector<pair<string, string> > functions;
    compute << OpenCLExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["EXTRA_ARGUMENTS"] = extraArguments;
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customAngleForce, replacements));
    kernel = cl::Kernel(program, "computeCustomAngleForces");
}

void OpenCLCalcCustomAngleForceKernel::executeForces(ContextImpl& context) {
    if (numAngles == 0)
        return;
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numAngles);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, indices->getDeviceBuffer());
        int nextIndex = 6;
        if (globals != NULL)
            kernel.setArg<cl::Buffer>(nextIndex++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            kernel.setArg<cl::Memory>(nextIndex++, buffer.getMemory());
        }
    }
    cl.executeKernel(kernel, numAngles);
}

double OpenCLCalcCustomAngleForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLPeriodicTorsionForceInfo : public OpenCLForceInfo {
public:
    OpenCLPeriodicTorsionForceInfo(int requiredBuffers, const PeriodicTorsionForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, periodicity2, phase2, k2);
        return (periodicity1 == periodicity2 && phase1 == phase2 && k1 == k2);
    }
private:
    const PeriodicTorsionForce& force;
};

OpenCLCalcPeriodicTorsionForceKernel::~OpenCLCalcPeriodicTorsionForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    if (numTorsions == 0)
        return;
    params = new OpenCLArray<mm_float4>(cl, numTorsions, "periodicTorsionParams");
    indices = new OpenCLArray<mm_int8>(cl, numTorsions, "periodicTorsionIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float4> paramVector(numTorsions);
    vector<mm_int8> indicesVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        paramVector[i] = mm_float4((cl_float) k, (cl_float) phase, (cl_float) periodicity, 0.0f);
        indicesVector[i] = mm_int8(particle1, particle2, particle3, particle4,
                forceBufferCounter[particle1]++, forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, forceBufferCounter[particle4]++);

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLPeriodicTorsionForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(OpenCLKernelSources::periodicTorsionForce);
    kernel = cl::Kernel(program, "calcPeriodicTorsionForce");
}

void OpenCLCalcPeriodicTorsionForceKernel::executeForces(ContextImpl& context) {
    if (numTorsions == 0)
        return;
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numTorsions);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
        kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    }
    cl.executeKernel(kernel, numTorsions);
}

double OpenCLCalcPeriodicTorsionForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLRBTorsionForceInfo : public OpenCLForceInfo {
public:
    OpenCLRBTorsionForceInfo(int requiredBuffers, const RBTorsionForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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
        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, c0b, c1b, c2b, c3b, c4b, c5b);
        return (c0a == c0b && c1a == c1b && c2a == c2b && c3a == c3b && c4a == c4b && c5a == c5b);
    }
private:
    const RBTorsionForce& force;
};

OpenCLCalcRBTorsionForceKernel::~OpenCLCalcRBTorsionForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
}

void OpenCLCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    if (numTorsions == 0)
        return;
    params = new OpenCLArray<mm_float8>(cl, numTorsions, "rbTorsionParams");
    indices = new OpenCLArray<mm_int8>(cl, numTorsions, "rbTorsionIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<mm_float8> paramVector(numTorsions);
    vector<mm_int8> indicesVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        paramVector[i] = mm_float8((cl_float) c0, (cl_float) c1, (cl_float) c2, (cl_float) c3, (cl_float) c4, (cl_float) c5, 0.0f, 0.0f);
        indicesVector[i] = mm_int8(particle1, particle2, particle3, particle4,
                forceBufferCounter[particle1]++, forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, forceBufferCounter[particle4]++);

    }
    params->upload(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLRBTorsionForceInfo(maxBuffers, force));
    cl::Program program = cl.createProgram(OpenCLKernelSources::rbTorsionForce);
    kernel = cl::Kernel(program, "calcRBTorsionForce");
}

void OpenCLCalcRBTorsionForceKernel::executeForces(ContextImpl& context) {
    if (numTorsions == 0)
        return;
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numTorsions);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, params->getDeviceBuffer());
        kernel.setArg<cl::Buffer>(6, indices->getDeviceBuffer());
    }
    cl.executeKernel(kernel, numTorsions);
}

double OpenCLCalcRBTorsionForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLCustomTorsionForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomTorsionForceInfo(int requiredBuffers, const CustomTorsionForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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

OpenCLCalcCustomTorsionForceKernel::~OpenCLCalcCustomTorsionForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
    if (globals != NULL)
        delete globals;
}

void OpenCLCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    if (numTorsions == 0)
        return;
    params = new OpenCLParameterSet(cl, force.getNumPerTorsionParameters(), numTorsions, "customTorsionParams");
    indices = new OpenCLArray<mm_int8>(cl, numTorsions, "customTorsionIndices");
    string extraArguments;
    if (force.getNumGlobalParameters() > 0) {
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customTorsionGlobals", false, CL_MEM_READ_ONLY);
        extraArguments += ", __constant float* globals";
    }
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    vector<vector<cl_float> > paramVector(numTorsions);
    vector<mm_int8> indicesVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int particle1, particle2, particle3, particle4;
        vector<double> parameters;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
        indicesVector[i] = mm_int8(particle1, particle2, particle3, particle4, forceBufferCounter[particle1]++,
                forceBufferCounter[particle2]++, forceBufferCounter[particle3]++, forceBufferCounter[particle4]++);
    }
    params->setParameterValues(paramVector);
    indices->upload(indicesVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    cl.addForce(new OpenCLCustomTorsionForceInfo(maxBuffers, force));

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("theta").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["float dEdAngle = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["theta"] = "theta";
    for (int i = 0; i < force.getNumPerTorsionParameters(); i++) {
        const string& name = force.getPerTorsionParameterName(i);
        variables[name] = "torsionParams"+params->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+intToString(i)+"]";
        variables[name] = value;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArguments += ", __global "+buffer.getType()+"* "+buffer.getName();
        compute<<buffer.getType()<<" torsionParams"<<(i+1)<<" = "<<buffer.getName()<<"[index];\n";
    }
    vector<pair<string, string> > functions;
    compute << OpenCLExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["EXTRA_ARGUMENTS"] = extraArguments;
    replacements["M_PI"] = doubleToString(M_PI);
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customTorsionForce, replacements));
    kernel = cl::Kernel(program, "computeCustomTorsionForces");
}

void OpenCLCalcCustomTorsionForceKernel::executeForces(ContextImpl& context) {
    if (numTorsions == 0)
        return;
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        kernel.setArg<cl_int>(1, numTorsions);
        kernel.setArg<cl::Buffer>(2, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(5, indices->getDeviceBuffer());
        int nextIndex = 6;
        if (globals != NULL)
            kernel.setArg<cl::Buffer>(nextIndex++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            kernel.setArg<cl::Memory>(nextIndex++, buffer.getMemory());
        }
    }
    cl.executeKernel(kernel, numTorsions);
}

double OpenCLCalcCustomTorsionForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLNonbondedForceInfo : public OpenCLForceInfo {
public:
    OpenCLNonbondedForceInfo(int requiredBuffers, const NonbondedForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, sigma1, sigma2, epsilon1, epsilon2;
        force.getParticleParameters(particle1, charge1, sigma1, epsilon1);
        force.getParticleParameters(particle2, charge2, sigma2, epsilon2);
        return (charge1 == charge2 && sigma1 == sigma2 && epsilon1 == epsilon2);
    }
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double chargeProd1, chargeProd2, sigma1, sigma2, epsilon1, epsilon2;
        force.getExceptionParameters(group1, particle1, particle2, chargeProd1, sigma1, epsilon1);
        force.getExceptionParameters(group2, particle1, particle2, chargeProd2, sigma2, epsilon2);
        return (chargeProd1 == chargeProd2 && sigma1 == sigma2 && epsilon1 == epsilon2);
    }
private:
    const NonbondedForce& force;
};

OpenCLCalcNonbondedForceKernel::~OpenCLCalcNonbondedForceKernel() {
    if (sigmaEpsilon != NULL)
        delete sigmaEpsilon;
    if (exceptionParams != NULL)
        delete exceptionParams;
    if (exceptionIndices != NULL)
        delete exceptionIndices;
    if (cosSinSums != NULL)
        delete cosSinSums;
    if (pmeGrid != NULL)
        delete pmeGrid;
    if (pmeBsplineModuliX != NULL)
        delete pmeBsplineModuliX;
    if (pmeBsplineModuliY != NULL)
        delete pmeBsplineModuliY;
    if (pmeBsplineModuliZ != NULL)
        delete pmeBsplineModuliZ;
    if (pmeBsplineTheta != NULL)
        delete pmeBsplineTheta;
    if (pmeBsplineDtheta != NULL)
        delete pmeBsplineDtheta;
    if (pmeAtomRange != NULL)
        delete pmeAtomRange;
    if (pmeAtomGridIndex != NULL)
        delete pmeAtomGridIndex;
    if (erfcTable != NULL)
        delete erfcTable;
    if (sort != NULL)
        delete sort;
    if (fft != NULL)
        delete fft;
}

void OpenCLCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {

    // Identify which exceptions are 1-4 interactions.

    vector<pair<int, int> > exclusions;
    vector<int> exceptions;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions.push_back(pair<int, int>(particle1, particle2));
        if (chargeProd != 0.0 || epsilon != 0.0)
            exceptions.push_back(i);
    }

    // Initialize nonbonded interactions.

    int numParticles = force.getNumParticles();
    sigmaEpsilon = new OpenCLArray<mm_float2>(cl, numParticles, "sigmaEpsilon");
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    vector<mm_float2> sigmaEpsilonVector(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    double sumSquaredCharges = 0.0;
    bool hasCoulomb = false;
    bool hasLJ = false;
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        posq[i].w = (float) charge;
        sigmaEpsilonVector[i] = mm_float2((float) (0.5*sigma), (float) (2.0*sqrt(epsilon)));
        exclusionList[i].push_back(i);
        sumSquaredCharges += charge*charge;
        if (charge != 0.0)
            hasCoulomb = true;
        if (epsilon != 0.0)
            hasLJ = true;
    }
    for (int i = 0; i < (int) exclusions.size(); i++) {
        exclusionList[exclusions[i].first].push_back(exclusions[i].second);
        exclusionList[exclusions[i].second].push_back(exclusions[i].first);
    }
    posq.upload();
    sigmaEpsilon->upload(sigmaEpsilonVector);
    bool useCutoff = (force.getNonbondedMethod() != NonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != NonbondedForce::NoCutoff && force.getNonbondedMethod() != NonbondedForce::CutoffNonPeriodic);
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    map<string, string> defines;
    defines["HAS_COULOMB"] = (hasCoulomb ? "1" : "0");
    defines["HAS_LENNARD_JONES"] = (hasLJ ? "1" : "0");
    if (useCutoff) {
        // Compute the reaction field constants.

        double reactionFieldK = pow(force.getCutoffDistance(), -3.0)*(force.getReactionFieldDielectric()-1.0)/(2.0*force.getReactionFieldDielectric()+1.0);
        double reactionFieldC = (1.0 / force.getCutoffDistance())*(3.0*force.getReactionFieldDielectric())/(2.0*force.getReactionFieldDielectric()+1.0);
        defines["REACTION_FIELD_K"] = doubleToString(reactionFieldK);
        defines["REACTION_FIELD_C"] = doubleToString(reactionFieldC);
    }
    double alpha = 0;
    if (force.getNonbondedMethod() == NonbondedForce::Ewald) {
        // Compute the Ewald parameters.

        int kmaxx, kmaxy, kmaxz;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
        defines["EWALD_ALPHA"] = doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        ewaldSelfEnergy = -ONE_4PI_EPS0*alpha*sumSquaredCharges/std::sqrt(M_PI);

        // Create the reciprocal space kernels.

        map<string, string> replacements;
        replacements["NUM_ATOMS"] = intToString(numParticles);
        replacements["KMAX_X"] = intToString(kmaxx);
        replacements["KMAX_Y"] = intToString(kmaxy);
        replacements["KMAX_Z"] = intToString(kmaxz);
        replacements["RECIPROCAL_BOX_SIZE_X"] = doubleToString(2.0*M_PI/boxVectors[0][0]);
        replacements["RECIPROCAL_BOX_SIZE_Y"] = doubleToString(2.0*M_PI/boxVectors[1][1]);
        replacements["RECIPROCAL_BOX_SIZE_Z"] = doubleToString(2.0*M_PI/boxVectors[2][2]);
        replacements["RECIPROCAL_COEFFICIENT"] = doubleToString(ONE_4PI_EPS0*4*M_PI/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]));
        replacements["EXP_COEFFICIENT"] = doubleToString(-1.0/(4.0*alpha*alpha));
        cl::Program program = cl.createProgram(OpenCLKernelSources::ewald, replacements);
        ewaldSumsKernel = cl::Kernel(program, "calculateEwaldCosSinSums");
        ewaldForcesKernel = cl::Kernel(program, "calculateEwaldForces");
        cosSinSums = new OpenCLArray<mm_float2>(cl, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), "cosSinSums");
    }
    else if (force.getNonbondedMethod() == NonbondedForce::PME) {
        // Compute the PME parameters.

        int gridSizeX, gridSizeY, gridSizeZ;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ);
        gridSizeX = OpenCLFFT3D::findLegalDimension(gridSizeX);
        gridSizeY = OpenCLFFT3D::findLegalDimension(gridSizeY);
        gridSizeZ = OpenCLFFT3D::findLegalDimension(gridSizeZ);
        defines["EWALD_ALPHA"] = doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        ewaldSelfEnergy = -ONE_4PI_EPS0*alpha*sumSquaredCharges/std::sqrt(M_PI);
        pmeDefines["PME_ORDER"] = intToString(PmeOrder);
        pmeDefines["NUM_ATOMS"] = intToString(numParticles);
        pmeDefines["RECIP_EXP_FACTOR"] = doubleToString(M_PI*M_PI/(alpha*alpha));
        pmeDefines["GRID_SIZE_X"] = intToString(gridSizeX);
        pmeDefines["GRID_SIZE_Y"] = intToString(gridSizeY);
        pmeDefines["GRID_SIZE_Z"] = intToString(gridSizeZ);
        pmeDefines["EPSILON_FACTOR"] = doubleToString(std::sqrt(ONE_4PI_EPS0));

        // Create required data structures.

        pmeGrid = new OpenCLArray<mm_float2>(cl, gridSizeX*gridSizeY*gridSizeZ, "pmeGrid");
        pmeBsplineModuliX = new OpenCLArray<cl_float>(cl, gridSizeX, "pmeBsplineModuliX");
        pmeBsplineModuliY = new OpenCLArray<cl_float>(cl, gridSizeY, "pmeBsplineModuliY");
        pmeBsplineModuliZ = new OpenCLArray<cl_float>(cl, gridSizeZ, "pmeBsplineModuliZ");
        pmeBsplineTheta = new OpenCLArray<mm_float4>(cl, PmeOrder*numParticles, "pmeBsplineTheta");
        pmeBsplineDtheta = new OpenCLArray<mm_float4>(cl, PmeOrder*numParticles, "pmeBsplineDtheta");
        pmeAtomRange = new OpenCLArray<cl_int>(cl, gridSizeX*gridSizeY*gridSizeZ+1, "pmeAtomRange");
        pmeAtomGridIndex = new OpenCLArray<mm_float2>(cl, numParticles, "pmeAtomGridIndex");
        sort = new OpenCLSort<mm_float2>(cl, cl.getNumAtoms(), "float2", "value.y");
        fft = new OpenCLFFT3D(cl, gridSizeX, gridSizeY, gridSizeZ);

        // Initialize the b-spline moduli.

        int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
        vector<double> data(PmeOrder);
        vector<double> ddata(PmeOrder);
        vector<double> bsplines_data(maxSize);
        data[PmeOrder-1] = 0.0;
        data[1] = 0.0;
        data[0] = 1.0;
        for (int i = 3; i < PmeOrder; i++) {
            double div = 1.0/(i-1.0);
            data[i-1] = 0.0;
            for (int j = 1; j < (i-1); j++)
                data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
            data[0] = div*data[0];
        }

        // Differentiate.

        ddata[0] = -data[0];
        for (int i = 1; i < PmeOrder; i++)
            ddata[i] = data[i-1]-data[i];
        double div = 1.0/(PmeOrder-1);
        data[PmeOrder-1] = 0.0;
        for (int i = 1; i < (PmeOrder-1); i++)
            data[PmeOrder-i-1] = div*(i*data[PmeOrder-i-2]+(PmeOrder-i)*data[PmeOrder-i-1]);
        data[0] = div*data[0];
        for (int i = 0; i < maxSize; i++)
            bsplines_data[i] = 0.0;
        for (int i = 1; i <= PmeOrder; i++)
            bsplines_data[i] = data[i-1];

        // Evaluate the actual bspline moduli for X/Y/Z.

        for(int dim = 0; dim < 3; dim++) {
            int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
            vector<cl_float> moduli(ndata);
            for (int i = 0; i < ndata; i++) {
                double sc = 0.0;
                double ss = 0.0;
                for (int j = 0; j < ndata; j++) {
                    double arg = (2.0*M_PI*i*j)/ndata;
                    sc += bsplines_data[j]*cos(arg);
                    ss += bsplines_data[j]*sin(arg);
                }
                moduli[i] = (float) (sc*sc+ss*ss);
            }
            for (int i = 0; i < ndata; i++)
            {
                if (moduli[i] < 1.0e-7)
                    moduli[i] = (moduli[i-1]+moduli[i+1])*0.5f;
            }
            if (dim == 0)
                pmeBsplineModuliX->upload(moduli);
            else if (dim == 1)
                pmeBsplineModuliY->upload(moduli);
            else
                pmeBsplineModuliZ->upload(moduli);
        }
    }
    else
        ewaldSelfEnergy = 0.0;
    
    // Tabulate values of erfc().

    if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
        const int tableSize = 2048;
        defines["ERFC_TABLE_SCALE"] = doubleToString((tableSize-1)/(alpha*force.getCutoffDistance()));
        erfcTable = new OpenCLArray<cl_float>(cl, tableSize, "ErfcTable", false, CL_MEM_READ_ONLY);
        vector<cl_float> erfcVector(tableSize);
        for (int i = 0; i < tableSize; ++i)
            erfcVector[i] = (float) erfc(i*(alpha*force.getCutoffDistance())/(tableSize-1));
        erfcTable->upload(erfcVector);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo("erfcTable", "float", 1, sizeof(cl_float), erfcTable->getDeviceBuffer()));
    }

    // Add the interaction to the default nonbonded kernel.
    
    string source = cl.replaceStrings(OpenCLKernelSources::coulombLennardJones, defines);
    cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source);
    if (hasLJ)
        cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo("sigmaEpsilon", "float", 2, sizeof(cl_float2), sigmaEpsilon->getDeviceBuffer()));

    // Initialize the exceptions.

    int numExceptions = exceptions.size();
    int maxBuffers = cl.getNonbondedUtilities().getNumForceBuffers();
    if (numExceptions > 0) {
        exceptionParams = new OpenCLArray<mm_float4>(cl, numExceptions, "exceptionParams");
        exceptionIndices = new OpenCLArray<mm_int4>(cl, numExceptions, "exceptionIndices");
        vector<mm_float4> exceptionParamsVector(numExceptions);
        vector<mm_int4> exceptionIndicesVector(numExceptions);
        vector<int> forceBufferCounter(system.getNumParticles(), 0);
        for (int i = 0; i < numExceptions; i++) {
            int particle1, particle2;
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[i], particle1, particle2, chargeProd, sigma, epsilon);
            exceptionParamsVector[i] = mm_float4((float) (ONE_4PI_EPS0*chargeProd), (float) sigma, (float) (4.0*epsilon), 0.0f);
            exceptionIndicesVector[i] = mm_int4(particle1, particle2, forceBufferCounter[particle1]++, forceBufferCounter[particle2]++);
        }
        exceptionParams->upload(exceptionParamsVector);
        exceptionIndices->upload(exceptionIndicesVector);
        for (int i = 0; i < (int) forceBufferCounter.size(); i++)
            maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    }
    cl.addForce(new OpenCLNonbondedForceInfo(maxBuffers, force));
    defines.clear();
    defines["NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
    defines["NUM_EXCEPTIONS"] = intToString(numExceptions);
    cl::Program program = cl.createProgram(OpenCLKernelSources::nonbondedExceptions, defines);
    exceptionsKernel = cl::Kernel(program, "computeNonbondedExceptions");
}

void OpenCLCalcNonbondedForceKernel::executeForces(ContextImpl& context) {
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        if (exceptionIndices != NULL) {
            exceptionsKernel.setArg<cl::Buffer>(0, cl.getForceBuffers().getDeviceBuffer());
            exceptionsKernel.setArg<cl::Buffer>(1, cl.getEnergyBuffer().getDeviceBuffer());
            exceptionsKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
            exceptionsKernel.setArg<cl::Buffer>(3, exceptionParams->getDeviceBuffer());
            exceptionsKernel.setArg<cl::Buffer>(4, exceptionIndices->getDeviceBuffer());
        }
        if (cosSinSums != NULL) {
            ewaldSumsKernel.setArg<cl::Buffer>(0, cl.getEnergyBuffer().getDeviceBuffer());
            ewaldSumsKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            ewaldSumsKernel.setArg<cl::Buffer>(2, cosSinSums->getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(0, cl.getForceBuffers().getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(2, cosSinSums->getDeviceBuffer());
        }
        if (pmeGrid != NULL) {
            mm_float4 boxSize = cl.getNonbondedUtilities().getPeriodicBoxSize();
            pmeDefines["PERIODIC_BOX_SIZE_X"] = doubleToString(boxSize.x);
            pmeDefines["PERIODIC_BOX_SIZE_Y"] = doubleToString(boxSize.y);
            pmeDefines["PERIODIC_BOX_SIZE_Z"] = doubleToString(boxSize.z);
            pmeDefines["INV_PERIODIC_BOX_SIZE_X"] = doubleToString(1.0/boxSize.x);
            pmeDefines["INV_PERIODIC_BOX_SIZE_Y"] = doubleToString(1.0/boxSize.y);
            pmeDefines["INV_PERIODIC_BOX_SIZE_Z"] = doubleToString(1.0/boxSize.z);
            pmeDefines["RECIP_SCALE_FACTOR"] = doubleToString(1.0/(M_PI*boxSize.x*boxSize.y*boxSize.z));
            cl::Program program = cl.createProgram(OpenCLKernelSources::pme, pmeDefines);
            pmeGridIndexKernel = cl::Kernel(program, "updateGridIndexAndFraction");
            pmeAtomRangeKernel = cl::Kernel(program, "findAtomRangeForGrid");
            pmeUpdateBsplinesKernel = cl::Kernel(program, "updateBsplines");
            pmeSpreadChargeKernel = cl::Kernel(program, "gridSpreadCharge");
            pmeConvolutionKernel = cl::Kernel(program, "reciprocalConvolution");
            pmeInterpolateForceKernel = cl::Kernel(program, "gridInterpolateForce");
            pmeGridIndexKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeGridIndexKernel.setArg<cl::Buffer>(1, pmeAtomGridIndex->getDeviceBuffer());
            pmeAtomRangeKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeAtomRangeKernel.setArg<cl::Buffer>(1, pmeAtomGridIndex->getDeviceBuffer());
            pmeAtomRangeKernel.setArg<cl::Buffer>(2, pmeAtomRange->getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(1, pmeBsplineTheta->getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(2, pmeBsplineDtheta->getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg(3, 2*OpenCLContext::ThreadBlockSize*PmeOrder*sizeof(mm_float4), NULL);
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(4, pmeAtomGridIndex->getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(0, pmeAtomGridIndex->getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(1, pmeAtomRange->getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(2, pmeGrid->getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(3, pmeBsplineTheta->getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(0, pmeGrid->getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(1, cl.getEnergyBuffer().getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(2, pmeBsplineModuliX->getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(3, pmeBsplineModuliY->getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(4, pmeBsplineModuliZ->getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(1, cl.getForceBuffers().getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(2, pmeBsplineTheta->getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(3, pmeBsplineDtheta->getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(4, pmeGrid->getDeviceBuffer());
       }
    }
    if (exceptionIndices != NULL)
        cl.executeKernel(exceptionsKernel, exceptionIndices->getSize());
    if (cosSinSums != NULL) {
        cl.executeKernel(ewaldSumsKernel, cosSinSums->getSize());
        cl.executeKernel(ewaldForcesKernel, cl.getNumAtoms());
    }
    if (pmeGrid != NULL) {
        cl.executeKernel(pmeGridIndexKernel, cl.getNumAtoms());
        sort->sort(*pmeAtomGridIndex);
        cl.executeKernel(pmeAtomRangeKernel, cl.getNumAtoms());
        cl.executeKernel(pmeUpdateBsplinesKernel, cl.getNumAtoms());
        cl.executeKernel(pmeSpreadChargeKernel, cl.getNumAtoms());
        fft->execFFT(*pmeGrid, true);
        cl.executeKernel(pmeConvolutionKernel, cl.getNumAtoms());
        fft->execFFT(*pmeGrid, false);
        cl.executeKernel(pmeInterpolateForceKernel, cl.getNumAtoms());
    }
}

double OpenCLCalcNonbondedForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return ewaldSelfEnergy;
}

class OpenCLCustomNonbondedForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomNonbondedForceInfo(int requiredBuffers, const CustomNonbondedForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        vector<double> params1;
        vector<double> params2;
        force.getParticleParameters(particle1, params1);
        force.getParticleParameters(particle2, params2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExclusions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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
    const CustomNonbondedForce& force;
};

OpenCLCalcCustomNonbondedForceKernel::~OpenCLCalcCustomNonbondedForceKernel() {
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
    if (tabulatedFunctionParams != NULL)
        delete tabulatedFunctionParams;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
}

void OpenCLCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "custom"+intToString(forceIndex)+"_";

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), numParticles, "customNonbondedParameters");
    if (force.getNumGlobalParameters() > 0)
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customNonbondedGlobals", false, CL_MEM_READ_ONLY);
    vector<vector<cl_float> > paramVector(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
        exclusionList[i].push_back(i);
    }
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusionList[particle1].push_back(particle2);
        exclusionList[particle2].push_back(particle1);
    }
    params->setParameterValues(paramVector);

    // Record the tabulated functions.

    OpenCLExpressionUtilities::FunctionPlaceholder fp;
    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        bool interpolating;
        force.getFunctionParameters(i, name, values, min, max, interpolating);
        string arrayName = prefix+"table"+intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = &fp;
        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), 0.0f);
        vector<mm_float4> f = OpenCLExpressionUtilities::computeFunctionCoefficients(values, interpolating);
        tabulatedFunctions.push_back(new OpenCLArray<mm_float4>(cl, values.size()-1, "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(arrayName, "float", 4, sizeof(cl_float4), tabulatedFunctions[tabulatedFunctions.size()-1]->getDeviceBuffer()));
    }
    if (force.getNumFunctions() > 0) {
        tabulatedFunctionParams = new OpenCLArray<mm_float4>(cl, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(prefix+"functionParams", "float", 4, sizeof(cl_float4), tabulatedFunctionParams->getDeviceBuffer()));
    }

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);
    bool useCutoff = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && force.getNonbondedMethod() != CustomNonbondedForce::CutoffNonPeriodic);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
    map<string, Lepton::ParsedExpression> forceExpressions;
    forceExpressions["tempEnergy += "] = energyExpression;
    forceExpressions["tempForce -= "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["r"] = "r";
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        variables[name+"1"] = prefix+"params"+params->getParameterSuffix(i, "1");
        variables[name+"2"] = prefix+"params"+params->getParameterSuffix(i, "2");
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+intToString(i)+"]";
        variables[name] = prefix+value;
    }
    stringstream compute;
    compute << OpenCLExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, prefix+"temp", prefix+"functionParams");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    string source = cl.replaceStrings(OpenCLKernelSources::customNonbonded, replacements);
    cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source);
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"params"+intToString(i+1), buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
    }
    if (globals != NULL) {
        globals->upload(globalParamValues);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(cl_float), globals->getDeviceBuffer()));
    }
    cl.addForce(new OpenCLCustomNonbondedForceInfo(cl.getNonbondedUtilities().getNumForceBuffers(), force));
}

void OpenCLCalcCustomNonbondedForceKernel::executeForces(ContextImpl& context) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
}

double OpenCLCalcCustomNonbondedForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLGBSAOBCForceInfo : public OpenCLForceInfo {
public:
    OpenCLGBSAOBCForceInfo(int requiredBuffers, const GBSAOBCForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, scale1, scale2;
        force.getParticleParameters(particle1, charge1, radius1, scale1);
        force.getParticleParameters(particle2, charge2, radius2, scale2);
        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2);
    }
private:
    const GBSAOBCForce& force;
};

OpenCLCalcGBSAOBCForceKernel::~OpenCLCalcGBSAOBCForceKernel() {
    if (params != NULL)
        delete params;
    if (bornSum != NULL)
        delete bornSum;
    if (bornRadii != NULL)
        delete bornRadii;
    if (bornForce != NULL)
        delete bornForce;
    if (obcChain != NULL)
        delete obcChain;
}

void OpenCLCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    params = new OpenCLArray<mm_float2>(cl, cl.getPaddedNumAtoms(), "gbsaObcParams");
    bornRadii = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms(), "bornRadii");
    obcChain = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms(), "obcChain");
    bornSum = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), "bornSum");
    bornForce = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), "bornForce");
    OpenCLArray<mm_float4>& posq = cl.getPosq();
    int numParticles = force.getNumParticles();
    vector<mm_float2> paramsVector(numParticles);
    const double dielectricOffset = 0.009;
    for (int i = 0; i < numParticles; i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
        posq[i].w = (float) charge;
    }
    posq.upload();
    params->upload(paramsVector);
    prefactor = -ONE_4PI_EPS0*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
    string source = OpenCLKernelSources::gbsaObc2;
    nb.addInteraction(useCutoff, usePeriodic, false, force.getCutoffDistance(), vector<vector<int> >(), source);
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo("obcParams", "float", 2, sizeof(cl_float2), params->getDeviceBuffer()));;
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo("bornForce", "float", 1, sizeof(cl_float), bornForce->getDeviceBuffer()));;
    cl.addForce(new OpenCLGBSAOBCForceInfo(nb.getNumForceBuffers(), force));
}

void OpenCLCalcGBSAOBCForceKernel::executeForces(ContextImpl& context) {
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    if (!hasCreatedKernels) {
        // These Kernels cannot be created in initialize(), because the OpenCLNonbondedUtilities has not been initialized yet then.

        hasCreatedKernels = true;
        map<string, string> defines;
        if (nb.getForceBufferPerAtomBlock())
            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
        if (nb.getUseCutoff())
            defines["USE_CUTOFF"] = "1";
        if (nb.getUsePeriodic())
            defines["USE_PERIODIC"] = "1";
        defines["PERIODIC_BOX_SIZE_X"] = doubleToString(nb.getPeriodicBoxSize().x);
        defines["PERIODIC_BOX_SIZE_Y"] = doubleToString(nb.getPeriodicBoxSize().y);
        defines["PERIODIC_BOX_SIZE_Z"] = doubleToString(nb.getPeriodicBoxSize().z);
        defines["INV_PERIODIC_BOX_SIZE_X"] = doubleToString(1.0/nb.getPeriodicBoxSize().x);
        defines["INV_PERIODIC_BOX_SIZE_Y"] = doubleToString(1.0/nb.getPeriodicBoxSize().y);
        defines["INV_PERIODIC_BOX_SIZE_Z"] = doubleToString(1.0/nb.getPeriodicBoxSize().z);
        defines["CUTOFF_SQUARED"] = doubleToString(nb.getCutoffDistance()*nb.getCutoffDistance());
        defines["PREFACTOR"] = doubleToString(prefactor);
        defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
        string file = (cl.getSIMDWidth() == 32 ? OpenCLKernelSources::gbsaObc_nvidia : OpenCLKernelSources::gbsaObc_default);
        cl::Program program = cl.createProgram(file, defines);
        int index = 0;
        computeBornSumKernel = cl::Kernel(program, "computeBornSum");
        computeBornSumKernel.setArg<cl::Buffer>(index++, bornSum->getDeviceBuffer());
        computeBornSumKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        computeBornSumKernel.setArg<cl::Buffer>(index++, params->getDeviceBuffer());
        computeBornSumKernel.setArg(index++, OpenCLContext::ThreadBlockSize*13*sizeof(cl_float), NULL);
        computeBornSumKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        if (nb.getUseCutoff()) {
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
        }
        else {
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getTiles().getDeviceBuffer());
            computeBornSumKernel.setArg<cl_uint>(index++, nb.getTiles().getSize());
        }
        force1Kernel = cl::Kernel(program, "computeGBSAForce1");
        index = 0;
        force1Kernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, bornRadii->getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, bornForce->getDeviceBuffer());
        force1Kernel.setArg(index++, OpenCLContext::ThreadBlockSize*13*sizeof(cl_float), NULL);
        force1Kernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
        if (nb.getUseCutoff()) {
            force1Kernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
        }
        else {
            force1Kernel.setArg<cl::Buffer>(index++, nb.getTiles().getDeviceBuffer());
            force1Kernel.setArg<cl_uint>(index++, nb.getTiles().getSize());
        }
        program = cl.createProgram(OpenCLKernelSources::gbsaObcReductions, defines);
        reduceBornSumKernel = cl::Kernel(program, "reduceBornSum");
        reduceBornSumKernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        reduceBornSumKernel.setArg<cl_int>(1, nb.getNumForceBuffers());
        reduceBornSumKernel.setArg<cl_float>(2, 1.0f);
        reduceBornSumKernel.setArg<cl_float>(3, 0.8f);
        reduceBornSumKernel.setArg<cl_float>(4, 4.85f);
        reduceBornSumKernel.setArg<cl::Buffer>(5, bornSum->getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(6, params->getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(7, bornRadii->getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(8, obcChain->getDeviceBuffer());
        reduceBornForceKernel = cl::Kernel(program, "reduceBornForce");
        reduceBornForceKernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        reduceBornForceKernel.setArg<cl_int>(1, nb.getNumForceBuffers());
        reduceBornForceKernel.setArg<cl::Buffer>(2, bornForce->getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(3, cl.getEnergyBuffer().getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(4, params->getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(5, bornRadii->getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(6, obcChain->getDeviceBuffer());
    }
    cl.clearBuffer(*bornSum);
    cl.clearBuffer(*bornForce);
    cl.executeKernel(computeBornSumKernel, nb.getTiles().getSize()*OpenCLContext::TileSize);
    cl.executeKernel(reduceBornSumKernel, cl.getPaddedNumAtoms());
    cl.executeKernel(force1Kernel, nb.getTiles().getSize()*OpenCLContext::TileSize);
    cl.executeKernel(reduceBornForceKernel, cl.getPaddedNumAtoms());
}

double OpenCLCalcGBSAOBCForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLCustomGBForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomGBForceInfo(int requiredBuffers, const CustomGBForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        vector<double> params1;
        vector<double> params2;
        force.getParticleParameters(particle1, params1);
        force.getParticleParameters(particle2, params2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExclusions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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
    const CustomGBForce& force;
};

OpenCLCalcCustomGBForceKernel::~OpenCLCalcCustomGBForceKernel() {
    if (params != NULL)
        delete params;
    if (computedValues != NULL)
        delete computedValues;
    if (energyDerivs != NULL)
        delete energyDerivs;
    if (globals != NULL)
        delete globals;
    if (valueBuffers != NULL)
        delete valueBuffers;
    if (tabulatedFunctionParams != NULL)
        delete tabulatedFunctionParams;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
}

void OpenCLCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    bool useExclusionsForValue = false;
    vector<string> computedValueNames(force.getNumComputedValues());
    vector<string> computedValueExpressions(force.getNumComputedValues());
    if (force.getNumComputedValues() > 0) {
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(0, computedValueNames[0], computedValueExpressions[0], type);
        if (type == CustomGBForce::SingleParticle)
            throw OpenMMException("OpenCLPlatform requires that the first computed value for a CustomGBForce be of type ParticlePair or ParticlePairNoExclusions.");
        useExclusionsForValue = (type == CustomGBForce::ParticlePair);
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            force.getComputedValueParameters(i, computedValueNames[i], computedValueExpressions[i], type);
            if (type != CustomGBForce::SingleParticle)
                throw OpenMMException("OpenCLPlatform requires that a CustomGBForce only have one computed value of type ParticlePair or ParticlePairNoExclusions.");
        }
    }
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "custom"+intToString(forceIndex)+"_";

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), numParticles, "customGBParameters");
    computedValues = new OpenCLParameterSet(cl, force.getNumComputedValues(), numParticles, "customGBComputedValues");
    if (force.getNumGlobalParameters() > 0)
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customGBGlobals", false, CL_MEM_READ_ONLY);
    vector<vector<cl_float> > paramVector(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
        exclusionList[i].push_back(i);
    }
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusionList[particle1].push_back(particle2);
        exclusionList[particle2].push_back(particle1);
    }
    params->setParameterValues(paramVector);

    // Record the tabulated functions.

    OpenCLExpressionUtilities::FunctionPlaceholder fp;
    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
    stringstream tableArgs;
    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        bool interpolating;
        force.getFunctionParameters(i, name, values, min, max, interpolating);
        string arrayName = prefix+"table"+intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = &fp;
        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), 0.0f);
        vector<mm_float4> f = OpenCLExpressionUtilities::computeFunctionCoefficients(values, interpolating);
        tabulatedFunctions.push_back(new OpenCLArray<mm_float4>(cl, values.size()-1, "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(arrayName, "float", 4, sizeof(cl_float4), tabulatedFunctions[tabulatedFunctions.size()-1]->getDeviceBuffer()));
        tableArgs << ", __global float4* " << arrayName;
    }
    if (force.getNumFunctions() > 0) {
        tabulatedFunctionParams = new OpenCLArray<mm_float4>(cl, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(prefix+"functionParams", "float", 4, sizeof(cl_float4), tabulatedFunctionParams->getDeviceBuffer()));
        tableArgs << ", __constant float4* " << prefix << "functionParams";
    }

    // Record the global parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);

    // Record derivatives of expressions needed for the chain rule terms.

    vector<vector<Lepton::ParsedExpression> > valueGradientExpressions(force.getNumComputedValues());
    bool needParameterGradient = false;
    for (int i = 1; i < force.getNumComputedValues(); i++) {
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
        valueGradientExpressions[i].push_back(ex.differentiate("x").optimize());
        valueGradientExpressions[i].push_back(ex.differentiate("y").optimize());
        valueGradientExpressions[i].push_back(ex.differentiate("z").optimize());
        if (!isZeroExpression(valueGradientExpressions[i][0]) || !isZeroExpression(valueGradientExpressions[i][1]) || !isZeroExpression(valueGradientExpressions[i][2]))
            needParameterGradient = true;
    }
    vector<vector<Lepton::ParsedExpression> > energyDerivExpressions(force.getNumEnergyTerms());
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        for (int j = 0; j < force.getNumComputedValues(); j++) {
            if (type == CustomGBForce::SingleParticle)
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"1").optimize());
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"2").optimize());
            }
        }
    }
    energyDerivs = new OpenCLParameterSet(cl, force.getNumComputedValues(), cl.getPaddedNumAtoms()*cl.getNonbondedUtilities().getNumForceBuffers(), "customGBEnergyDerivatives");

    // Create the kernels.

    bool useCutoff = (force.getNonbondedMethod() != CustomGBForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomGBForce::NoCutoff && force.getNonbondedMethod() != CustomGBForce::CutoffNonPeriodic);
    {
        // Create the N2 value kernel.

        map<string, string> variables;
        map<string, string> rename;
        variables["r"] = "r";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables[name+"1"] = "params"+params->getParameterSuffix(i, "1");
            variables[name+"2"] = "params"+params->getParameterSuffix(i, "2");
            rename[name+"1"] = name+"2";
            rename[name+"2"] = name+"1";
        }
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+intToString(i)+"]";
            variables[name] = value;
        }
        map<string, Lepton::ParsedExpression> n2ValueExpressions;
        stringstream n2ValueSource;
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[0], functions).optimize();
        n2ValueExpressions["tempValue1 = "] = ex;
        n2ValueExpressions["tempValue2 = "] = ex.renameVariables(rename);
        n2ValueSource << OpenCLExpressionUtilities::createExpressions(n2ValueExpressions, variables, functionDefinitions, "temp", prefix+"functionParams");
        map<string, string> replacements;
        replacements["COMPUTE_VALUE"] = n2ValueSource.str();
        stringstream extraArgs, loadLocal1, loadLocal2, load1, load2;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __constant float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* global_" << paramName << ", __local " << buffer.getType() << "* local_" << paramName;
            loadLocal1 << "local_" << paramName << "[get_local_id(0)] = " << paramName << "1;\n";
            loadLocal2 << "local_" << paramName << "[get_local_id(0)] = global_" << paramName << "[j];\n";
            load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
            load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        map<string, string> defines;
        if (cl.getNonbondedUtilities().getForceBufferPerAtomBlock())
            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
        if (useCutoff)
            defines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        if (useExclusionsForValue)
            defines["USE_EXCLUSIONS"] = "1";
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        defines["PERIODIC_BOX_SIZE_X"] = doubleToString(boxVectors[0][0]);
        defines["PERIODIC_BOX_SIZE_Y"] = doubleToString(boxVectors[1][1]);
        defines["PERIODIC_BOX_SIZE_Z"] = doubleToString(boxVectors[2][2]);
        defines["INV_PERIODIC_BOX_SIZE_X"] = doubleToString(1.0/boxVectors[0][0]);
        defines["INV_PERIODIC_BOX_SIZE_Y"] = doubleToString(1.0/boxVectors[1][1]);
        defines["INV_PERIODIC_BOX_SIZE_Z"] = doubleToString(1.0/boxVectors[2][2]);
        defines["CUTOFF_SQUARED"] = doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
        defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
        string file = (cl.getSIMDWidth() == 32 ? OpenCLKernelSources::customGBValueN2_nvidia : OpenCLKernelSources::customGBValueN2_default);
        cl::Program program = cl.createProgram(cl.replaceStrings(file, replacements), defines);
        pairValueKernel = cl::Kernel(program, "computeN2Value");
    }
    {
        // Create the kernel to reduce the N2 value and calculate other values.

        stringstream reductionSource, extraArgs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __constant float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* global_" << valueName;
            reductionSource << buffer.getType() << " local_" << valueName << ";\n";
        }
        reductionSource << "local_values" << computedValues->getParameterSuffix(0) << " = sum;\n";
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+intToString(i)+"]";
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            variables[computedValueNames[i-1]] = "local_values"+computedValues->getParameterSuffix(i-1);
            map<string, Lepton::ParsedExpression> valueExpressions;
            valueExpressions["local_values"+computedValues->getParameterSuffix(i)+" = "] = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
            reductionSource << OpenCLExpressionUtilities::createExpressions(valueExpressions, variables, functionDefinitions, "value"+intToString(i)+"_temp", "functionParams");
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            string valueName = "values"+intToString(i+1);
            reductionSource << "global_" << valueName << "[index] = local_" << valueName << ";\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["COMPUTE_VALUES"] = reductionSource.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
        cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customGBValuePerParticle, replacements), defines);
        perParticleValueKernel = cl::Kernel(program, "computePerParticleValues");
    }
    {
        // Create the N2 energy kernel.

        map<string, string> variables;
        variables["r"] = "r";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables[name+"1"] = "params"+params->getParameterSuffix(i, "1");
            variables[name+"2"] = "params"+params->getParameterSuffix(i, "2");
        }
        for (int i = 0; i < force.getNumComputedValues(); i++) {
            variables[computedValueNames[i]+"1"] = "values"+computedValues->getParameterSuffix(i, "1");
            variables[computedValueNames[i]+"2"] = "values"+computedValues->getParameterSuffix(i, "2");
        }
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+intToString(i)+"]";
        map<string, Lepton::ParsedExpression> n2EnergyExpressions;
        stringstream n2EnergySource;
        bool anyExclusions = false;
        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
            string expression;
            CustomGBForce::ComputationType type;
            force.getEnergyTermParameters(i, expression, type);
            if (type == CustomGBForce::SingleParticle)
                continue;
            bool exclude = (type == CustomGBForce::ParticlePair);
            anyExclusions |= exclude;
            n2EnergyExpressions["tempEnergy += "] = Lepton::Parser::parse(expression, functions).optimize();
            n2EnergyExpressions["dEdR += "] = Lepton::Parser::parse(expression, functions).differentiate("r").optimize();
            for (int j = 0; j < force.getNumComputedValues(); j++) {
                n2EnergyExpressions["/*"+intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j, "_1")+" += "] = energyDerivExpressions[i][2*j];
                n2EnergyExpressions["/*"+intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j, "_2")+" += "] = energyDerivExpressions[i][2*j+1];
            }
            if (exclude)
                n2EnergySource << "if (!isExcluded) {\n";
            n2EnergySource << OpenCLExpressionUtilities::createExpressions(n2EnergyExpressions, variables, functionDefinitions, "temp", prefix+"functionParams");
            if (exclude)
                n2EnergySource << "}\n";
        }
        map<string, string> replacements;
        replacements["COMPUTE_INTERACTION"] = n2EnergySource.str();
        stringstream extraArgs, loadLocal1, loadLocal2, clearLocal, load1, load2, recordDeriv, storeDerivs1, storeDerivs2, declareTemps, setTemps;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __constant float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* global_" << paramName << ", __local " << buffer.getType() << "* local_" << paramName;
            loadLocal1 << "local_" << paramName << "[get_local_id(0)] = " << paramName << "1;\n";
            loadLocal2 << "local_" << paramName << "[get_local_id(0)] = global_" << paramName << "[j];\n";
            load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
            load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* global_" << valueName << ", __local " << buffer.getType() << "* local_" << valueName;
            loadLocal1 << "local_" << valueName << "[get_local_id(0)] = " << valueName << "1;\n";
            loadLocal2 << "local_" << valueName << "[get_local_id(0)] = global_" << valueName << "[j];\n";
            load1 << buffer.getType() << " " << valueName << "1 = global_" << valueName << "[atom1];\n";
            load2 << buffer.getType() << " " << valueName << "2 = local_" << valueName << "[atom2];\n";
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string index = intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* derivBuffers" << index << ", __local " << buffer.getType() << "* local_deriv" << index;
            clearLocal << "local_deriv" << index << "[get_local_id(0)] = 0.0f;\n";
            load1 << buffer.getType() << " deriv" << index << "_1 = 0.0f;\n";
            load2 << buffer.getType() << " deriv" << index << "_2 = 0.0f;\n";
            recordDeriv << "local_deriv" << index << "[atom2] += deriv" << index << "_2;\n";
            storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")";
            storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")";
            declareTemps << "__local " << buffer.getType() << " tempDerivBuffer" << index << "[64];\n";
            setTemps << "tempDerivBuffer" << index << "[get_local_id(0)] = deriv" << index << "_1;\n";
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["CLEAR_LOCAL_DERIVATIVES"] = clearLocal.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        replacements["RECORD_DERIVATIVE_2"] = recordDeriv.str();
        replacements["STORE_DERIVATIVES_1"] = storeDerivs1.str();
        replacements["STORE_DERIVATIVES_2"] = storeDerivs2.str();
        replacements["DECLARE_TEMP_BUFFERS"] = declareTemps.str();
        replacements["SET_TEMP_BUFFERS"] = setTemps.str();
        map<string, string> defines;
        if (cl.getNonbondedUtilities().getForceBufferPerAtomBlock())
            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
        if (useCutoff)
            defines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        if (anyExclusions)
            defines["USE_EXCLUSIONS"] = "1";
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        defines["PERIODIC_BOX_SIZE_X"] = doubleToString(boxVectors[0][0]);
        defines["PERIODIC_BOX_SIZE_Y"] = doubleToString(boxVectors[1][1]);
        defines["PERIODIC_BOX_SIZE_Z"] = doubleToString(boxVectors[2][2]);
        defines["INV_PERIODIC_BOX_SIZE_X"] = doubleToString(1.0/boxVectors[0][0]);
        defines["INV_PERIODIC_BOX_SIZE_Y"] = doubleToString(1.0/boxVectors[1][1]);
        defines["INV_PERIODIC_BOX_SIZE_Z"] = doubleToString(1.0/boxVectors[2][2]);
        defines["CUTOFF_SQUARED"] = doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
        defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
        string file = (cl.getSIMDWidth() == 32 ? OpenCLKernelSources::customGBEnergyN2_nvidia : OpenCLKernelSources::customGBEnergyN2_default);
        cl::Program program = cl.createProgram(cl.replaceStrings(file, replacements), defines);
        pairEnergyKernel = cl::Kernel(program, "computeN2Energy");
    }
    {
        // Create the kernel to reduce the derivatives and calculate per-particle energy terms.

        stringstream compute, extraArgs, reduce;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __constant float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string index = intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* derivBuffers" << index;
            reduce << "REDUCE_VALUE(derivBuffers" << index << ", " << buffer.getType() << ")\n";
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+intToString(i)+"]";
        for (int i = 0; i < force.getNumComputedValues(); i++)
            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
        map<string, Lepton::ParsedExpression> energyExpressions;
        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
            string expression;
            CustomGBForce::ComputationType type;
            force.getEnergyTermParameters(i, expression, type);
            if (type != CustomGBForce::SingleParticle)
                continue;
            Lepton::ParsedExpression parsed = Lepton::Parser::parse(expression, functions).optimize();
            energyExpressions["/*"+intToString(i+1)+"*/ energy += "] = parsed;
            for (int j = 0; j < force.getNumComputedValues(); j++)
                energyExpressions["/*"+intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j)+" += "] = energyDerivExpressions[i][j];
            Lepton::ParsedExpression gradx = parsed.differentiate("x").optimize();
            Lepton::ParsedExpression grady = parsed.differentiate("y").optimize();
            Lepton::ParsedExpression gradz = parsed.differentiate("z").optimize();
            if (!isZeroExpression(gradx))
                energyExpressions["/*"+intToString(i+1)+"*/ force.x -= "] = gradx;
            if (!isZeroExpression(grady))
                energyExpressions["/*"+intToString(i+1)+"*/ force.y -= "] = grady;
            if (!isZeroExpression(gradz))
                energyExpressions["/*"+intToString(i+1)+"*/ force.z -= "] = gradz;
        }
        compute << OpenCLExpressionUtilities::createExpressions(energyExpressions, variables, functionDefinitions, "temp", prefix+"functionParams");
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            string index = intToString(i+1);
            compute << "derivBuffers" << index << "[index] = deriv" << index << ";\n";
        }
        compute << "forceBuffers[index] = forceBuffers[index]+force;\n";
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["REDUCE_DERIVATIVES"] = reduce.str();
        replacements["COMPUTE_ENERGY"] = compute.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
        cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customGBEnergyPerParticle, replacements), defines);
        perParticleEnergyKernel = cl::Kernel(program, "computePerParticleEnergy");
    }
    {
        // Create the code to calculate chain rules terms (possibly as part of the default nonbonded kernel).

        map<string, string> globalVariables;
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+intToString(i)+"]";
            globalVariables[name] = prefix+value;
        }
        map<string, string> variables = globalVariables;
        map<string, string> rename;
        variables["r"] = "r";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables[name+"1"] = prefix+"params"+params->getParameterSuffix(i, "1");
            variables[name+"2"] = prefix+"params"+params->getParameterSuffix(i, "2");
            rename[name+"1"] =  name+"2";
            rename[name+"2"] =  name+"1";
        }
        map<string, Lepton::ParsedExpression> derivExpressions;
        stringstream chainSource;
        Lepton::ParsedExpression dVdR = Lepton::Parser::parse(computedValueExpressions[0], functions).differentiate("r").optimize();
        derivExpressions["float dV0dR1 = "] = dVdR;
        derivExpressions["float dV0dR2 = "] = dVdR.renameVariables(rename);
        chainSource << OpenCLExpressionUtilities::createExpressions(derivExpressions, variables, functionDefinitions, prefix+"temp0_", prefix+"functionParams");
        if (needParameterGradient) {
            chainSource << "float4 grad1_0_1 = dV0dR1*delta*invR;\n";
            chainSource << "float4 grad1_0_2 = dV0dR2*delta*invR;\n";
            chainSource << "float4 grad2_0_1 = -grad1_0_1;\n";
            chainSource << "float4 grad2_0_2 = -grad1_0_2;\n";
            chainSource << "tempForce1 -= grad1_0_1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "1") << ";\n";
            chainSource << "tempForce1 -= grad1_0_2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "2") << ";\n";
            chainSource << "tempForce2 -= grad2_0_1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "1") << ";\n";
            chainSource << "tempForce2 -= grad2_0_2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "2") << ";\n";
        }
        else {
            chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "1") << ";\n";
            chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "2") << ";\n";
        }
        variables = globalVariables;
        map<string, string> rename1;
        map<string, string> rename2;
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables[name+"1"] = prefix+"params"+params->getParameterSuffix(i, "1");
            variables[name+"2"] = prefix+"params"+params->getParameterSuffix(i, "2");
            rename1[name] = name+"1";
            rename2[name] = name+"2";
        }
        for (int i = 0; i < force.getNumComputedValues(); i++) {
            const string& name = computedValueNames[i];
            variables[name+"1"] = prefix+"values"+computedValues->getParameterSuffix(i, "1");
            variables[name+"2"] = prefix+"values"+computedValues->getParameterSuffix(i, "2");
            rename1[name] = name+"1";
            rename2[name] = name+"2";
            if (i == 0)
                continue;
            string is = intToString(i);
            if (needParameterGradient) {
                chainSource << "float4 grad1_"+is+"_1 = 0;\n";
                chainSource << "float4 grad1_"+is+"_2 = 0;\n";
                chainSource << "float4 grad2_"+is+"_1 = 0;\n";
                chainSource << "float4 grad2_"+is+"_2 = 0;\n";
                for (int j = 0; j < i; j++) {
                    string js = intToString(j);
                    Lepton::ParsedExpression dVdV = Lepton::Parser::parse(computedValueExpressions[i], functions).differentiate(computedValueNames[j]).optimize();
                    derivExpressions.clear();
                    derivExpressions["float dV"+is+"dV"+js+"_1 = "] = dVdV.renameVariables(rename1);
                    derivExpressions["float dV"+is+"dV"+js+"_2 = "] = dVdV.renameVariables(rename2);
                    chainSource << OpenCLExpressionUtilities::createExpressions(derivExpressions, variables, functionDefinitions, prefix+"temp"+is+"_"+js+"_", prefix+"functionParams");
                    chainSource << "grad1_"+is+"_1 += dV"+is+"dV"+js+"_1*grad1_"+js+"_1;\n";
                    chainSource << "grad2_"+is+"_1 += dV"+is+"dV"+js+"_1*grad2_"+js+"_1;\n";
                    chainSource << "grad1_"+is+"_2 += dV"+is+"dV"+js+"_2*grad1_"+js+"_2;\n";
                    chainSource << "grad2_"+is+"_2 += dV"+is+"dV"+js+"_2*grad2_"+js+"_2;\n";
                }
                derivExpressions.clear();
                if (!isZeroExpression(valueGradientExpressions[i][0])) {
                    derivExpressions["grad1_"+is+"_1.x -= "] = valueGradientExpressions[i][0].renameVariables(rename1);
                    derivExpressions["grad2_"+is+"_2.x -= "] = valueGradientExpressions[i][0].renameVariables(rename2);
                }
                if (!isZeroExpression(valueGradientExpressions[i][1])) {
                    derivExpressions["grad1_"+is+"_1.y -= "] = valueGradientExpressions[i][1].renameVariables(rename1);
                    derivExpressions["grad2_"+is+"_2.y -= "] = valueGradientExpressions[i][1].renameVariables(rename2);
                }
                if (!isZeroExpression(valueGradientExpressions[i][2])) {
                    derivExpressions["grad1_"+is+"_1.z -= "] = valueGradientExpressions[i][2].renameVariables(rename1);
                    derivExpressions["grad2_"+is+"_2.z -= "] = valueGradientExpressions[i][2].renameVariables(rename2);
                }
                chainSource << OpenCLExpressionUtilities::createExpressions(derivExpressions, variables, functionDefinitions, prefix+"temp"+is+"_", prefix+"functionParams");
                chainSource << "tempForce1 -= grad1_"<<is<<"_1*"<<prefix<<"dEdV"<<energyDerivs->getParameterSuffix(i, "1")<<";\n";
                chainSource << "tempForce2 -= grad2_"<<is<<"_1*"<<prefix<<"dEdV"<<energyDerivs->getParameterSuffix(i, "1")<<";\n";
                chainSource << "tempForce1 -= grad1_"<<is<<"_2*"<<prefix<<"dEdV"<<energyDerivs->getParameterSuffix(i, "2")<<";\n";
                chainSource << "tempForce2 -= grad2_"<<is<<"_2*"<<prefix<<"dEdV"<<energyDerivs->getParameterSuffix(i, "2")<<";\n";
            }
            else {
                chainSource << "float dV"+is+"dR1 = 0;\n";
                chainSource << "float dV"+is+"dR2 = 0;\n";
                for (int j = 0; j < i; j++) {
                    string js = intToString(j);
                    Lepton::ParsedExpression dVdV = Lepton::Parser::parse(computedValueExpressions[i], functions).differentiate(computedValueNames[j]).optimize();
                    derivExpressions.clear();
                    derivExpressions["dV"+is+"dR1 += dV"+js+"dR1*"] = dVdV.renameVariables(rename1);
                    derivExpressions["dV"+is+"dR2 += dV"+js+"dR2*"] = dVdV.renameVariables(rename2);
                    chainSource << OpenCLExpressionUtilities::createExpressions(derivExpressions, variables, functionDefinitions, prefix+"temp"+is+"_"+js+"_", prefix+"functionParams");
                }
                chainSource << "tempForce -= dV"<< is << "dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "1") << ";\n";
                chainSource << "tempForce -= dV"<< is << "dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "2") << ";\n";
            }
        }
        map<string, string> replacements;
        replacements["COMPUTE_FORCE"] = chainSource.str();
        string source = cl.replaceStrings(OpenCLKernelSources::customGBChainRule, replacements);
        vector<OpenCLNonbondedUtilities::ParameterInfo> parameters;
        vector<OpenCLNonbondedUtilities::ParameterInfo> arguments;
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = prefix+"params"+intToString(i+1);
            parameters.push_back(OpenCLNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string paramName = prefix+"values"+intToString(i+1);
            parameters.push_back(OpenCLNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string paramName = prefix+"dEdV"+intToString(i+1);
            parameters.push_back(OpenCLNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        if (globals != NULL) {
            globals->upload(globalParamValues);
            arguments.push_back(OpenCLNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(cl_float), globals->getDeviceBuffer()));
        }
        if (needParameterGradient) {
            chainRuleParameters = parameters;
            chainRuleArguments = arguments;
            chainRuleSource = source;
            separateChainRuleKernel = true;
            cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, "");
        }
        else {
            cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source);
            for (int i = 0; i < (int) parameters.size(); i++)
                cl.getNonbondedUtilities().addParameter(parameters[i]);
            for (int i = 0; i < (int) arguments.size(); i++)
                cl.getNonbondedUtilities().addArgument(arguments[i]);
            separateChainRuleKernel = false;
        }
    }
    cl.addForce(new OpenCLCustomGBForceInfo(cl.getNonbondedUtilities().getNumForceBuffers(), force));
}

void OpenCLCalcCustomGBForceKernel::executeForces(ContextImpl& context) {
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        valueBuffers = new OpenCLArray<cl_float>(cl, cl.getPaddedNumAtoms()*cl.getNumForceBuffers(), "customGBValueBuffers");
        int index = 0;
        pairValueKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        pairValueKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        pairValueKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusions().getDeviceBuffer());
        pairValueKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusionIndices().getDeviceBuffer());
        pairValueKernel.setArg<cl::Buffer>(index++, valueBuffers->getDeviceBuffer());
        pairValueKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        pairValueKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float), NULL);
        if (nb.getUseCutoff()) {
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
        }
        else {
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getTiles().getDeviceBuffer());
            pairValueKernel.setArg<cl_uint>(index++, nb.getTiles().getSize());
        }
        if (globals != NULL)
            pairValueKernel.setArg<cl::Buffer>(index++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            pairValueKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            pairValueKernel.setArg(index++, OpenCLContext::ThreadBlockSize*buffer.getSize(), NULL);
        }
        if (tabulatedFunctionParams != NULL) {
            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
                pairValueKernel.setArg<cl::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
            pairValueKernel.setArg<cl::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
        }
        index = 0;
        perParticleValueKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        perParticleValueKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
        perParticleValueKernel.setArg<cl::Buffer>(index++, valueBuffers->getDeviceBuffer());
        perParticleValueKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        if (globals != NULL)
            perParticleValueKernel.setArg<cl::Buffer>(index++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++)
            perParticleValueKernel.setArg<cl::Memory>(index++, params->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
            perParticleValueKernel.setArg<cl::Memory>(index++, computedValues->getBuffers()[i].getMemory());
        if (tabulatedFunctionParams != NULL) {
            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
                perParticleValueKernel.setArg<cl::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
            perParticleValueKernel.setArg<cl::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
        }
        index = 0;
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        pairEnergyKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        pairEnergyKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusions().getDeviceBuffer());
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusionIndices().getDeviceBuffer());
        pairEnergyKernel.setArg(index++, OpenCLContext::ThreadBlockSize*sizeof(cl_float4), NULL);
        if (nb.getUseCutoff()) {
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
        }
        else {
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getTiles().getDeviceBuffer());
            pairEnergyKernel.setArg<cl_uint>(index++, nb.getTiles().getSize());
        }
        if (globals != NULL)
            pairEnergyKernel.setArg<cl::Buffer>(index++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            pairEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            pairEnergyKernel.setArg(index++, OpenCLContext::ThreadBlockSize*buffer.getSize(), NULL);
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            pairEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            pairEnergyKernel.setArg(index++, OpenCLContext::ThreadBlockSize*buffer.getSize(), NULL);
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            pairEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            pairEnergyKernel.setArg(index++, OpenCLContext::ThreadBlockSize*buffer.getSize(), NULL);
        }
        if (tabulatedFunctionParams != NULL) {
            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
                pairEnergyKernel.setArg<cl::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
        }
        index = 0;
        perParticleEnergyKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        perParticleEnergyKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
        perParticleEnergyKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        perParticleEnergyKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        perParticleEnergyKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        if (globals != NULL)
            perParticleEnergyKernel.setArg<cl::Buffer>(index++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++)
            perParticleEnergyKernel.setArg<cl::Memory>(index++, params->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
            perParticleEnergyKernel.setArg<cl::Memory>(index++, computedValues->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
            perParticleEnergyKernel.setArg<cl::Memory>(index++, energyDerivs->getBuffers()[i].getMemory());
        if (tabulatedFunctionParams != NULL) {
            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
                perParticleEnergyKernel.setArg<cl::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
            perParticleEnergyKernel.setArg<cl::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
        }
        chainRuleKernel = nb.createInteractionKernel(chainRuleSource, chainRuleParameters, chainRuleArguments, true, false);
    }
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    cl.clearBuffer(*valueBuffers);
    for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
        cl.clearBuffer(buffer.getMemory(), buffer.getSize()*energyDerivs->getNumObjects()/sizeof(cl_float));
    }
    cl.executeKernel(pairValueKernel, nb.getTiles().getSize()*OpenCLContext::TileSize);
    cl.executeKernel(perParticleValueKernel, cl.getPaddedNumAtoms());
    cl.executeKernel(pairEnergyKernel, nb.getTiles().getSize()*OpenCLContext::TileSize);
    cl.executeKernel(perParticleEnergyKernel, cl.getPaddedNumAtoms());
    if (separateChainRuleKernel)
        cl.executeKernel(chainRuleKernel, nb.getTiles().getSize()*OpenCLContext::TileSize);
}

double OpenCLCalcCustomGBForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLCustomExternalForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomExternalForceInfo(const CustomExternalForce& force, int numParticles) : OpenCLForceInfo(1), force(force), indices(numParticles, -1) {
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

OpenCLCalcCustomExternalForceKernel::~OpenCLCalcCustomExternalForceKernel() {
    if (params != NULL)
        delete params;
    if (indices != NULL)
        delete indices;
    if (globals != NULL)
        delete globals;
}

void OpenCLCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    numParticles = force.getNumParticles();
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), numParticles, "customExternalParams");
    indices = new OpenCLArray<cl_int>(cl, numParticles, "customExternalIndices");
    string extraArguments;
    if (force.getNumGlobalParameters() > 0) {
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customExternalGlobals", false, CL_MEM_READ_ONLY);
        extraArguments += ", __constant float* globals";
    }
    vector<vector<cl_float> > paramVector(numParticles);
    vector<cl_int> indicesVector(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, indicesVector[i], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    indices->upload(indicesVector);
    cl.addForce(new OpenCLCustomExternalForceInfo(force, system.getNumParticles()));

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x").optimize();
    Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y").optimize();
    Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z").optimize();
    map<string, Lepton::ParsedExpression> expressions;
    expressions["energy += "] = energyExpression;
    expressions["float dEdX = "] = forceExpressionX;
    expressions["float dEdY = "] = forceExpressionY;
    expressions["float dEdZ = "] = forceExpressionZ;

    // Create the kernels.

    map<string, string> variables;
    variables["x"] = "pos.x";
    variables["y"] = "pos.y";
    variables["z"] = "pos.z";
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        variables[name] = "particleParams"+params->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+intToString(i)+"]";
        variables[name] = value;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArguments += ", __global "+buffer.getType()+"* "+buffer.getName();
        compute<<buffer.getType()<<" particleParams"<<(i+1)<<" = "<<buffer.getName()<<"[index];\n";
    }
    vector<pair<string, string> > functions;
    compute << OpenCLExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["EXTRA_ARGUMENTS"] = extraArguments;
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customExternalForce, replacements));
    kernel = cl::Kernel(program, "computeCustomExternalForces");
}

void OpenCLCalcCustomExternalForceKernel::executeForces(ContextImpl& context) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        kernel.setArg<cl_int>(0, numParticles);
        kernel.setArg<cl::Buffer>(1, cl.getForceBuffers().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(2, cl.getEnergyBuffer().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, indices->getDeviceBuffer());
        int nextIndex = 5;
        if (globals != NULL)
            kernel.setArg<cl::Buffer>(nextIndex++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            kernel.setArg<cl::Memory>(nextIndex++, buffer.getMemory());
        }
    }
    cl.executeKernel(kernel, numParticles);
}

double OpenCLCalcCustomExternalForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

class OpenCLCustomHbondForceInfo : public OpenCLForceInfo {
public:
    OpenCLCustomHbondForceInfo(int requiredBuffers, const CustomHbondForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumDonors()+force.getNumAcceptors()+force.getNumExclusions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int p1, p2, p3;
        vector<double> parameters;
        if (index < force.getNumDonors()) {
            force.getDonorParameters(index, p1, p2, p3, parameters);
            particles.clear();
            particles.push_back(p1);
            if (p2 > -1)
                particles.push_back(p2);
            if (p3 > -1)
                particles.push_back(p3);
            return;
        }
        index -= force.getNumDonors();
        if (index < force.getNumAcceptors()) {
            force.getAcceptorParameters(index, p1, p2, p3, parameters);
            particles.clear();
            particles.push_back(p1);
            if (p2 > -1)
                particles.push_back(p2);
            if (p3 > -1)
                particles.push_back(p3);
            return;
        }
        index -= force.getNumAcceptors();
        int donor, acceptor;
        force.getExclusionParticles(index, donor, acceptor);
        particles.clear();
        force.getDonorParameters(donor, p1, p2, p3, parameters);
        particles.push_back(p1);
        if (p2 > -1)
            particles.push_back(p2);
        if (p3 > -1)
            particles.push_back(p3);
        force.getAcceptorParameters(acceptor, p1, p2, p3, parameters);
        particles.push_back(p1);
        if (p2 > -1)
            particles.push_back(p2);
        if (p3 > -1)
            particles.push_back(p3);
    }
    bool areGroupsIdentical(int group1, int group2) {
        int p1, p2, p3;
        vector<double> params1, params2;
        if (group1 < force.getNumDonors() && group2 < force.getNumDonors()) {
            force.getDonorParameters(group1, p1, p2, p3, params1);
            force.getDonorParameters(group2, p1, p2, p3, params2);
            return (params1 == params2 && params1 == params2);
        }
        if (group1 < force.getNumDonors() || group2 < force.getNumDonors())
            return false;
        group1 -= force.getNumDonors();
        group2 -= force.getNumDonors();
        if (group1 < force.getNumAcceptors() && group2 < force.getNumAcceptors()) {
            force.getAcceptorParameters(group1, p1, p2, p3, params1);
            force.getAcceptorParameters(group2, p1, p2, p3, params2);
            return (params1 == params2 && params1 == params2);
        }
        if (group1 < force.getNumAcceptors() || group2 < force.getNumAcceptors())
            return false;
        return true;
    }
private:
    const CustomHbondForce& force;
};

OpenCLCalcCustomHbondForceKernel::~OpenCLCalcCustomHbondForceKernel() {
    if (donorParams != NULL)
        delete donorParams;
    if (acceptorParams != NULL)
        delete acceptorParams;
    if (donors != NULL)
        delete donors;
    if (acceptors != NULL)
        delete acceptors;
    if (donorBufferIndices != NULL)
        delete donorBufferIndices;
    if (acceptorBufferIndices != NULL)
        delete acceptorBufferIndices;
    if (globals != NULL)
        delete globals;
    if (donorExclusions != NULL)
        delete donorExclusions;
    if (acceptorExclusions != NULL)
        delete acceptorExclusions;
    if (tabulatedFunctionParams != NULL)
        delete tabulatedFunctionParams;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
}

static void addDonorAndAcceptorCode(stringstream& computeDonor, stringstream& computeAcceptor, const string& value) {
    computeDonor << value;
    computeAcceptor << value;
}

static void applyDonorAndAcceptorForces(stringstream& applyToDonor, stringstream& applyToAcceptor, int atom, const string& value) {
    string forceNames[] = {"f1", "f2", "f3"};
    if (atom < 3)
        applyToAcceptor << forceNames[atom]<<".xyz += "<<value<<";\n";
    else
        applyToDonor << forceNames[atom-3]<<".xyz += "<<value<<";\n";
}

void OpenCLCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    // Record the lists of donors and acceptors, and the parameters for each one.

    numDonors = force.getNumDonors();
    numAcceptors = force.getNumAcceptors();
    int numParticles = system.getNumParticles();
    donors = new OpenCLArray<mm_int4>(cl, numDonors, "customHbondDonors");
    acceptors = new OpenCLArray<mm_int4>(cl, numAcceptors, "customHbondAcceptors");
    donorParams = new OpenCLParameterSet(cl, force.getNumPerDonorParameters(), numDonors, "customHbondDonorParameters");
    acceptorParams = new OpenCLParameterSet(cl, force.getNumPerAcceptorParameters(), numAcceptors, "customHbondAcceptorParameters");
    if (force.getNumGlobalParameters() > 0)
        globals = new OpenCLArray<cl_float>(cl, force.getNumGlobalParameters(), "customHbondGlobals", false, CL_MEM_READ_ONLY);
    vector<vector<cl_float> > donorParamVector(numDonors);
    vector<mm_int4> donorVector(numDonors);
    for (int i = 0; i < numDonors; i++) {
        vector<double> parameters;
        force.getDonorParameters(i, donorVector[i].x, donorVector[i].y, donorVector[i].z, parameters);
        donorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            donorParamVector[i][j] = (cl_float) parameters[j];
    }
    donors->upload(donorVector);
    donorParams->setParameterValues(donorParamVector);
    vector<vector<cl_float> > acceptorParamVector(numAcceptors);
    vector<mm_int4> acceptorVector(numAcceptors);
    for (int i = 0; i < numAcceptors; i++) {
        vector<double> parameters;
        force.getAcceptorParameters(i, acceptorVector[i].x, acceptorVector[i].y, acceptorVector[i].z, parameters);
        acceptorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            acceptorParamVector[i][j] = (cl_float) parameters[j];
    }
    acceptors->upload(acceptorVector);
    acceptorParams->setParameterValues(acceptorParamVector);

    // Select an output buffer index for each donor and acceptor.

    donorBufferIndices = new OpenCLArray<mm_int4>(cl, numDonors, "customHbondDonorBuffers");
    acceptorBufferIndices = new OpenCLArray<mm_int4>(cl, numAcceptors, "customHbondAcceptorBuffers");
    vector<mm_int4> donorBufferVector(numDonors);
    vector<mm_int4> acceptorBufferVector(numAcceptors);
    vector<int> donorBufferCounter(numParticles, 0);
    for (int i = 0; i < numDonors; i++)
        donorBufferVector[i] = mm_int4(donorVector[i].x > -1 ? donorBufferCounter[donorVector[i].x]++ : 0,
                                       donorVector[i].y > -1 ? donorBufferCounter[donorVector[i].y]++ : 0,
                                       donorVector[i].z > -1 ? donorBufferCounter[donorVector[i].z]++ : 0, 0);
    vector<int> acceptorBufferCounter(numParticles, 0);
    for (int i = 0; i < numAcceptors; i++)
        acceptorBufferVector[i] = mm_int4(acceptorVector[i].x > -1 ? acceptorBufferCounter[acceptorVector[i].x]++ : 0,
                                       acceptorVector[i].y > -1 ? acceptorBufferCounter[acceptorVector[i].y]++ : 0,
                                       acceptorVector[i].z > -1 ? acceptorBufferCounter[acceptorVector[i].z]++ : 0, 0);
    donorBufferIndices->upload(donorBufferVector);
    acceptorBufferIndices->upload(acceptorBufferVector);
    int maxBuffers = 1;
    for (int i = 0; i < (int) donorBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, donorBufferCounter[i]);
    for (int i = 0; i < (int) acceptorBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, acceptorBufferCounter[i]);
    cl.addForce(new OpenCLCustomHbondForceInfo(maxBuffers, force));

    // Record exclusions.

    vector<mm_int4> donorExclusionVector(numDonors, mm_int4(-1, -1, -1, -1));
    vector<mm_int4> acceptorExclusionVector(numAcceptors, mm_int4(-1, -1, -1, -1));
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int donor, acceptor;
        force.getExclusionParticles(i, donor, acceptor);
        if (donorExclusionVector[donor].x == -1)
            donorExclusionVector[donor].x = acceptor;
        else if (donorExclusionVector[donor].y == -1)
            donorExclusionVector[donor].y = acceptor;
        else if (donorExclusionVector[donor].z == -1)
            donorExclusionVector[donor].z = acceptor;
        else if (donorExclusionVector[donor].w == -1)
            donorExclusionVector[donor].w = acceptor;
        else
            throw OpenMMException("CustomHbondForce: OpenCLPlatform does not support more than four exclusions per donor");
        if (acceptorExclusionVector[acceptor].x == -1)
            acceptorExclusionVector[acceptor].x = donor;
        else if (acceptorExclusionVector[acceptor].y == -1)
            acceptorExclusionVector[acceptor].y = donor;
        else if (acceptorExclusionVector[acceptor].z == -1)
            acceptorExclusionVector[acceptor].z = donor;
        else if (acceptorExclusionVector[acceptor].w == -1)
            acceptorExclusionVector[acceptor].w = donor;
        else
            throw OpenMMException("CustomHbondForce: OpenCLPlatform does not support more than four exclusions per acceptor");
    }
    donorExclusions = new OpenCLArray<mm_int4>(cl, numDonors, "customHbondDonorExclusions");
    acceptorExclusions = new OpenCLArray<mm_int4>(cl, numDonors, "customHbondAcceptorExclusions");
    donorExclusions->upload(donorExclusionVector);
    acceptorExclusions->upload(acceptorExclusionVector);

    // Record the tabulated functions.

    OpenCLExpressionUtilities::FunctionPlaceholder fp;
    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
    stringstream tableArgs;
    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        bool interpolating;
        force.getFunctionParameters(i, name, values, min, max, interpolating);
        string arrayName = "table"+intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = &fp;
        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), 0.0f);
        vector<mm_float4> f = OpenCLExpressionUtilities::computeFunctionCoefficients(values, interpolating);
        tabulatedFunctions.push_back(new OpenCLArray<mm_float4>(cl, values.size()-1, "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        tableArgs << ", __global float4* " << arrayName;
    }
    if (force.getNumFunctions() > 0) {
        tabulatedFunctionParams = new OpenCLArray<mm_float4>(cl, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
        tableArgs << ", __constant float4* functionParams";
    }

    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);
    map<string, string> variables;
    for (int i = 0; i < force.getNumPerDonorParameters(); i++) {
        const string& name = force.getPerDonorParameterName(i);
        variables[name] = "donorParams"+donorParams->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumPerAcceptorParameters(); i++) {
        const string& name = force.getPerAcceptorParameterName(i);
        variables[name] = "acceptorParams"+acceptorParams->getParameterSuffix(i);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        variables[name] = "globals["+intToString(i)+"]";
    }

    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
    // and dihedrals the expression depends on.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomHbondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    map<string, Lepton::ParsedExpression> forceExpressions;
    set<string> computedDeltas;
    computedDeltas.insert("D1A1");
    string atomNames[] = {"A1", "A2", "A3", "D1", "D2", "D3"};
    string atomNamesLower[] = {"a1", "a2", "a3", "d1", "d2", "d3"};
    stringstream computeDonor, computeAcceptor, extraArgs;
    int index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+");\n");
            computedDeltas.insert(deltaName);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float r_"+deltaName+" = sqrt(delta"+deltaName+".w);\n");
        variables[iter->first] = "r_"+deltaName;
        forceExpressions["float dEdDistance"+intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[0]]+");\n");
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[2]]+");\n");
            computedDeltas.insert(deltaName2);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float "+angleName+" = computeAngle(delta"+deltaName1+", delta"+deltaName2+");\n");
        variables[iter->first] = angleName;
        forceExpressions["float dEdAngle"+intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        string dihedralName = "dihedral_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]]+atomNames[atoms[3]];
        if (computedDeltas.count(deltaName1) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+");\n");
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[1]]+");\n");
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName3+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[3]]+");\n");
            computedDeltas.insert(deltaName3);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 "+crossName1+" = computeCross(delta"+deltaName1+", delta"+deltaName2+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 "+crossName2+" = computeCross(delta"+deltaName2+", delta"+deltaName3+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float "+dihedralName+" = computeAngle("+crossName1+", "+crossName2+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, dihedralName+" *= (delta"+deltaName1+".x*"+crossName2+".x + delta"+deltaName1+".y*"+crossName2+".y + delta"+deltaName1+".z*"+crossName2+".z < 0 ? -1 : 1);\n");
        variables[iter->first] = dihedralName;
        forceExpressions["float dEdDihedral"+intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }

    // Next it needs to load parameters from global memory.

    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", __constant float* globals";
    for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
        extraArgs << ", __global "+buffer.getType()+"* donor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" donorParams"+intToString(i+1)+" = donor"+buffer.getName()+"[index];\n");
    }
    for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
        extraArgs << ", __global "+buffer.getType()+"* acceptor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" acceptorParams"+intToString(i+1)+" = acceptor"+buffer.getName()+"[index];\n");
    }

    // Now evaluate the expressions.

    computeAcceptor << OpenCLExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, "temp", "functionParams");
    forceExpressions["energy += "] = energyExpression;
    computeDonor << OpenCLExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, "temp", "functionParams");

    // Finally, apply forces to atoms.

    index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "-"+value);
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], value);
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 crossProd = cross(delta"+deltaName2+", delta"+deltaName1+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float lengthCross = max(length(crossProd), 1e-6f);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 deltaCross0 = -cross(delta"+deltaName1+", crossProd)*dEdAngle"+intToString(index)+"/(delta"+deltaName1+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 deltaCross2 = cross(delta"+deltaName2+", crossProd)*dEdAngle"+intToString(index)+"/(delta"+deltaName2+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 deltaCross1 = -(deltaCross0+deltaCross2);\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "deltaCross0.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "deltaCross1.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "deltaCross2.xyz");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float r = sqrt(delta"+deltaName2+".w);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 ff;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.x = (-dEdDihedral"+intToString(index)+"*r)/"+crossName1+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.y = (delta"+deltaName1+".x*delta"+deltaName2+".x + delta"+deltaName1+".y*delta"+deltaName2+".y + delta"+deltaName1+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.z = (delta"+deltaName3+".x*delta"+deltaName2+".x + delta"+deltaName3+".y*delta"+deltaName2+".y + delta"+deltaName3+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.w = (dEdDihedral"+intToString(index)+"*r)/"+crossName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 internalF0 = ff.x*"+crossName1+";\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 internalF3 = ff.w*"+crossName2+";\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 s = ff.y*internalF0 - ff.z*internalF3;\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "internalF0.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "s.xyz-internalF0.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "-s.xyz-internalF3.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[3], "internalF3.xyz");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
    }

    // Generate the kernels.

    map<string, string> replacements;
    replacements["COMPUTE_DONOR_FORCE"] = computeDonor.str();
    replacements["COMPUTE_ACCEPTOR_FORCE"] = computeAcceptor.str();
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
    defines["NUM_DONORS"] = intToString(force.getNumDonors());
    defines["NUM_ACCEPTORS"] = intToString(force.getNumAcceptors());
    defines["M_PI"] = doubleToString(M_PI);
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff) {
        defines["USE_CUTOFF"] = "1";
        defines["CUTOFF_SQUARED"] = doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff && force.getNonbondedMethod() != CustomHbondForce::CutoffNonPeriodic)
        defines["USE_PERIODIC"] = "1";
    if (force.getNumExclusions() > 0)
        defines["USE_EXCLUSIONS"] = "1";
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customHbondForce, replacements), defines);
    donorKernel = cl::Kernel(program, "computeDonorForces");
    acceptorKernel = cl::Kernel(program, "computeAcceptorForces");
}

void OpenCLCalcCustomHbondForceKernel::executeForces(ContextImpl& context) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        int index = 0;
        donorKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, donorExclusions->getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, donors->getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, acceptors->getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, donorBufferIndices->getDeviceBuffer());
        donorKernel.setArg(index++, 3*OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
        if (globals != NULL)
            donorKernel.setArg<cl::Buffer>(index++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
            donorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        }
        for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
            donorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        }
        if (tabulatedFunctionParams != NULL) {
            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
                donorKernel.setArg<cl::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
            donorKernel.setArg<cl::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
        }
        index = 0;
        acceptorKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, acceptorExclusions->getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, donors->getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, acceptors->getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, acceptorBufferIndices->getDeviceBuffer());
        acceptorKernel.setArg(index++, 3*OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
        if (globals != NULL)
            acceptorKernel.setArg<cl::Buffer>(index++, globals->getDeviceBuffer());
        for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
            acceptorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        }
        for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
            acceptorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        }
        if (tabulatedFunctionParams != NULL) {
            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
                acceptorKernel.setArg<cl::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
            acceptorKernel.setArg<cl::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
        }
    }
    cl.executeKernel(donorKernel, std::max(numDonors, numAcceptors));
    cl.executeKernel(acceptorKernel, std::max(numDonors, numAcceptors));
}

double OpenCLCalcCustomHbondForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
    return 0.0;
}

OpenCLIntegrateVerletStepKernel::~OpenCLIntegrateVerletStepKernel() {
}

void OpenCLIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cl.initialize(system);
    cl::Program program = cl.createProgram(OpenCLKernelSources::verlet);
    kernel1 = cl::Kernel(program, "integrateVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVerletPart2");
    prevStepSize = -1.0;
}

void OpenCLIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    double dt = integrator.getStepSize();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl_int>(0, numAtoms);
        kernel1.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl_int>(0, numAtoms);
        kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getPosDelta().getDeviceBuffer());
    }
    if (dt != prevStepSize) {
        vector<mm_float2> stepSizeVec(1);
        stepSizeVec[0] = mm_float2((cl_float) dt, (cl_float) dt);
        cl.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        prevStepSize = dt;
    }

    // Call the first integration kernel.

    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
}

OpenCLIntegrateLangevinStepKernel::~OpenCLIntegrateLangevinStepKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    cl.initialize(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::langevin, defines);
    kernel1 = cl::Kernel(program, "integrateLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateLangevinPart2");
    params = new OpenCLArray<cl_float>(cl, 3, "langevinParams");
    prevStepSize = -1.0;
}

void OpenCLIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl::Buffer>(0, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, params->getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(1, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, integration.getStepSize().getDeviceBuffer());
    }
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        double kT = BOLTZ*temperature;
        double vscale = std::exp(-stepSize/tau);
        double fscale = (1-vscale)*tau;
        double noisescale = std::sqrt(2*kT/tau)*std::sqrt(0.5*(1-vscale*vscale)*tau);
        vector<cl_float> p(params->getSize());
        p[0] = (cl_float) vscale;
        p[1] = (cl_float) fscale;
        p[2] = (cl_float) noisescale;
        params->upload(p);
        integration.getStepSize()[0].y = stepSize;
        integration.getStepSize().upload();
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Call the first integration kernel.

    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
}

OpenCLIntegrateBrownianStepKernel::~OpenCLIntegrateBrownianStepKernel() {
}

void OpenCLIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    cl.initialize(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::brownian, defines);
    kernel1 = cl::Kernel(program, "integrateBrownianPart1");
    kernel2 = cl::Kernel(program, "integrateBrownianPart2");
    prevStepSize = -1.0;
}

void OpenCLIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl::Buffer>(2, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, integration.getPosDelta().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, integration.getPosDelta().getDeviceBuffer());
    }
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        kernel1.setArg<cl_float>(0, (cl_float) (tau*stepSize));
        kernel1.setArg<cl_float>(1, (cl_float) (sqrt(2.0f*BOLTZ*temperature*stepSize*tau)));
        kernel2.setArg<cl_float>(0, (cl_float) (1.0/stepSize));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Call the first integration kernel.

    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
}

OpenCLIntegrateVariableVerletStepKernel::~OpenCLIntegrateVariableVerletStepKernel() {
}

void OpenCLIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    cl.initialize(system);
    cl::Program program = cl.createProgram(OpenCLKernelSources::verlet);
    kernel1 = cl::Kernel(program, "integrateVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVerletPart2");
    selectSizeKernel = cl::Kernel(program, "selectVerletStepSize");
    blockSize = std::min(std::min(256, system.getNumParticles()), (int) cl.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
}

void OpenCLIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl_int>(0, numAtoms);
        kernel1.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl_int>(0, numAtoms);
        kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getPosDelta().getDeviceBuffer());
        selectSizeKernel.setArg<cl_int>(0, numAtoms);
        selectSizeKernel.setArg<cl::Buffer>(3, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(5, cl.getForce().getDeviceBuffer());
        selectSizeKernel.setArg(6, blockSize*sizeof(cl_float), NULL);
    }

    // Select the step size to use.

    float maxStepSize = (float)(maxTime-cl.getTime());
    selectSizeKernel.setArg<cl_float>(1, maxStepSize);
    selectSizeKernel.setArg<cl_float>(2, (cl_float) integrator.getErrorTolerance());
    cl.executeKernel(selectSizeKernel, blockSize, blockSize);

    // Call the first integration kernel.

    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);

    // Update the time and step count.

    cl.getIntegrationUtilities().getStepSize().download();
    double dt = cl.getIntegrationUtilities().getStepSize()[0].y;
    double time = cl.getTime()+dt;
    if (dt == maxStepSize)
        time = maxTime; // Avoid round-off error
    cl.setTime(time);
    cl.setStepCount(cl.getStepCount()+1);
}

OpenCLIntegrateVariableLangevinStepKernel::~OpenCLIntegrateVariableLangevinStepKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    cl.initialize(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::langevin, defines);
    kernel1 = cl::Kernel(program, "integrateLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateLangevinPart2");
    selectSizeKernel = cl::Kernel(program, "selectLangevinStepSize");
    params = new OpenCLArray<cl_float>(cl, 3, "langevinParams");
    blockSize = std::min(256, system.getNumParticles());
    blockSize = std::max(blockSize, params->getSize());
    blockSize = std::min(blockSize, (int) cl.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
}

void OpenCLIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl::Buffer>(0, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, params->getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(1, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, integration.getStepSize().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(5, cl.getVelm().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(6, cl.getForce().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(7, params->getDeviceBuffer());
        selectSizeKernel.setArg(8, params->getSize()*sizeof(cl_float), NULL);
        selectSizeKernel.setArg(9, blockSize*sizeof(cl_float), NULL);
    }

    // Select the step size to use.

    float maxStepSize = (float)(maxTime-cl.getTime());
    selectSizeKernel.setArg<cl_float>(0, maxStepSize);
    selectSizeKernel.setArg<cl_float>(1, (cl_float) integrator.getErrorTolerance());
    selectSizeKernel.setArg<cl_float>(2, (cl_float) (integrator.getFriction() == 0.0 ? 0.0 : 1.0/integrator.getFriction()));
    selectSizeKernel.setArg<cl_float>(3, (cl_float) (BOLTZ*integrator.getTemperature()));
    cl.executeKernel(selectSizeKernel, blockSize, blockSize);

    // Call the first integration kernel.

    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);

    // Update the time and step count.

    cl.getIntegrationUtilities().getStepSize().download();
    double dt = cl.getIntegrationUtilities().getStepSize()[0].y;
    double time = cl.getTime()+dt;
    if (dt == maxStepSize)
        time = maxTime; // Avoid round-off error
    cl.setTime(time);
    cl.setStepCount(cl.getStepCount()+1);
}

OpenCLApplyAndersenThermostatKernel::~OpenCLApplyAndersenThermostatKernel() {
}

void OpenCLApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    randomSeed = thermostat.getRandomNumberSeed();
    map<string, string> defines;
    defines["NUM_ATOMS"] = intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::andersenThermostat, defines);
    kernel = cl::Kernel(program, "applyAndersenThermostat");
}

void OpenCLApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        cl.getIntegrationUtilities().initRandomNumberGenerator(randomSeed);
        kernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getIntegrationUtilities().getRandom().getDeviceBuffer());
    }
    kernel.setArg<cl_float>(0, (cl_float) context.getParameter(AndersenThermostat::CollisionFrequency()));
    kernel.setArg<cl_float>(1, (cl_float) (BOLTZ*context.getParameter(AndersenThermostat::Temperature())));
    kernel.setArg<cl_uint>(5, cl.getIntegrationUtilities().prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel, cl.getNumAtoms());
}

void OpenCLCalcKineticEnergyKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

double OpenCLCalcKineticEnergyKernel::execute(ContextImpl& context) {
    // We don't currently have a GPU kernel to do this, so we retrieve the velocities and calculate the energy
    // on the CPU.

    OpenCLArray<mm_float4>& velm = cl.getVelm();
    velm.download();
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i) {
        mm_float4 v = velm[i];
        energy += masses[i]*(v.x*v.x+v.y*v.y+v.z*v.z);
    }
    return 0.5*energy;
}

OpenCLRemoveCMMotionKernel::~OpenCLRemoveCMMotionKernel() {
    if (cmMomentum != NULL)
        delete cmMomentum;
}

void OpenCLRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    frequency = force.getFrequency();
    int numAtoms = cl.getNumAtoms();
    cmMomentum = new OpenCLArray<mm_float4>(cl, (numAtoms+OpenCLContext::ThreadBlockSize-1)/OpenCLContext::ThreadBlockSize, "cmMomentum");
    double totalMass = 0.0;
    for (int i = 0; i < numAtoms; i++)
        totalMass += system.getParticleMass(i);
    map<string, string> defines;
    defines["INVERSE_TOTAL_MASS"] = doubleToString(1.0/totalMass);
    cl::Program program = cl.createProgram(OpenCLKernelSources::removeCM, defines);
    kernel1 = cl::Kernel(program, "calcCenterOfMassMomentum");
    kernel1.setArg<cl_int>(0, numAtoms);
    kernel1.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(2, cmMomentum->getDeviceBuffer());
    kernel1.setArg(3, OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
    kernel2 = cl::Kernel(program, "removeCenterOfMassMomentum");
    kernel2.setArg<cl_int>(0, numAtoms);
    kernel2.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(2, cmMomentum->getDeviceBuffer());
    kernel2.setArg(3, OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
}

void OpenCLRemoveCMMotionKernel::execute(ContextImpl& context) {
    cl.executeKernel(kernel1, cl.getNumAtoms());
    cl.executeKernel(kernel2, cl.getNumAtoms());
}
