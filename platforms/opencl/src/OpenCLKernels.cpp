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

#include "OpenCLKernels.h"
#include "OpenCLForceInfo.h"
#include "openmm/Context.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/internal/CustomHbondForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "OpenCLBondedUtilities.h"
#include "OpenCLExpressionUtilities.h"
#include "OpenCLIntegrationUtilities.h"
#include "OpenCLNonbondedUtilities.h"
#include "OpenCLKernelSources.h"
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "ReferenceTabulatedFunction.h"
#include "SimTKOpenMMRealType.h"
#include "SimTKOpenMMUtilities.h"
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iterator>
#include <set>

using namespace OpenMM;
using namespace std;
using namespace Lepton;

static void setPosqCorrectionArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseMixedPrecision())
        kernel.setArg<cl::Buffer>(index, cl.getPosqCorrection().getDeviceBuffer());
    else
        kernel.setArg<void*>(index, NULL);
}

static void setPeriodicBoxSizeArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseDoublePrecision())
        kernel.setArg<mm_double4>(index, cl.getPeriodicBoxSizeDouble());
    else
        kernel.setArg<mm_float4>(index, cl.getPeriodicBoxSize());
}

static void setPeriodicBoxArgs(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseDoublePrecision()) {
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(index++, cl.getInvPeriodicBoxSizeDouble());
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxVecXDouble());
        kernel.setArg<mm_double4>(index++, cl.getPeriodicBoxVecYDouble());
        kernel.setArg<mm_double4>(index, cl.getPeriodicBoxVecZDouble());
    }
    else {
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxSize());
        kernel.setArg<mm_float4>(index++, cl.getInvPeriodicBoxSize());
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxVecX());
        kernel.setArg<mm_float4>(index++, cl.getPeriodicBoxVecY());
        kernel.setArg<mm_float4>(index, cl.getPeriodicBoxVecZ());
    }
}

static bool isZeroExpression(const Lepton::ParsedExpression& expression) {
    const Lepton::Operation& op = expression.getRootNode().getOperation();
    if (op.getId() != Lepton::Operation::CONSTANT)
        return false;
    return (dynamic_cast<const Lepton::Operation::Constant&>(op).getValue() == 0.0);
}

static bool usesVariable(const Lepton::ExpressionTreeNode& node, const string& variable) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::VARIABLE && op.getName() == variable)
        return true;
    for (auto& child : node.getChildren())
        if (usesVariable(child, variable))
            return true;
    return false;
}

static bool usesVariable(const Lepton::ParsedExpression& expression, const string& variable) {
    return usesVariable(expression.getRootNode(), variable);
}

static pair<ExpressionTreeNode, string> makeVariable(const string& name, const string& value) {
    return make_pair(ExpressionTreeNode(new Operation::Variable(name)), value);
}

static void replaceFunctionsInExpression(map<string, CustomFunction*>& functions, ExpressionProgram& expression) {
    for (int i = 0; i < expression.getNumOperations(); i++) {
        if (expression.getOperation(i).getId() == Operation::CUSTOM) {
            const Operation::Custom& op = dynamic_cast<const Operation::Custom&>(expression.getOperation(i));
            expression.setOperation(i, new Operation::Custom(op.getName(), functions[op.getName()]->clone(), op.getDerivOrder()));
        }
    }
}

void OpenCLCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void OpenCLCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cl.setForcesValid(true);
    cl.clearAutoclearBuffers();
    for (auto computation : cl.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    cl.setComputeForceCount(cl.getComputeForceCount()+1);
    nb.prepareInteractions(groups);
    map<string, double>& derivs = cl.getEnergyParamDerivWorkspace();
    for (auto& param : context.getParameters())
        derivs[param.first] = 0;
}

double OpenCLCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    cl.getBondedUtilities().computeInteractions(groups);
    cl.getNonbondedUtilities().computeInteractions(groups, includeForces, includeEnergy);
    double sum = 0.0;
    for (auto computation : cl.getPostComputations())
        sum += computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    cl.reduceForces();
    cl.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy)
        sum += cl.reduceEnergy();
    if (!cl.getForcesValid())
        valid = false;
    return sum;
}

void OpenCLUpdateStateDataKernel::initialize(const System& system) {
}

double OpenCLUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return cl.getTime();
}

void OpenCLUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    vector<OpenCLContext*>& contexts = cl.getPlatformData().contexts;
    for (auto ctx : contexts)
        ctx->setTime(time);
}

void OpenCLUpdateStateDataKernel::getPositions(ContextImpl& context, vector<Vec3>& positions) {
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    vector<mm_float4> posCorrection;
    if (cl.getUseDoublePrecision()) {
        mm_double4* posq = (mm_double4*) cl.getPinnedBuffer();
        cl.getPosq().download(posq);
    }
    else if (cl.getUseMixedPrecision()) {
        mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
        cl.getPosq().download(posq, false);
        posCorrection.resize(numParticles);
        cl.getPosqCorrection().download(posCorrection);
    }
    else {
        mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
        cl.getPosq().download(posq);
    }
    
    // Filling in the output array is done in parallel for speed.
    
    cl.getPlatformData().threads.execute([&] (ThreadPool& threads, int threadIndex) {
        // Compute the position of each particle to return to the user.  This is done in parallel for speed.
        
        const vector<int>& order = cl.getAtomIndex();
        int numParticles = cl.getNumAtoms();
        Vec3 boxVectors[3];
        cl.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (cl.getUseDoublePrecision()) {
            mm_double4* posq = (mm_double4*) cl.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                mm_double4 pos = posq[i];
                mm_int4 offset = cl.getPosCellOffsets()[i];
                positions[order[i]] = Vec3(pos.x, pos.y, pos.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
        else if (cl.getUseMixedPrecision()) {
            mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                mm_float4 pos1 = posq[i];
                mm_float4 pos2 = posCorrection[i];
                mm_int4 offset = cl.getPosCellOffsets()[i];
                positions[order[i]] = Vec3((double)pos1.x+(double)pos2.x, (double)pos1.y+(double)pos2.y, (double)pos1.z+(double)pos2.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
        else {
            mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                mm_float4 pos = posq[i];
                mm_int4 offset = cl.getPosCellOffsets()[i];
                positions[order[i]] = Vec3(pos.x, pos.y, pos.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
    });
    cl.getPlatformData().threads.waitForThreads();
}

void OpenCLUpdateStateDataKernel::setPositions(ContextImpl& context, const vector<Vec3>& positions) {
    const vector<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cl.getUseDoublePrecision()) {
        mm_double4* posq = (mm_double4*) cl.getPinnedBuffer();
        cl.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4& pos = posq[i];
            const Vec3& p = positions[order[i]];
            pos.x = p[0];
            pos.y = p[1];
            pos.z = p[2];
        }
        for (int i = numParticles; i < cl.getPaddedNumAtoms(); i++)
            posq[i] = mm_double4(0.0, 0.0, 0.0, 0.0);
        cl.getPosq().upload(posq);
    }
    else {
        mm_float4* posq = (mm_float4*) cl.getPinnedBuffer();
        cl.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4& pos = posq[i];
            const Vec3& p = positions[order[i]];
            pos.x = (cl_float) p[0];
            pos.y = (cl_float) p[1];
            pos.z = (cl_float) p[2];
        }
        for (int i = numParticles; i < cl.getPaddedNumAtoms(); i++)
            posq[i] = mm_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cl.getPosq().upload(posq);
    }
    if (cl.getUseMixedPrecision()) {
        mm_float4* posCorrection = (mm_float4*) cl.getPinnedBuffer();
        for (int i = 0; i < numParticles; ++i) {
            mm_float4& c = posCorrection[i];
            const Vec3& p = positions[order[i]];
            c.x = (cl_float) (p[0]-(cl_float)p[0]);
            c.y = (cl_float) (p[1]-(cl_float)p[1]);
            c.z = (cl_float) (p[2]-(cl_float)p[2]);
            c.w = 0;
        }
        for (int i = numParticles; i < cl.getPaddedNumAtoms(); i++)
            posCorrection[i] = mm_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cl.getPosqCorrection().upload(posCorrection);
    }
    for (auto& offset : cl.getPosCellOffsets())
        offset = mm_int4(0, 0, 0, 0);
    cl.reorderAtoms();
}

void OpenCLUpdateStateDataKernel::getVelocities(ContextImpl& context, vector<Vec3>& velocities) {
    const vector<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        mm_double4* velm = (mm_double4*) cl.getPinnedBuffer();
        cl.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4 vel = velm[i];
            mm_int4 offset = cl.getPosCellOffsets()[i];
            velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
        }
    }
    else {
        mm_float4* velm = (mm_float4*) cl.getPinnedBuffer();
        cl.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4 vel = velm[i];
            mm_int4 offset = cl.getPosCellOffsets()[i];
            velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
        }
    }
}

void OpenCLUpdateStateDataKernel::setVelocities(ContextImpl& context, const vector<Vec3>& velocities) {
    const vector<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        mm_double4* velm = (mm_double4*) cl.getPinnedBuffer();
        cl.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4& vel = velm[i];
            const Vec3& p = velocities[order[i]];
            vel.x = p[0];
            vel.y = p[1];
            vel.z = p[2];
        }
        for (int i = numParticles; i < cl.getPaddedNumAtoms(); i++)
            velm[i] = mm_double4(0.0, 0.0, 0.0, 0.0);
        cl.getVelm().upload(velm);
    }
    else {
        mm_float4* velm = (mm_float4*) cl.getPinnedBuffer();
        cl.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4& vel = velm[i];
            const Vec3& p = velocities[order[i]];
            vel.x = p[0];
            vel.y = p[1];
            vel.z = p[2];
        }
        for (int i = numParticles; i < cl.getPaddedNumAtoms(); i++)
            velm[i] = mm_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cl.getVelm().upload(velm);
    }
}

void OpenCLUpdateStateDataKernel::getForces(ContextImpl& context, vector<Vec3>& forces) {
    const vector<cl_int>& order = cl.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    forces.resize(numParticles);
    if (cl.getUseDoublePrecision()) {
        mm_double4* force = (mm_double4*) cl.getPinnedBuffer();
        cl.getForce().download(force);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4 f = force[i];
            forces[order[i]] = Vec3(f.x, f.y, f.z);
        }
    }
    else {
        mm_float4* force = (mm_float4*) cl.getPinnedBuffer();
        cl.getForce().download(force);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4 f = force[i];
            forces[order[i]] = Vec3(f.x, f.y, f.z);
        }
    }
}

void OpenCLUpdateStateDataKernel::getEnergyParameterDerivatives(ContextImpl& context, map<string, double>& derivs) {
    const vector<string>& paramDerivNames = cl.getEnergyParamDerivNames();
    int numDerivs = paramDerivNames.size();
    if (numDerivs == 0)
        return;
    derivs = cl.getEnergyParamDerivWorkspace();
    OpenCLArray& derivArray = cl.getEnergyParamDerivBuffer();
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        vector<double> derivBuffers;
        derivArray.download(derivBuffers);
        for (int i = numDerivs; i < derivArray.getSize(); i += numDerivs)
            for (int j = 0; j < numDerivs; j++)
                derivBuffers[j] += derivBuffers[i+j];
        for (int i = 0; i < numDerivs; i++)
            derivs[paramDerivNames[i]] += derivBuffers[i];
    }
    else {
        vector<float> derivBuffers;
        derivArray.download(derivBuffers);
        for (int i = numDerivs; i < derivArray.getSize(); i += numDerivs)
            for (int j = 0; j < numDerivs; j++)
                derivBuffers[j] += derivBuffers[i+j];
        for (int i = 0; i < numDerivs; i++)
            derivs[paramDerivNames[i]] += derivBuffers[i];
    }
}

void OpenCLUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    cl.getPeriodicBoxVectors(a, b, c);
}

void OpenCLUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) {
    vector<OpenCLContext*>& contexts = cl.getPlatformData().contexts;

    // If any particles have been wrapped to the first periodic box, we need to unwrap them
    // to avoid changing their positions.

    vector<Vec3> positions;
    for (auto offset : cl.getPosCellOffsets()) {
        if (offset.x != 0 || offset.y != 0 || offset.z != 0) {
            getPositions(context, positions);
            break;
        }
    }
    
    // Update the vectors.

    for (auto ctx : contexts)
        ctx->setPeriodicBoxVectors(a, b, c);
    if (positions.size() > 0)
        setPositions(context, positions);
}

void OpenCLUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    int version = 3;
    stream.write((char*) &version, sizeof(int));
    int precision = (cl.getUseDoublePrecision() ? 2 : cl.getUseMixedPrecision() ? 1 : 0);
    stream.write((char*) &precision, sizeof(int));
    double time = cl.getTime();
    stream.write((char*) &time, sizeof(double));
    int stepCount = cl.getStepCount();
    stream.write((char*) &stepCount, sizeof(int));
    int stepsSinceReorder = cl.getStepsSinceReorder();
    stream.write((char*) &stepsSinceReorder, sizeof(int));
    char* buffer = (char*) cl.getPinnedBuffer();
    cl.getPosq().download(buffer);
    stream.write(buffer, cl.getPosq().getSize()*cl.getPosq().getElementSize());
    if (cl.getUseMixedPrecision()) {
        cl.getPosqCorrection().download(buffer);
        stream.write(buffer, cl.getPosqCorrection().getSize()*cl.getPosqCorrection().getElementSize());
    }
    cl.getVelm().download(buffer);
    stream.write(buffer, cl.getVelm().getSize()*cl.getVelm().getElementSize());
    stream.write((char*) &cl.getAtomIndex()[0], sizeof(cl_int)*cl.getAtomIndex().size());
    stream.write((char*) &cl.getPosCellOffsets()[0], sizeof(mm_int4)*cl.getPosCellOffsets().size());
    Vec3 boxVectors[3];
    cl.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    stream.write((char*) boxVectors, 3*sizeof(Vec3));
    cl.getIntegrationUtilities().createCheckpoint(stream);
    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void OpenCLUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    int version;
    stream.read((char*) &version, sizeof(int));
    if (version != 3)
        throw OpenMMException("Checkpoint was created with a different version of OpenMM");
    int precision;
    stream.read((char*) &precision, sizeof(int));
    int expectedPrecision = (cl.getUseDoublePrecision() ? 2 : cl.getUseMixedPrecision() ? 1 : 0);
    if (precision != expectedPrecision)
        throw OpenMMException("Checkpoint was created with a different numeric precision");
    double time;
    stream.read((char*) &time, sizeof(double));
    int stepCount, stepsSinceReorder;
    stream.read((char*) &stepCount, sizeof(int));
    stream.read((char*) &stepsSinceReorder, sizeof(int));
    vector<OpenCLContext*>& contexts = cl.getPlatformData().contexts;
    for (auto ctx : contexts) {
        ctx->setTime(time);
        ctx->setStepCount(stepCount);
        ctx->setStepsSinceReorder(stepsSinceReorder);
    }
    char* buffer = (char*) cl.getPinnedBuffer();
    stream.read(buffer, cl.getPosq().getSize()*cl.getPosq().getElementSize());
    cl.getPosq().upload(buffer);
    if (cl.getUseMixedPrecision()) {
        stream.read(buffer, cl.getPosqCorrection().getSize()*cl.getPosqCorrection().getElementSize());
        cl.getPosqCorrection().upload(buffer);
    }
    stream.read(buffer, cl.getVelm().getSize()*cl.getVelm().getElementSize());
    cl.getVelm().upload(buffer);
    stream.read((char*) &cl.getAtomIndex()[0], sizeof(cl_int)*cl.getAtomIndex().size());
    cl.getAtomIndexArray().upload(cl.getAtomIndex());
    stream.read((char*) &cl.getPosCellOffsets()[0], sizeof(mm_int4)*cl.getPosCellOffsets().size());
    Vec3 boxVectors[3];
    stream.read((char*) &boxVectors, 3*sizeof(Vec3));
    for (auto ctx : contexts)
        ctx->setPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    cl.getIntegrationUtilities().loadCheckpoint(stream);
    SimTKOpenMMUtilities::loadCheckpoint(stream);
    for (auto listener : cl.getReorderListeners())
        listener->execute();
}

void OpenCLApplyConstraintsKernel::initialize(const System& system) {
}

void OpenCLApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        map<string, string> defines;
        defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        cl::Program program = cl.createProgram(OpenCLKernelSources::constraints, defines);
        applyDeltasKernel = cl::Kernel(program, "applyPositionDeltas");
        applyDeltasKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, applyDeltasKernel, 1);
        applyDeltasKernel.setArg<cl::Buffer>(2, cl.getIntegrationUtilities().getPosDelta().getDeviceBuffer());
    }
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    cl.clearBuffer(integration.getPosDelta());
    integration.applyConstraints(tol);
    cl.executeKernel(applyDeltasKernel, cl.getNumAtoms());
    integration.computeVirtualSites();
}

void OpenCLApplyConstraintsKernel::applyToVelocities(ContextImpl& context, double tol) {
    cl.getIntegrationUtilities().applyVelocityConstraints(tol);
}

void OpenCLVirtualSitesKernel::initialize(const System& system) {
}

void OpenCLVirtualSitesKernel::computePositions(ContextImpl& context) {
    cl.getIntegrationUtilities().computeVirtualSites();
}

class OpenCLCalcNonbondedForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const NonbondedForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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
    void getParticlesInGroup(int index, vector<int>& particles) {
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

class OpenCLCalcNonbondedForceKernel::PmeIO : public CalcPmeReciprocalForceKernel::IO {
public:
    PmeIO(OpenCLContext& cl, cl::Kernel addForcesKernel) : cl(cl), addForcesKernel(addForcesKernel) {
        forceTemp.initialize<mm_float4>(cl, cl.getNumAtoms(), "PmeForce");
        addForcesKernel.setArg<cl::Buffer>(0, forceTemp.getDeviceBuffer());
    }
    float* getPosq() {
        cl.getPosq().download(posq);
        return (float*) &posq[0];
    }
    void setForce(float* force) {
        forceTemp.upload(force);
        addForcesKernel.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        cl.executeKernel(addForcesKernel, cl.getNumAtoms());
    }
private:
    OpenCLContext& cl;
    vector<mm_float4> posq;
    OpenCLArray forceTemp;
    cl::Kernel addForcesKernel;
};

class OpenCLCalcNonbondedForceKernel::PmePreComputation : public OpenCLContext::ForcePreComputation {
public:
    PmePreComputation(OpenCLContext& cl, Kernel& pme, CalcPmeReciprocalForceKernel::IO& io) : cl(cl), pme(pme), io(io) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        Vec3 boxVectors[3] = {Vec3(cl.getPeriodicBoxSize().x, 0, 0), Vec3(0, cl.getPeriodicBoxSize().y, 0), Vec3(0, 0, cl.getPeriodicBoxSize().z)};
        pme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxVectors, includeEnergy);
    }
private:
    OpenCLContext& cl;
    Kernel pme;
    CalcPmeReciprocalForceKernel::IO& io;
};

class OpenCLCalcNonbondedForceKernel::PmePostComputation : public OpenCLContext::ForcePostComputation {
public:
    PmePostComputation(Kernel& pme, CalcPmeReciprocalForceKernel::IO& io) : pme(pme), io(io) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return pme.getAs<CalcPmeReciprocalForceKernel>().finishComputation(io);
    }
private:
    Kernel pme;
    CalcPmeReciprocalForceKernel::IO& io;
};

class OpenCLCalcNonbondedForceKernel::SyncQueuePreComputation : public OpenCLContext::ForcePreComputation {
public:
    SyncQueuePreComputation(OpenCLContext& cl, cl::CommandQueue queue, int forceGroup) : cl(cl), queue(queue), forceGroup(forceGroup) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0) {
            vector<cl::Event> events(1);
            cl.getQueue().enqueueMarker(&events[0]);
            queue.enqueueWaitForEvents(events);
        }
    }
private:
    OpenCLContext& cl;
    cl::CommandQueue queue;
    int forceGroup;
};

class OpenCLCalcNonbondedForceKernel::SyncQueuePostComputation : public OpenCLContext::ForcePostComputation {
public:
    SyncQueuePostComputation(OpenCLContext& cl, cl::Event& event, OpenCLArray& pmeEnergyBuffer, int forceGroup) : cl(cl), event(event),
            pmeEnergyBuffer(pmeEnergyBuffer), forceGroup(forceGroup) {
    }
    void setKernel(cl::Kernel kernel) {
        addEnergyKernel = kernel;
        addEnergyKernel.setArg<cl::Buffer>(0, pmeEnergyBuffer.getDeviceBuffer());
        addEnergyKernel.setArg<cl::Buffer>(1, cl.getEnergyBuffer().getDeviceBuffer());
        addEnergyKernel.setArg<cl_int>(2, pmeEnergyBuffer.getSize());
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0) {
            vector<cl::Event> events(1);
            events[0] = event;
            event = cl::Event();
            cl.getQueue().enqueueWaitForEvents(events);
            if (includeEnergy)
                cl.executeKernel(addEnergyKernel, pmeEnergyBuffer.getSize());
        }
        return 0.0;
    }
private:
    OpenCLContext& cl;
    cl::Event& event;
    cl::Kernel addEnergyKernel;
    OpenCLArray& pmeEnergyBuffer;
    int forceGroup;
};

OpenCLCalcNonbondedForceKernel::~OpenCLCalcNonbondedForceKernel() {
    if (sort != NULL)
        delete sort;
    if (fft != NULL)
        delete fft;
    if (dispersionFft != NULL)
        delete dispersionFft;
    if (pmeio != NULL)
        delete pmeio;
}

void OpenCLCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "nonbonded"+cl.intToString(forceIndex)+"_";

    // Identify which exceptions are 1-4 interactions.

    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    vector<pair<int, int> > exclusions;
    vector<int> exceptions;
    map<int, int> exceptionIndex;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions.push_back(pair<int, int>(particle1, particle2));
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end()) {
            exceptionIndex[i] = exceptions.size();
            exceptions.push_back(i);
        }
    }

    // Initialize nonbonded interactions.

    int numParticles = force.getNumParticles();
    vector<mm_float4> baseParticleParamVec(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<vector<int> > exclusionList(numParticles);
    hasCoulomb = false;
    hasLJ = false;
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        baseParticleParamVec[i] = mm_float4(charge, sigma, epsilon, 0);
        exclusionList[i].push_back(i);
        if (charge != 0.0)
            hasCoulomb = true;
        if (epsilon != 0.0)
            hasLJ = true;
    }
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string param;
        int particle;
        double charge, sigma, epsilon;
        force.getParticleParameterOffset(i, param, particle, charge, sigma, epsilon);
        if (charge != 0.0)
            hasCoulomb = true;
        if (epsilon != 0.0)
            hasLJ = true;
    }
    for (auto exclusion : exclusions) {
        exclusionList[exclusion.first].push_back(exclusion.second);
        exclusionList[exclusion.second].push_back(exclusion.first);
    }
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    bool useCutoff = (nonbondedMethod != NoCutoff);
    bool usePeriodic = (nonbondedMethod != NoCutoff && nonbondedMethod != CutoffNonPeriodic);
    doLJPME = (nonbondedMethod == LJPME && hasLJ);
    usePosqCharges = hasCoulomb ? cl.requestPosqCharges() : false;
    map<string, string> defines;
    defines["HAS_COULOMB"] = (hasCoulomb ? "1" : "0");
    defines["HAS_LENNARD_JONES"] = (hasLJ ? "1" : "0");
    defines["USE_LJ_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    if (useCutoff) {
        // Compute the reaction field constants.

        double reactionFieldK = pow(force.getCutoffDistance(), -3.0)*(force.getReactionFieldDielectric()-1.0)/(2.0*force.getReactionFieldDielectric()+1.0);
        double reactionFieldC = (1.0 / force.getCutoffDistance())*(3.0*force.getReactionFieldDielectric())/(2.0*force.getReactionFieldDielectric()+1.0);
        defines["REACTION_FIELD_K"] = cl.doubleToString(reactionFieldK);
        defines["REACTION_FIELD_C"] = cl.doubleToString(reactionFieldC);
        
        // Compute the switching coefficients.
        
        if (force.getUseSwitchingFunction()) {
            defines["LJ_SWITCH_CUTOFF"] = cl.doubleToString(force.getSwitchingDistance());
            defines["LJ_SWITCH_C3"] = cl.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
            defines["LJ_SWITCH_C4"] = cl.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
            defines["LJ_SWITCH_C5"] = cl.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
        }
    }
    if (force.getUseDispersionCorrection() && cl.getContextIndex() == 0 && !doLJPME)
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
    alpha = 0;
    ewaldSelfEnergy = 0.0;
    map<string, string> paramsDefines;
    hasOffsets = (force.getNumParticleParameterOffsets() > 0 || force.getNumExceptionParameterOffsets() > 0);
    if (hasOffsets)
        paramsDefines["HAS_OFFSETS"] = "1";
    if (usePosqCharges)
        paramsDefines["USE_POSQ_CHARGES"] = "1";
    if (nonbondedMethod == Ewald) {
        // Compute the Ewald parameters.

        int kmaxx, kmaxy, kmaxz;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
        defines["EWALD_ALPHA"] = cl.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cl.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        if (cl.getContextIndex() == 0) {
            paramsDefines["INCLUDE_EWALD"] = "1";
            paramsDefines["EWALD_SELF_ENERGY_SCALE"] = cl.doubleToString(ONE_4PI_EPS0*alpha/sqrt(M_PI));
            for (int i = 0; i < numParticles; i++)
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);

            // Create the reciprocal space kernels.

            map<string, string> replacements;
            replacements["NUM_ATOMS"] = cl.intToString(numParticles);
            replacements["KMAX_X"] = cl.intToString(kmaxx);
            replacements["KMAX_Y"] = cl.intToString(kmaxy);
            replacements["KMAX_Z"] = cl.intToString(kmaxz);
            replacements["EXP_COEFFICIENT"] = cl.doubleToString(-1.0/(4.0*alpha*alpha));
            cl::Program program = cl.createProgram(OpenCLKernelSources::ewald, replacements);
            ewaldSumsKernel = cl::Kernel(program, "calculateEwaldCosSinSums");
            ewaldForcesKernel = cl::Kernel(program, "calculateEwaldForces");
            int elementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double2) : sizeof(mm_float2));
            cosSinSums.initialize(cl, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), elementSize, "cosSinSums");
        }
    }
    else if (((nonbondedMethod == PME || nonbondedMethod == LJPME) && hasCoulomb) || doLJPME) {
        // Compute the PME parameters.

        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ, false);
        gridSizeX = OpenCLFFT3D::findLegalDimension(gridSizeX);
        gridSizeY = OpenCLFFT3D::findLegalDimension(gridSizeY);
        gridSizeZ = OpenCLFFT3D::findLegalDimension(gridSizeZ);
        if (doLJPME) {
            NonbondedForceImpl::calcPMEParameters(system, force, dispersionAlpha, dispersionGridSizeX,
                                                  dispersionGridSizeY, dispersionGridSizeZ, true);
            dispersionGridSizeX = OpenCLFFT3D::findLegalDimension(dispersionGridSizeX);
            dispersionGridSizeY = OpenCLFFT3D::findLegalDimension(dispersionGridSizeY);
            dispersionGridSizeZ = OpenCLFFT3D::findLegalDimension(dispersionGridSizeZ);
        }
        defines["EWALD_ALPHA"] = cl.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cl.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        defines["DO_LJPME"] = doLJPME ? "1" : "0";
        if (doLJPME)
            defines["EWALD_DISPERSION_ALPHA"] = cl.doubleToString(dispersionAlpha);
        if (cl.getContextIndex() == 0) {
            paramsDefines["INCLUDE_EWALD"] = "1";
            paramsDefines["EWALD_SELF_ENERGY_SCALE"] = cl.doubleToString(ONE_4PI_EPS0*alpha/sqrt(M_PI));
            for (int i = 0; i < numParticles; i++)
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);
            if (doLJPME) {
                paramsDefines["INCLUDE_LJPME"] = "1";
                paramsDefines["LJPME_SELF_ENERGY_SCALE"] = cl.doubleToString(pow(dispersionAlpha, 6)/3.0);
                for (int i = 0; i < numParticles; i++)
                    ewaldSelfEnergy += baseParticleParamVec[i].z*pow(baseParticleParamVec[i].y*dispersionAlpha, 6)/3.0;
            }
            pmeDefines["PME_ORDER"] = cl.intToString(PmeOrder);
            pmeDefines["NUM_ATOMS"] = cl.intToString(numParticles);
            pmeDefines["RECIP_EXP_FACTOR"] = cl.doubleToString(M_PI*M_PI/(alpha*alpha));
            pmeDefines["GRID_SIZE_X"] = cl.intToString(gridSizeX);
            pmeDefines["GRID_SIZE_Y"] = cl.intToString(gridSizeY);
            pmeDefines["GRID_SIZE_Z"] = cl.intToString(gridSizeZ);
            pmeDefines["EPSILON_FACTOR"] = cl.doubleToString(sqrt(ONE_4PI_EPS0));
            pmeDefines["M_PI"] = cl.doubleToString(M_PI);
            bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
            if (deviceIsCpu)
                pmeDefines["DEVICE_IS_CPU"] = "1";
            if (cl.getPlatformData().useCpuPme && !doLJPME && usePosqCharges) {
                // Create the CPU PME kernel.

                try {
                    cpuPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), *cl.getPlatformData().context);
                    cpuPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSizeX, gridSizeY, gridSizeZ, numParticles, alpha, false);
                    cl::Program program = cl.createProgram(OpenCLKernelSources::pme, pmeDefines);
                    cl::Kernel addForcesKernel = cl::Kernel(program, "addForces");
                    pmeio = new PmeIO(cl, addForcesKernel);
                    cl.addPreComputation(new PmePreComputation(cl, cpuPme, *pmeio));
                    cl.addPostComputation(new PmePostComputation(cpuPme, *pmeio));
                }
                catch (OpenMMException& ex) {
                    // The CPU PME plugin isn't available.
                }
            }
            if (pmeio == NULL) {
                // Create required data structures.

                if (doLJPME) {
                    double invRCut6 = pow(force.getCutoffDistance(), -6);
                    double dalphaR = dispersionAlpha * force.getCutoffDistance();
                    double dar2 = dalphaR*dalphaR;
                    double dar4 = dar2*dar2;
                    double multShift6 = -invRCut6*(1.0 - exp(-dar2) * (1.0 + dar2 + 0.5*dar4));
                    defines["INVCUT6"] = cl.doubleToString(invRCut6);
                    defines["MULTSHIFT6"] = cl.doubleToString(multShift6);
                }
                int elementSize = (cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
                int roundedZSize = PmeOrder*(int) ceil(gridSizeZ/(double) PmeOrder);
                int gridElements = gridSizeX*gridSizeY*roundedZSize;
                if (doLJPME) {
                    roundedZSize = PmeOrder*(int) ceil(dispersionGridSizeZ/(double) PmeOrder);
                    gridElements = max(gridElements, dispersionGridSizeX*dispersionGridSizeY*roundedZSize);
                }
                pmeGrid1.initialize(cl, gridElements, 2*elementSize, "pmeGrid1");
                pmeGrid2.initialize(cl, gridElements, 2*elementSize, "pmeGrid2");
                if (cl.getSupports64BitGlobalAtomics())
                    cl.addAutoclearBuffer(pmeGrid2);
                else
                    cl.addAutoclearBuffer(pmeGrid1);
                pmeBsplineModuliX.initialize(cl, gridSizeX, elementSize, "pmeBsplineModuliX");
                pmeBsplineModuliY.initialize(cl, gridSizeY, elementSize, "pmeBsplineModuliY");
                pmeBsplineModuliZ.initialize(cl, gridSizeZ, elementSize, "pmeBsplineModuliZ");
                if (doLJPME) {
                    pmeDispersionBsplineModuliX.initialize(cl, dispersionGridSizeX, elementSize, "pmeDispersionBsplineModuliX");
                    pmeDispersionBsplineModuliY.initialize(cl, dispersionGridSizeY, elementSize, "pmeDispersionBsplineModuliY");
                    pmeDispersionBsplineModuliZ.initialize(cl, dispersionGridSizeZ, elementSize, "pmeDispersionBsplineModuliZ");
                }
                pmeBsplineTheta.initialize(cl, PmeOrder*numParticles, 4*elementSize, "pmeBsplineTheta");
                pmeAtomRange.initialize<cl_int>(cl, gridSizeX*gridSizeY*gridSizeZ+1, "pmeAtomRange");
                pmeAtomGridIndex.initialize<mm_int2>(cl, numParticles, "pmeAtomGridIndex");
                int energyElementSize = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
                pmeEnergyBuffer.initialize(cl, cl.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize, energyElementSize, "pmeEnergyBuffer");
                cl.clearBuffer(pmeEnergyBuffer);
                sort = new OpenCLSort(cl, new SortTrait(), cl.getNumAtoms());
                fft = new OpenCLFFT3D(cl, gridSizeX, gridSizeY, gridSizeZ, true);
                if (doLJPME)
                    dispersionFft = new OpenCLFFT3D(cl, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, true);
                string vendor = cl.getDevice().getInfo<CL_DEVICE_VENDOR>();
                bool isNvidia = (vendor.size() >= 6 && vendor.substr(0, 6) == "NVIDIA");
                usePmeQueue = (!cl.getPlatformData().disablePmeStream && isNvidia);
                if (usePmeQueue) {
                    pmeDefines["USE_PME_STREAM"] = "1";
                    pmeQueue = cl::CommandQueue(cl.getContext(), cl.getDevice());
                    int recipForceGroup = force.getReciprocalSpaceForceGroup();
                    if (recipForceGroup < 0)
                        recipForceGroup = force.getForceGroup();
                    cl.addPreComputation(new SyncQueuePreComputation(cl, pmeQueue, recipForceGroup));
                    cl.addPostComputation(syncQueue = new SyncQueuePostComputation(cl, pmeSyncEvent, pmeEnergyBuffer, recipForceGroup));
                }

                // Initialize the b-spline moduli.

                for (int grid = 0; grid < 2; grid++) {
                    int xsize, ysize, zsize;
                    OpenCLArray *xmoduli, *ymoduli, *zmoduli;
                    if (grid == 0) {
                        xsize = gridSizeX;
                        ysize = gridSizeY;
                        zsize = gridSizeZ;
                        xmoduli = &pmeBsplineModuliX;
                        ymoduli = &pmeBsplineModuliY;
                        zmoduli = &pmeBsplineModuliZ;
                    }
                    else {
                        if (!doLJPME)
                            continue;
                        xsize = dispersionGridSizeX;
                        ysize = dispersionGridSizeY;
                        zsize = dispersionGridSizeZ;
                        xmoduli = &pmeDispersionBsplineModuliX;
                        ymoduli = &pmeDispersionBsplineModuliY;
                        zmoduli = &pmeDispersionBsplineModuliZ;
                    }
                    int maxSize = max(max(xsize, ysize), zsize);
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

                    for (int dim = 0; dim < 3; dim++) {
                        int ndata = (dim == 0 ? xsize : dim == 1 ? ysize : zsize);
                        vector<cl_double> moduli(ndata);
                        for (int i = 0; i < ndata; i++) {
                            double sc = 0.0;
                            double ss = 0.0;
                            for (int j = 0; j < ndata; j++) {
                                double arg = (2.0*M_PI*i*j)/ndata;
                                sc += bsplines_data[j]*cos(arg);
                                ss += bsplines_data[j]*sin(arg);
                            }
                            moduli[i] = sc*sc+ss*ss;
                        }
                        for (int i = 0; i < ndata; i++)
                        {
                            if (moduli[i] < 1.0e-7)
                                moduli[i] = (moduli[i-1]+moduli[i+1])*0.5f;
                        }
                        if (dim == 0)
                            xmoduli->upload(moduli, true);
                        else if (dim == 1)
                            ymoduli->upload(moduli, true);
                        else
                            zmoduli->upload(moduli, true);
                    }
                }
            }
        }
    }

    // Add code to subtract off the reciprocal part of excluded interactions.

    if ((nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME) && pmeio == NULL) {
        int numContexts = cl.getPlatformData().contexts.size();
        int startIndex = cl.getContextIndex()*force.getNumExceptions()/numContexts;
        int endIndex = (cl.getContextIndex()+1)*force.getNumExceptions()/numContexts;
        int numExclusions = endIndex-startIndex;
        if (numExclusions > 0) {
            paramsDefines["HAS_EXCLUSIONS"] = "1";
            vector<vector<int> > atoms(numExclusions, vector<int>(2));
            exclusionAtoms.initialize<mm_int2>(cl, numExclusions, "exclusionAtoms");
            exclusionParams.initialize<mm_float4>(cl, numExclusions, "exclusionParams");
            vector<mm_int2> exclusionAtomsVec(numExclusions);
            for (int i = 0; i < numExclusions; i++) {
                int j = i+startIndex;
                exclusionAtomsVec[i] = mm_int2(exclusions[j].first, exclusions[j].second);
                atoms[i][0] = exclusions[j].first;
                atoms[i][1] = exclusions[j].second;
            }
            exclusionAtoms.upload(exclusionAtomsVec);
            map<string, string> replacements;
            replacements["PARAMS"] = cl.getBondedUtilities().addArgument(exclusionParams.getDeviceBuffer(), "float4");
            replacements["EWALD_ALPHA"] = cl.doubleToString(alpha);
            replacements["TWO_OVER_SQRT_PI"] = cl.doubleToString(2.0/sqrt(M_PI));
            replacements["DO_LJPME"] = doLJPME ? "1" : "0";
            if (doLJPME)
                replacements["EWALD_DISPERSION_ALPHA"] = cl.doubleToString(dispersionAlpha);
            cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::pmeExclusions, replacements), force.getForceGroup());
        }
    }

    // Add the interaction to the default nonbonded kernel.
    
    string source = cl.replaceStrings(OpenCLKernelSources::coulombLennardJones, defines);
    charges.initialize(cl, cl.getPaddedNumAtoms(), cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "charges");
    baseParticleParams.initialize<mm_float4>(cl, cl.getPaddedNumAtoms(), "baseParticleParams");
    baseParticleParams.upload(baseParticleParamVec);
    map<string, string> replacements;
    if (usePosqCharges) {
        replacements["CHARGE1"] = "posq1.w";
        replacements["CHARGE2"] = "posq2.w";
    }
    else {
        replacements["CHARGE1"] = prefix+"charge1";
        replacements["CHARGE2"] = prefix+"charge2";
    }
    if (hasCoulomb)
        cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"charge", "real", 1, charges.getElementSize(), charges.getDeviceBuffer()));
    sigmaEpsilon.initialize<mm_float2>(cl, cl.getPaddedNumAtoms(), "sigmaEpsilon");
    if (hasLJ) {
        replacements["SIGMA_EPSILON1"] = prefix+"sigmaEpsilon1";
        replacements["SIGMA_EPSILON2"] = prefix+"sigmaEpsilon2";
        cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"sigmaEpsilon", "float", 2, sizeof(cl_float2), sigmaEpsilon.getDeviceBuffer()));
    }
    source = cl.replaceStrings(source, replacements);
    cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());

    // Initialize the exceptions.

    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    if (numExceptions > 0) {
        paramsDefines["HAS_EXCEPTIONS"] = "1";
        exceptionAtoms.resize(numExceptions);
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        exceptionParams.initialize<mm_float4>(cl, numExceptions, "exceptionParams");
        baseExceptionParams.initialize<mm_float4>(cl, numExceptions, "baseExceptionParams");
        vector<mm_float4> baseExceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            baseExceptionParamsVec[i] = mm_float4(chargeProd, sigma, epsilon, 0);
            exceptionAtoms[i] = make_pair(atoms[i][0], atoms[i][1]);
        }
        baseExceptionParams.upload(baseExceptionParamsVec);
        map<string, string> replacements;
        replacements["APPLY_PERIODIC"] = (usePeriodic && force.getExceptionsUsePeriodicBoundaryConditions() ? "1" : "0");
        replacements["PARAMS"] = cl.getBondedUtilities().addArgument(exceptionParams.getDeviceBuffer(), "float4");
        cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::nonbondedExceptions, replacements), force.getForceGroup());
    }
    
    // Initialize parameter offsets.

    vector<vector<mm_float4> > particleOffsetVec(force.getNumParticles());
    vector<vector<mm_float4> > exceptionOffsetVec(force.getNumExceptions());
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string param;
        int particle;
        double charge, sigma, epsilon;
        force.getParticleParameterOffset(i, param, particle, charge, sigma, epsilon);
        auto paramPos = find(paramNames.begin(), paramNames.end(), param);
        int paramIndex;
        if (paramPos == paramNames.end()) {
            paramIndex = paramNames.size();
            paramNames.push_back(param);
        }
        else
            paramIndex = paramPos-paramNames.begin();
        particleOffsetVec[particle].push_back(mm_float4(charge, sigma, epsilon, paramIndex));
    }
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        auto paramPos = find(paramNames.begin(), paramNames.end(), param);
        int paramIndex;
        if (paramPos == paramNames.end()) {
            paramIndex = paramNames.size();
            paramNames.push_back(param);
        }
        else
            paramIndex = paramPos-paramNames.begin();
        exceptionOffsetVec[exceptionIndex[exception]].push_back(mm_float4(charge, sigma, epsilon, paramIndex));
    }
    paramValues.resize(paramNames.size(), 0.0);
    particleParamOffsets.initialize<mm_float4>(cl, max(force.getNumParticleParameterOffsets(), 1), "particleParamOffsets");
    exceptionParamOffsets.initialize<mm_float4>(cl, max(force.getNumExceptionParameterOffsets(), 1), "exceptionParamOffsets");
    particleOffsetIndices.initialize<cl_int>(cl, cl.getPaddedNumAtoms()+1, "particleOffsetIndices");
    exceptionOffsetIndices.initialize<cl_int>(cl, force.getNumExceptions()+1, "exceptionOffsetIndices");
    vector<cl_int> particleOffsetIndicesVec, exceptionOffsetIndicesVec;
    vector<mm_float4> p, e;
    for (int i = 0; i < particleOffsetVec.size(); i++) {
        particleOffsetIndicesVec.push_back(p.size());
        for (int j = 0; j < particleOffsetVec[i].size(); j++)
            p.push_back(particleOffsetVec[i][j]);
    }
    while (particleOffsetIndicesVec.size() < particleOffsetIndices.getSize())
        particleOffsetIndicesVec.push_back(p.size());
    for (int i = 0; i < exceptionOffsetVec.size(); i++) {
        exceptionOffsetIndicesVec.push_back(e.size());
        for (int j = 0; j < exceptionOffsetVec[i].size(); j++)
            e.push_back(exceptionOffsetVec[i][j]);
    }
    exceptionOffsetIndicesVec.push_back(e.size());
    if (force.getNumParticleParameterOffsets() > 0) {
        particleParamOffsets.upload(p);
        particleOffsetIndices.upload(particleOffsetIndicesVec);
    }
    if (force.getNumExceptionParameterOffsets() > 0) {
        exceptionParamOffsets.upload(e);
        exceptionOffsetIndices.upload(exceptionOffsetIndicesVec);
    }
    globalParams.initialize(cl, max((int) paramValues.size(), 1), cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "globalParams");
    recomputeParams = true;
    
    // Initialize the kernel for updating parameters.
    
    cl::Program program = cl.createProgram(OpenCLKernelSources::nonbondedParameters, paramsDefines);
    computeParamsKernel = cl::Kernel(program, "computeParameters");
    computeExclusionParamsKernel = cl::Kernel(program, "computeExclusionParameters");
    info = new ForceInfo(cl.getNonbondedUtilities().getNumForceBuffers(), force);
    cl.addForce(info);
}

double OpenCLCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        int index = 0;
        computeParamsKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        index++;
        computeParamsKernel.setArg<cl::Buffer>(index++, globalParams.getDeviceBuffer());
        computeParamsKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        computeParamsKernel.setArg<cl::Buffer>(index++, baseParticleParams.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, charges.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, sigmaEpsilon.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, particleParamOffsets.getDeviceBuffer());
        computeParamsKernel.setArg<cl::Buffer>(index++, particleOffsetIndices.getDeviceBuffer());
        if (exceptionParams.isInitialized()) {
            computeParamsKernel.setArg<cl_int>(index++, exceptionParams.getSize());
            computeParamsKernel.setArg<cl::Buffer>(index++, baseExceptionParams.getDeviceBuffer());
            computeParamsKernel.setArg<cl::Buffer>(index++, exceptionParams.getDeviceBuffer());
            computeParamsKernel.setArg<cl::Buffer>(index++, exceptionParamOffsets.getDeviceBuffer());
            computeParamsKernel.setArg<cl::Buffer>(index++, exceptionOffsetIndices.getDeviceBuffer());
        }
        if (exclusionParams.isInitialized()) {
            computeExclusionParamsKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl::Buffer>(1, charges.getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl::Buffer>(2, sigmaEpsilon.getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl_int>(3, exclusionParams.getSize());
            computeExclusionParamsKernel.setArg<cl::Buffer>(4, exclusionAtoms.getDeviceBuffer());
            computeExclusionParamsKernel.setArg<cl::Buffer>(5, exclusionParams.getDeviceBuffer());
        }
        if (cosSinSums.isInitialized()) {
            ewaldSumsKernel.setArg<cl::Buffer>(0, cl.getEnergyBuffer().getDeviceBuffer());
            ewaldSumsKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            ewaldSumsKernel.setArg<cl::Buffer>(2, cosSinSums.getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(0, cl.getForceBuffers().getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            ewaldForcesKernel.setArg<cl::Buffer>(2, cosSinSums.getDeviceBuffer());
        }
        if (pmeGrid1.isInitialized()) {
            // Create kernels for Coulomb PME.
            
            map<string, string> replacements;
            replacements["CHARGE"] = (usePosqCharges ? "pos.w" : "charges[atom]");
            cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::pme, replacements), pmeDefines);
            pmeUpdateBsplinesKernel = cl::Kernel(program, "updateBsplines");
            pmeAtomRangeKernel = cl::Kernel(program, "findAtomRangeForGrid");
            pmeZIndexKernel = cl::Kernel(program, "recordZIndex");
            pmeSpreadChargeKernel = cl::Kernel(program, "gridSpreadCharge");
            pmeConvolutionKernel = cl::Kernel(program, "reciprocalConvolution");
            pmeEvalEnergyKernel = cl::Kernel(program, "gridEvaluateEnergy");
            pmeInterpolateForceKernel = cl::Kernel(program, "gridInterpolateForce");
            int elementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(1, pmeBsplineTheta.getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg(2, OpenCLContext::ThreadBlockSize*PmeOrder*elementSize, NULL);
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(3, pmeAtomGridIndex.getDeviceBuffer());
            pmeUpdateBsplinesKernel.setArg<cl::Buffer>(12, charges.getDeviceBuffer());
            pmeAtomRangeKernel.setArg<cl::Buffer>(0, pmeAtomGridIndex.getDeviceBuffer());
            pmeAtomRangeKernel.setArg<cl::Buffer>(1, pmeAtomRange.getDeviceBuffer());
            pmeAtomRangeKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
            pmeZIndexKernel.setArg<cl::Buffer>(0, pmeAtomGridIndex.getDeviceBuffer());
            pmeZIndexKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(1, pmeAtomGridIndex.getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(2, pmeAtomRange.getDeviceBuffer());
            if (cl.getSupports64BitGlobalAtomics())
                pmeSpreadChargeKernel.setArg<cl::Buffer>(3, pmeGrid2.getDeviceBuffer());
            else
                pmeSpreadChargeKernel.setArg<cl::Buffer>(3, pmeGrid1.getDeviceBuffer());
            pmeSpreadChargeKernel.setArg<cl::Buffer>(4, pmeBsplineTheta.getDeviceBuffer());
            if (deviceIsCpu || cl.getSupports64BitGlobalAtomics())
                pmeSpreadChargeKernel.setArg<cl::Buffer>(13, charges.getDeviceBuffer());
            else
                pmeSpreadChargeKernel.setArg<cl::Buffer>(5, charges.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(1, pmeBsplineModuliX.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(2, pmeBsplineModuliY.getDeviceBuffer());
            pmeConvolutionKernel.setArg<cl::Buffer>(3, pmeBsplineModuliZ.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(1, usePmeQueue ? pmeEnergyBuffer.getDeviceBuffer() : cl.getEnergyBuffer().getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(2, pmeBsplineModuliX.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(3, pmeBsplineModuliY.getDeviceBuffer());
            pmeEvalEnergyKernel.setArg<cl::Buffer>(4, pmeBsplineModuliZ.getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(1, cl.getForceBuffers().getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(2, pmeGrid1.getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(11, pmeAtomGridIndex.getDeviceBuffer());
            pmeInterpolateForceKernel.setArg<cl::Buffer>(12, charges.getDeviceBuffer());
            if (cl.getSupports64BitGlobalAtomics()) {
                pmeFinishSpreadChargeKernel = cl::Kernel(program, "finishSpreadCharge");
                pmeFinishSpreadChargeKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                pmeFinishSpreadChargeKernel.setArg<cl::Buffer>(1, pmeGrid1.getDeviceBuffer());
            }
            if (usePmeQueue)
                syncQueue->setKernel(cl::Kernel(program, "addEnergy"));

            if (doLJPME) {
                // Create kernels for LJ PME.

                pmeDefines["EWALD_ALPHA"] = cl.doubleToString(dispersionAlpha);
                pmeDefines["GRID_SIZE_X"] = cl.intToString(dispersionGridSizeX);
                pmeDefines["GRID_SIZE_Y"] = cl.intToString(dispersionGridSizeY);
                pmeDefines["GRID_SIZE_Z"] = cl.intToString(dispersionGridSizeZ);
                pmeDefines["EPSILON_FACTOR"] = "1";
                pmeDefines["RECIP_EXP_FACTOR"] = cl.doubleToString(M_PI*M_PI/(dispersionAlpha*dispersionAlpha));
                pmeDefines["USE_LJPME"] = "1";
                program = cl.createProgram(OpenCLKernelSources::pme, pmeDefines);
                pmeDispersionUpdateBsplinesKernel = cl::Kernel(program, "updateBsplines");
                pmeDispersionAtomRangeKernel = cl::Kernel(program, "findAtomRangeForGrid");
                pmeDispersionZIndexKernel = cl::Kernel(program, "recordZIndex");
                pmeDispersionSpreadChargeKernel = cl::Kernel(program, "gridSpreadCharge");
                pmeDispersionConvolutionKernel = cl::Kernel(program, "reciprocalConvolution");
                pmeDispersionEvalEnergyKernel = cl::Kernel(program, "gridEvaluateEnergy");
                pmeDispersionInterpolateForceKernel = cl::Kernel(program, "gridInterpolateForce");
                int elementSize = (cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
                pmeDispersionUpdateBsplinesKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
                pmeDispersionUpdateBsplinesKernel.setArg<cl::Buffer>(1, pmeBsplineTheta.getDeviceBuffer());
                pmeDispersionUpdateBsplinesKernel.setArg(2, OpenCLContext::ThreadBlockSize*PmeOrder*elementSize, NULL);
                pmeDispersionUpdateBsplinesKernel.setArg<cl::Buffer>(3, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionUpdateBsplinesKernel.setArg<cl::Buffer>(12, sigmaEpsilon.getDeviceBuffer());
                pmeDispersionAtomRangeKernel.setArg<cl::Buffer>(0, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionAtomRangeKernel.setArg<cl::Buffer>(1, pmeAtomRange.getDeviceBuffer());
                pmeDispersionAtomRangeKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
                pmeDispersionZIndexKernel.setArg<cl::Buffer>(0, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionZIndexKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(1, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(2, pmeAtomRange.getDeviceBuffer());
                if (cl.getSupports64BitGlobalAtomics())
                    pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(3, pmeGrid2.getDeviceBuffer());
                else
                    pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(3, pmeGrid1.getDeviceBuffer());
                pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(4, pmeBsplineTheta.getDeviceBuffer());
                if (deviceIsCpu || cl.getSupports64BitGlobalAtomics())
                    pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(13, sigmaEpsilon.getDeviceBuffer());
                else
                    pmeDispersionSpreadChargeKernel.setArg<cl::Buffer>(5, sigmaEpsilon.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(1, pmeDispersionBsplineModuliX.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(2, pmeDispersionBsplineModuliY.getDeviceBuffer());
                pmeDispersionConvolutionKernel.setArg<cl::Buffer>(3, pmeDispersionBsplineModuliZ.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(1, usePmeQueue ? pmeEnergyBuffer.getDeviceBuffer() : cl.getEnergyBuffer().getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(2, pmeDispersionBsplineModuliX.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(3, pmeDispersionBsplineModuliY.getDeviceBuffer());
                pmeDispersionEvalEnergyKernel.setArg<cl::Buffer>(4, pmeDispersionBsplineModuliZ.getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(1, cl.getForceBuffers().getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(2, pmeGrid1.getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(11, pmeAtomGridIndex.getDeviceBuffer());
                pmeDispersionInterpolateForceKernel.setArg<cl::Buffer>(12, sigmaEpsilon.getDeviceBuffer());
                if (cl.getSupports64BitGlobalAtomics()) {
                    pmeDispersionFinishSpreadChargeKernel = cl::Kernel(program, "finishSpreadCharge");
                    pmeDispersionFinishSpreadChargeKernel.setArg<cl::Buffer>(0, pmeGrid2.getDeviceBuffer());
                    pmeDispersionFinishSpreadChargeKernel.setArg<cl::Buffer>(1, pmeGrid1.getDeviceBuffer());
                }
            }
       }
    }
    
    // Update particle and exception parameters.

    bool paramChanged = false;
    for (int i = 0; i < paramNames.size(); i++) {
        double value = context.getParameter(paramNames[i]);
        if (value != paramValues[i]) {
            paramValues[i] = value;;
            paramChanged = true;
        }
    }
    if (paramChanged) {
        recomputeParams = true;
        globalParams.upload(paramValues, true);
    }
    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
    if (recomputeParams || hasOffsets) {
        computeParamsKernel.setArg<cl_int>(1, includeEnergy && includeReciprocal);
        cl.executeKernel(computeParamsKernel, cl.getPaddedNumAtoms());
        if (exclusionParams.isInitialized())
            cl.executeKernel(computeExclusionParamsKernel, exclusionParams.getSize());
        if (usePmeQueue) {
            vector<cl::Event> events(1);
            cl.getQueue().enqueueMarker(&events[0]);
            pmeQueue.enqueueWaitForEvents(events);
        }
        if (hasOffsets)
            energy = 0.0; // The Ewald self energy was computed in the kernel.
        recomputeParams = false;
    }
    
    // Do reciprocal space calculations.
    
    if (cosSinSums.isInitialized() && includeReciprocal) {
        mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
        mm_double4 recipBoxSize = mm_double4(2*M_PI/boxSize.x, 2*M_PI/boxSize.y, 2*M_PI/boxSize.z, 0.0);
        double recipCoefficient = ONE_4PI_EPS0*4*M_PI/(boxSize.x*boxSize.y*boxSize.z);
        if (cl.getUseDoublePrecision()) {
            ewaldSumsKernel.setArg<mm_double4>(3, recipBoxSize);
            ewaldSumsKernel.setArg<cl_double>(4, recipCoefficient);
            ewaldForcesKernel.setArg<mm_double4>(3, recipBoxSize);
            ewaldForcesKernel.setArg<cl_double>(4, recipCoefficient);
        }
        else {
            ewaldSumsKernel.setArg<mm_float4>(3, mm_float4((float) recipBoxSize.x, (float) recipBoxSize.y, (float) recipBoxSize.z, 0));
            ewaldSumsKernel.setArg<cl_float>(4, (cl_float) recipCoefficient);
            ewaldForcesKernel.setArg<mm_float4>(3, mm_float4((float) recipBoxSize.x, (float) recipBoxSize.y, (float) recipBoxSize.z, 0));
            ewaldForcesKernel.setArg<cl_float>(4, (cl_float) recipCoefficient);
        }
        cl.executeKernel(ewaldSumsKernel, cosSinSums.getSize());
        cl.executeKernel(ewaldForcesKernel, cl.getNumAtoms());
    }
    if (pmeGrid1.isInitialized() && includeReciprocal) {
        if (usePmeQueue && !includeEnergy)
            cl.setQueue(pmeQueue);
        
        // Invert the periodic box vectors.
        
        Vec3 boxVectors[3];
        cl.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double determinant = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
        double scale = 1.0/determinant;
        mm_double4 recipBoxVectors[3];
        recipBoxVectors[0] = mm_double4(boxVectors[1][1]*boxVectors[2][2]*scale, 0, 0, 0);
        recipBoxVectors[1] = mm_double4(-boxVectors[1][0]*boxVectors[2][2]*scale, boxVectors[0][0]*boxVectors[2][2]*scale, 0, 0);
        recipBoxVectors[2] = mm_double4((boxVectors[1][0]*boxVectors[2][1]-boxVectors[1][1]*boxVectors[2][0])*scale, -boxVectors[0][0]*boxVectors[2][1]*scale, boxVectors[0][0]*boxVectors[1][1]*scale, 0);
        mm_float4 recipBoxVectorsFloat[3];
        for (int i = 0; i < 3; i++)
            recipBoxVectorsFloat[i] = mm_float4((float) recipBoxVectors[i].x, (float) recipBoxVectors[i].y, (float) recipBoxVectors[i].z, 0);
        
        // Execute the reciprocal space kernels.

        if (hasCoulomb) {
            setPeriodicBoxArgs(cl, pmeUpdateBsplinesKernel, 4);
            if (cl.getUseDoublePrecision()) {
                pmeUpdateBsplinesKernel.setArg<mm_double4>(9, recipBoxVectors[0]);
                pmeUpdateBsplinesKernel.setArg<mm_double4>(10, recipBoxVectors[1]);
                pmeUpdateBsplinesKernel.setArg<mm_double4>(11, recipBoxVectors[2]);
            }
            else {
                pmeUpdateBsplinesKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[0]);
                pmeUpdateBsplinesKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[1]);
                pmeUpdateBsplinesKernel.setArg<mm_float4>(11, recipBoxVectorsFloat[2]);
            }
            cl.executeKernel(pmeUpdateBsplinesKernel, cl.getNumAtoms());
            if (deviceIsCpu && !cl.getSupports64BitGlobalAtomics()) {
                setPeriodicBoxArgs(cl, pmeSpreadChargeKernel, 5);
                if (cl.getUseDoublePrecision()) {
                    pmeSpreadChargeKernel.setArg<mm_double4>(10, recipBoxVectors[0]);
                    pmeSpreadChargeKernel.setArg<mm_double4>(11, recipBoxVectors[1]);
                    pmeSpreadChargeKernel.setArg<mm_double4>(12, recipBoxVectors[2]);
                }
                else {
                    pmeSpreadChargeKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[0]);
                    pmeSpreadChargeKernel.setArg<mm_float4>(11, recipBoxVectorsFloat[1]);
                    pmeSpreadChargeKernel.setArg<mm_float4>(12, recipBoxVectorsFloat[2]);
                }
                cl.executeKernel(pmeSpreadChargeKernel, 2*cl.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
            }
            else {
                sort->sort(pmeAtomGridIndex);
                if (cl.getSupports64BitGlobalAtomics()) {
                    setPeriodicBoxArgs(cl, pmeSpreadChargeKernel, 5);
                    if (cl.getUseDoublePrecision()) {
                        pmeSpreadChargeKernel.setArg<mm_double4>(10, recipBoxVectors[0]);
                        pmeSpreadChargeKernel.setArg<mm_double4>(11, recipBoxVectors[1]);
                        pmeSpreadChargeKernel.setArg<mm_double4>(12, recipBoxVectors[2]);
                    }
                    else {
                        pmeSpreadChargeKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[0]);
                        pmeSpreadChargeKernel.setArg<mm_float4>(11, recipBoxVectorsFloat[1]);
                        pmeSpreadChargeKernel.setArg<mm_float4>(12, recipBoxVectorsFloat[2]);
                    }
                    cl.executeKernel(pmeSpreadChargeKernel, cl.getNumAtoms());
                    cl.executeKernel(pmeFinishSpreadChargeKernel, gridSizeX*gridSizeY*gridSizeZ);
                }
                else {
                    cl.executeKernel(pmeAtomRangeKernel, cl.getNumAtoms());
                    setPeriodicBoxSizeArg(cl, pmeZIndexKernel, 2);
                    if (cl.getUseDoublePrecision())
                        pmeZIndexKernel.setArg<mm_double4>(3, recipBoxVectors[2]);
                    else
                        pmeZIndexKernel.setArg<mm_float4>(3, recipBoxVectorsFloat[2]);
                    cl.executeKernel(pmeZIndexKernel, cl.getNumAtoms());
                    cl.executeKernel(pmeSpreadChargeKernel, cl.getNumAtoms());
                }
            }
            fft->execFFT(pmeGrid1, pmeGrid2, true);
            mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
            if (cl.getUseDoublePrecision()) {
                pmeConvolutionKernel.setArg<mm_double4>(4, recipBoxVectors[0]);
                pmeConvolutionKernel.setArg<mm_double4>(5, recipBoxVectors[1]);
                pmeConvolutionKernel.setArg<mm_double4>(6, recipBoxVectors[2]);
                pmeEvalEnergyKernel.setArg<mm_double4>(5, recipBoxVectors[0]);
                pmeEvalEnergyKernel.setArg<mm_double4>(6, recipBoxVectors[1]);
                pmeEvalEnergyKernel.setArg<mm_double4>(7, recipBoxVectors[2]);
            }
            else {
                pmeConvolutionKernel.setArg<mm_float4>(4, recipBoxVectorsFloat[0]);
                pmeConvolutionKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[1]);
                pmeConvolutionKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[2]);
                pmeEvalEnergyKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[0]);
                pmeEvalEnergyKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[1]);
                pmeEvalEnergyKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[2]);
            }
            if (includeEnergy)
                cl.executeKernel(pmeEvalEnergyKernel, gridSizeX*gridSizeY*gridSizeZ);
            cl.executeKernel(pmeConvolutionKernel, gridSizeX*gridSizeY*gridSizeZ);
            fft->execFFT(pmeGrid2, pmeGrid1, false);
            setPeriodicBoxArgs(cl, pmeInterpolateForceKernel, 3);
            if (cl.getUseDoublePrecision()) {
                pmeInterpolateForceKernel.setArg<mm_double4>(8, recipBoxVectors[0]);
                pmeInterpolateForceKernel.setArg<mm_double4>(9, recipBoxVectors[1]);
                pmeInterpolateForceKernel.setArg<mm_double4>(10, recipBoxVectors[2]);
            }
            else {
                pmeInterpolateForceKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[0]);
                pmeInterpolateForceKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[1]);
                pmeInterpolateForceKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[2]);
            }
            if (deviceIsCpu)
                cl.executeKernel(pmeInterpolateForceKernel, 2*cl.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
            else
                cl.executeKernel(pmeInterpolateForceKernel, cl.getNumAtoms());
        }
        
        if (doLJPME && hasLJ) {
            setPeriodicBoxArgs(cl, pmeDispersionUpdateBsplinesKernel, 4);
            if (cl.getUseDoublePrecision()) {
                pmeDispersionUpdateBsplinesKernel.setArg<mm_double4>(9, recipBoxVectors[0]);
                pmeDispersionUpdateBsplinesKernel.setArg<mm_double4>(10, recipBoxVectors[1]);
                pmeDispersionUpdateBsplinesKernel.setArg<mm_double4>(11, recipBoxVectors[2]);
            }
            else {
                pmeDispersionUpdateBsplinesKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[0]);
                pmeDispersionUpdateBsplinesKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[1]);
                pmeDispersionUpdateBsplinesKernel.setArg<mm_float4>(11, recipBoxVectorsFloat[2]);
            }
            cl.executeKernel(pmeDispersionUpdateBsplinesKernel, cl.getNumAtoms());
            if (deviceIsCpu && !cl.getSupports64BitGlobalAtomics()) {
                cl.clearBuffer(pmeGrid1);
                setPeriodicBoxArgs(cl, pmeDispersionSpreadChargeKernel, 5);
                if (cl.getUseDoublePrecision()) {
                    pmeDispersionSpreadChargeKernel.setArg<mm_double4>(10, recipBoxVectors[0]);
                    pmeDispersionSpreadChargeKernel.setArg<mm_double4>(11, recipBoxVectors[1]);
                    pmeDispersionSpreadChargeKernel.setArg<mm_double4>(12, recipBoxVectors[2]);
                }
                else {
                    pmeDispersionSpreadChargeKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[0]);
                    pmeDispersionSpreadChargeKernel.setArg<mm_float4>(11, recipBoxVectorsFloat[1]);
                    pmeDispersionSpreadChargeKernel.setArg<mm_float4>(12, recipBoxVectorsFloat[2]);
                }
                cl.executeKernel(pmeDispersionSpreadChargeKernel, 2*cl.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
            }
            else {
                if (cl.getSupports64BitGlobalAtomics()) {
                    if (!hasCoulomb)
                        sort->sort(pmeAtomGridIndex);
                    cl.clearBuffer(pmeGrid2);
                    setPeriodicBoxArgs(cl, pmeDispersionSpreadChargeKernel, 5);
                    if (cl.getUseDoublePrecision()) {
                        pmeDispersionSpreadChargeKernel.setArg<mm_double4>(10, recipBoxVectors[0]);
                        pmeDispersionSpreadChargeKernel.setArg<mm_double4>(11, recipBoxVectors[1]);
                        pmeDispersionSpreadChargeKernel.setArg<mm_double4>(12, recipBoxVectors[2]);
                    }
                    else {
                        pmeDispersionSpreadChargeKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[0]);
                        pmeDispersionSpreadChargeKernel.setArg<mm_float4>(11, recipBoxVectorsFloat[1]);
                        pmeDispersionSpreadChargeKernel.setArg<mm_float4>(12, recipBoxVectorsFloat[2]);
                    }
                    cl.executeKernel(pmeDispersionSpreadChargeKernel, cl.getNumAtoms());
                    cl.executeKernel(pmeDispersionFinishSpreadChargeKernel, gridSizeX*gridSizeY*gridSizeZ);
                }
                else {
                    sort->sort(pmeAtomGridIndex);
                    cl.clearBuffer(pmeGrid1);
                    cl.executeKernel(pmeDispersionAtomRangeKernel, cl.getNumAtoms());
                    setPeriodicBoxSizeArg(cl, pmeDispersionZIndexKernel, 2);
                    if (cl.getUseDoublePrecision())
                        pmeDispersionZIndexKernel.setArg<mm_double4>(3, recipBoxVectors[2]);
                    else
                        pmeDispersionZIndexKernel.setArg<mm_float4>(3, recipBoxVectorsFloat[2]);
                    cl.executeKernel(pmeDispersionZIndexKernel, cl.getNumAtoms());
                    cl.executeKernel(pmeDispersionSpreadChargeKernel, cl.getNumAtoms());
                }
            }
            dispersionFft->execFFT(pmeGrid1, pmeGrid2, true);
            mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
            if (cl.getUseDoublePrecision()) {
                pmeDispersionConvolutionKernel.setArg<mm_double4>(4, recipBoxVectors[0]);
                pmeDispersionConvolutionKernel.setArg<mm_double4>(5, recipBoxVectors[1]);
                pmeDispersionConvolutionKernel.setArg<mm_double4>(6, recipBoxVectors[2]);
                pmeDispersionEvalEnergyKernel.setArg<mm_double4>(5, recipBoxVectors[0]);
                pmeDispersionEvalEnergyKernel.setArg<mm_double4>(6, recipBoxVectors[1]);
                pmeDispersionEvalEnergyKernel.setArg<mm_double4>(7, recipBoxVectors[2]);
            }
            else {
                pmeDispersionConvolutionKernel.setArg<mm_float4>(4, recipBoxVectorsFloat[0]);
                pmeDispersionConvolutionKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[1]);
                pmeDispersionConvolutionKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[2]);
                pmeDispersionEvalEnergyKernel.setArg<mm_float4>(5, recipBoxVectorsFloat[0]);
                pmeDispersionEvalEnergyKernel.setArg<mm_float4>(6, recipBoxVectorsFloat[1]);
                pmeDispersionEvalEnergyKernel.setArg<mm_float4>(7, recipBoxVectorsFloat[2]);
            }
            if (!hasCoulomb) cl.clearBuffer(pmeEnergyBuffer);
            if (includeEnergy)
                cl.executeKernel(pmeDispersionEvalEnergyKernel, gridSizeX*gridSizeY*gridSizeZ);
            cl.executeKernel(pmeDispersionConvolutionKernel, gridSizeX*gridSizeY*gridSizeZ);
            dispersionFft->execFFT(pmeGrid2, pmeGrid1, false);
            setPeriodicBoxArgs(cl, pmeDispersionInterpolateForceKernel, 3);
            if (cl.getUseDoublePrecision()) {
                pmeDispersionInterpolateForceKernel.setArg<mm_double4>(8, recipBoxVectors[0]);
                pmeDispersionInterpolateForceKernel.setArg<mm_double4>(9, recipBoxVectors[1]);
                pmeDispersionInterpolateForceKernel.setArg<mm_double4>(10, recipBoxVectors[2]);
            }
            else {
                pmeDispersionInterpolateForceKernel.setArg<mm_float4>(8, recipBoxVectorsFloat[0]);
                pmeDispersionInterpolateForceKernel.setArg<mm_float4>(9, recipBoxVectorsFloat[1]);
                pmeDispersionInterpolateForceKernel.setArg<mm_float4>(10, recipBoxVectorsFloat[2]);
            }
            if (deviceIsCpu)
                cl.executeKernel(pmeDispersionInterpolateForceKernel, 2*cl.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
            else
                cl.executeKernel(pmeDispersionInterpolateForceKernel, cl.getNumAtoms());
        }
        if (usePmeQueue) {
            pmeQueue.enqueueMarker(&pmeSyncEvent);
            cl.restoreDefaultQueue();
        }
    }
    if (dispersionCoefficient != 0.0 && includeDirect) {
        mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
        energy += dispersionCoefficient/(boxSize.x*boxSize.y*boxSize.z);
    }
    return energy;
}

void OpenCLCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
    // Make sure the new parameters are acceptable.
    
    if (force.getNumParticles() != cl.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    if (!hasCoulomb || !hasLJ) {
        for (int i = 0; i < force.getNumParticles(); i++) {
            double charge, sigma, epsilon;
            force.getParticleParameters(i, charge, sigma, epsilon);
            if (!hasCoulomb && charge != 0.0)
                throw OpenMMException("updateParametersInContext: The nonbonded force kernel does not include Coulomb interactions, because all charges were originally 0");
            if (!hasLJ && epsilon != 0.0)
                throw OpenMMException("updateParametersInContext: The nonbonded force kernel does not include Lennard-Jones interactions, because all epsilons were originally 0");
        }
    }
    vector<int> exceptions;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (exceptionAtoms.size() > exceptions.size() && make_pair(particle1, particle2) == exceptionAtoms[exceptions.size()])
            exceptions.push_back(i);
        else if (chargeProd != 0.0 || epsilon != 0.0)
            throw OpenMMException("updateParametersInContext: The set of non-excluded exceptions has changed");
    }
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    
    // Record the per-particle parameters.
    
    vector<mm_float4> baseParticleParamVec(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        baseParticleParamVec[i] = mm_float4(charge, sigma, epsilon, 0);
    }
    baseParticleParams.upload(baseParticleParamVec);
    
    // Record the exceptions.
    
    if (numExceptions > 0) {
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        vector<mm_float4> baseExceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            baseExceptionParamsVec[i] = mm_float4(chargeProd, sigma, epsilon, 0);
        }
        baseExceptionParams.upload(baseExceptionParamsVec);
    }
    
    // Compute other values.
    
    ewaldSelfEnergy = 0.0;
    if (nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME) {
        if (cl.getContextIndex() == 0) {
            for (int i = 0; i < force.getNumParticles(); i++) {
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);
                if (doLJPME)
                    ewaldSelfEnergy += baseParticleParamVec[i].z*pow(baseParticleParamVec[i].y*dispersionAlpha, 6)/3.0;
            }
        }
    }
    if (force.getUseDispersionCorrection() && cl.getContextIndex() == 0 && (nonbondedMethod == CutoffPeriodic || nonbondedMethod == Ewald || nonbondedMethod == PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
    cl.invalidateMolecules(info);
    recomputeParams = true;
}

void OpenCLCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (cl.getPlatformData().useCpuPme)
        cpuPme.getAs<CalcPmeReciprocalForceKernel>().getPMEParameters(alpha, nx, ny, nz);
    else {
        alpha = this->alpha;
        nx = gridSizeX;
        ny = gridSizeY;
        nz = gridSizeZ;
    }
}

void OpenCLCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != LJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (cl.getPlatformData().useCpuPme)
        //cpuPme.getAs<CalcPmeReciprocalForceKernel>().getLJPMEParameters(alpha, nx, ny, nz);
        throw OpenMMException("getPMEParametersInContext: CPUPME has not been implemented for LJPME yet.");
    else {
        alpha = this->dispersionAlpha;
        nx = dispersionGridSizeX;
        ny = dispersionGridSizeY;
        nz = dispersionGridSizeZ;
    }
}

void OpenCLIntegrateVelocityVerletStepKernel::initialize(const System& system, const NoseHooverIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    map<string, string> defines;
    defines["BOLTZ"] = cl.doubleToString(BOLTZ);
    cl::Program program = cl.createProgram(OpenCLKernelSources::velocityVerlet, defines, "");
    kernel1 = cl::Kernel(program, "integrateVelocityVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVelocityVerletPart2");
    kernel3 = cl::Kernel(program, "integrateVelocityVerletPart3");
    kernelHardWall = cl::Kernel(program, "integrateVelocityVerletHardWall");
    prevMaxPairDistance = (cl_float) -1.0;
    maxPairDistanceBuffer.initialize<cl_float>(cl, 1, "maxPairDistanceBuffer");
}

void OpenCLIntegrateVelocityVerletStepKernel::execute(ContextImpl& context, const NoseHooverIntegrator& integrator, bool &forcesAreValid) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int paddedNumAtoms = cl.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    cl.getIntegrationUtilities().setNextStepSize(dt);

    if (!forcesAreValid) context.calcForcesAndEnergy(true, false);

    const auto& atomList = integrator.getAllThermostatedIndividualParticles();
    const auto& pairList = integrator.getAllThermostatedPairs();
    int numAtoms = atomList.size();
    int numPairs = pairList.size();
    int numParticles = numAtoms + 2*numPairs;
    float maxPairDistance = integrator.getMaximumPairDistance();
    // Make sure atom and pair metadata is uploaded and has the correct dimensions
    if (prevMaxPairDistance != maxPairDistance) {
        std::vector<float> tmp(1, maxPairDistance);
        maxPairDistanceBuffer.upload(tmp);
        prevMaxPairDistance = maxPairDistance;
    }
    if (numAtoms !=0 && (!atomListBuffer.isInitialized() || atomListBuffer.getSize() != numAtoms)) {
        if (atomListBuffer.isInitialized())
            atomListBuffer.resize(atomList.size());
        else
            atomListBuffer.initialize<cl_int>(cl, atomList.size(), "atomListBuffer");
        atomListBuffer.upload(atomList);
    }
    if (numPairs !=0 && (!pairListBuffer.isInitialized() || pairListBuffer.getSize() != numPairs)) {
        if (pairListBuffer.isInitialized()) {
            pairListBuffer.resize(pairList.size());
            pairTemperatureBuffer.resize(pairList.size());
        }
        else {
            pairListBuffer.initialize<mm_int2>(cl, pairList.size(), "pairListBuffer");
            pairTemperatureBuffer.initialize<cl_float>(cl, pairList.size(), "pairTemperatureBuffer");
        }
        std::vector<mm_int2> tmp;
        std::vector<float> tmp2;
        for(const auto &pair : pairList) {
            tmp.push_back(mm_int2(std::get<0>(pair), std::get<1>(pair)));
            tmp2.push_back(std::get<2>(pair));
        }
        pairListBuffer.upload(tmp);
        pairTemperatureBuffer.upload(tmp2);
    }

//// Call the first integration kernel.
    kernel1.setArg<cl_int>(0, numAtoms);
    kernel1.setArg<cl_int>(1, numPairs);
    kernel1.setArg<cl_int>(2, paddedNumAtoms);
    kernel1.setArg<cl::Buffer>(3, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    setPosqCorrectionArg(cl, kernel1, 5);
    kernel1.setArg<cl::Buffer>(6, cl.getVelm().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(7, cl.getForce().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(8, integration.getPosDelta().getDeviceBuffer());
    if (numAtoms > 0)
        kernel1.setArg<cl::Buffer>(9, atomListBuffer.getDeviceBuffer());
    else
        kernel1.setArg<void*>(9, NULL);
    if (numPairs > 0)
        kernel1.setArg<cl::Buffer>(10, pairListBuffer.getDeviceBuffer());
    else
        kernel1.setArg<void*>(10, NULL);
    cl.executeKernel(kernel1, std::max(numAtoms, numPairs));

    //// Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    //// Call the second integration kernel.
    kernel2.setArg<cl_int>(0, numParticles);
    kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
    setPosqCorrectionArg(cl, kernel2, 3);
    kernel2.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
    cl.executeKernel(kernel2, numParticles);

    if (numPairs > 0) {
        //// Enforce hard wall constraint
        kernelHardWall.setArg<cl_int>(0, numPairs);
        kernelHardWall.setArg<cl::Buffer>(1, maxPairDistanceBuffer.getDeviceBuffer());
        kernelHardWall.setArg<cl::Buffer>(2, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernelHardWall.setArg<cl::Buffer>(3, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernelHardWall, 4);
        kernelHardWall.setArg<cl::Buffer>(5, cl.getVelm().getDeviceBuffer());
        kernelHardWall.setArg<cl::Buffer>(6, pairListBuffer.getDeviceBuffer());
        kernelHardWall.setArg<cl::Buffer>(7, pairTemperatureBuffer.getDeviceBuffer());
        cl.executeKernel(kernelHardWall, numPairs);
    }


    integration.computeVirtualSites();

    //// Update forces
    context.calcForcesAndEnergy(true, false);
    forcesAreValid = true;

    //// Call the third integration kernel.
    kernel3.setArg<cl_int>(0, numAtoms);
    kernel3.setArg<cl_int>(1, numPairs);
    kernel3.setArg<cl_int>(2, paddedNumAtoms);
    kernel3.setArg<cl::Buffer>(3, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
    kernel3.setArg<cl::Buffer>(4, cl.getPosq().getDeviceBuffer());
    setPosqCorrectionArg(cl, kernel3, 5);
    kernel3.setArg<cl::Buffer>(6, cl.getVelm().getDeviceBuffer());
    kernel3.setArg<cl::Buffer>(7, cl.getForce().getDeviceBuffer());
    kernel3.setArg<cl::Buffer>(8, integration.getPosDelta().getDeviceBuffer());
    if (numAtoms > 0)
        kernel3.setArg<cl::Buffer>(9, atomListBuffer.getDeviceBuffer());
    else
        kernel3.setArg<void*>(9, NULL);
    if (numPairs > 0)
        kernel3.setArg<cl::Buffer>(10, pairListBuffer.getDeviceBuffer());
    else
        kernel3.setArg<void*>(10, NULL);
    cl.executeKernel(kernel3, std::max(numAtoms, numPairs));

    integration.applyVelocityConstraints(integrator.getConstraintTolerance());

    //// Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
}

double OpenCLIntegrateVelocityVerletStepKernel::computeKineticEnergy(ContextImpl& context, const NoseHooverIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0);
}

class OpenCLCalcCustomCVForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(ComputeForceInfo& force) : OpenCLForceInfo(0), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        return force.areParticlesIdentical(particle1, particle2);
    }
    int getNumParticleGroups() {
        return force.getNumParticleGroups();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        force.getParticlesInGroup(index, particles);
    }
    bool areGroupsIdentical(int group1, int group2) {
        return force.areGroupsIdentical(group1, group2);
    }
private:
    ComputeForceInfo& force;
};

class OpenCLCalcCustomCVForceKernel::ReorderListener : public OpenCLContext::ReorderListener {
public:
    ReorderListener(OpenCLContext& cl, OpenCLArray& invAtomOrder) : cl(cl), invAtomOrder(invAtomOrder) {
    }
    void execute() {
        vector<cl_int> invOrder(cl.getPaddedNumAtoms());
        const vector<int>& order = cl.getAtomIndex();
        for (int i = 0; i < order.size(); i++)
            invOrder[order[i]] = i;
        invAtomOrder.upload(invOrder);
    }
private:
    OpenCLContext& cl;
    OpenCLArray& invAtomOrder;
};

void OpenCLCalcCustomCVForceKernel::initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext) {
    int numCVs = force.getNumCollectiveVariables();
    cl.addForce(new OpenCLForceInfo(1));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < numCVs; i++)
        variableNames.push_back(force.getCollectiveVariableName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string name = force.getEnergyParameterDerivativeName(i);
        paramDerivNames.push_back(name);
        cl.addEnergyParameterDerivative(name);
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Create the expressions.

    Lepton::ParsedExpression energyExpr = Lepton::Parser::parse(force.getEnergyFunction(), functions);
    energyExpression = energyExpr.createProgram();
    variableDerivExpressions.clear();
    for (auto& name : variableNames)
        variableDerivExpressions.push_back(energyExpr.differentiate(name).optimize().createProgram());
    paramDerivExpressions.clear();
    for (auto& name : paramDerivNames)
        paramDerivExpressions.push_back(energyExpr.differentiate(name).optimize().createProgram());

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;

    // Copy parameter derivatives from the inner context.

    OpenCLContext& cl2 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    for (auto& param : cl2.getEnergyParamDerivNames())
        cl.addEnergyParameterDerivative(param);
    
    // Create arrays for storing information.
    
    int elementSize = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    cvForces.resize(numCVs);
    for (int i = 0; i < numCVs; i++)
        cvForces[i].initialize(cl, cl.getNumAtoms(), 4*elementSize, "cvForce");
    invAtomOrder.initialize<cl_int>(cl, cl.getPaddedNumAtoms(), "invAtomOrder");
    innerInvAtomOrder.initialize<cl_int>(cl, cl.getPaddedNumAtoms(), "innerInvAtomOrder");
    
    // Create the kernels.
    
    stringstream args, add;
    for (int i = 0; i < numCVs; i++) {
        args << ", __global real4* restrict force" << i << ", real dEdV" << i;
        add << "f += force" << i << "[i]*dEdV" << i << ";\n";
    }
    map<string, string> replacements;
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    replacements["ADD_FORCES"] = add.str();
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customCVForce, replacements));
    copyStateKernel = cl::Kernel(program, "copyState");
    copyForcesKernel = cl::Kernel(program, "copyForces");
    addForcesKernel = cl::Kernel(program, "addForces");

    // This context needs to respect all forces in the inner context when reordering atoms.

    for (auto* info : cl2.getForceInfos())
        cl.addForce(new ForceInfo(*info));
}

double OpenCLCalcCustomCVForceKernel::execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy) {
    copyState(context, innerContext);
    int numCVs = variableNames.size();
    int numAtoms = cl.getNumAtoms();
    OpenCLContext& cl2 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    vector<double> cvValues;
    vector<map<string, double> > cvDerivs(numCVs);
    for (int i = 0; i < numCVs; i++) {
        cvValues.push_back(innerContext.calcForcesAndEnergy(true, true, 1<<i));
        copyForcesKernel.setArg<cl::Buffer>(0, cvForces[i].getDeviceBuffer());
        cl.executeKernel(copyForcesKernel, numAtoms);
        innerContext.getEnergyParameterDerivatives(cvDerivs[i]);
    }
    
    // Compute the energy and forces.
    
    map<string, double> variables;
    for (auto& name : globalParameterNames)
        variables[name] = context.getParameter(name);
    for (int i = 0; i < numCVs; i++)
        variables[variableNames[i]] = cvValues[i];
    double energy = energyExpression.evaluate(variables);
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate(variables);
        if (cl.getUseDoublePrecision())
            addForcesKernel.setArg<cl_double>(2*i+3, dEdV);
        else
            addForcesKernel.setArg<cl_float>(2*i+3, dEdV);
    }
    cl.executeKernel(addForcesKernel, numAtoms);
    
    // Compute the energy parameter derivatives.
    
    map<string, double>& energyParamDerivs = cl.getEnergyParamDerivWorkspace();
    for (int i = 0; i < paramDerivExpressions.size(); i++)
        energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate(variables);
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate(variables);
        for (auto& deriv : cvDerivs[i])
            energyParamDerivs[deriv.first] += dEdV*deriv.second;
    }
    return energy;
}

void OpenCLCalcCustomCVForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext) {
    int numAtoms = cl.getNumAtoms();
    OpenCLContext& cl2 = *reinterpret_cast<OpenCLPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // Initialize the listeners.
        
        ReorderListener* listener1 = new ReorderListener(cl, invAtomOrder);
        ReorderListener* listener2 = new ReorderListener(cl2, innerInvAtomOrder);
        cl.addReorderListener(listener1);
        cl2.addReorderListener(listener2);
        listener1->execute();
        listener2->execute();
        
        // Initialize the kernels.
        
        copyStateKernel.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        copyStateKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        copyStateKernel.setArg<cl::Buffer>(3, cl.getAtomIndexArray().getDeviceBuffer());
        copyStateKernel.setArg<cl::Buffer>(4, cl2.getPosq().getDeviceBuffer());
        copyStateKernel.setArg<cl::Buffer>(6, cl2.getVelm().getDeviceBuffer());
        copyStateKernel.setArg<cl::Buffer>(7, innerInvAtomOrder.getDeviceBuffer());
        copyStateKernel.setArg<cl_int>(8, numAtoms);
        if (cl.getUseMixedPrecision()) {
            copyStateKernel.setArg<cl::Buffer>(1, cl.getPosqCorrection().getDeviceBuffer());
            copyStateKernel.setArg<cl::Buffer>(5, cl2.getPosqCorrection().getDeviceBuffer());
        }
        else {
            copyStateKernel.setArg<void*>(1, NULL);
            copyStateKernel.setArg<void*>(5, NULL);
        }

        copyForcesKernel.setArg<cl::Buffer>(1, invAtomOrder.getDeviceBuffer());
        copyForcesKernel.setArg<cl::Buffer>(2, cl2.getForce().getDeviceBuffer());
        copyForcesKernel.setArg<cl::Buffer>(3, cl2.getAtomIndexArray().getDeviceBuffer());
        copyForcesKernel.setArg<cl_int>(4, numAtoms);

        addForcesKernel.setArg<cl::Buffer>(0, cl.getForce().getDeviceBuffer());
        addForcesKernel.setArg<cl_int>(1, numAtoms);
        for (int i = 0; i < cvForces.size(); i++)
            addForcesKernel.setArg<cl::Buffer>(2*i+2, cvForces[i].getDeviceBuffer());
    }
    cl.executeKernel(copyStateKernel, numAtoms);
    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext.setPeriodicBoxVectors(a, b, c);
    innerContext.setTime(context.getTime());
    map<string, double> innerParameters = innerContext.getParameters();
    for (auto& param : innerParameters)
        innerContext.setParameter(param.first, context.getParameter(param.first));
}

void OpenCLCalcCustomCVForceKernel::copyParametersToContext(ContextImpl& context, const CustomCVForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, CustomFunction*> functions;
    for (int i = 0; i < (int) force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Replace tabulated functions in the expressions.

    replaceFunctionsInExpression(functions, energyExpression);
    for (auto& expression : variableDerivExpressions)
        replaceFunctionsInExpression(functions, expression);
    for (auto& expression : paramDerivExpressions)
        replaceFunctionsInExpression(functions, expression);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

void OpenCLNoseHooverChainKernel::initialize() {
    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();
    map<string, string> defines;
    defines["BEGIN_YS_LOOP"] = "const real arr[1] = {1.0}; for(int i=0;i<1;++i) { const real ys = arr[i];";
    defines["END_YS_LOOP"] = "}";
    cl::Program program = cl.createProgram(OpenCLKernelSources::noseHooverChain, defines);
    propagateKernels[1] = cl::Kernel(program, "propagateNoseHooverChain");
    defines["BEGIN_YS_LOOP"] = "const real arr[3] = {0.828981543588751, -0.657963087177502, 0.828981543588751}; for(int i=0;i<3;++i) { const real ys = arr[i];";
    program = cl.createProgram(OpenCLKernelSources::noseHooverChain, defines);
    propagateKernels[3] = cl::Kernel(program, "propagateNoseHooverChain");
    defines["BEGIN_YS_LOOP"] = "const real arr[5] = {0.2967324292201065, 0.2967324292201065, -0.186929716880426, 0.2967324292201065, 0.2967324292201065}; for(int i=0;i<5;++i) { const real ys = arr[i];";
    program = cl.createProgram(OpenCLKernelSources::noseHooverChain, defines);
    propagateKernels[5] = cl::Kernel(program, "propagateNoseHooverChain");
    program = cl.createProgram(OpenCLKernelSources::noseHooverChain, defines);
    reduceEnergyKernel = cl::Kernel(program, "reduceEnergyPair");

    computeHeatBathEnergyKernel = cl::Kernel(program, "computeHeatBathEnergy");
    computeAtomsKineticEnergyKernel = cl::Kernel(program, "computeAtomsKineticEnergy");
    computePairsKineticEnergyKernel = cl::Kernel(program, "computePairsKineticEnergy");
    scaleAtomsVelocitiesKernel = cl::Kernel(program, "scaleAtomsVelocities");
    scalePairsVelocitiesKernel = cl::Kernel(program, "scalePairsVelocities");
    int energyBufferSize = cl.getEnergyBuffer().getSize();
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision())
        energyBuffer.initialize<mm_double2>(cl, energyBufferSize, "energyBuffer");
    else
        energyBuffer.initialize<mm_float2>(cl, energyBufferSize, "energyBuffer");
}

std::pair<double, double> OpenCLNoseHooverChainKernel::propagateChain(ContextImpl& context, const NoseHooverChain &nhc, std::pair<double, double> kineticEnergies, double timeStep) {
    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();
    int chainID = nhc.getChainID();
    int nAtoms = nhc.getThermostatedAtoms().size();
    int nPairs = nhc.getThermostatedPairs().size();
    int chainLength = nhc.getChainLength();
    int numYS = nhc.getNumYoshidaSuzukiTimeSteps();
    int numMTS = nhc.getNumMultiTimeSteps();
    int numDOFs = nhc.getNumDegreesOfFreedom();
    double temperature = nhc.getTemperature();
    double frequency = nhc.getCollisionFrequency();
    double relativeTemperature = nhc.getRelativeTemperature();
    double relativeFrequency = nhc.getRelativeCollisionFrequency();

    if (numYS != 1 && numYS != 3 && numYS != 5) {
        throw OpenMMException("Number of Yoshida Suzuki time steps has to be 1, 3, or 5.");
    }

    auto & chainState = cl.getIntegrationUtilities().getNoseHooverChainState();

    if (!scaleFactorBuffer.isInitialized() || scaleFactorBuffer.getSize() == 0) {
        if (useDouble) {
            std::vector<mm_double2> zeros{{0,0}};
            if (scaleFactorBuffer.isInitialized())
                scaleFactorBuffer.resize(1);
            else
                scaleFactorBuffer.initialize<mm_double2>(cl, 1, "scaleFactorBuffer");
            scaleFactorBuffer.upload(zeros);
        }
        else {
            std::vector<mm_float2> zeros{{0,0}};
            if (scaleFactorBuffer.isInitialized())
                scaleFactorBuffer.resize(1);
            else
                scaleFactorBuffer.initialize<mm_float2>(cl, 1, "scaleFactorBuffer");
            scaleFactorBuffer.upload(zeros);
        }
    }
    if (!chainForces.isInitialized() || !chainMasses.isInitialized()) {
        if (useDouble) {
            std::vector<cl_double> zeros(chainLength,0);
            if (chainForces.isInitialized()) {
                chainMasses.resize(chainLength);
                chainForces.resize(chainLength);
            }
            else {
                chainMasses.initialize<cl_double>(cl, chainLength, "chainMasses");
                chainForces.initialize<cl_double>(cl, chainLength, "chainForces");
            }
            chainMasses.upload(zeros);
            chainForces.upload(zeros);
        }
        else {
            std::vector<cl_float> zeros(chainLength,0);
            if (chainForces.isInitialized()) {
                chainMasses.resize(chainLength);
                chainForces.resize(chainLength);
            }
            else {
                chainMasses.initialize<cl_float>(cl, chainLength, "chainMasses");
                chainForces.initialize<cl_float>(cl, chainLength, "chainForces");
            }
            chainMasses.upload(zeros);
            chainForces.upload(zeros);
        }
    }
    if (chainForces.getSize() < chainLength)
        chainMasses.resize(chainLength);
    if (chainMasses.getSize() < chainLength)
        chainMasses.resize(chainLength);

    float timeStepFloat = (float) timeStep;
    // N.B. We ignore the incoming kineticEnergy and grab it from the device buffer instead
    if (nAtoms) {
        if (!chainState.count(2*chainID))
            chainState[2*chainID] = ComputeArray();
        if (!chainState.at(2*chainID).isInitialized() || chainState.at(2*chainID).getSize() != chainLength) {
            // We need to upload the OpenCL array
            if (useDouble) {
                if (chainState.at(2*chainID).isInitialized())
                    chainState.at(2*chainID).resize(chainLength);
                else
                    chainState.at(2*chainID).initialize<mm_double2>(cl, chainLength, "chainState" + std::to_string(2*chainID));
                std::vector<mm_double2> zeros(chainLength, mm_double2(0.0, 0.0));
                chainState.at(2*chainID).upload(zeros.data());
            }
            else {
                if (chainState.at(2*chainID).isInitialized())
                    chainState.at(2*chainID).resize(chainLength);
                else
                    chainState.at(2*chainID).initialize<mm_float2>(cl, chainLength, "chainState" + std::to_string(2*chainID));
                std::vector<mm_float2> zeros(chainLength, mm_float2(0.0f, 0.0f));
                chainState.at(2*chainID).upload(zeros.data());
            }
        }
        int chainType = 0;
        double kT = BOLTZ * temperature;
        float kTfloat = (float) kT;
        float frequencyFloat = (float) frequency;
        propagateKernels[numYS].setArg<cl::Buffer>(0, cl.unwrap(chainState[2*chainID]).getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(1, kineticEnergyBuffer.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(2, scaleFactorBuffer.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(3, chainMasses.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(4, chainForces.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl_int>(5, chainType);
        propagateKernels[numYS].setArg<cl_int>(6, chainLength);
        propagateKernels[numYS].setArg<cl_int>(7, numMTS);
        propagateKernels[numYS].setArg<cl_int>(8, numDOFs);
        propagateKernels[numYS].setArg<cl_float>(9, timeStepFloat);
        if (useDouble) 
            propagateKernels[numYS].setArg<cl_double>(10, kT);
        else
            propagateKernels[numYS].setArg<cl_float>(10, kTfloat);
        propagateKernels[numYS].setArg<cl_float>(11, frequencyFloat);
        cl.executeKernel(propagateKernels[numYS], 1, 1);
    }
    if (nPairs) {
        if (!chainState.count(2*chainID+1))
            chainState[2*chainID+1] = ComputeArray();
        if (!chainState.at(2*chainID+1).isInitialized() || chainState.at(2*chainID+1).getSize() != chainLength) {
            // We need to upload the OpenCL array
            if (useDouble) {
                if (chainState.at(2*chainID+1).isInitialized())
                    chainState.at(2*chainID+1).resize(chainLength);
                else
                    chainState.at(2*chainID+1).initialize<mm_double2>(cl, chainLength, "chainState" + std::to_string(2*chainID+1));
                std::vector<mm_double2> zeros(chainLength, mm_double2(0.0, 0.0));
                chainState.at(2*chainID+1).upload(zeros.data());
            }
            else {
                if (chainState.at(2*chainID+1).isInitialized())
                    chainState.at(2*chainID+1).resize(chainLength);
                else
                    chainState.at(2*chainID+1).initialize<mm_float2>(cl, chainLength, "chainState" + std::to_string(2*chainID+1));
                std::vector<mm_float2> zeros(chainLength, mm_float2(0.0f, 0.0f));
                chainState.at(2*chainID+1).upload(zeros.data());
            }
        }
        int chainType = 1;
        double kT = BOLTZ * relativeTemperature;
        int ndf = 3*nPairs;
        float kTfloat = (float) kT;
        float frequencyFloat = (float) relativeFrequency;
        propagateKernels[numYS].setArg<cl::Buffer>(0, cl.unwrap(chainState[2*chainID+1]).getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(1, kineticEnergyBuffer.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(2, scaleFactorBuffer.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(3, chainMasses.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl::Buffer>(4, chainForces.getDeviceBuffer());
        propagateKernels[numYS].setArg<cl_int>(5, chainType);
        propagateKernels[numYS].setArg<cl_int>(6, chainLength);
        propagateKernels[numYS].setArg<cl_int>(7, numMTS);
        propagateKernels[numYS].setArg<cl_int>(8, ndf);
        propagateKernels[numYS].setArg<cl_float>(9, timeStepFloat);
        if (useDouble) 
            propagateKernels[numYS].setArg<cl_double>(10, kT);
        else
            propagateKernels[numYS].setArg<cl_float>(10, kTfloat);
        propagateKernels[numYS].setArg<cl_float>(11, frequencyFloat);
        cl.executeKernel(propagateKernels[numYS], 1, 1);
    }
    return {0, 0};
}

double OpenCLNoseHooverChainKernel::computeHeatBathEnergy(ContextImpl& context, const NoseHooverChain &nhc) {

    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();

    int chainID = nhc.getChainID();
    int chainLength = nhc.getChainLength();

    auto & chainState = cl.getIntegrationUtilities().getNoseHooverChainState();

    bool absChainIsValid = chainState.count(2*chainID) != 0 &&
                           chainState[2*chainID].isInitialized() &&
                           chainState[2*chainID].getSize() == chainLength;
    bool relChainIsValid = chainState.count(2*chainID+1) != 0 &&
                           chainState[2*chainID+1].isInitialized() &&
                           chainState[2*chainID+1].getSize() == chainLength;

    if (!absChainIsValid && !relChainIsValid) return 0.0;

    if (!heatBathEnergy.isInitialized() || heatBathEnergy.getSize() == 0) {
        if (useDouble) {
            std::vector<cl_double> one(1);
            heatBathEnergy.initialize<cl_double>(cl, 1, "heatBathEnergy");
            heatBathEnergy.upload(one);
        }
        else {
            std::vector<cl_float> one(1);
            heatBathEnergy.initialize<cl_float>(cl, 1, "heatBathEnergy");
            heatBathEnergy.upload(one);
        }
    }

    cl.clearBuffer(heatBathEnergy);

    if (absChainIsValid) {
        int numDOFs = nhc.getNumDegreesOfFreedom();
        double temperature = nhc.getTemperature();
        double frequency = nhc.getCollisionFrequency();
        double kT = BOLTZ * temperature;
        float kTfloat = (float) kT;
        float frequencyFloat = (float) frequency;

        computeHeatBathEnergyKernel.setArg<cl::Buffer>(0, heatBathEnergy.getDeviceBuffer());
        computeHeatBathEnergyKernel.setArg<cl_int>(1, chainLength);
        computeHeatBathEnergyKernel.setArg<cl_int>(2, numDOFs); 
        if (useDouble)
            computeHeatBathEnergyKernel.setArg<cl_double>(3, kT);
        else
            computeHeatBathEnergyKernel.setArg<cl_float>(3, kTfloat);
        computeHeatBathEnergyKernel.setArg<cl_float>(4, frequencyFloat);
        computeHeatBathEnergyKernel.setArg<cl::Buffer>(5, cl.unwrap(chainState[2*chainID]).getDeviceBuffer());
        cl.executeKernel(computeHeatBathEnergyKernel, 1, 1);
    }
    if (relChainIsValid) {
        int numDOFs = 3 * nhc.getThermostatedPairs().size();
        double temperature = nhc.getRelativeTemperature();
        double frequency = nhc.getRelativeCollisionFrequency();
        double kT = BOLTZ * temperature;
        float kTfloat = (float) kT;
        float frequencyFloat = (float) frequency;

        computeHeatBathEnergyKernel.setArg<cl::Buffer>(0, heatBathEnergy.getDeviceBuffer());
        computeHeatBathEnergyKernel.setArg<cl_int>(1, chainLength);
        computeHeatBathEnergyKernel.setArg<cl_int>(2, numDOFs); 
        if (useDouble)
            computeHeatBathEnergyKernel.setArg<cl_double>(3, kT);
        else
            computeHeatBathEnergyKernel.setArg<cl_float>(3, kTfloat);
        computeHeatBathEnergyKernel.setArg<cl_float>(4, frequencyFloat);
        computeHeatBathEnergyKernel.setArg<cl::Buffer>(5, cl.unwrap(chainState[2*chainID+1]).getDeviceBuffer());
        cl.executeKernel(computeHeatBathEnergyKernel, 1, 1);
    }


    void * pinnedBuffer = cl.getPinnedBuffer();
    heatBathEnergy.download(pinnedBuffer);
    if (useDouble)
        return *((double*) pinnedBuffer);
    else
        return *((float*) pinnedBuffer);
}

std::pair<double, double> OpenCLNoseHooverChainKernel::computeMaskedKineticEnergy(ContextImpl& context, const NoseHooverChain &nhc, bool downloadValue) {

    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();

    int chainID = nhc.getChainID();
    const auto & nhcAtoms = nhc.getThermostatedAtoms();
    const auto & nhcPairs = nhc.getThermostatedPairs();
    auto nAtoms = nhcAtoms.size();
    auto nPairs = nhcPairs.size();
    if (nAtoms) {
        if (!atomlists.count(chainID)) { 
            // We need to upload the OpenCL array
            atomlists[chainID] = OpenCLArray();
            atomlists[chainID].initialize<int>(cl, nAtoms, "atomlist" + std::to_string(chainID));
            atomlists[chainID].upload(nhcAtoms);
        }
        if (atomlists[chainID].getSize() != nAtoms) {
            throw OpenMMException("Number of atoms changed. Cannot be handled by the same Nose-Hoover thermostat.");
        }
    }
    if (nPairs) {
        if (!pairlists.count(chainID)) { 
            // We need to upload the OpenCL array
            pairlists[chainID] = OpenCLArray();
            pairlists[chainID].initialize<mm_int2>(cl, nPairs, "pairlist" + std::to_string(chainID));
            std::vector<mm_int2> int2vec;
            for(const auto &p : nhcPairs) int2vec.push_back(mm_int2(p.first, p.second));
            pairlists[chainID].upload(int2vec);
        }
        if (pairlists[chainID].getSize() != nPairs) {
            throw OpenMMException("Number of thermostated pairs changed. Cannot be handled by the same Nose-Hoover thermostat.");
        }
    }
    if (!kineticEnergyBuffer.isInitialized() || kineticEnergyBuffer.getSize() == 0) {
        if (useDouble) {
            std::vector<mm_double2> zeros{{0,0}};
            kineticEnergyBuffer.initialize<mm_double2>(cl, 1, "kineticEnergyBuffer");
            kineticEnergyBuffer.upload(zeros);
        }
        else {
            std::vector<mm_float2> zeros{{0,0}};
            kineticEnergyBuffer.initialize<mm_float2>(cl, 1, "kineticEnergyBuffer");
            kineticEnergyBuffer.upload(zeros);
        }
    }
    cl.clearBuffer(cl.getEnergyBuffer());
    if (nAtoms) {
        computeAtomsKineticEnergyKernel.setArg<cl::Buffer>(0, energyBuffer.getDeviceBuffer());
        computeAtomsKineticEnergyKernel.setArg<cl_int>(1, nAtoms);
        computeAtomsKineticEnergyKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        computeAtomsKineticEnergyKernel.setArg<cl::Buffer>(3, atomlists[chainID].getDeviceBuffer());
        cl.executeKernel(computeAtomsKineticEnergyKernel, nAtoms);
    }
    if (nPairs) {
        computePairsKineticEnergyKernel.setArg<cl::Buffer>(0, energyBuffer.getDeviceBuffer());
        computePairsKineticEnergyKernel.setArg<cl_int>(1, nPairs);
        computePairsKineticEnergyKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        computePairsKineticEnergyKernel.setArg<cl::Buffer>(3, pairlists[chainID].getDeviceBuffer());
        cl.executeKernel(computePairsKineticEnergyKernel, nPairs);
    }
    int bufferSize = energyBuffer.getSize();
    int workGroupSize  = cl.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
    if (workGroupSize > 512)
        workGroupSize = 512;
    reduceEnergyKernel.setArg<cl::Buffer>(0, energyBuffer.getDeviceBuffer());
    reduceEnergyKernel.setArg<cl::Buffer>(1, kineticEnergyBuffer.getDeviceBuffer());
    reduceEnergyKernel.setArg<cl_int>(2, bufferSize);
    reduceEnergyKernel.setArg<cl_int>(3, workGroupSize);
    reduceEnergyKernel.setArg(4, workGroupSize*energyBuffer.getElementSize(), NULL);
    cl.executeKernel(reduceEnergyKernel, workGroupSize, workGroupSize);

    std::pair<double, double> KEs = {0, 0};
    if (downloadValue) {
        if (useDouble) {
            mm_double2 tmp;
            kineticEnergyBuffer.download(&tmp);
            KEs.first = tmp.x;
            KEs.second = tmp.y;
        }
        else {
            mm_float2 tmp;
            kineticEnergyBuffer.download(&tmp);
            KEs.first = tmp.x;
            KEs.second = tmp.y;
        }
    }
    return KEs;
}

void OpenCLNoseHooverChainKernel::scaleVelocities(ContextImpl& context, const NoseHooverChain &nhc, std::pair<double, double> scaleFactor) {
    // For now we assume that the atoms and pairs info is valid, because compute{Atoms|Pairs}KineticEnergy must have been
    // called before this kernel.  If that ever ceases to be true, some sanity checks are needed here.

    int chainID = nhc.getChainID();
    auto nAtoms = nhc.getThermostatedAtoms().size();
    auto nPairs = nhc.getThermostatedPairs().size();
    if (nAtoms) {
        scaleAtomsVelocitiesKernel.setArg<cl::Buffer>(0, scaleFactorBuffer.getDeviceBuffer());
        scaleAtomsVelocitiesKernel.setArg<cl_int>(1, nAtoms);
        scaleAtomsVelocitiesKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        scaleAtomsVelocitiesKernel.setArg<cl::Buffer>(3, atomlists[chainID].getDeviceBuffer());
        cl.executeKernel(scaleAtomsVelocitiesKernel, nAtoms);
    }
    if (nPairs) {
        scalePairsVelocitiesKernel.setArg<cl::Buffer>(0, scaleFactorBuffer.getDeviceBuffer());
        scalePairsVelocitiesKernel.setArg<cl_int>(1, nPairs);
        scalePairsVelocitiesKernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        scalePairsVelocitiesKernel.setArg<cl::Buffer>(3, pairlists[chainID].getDeviceBuffer());
        cl.executeKernel(scalePairsVelocitiesKernel, nPairs);
    }
}

void OpenCLApplyMonteCarloBarostatKernel::initialize(const System& system, const Force& thermostat) {
    savedPositions.initialize(cl, cl.getPaddedNumAtoms(), cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4), "savedPositions");
    savedForces.initialize(cl, cl.getPaddedNumAtoms(), cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4), "savedForces");
    cl::Program program = cl.createProgram(OpenCLKernelSources::monteCarloBarostat);
    kernel = cl::Kernel(program, "scalePositions");
}

void OpenCLApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;

        // Create the arrays with the molecule definitions.

        vector<vector<int> > molecules = context.getMolecules();
        numMolecules = molecules.size();
        moleculeAtoms.initialize<int>(cl, cl.getNumAtoms(), "moleculeAtoms");
        moleculeStartIndex.initialize<int>(cl, numMolecules+1, "moleculeStartIndex");
        vector<int> atoms(moleculeAtoms.getSize());
        vector<int> startIndex(moleculeStartIndex.getSize());
        int index = 0;
        for (int i = 0; i < numMolecules; i++) {
            startIndex[i] = index;
            for (int molecule : molecules[i])
                atoms[index++] = molecule;
        }
        startIndex[numMolecules] = index;
        moleculeAtoms.upload(atoms);
        moleculeStartIndex.upload(startIndex);

        // Initialize the kernel arguments.
        
        kernel.setArg<cl_int>(3, numMolecules);
        kernel.setArg<cl::Buffer>(9, cl.getPosq().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(10, moleculeAtoms.getDeviceBuffer());
        kernel.setArg<cl::Buffer>(11, moleculeStartIndex.getDeviceBuffer());
    }
    int bytesToCopy = cl.getPosq().getSize()*(cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
    cl.getQueue().enqueueCopyBuffer(cl.getPosq().getDeviceBuffer(), savedPositions.getDeviceBuffer(), 0, 0, bytesToCopy);
    cl.getQueue().enqueueCopyBuffer(cl.getForce().getDeviceBuffer(), savedForces.getDeviceBuffer(), 0, 0, bytesToCopy);
    kernel.setArg<cl_float>(0, (cl_float) scaleX);
    kernel.setArg<cl_float>(1, (cl_float) scaleY);
    kernel.setArg<cl_float>(2, (cl_float) scaleZ);
    setPeriodicBoxArgs(cl, kernel, 4);
    cl.executeKernel(kernel, cl.getNumAtoms());
    for (auto& offset : cl.getPosCellOffsets())
        offset = mm_int4(0, 0, 0, 0);
    lastAtomOrder = cl.getAtomIndex();
}

void OpenCLApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    int bytesToCopy = cl.getPosq().getSize()*(cl.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4));
    cl.getQueue().enqueueCopyBuffer(savedPositions.getDeviceBuffer(), cl.getPosq().getDeviceBuffer(), 0, 0, bytesToCopy);
    cl.getQueue().enqueueCopyBuffer(savedForces.getDeviceBuffer(), cl.getForce().getDeviceBuffer(), 0, 0, bytesToCopy);
}
