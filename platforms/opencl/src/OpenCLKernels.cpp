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
#include "jama_eig.h"
#include <algorithm>
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
    int version = 2;
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
    if (version != 2)
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

class OpenCLCalcCustomHbondForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const CustomHbondForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumDonors()+force.getNumAcceptors()+force.getNumExclusions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
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

    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumDonors()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumDonors()/numContexts;
    numDonors = endIndex-startIndex;
    numAcceptors = force.getNumAcceptors();
    if (numDonors == 0 || numAcceptors == 0)
        return;
    int numParticles = system.getNumParticles();
    donors.initialize<mm_int4>(cl, numDonors, "customHbondDonors");
    acceptors.initialize<mm_int4>(cl, numAcceptors, "customHbondAcceptors");
    donorParams = new OpenCLParameterSet(cl, force.getNumPerDonorParameters(), numDonors, "customHbondDonorParameters");
    acceptorParams = new OpenCLParameterSet(cl, force.getNumPerAcceptorParameters(), numAcceptors, "customHbondAcceptorParameters");
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customHbondGlobals", CL_MEM_READ_ONLY);
    vector<vector<cl_float> > donorParamVector(numDonors);
    vector<mm_int4> donorVector(numDonors);
    for (int i = 0; i < numDonors; i++) {
        vector<double> parameters;
        force.getDonorParameters(startIndex+i, donorVector[i].x, donorVector[i].y, donorVector[i].z, parameters);
        donorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            donorParamVector[i][j] = (cl_float) parameters[j];
    }
    donors.upload(donorVector);
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
    acceptors.upload(acceptorVector);
    acceptorParams->setParameterValues(acceptorParamVector);

    // Select an output buffer index for each donor and acceptor.

    donorBufferIndices.initialize<mm_int4>(cl, numDonors, "customHbondDonorBuffers");
    acceptorBufferIndices.initialize<mm_int4>(cl, numAcceptors, "customHbondAcceptorBuffers");
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
    donorBufferIndices.upload(donorBufferVector);
    acceptorBufferIndices.upload(acceptorBufferVector);
    int maxBuffers = 1;
    for (int i : donorBufferCounter)
        maxBuffers = max(maxBuffers, i);
    for (int i : acceptorBufferCounter)
        maxBuffers = max(maxBuffers, i);
    info = new ForceInfo(maxBuffers, force);
    cl.addForce(info);

    // Record exclusions.

    vector<mm_int4> donorExclusionVector(numDonors, mm_int4(-1, -1, -1, -1));
    vector<mm_int4> acceptorExclusionVector(numAcceptors, mm_int4(-1, -1, -1, -1));
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int donor, acceptor;
        force.getExclusionParticles(i, donor, acceptor);
        if (donor < startIndex || donor >= endIndex)
            continue;
        donor -= startIndex;
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
    donorExclusions.initialize<mm_int4>(cl, numDonors, "customHbondDonorExclusions");
    acceptorExclusions.initialize<mm_int4>(cl, numAcceptors, "customHbondAcceptorExclusions");
    donorExclusions.upload(donorExclusionVector);
    acceptorExclusions.upload(acceptorExclusionVector);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumFunctions());
    for (int i = 0; i < force.getNumFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cl.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cl.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cl.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cl, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        tableArgs << ", __global const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* restrict " << arrayName;
    }

    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized())
        globals.upload(globalParamValues);
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
        variables[name] = "globals["+cl.intToString(i)+"]";
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
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real r_"+deltaName+" = SQRT(delta"+deltaName+".w);\n");
        variables[distance.first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cl.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[0]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[2]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName2);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real "+angleName+" = computeAngle(delta"+deltaName1+", delta"+deltaName2+");\n");
        variables[angle.first] = angleName;
        forceExpressions["real dEdAngle"+cl.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
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
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName3+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[3]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName3);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 "+crossName1+" = computeCross(delta"+deltaName1+", delta"+deltaName2+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 "+crossName2+" = computeCross(delta"+deltaName2+", delta"+deltaName3+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real "+dihedralName+" = computeAngle("+crossName1+", "+crossName2+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, dihedralName+" *= (delta"+deltaName1+".x*"+crossName2+".x + delta"+deltaName1+".y*"+crossName2+".y + delta"+deltaName1+".z*"+crossName2+".z < 0 ? -1 : 1);\n");
        variables[dihedral.first] = dihedralName;
        forceExpressions["real dEdDihedral"+cl.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Next it needs to load parameters from global memory.

    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", __global const float* restrict globals";
    for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
        extraArgs << ", __global const "+buffer.getType()+"* restrict donor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" donorParams"+cl.intToString(i+1)+" = donor"+buffer.getName()+"[donorIndex];\n");
    }
    for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
        extraArgs << ", __global const "+buffer.getType()+"* restrict acceptor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" acceptorParams"+cl.intToString(i+1)+" = acceptor"+buffer.getName()+"[acceptorIndex];\n");
    }

    // Now evaluate the expressions.

    computeAcceptor << cl.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");
    forceExpressions["energy += "] = energyExpression;
    computeDonor << cl.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cl.intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "-"+value);
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], value);
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 crossProd = cross(delta"+deltaName2+", delta"+deltaName1+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real lengthCross = max(length(crossProd), (real) 1e-6f);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 deltaCross0 = -cross(delta"+deltaName1+", crossProd)*dEdAngle"+cl.intToString(index)+"/(delta"+deltaName1+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 deltaCross2 = cross(delta"+deltaName2+", crossProd)*dEdAngle"+cl.intToString(index)+"/(delta"+deltaName2+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 deltaCross1 = -(deltaCross0+deltaCross2);\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "deltaCross0.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "deltaCross1.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "deltaCross2.xyz");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
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
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real r = SQRT(delta"+deltaName2+".w);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 ff;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.x = (-dEdDihedral"+cl.intToString(index)+"*r)/"+crossName1+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.y = (delta"+deltaName1+".x*delta"+deltaName2+".x + delta"+deltaName1+".y*delta"+deltaName2+".y + delta"+deltaName1+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.z = (delta"+deltaName3+".x*delta"+deltaName2+".x + delta"+deltaName3+".y*delta"+deltaName2+".y + delta"+deltaName3+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.w = (dEdDihedral"+cl.intToString(index)+"*r)/"+crossName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 internalF0 = ff.x*"+crossName1+";\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 internalF3 = ff.w*"+crossName2+";\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 s = ff.y*internalF0 - ff.z*internalF3;\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "internalF0.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "s.xyz-internalF0.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "-s.xyz-internalF3.xyz");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[3], "internalF3.xyz");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
        index++;
    }

    // Generate the kernels.

    map<string, string> replacements;
    replacements["COMPUTE_DONOR_FORCE"] = computeDonor.str();
    replacements["COMPUTE_ACCEPTOR_FORCE"] = computeAcceptor.str();
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    defines["NUM_DONORS"] = cl.intToString(numDonors);
    defines["NUM_ACCEPTORS"] = cl.intToString(numAcceptors);
    defines["PI"] = cl.doubleToString(M_PI);
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff) {
        defines["USE_CUTOFF"] = "1";
        defines["CUTOFF_SQUARED"] = cl.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff && force.getNonbondedMethod() != CustomHbondForce::CutoffNonPeriodic)
        defines["USE_PERIODIC"] = "1";
    if (force.getNumExclusions() > 0)
        defines["USE_EXCLUSIONS"] = "1";
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customHbondForce, replacements), defines);
    donorKernel = cl::Kernel(program, "computeDonorForces");
    acceptorKernel = cl::Kernel(program, "computeAcceptorForces");
}

double OpenCLCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (numDonors == 0 || numAcceptors == 0)
        return 0.0;
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals.upload(globalParamValues);
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        int index = 0;
        donorKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, donorExclusions.getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, donors.getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, acceptors.getDeviceBuffer());
        donorKernel.setArg<cl::Buffer>(index++, donorBufferIndices.getDeviceBuffer());
        donorKernel.setArg(index++, 3*OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
        index += 5; // Periodic box size arguments are set when the kernel is executed.
        if (globals.isInitialized())
            donorKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (auto& buffer : donorParams->getBuffers())
            donorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& buffer : acceptorParams->getBuffers())
            donorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& function : tabulatedFunctions)
            donorKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        index = 0;
        acceptorKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, acceptorExclusions.getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, donors.getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, acceptors.getDeviceBuffer());
        acceptorKernel.setArg<cl::Buffer>(index++, acceptorBufferIndices.getDeviceBuffer());
        acceptorKernel.setArg(index++, 3*OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
        index += 5; // Periodic box size arguments are set when the kernel is executed.
        if (globals.isInitialized())
            acceptorKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (auto& buffer : donorParams->getBuffers())
            acceptorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& buffer : acceptorParams->getBuffers())
            acceptorKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& function : tabulatedFunctions)
            acceptorKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
    }
    setPeriodicBoxArgs(cl, donorKernel, 8);
    cl.executeKernel(donorKernel, max(numDonors, numAcceptors));
    setPeriodicBoxArgs(cl, acceptorKernel, 8);
    cl.executeKernel(acceptorKernel, max(numDonors, numAcceptors));
    return 0.0;
}

void OpenCLCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumDonors()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumDonors()/numContexts;
    if (numDonors != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of donors has changed");
    if (numAcceptors != force.getNumAcceptors())
        throw OpenMMException("updateParametersInContext: The number of acceptors has changed");
    
    // Record the per-donor parameters.
    
    if (numDonors > 0) {
        vector<vector<cl_float> > donorParamVector(numDonors);
        vector<double> parameters;
        for (int i = 0; i < numDonors; i++) {
            int d1, d2, d3;
            force.getDonorParameters(startIndex+i, d1, d2, d3, parameters);
            donorParamVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                donorParamVector[i][j] = (cl_float) parameters[j];
        }
        donorParams->setParameterValues(donorParamVector);
    }
    
    // Record the per-acceptor parameters.
    
    if (numAcceptors > 0) {
        vector<vector<cl_float> > acceptorParamVector(numAcceptors);
        vector<double> parameters;
        for (int i = 0; i < numAcceptors; i++) {
            int a1, a2, a3;
            force.getAcceptorParameters(i, a1, a2, a3, parameters);
            acceptorParamVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                acceptorParamVector[i][j] = (cl_float) parameters[j];
        }
        acceptorParams->setParameterValues(acceptorParamVector);
    }
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
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

class OpenCLCalcRMSDForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const RMSDForce& force) : OpenCLForceInfo(0), force(force) {
        updateParticles();
    }
    void updateParticles() {
        particles.clear();
        for (int i : force.getParticles())
            particles.insert(i);
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        bool include1 = (particles.find(particle1) != particles.end());
        bool include2 = (particles.find(particle2) != particles.end());
        return (include1 == include2);
    }
private:
    const RMSDForce& force;
    set<int> particles;
};

void OpenCLCalcRMSDForceKernel::initialize(const System& system, const RMSDForce& force) {
    // Create data structures.
    
    bool useDouble = cl.getUseDoublePrecision();
    int elementSize = (useDouble ? sizeof(cl_double) : sizeof(cl_float));
    int numParticles = force.getParticles().size();
    if (numParticles == 0)
        numParticles = system.getNumParticles();
    referencePos.initialize(cl, system.getNumParticles(), 4*elementSize, "referencePos");
    particles.initialize<cl_int>(cl, numParticles, "particles");
    buffer.initialize(cl, 13, elementSize, "buffer");
    recordParameters(force);
    info = new ForceInfo(force);
    cl.addForce(info);
    
    // Create the kernels.

    cl::Program program = cl.createProgram(OpenCLKernelSources::rmsd);
    kernel1 = cl::Kernel(program, "computeRMSDPart1");
    kernel2 = cl::Kernel(program, "computeRMSDForces");
}

void OpenCLCalcRMSDForceKernel::recordParameters(const RMSDForce& force) {
    // Record the parameters and center the reference positions.
    
    vector<int> particleVec = force.getParticles();
    if (particleVec.size() == 0)
        for (int i = 0; i < cl.getNumAtoms(); i++)
            particleVec.push_back(i);
    vector<Vec3> centeredPositions = force.getReferencePositions();
    Vec3 center;
    for (int i : particleVec)
        center += centeredPositions[i];
    center /= particleVec.size();
    for (Vec3& p : centeredPositions)
        p -= center;

    // Upload them to the device.

    particles.upload(particleVec);
    vector<mm_double4> pos;
    for (Vec3 p : centeredPositions)
        pos.push_back(mm_double4(p[0], p[1], p[2], 0));
    referencePos.upload(pos, true);

    // Record the sum of the norms of the reference positions.

    sumNormRef = 0.0;
    for (int i : particleVec) {
        Vec3 p = centeredPositions[i];
        sumNormRef += p.dot(p);
    }
}

double OpenCLCalcRMSDForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (cl.getUseDoublePrecision())
        return executeImpl<double>(context);
    return executeImpl<float>(context);
}

template <class REAL>
double OpenCLCalcRMSDForceKernel::executeImpl(ContextImpl& context) {
    // Execute the first kernel.

    int numParticles = particles.getSize();
    int blockSize = min(256, (int) kernel1.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(cl.getDevice()));
    kernel1.setArg<cl_int>(0, numParticles);
    kernel1.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(2, referencePos.getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(3, particles.getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(4, buffer.getDeviceBuffer());
    kernel1.setArg(5, blockSize*sizeof(REAL), NULL);
    cl.executeKernel(kernel1, blockSize, blockSize);
    
    // Download the results, build the F matrix, and find the maximum eigenvalue
    // and eigenvector.

    vector<REAL> b;
    buffer.download(b);
    Array2D<double> F(4, 4);
    F[0][0] =  b[0*3+0] + b[1*3+1] + b[2*3+2];
    F[1][0] =  b[1*3+2] - b[2*3+1];
    F[2][0] =  b[2*3+0] - b[0*3+2];
    F[3][0] =  b[0*3+1] - b[1*3+0];
    F[0][1] =  b[1*3+2] - b[2*3+1];
    F[1][1] =  b[0*3+0] - b[1*3+1] - b[2*3+2];
    F[2][1] =  b[0*3+1] + b[1*3+0];
    F[3][1] =  b[0*3+2] + b[2*3+0];
    F[0][2] =  b[2*3+0] - b[0*3+2];
    F[1][2] =  b[0*3+1] + b[1*3+0];
    F[2][2] = -b[0*3+0] + b[1*3+1] - b[2*3+2];
    F[3][2] =  b[1*3+2] + b[2*3+1];
    F[0][3] =  b[0*3+1] - b[1*3+0];
    F[1][3] =  b[0*3+2] + b[2*3+0];
    F[2][3] =  b[1*3+2] + b[2*3+1];
    F[3][3] = -b[0*3+0] - b[1*3+1] + b[2*3+2];
    JAMA::Eigenvalue<double> eigen(F);
    Array1D<double> values;
    eigen.getRealEigenvalues(values);
    Array2D<double> vectors;
    eigen.getV(vectors);

    // Compute the RMSD.

    double msd = (sumNormRef+b[9]-2*values[3])/numParticles;
    if (msd < 1e-20) {
        // The particles are perfectly aligned, so all the forces should be zero.
        // Numerical error can lead to NaNs, so just return 0 now.
        return 0.0;
    }
    double rmsd = sqrt(msd);
    b[9] = rmsd;

    // Compute the rotation matrix.

    double q[] = {vectors[0][3], vectors[1][3], vectors[2][3], vectors[3][3]};
    double q00 = q[0]*q[0], q01 = q[0]*q[1], q02 = q[0]*q[2], q03 = q[0]*q[3];
    double q11 = q[1]*q[1], q12 = q[1]*q[2], q13 = q[1]*q[3];
    double q22 = q[2]*q[2], q23 = q[2]*q[3];
    double q33 = q[3]*q[3];
    b[0] = q00+q11-q22-q33;
    b[1] = 2*(q12-q03);
    b[2] = 2*(q13+q02);
    b[3] = 2*(q12+q03);
    b[4] = q00-q11+q22-q33;
    b[5] = 2*(q23-q01);
    b[6] = 2*(q13-q02);
    b[7] = 2*(q23+q01);
    b[8] = q00-q11-q22+q33;

    // Upload it to the device and invoke the kernel to apply forces.
    
    buffer.upload(b);
    kernel2.setArg<cl_int>(0, numParticles);
    kernel2.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(2, referencePos.getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(3, particles.getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(4, buffer.getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(5, cl.getForceBuffers().getDeviceBuffer());
    cl.executeKernel(kernel2, numParticles);
    return rmsd;
}

void OpenCLCalcRMSDForceKernel::copyParametersToContext(ContextImpl& context, const RMSDForce& force) {
    if (referencePos.getSize() != force.getReferencePositions().size())
        throw OpenMMException("updateParametersInContext: The number of reference positions has changed");
    int numParticles = force.getParticles().size();
    if (numParticles == 0)
        numParticles = context.getSystem().getNumParticles();
    if (numParticles != particles.getSize())
        particles.resize(numParticles);
    recordParameters(force);
    
    // Mark that the current reordering may be invalid.
    
    info->updateParticles();
    cl.invalidateMolecules(info);
}

void OpenCLApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    randomSeed = thermostat.getRandomNumberSeed();
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::andersenThermostat, defines);
    kernel = cl::Kernel(program, "applyAndersenThermostat");
    cl.getIntegrationUtilities().initRandomNumberGenerator(randomSeed);

    // Create the arrays with the group definitions.

    vector<vector<int> > groups = AndersenThermostatImpl::calcParticleGroups(system);
    atomGroups.initialize<int>(cl, cl.getNumAtoms(), "atomGroups");
    vector<int> atoms(atomGroups.getSize());
    for (int i = 0; i < (int) groups.size(); i++) {
        for (int j = 0; j < (int) groups[i].size(); j++)
            atoms[groups[i][j]] = i;
    }
    atomGroups.upload(atoms);
}

void OpenCLApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel.setArg<cl::Buffer>(2, cl.getVelm().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(3, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(4, cl.getIntegrationUtilities().getRandom().getDeviceBuffer());
        kernel.setArg<cl::Buffer>(6, atomGroups.getDeviceBuffer());
    }
    kernel.setArg<cl_float>(0, (cl_float) context.getParameter(AndersenThermostat::CollisionFrequency()));
    kernel.setArg<cl_float>(1, (cl_float) (BOLTZ*context.getParameter(AndersenThermostat::Temperature())));
    kernel.setArg<cl_uint>(5, cl.getIntegrationUtilities().prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel, cl.getNumAtoms());
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
