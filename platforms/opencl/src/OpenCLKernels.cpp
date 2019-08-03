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
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCentroidBondForceImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/internal/CustomHbondForceImpl.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
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

class OpenCLCalcHarmonicBondForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const HarmonicBondForce& force) : OpenCLForceInfo(0), force(force) {
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

void OpenCLCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params.initialize<mm_float2>(cl, numBonds, "bondParams");
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = OpenCLKernelSources::harmonicBondForce;
    replacements["PARAMS"] = cl.getBondedUtilities().addArgument(params.getDeviceBuffer(), "float2");
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::bondForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cl.addForce(info);
}

double OpenCLCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
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
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCustomBondForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomBondForce& force) : OpenCLForceInfo(0), force(force) {
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

OpenCLCalcCustomBondForceKernel::~OpenCLCalcCustomBondForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = new OpenCLParameterSet(cl, force.getNumPerBondParameters(), numBonds, "customBondParams");
    vector<vector<cl_float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cl.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
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
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customBondGlobals", CL_MEM_READ_ONLY);
        globals.upload(globalParamValues);
        string argName = cl.getBondedUtilities().addArgument(globals.getDeviceBuffer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cl.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cl.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cl.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cl.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::bondForce, replacements), force.getForceGroup());
}

double OpenCLCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    return 0.0;
}

void OpenCLCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<vector<cl_float> > paramVector(numBonds);
    vector<double> parameters;
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        force.getBondParameters(startIndex+i, atom1, atom2, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcHarmonicAngleForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const HarmonicAngleForce& force) : OpenCLForceInfo(0), force(force) {
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

void OpenCLCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params.initialize<mm_float2>(cl, numAngles, "angleParams");
    vector<mm_float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], angle, k);
        paramVector[i] = mm_float2((cl_float) angle, (cl_float) k);

    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = OpenCLKernelSources::harmonicAngleForce;
    replacements["PARAMS"] = cl.getBondedUtilities().addArgument(params.getDeviceBuffer(), "float2");
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::angleForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cl.addForce(info);
}

double OpenCLCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumAngles()/numContexts;
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
        paramVector[i] = mm_float2((cl_float) angle, (cl_float) k);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCustomAngleForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomAngleForce& force) : OpenCLForceInfo(0), force(force) {
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

OpenCLCalcCustomAngleForceKernel::~OpenCLCalcCustomAngleForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = new OpenCLParameterSet(cl, force.getNumPerAngleParameters(), numAngles, "customAngleParams");
    vector<vector<cl_float> > paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        vector<double> parameters;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cl.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
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
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customAngleGlobals", CL_MEM_READ_ONLY);
        globals.upload(globalParamValues);
        string argName = cl.getBondedUtilities().addArgument(globals.getDeviceBuffer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cl.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cl.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cl.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" angleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cl.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::angleForce, replacements), force.getForceGroup());
}

double OpenCLCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    return 0.0;
}

void OpenCLCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0)
        return;
    
    // Record the per-angle parameters.
    
    vector<vector<cl_float> > paramVector(numAngles);
    vector<double> parameters;
    for (int i = 0; i < numAngles; i++) {
        int atom1, atom2, atom3;
        force.getAngleParameters(startIndex+i, atom1, atom2, atom3, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcPeriodicTorsionForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const PeriodicTorsionForce& force) : OpenCLForceInfo(0), force(force) {
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

void OpenCLCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params.initialize<mm_float4>(cl, numTorsions, "periodicTorsionParams");
    vector<mm_float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], periodicity, phase, k);
        paramVector[i] = mm_float4((cl_float) k, (cl_float) phase, (cl_float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = OpenCLKernelSources::periodicTorsionForce;
    replacements["PARAMS"] = cl.getBondedUtilities().addArgument(params.getDeviceBuffer(), "float4");
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cl.addForce(info);
}

double OpenCLCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
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
        paramVector[i] = mm_float4((cl_float) k, (cl_float) phase, (cl_float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcRBTorsionForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const RBTorsionForce& force) : OpenCLForceInfo(0), force(force) {
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

void OpenCLCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params.initialize<mm_float8>(cl, numTorsions, "rbTorsionParams");
    vector<mm_float8> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], c0, c1, c2, c3, c4, c5);
        paramVector[i] = mm_float8((cl_float) c0, (cl_float) c1, (cl_float) c2, (cl_float) c3, (cl_float) c4, (cl_float) c5, 0.0f, 0.0f);

    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = OpenCLKernelSources::rbTorsionForce;
    replacements["PARAMS"] = cl.getBondedUtilities().addArgument(params.getDeviceBuffer(), "float8");
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cl.addForce(info);
}

double OpenCLCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<mm_float8> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        paramVector[i] = mm_float8((cl_float) c0, (cl_float) c1, (cl_float) c2, (cl_float) c3, (cl_float) c4, (cl_float) c5, 0.0f, 0.0f);
    }
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCMAPTorsionForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CMAPTorsionForce& force) : OpenCLForceInfo(0), force(force) {
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

void OpenCLCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
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
    vector<cl_int> torsionMapsVec(numTorsions);
    for (int i = 0; i < numTorsions; i++)
        force.getTorsionParameters(startIndex+i, torsionMapsVec[i], atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], atoms[i][6], atoms[i][7]);
    coefficients.initialize<mm_float4>(cl, coeffVec.size(), "cmapTorsionCoefficients");
    mapPositions.initialize<mm_int2>(cl, numMaps, "cmapTorsionMapPositions");
    torsionMaps.initialize<cl_int>(cl, numTorsions, "cmapTorsionMaps");
    coefficients.upload(coeffVec);
    mapPositions.upload(mapPositionsVec);
    torsionMaps.upload(torsionMapsVec);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COEFF"] = cl.getBondedUtilities().addArgument(coefficients.getDeviceBuffer(), "float4");
    replacements["MAP_POS"] = cl.getBondedUtilities().addArgument(mapPositions.getDeviceBuffer(), "int2");
    replacements["MAPS"] = cl.getBondedUtilities().addArgument(torsionMaps.getDeviceBuffer(), "int");
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::cmapTorsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cl.addForce(info);
}

double OpenCLCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    int numMaps = force.getNumMaps();
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
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

class OpenCLCalcCustomTorsionForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomTorsionForce& force) : OpenCLForceInfo(0), force(force) {
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

OpenCLCalcCustomTorsionForceKernel::~OpenCLCalcCustomTorsionForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params = new OpenCLParameterSet(cl, force.getNumPerTorsionParameters(), numTorsions, "customTorsionParams");
    vector<vector<cl_float> > paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        vector<double> parameters;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cl.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
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
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customTorsionGlobals", CL_MEM_READ_ONLY);
        globals.upload(globalParamValues);
        string argName = cl.getBondedUtilities().addArgument(globals.getDeviceBuffer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cl.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cl.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cl.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" torsionParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cl.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::torsionForce, replacements), force.getForceGroup());
}

double OpenCLCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    return 0.0;
}

void OpenCLCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<vector<cl_float> > paramVector(numTorsions);
    vector<double> parameters;
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
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
                            xmoduli->upload(moduli, true, true);
                        else if (dim == 1)
                            ymoduli->upload(moduli, true, true);
                        else
                            zmoduli->upload(moduli, true, true);
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
        globalParams.upload(paramValues, true, true);
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

class OpenCLCalcCustomNonbondedForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const CustomNonbondedForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
        if (force.getNumInteractionGroups() > 0) {
            groupsForParticle.resize(force.getNumParticles());
            for (int i = 0; i < force.getNumInteractionGroups(); i++) {
                set<int> set1, set2;
                force.getInteractionGroupParameters(i, set1, set2);
                for (int p : set1)
                    groupsForParticle[p].insert(2*i);
                for (int p : set2)
                    groupsForParticle[p].insert(2*i+1);
            }
        }
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        vector<double> params1;
        vector<double> params2;
        force.getParticleParameters(particle1, params1);
        force.getParticleParameters(particle2, params2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i])
                return false;
        if (groupsForParticle.size() > 0 && groupsForParticle[particle1] != groupsForParticle[particle2])
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
    const CustomNonbondedForce& force;
    vector<set<int> > groupsForParticle;
};

OpenCLCalcCustomNonbondedForceKernel::~OpenCLCalcCustomNonbondedForceKernel() {
    if (params != NULL)
        delete params;
    if (forceCopy != NULL)
        delete forceCopy;
}

void OpenCLCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = (force.getNumInteractionGroups() == 0 ? "custom"+cl.intToString(forceIndex)+"_" : "");

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), numParticles, "customNonbondedParameters");
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customNonbondedGlobals", CL_MEM_READ_ONLY);
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

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    vector<string> tableTypes;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cl.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cl.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cl.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cl, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(arrayName, "float", width, width*sizeof(float), tabulatedFunctions[i].getDeviceBuffer()));
        if (width == 1)
            tableTypes.push_back("float");
        else
            tableTypes.push_back("float"+cl.intToString(width));
    }

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized())
        globals.upload(globalParamValues);
    bool useCutoff = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && force.getNonbondedMethod() != CustomNonbondedForce::CutoffNonPeriodic);
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
    map<string, Lepton::ParsedExpression> forceExpressions;
    forceExpressions["real customEnergy = "] = energyExpression;
    forceExpressions["tempForce -= "] = forceExpression;

    // Create the kernels.

    vector<pair<ExpressionTreeNode, string> > variables;
    ExpressionTreeNode rnode(new Operation::Variable("r"));
    variables.push_back(make_pair(rnode, "r"));
    variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
    variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        variables.push_back(makeVariable(name+"1", prefix+"params"+params->getParameterSuffix(i, "1")));
        variables.push_back(makeVariable(name+"2", prefix+"params"+params->getParameterSuffix(i, "2")));
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        string value = "globals["+cl.intToString(i)+"]";
        variables.push_back(makeVariable(name, prefix+value));
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cl.getNonbondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        forceExpressions[derivVariable+" += interactionScale*switchValue*"] = derivExpression;
    }
    stringstream compute;
    compute << cl.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, prefix+"temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["USE_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    if (force.getUseSwitchingFunction()) {
        // Compute the switching coefficients.
        
        replacements["SWITCH_CUTOFF"] = cl.doubleToString(force.getSwitchingDistance());
        replacements["SWITCH_C3"] = cl.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
        replacements["SWITCH_C4"] = cl.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
        replacements["SWITCH_C5"] = cl.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
    }
    string source = cl.replaceStrings(OpenCLKernelSources::customNonbonded, replacements);
    if (force.getNumInteractionGroups() > 0)
        initInteractionGroups(force, source, tableTypes);
    else {
        cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            cl.getNonbondedUtilities().addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"params"+cl.intToString(i+1), buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        if (globals.isInitialized()) {
            globals.upload(globalParamValues);
            cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(cl_float), globals.getDeviceBuffer()));
        }
    }
    info = new ForceInfo(cl.getNonbondedUtilities().getNumForceBuffers(), force);
    cl.addForce(info);
    
    // Record information for the long range correction.
    
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic && force.getUseLongRangeCorrection() && cl.getContextIndex() == 0) {
        forceCopy = new CustomNonbondedForce(force);
        hasInitializedLongRangeCorrection = false;
    }
    else {
        longRangeCoefficient = 0.0;
        hasInitializedLongRangeCorrection = true;
    }
}

void OpenCLCalcCustomNonbondedForceKernel::initInteractionGroups(const CustomNonbondedForce& force, const string& interactionSource, const vector<string>& tableTypes) {
    // Process groups to form tiles.
    
    vector<vector<int> > atomLists;
    vector<pair<int, int> > tiles;
    vector<int> tileGroup;
    vector<vector<int> > duplicateAtomsForGroup;
    for (int group = 0; group < force.getNumInteractionGroups(); group++) {
        // Get the list of atoms in this group and sort them.
        
        set<int> set1, set2;
        force.getInteractionGroupParameters(group, set1, set2);
        vector<int> atoms1, atoms2;
        atoms1.insert(atoms1.begin(), set1.begin(), set1.end());
        atoms2.insert(atoms2.begin(), set2.begin(), set2.end());
        sort(atoms1.begin(), atoms1.end());
        sort(atoms2.begin(), atoms2.end());
        duplicateAtomsForGroup.push_back(vector<int>());
        set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
                inserter(duplicateAtomsForGroup[group], duplicateAtomsForGroup[group].begin()));
        sort(duplicateAtomsForGroup[group].begin(), duplicateAtomsForGroup[group].end());
        
        // Find how many tiles we will create for this group.
        
        int tileWidth = min(min(32, (int) atoms1.size()), (int) atoms2.size());
        if (tileWidth == 0)
            continue;
        int numBlocks1 = (atoms1.size()+tileWidth-1)/tileWidth;
        int numBlocks2 = (atoms2.size()+tileWidth-1)/tileWidth;
        
        // Add the tiles.
        
        int firstTile = tiles.size();
        for (int i = 0; i < numBlocks1; i++)
            for (int j = 0; j < numBlocks2; j++) {
                tiles.push_back(make_pair(atomLists.size()+i, atomLists.size()+numBlocks1+j));
                tileGroup.push_back(group);
            }
        
        // Add the atom lists.
        
        for (int i = 0; i < numBlocks1; i++) {
            vector<int> atoms;
            int first = i*tileWidth;
            int last = min((i+1)*tileWidth, (int) atoms1.size());
            for (int j = first; j < last; j++)
                atoms.push_back(atoms1[j]);
            atomLists.push_back(atoms);
        }
        for (int i = 0; i < numBlocks2; i++) {
            vector<int> atoms;
            int first = i*tileWidth;
            int last = min((i+1)*tileWidth, (int) atoms2.size());
            for (int j = first; j < last; j++)
                atoms.push_back(atoms2[j]);
            atomLists.push_back(atoms);
        }
    }
    
    // Build a lookup table for quickly identifying excluded interactions.
    
    vector<set<int> > exclusions(force.getNumParticles());
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int p1, p2;
        force.getExclusionParticles(i, p1, p2);
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
    }
    
    // Build the exclusion flags for each tile.  While we're at it, filter out tiles
    // where all interactions are excluded, and sort the tiles by size.

    vector<vector<int> > exclusionFlags(tiles.size());
    vector<pair<int, int> > tileOrder;
    for (int tile = 0; tile < tiles.size(); tile++) {
        bool swapped = false;
        if (atomLists[tiles[tile].first].size() < atomLists[tiles[tile].second].size()) {
            // For efficiency, we want the first axis to be the larger one.
            
            int swap = tiles[tile].first;
            tiles[tile].first = tiles[tile].second;
            tiles[tile].second = swap;
            swapped = true;
        }
        vector<int>& atoms1 = atomLists[tiles[tile].first];
        vector<int>& atoms2 = atomLists[tiles[tile].second];
        vector<int>& duplicateAtoms = duplicateAtomsForGroup[tileGroup[tile]];
        vector<int>& flags = exclusionFlags[tile];
        flags.resize(atoms1.size(), (int) (1LL<<atoms2.size())-1);
        int numExcluded = 0;
        for (int i = 0; i < (int) atoms1.size(); i++) {
            int a1 = atoms1[i];
            bool a1IsDuplicate = binary_search(duplicateAtoms.begin(), duplicateAtoms.end(), a1);
            for (int j = 0; j < (int) atoms2.size(); j++) {
                int a2 = atoms2[j];
                bool isExcluded = false;
                if (a1 == a2 || exclusions[a1].find(a2) != exclusions[a1].end())
                    isExcluded = true; // This is an excluded interaction.
                else if ((a1 > a2) == swapped && a1IsDuplicate && binary_search(duplicateAtoms.begin(), duplicateAtoms.end(), a2))
                    isExcluded = true; // Both atoms are in both sets, so skip duplicate interactions.
                if (isExcluded) {
                    flags[i] &= -1-(1<<j);
                    numExcluded++;
                }
            }
        }
        if (numExcluded == atoms1.size()*atoms2.size())
            continue; // All interactions are excluded.
        tileOrder.push_back(make_pair((int) -atoms2.size(), tile));
    }
    sort(tileOrder.begin(), tileOrder.end());
    
    // Merge tiles to get as close as possible to 32 along the first axis of each one.
    
    vector<int> tileSetStart;
    tileSetStart.push_back(0);
    int tileSetSize = 0;
    for (int i = 0; i < tileOrder.size(); i++) {
        int tile = tileOrder[i].second;
        int size = atomLists[tiles[tile].first].size();
        if (tileSetSize+size > 32) {
            tileSetStart.push_back(i);
            tileSetSize = 0;
        }
        tileSetSize += size;
    }
    tileSetStart.push_back(tileOrder.size());
    
    // Build the data structures.
    
    int numTileSets = tileSetStart.size()-1;
    vector<mm_int4> groupData;
    for (int tileSet = 0; tileSet < numTileSets; tileSet++) {
        int indexInTileSet = 0;
        int minSize = 0;
        if (cl.getSIMDWidth() < 32) {
            // We need to include a barrier inside the inner loop, so ensure that all
            // threads will loop the same number of times.
            
            for (int i = tileSetStart[tileSet]; i < tileSetStart[tileSet+1]; i++)
                minSize = max(minSize, (int) atomLists[tiles[tileOrder[i].second].first].size());
        }
        for (int i = tileSetStart[tileSet]; i < tileSetStart[tileSet+1]; i++) {
            int tile = tileOrder[i].second;
            vector<int>& atoms1 = atomLists[tiles[tile].first];
            vector<int>& atoms2 = atomLists[tiles[tile].second];
            int range = indexInTileSet + ((indexInTileSet+max(minSize, (int) atoms1.size()))<<16);
            int allFlags = (1<<atoms2.size())-1;
            for (int j = 0; j < (int) atoms1.size(); j++) {
                int a1 = atoms1[j];
                int a2 = (j < atoms2.size() ? atoms2[j] : 0);
                int flags = (exclusionFlags[tile].size() > 0 ? exclusionFlags[tile][j] : allFlags);
                groupData.push_back(mm_int4(a1, a2, range, flags<<indexInTileSet));
            }
            indexInTileSet += atoms1.size();
        }
        for (; indexInTileSet < 32; indexInTileSet++)
            groupData.push_back(mm_int4(0, 0, minSize<<16, 0));
    }
    interactionGroupData.initialize<mm_int4>(cl, groupData.size(), "interactionGroupData");
    interactionGroupData.upload(groupData);
    numGroupTiles.initialize<cl_int>(cl, 1, "numGroupTiles");

    // Allocate space for a neighbor list, if necessary.

    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && groupData.size() > cl.getNumThreadBlocks()) {
        filteredGroupData.initialize<mm_int4>(cl, groupData.size(), "filteredGroupData");
        interactionGroupData.copyTo(filteredGroupData);
        int numTiles = groupData.size()/32;
        numGroupTiles.upload(&numTiles);
    }
    
    // Create the kernel.
    
    hasParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = interactionSource;
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream localData;
    int localDataSize = 0;
    vector<OpenCLNonbondedUtilities::ParameterInfo>& buffers = params->getBuffers(); 
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getNumComponents() == 1)
            localData<<buffers[i].getComponentType()<<" params"<<(i+1)<<";\n";
        else {
            for (int j = 0; j < buffers[i].getNumComponents(); ++j)
                localData<<buffers[i].getComponentType()<<" params"<<(i+1)<<"_"<<suffixes[j]<<";\n";
        }
        localDataSize += buffers[i].getSize();
    }
    replacements["ATOM_PARAMETER_DATA"] = localData.str();
    stringstream args;
    for (int i = 0; i < (int) buffers.size(); i++)
        args<<", __global const "<<buffers[i].getType()<<"* restrict global_params"<<(i+1);
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        args << ", __global const " << tableTypes[i]<< "* restrict table" << i;
    if (globals.isInitialized())
        args<<", __global const float* restrict globals";
    if (hasParamDerivs)
        args << ", __global mixed* restrict energyParamDerivs";
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    stringstream load1;
    for (int i = 0; i < (int) buffers.size(); i++)
        load1<<buffers[i].getType()<<" params"<<(i+1)<<"1 = global_params"<<(i+1)<<"[atom1];\n";
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    stringstream loadLocal2;
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getNumComponents() == 1)
            loadLocal2<<"localData[get_local_id(0)].params"<<(i+1)<<" = global_params"<<(i+1)<<"[atom2];\n";
        else {
            loadLocal2<<buffers[i].getType()<<" temp_params"<<(i+1)<<" = global_params"<<(i+1)<<"[atom2];\n";
            for (int j = 0; j < buffers[i].getNumComponents(); ++j)
                loadLocal2<<"localData[get_local_id(0)].params"<<(i+1)<<"_"<<suffixes[j]<<" = temp_params"<<(i+1)<<"."<<suffixes[j]<<";\n";
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS"] = loadLocal2.str();
    stringstream load2;
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getNumComponents() == 1)
            load2<<buffers[i].getType()<<" params"<<(i+1)<<"2 = localData[localIndex].params"<<(i+1)<<";\n";
        else {
            load2<<buffers[i].getType()<<" params"<<(i+1)<<"2 = ("<<buffers[i].getType()<<") (";
            for (int j = 0; j < buffers[i].getNumComponents(); ++j) {
                if (j > 0)
                    load2<<", ";
                load2<<"localData[localIndex].params"<<(i+1)<<"_"<<suffixes[j];
            }
            load2<<");\n";
        }
    }
    replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
    stringstream initDerivs, saveDerivs;
    const vector<string>& allParamDerivNames = cl.getEnergyParamDerivNames();
    int numDerivs = allParamDerivNames.size();
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cl.getNonbondedUtilities().addEnergyParameterDerivative(paramName);
        initDerivs<<"mixed "<<derivVariable<<" = 0;\n";
        for (int index = 0; index < numDerivs; index++)
            if (allParamDerivNames[index] == paramName)
                saveDerivs<<"energyParamDerivs[get_global_id(0)*"<<numDerivs<<"+"<<index<<"] += "<<derivVariable<<";\n";
    }
    replacements["INIT_DERIVATIVES"] = initDerivs.str();
    replacements["SAVE_DERIVATIVES"] = saveDerivs.str();
    map<string, string> defines;
    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
        defines["USE_CUTOFF"] = "1";
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic)
        defines["USE_PERIODIC"] = "1";
    int localMemorySize = max(32, cl.getNonbondedUtilities().getForceThreadBlockSize());
    defines["LOCAL_MEMORY_SIZE"] = cl.intToString(localMemorySize);
    defines["WARPS_IN_BLOCK"] = cl.intToString(localMemorySize/32);
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cl.doubleToString(cutoff*cutoff);
    double paddedCutoff = cl.getNonbondedUtilities().padCutoff(cutoff);
    defines["PADDED_CUTOFF_SQUARED"] = cl.doubleToString(paddedCutoff*paddedCutoff);
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    defines["TILE_SIZE"] = "32";
    defines["NUM_TILES"] = cl.intToString(numTileSets);
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*numTileSets/numContexts;
    int endIndex = (cl.getContextIndex()+1)*numTileSets/numContexts;
    defines["FIRST_TILE"] = cl.intToString(startIndex);
    defines["LAST_TILE"] = cl.intToString(endIndex);
    if ((localDataSize/4)%2 == 0 && !cl.getUseDoublePrecision())
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customNonbondedGroups, replacements), defines);
    interactionGroupKernel = cl::Kernel(program, "computeInteractionGroups");
    prepareNeighborListKernel = cl::Kernel(program, "prepareToBuildNeighborList");
    buildNeighborListKernel = cl::Kernel(program, "buildNeighborList");
    numGroupThreadBlocks = cl.getNonbondedUtilities().getNumForceThreadBlocks();
}

double OpenCLCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    useNeighborList = (filteredGroupData.isInitialized() && cl.getNonbondedUtilities().getUseCutoff());
    if (useNeighborList && cl.getContextIndex() > 0) {
        // When using a neighbor list, run the whole calculation on a single device.
        return 0.0;
    }
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed) {
            globals.upload(globalParamValues);
            if (forceCopy != NULL) {
                CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
                hasInitializedLongRangeCorrection = true;
            }
        }
    }
    if (!hasInitializedLongRangeCorrection) {
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
    }
    if (interactionGroupData.isInitialized()) {
        if (!hasInitializedKernel) {
            hasInitializedKernel = true;
            int index = 0;
            bool useLong = cl.getSupports64BitGlobalAtomics();
            interactionGroupKernel.setArg<cl::Buffer>(index++, (useLong ? cl.getLongForceBuffer() : cl.getForceBuffers()).getDeviceBuffer());
            interactionGroupKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
            interactionGroupKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
            interactionGroupKernel.setArg<cl::Buffer>(index++, (useNeighborList ? filteredGroupData : interactionGroupData).getDeviceBuffer());
            interactionGroupKernel.setArg<cl::Buffer>(index++, numGroupTiles.getDeviceBuffer());
            interactionGroupKernel.setArg<cl_int>(index++, useNeighborList);
            index += 5;
            for (auto& buffer : params->getBuffers())
                interactionGroupKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            for (auto& function : tabulatedFunctions)
                interactionGroupKernel.setArg<cl::Memory>(index++, function.getDeviceBuffer());
            if (globals.isInitialized())
                interactionGroupKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
            if (hasParamDerivs)
                interactionGroupKernel.setArg<cl::Memory>(index++, cl.getEnergyParamDerivBuffer().getDeviceBuffer());
            if (useNeighborList) {
                // Initialize kernels for building the interaction group neighbor list.
                
                prepareNeighborListKernel.setArg<cl::Buffer>(0, cl.getNonbondedUtilities().getRebuildNeighborList().getDeviceBuffer());
                prepareNeighborListKernel.setArg<cl::Buffer>(1, numGroupTiles.getDeviceBuffer());
                buildNeighborListKernel.setArg<cl::Buffer>(0, cl.getNonbondedUtilities().getRebuildNeighborList().getDeviceBuffer());
                buildNeighborListKernel.setArg<cl::Buffer>(1, numGroupTiles.getDeviceBuffer());
                buildNeighborListKernel.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
                buildNeighborListKernel.setArg<cl::Buffer>(3, interactionGroupData.getDeviceBuffer());
                buildNeighborListKernel.setArg<cl::Buffer>(4, filteredGroupData.getDeviceBuffer());
            }
        }
        int forceThreadBlockSize = max(32, cl.getNonbondedUtilities().getForceThreadBlockSize());
        if (useNeighborList) {
            // Rebuild the neighbor list, if necessary.

            setPeriodicBoxArgs(cl, buildNeighborListKernel, 5);
            cl.executeKernel(prepareNeighborListKernel, 1, 1);
            cl.executeKernel(buildNeighborListKernel, numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
        }
        setPeriodicBoxArgs(cl, interactionGroupKernel, 6);
        cl.executeKernel(interactionGroupKernel, numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    mm_double4 boxSize = cl.getPeriodicBoxSizeDouble();
    double volume = boxSize.x*boxSize.y*boxSize.z;
    map<string, double>& derivs = cl.getEnergyParamDerivWorkspace();
    for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
        derivs[forceCopy->getEnergyParameterDerivativeName(i)] += longRangeCoefficientDerivs[i]/volume;
    return longRangeCoefficient/volume;
}

void OpenCLCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
    int numParticles = force.getNumParticles();
    if (numParticles != cl.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<vector<cl_float> > paramVector(numParticles);
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // If necessary, recompute the long range correction.
    
    if (forceCopy != NULL) {
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcGBSAOBCForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const GBSAOBCForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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

void OpenCLCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    if (cl.getPlatformData().contexts.size() > 1)
        throw OpenMMException("GBSAOBCForce does not support using multiple OpenCL devices");
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "obc"+cl.intToString(forceIndex)+"_";
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    params.initialize<mm_float2>(cl, cl.getPaddedNumAtoms(), "gbsaObcParams");
    int elementSize = (cl.getUseDoublePrecision() ? sizeof(cl_double) : sizeof(cl_float));
    charges.initialize(cl, cl.getPaddedNumAtoms(), elementSize, "gbsaObcCharges");
    bornRadii.initialize(cl, cl.getPaddedNumAtoms(), elementSize, "bornRadii");
    obcChain.initialize(cl, cl.getPaddedNumAtoms(), elementSize, "obcChain");
    if (cl.getSupports64BitGlobalAtomics()) {
        longBornSum.initialize<cl_long>(cl, cl.getPaddedNumAtoms(), "longBornSum");
        longBornForce.initialize<cl_long>(cl, cl.getPaddedNumAtoms(), "longBornForce");
        bornForce.initialize(cl, cl.getPaddedNumAtoms(), elementSize, "bornForce");
        cl.addAutoclearBuffer(longBornSum);
        cl.addAutoclearBuffer(longBornForce);
    }
    else {
        bornSum.initialize(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), elementSize, "bornSum");
        bornForce.initialize(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), elementSize, "bornForce");
        cl.addAutoclearBuffer(bornSum);
        cl.addAutoclearBuffer(bornForce);
    }
    vector<double> chargeVec(cl.getPaddedNumAtoms());
    vector<mm_float2> paramsVector(cl.getPaddedNumAtoms(), mm_float2(1,1));
    const double dielectricOffset = 0.009;
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        chargeVec[i] = charge;
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
    }
    charges.upload(chargeVec, true, true);
    params.upload(paramsVector);
    prefactor = -ONE_4PI_EPS0*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
    surfaceAreaFactor = -6.0*4*M_PI*force.getSurfaceAreaEnergy();
    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
    cutoff = force.getCutoffDistance();
    string source = OpenCLKernelSources::gbsaObc2;
    map<string, string> replacements;
    replacements["CHARGE1"] = prefix+"charge1";
    replacements["CHARGE2"] = prefix+"charge2";
    replacements["OBC_PARAMS1"] = prefix+"obcParams1";
    replacements["OBC_PARAMS2"] = prefix+"obcParams2";
    replacements["BORN_FORCE1"] = prefix+"bornForce1";
    replacements["BORN_FORCE2"] = prefix+"bornForce2";
    source = cl.replaceStrings(source, replacements);
    nb.addInteraction(useCutoff, usePeriodic, false, cutoff, vector<vector<int> >(), source, force.getForceGroup());
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"charge", "float", 1, sizeof(cl_float), charges.getDeviceBuffer()));;
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"obcParams", "float", 2, sizeof(cl_float2), params.getDeviceBuffer()));;
    nb.addParameter(OpenCLNonbondedUtilities::ParameterInfo(prefix+"bornForce", "real", 1, elementSize, bornForce.getDeviceBuffer()));;
    info = new ForceInfo(nb.getNumForceBuffers(), force);
    cl.addForce(info);
}

double OpenCLCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    if (!hasCreatedKernels) {
        // These Kernels cannot be created in initialize(), because the OpenCLNonbondedUtilities has not been initialized yet then.

        hasCreatedKernels = true;
        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : 0);
        map<string, string> defines;
        if (nb.getUseCutoff())
            defines["USE_CUTOFF"] = "1";
        if (nb.getUsePeriodic())
            defines["USE_PERIODIC"] = "1";
        defines["CUTOFF_SQUARED"] = cl.doubleToString(cutoff*cutoff);
        defines["CUTOFF"] = cl.doubleToString(cutoff);
        defines["PREFACTOR"] = cl.doubleToString(prefactor);
        defines["SURFACE_AREA_FACTOR"] = cl.doubleToString(surfaceAreaFactor);
        defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
        defines["NUM_BLOCKS"] = cl.intToString(cl.getNumAtomBlocks());
        defines["FORCE_WORK_GROUP_SIZE"] = cl.intToString(nb.getForceThreadBlockSize());
        defines["TILE_SIZE"] = cl.intToString(OpenCLContext::TileSize);
        int numExclusionTiles = nb.getExclusionTiles().getSize();
        defines["NUM_TILES_WITH_EXCLUSIONS"] = cl.intToString(numExclusionTiles);
        int numContexts = cl.getPlatformData().contexts.size();
        int startExclusionIndex = cl.getContextIndex()*numExclusionTiles/numContexts;
        int endExclusionIndex = (cl.getContextIndex()+1)*numExclusionTiles/numContexts;
        defines["FIRST_EXCLUSION_TILE"] = cl.intToString(startExclusionIndex);
        defines["LAST_EXCLUSION_TILE"] = cl.intToString(endExclusionIndex);
        string platformVendor = cl::Platform(cl.getDevice().getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_VENDOR>();
        if (platformVendor == "Apple")
            defines["USE_APPLE_WORKAROUND"] = "1";
        string file;
        if (deviceIsCpu)
            file = OpenCLKernelSources::gbsaObc_cpu;
        else
            file = OpenCLKernelSources::gbsaObc;
        cl::Program program = cl.createProgram(file, defines);
        bool useLong = cl.getSupports64BitGlobalAtomics();
        int index = 0;
        computeBornSumKernel = cl::Kernel(program, "computeBornSum");
        computeBornSumKernel.setArg<cl::Buffer>(index++, (useLong ? longBornSum.getDeviceBuffer() : bornSum.getDeviceBuffer()));
        computeBornSumKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        computeBornSumKernel.setArg<cl::Buffer>(index++, charges.getDeviceBuffer());
        computeBornSumKernel.setArg<cl::Buffer>(index++, params.getDeviceBuffer());
        if (nb.getUseCutoff()) {
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
            index += 5; // The periodic box size arguments are set when the kernel is executed.
            computeBornSumKernel.setArg<cl_uint>(index++, maxTiles);
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getBlockCenters().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getBlockBoundingBoxes().getDeviceBuffer());
            computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getInteractingAtoms().getDeviceBuffer());
        }
        else
            computeBornSumKernel.setArg<cl_uint>(index++, cl.getNumAtomBlocks()*(cl.getNumAtomBlocks()+1)/2);
        computeBornSumKernel.setArg<cl::Buffer>(index++, nb.getExclusionTiles().getDeviceBuffer());
        force1Kernel = cl::Kernel(program, "computeGBSAForce1");
        index = 0;
        force1Kernel.setArg<cl::Buffer>(index++, (useLong ? cl.getLongForceBuffer().getDeviceBuffer() : cl.getForceBuffers().getDeviceBuffer()));
        force1Kernel.setArg<cl::Buffer>(index++, (useLong ? longBornForce.getDeviceBuffer() : bornForce.getDeviceBuffer()));
        force1Kernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, charges.getDeviceBuffer());
        force1Kernel.setArg<cl::Buffer>(index++, bornRadii.getDeviceBuffer());
        index++; // Whether to include energy.
        if (nb.getUseCutoff()) {
            force1Kernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
            index += 5; // The periodic box size arguments are set when the kernel is executed.
            force1Kernel.setArg<cl_uint>(index++, maxTiles);
            force1Kernel.setArg<cl::Buffer>(index++, nb.getBlockCenters().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(index++, nb.getBlockBoundingBoxes().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(index++, nb.getInteractingAtoms().getDeviceBuffer());
        }
        else
            force1Kernel.setArg<cl_uint>(index++, cl.getNumAtomBlocks()*(cl.getNumAtomBlocks()+1)/2);
        force1Kernel.setArg<cl::Buffer>(index++, nb.getExclusionTiles().getDeviceBuffer());
        program = cl.createProgram(OpenCLKernelSources::gbsaObcReductions, defines);
        reduceBornSumKernel = cl::Kernel(program, "reduceBornSum");
        reduceBornSumKernel.setArg<cl_int>(0, cl.getPaddedNumAtoms());
        reduceBornSumKernel.setArg<cl_int>(1, nb.getNumForceBuffers());
        reduceBornSumKernel.setArg<cl_float>(2, 1.0f);
        reduceBornSumKernel.setArg<cl_float>(3, 0.8f);
        reduceBornSumKernel.setArg<cl_float>(4, 4.85f);
        reduceBornSumKernel.setArg<cl::Buffer>(5, (useLong ? longBornSum.getDeviceBuffer() : bornSum.getDeviceBuffer()));
        reduceBornSumKernel.setArg<cl::Buffer>(6, params.getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(7, bornRadii.getDeviceBuffer());
        reduceBornSumKernel.setArg<cl::Buffer>(8, obcChain.getDeviceBuffer());
        reduceBornForceKernel = cl::Kernel(program, "reduceBornForce");
        index = 0;
        reduceBornForceKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        reduceBornForceKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
        reduceBornForceKernel.setArg<cl::Buffer>(index++, bornForce.getDeviceBuffer());
        if (useLong)
            reduceBornForceKernel.setArg<cl::Buffer>(index++, longBornForce.getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(index++, params.getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(index++, bornRadii.getDeviceBuffer());
        reduceBornForceKernel.setArg<cl::Buffer>(index++, obcChain.getDeviceBuffer());
    }
    force1Kernel.setArg<cl_int>(6, includeEnergy);
    if (nb.getUseCutoff()) {
        setPeriodicBoxArgs(cl, computeBornSumKernel, 6);
        setPeriodicBoxArgs(cl, force1Kernel, 9);
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            computeBornSumKernel.setArg<cl::Buffer>(4, nb.getInteractingTiles().getDeviceBuffer());
            computeBornSumKernel.setArg<cl_uint>(11, maxTiles);
            computeBornSumKernel.setArg<cl::Buffer>(14, nb.getInteractingAtoms().getDeviceBuffer());
            force1Kernel.setArg<cl::Buffer>(7, nb.getInteractingTiles().getDeviceBuffer());
            force1Kernel.setArg<cl_uint>(14, maxTiles);
            force1Kernel.setArg<cl::Buffer>(17, nb.getInteractingAtoms().getDeviceBuffer());
        }
    }
    cl.executeKernel(computeBornSumKernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cl.executeKernel(reduceBornSumKernel, cl.getPaddedNumAtoms());
    cl.executeKernel(force1Kernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cl.executeKernel(reduceBornForceKernel, cl.getPaddedNumAtoms());
    return 0.0;
}

void OpenCLCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
    // Make sure the new parameters are acceptable.
    
    int numParticles = force.getNumParticles();
    if (numParticles != cl.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<double> chargeVector(cl.getPaddedNumAtoms(), 0.0);
    vector<mm_float2> paramsVector(cl.getPaddedNumAtoms());
    const double dielectricOffset = 0.009;
    for (int i = 0; i < numParticles; i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        chargeVector[i] = charge;
        radius -= dielectricOffset;
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
    }
    for (int i = numParticles; i < cl.getPaddedNumAtoms(); i++)
        paramsVector[i] = mm_float2(1,1);
    charges.upload(chargeVector, true, true);
    params.upload(paramsVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCustomGBForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const CustomGBForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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
    const CustomGBForce& force;
};

OpenCLCalcCustomGBForceKernel::~OpenCLCalcCustomGBForceKernel() {
    if (params != NULL)
        delete params;
    if (computedValues != NULL)
        delete computedValues;
    if (energyDerivs != NULL)
        delete energyDerivs;
    if (energyDerivChain != NULL)
        delete energyDerivChain;
    for (auto d : dValuedParam)
        delete d;
}

void OpenCLCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    if (cl.getPlatformData().contexts.size() > 1)
        throw OpenMMException("CustomGBForce does not support using multiple OpenCL devices");
    cutoff = force.getCutoffDistance();
    bool useExclusionsForValue = false;
    numComputedValues = force.getNumComputedValues();
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
    string prefix = "custom"+cl.intToString(forceIndex)+"_";

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    int paddedNumParticles = cl.getPaddedNumAtoms();
    int numParams = force.getNumPerParticleParameters();
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), paddedNumParticles, "customGBParameters", true);
    computedValues = new OpenCLParameterSet(cl, force.getNumComputedValues(), paddedNumParticles, "customGBComputedValues", true, cl.getUseDoublePrecision());
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customGBGlobals", CL_MEM_READ_ONLY);
    vector<vector<cl_float> > paramVector(paddedNumParticles, vector<cl_float>(numParams, 0));
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
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

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cl.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cl.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cl.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cl, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        cl.getNonbondedUtilities().addArgument(OpenCLNonbondedUtilities::ParameterInfo(arrayName, "float", width, width*sizeof(float), tabulatedFunctions[i].getDeviceBuffer()));
        tableArgs << ", __global const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* restrict " << arrayName;
    }

    // Record the global parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized())
        globals.upload(globalParamValues);

    // Record derivatives of expressions needed for the chain rule terms.

    vector<vector<Lepton::ParsedExpression> > valueGradientExpressions(force.getNumComputedValues());
    vector<vector<Lepton::ParsedExpression> > valueDerivExpressions(force.getNumComputedValues());
    vector<vector<Lepton::ParsedExpression> > valueParamDerivExpressions(force.getNumComputedValues());
    needParameterGradient = false;
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
        if (i > 0) {
            valueGradientExpressions[i].push_back(ex.differentiate("x").optimize());
            valueGradientExpressions[i].push_back(ex.differentiate("y").optimize());
            valueGradientExpressions[i].push_back(ex.differentiate("z").optimize());
            if (!isZeroExpression(valueGradientExpressions[i][0]) || !isZeroExpression(valueGradientExpressions[i][1]) || !isZeroExpression(valueGradientExpressions[i][2]))
                needParameterGradient = true;
             for (int j = 0; j < i; j++)
                valueDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
            valueParamDerivExpressions[i].push_back(ex.differentiate(force.getEnergyParameterDerivativeName(j)).optimize());
    }
    vector<vector<Lepton::ParsedExpression> > energyDerivExpressions(force.getNumEnergyTerms());
    vector<vector<Lepton::ParsedExpression> > energyParamDerivExpressions(force.getNumEnergyTerms());
    vector<bool> needChainForValue(force.getNumComputedValues(), false);
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        for (int j = 0; j < force.getNumComputedValues(); j++) {
            if (type == CustomGBForce::SingleParticle) {
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
                if (!isZeroExpression(energyDerivExpressions[i].back()))
                    needChainForValue[j] = true;
            }
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"1").optimize());
                if (!isZeroExpression(energyDerivExpressions[i].back()))
                    needChainForValue[j] = true;
                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"2").optimize());
                if (!isZeroExpression(energyDerivExpressions[i].back()))
                    needChainForValue[j] = true;
            }
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
            energyParamDerivExpressions[i].push_back(ex.differentiate(force.getEnergyParameterDerivativeName(j)).optimize());
    }
    bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    bool useLong = cl.getSupports64BitGlobalAtomics();
    if (useLong) {
        longEnergyDerivs.initialize<cl_long>(cl, force.getNumComputedValues()*cl.getPaddedNumAtoms(), "customGBLongEnergyDerivatives");
        energyDerivs = new OpenCLParameterSet(cl, force.getNumComputedValues(), cl.getPaddedNumAtoms(), "customGBEnergyDerivatives", true);
    }
    else
        energyDerivs = new OpenCLParameterSet(cl, force.getNumComputedValues(), cl.getPaddedNumAtoms()*cl.getNonbondedUtilities().getNumForceBuffers(), "customGBEnergyDerivatives", true);
    energyDerivChain = new OpenCLParameterSet(cl, force.getNumComputedValues(), cl.getPaddedNumAtoms(), "customGBEnergyDerivativeChain", true);
    int elementSize = (cl.getUseDoublePrecision() ? sizeof(cl_double) : sizeof(cl_float));
    needEnergyParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    dValue0dParam.resize(force.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        dValuedParam.push_back(new OpenCLParameterSet(cl, force.getNumComputedValues(), cl.getPaddedNumAtoms(), "dValuedParam", true, cl.getUseDoublePrecision()));
        if (useLong)
            dValue0dParam[i].initialize<cl_long>(cl, cl.getPaddedNumAtoms(), "dValue0dParam");
        else
            dValue0dParam[i].initialize(cl, cl.getPaddedNumAtoms()*cl.getNonbondedUtilities().getNumForceBuffers(), elementSize, "dValue0dParam");
        cl.addAutoclearBuffer(dValue0dParam[i]);
        string name = force.getEnergyParameterDerivativeName(i);
        cl.addEnergyParameterDerivative(name);
    }

    // Create the kernels.

    bool useCutoff = (force.getNonbondedMethod() != CustomGBForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomGBForce::NoCutoff && force.getNonbondedMethod() != CustomGBForce::CutoffNonPeriodic);
    {
        // Create the N2 value kernel.

        vector<pair<ExpressionTreeNode, string> > variables;
        map<string, string> rename;
        ExpressionTreeNode rnode(new Operation::Variable("r"));
        variables.push_back(make_pair(rnode, "r"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables.push_back(makeVariable(name+"1", "params"+params->getParameterSuffix(i, "1")));
            variables.push_back(makeVariable(name+"2", "params"+params->getParameterSuffix(i, "2")));
            rename[name+"1"] = name+"2";
            rename[name+"2"] = name+"1";
        }
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cl.intToString(i)+"]";
            variables.push_back(makeVariable(name, value));
        }
        map<string, Lepton::ParsedExpression> n2ValueExpressions;
        stringstream n2ValueSource;
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[0], functions).optimize();
        n2ValueExpressions["tempValue1 = "] = ex;
        n2ValueExpressions["tempValue2 = "] = ex.renameVariables(rename);
        for (int i = 0; i < valueParamDerivExpressions[0].size(); i++) {
            string variableBase = "temp_dValue0dParam"+cl.intToString(i+1);
            if (!isZeroExpression(valueParamDerivExpressions[0][i])) {
                n2ValueExpressions[variableBase+"_1 = "] = valueParamDerivExpressions[0][i];
                n2ValueExpressions[variableBase+"_2 = "] = valueParamDerivExpressions[0][i].renameVariables(rename);
            }
        }
        n2ValueSource << cl.getExpressionUtilities().createExpressions(n2ValueExpressions, variables, functionList, functionDefinitions, "temp");
        map<string, string> replacements;
        string n2ValueStr = n2ValueSource.str();
        replacements["COMPUTE_VALUE"] = n2ValueStr;
        stringstream extraArgs, loadLocal1, loadLocal2, load1, load2, tempDerivs1, tempDerivs2, storeDeriv1, storeDeriv2;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __global const float* globals";
        pairValueUsesParam.resize(params->getBuffers().size(), false);
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cl.intToString(i+1);
            if (n2ValueStr.find(paramName+"1") != n2ValueStr.npos || n2ValueStr.find(paramName+"2") != n2ValueStr.npos) {
                extraArgs << ", __global const " << buffer.getType() << "* restrict global_" << paramName << ", __local " << buffer.getType() << "* restrict local_" << paramName;
                loadLocal1 << "local_" << paramName << "[localAtomIndex] = " << paramName << "1;\n";
                loadLocal2 << "local_" << paramName << "[localAtomIndex] = global_" << paramName << "[j];\n";
                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
                load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
                pairValueUsesParam[i] = true;
            }
        }
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string derivName = "dValue0dParam"+cl.intToString(i+1);
            if (useLong)
                extraArgs << ", __global long* restrict global_" << derivName;
            else
                extraArgs << ", __global real* restrict global_" << derivName;
            extraArgs << ", __local real* restrict local_" << derivName;
            loadLocal2 << "local_" << derivName << "[localAtomIndex] = 0;\n";
            load1 << "real " << derivName << " = 0;\n";
            if (!isZeroExpression(valueParamDerivExpressions[0][i])) {
                load2 << "real temp_" << derivName << "_1 = 0;\n";
                load2 << "real temp_" << derivName << "_2 = 0;\n";
                tempDerivs1 << derivName << " += temp_" << derivName << "_1;\n";
                if (deviceIsCpu)
                    tempDerivs2 << "local_" << derivName << "[j] += temp_" << derivName << "_2;\n";
                else
                    tempDerivs2 << "local_" << derivName << "[tbx+tj] += temp_" << derivName << "_2;\n";
                if (useLong) {
                    storeDeriv1 << "atom_add(&global_" << derivName << "[offset1], (long) (" << derivName << "*0x100000000));\n";
                    if (deviceIsCpu)
                        storeDeriv2 << "atom_add(&global_" << derivName << "[offset2], (long) (local_" << derivName << "[tgx]*0x100000000));\n";
                    else
                        storeDeriv2 << "atom_add(&global_" << derivName << "[offset2], (long) (local_" << derivName << "[get_local_id(0)]*0x100000000));\n";
                }
                else {
                    storeDeriv1 << "global_" << derivName << "[offset1] += " << derivName << ";\n";
                    if (deviceIsCpu)
                        storeDeriv2 << "global_" << derivName << "[offset2] += local_" << derivName << "[tgx];\n";
                    else
                        storeDeriv2 << "global_" << derivName << "[offset2] += local_" << derivName << "[get_local_id(0)];\n";
                }
            }
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        replacements["ADD_TEMP_DERIVS1"] = tempDerivs1.str();
        replacements["ADD_TEMP_DERIVS2"] = tempDerivs2.str();
        replacements["STORE_PARAM_DERIVS1"] = storeDeriv1.str();
        replacements["STORE_PARAM_DERIVS2"] = storeDeriv2.str();
        if (useCutoff)
            pairValueDefines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            pairValueDefines["USE_PERIODIC"] = "1";
        if (useExclusionsForValue)
            pairValueDefines["USE_EXCLUSIONS"] = "1";
        pairValueDefines["FORCE_WORK_GROUP_SIZE"] = cl.intToString(cl.getNonbondedUtilities().getForceThreadBlockSize());
        pairValueDefines["CUTOFF_SQUARED"] = cl.doubleToString(cutoff*cutoff);
        pairValueDefines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        pairValueDefines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
        pairValueDefines["NUM_BLOCKS"] = cl.intToString(cl.getNumAtomBlocks());
        pairValueDefines["TILE_SIZE"] = cl.intToString(OpenCLContext::TileSize);
        string file;
        if (deviceIsCpu)
            file = OpenCLKernelSources::customGBValueN2_cpu;
        else
            file = OpenCLKernelSources::customGBValueN2;
        pairValueSrc = cl.replaceStrings(file, replacements);
        if (useExclusionsForValue)
            cl.getNonbondedUtilities().requestExclusions(exclusionList);
    }
    {
        // Create the kernel to reduce the N2 value and calculate other values.

        stringstream reductionSource, extraArgs, deriv0;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __global const float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cl.intToString(i+1);
            extraArgs << ", __global const " << buffer.getType() << "* restrict " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cl.intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* restrict global_" << valueName;
            reductionSource << buffer.getType() << " local_" << valueName << ";\n";
        }
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string variableName = "dValuedParam_0_"+cl.intToString(i);
            if (useLong) {
                extraArgs << ", __global const long* restrict dValue0dParam" << i;
                deriv0 << "real " << variableName << " = (1.0f/0x100000000)*dValue0dParam" << i << "[index];\n";
            }
            else {
                extraArgs << ", __global const real* restrict dValue0dParam" << i;
                deriv0 << "real " << variableName << " = dValue0dParam" << i << "[index];\n";
                deriv0 << "for (int i = index+bufferSize; i < totalSize; i += bufferSize)\n";
                deriv0 << "    " << variableName << " += dValue0dParam" << i << "[i];\n";
            }
            for (int j = 0; j < dValuedParam[i]->getBuffers().size(); j++)
                extraArgs << ", __global real* restrict global_dValuedParam_" << j << "_" << i;
            deriv0 << "global_dValuedParam_0_" << i << "[index] = dValuedParam_0_" << i << ";\n";
        }
        reductionSource << "local_values" << computedValues->getParameterSuffix(0) << " = sum;\n";
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cl.intToString(i)+"]";
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            variables[computedValueNames[i-1]] = "local_values"+computedValues->getParameterSuffix(i-1);
            map<string, Lepton::ParsedExpression> valueExpressions;
            valueExpressions["local_values"+computedValues->getParameterSuffix(i)+" = "] = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
            reductionSource << cl.getExpressionUtilities().createExpressions(valueExpressions, variables, functionList, functionDefinitions, "value"+cl.intToString(i)+"_temp");
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            string valueName = "values"+cl.intToString(i+1);
            reductionSource << "global_" << valueName << "[index] = local_" << valueName << ";\n";
        }
        if (needEnergyParamDerivs) {
            map<string, Lepton::ParsedExpression> derivExpressions;
            for (int i = 1; i < force.getNumComputedValues(); i++) {
                for (int j = 0; j < valueParamDerivExpressions[i].size(); j++)
                    derivExpressions["real dValuedParam_"+cl.intToString(i)+"_"+cl.intToString(j)+" = "] = valueParamDerivExpressions[i][j];
                for (int j = 0; j < i; j++)
                    derivExpressions["real dVdV_"+cl.intToString(i)+"_"+cl.intToString(j)+" = "] = valueDerivExpressions[i][j];
            }
            reductionSource << cl.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, "derivChain_temp");
            for (int i = 1; i < force.getNumComputedValues(); i++) {
                for (int j = 0; j < i; j++)
                    for (int k = 0; k < valueParamDerivExpressions[i].size(); k++)
                        reductionSource << "dValuedParam_" << i << "_" << k << " += dVdV_" << i << "_" << j << "*dValuedParam_" << j <<"_" << k << ";\n";
                for (int j = 0; j < valueParamDerivExpressions[i].size(); j++)
                    reductionSource << "global_dValuedParam_" << i << "_" << j << "[index] = dValuedParam_" << i << "_" << j << ";\n";
            }
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["REDUCE_PARAM0_DERIV"] = deriv0.str();
        replacements["COMPUTE_VALUES"] = reductionSource.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customGBValuePerParticle, replacements), defines);
        perParticleValueKernel = cl::Kernel(program, "computePerParticleValues");
    }
    {
        // Create the N2 energy kernel.

        vector<pair<ExpressionTreeNode, string> > variables;
        ExpressionTreeNode rnode(new Operation::Variable("r"));
        variables.push_back(make_pair(rnode, "r"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables.push_back(makeVariable(name+"1", "params"+params->getParameterSuffix(i, "1")));
            variables.push_back(makeVariable(name+"2", "params"+params->getParameterSuffix(i, "2")));
        }
        for (int i = 0; i < force.getNumComputedValues(); i++) {
            variables.push_back(makeVariable(computedValueNames[i]+"1", "values"+computedValues->getParameterSuffix(i, "1")));
            variables.push_back(makeVariable(computedValueNames[i]+"2", "values"+computedValues->getParameterSuffix(i, "2")));
        }
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables.push_back(makeVariable(force.getGlobalParameterName(i), "globals["+cl.intToString(i)+"]"));
        stringstream n2EnergySource;
        bool anyExclusions = (force.getNumExclusions() > 0);
        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
            string expression;
            CustomGBForce::ComputationType type;
            force.getEnergyTermParameters(i, expression, type);
            if (type == CustomGBForce::SingleParticle)
                continue;
            bool exclude = (anyExclusions && type == CustomGBForce::ParticlePair);
            map<string, Lepton::ParsedExpression> n2EnergyExpressions;
            n2EnergyExpressions["tempEnergy += "] = Lepton::Parser::parse(expression, functions).optimize();
            n2EnergyExpressions["dEdR += "] = Lepton::Parser::parse(expression, functions).differentiate("r").optimize();
            if (useLong) {
                for (int j = 0; j < force.getNumComputedValues(); j++) {
                    if (needChainForValue[j]) {
                        string index = cl.intToString(j+1);
                        n2EnergyExpressions["/*"+cl.intToString(i+1)+"*/ deriv"+index+"_1 += "] = energyDerivExpressions[i][2*j];
                        n2EnergyExpressions["/*"+cl.intToString(i+1)+"*/ deriv"+index+"_2 += "] = energyDerivExpressions[i][2*j+1];
                    }
                }
            }
            else {
                for (int j = 0; j < force.getNumComputedValues(); j++) {
                    if (needChainForValue[j]) {
                        n2EnergyExpressions["/*"+cl.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j, "_1")+" += "] = energyDerivExpressions[i][2*j];
                        n2EnergyExpressions["/*"+cl.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j, "_2")+" += "] = energyDerivExpressions[i][2*j+1];
                    }
                }
            }
            for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
                n2EnergyExpressions["energyParamDeriv"+cl.intToString(j)+" += interactionScale*"] = energyParamDerivExpressions[i][j];
            if (exclude)
                n2EnergySource << "if (!isExcluded) {\n";
            n2EnergySource << cl.getExpressionUtilities().createExpressions(n2EnergyExpressions, variables, functionList, functionDefinitions, "temp");
            if (exclude)
                n2EnergySource << "}\n";
        }
        map<string, string> replacements;
        string n2EnergyStr = n2EnergySource.str();
        replacements["COMPUTE_INTERACTION"] = n2EnergyStr;
        stringstream extraArgs, loadLocal1, loadLocal2, clearLocal, load1, load2, declare1, recordDeriv, storeDerivs1, storeDerivs2, declareTemps, setTemps, initParamDerivs, saveParamDerivs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __global const float* globals";
        pairEnergyUsesParam.resize(params->getBuffers().size(), false);
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cl.intToString(i+1);
            if (n2EnergyStr.find(paramName+"1") != n2EnergyStr.npos || n2EnergyStr.find(paramName+"2") != n2EnergyStr.npos) {
                extraArgs << ", __global const " << buffer.getType() << "* restrict global_" << paramName << ", __local " << buffer.getType() << "* restrict local_" << paramName;
                loadLocal1 << "local_" << paramName << "[localAtomIndex] = " << paramName << "1;\n";
                loadLocal2 << "local_" << paramName << "[localAtomIndex] = global_" << paramName << "[j];\n";
                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
                load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
                pairEnergyUsesParam[i] = true;
            }
        }
        pairEnergyUsesValue.resize(computedValues->getBuffers().size(), false);
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cl.intToString(i+1);
            if (n2EnergyStr.find(valueName+"1") != n2EnergyStr.npos || n2EnergyStr.find(valueName+"2") != n2EnergyStr.npos) {
                extraArgs << ", __global const " << buffer.getType() << "* restrict global_" << valueName << ", __local " << buffer.getType() << "* restrict local_" << valueName;
                loadLocal1 << "local_" << valueName << "[localAtomIndex] = " << valueName << "1;\n";
                loadLocal2 << "local_" << valueName << "[localAtomIndex] = global_" << valueName << "[j];\n";
                load1 << buffer.getType() << " " << valueName << "1 = global_" << valueName << "[atom1];\n";
                load2 << buffer.getType() << " " << valueName << "2 = local_" << valueName << "[atom2];\n";
                pairEnergyUsesValue[i] = true;
            }
        }
        if (useLong) {
            extraArgs << ", __global long* restrict derivBuffers";
            for (int i = 0; i < force.getNumComputedValues(); i++) {
                string index = cl.intToString(i+1);
                extraArgs << ", __local real* restrict local_deriv" << index;
                clearLocal << "local_deriv" << index << "[localAtomIndex] = 0.0f;\n";
                declare1 << "real deriv" << index << "_1 = 0;\n";
                load2 << "real deriv" << index << "_2 = 0;\n";
                recordDeriv << "local_deriv" << index << "[atom2] += deriv" << index << "_2;\n";
                storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")\n";
                storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")\n";
                declareTemps << "__local real tempDerivBuffer" << index << "[64];\n";
                setTemps << "tempDerivBuffer" << index << "[get_local_id(0)] = deriv" << index << "_1;\n";
            }
        }
        else {
            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
                const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
                string index = cl.intToString(i+1);
                extraArgs << ", __global " << buffer.getType() << "* restrict derivBuffers" << index << ", __local " << buffer.getType() << "* restrict local_deriv" << index;
                clearLocal << "local_deriv" << index << "[localAtomIndex] = 0.0f;\n";
                declare1 << buffer.getType() << " deriv" << index << "_1 = 0.0f;\n";
                load2 << buffer.getType() << " deriv" << index << "_2 = 0.0f;\n";
                recordDeriv << "local_deriv" << index << "[atom2] += deriv" << index << "_2;\n";
                storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")\n";
                storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")\n";
                declareTemps << "__local " << buffer.getType() << " tempDerivBuffer" << index << "[64];\n";
                setTemps << "tempDerivBuffer" << index << "[get_local_id(0)] = deriv" << index << "_1;\n";
            }
        }
        if (needEnergyParamDerivs) {
            extraArgs << ", __global mixed* restrict energyParamDerivs";
            const vector<string>& allParamDerivNames = cl.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[get_global_id(0)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["CLEAR_LOCAL_DERIVATIVES"] = clearLocal.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        replacements["DECLARE_ATOM1_DERIVATIVES"] = declare1.str();
        replacements["RECORD_DERIVATIVE_2"] = recordDeriv.str();
        replacements["STORE_DERIVATIVES_1"] = storeDerivs1.str();
        replacements["STORE_DERIVATIVES_2"] = storeDerivs2.str();
        replacements["DECLARE_TEMP_BUFFERS"] = declareTemps.str();
        replacements["SET_TEMP_BUFFERS"] = setTemps.str();
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        if (useCutoff)
            pairEnergyDefines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            pairEnergyDefines["USE_PERIODIC"] = "1";
        if (anyExclusions)
            pairEnergyDefines["USE_EXCLUSIONS"] = "1";
        pairEnergyDefines["FORCE_WORK_GROUP_SIZE"] = cl.intToString(cl.getNonbondedUtilities().getForceThreadBlockSize());
        pairEnergyDefines["CUTOFF_SQUARED"] = cl.doubleToString(cutoff*cutoff);
        pairEnergyDefines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        pairEnergyDefines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
        pairEnergyDefines["NUM_BLOCKS"] = cl.intToString(cl.getNumAtomBlocks());
        pairEnergyDefines["TILE_SIZE"] = cl.intToString(OpenCLContext::TileSize);
        string file;
        if (deviceIsCpu)
            file = OpenCLKernelSources::customGBEnergyN2_cpu;
        else
            file = OpenCLKernelSources::customGBEnergyN2;
        pairEnergySrc = cl.replaceStrings(file, replacements);
    }
    {
        // Create the kernel to reduce the derivatives and calculate per-particle energy terms.

        stringstream compute, extraArgs, reduce, initParamDerivs, saveParamDerivs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __global const float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cl.intToString(i+1);
            extraArgs << ", __global const " << buffer.getType() << "* restrict " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cl.intToString(i+1);
            extraArgs << ", __global const " << buffer.getType() << "* restrict " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string index = cl.intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* restrict derivBuffers" << index;
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        for (int i = 0; i < (int) energyDerivChain->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivChain->getBuffers()[i];
            string index = cl.intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* restrict derivChain" << index;
        }
        if (useLong) {
            extraArgs << ", __global const long* restrict derivBuffersIn";
            for (int i = 0; i < energyDerivs->getNumParameters(); ++i)
                reduce << "derivBuffers" << energyDerivs->getParameterSuffix(i, "[index]") <<
                        " = (1.0f/0x100000000)*derivBuffersIn[index+PADDED_NUM_ATOMS*" << cl.intToString(i) << "];\n";
        }
        else {
            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
                reduce << "REDUCE_VALUE(derivBuffers" << cl.intToString(i+1) << ", " << energyDerivs->getBuffers()[i].getType() << ")\n";
        }
        if (needEnergyParamDerivs) {
            extraArgs << ", __global mixed* restrict energyParamDerivs";
            const vector<string>& allParamDerivNames = cl.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[get_global_id(0)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
        }
        
        // Compute the various expressions.
        
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cl.intToString(i)+"]";
        for (int i = 0; i < force.getNumComputedValues(); i++)
            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
        map<string, Lepton::ParsedExpression> expressions;
        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
            string expression;
            CustomGBForce::ComputationType type;
            force.getEnergyTermParameters(i, expression, type);
            if (type != CustomGBForce::SingleParticle)
                continue;
            Lepton::ParsedExpression parsed = Lepton::Parser::parse(expression, functions).optimize();
            expressions["/*"+cl.intToString(i+1)+"*/ energy += "] = parsed;
            for (int j = 0; j < force.getNumComputedValues(); j++)
                expressions["/*"+cl.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j)+" += "] = energyDerivExpressions[i][j];
            Lepton::ParsedExpression gradx = parsed.differentiate("x").optimize();
            Lepton::ParsedExpression grady = parsed.differentiate("y").optimize();
            Lepton::ParsedExpression gradz = parsed.differentiate("z").optimize();
            if (!isZeroExpression(gradx))
                expressions["/*"+cl.intToString(i+1)+"*/ force.x -= "] = gradx;
            if (!isZeroExpression(grady))
                expressions["/*"+cl.intToString(i+1)+"*/ force.y -= "] = grady;
            if (!isZeroExpression(gradz))
                expressions["/*"+cl.intToString(i+1)+"*/ force.z -= "] = gradz;
            for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
                expressions["/*"+cl.intToString(i+1)+"*/ energyParamDeriv"+cl.intToString(j)+" += "] = energyParamDerivExpressions[i][j];
        }
        for (int i = 1; i < force.getNumComputedValues(); i++)
            for (int j = 0; j < i; j++)
                expressions["real dV"+cl.intToString(i)+"dV"+cl.intToString(j)+" = "] = valueDerivExpressions[i][j];
        compute << cl.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "temp");
        
        // Record values.
        
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            string index = cl.intToString(i+1);
            compute << "derivBuffers" << index << "[index] = deriv" << index << ";\n";
        }
        compute << "forceBuffers[index] = forceBuffers[index]+force;\n";
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            compute << "real totalDeriv"<<i<<" = dV"<<i<<"dV0";
            for (int j = 1; j < i; j++)
                compute << " + totalDeriv"<<j<<"*dV"<<i<<"dV"<<j;
            compute << ";\n";
            compute << "deriv"<<(i+1)<<" *= totalDeriv"<<i<<";\n";
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            string index = cl.intToString(i+1);
            compute << "derivChain" << index << "[index] = deriv" << index << ";\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["REDUCE_DERIVATIVES"] = reduce.str();
        replacements["COMPUTE_ENERGY"] = compute.str();
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
        cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customGBEnergyPerParticle, replacements), defines);
        perParticleEnergyKernel = cl::Kernel(program, "computePerParticleEnergy");
    }
    if (needParameterGradient || needEnergyParamDerivs) {
        // Create the kernel to compute chain rule terms for computed values that depend explicitly on particle coordinates, and for
        // derivatives with respect to global parameters.

        stringstream compute, extraArgs, initParamDerivs, saveParamDerivs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", __global const float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cl.intToString(i+1);
            extraArgs << ", __global const " << buffer.getType() << "* restrict " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cl.intToString(i+1);
            extraArgs << ", __global const " << buffer.getType() << "* restrict " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string index = cl.intToString(i+1);
            extraArgs << ", __global " << buffer.getType() << "* restrict derivBuffers" << index;
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        if (needEnergyParamDerivs) {
            extraArgs << ", __global mixed* restrict energyParamDerivs";
            const vector<string>& allParamDerivNames = cl.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                for (int j = 0; j < dValuedParam[i]->getBuffers().size(); j++)
                    extraArgs << ", __global real* restrict dValuedParam_" << j << "_" << i;
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[get_global_id(0)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
        }
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cl.intToString(i)+"]";
        for (int i = 0; i < force.getNumComputedValues(); i++)
            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
        if (needParameterGradient) {
            for (int i = 1; i < force.getNumComputedValues(); i++) {
                string is = cl.intToString(i);
                compute << "real4 dV"<<is<<"dR = (real4) 0;\n";
                for (int j = 1; j < i; j++) {
                    if (!isZeroExpression(valueDerivExpressions[i][j])) {
                        map<string, Lepton::ParsedExpression> derivExpressions;
                        string js = cl.intToString(j);
                        derivExpressions["real dV"+is+"dV"+js+" = "] = valueDerivExpressions[i][j];
                        compute << cl.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, "temp_"+is+"_"+js);
                        compute << "dV"<<is<<"dR += dV"<<is<<"dV"<<js<<"*dV"<<js<<"dR;\n";
                    }
                }
                map<string, Lepton::ParsedExpression> gradientExpressions;
                if (!isZeroExpression(valueGradientExpressions[i][0]))
                    gradientExpressions["dV"+is+"dR.x += "] = valueGradientExpressions[i][0];
                if (!isZeroExpression(valueGradientExpressions[i][1]))
                    gradientExpressions["dV"+is+"dR.y += "] = valueGradientExpressions[i][1];
                if (!isZeroExpression(valueGradientExpressions[i][2]))
                    gradientExpressions["dV"+is+"dR.z += "] = valueGradientExpressions[i][2];
                compute << cl.getExpressionUtilities().createExpressions(gradientExpressions, variables, functionList, functionDefinitions, "temp");
            }
            for (int i = 1; i < force.getNumComputedValues(); i++)
                compute << "force -= deriv"<<energyDerivs->getParameterSuffix(i)<<"*dV"<<i<<"dR;\n";
        }
        if (needEnergyParamDerivs)
            for (int i = 0; i < force.getNumComputedValues(); i++)
                for (int j = 0; j < dValuedParam.size(); j++)
                    compute << "energyParamDeriv"<<j<<" += deriv"<<energyDerivs->getParameterSuffix(i)<<"*dValuedParam_"<<i<<"_"<<j<<"[index];\n";
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["COMPUTE_FORCES"] = compute.str();
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customGBGradientChainRule, replacements), defines);
        gradientChainRuleKernel = cl::Kernel(program, "computeGradientChainRuleTerms");
    }
    {
        // Create the code to calculate chain rule terms as part of the default nonbonded kernel.

        vector<pair<ExpressionTreeNode, string> > globalVariables;
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cl.intToString(i)+"]";
            globalVariables.push_back(makeVariable(name, prefix+value));
        }
        vector<pair<ExpressionTreeNode, string> > variables = globalVariables;
        map<string, string> rename;
        ExpressionTreeNode rnode(new Operation::Variable("r"));
        variables.push_back(make_pair(rnode, "r"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
            const string& name = force.getPerParticleParameterName(i);
            variables.push_back(makeVariable(name+"1", prefix+"params"+params->getParameterSuffix(i, "1")));
            variables.push_back(makeVariable(name+"2", prefix+"params"+params->getParameterSuffix(i, "2")));
            rename[name+"1"] = name+"2";
            rename[name+"2"] = name+"1";
        }
        map<string, Lepton::ParsedExpression> derivExpressions;
        stringstream chainSource;
        Lepton::ParsedExpression dVdR = Lepton::Parser::parse(computedValueExpressions[0], functions).differentiate("r").optimize();
        derivExpressions["real dV0dR1 = "] = dVdR;
        derivExpressions["real dV0dR2 = "] = dVdR.renameVariables(rename);
        chainSource << cl.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, prefix+"temp0_");
        if (needChainForValue[0]) {
            if (useExclusionsForValue)
                chainSource << "if (!isExcluded) {\n";
            chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "1") << ";\n";
            chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "2") << ";\n";
            if (useExclusionsForValue)
                chainSource << "}\n";
        }
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            if (needChainForValue[i]) {
                chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "1") << ";\n";
                chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "2") << ";\n";
            }
        }
        map<string, string> replacements;
        string chainStr = chainSource.str();
        replacements["COMPUTE_FORCE"] = chainStr;
        string source = cl.replaceStrings(OpenCLKernelSources::customGBChainRule, replacements);
        vector<OpenCLNonbondedUtilities::ParameterInfo> parameters;
        vector<OpenCLNonbondedUtilities::ParameterInfo> arguments;
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = prefix+"params"+cl.intToString(i+1);
            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
                parameters.push_back(OpenCLNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string paramName = prefix+"values"+cl.intToString(i+1);
            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
                parameters.push_back(OpenCLNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        for (int i = 0; i < (int) energyDerivChain->getBuffers().size(); i++) {
            if (needChainForValue[i]) { 
                const OpenCLNonbondedUtilities::ParameterInfo& buffer = energyDerivChain->getBuffers()[i];
                string paramName = prefix+"dEdV"+cl.intToString(i+1);
                parameters.push_back(OpenCLNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
            }
        }
        if (globals.isInitialized()) {
            globals.upload(globalParamValues);
            arguments.push_back(OpenCLNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(cl_float), globals.getDeviceBuffer()));
        }
        cl.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, force.getNumExclusions() > 0, cutoff, exclusionList, source, force.getForceGroup());
        for (auto param : parameters)
            cl.getNonbondedUtilities().addParameter(param);
        for (auto arg : arguments)
            cl.getNonbondedUtilities().addArgument(arg);
    }
    info = new ForceInfo(cl.getNonbondedUtilities().getNumForceBuffers(), force);
    cl.addForce(info);
    if (useLong)
        cl.addAutoclearBuffer(longEnergyDerivs);
    else {
        for (auto& buffer : energyDerivs->getBuffers())
            cl.addAutoclearBuffer(buffer.getMemory(), buffer.getSize()*energyDerivs->getNumObjects());
    }
}

double OpenCLCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    bool deviceIsCpu = (cl.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    OpenCLNonbondedUtilities& nb = cl.getNonbondedUtilities();
    int elementSize = (cl.getUseDoublePrecision() ? sizeof(cl_double) : sizeof(cl_float));
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // These two kernels can't be compiled in initialize(), because the nonbonded utilities object
        // has not yet been initialized then.

        {
            int numExclusionTiles = nb.getExclusionTiles().getSize();
            pairValueDefines["NUM_TILES_WITH_EXCLUSIONS"] = cl.intToString(numExclusionTiles);
            int numContexts = cl.getPlatformData().contexts.size();
            int startExclusionIndex = cl.getContextIndex()*numExclusionTiles/numContexts;
            int endExclusionIndex = (cl.getContextIndex()+1)*numExclusionTiles/numContexts;
            pairValueDefines["FIRST_EXCLUSION_TILE"] = cl.intToString(startExclusionIndex);
            pairValueDefines["LAST_EXCLUSION_TILE"] = cl.intToString(endExclusionIndex);
            pairValueDefines["CUTOFF"] = cl.doubleToString(cutoff);
            cl::Program program = cl.createProgram(pairValueSrc, pairValueDefines);
            pairValueKernel = cl::Kernel(program, "computeN2Value");
            pairValueSrc = "";
            pairValueDefines.clear();
        }
        {
            int numExclusionTiles = nb.getExclusionTiles().getSize();
            pairEnergyDefines["NUM_TILES_WITH_EXCLUSIONS"] = cl.intToString(numExclusionTiles);
            int numContexts = cl.getPlatformData().contexts.size();
            int startExclusionIndex = cl.getContextIndex()*numExclusionTiles/numContexts;
            int endExclusionIndex = (cl.getContextIndex()+1)*numExclusionTiles/numContexts;
            pairEnergyDefines["FIRST_EXCLUSION_TILE"] = cl.intToString(startExclusionIndex);
            pairEnergyDefines["LAST_EXCLUSION_TILE"] = cl.intToString(endExclusionIndex);
            pairEnergyDefines["CUTOFF"] = cl.doubleToString(cutoff);
            cl::Program program = cl.createProgram(pairEnergySrc, pairEnergyDefines);
            pairEnergyKernel = cl::Kernel(program, "computeN2Energy");
            pairEnergySrc = "";
            pairEnergyDefines.clear();
        }

        // Set arguments for kernels.
        
        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : 0);
        bool useLong = cl.getSupports64BitGlobalAtomics();
        if (useLong) {
            longValueBuffers.initialize<cl_long>(cl, cl.getPaddedNumAtoms(), "customGBLongValueBuffers");
            cl.addAutoclearBuffer(longValueBuffers);
            cl.clearBuffer(longValueBuffers);
        }
        else {
            valueBuffers.initialize(cl, cl.getPaddedNumAtoms()*nb.getNumForceBuffers(), elementSize, "customGBValueBuffers");
            cl.addAutoclearBuffer(valueBuffers);
            cl.clearBuffer(valueBuffers);
        }
        int index = 0;
        pairValueKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        pairValueKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*4*elementSize, NULL);
        pairValueKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusions().getDeviceBuffer());
        pairValueKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusionTiles().getDeviceBuffer());
        pairValueKernel.setArg<cl::Buffer>(index++, useLong ? longValueBuffers.getDeviceBuffer() : valueBuffers.getDeviceBuffer());
        pairValueKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*elementSize, NULL);
        if (nb.getUseCutoff()) {
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
            index += 5; // Periodic box size arguments are set when the kernel is executed.
            pairValueKernel.setArg<cl_uint>(index++, maxTiles);
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getBlockCenters().getDeviceBuffer());
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getBlockBoundingBoxes().getDeviceBuffer());
            pairValueKernel.setArg<cl::Buffer>(index++, nb.getInteractingAtoms().getDeviceBuffer());
        }
        else
            pairValueKernel.setArg<cl_uint>(index++, cl.getNumAtomBlocks()*(cl.getNumAtomBlocks()+1)/2);
        if (globals.isInitialized())
            pairValueKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            if (pairValueUsesParam[i]) {
                const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
                pairValueKernel.setArg<cl::Memory>(index++, buffer.getMemory());
                pairValueKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
            }
        }
        for (auto& d : dValue0dParam) {
            pairValueKernel.setArg<cl::Buffer>(index++, d.getDeviceBuffer());
            pairValueKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*d.getElementSize(), NULL);
        }
        for (auto& function : tabulatedFunctions)
            pairValueKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        index = 0;
        perParticleValueKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        perParticleValueKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
        perParticleValueKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        perParticleValueKernel.setArg<cl::Buffer>(index++, useLong ? longValueBuffers.getDeviceBuffer() : valueBuffers.getDeviceBuffer());
        if (globals.isInitialized())
            perParticleValueKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (auto& buffer : params->getBuffers())
            perParticleValueKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& buffer : computedValues->getBuffers())
            perParticleValueKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (int i = 0; i < dValuedParam.size(); i++) {
            perParticleValueKernel.setArg<cl::Memory>(index++, dValue0dParam[i].getDeviceBuffer());
            for (int j = 0; j < dValuedParam[i]->getBuffers().size(); j++)
                perParticleValueKernel.setArg<cl::Memory>(index++, dValuedParam[i]->getBuffers()[j].getMemory());
        }
        for (auto& function : tabulatedFunctions)
            perParticleValueKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        index = 0;
        pairEnergyKernel.setArg<cl::Buffer>(index++, useLong ? cl.getLongForceBuffer().getDeviceBuffer() : cl.getForceBuffers().getDeviceBuffer());
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        pairEnergyKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*4*elementSize, NULL);
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        pairEnergyKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*4*elementSize, NULL);
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusions().getDeviceBuffer());
        pairEnergyKernel.setArg<cl::Buffer>(index++, cl.getNonbondedUtilities().getExclusionTiles().getDeviceBuffer());
        index++; // Whether to include energy.
        if (nb.getUseCutoff()) {
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
            index += 5; // Periodic box size arguments are set when the kernel is executed.
            pairEnergyKernel.setArg<cl_uint>(index++, maxTiles);
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getBlockCenters().getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getBlockBoundingBoxes().getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(index++, nb.getInteractingAtoms().getDeviceBuffer());
        }
        else
            pairEnergyKernel.setArg<cl_uint>(index++, cl.getNumAtomBlocks()*(cl.getNumAtomBlocks()+1)/2);
        if (globals.isInitialized())
            pairEnergyKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            if (pairEnergyUsesParam[i]) {
                const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
                pairEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
                pairEnergyKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
            }
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            if (pairEnergyUsesValue[i]) {
                const OpenCLNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
                pairEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
                pairEnergyKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
            }
        }
        if (useLong) {
            pairEnergyKernel.setArg<cl::Memory>(index++, longEnergyDerivs.getDeviceBuffer());
            for (int i = 0; i < numComputedValues; ++i)
                pairEnergyKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*elementSize, NULL);
        }
        else {
            for (auto& buffer : energyDerivs->getBuffers()) {
                pairEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
                pairEnergyKernel.setArg(index++, (deviceIsCpu ? OpenCLContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
            }
        }
        if (needEnergyParamDerivs)
            pairEnergyKernel.setArg<cl::Memory>(index++, cl.getEnergyParamDerivBuffer().getDeviceBuffer());
        for (auto& function : tabulatedFunctions)
            pairEnergyKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        index = 0;
        perParticleEnergyKernel.setArg<cl_int>(index++, cl.getPaddedNumAtoms());
        perParticleEnergyKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
        perParticleEnergyKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
        perParticleEnergyKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        perParticleEnergyKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        if (globals.isInitialized())
            perParticleEnergyKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (auto& buffer : params->getBuffers())
            perParticleEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& buffer : computedValues->getBuffers())
            perParticleEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& buffer : energyDerivs->getBuffers())
            perParticleEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& buffer : energyDerivChain->getBuffers())
            perParticleEnergyKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        if (useLong)
            perParticleEnergyKernel.setArg<cl::Memory>(index++, longEnergyDerivs.getDeviceBuffer());
        if (needEnergyParamDerivs)
            perParticleEnergyKernel.setArg<cl::Memory>(index++, cl.getEnergyParamDerivBuffer().getDeviceBuffer());
        for (auto& function : tabulatedFunctions)
            perParticleEnergyKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        if (needParameterGradient || needEnergyParamDerivs) {
            index = 0;
            gradientChainRuleKernel.setArg<cl::Buffer>(index++, cl.getForceBuffers().getDeviceBuffer());
            gradientChainRuleKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
            if (globals.isInitialized())
                gradientChainRuleKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
            for (auto& buffer : params->getBuffers())
                gradientChainRuleKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            for (auto& buffer : computedValues->getBuffers())
                gradientChainRuleKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            for (auto& buffer : energyDerivs->getBuffers())
                gradientChainRuleKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            if (needEnergyParamDerivs) {
                gradientChainRuleKernel.setArg<cl::Buffer>(index++, cl.getEnergyParamDerivBuffer().getDeviceBuffer());
                for (auto d : dValuedParam)
                    for (auto& buffer : d->getBuffers())
                        gradientChainRuleKernel.setArg<cl::Memory>(index++, buffer.getMemory());
            }
            for (auto& function : tabulatedFunctions)
                gradientChainRuleKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        }
    }
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
    pairEnergyKernel.setArg<cl_int>(7, includeEnergy);
    if (nb.getUseCutoff()) {
        setPeriodicBoxArgs(cl, pairValueKernel, 8);
        setPeriodicBoxArgs(cl, pairEnergyKernel, 10);
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            pairValueKernel.setArg<cl::Buffer>(6, nb.getInteractingTiles().getDeviceBuffer());
            pairValueKernel.setArg<cl_uint>(13, maxTiles);
            pairValueKernel.setArg<cl::Buffer>(16, nb.getInteractingAtoms().getDeviceBuffer());
            pairEnergyKernel.setArg<cl::Buffer>(8, nb.getInteractingTiles().getDeviceBuffer());
            pairEnergyKernel.setArg<cl_uint>(15, maxTiles);
            pairEnergyKernel.setArg<cl::Buffer>(18, nb.getInteractingAtoms().getDeviceBuffer());
        }
    }
    cl.executeKernel(pairValueKernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cl.executeKernel(perParticleValueKernel, cl.getPaddedNumAtoms());
    cl.executeKernel(pairEnergyKernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cl.executeKernel(perParticleEnergyKernel, cl.getPaddedNumAtoms());
    if (needParameterGradient || needEnergyParamDerivs)
        cl.executeKernel(gradientChainRuleKernel, cl.getPaddedNumAtoms());
    return 0.0;
}

void OpenCLCalcCustomGBForceKernel::copyParametersToContext(ContextImpl& context, const CustomGBForce& force) {
    int numParticles = force.getNumParticles();
    if (numParticles != cl.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<vector<cl_float> > paramVector(cl.getPaddedNumAtoms(), vector<cl_float>(force.getNumPerParticleParameters(), 0));
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCustomExternalForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomExternalForce& force, int numParticles) : OpenCLForceInfo(0), force(force), indices(numParticles, -1) {
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
}

void OpenCLCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumParticles()/numContexts;
    numParticles = endIndex-startIndex;
    if (numParticles == 0)
        return;
    vector<vector<int> > atoms(numParticles, vector<int>(1));
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), numParticles, "customExternalParams");
    vector<vector<cl_float> > paramVector(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(startIndex+i, atoms[i][0], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force, system.getNumParticles());
    cl.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    map<string, Lepton::CustomFunction*> customFunctions;
    customFunctions["periodicdistance"] = cl.getExpressionUtilities().getPeriodicDistancePlaceholder();
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
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customExternalGlobals", CL_MEM_READ_ONLY);
        globals.upload(globalParamValues);
        string argName = cl.getBondedUtilities().addArgument(globals.getDeviceBuffer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cl.intToString(i)+"]";
            variables[name] = value;
        }
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cl.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" particleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cl.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLKernelSources::customExternalForce, replacements), force.getForceGroup());
}

double OpenCLCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    return 0.0;
}

void OpenCLCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumParticles()/numContexts;
    if (numParticles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    if (numParticles == 0)
        return;
    
    // Record the per-particle parameters.
    
    vector<vector<cl_float> > paramVector(numParticles);
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        int particle;
        force.getParticleParameters(startIndex+i, particle, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
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

class OpenCLCalcCustomCentroidBondForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomCentroidBondForce& force) : OpenCLForceInfo(0), force(force) {
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

OpenCLCalcCustomCentroidBondForceKernel::~OpenCLCalcCustomCentroidBondForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcCustomCentroidBondForceKernel::initialize(const System& system, const CustomCentroidBondForce& force) {
    numBonds = force.getNumBonds();
    if (numBonds == 0)
        return;
    if (!cl.getSupports64BitGlobalAtomics())
        throw OpenMMException("CustomCentroidBondForce requires a device that supports 64 bit atomic operations");
    info = new ForceInfo(force);
    cl.addForce(info);
    
    // Record the groups.
    
    numGroups = force.getNumGroups();
    vector<cl_int> groupParticleVec;
    vector<cl_double> groupWeightVec;
    vector<cl_int> groupOffsetVec;
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
    groupParticles.initialize<int>(cl, groupParticleVec.size(), "groupParticles");
    groupParticles.upload(groupParticleVec);
    if (cl.getUseDoublePrecision()) {
        groupWeights.initialize<double>(cl, groupParticleVec.size(), "groupWeights");
        centerPositions.initialize<mm_double4>(cl, numGroups, "centerPositions");
    }
    else {
        groupWeights.initialize<float>(cl, groupParticleVec.size(), "groupWeights");
        centerPositions.initialize<mm_float4>(cl, numGroups, "centerPositions");
    }
    groupWeights.upload(groupWeightVec, true, true);
    groupOffsets.initialize<int>(cl, groupOffsetVec.size(), "groupOffsets");
    groupOffsets.upload(groupOffsetVec);
    groupForces.initialize<long long>(cl, numGroups*3, "groupForces");
    cl.addAutoclearBuffer(groupForces);
    
    // Record the bonds.
    
    int groupsPerBond = force.getNumGroupsPerBond();
    vector<cl_int> bondGroupVec(numBonds*groupsPerBond);
    params = new OpenCLParameterSet(cl, force.getNumPerBondParameters(), numBonds, "customCentroidBondParams");
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
    bondGroups.initialize<int>(cl, bondGroupVec.size(), "bondGroups");
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
        string arrayName = "table"+cl.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cl.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cl.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cl, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        extraArgs << ", __global const float";
        if (width > 1)
            extraArgs << width;
        extraArgs << "* restrict " << arrayName;
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
        string index = cl.intToString(i+1);
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
        extraArgs << ", __global mixed* restrict energyParamDerivs";
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cl, force.getNumGlobalParameters(), "customCentroidBondGlobals");
        globals.upload(globalParamValues);
        extraArgs << ", __global const float* restrict globals";
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cl.intToString(i)+"]";
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
        string index = cl.intToString(i+1);
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
        forceExpressions["real dEdDistance"+cl.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
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
        forceExpressions["real dEdAngle"+cl.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
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
        forceExpressions["real dEdDihedral"+cl.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArgs<<", __global const "<<buffer.getType()<<"* restrict globalParams"<<i;
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = globalParams"<<i<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    if (needEnergyParamDerivs) {
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string paramName = force.getEnergyParameterDerivativeName(i);
            cl.addEnergyParameterDerivative(paramName);
            Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
            forceExpressions[string("energyParamDeriv")+cl.intToString(i)+" += "] = derivExpression;
            initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
        }
        const vector<string>& allParamDerivNames = cl.getEnergyParamDerivNames();
        int numDerivs = allParamDerivNames.size();
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
            for (int index = 0; index < numDerivs; index++)
                if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                    saveParamDerivs << "energyParamDerivs[get_global_id(0)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
    }
    compute << cl.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to groups.

    vector<string> forceNames;
    for (int i = 0; i < groupsPerBond; i++) {
        string istr = cl.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = (real3) 0;\n";
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
            compute<<cl.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (auto& distance : distances) {
        const vector<int>& groups = distance.second;
        string deltaName = atomNames[groups[0]]+atomNames[groups[1]];
        string value = "(dEdDistance"+cl.intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
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
        compute<<"real4 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(length(crossProd), (real) 1e-6f);\n";
        compute<<"real4 deltaCross0 = -cross(delta"<<deltaName1<<", crossProd)*dEdAngle"<<cl.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real4 deltaCross2 = cross(delta"<<deltaName2<<", crossProd)*dEdAngle"<<cl.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real4 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[groups[0]]<<".xyz += deltaCross0.xyz;\n";
        compute<<forceNames[groups[1]]<<".xyz += deltaCross1.xyz;\n";
        compute<<forceNames[groups[2]]<<".xyz += deltaCross2.xyz;\n";
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
        compute<<"ff.x = (-dEdDihedral"<<cl.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cl.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real4 internalF0 = ff.x*"<<crossName1<<";\n";
        compute<<"real4 internalF3 = ff.w*"<<crossName2<<";\n";
        compute<<"real4 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[groups[0]]<<".xyz += internalF0.xyz;\n";
        compute<<forceNames[groups[1]]<<".xyz += s.xyz-internalF0.xyz;\n";
        compute<<forceNames[groups[2]]<<".xyz += -s.xyz-internalF3.xyz;\n";
        compute<<forceNames[groups[3]]<<".xyz += internalF3.xyz;\n";
        compute<<"}\n";
        index++;
    }
    
    // Save the forces to global memory.
    
    for (int i = 0; i < groupsPerBond; i++) {
        compute<<"atom_add(&groupForce[group"<<(i+1)<<"], (long) (force"<<(i+1)<<".x*0x100000000));\n";
        compute<<"atom_add(&groupForce[group"<<(i+1)<<"+NUM_GROUPS], (long) (force"<<(i+1)<<".y*0x100000000));\n";
        compute<<"atom_add(&groupForce[group"<<(i+1)<<"+NUM_GROUPS*2], (long) (force"<<(i+1)<<".z*0x100000000));\n";
    }
    map<string, string> replacements;
    replacements["M_PI"] = cl.doubleToString(M_PI);
    replacements["NUM_GROUPS"] = cl.intToString(numGroups);
    replacements["NUM_BONDS"] = cl.intToString(numBonds);
    replacements["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    replacements["EXTRA_ARGS"] = extraArgs.str();
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
    replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customCentroidBond, replacements));
    index = 0;
    computeCentersKernel = cl::Kernel(program, "computeGroupCenters");
    computeCentersKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
    computeCentersKernel.setArg<cl::Buffer>(index++, groupParticles.getDeviceBuffer());
    computeCentersKernel.setArg<cl::Buffer>(index++, groupWeights.getDeviceBuffer());
    computeCentersKernel.setArg<cl::Buffer>(index++, groupOffsets.getDeviceBuffer());
    computeCentersKernel.setArg<cl::Buffer>(index++, centerPositions.getDeviceBuffer());
    index = 0;
    groupForcesKernel = cl::Kernel(program, "computeGroupForces");
    groupForcesKernel.setArg<cl::Buffer>(index++, groupForces.getDeviceBuffer());
    index++; // Energy buffer hasn't been created yet
    groupForcesKernel.setArg<cl::Buffer>(index++, centerPositions.getDeviceBuffer());
    groupForcesKernel.setArg<cl::Buffer>(index++, bondGroups.getDeviceBuffer());
    index += 5; // Periodic box information
    if (needEnergyParamDerivs)
        index++; // Deriv buffer hasn't been created yet.
    for (auto& function : tabulatedFunctions)
        groupForcesKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
    if (globals.isInitialized())
        groupForcesKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
    for (auto& buffer : params->getBuffers())
        groupForcesKernel.setArg<cl::Memory>(index++, buffer.getMemory());
    index = 0;
    applyForcesKernel = cl::Kernel(program, "applyForcesToAtoms");
    applyForcesKernel.setArg<cl::Buffer>(index++, groupParticles.getDeviceBuffer());
    applyForcesKernel.setArg<cl::Buffer>(index++, groupWeights.getDeviceBuffer());
    applyForcesKernel.setArg<cl::Buffer>(index++, groupOffsets.getDeviceBuffer());
    applyForcesKernel.setArg<cl::Buffer>(index++, groupForces.getDeviceBuffer());
}

double OpenCLCalcCustomCentroidBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    cl.executeKernel(computeCentersKernel, OpenCLContext::TileSize*numGroups);
    groupForcesKernel.setArg<cl::Buffer>(1, cl.getEnergyBuffer().getDeviceBuffer());
    setPeriodicBoxArgs(cl, groupForcesKernel, 4);
    if (needEnergyParamDerivs)
        groupForcesKernel.setArg<cl::Memory>(9, cl.getEnergyParamDerivBuffer().getDeviceBuffer());
    cl.executeKernel(groupForcesKernel, numBonds);
    applyForcesKernel.setArg<cl::Buffer>(4, cl.getLongForceBuffer().getDeviceBuffer());
    cl.executeKernel(applyForcesKernel, OpenCLContext::TileSize*numGroups);
    return 0.0;
}

void OpenCLCalcCustomCentroidBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force) {
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
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCustomCompoundBondForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomCompoundBondForce& force) : OpenCLForceInfo(0), force(force) {
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

OpenCLCalcCustomCompoundBondForceKernel::~OpenCLCalcCustomCompoundBondForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    int particlesPerBond = force.getNumParticlesPerBond();
    vector<vector<int> > atoms(numBonds, vector<int>(particlesPerBond));
    params = new OpenCLParameterSet(cl, force.getNumPerBondParameters(), numBonds, "customCompoundBondParams");
    vector<vector<cl_float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        force.getBondParameters(startIndex+i, atoms[i], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    info = new ForceInfo(force);
    cl.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        functions[name] = cl.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cl.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cl, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        string arrayName = cl.getBondedUtilities().addArgument(tabulatedFunctions[i].getDeviceBuffer(), width == 1 ? "float" : "float"+cl.intToString(width));
        functionDefinitions.push_back(make_pair(name, arrayName));
    }
    
    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
    }
    map<string, string> variables;
    for (int i = 0; i < particlesPerBond; i++) {
        string index = cl.intToString(i+1);
        variables["x"+index] = "pos"+index+".x";
        variables["y"+index] = "pos"+index+".y";
        variables["z"+index] = "pos"+index+".z";
    }
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customCompoundBondGlobals", CL_MEM_READ_ONLY);
        globals.upload(globalParamValues);
        string argName = cl.getBondedUtilities().addArgument(globals.getDeviceBuffer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cl.intToString(i)+"]";
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
        string index = cl.intToString(i+1);
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
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
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
            compute<<"real4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = ccb_computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
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
        forceExpressions["real dEdDihedral"+cl.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        const OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cl.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cl.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        forceExpressions[derivVariable+" += "] = derivExpression;
    }
    compute << cl.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    vector<string> forceNames;
    for (int i = 0; i < particlesPerBond; i++) {
        string istr = cl.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real4 "<<forceName<<" = (real4) 0;\n";
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
            compute<<cl.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cl.intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
        compute<<forceNames[atoms[0]]<<".xyz += "<<"-"<<value<<";\n";
        compute<<forceNames[atoms[1]]<<".xyz += "<<value<<";\n";
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute<<"{\n";
        compute<<"real4 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(length(crossProd), (real) 1e-6f);\n";
        compute<<"real4 deltaCross0 = -cross(delta"<<deltaName1<<", crossProd)*dEdAngle"<<cl.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real4 deltaCross2 = cross(delta"<<deltaName2<<", crossProd)*dEdAngle"<<cl.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real4 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[atoms[0]]<<".xyz += deltaCross0.xyz;\n";
        compute<<forceNames[atoms[1]]<<".xyz += deltaCross1.xyz;\n";
        compute<<forceNames[atoms[2]]<<".xyz += deltaCross2.xyz;\n";
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
        compute<<"ff.x = (-dEdDihedral"<<cl.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cl.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real4 internalF0 = ff.x*"<<crossName1<<";\n";
        compute<<"real4 internalF3 = ff.w*"<<crossName2<<";\n";
        compute<<"real4 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[atoms[0]]<<".xyz += internalF0.xyz;\n";
        compute<<forceNames[atoms[1]]<<".xyz += s.xyz-internalF0.xyz;\n";
        compute<<forceNames[atoms[2]]<<".xyz += -s.xyz-internalF3.xyz;\n";
        compute<<forceNames[atoms[3]]<<".xyz += internalF3.xyz;\n";
        compute<<"}\n";
        index++;
    }
    cl.getBondedUtilities().addInteraction(atoms, compute.str(), force.getForceGroup());
    map<string, string> replacements;
    replacements["M_PI"] = cl.doubleToString(M_PI);
    cl.getBondedUtilities().addPrefixCode(cl.replaceStrings(OpenCLKernelSources::customCompoundBond, replacements));;
}

double OpenCLCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    return 0.0;
}

void OpenCLCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<vector<cl_float> > paramVector(numBonds);
    vector<int> particles;
    vector<double> parameters;
    for (int i = 0; i < numBonds; i++) {
        force.getBondParameters(startIndex+i, particles, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (cl_float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcCustomManyParticleForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(const CustomManyParticleForce& force) : OpenCLForceInfo(0), force(force) {
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

OpenCLCalcCustomManyParticleForceKernel::~OpenCLCalcCustomManyParticleForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {
    if (!cl.getSupports64BitGlobalAtomics())
        throw OpenMMException("CustomManyParticleForce requires a device that supports 64 bit atomic operations");
    int numParticles = force.getNumParticles();
    int particlesPerSet = force.getNumParticlesPerSet();
    bool centralParticleMode = (force.getPermutationMode() == CustomManyParticleForce::UniqueCentralParticle);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    forceWorkgroupSize = 128;
    findNeighborsWorkgroupSize = (cl.getSIMDWidth() >= 32 ? 128 : 32);
    
    // Record parameter values.
    
    params = new OpenCLParameterSet(cl, force.getNumPerParticleParameters(), numParticles, "customManyParticleParameters");
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
    cl.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
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
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    vector<pair<ExpressionTreeNode, string> > variables;
    for (int i = 0; i < particlesPerSet; i++) {
        string index = cl.intToString(i+1);
        variables.push_back(makeVariable("x"+index, "pos"+index+".x"));
        variables.push_back(makeVariable("y"+index, "pos"+index+".y"));
        variables.push_back(makeVariable("z"+index, "pos"+index+".z"));
    }
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        for (int j = 0; j < particlesPerSet; j++) {
            string index = cl.intToString(j+1);
            variables.push_back(makeVariable(name+index, "params"+params->getParameterSuffix(i, index)));
        }
    }
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<cl_float>(cl, force.getNumGlobalParameters(), "customManyParticleGlobals", CL_MEM_READ_ONLY);
        globals.upload(globalParamValues);
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cl.intToString(i)+"]";
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
        particleTypes.initialize<int>(cl, particleTypesVec.size(), "customManyParticleTypes");
        orderIndex.initialize<int>(cl, orderIndexVec.size(), "customManyParticleOrderIndex");
        particleOrder.initialize<int>(cl, particleOrderVec.size()*particlesPerSet, "customManyParticleOrder");
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
        exclusions.initialize<int>(cl, exclusionsVec.size(), "customManyParticleExclusions");
        exclusionStartIndex.initialize<int>(cl, exclusionStartIndexVec.size(), "customManyParticleExclusionStart");
        exclusions.upload(exclusionsVec);
        exclusionStartIndex.upload(exclusionStartIndexVec);
    }
    
    // Build data structures for the neighbor list.
    
    if (nonbondedMethod != NoCutoff) {
        int numAtomBlocks = cl.getNumAtomBlocks();
        int elementSize = (cl.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter.initialize(cl, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(cl, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        numNeighborPairs.initialize<int>(cl, 1, "customManyParticleNumNeighborPairs");
        neighborStartIndex.initialize<int>(cl, numParticles+1, "customManyParticleNeighborStartIndex");
        numNeighborsForAtom.initialize<int>(cl, numParticles, "customManyParticleNumNeighborsForAtom");

        // Select a size for the array that holds the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxNeighborPairs = 150*numParticles;
        neighborPairs.initialize<mm_int2>(cl, maxNeighborPairs, "customManyParticleNeighborPairs");
        neighbors.initialize<int>(cl, maxNeighborPairs, "customManyParticleNeighbors");
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
        string index = cl.intToString(i+1);
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
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        variables.push_back(makeVariable(angle.first, angleName));
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
        forceExpressions["real dEdDihedral"+cl.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        compute<<buffer.getType()<<" params"<<(i+1)<<" = global_params"<<(i+1)<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    compute << cl.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Apply forces to atoms.

    vector<string> forceNames;
    for (int i = 0; i < particlesPerSet; i++) {
        string istr = cl.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real4 "<<forceName<<" = (real4) 0;\n";
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
            compute<<cl.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cl.intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
        compute<<forceNames[atoms[0]]<<".xyz += "<<"-"<<value<<";\n";
        compute<<forceNames[atoms[1]]<<".xyz += "<<value<<";\n";
        index++;
    }
    index = 0;
    for (auto& angle : angles) {
        const vector<int>& atoms = angle.second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute<<"{\n";
        compute<<"real4 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), (real) 1e-6f);\n";
        compute<<"real4 deltaCross0 = -cross(delta"<<deltaName1<<", crossProd)*dEdAngle"<<cl.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real4 deltaCross2 = cross(delta"<<deltaName2<<", crossProd)*dEdAngle"<<cl.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real4 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[atoms[0]]<<".xyz += deltaCross0.xyz;\n";
        compute<<forceNames[atoms[1]]<<".xyz += deltaCross1.xyz;\n";
        compute<<forceNames[atoms[2]]<<".xyz += deltaCross2.xyz;\n";
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
        compute<<"ff.x = (-dEdDihedral"<<cl.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cl.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real4 internalF0 = ff.x*"<<crossName1<<";\n";
        compute<<"real4 internalF3 = ff.w*"<<crossName2<<";\n";
        compute<<"real4 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[atoms[0]]<<".xyz += internalF0.xyz;\n";
        compute<<forceNames[atoms[1]]<<".xyz += s.xyz-internalF0.xyz;\n";
        compute<<forceNames[atoms[2]]<<".xyz += -s.xyz-internalF3.xyz;\n";
        compute<<forceNames[atoms[3]]<<".xyz += internalF3.xyz;\n";
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
        loadData<<"real4 pos"<<(i+1)<<" = posq[atom"<<(i+1)<<"];\n";
        for (int j = 0; j < (int) params->getBuffers().size(); j++)
            loadData<<params->getBuffers()[j].getType()<<" params"<<(j+1)<<(i+1)<<" = global_params"<<(j+1)<<"[atom"<<(i+1)<<"];\n";
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
            verifyCutoff<<"real4 pos"<<(i+1)<<" = posq[p"<<(i+1)<<"];\n";
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
    string computeTypeIndex = "particleTypes[p"+cl.intToString(particlesPerSet)+"]";
    for (int i = particlesPerSet-2; i >= 0; i--)
        computeTypeIndex = "particleTypes[p"+cl.intToString(i+1)+"]+"+cl.intToString(numTypes)+"*("+computeTypeIndex+")";
    
    // Create replacements for extra arguments.
    
    stringstream extraArgs;
    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", __global const float* globals";
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        OpenCLNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArgs<<", __global const "<<buffer.getType()<<"* restrict global_params"<<(i+1);
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
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    defines["M_PI"] = cl.doubleToString(M_PI);
    defines["CUTOFF_SQUARED"] = cl.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    defines["TILE_SIZE"] = cl.intToString(OpenCLContext::TileSize);
    defines["NUM_BLOCKS"] = cl.intToString(cl.getNumAtomBlocks());
    defines["FIND_NEIGHBORS_WORKGROUP_SIZE"] = cl.intToString(findNeighborsWorkgroupSize);
    cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customManyParticle, replacements), defines);
    forceKernel = cl::Kernel(program, "computeInteraction");
    blockBoundsKernel = cl::Kernel(program, "findBlockBounds");
    neighborsKernel = cl::Kernel(program, "findNeighbors");
    startIndicesKernel = cl::Kernel(program, "computeNeighborStartIndices");
    copyPairsKernel = cl::Kernel(program, "copyPairsToNeighborList");
}

double OpenCLCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        
        // Set arguments for the force kernel.
        
        int index = 0;
        forceKernel.setArg<cl::Buffer>(index++, cl.getLongForceBuffer().getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        setPeriodicBoxArgs(cl, forceKernel, index);
        index += 5;
        if (nonbondedMethod != NoCutoff) {
            forceKernel.setArg<cl::Buffer>(index++, neighbors.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(index++, neighborStartIndex.getDeviceBuffer());
        }
        if (particleTypes.isInitialized()) {
            forceKernel.setArg<cl::Buffer>(index++, particleTypes.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(index++, orderIndex.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(index++, particleOrder.getDeviceBuffer());
        }
        if (exclusions.isInitialized()) {
            forceKernel.setArg<cl::Buffer>(index++, exclusions.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(index++, exclusionStartIndex.getDeviceBuffer());
        }
        if (globals.isInitialized())
            forceKernel.setArg<cl::Buffer>(index++, globals.getDeviceBuffer());
        for (auto& buffer : params->getBuffers())
            forceKernel.setArg<cl::Memory>(index++, buffer.getMemory());
        for (auto& function : tabulatedFunctions)
            forceKernel.setArg<cl::Buffer>(index++, function.getDeviceBuffer());
        
        if (nonbondedMethod != NoCutoff) {
            // Set arguments for the block bounds kernel.

            index = 0;
            setPeriodicBoxArgs(cl, blockBoundsKernel, index);
            index += 5;
            blockBoundsKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
            blockBoundsKernel.setArg<cl::Buffer>(index++, blockCenter.getDeviceBuffer());
            blockBoundsKernel.setArg<cl::Buffer>(index++, blockBoundingBox.getDeviceBuffer());
            blockBoundsKernel.setArg<cl::Buffer>(index++, numNeighborPairs.getDeviceBuffer());

            // Set arguments for the neighbor list kernel.

            index = 0;
            setPeriodicBoxArgs(cl, neighborsKernel, index);
            index += 5;
            neighborsKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
            neighborsKernel.setArg<cl::Buffer>(index++, blockCenter.getDeviceBuffer());
            neighborsKernel.setArg<cl::Buffer>(index++, blockBoundingBox.getDeviceBuffer());
            neighborsKernel.setArg<cl::Buffer>(index++, neighborPairs.getDeviceBuffer());
            neighborsKernel.setArg<cl::Buffer>(index++, numNeighborPairs.getDeviceBuffer());
            neighborsKernel.setArg<cl::Buffer>(index++, numNeighborsForAtom.getDeviceBuffer());
            index++;
            if (exclusions.isInitialized()) {
                neighborsKernel.setArg<cl::Buffer>(index++, exclusions.getDeviceBuffer());
                neighborsKernel.setArg<cl::Buffer>(index++, exclusionStartIndex.getDeviceBuffer());
            }
            
            // Set arguments for the kernel to find neighbor list start indices.
            
            index = 0;
            startIndicesKernel.setArg<cl::Buffer>(index++, numNeighborsForAtom.getDeviceBuffer());
            startIndicesKernel.setArg<cl::Buffer>(index++, neighborStartIndex.getDeviceBuffer());
            startIndicesKernel.setArg<cl::Buffer>(index++, numNeighborPairs.getDeviceBuffer());

            // Set arguments for the kernel to assemble the final neighbor list.
            
            index = 0;
            copyPairsKernel.setArg<cl::Buffer>(index++, neighborPairs.getDeviceBuffer());
            copyPairsKernel.setArg<cl::Buffer>(index++, neighbors.getDeviceBuffer());
            copyPairsKernel.setArg<cl::Buffer>(index++, numNeighborPairs.getDeviceBuffer());
            index++;
            copyPairsKernel.setArg<cl::Buffer>(index++, numNeighborsForAtom.getDeviceBuffer());
            copyPairsKernel.setArg<cl::Buffer>(index++, neighborStartIndex.getDeviceBuffer());
       }
    }
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
    while (true) {
        int* numPairs = (int*) cl.getPinnedBuffer();
        cl::Event event;
        if (nonbondedMethod != NoCutoff) {
            neighborsKernel.setArg<int>(11, maxNeighborPairs);
            startIndicesKernel.setArg<int>(3, maxNeighborPairs);
            copyPairsKernel.setArg<int>(3, maxNeighborPairs);
            cl.executeKernel(blockBoundsKernel, cl.getNumAtomBlocks());
            cl.executeKernel(neighborsKernel, cl.getNumAtoms(), findNeighborsWorkgroupSize);

            // We need to make sure there was enough memory for the neighbor list.  Download the
            // information asynchronously so kernels can be running at the same time.

            numNeighborPairs.download(numPairs, false);
            cl.getQueue().enqueueMarker(&event);
            cl.executeKernel(startIndicesKernel, 256, 256);
            cl.executeKernel(copyPairsKernel, maxNeighborPairs);
        }
        int maxThreads = min(cl.getNumAtoms()*forceWorkgroupSize, cl.getEnergyBuffer().getSize());
        cl.executeKernel(forceKernel, maxThreads, forceWorkgroupSize);
        if (nonbondedMethod != NoCutoff) {
            // Make sure there was enough memory for the neighbor list.

            event.wait();
            if (*numPairs > maxNeighborPairs) {
                // Resize the arrays and run the calculation again.

                maxNeighborPairs = (int) (1.1*(*numPairs));
                neighborPairs.resize(maxNeighborPairs);
                neighbors.resize(maxNeighborPairs);
                forceKernel.setArg<cl::Buffer>(8, neighbors.getDeviceBuffer());
                neighborsKernel.setArg<cl::Buffer>(8, neighborPairs.getDeviceBuffer());
                copyPairsKernel.setArg<cl::Buffer>(0, neighborPairs.getDeviceBuffer());
                copyPairsKernel.setArg<cl::Buffer>(1, neighbors.getDeviceBuffer());
                continue;
            }
        }
        break;
    }
    return 0.0;
}

void OpenCLCalcCustomManyParticleForceKernel::copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) {
    int numParticles = force.getNumParticles();
    if (numParticles != cl.getNumAtoms())
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
    
    cl.invalidateMolecules(info);
}

class OpenCLCalcGayBerneForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(int requiredBuffers, const GayBerneForce& force) : OpenCLForceInfo(requiredBuffers), force(force) {
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

class OpenCLCalcGayBerneForceKernel::ReorderListener : public OpenCLContext::ReorderListener {
public:
    ReorderListener(OpenCLCalcGayBerneForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.sortAtoms();
    }
private:
    OpenCLCalcGayBerneForceKernel& owner;
};

void OpenCLCalcGayBerneForceKernel::initialize(const System& system, const GayBerneForce& force) {
    if (!cl.getSupports64BitGlobalAtomics())
        throw OpenMMException("GayBerneForce requires a device that supports 64 bit atomic operations");

    // Initialize interactions.

    int numParticles = force.getNumParticles();
    sigParams.initialize<mm_float4>(cl, cl.getPaddedNumAtoms(), "sigParams");
    epsParams.initialize<mm_float2>(cl, cl.getPaddedNumAtoms(), "epsParams");
    scale.initialize<mm_float4>(cl, cl.getPaddedNumAtoms(), "scale");
    axisParticleIndices.initialize<mm_int2>(cl, cl.getPaddedNumAtoms(), "axisParticleIndices");
    sortedParticles.initialize<cl_int>(cl, cl.getPaddedNumAtoms(), "sortedParticles");
    aMatrix.initialize<cl_float>(cl, 9*cl.getPaddedNumAtoms(), "aMatrix");
    bMatrix.initialize<cl_float>(cl, 9*cl.getPaddedNumAtoms(), "bMatrix");
    gMatrix.initialize<cl_float>(cl, 9*cl.getPaddedNumAtoms(), "gMatrix");
    vector<mm_float4> sigParamsVector(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_float2> epsParamsVector(cl.getPaddedNumAtoms(), mm_float2(0, 0));
    vector<mm_float4> scaleVector(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_int2> axisParticleVector(cl.getPaddedNumAtoms(), mm_int2(0, 0));
    isRealParticle.resize(cl.getPaddedNumAtoms());
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
    exclusions.initialize<cl_int>(cl, max(1, (int) excludedPairs.size()), "exclusions");
    exclusionStartIndex.initialize<cl_int>(cl, numRealParticles+1, "exclusionStartIndex");
    exceptionParticles.initialize<mm_int4>(cl, max(1, numExceptions), "exceptionParticles");
    exceptionParams.initialize<mm_float2>(cl, max(1, numExceptions), "exceptionParams");
    if (numExceptions > 0)
        exceptionParams.upload(exceptionParamsVec);
    
    // Create data structures used for the neighbor list.

    int numAtomBlocks = (numRealParticles+31)/32;
    int elementSize = (cl.getUseDoublePrecision() ? sizeof(cl_double) : sizeof(cl_float));
    blockCenter.initialize(cl, numAtomBlocks, 4*elementSize, "blockCenter");
    blockBoundingBox.initialize(cl, numAtomBlocks, 4*elementSize, "blockBoundingBox");
    sortedPos.initialize(cl, numRealParticles, 4*elementSize, "sortedPos");
    maxNeighborBlocks = numRealParticles*2;
    neighbors.initialize<cl_int>(cl, maxNeighborBlocks*32, "neighbors");
    neighborIndex.initialize<cl_int>(cl, maxNeighborBlocks, "neighborIndex");
    neighborBlockCount.initialize<cl_int>(cl, 1, "neighborBlockCount");

    // Create array for accumulating torques.
    
    torque.initialize<cl_long>(cl, 3*cl.getPaddedNumAtoms(), "torque");
    cl.addAutoclearBuffer(torque);

    // Create the kernels.
    
    nonbondedMethod = force.getNonbondedMethod();
    bool useCutoff = (nonbondedMethod != GayBerneForce::NoCutoff);
    bool usePeriodic = (nonbondedMethod == GayBerneForce::CutoffPeriodic);
    map<string, string> defines;
    defines["USE_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cl.doubleToString(cutoff*cutoff);
    if (useCutoff) {
        defines["USE_CUTOFF"] = 1;
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        
        // Compute the switching coefficients.
        
        if (force.getUseSwitchingFunction()) {
            defines["SWITCH_CUTOFF"] = cl.doubleToString(force.getSwitchingDistance());
            defines["SWITCH_C3"] = cl.doubleToString(10/pow(force.getSwitchingDistance()-cutoff, 3.0));
            defines["SWITCH_C4"] = cl.doubleToString(15/pow(force.getSwitchingDistance()-cutoff, 4.0));
            defines["SWITCH_C5"] = cl.doubleToString(6/pow(force.getSwitchingDistance()-cutoff, 5.0));
        }
    }
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::gayBerne, defines);
    framesKernel = cl::Kernel(program, "computeEllipsoidFrames");
    blockBoundsKernel = cl::Kernel(program, "findBlockBounds");
    neighborsKernel = cl::Kernel(program, "findNeighbors");
    forceKernel = cl::Kernel(program, "computeForce");
    torqueKernel = cl::Kernel(program, "applyTorques");
    info = new ForceInfo(cl.getNonbondedUtilities().getNumForceBuffers(), force);
    cl.addForce(info);
    cl.addReorderListener(new ReorderListener(*this));
}

double OpenCLCalcGayBerneForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        sortAtoms();
        framesKernel.setArg<cl_int>(0, numRealParticles);
        framesKernel.setArg<cl::Buffer>(1, cl.getPosq().getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(2, axisParticleIndices.getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(3, sigParams.getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(4, scale.getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(5, aMatrix.getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(6, bMatrix.getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(7, gMatrix.getDeviceBuffer());
        framesKernel.setArg<cl::Buffer>(8, sortedParticles.getDeviceBuffer());
        blockBoundsKernel.setArg<cl_int>(0, numRealParticles);
        blockBoundsKernel.setArg<cl::Buffer>(6, sortedParticles.getDeviceBuffer());
        blockBoundsKernel.setArg<cl::Buffer>(7, cl.getPosq().getDeviceBuffer());
        blockBoundsKernel.setArg<cl::Buffer>(8, sortedPos.getDeviceBuffer());
        blockBoundsKernel.setArg<cl::Buffer>(9, blockCenter.getDeviceBuffer());
        blockBoundsKernel.setArg<cl::Buffer>(10, blockBoundingBox.getDeviceBuffer());
        blockBoundsKernel.setArg<cl::Buffer>(11, neighborBlockCount.getDeviceBuffer());
        neighborsKernel.setArg<cl_int>(0, numRealParticles);
        neighborsKernel.setArg<cl_int>(1, maxNeighborBlocks);
        neighborsKernel.setArg<cl::Buffer>(7, sortedPos.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(8, blockCenter.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(9, blockBoundingBox.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(10, neighbors.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(11, neighborIndex.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(12, neighborBlockCount.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(13, exclusions.getDeviceBuffer());
        neighborsKernel.setArg<cl::Buffer>(14, exclusionStartIndex.getDeviceBuffer());
        int index = 0;
        forceKernel.setArg<cl::Buffer>(index++, cl.getLongForceBuffer().getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, torque.getDeviceBuffer());
        forceKernel.setArg<cl_int>(index++, numRealParticles);
        forceKernel.setArg<cl_int>(index++, exceptionAtoms.size());
        forceKernel.setArg<cl::Buffer>(index++, cl.getEnergyBuffer().getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, sortedPos.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, sigParams.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, epsParams.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, sortedParticles.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, aMatrix.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, bMatrix.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, gMatrix.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, exclusions.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, exclusionStartIndex.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, exceptionParticles.getDeviceBuffer());
        forceKernel.setArg<cl::Buffer>(index++, exceptionParams.getDeviceBuffer());
        if (nonbondedMethod != GayBerneForce::NoCutoff) {
            forceKernel.setArg<cl_int>(index++, maxNeighborBlocks);
            forceKernel.setArg<cl::Buffer>(index++, neighbors.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(index++, neighborIndex.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(index++, neighborBlockCount.getDeviceBuffer());
        }
        index = 0;
        torqueKernel.setArg<cl::Buffer>(index++, cl.getLongForceBuffer().getDeviceBuffer());
        torqueKernel.setArg<cl::Buffer>(index++, torque.getDeviceBuffer());
        torqueKernel.setArg<cl_int>(index++, numRealParticles);
        torqueKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        torqueKernel.setArg<cl::Buffer>(index++, axisParticleIndices.getDeviceBuffer());
        torqueKernel.setArg<cl::Buffer>(index++, sortedParticles.getDeviceBuffer());
    }
    cl.executeKernel(framesKernel, numRealParticles);
    setPeriodicBoxArgs(cl, blockBoundsKernel, 1);
    cl.executeKernel(blockBoundsKernel, (numRealParticles+31)/32);
    if (nonbondedMethod == GayBerneForce::NoCutoff) {
        cl.executeKernel(forceKernel, cl.getNonbondedUtilities().getNumForceThreadBlocks()*cl.getNonbondedUtilities().getForceThreadBlockSize());
    }
    else {
        while (true) {
            setPeriodicBoxArgs(cl, neighborsKernel, 2);
            cl.executeKernel(neighborsKernel, numRealParticles);
            cl_int* count = (cl_int*) cl.getPinnedBuffer();
            cl::Event event;
            cl.getQueue().enqueueReadBuffer(neighborBlockCount.getDeviceBuffer(), CL_FALSE, 0, neighborBlockCount.getSize()*neighborBlockCount.getElementSize(), count, NULL, &event);
            setPeriodicBoxArgs(cl, forceKernel, 20);
            cl.executeKernel(forceKernel, cl.getNonbondedUtilities().getNumForceThreadBlocks()*cl.getNonbondedUtilities().getForceThreadBlockSize());
            event.wait();
            if (*count <= maxNeighborBlocks)
                break;
            
            // There wasn't enough room for the neighbor list, so we need to recreate it.

            maxNeighborBlocks = (int) ceil((*count)*1.1);
            neighbors.resize(maxNeighborBlocks*32);
            neighborIndex.resize(maxNeighborBlocks);
            neighborsKernel.setArg<cl::Buffer>(10, neighbors.getDeviceBuffer());
            neighborsKernel.setArg<cl::Buffer>(11, neighborIndex.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(17, neighbors.getDeviceBuffer());
            forceKernel.setArg<cl::Buffer>(18, neighborIndex.getDeviceBuffer());
        }
    }
    cl.executeKernel(torqueKernel, numRealParticles);
    return 0.0;
}

void OpenCLCalcGayBerneForceKernel::copyParametersToContext(ContextImpl& context, const GayBerneForce& force) {
    // Make sure the new parameters are acceptable.
    
    if (force.getNumParticles() != cl.getNumAtoms())
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
    
    vector<mm_float4> sigParamsVector(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_float2> epsParamsVector(cl.getPaddedNumAtoms(), mm_float2(0, 0));
    vector<mm_float4> scaleVector(cl.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
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
    cl.invalidateMolecules(info);
    sortAtoms();
}

void OpenCLCalcGayBerneForceKernel::sortAtoms() {
    // Sort the list of atoms by type to avoid thread divergence.  This is executed every time
    // the atoms are reordered.
    
    int nextIndex = 0;
    vector<cl_int> particles(cl.getPaddedNumAtoms(), 0);
    const vector<int>& order = cl.getAtomIndex();
    vector<int> inverseOrder(order.size(), -1);
    for (int i = 0; i < cl.getNumAtoms(); i++) {
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

class OpenCLCalcCustomCVForceKernel::ForceInfo : public OpenCLForceInfo {
public:
    ForceInfo(OpenCLForceInfo& force) : OpenCLForceInfo(0), force(force) {
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
    OpenCLForceInfo& force;
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

    for (OpenCLForceInfo* info : cl2.getForceInfos())
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
    referencePos.upload(pos, true, true);

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

OpenCLIntegrateVerletStepKernel::~OpenCLIntegrateVerletStepKernel() {
}

void OpenCLIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl::Program program = cl.createProgram(OpenCLKernelSources::verlet, "");
    kernel1 = cl::Kernel(program, "integrateVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVerletPart2");
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
        setPosqCorrectionArg(cl, kernel1, 3);
        kernel1.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(6, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl_int>(0, numAtoms);
        kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 3);
        kernel2.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
    }
    cl.getIntegrationUtilities().setNextStepSize(dt);

    // Call the first integration kernel.

    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void OpenCLIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::langevin, defines, "");
    kernel1 = cl::Kernel(program, "integrateLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateLangevinPart2");
    params.initialize(cl, 3, cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(cl_double) : sizeof(cl_float), "langevinParams");
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
        kernel1.setArg<cl::Buffer>(3, params.getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 1);
        kernel2.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
    }
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    cl.getIntegrationUtilities().setNextStepSize(stepSize);
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize*friction);
        double fscale = (friction == 0 ? stepSize : (1-vscale)/friction);
        double noisescale = sqrt(kT*(1-vscale*vscale));
        vector<cl_double> p(params.getSize());
        p[0] = vscale;
        p[1] = fscale;
        p[2] = noisescale;
        params.upload(p, true, true);
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
    integration.computeVirtualSites();

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

OpenCLIntegrateBrownianStepKernel::~OpenCLIntegrateBrownianStepKernel() {
}

void OpenCLIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::brownian, defines, "");
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
        setPosqCorrectionArg(cl, kernel2, 2);
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getPosDelta().getDeviceBuffer());
    }
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
            kernel1.setArg<cl_double>(0, tau*stepSize);
            kernel1.setArg<cl_double>(1, sqrt(2.0f*BOLTZ*temperature*stepSize*tau));
            kernel2.setArg<cl_double>(0, 1.0/stepSize);
        }
        else {
            kernel1.setArg<cl_float>(0, (cl_float) (tau*stepSize));
            kernel1.setArg<cl_float>(1, (cl_float) (sqrt(2.0f*BOLTZ*temperature*stepSize*tau)));
            kernel2.setArg<cl_float>(0, (cl_float) (1.0/stepSize));
        }
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
    integration.computeVirtualSites();

    // Update the time and step count.

    cl.setTime(cl.getTime()+stepSize);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateBrownianStepKernel::computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0);
}

OpenCLIntegrateVariableVerletStepKernel::~OpenCLIntegrateVariableVerletStepKernel() {
}

void OpenCLIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl::Program program = cl.createProgram(OpenCLKernelSources::verlet, "");
    kernel1 = cl::Kernel(program, "integrateVerletPart1");
    kernel2 = cl::Kernel(program, "integrateVerletPart2");
    selectSizeKernel = cl::Kernel(program, "selectVerletStepSize");
    blockSize = min(min(256, system.getNumParticles()), (int) selectSizeKernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(cl.getDevice()));
}

double OpenCLIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl_int>(0, numAtoms);
        kernel1.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel1, 3);
        kernel1.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(6, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl_int>(0, numAtoms);
        kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 3);
        kernel2.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
        selectSizeKernel.setArg<cl_int>(0, numAtoms);
        selectSizeKernel.setArg<cl::Buffer>(3, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(5, cl.getForce().getDeviceBuffer());
        int elementSize = (useDouble ? sizeof(cl_double) : sizeof(cl_float));
        selectSizeKernel.setArg(6, blockSize*elementSize, NULL);
    }

    // Select the step size to use.

    double maxStepSize = maxTime-cl.getTime();
    float maxStepSizeFloat = (float) maxStepSize;
    if (useDouble) {
        selectSizeKernel.setArg<cl_double>(1, maxStepSize);
        selectSizeKernel.setArg<cl_double>(2, integrator.getErrorTolerance());
    }
    else {
        selectSizeKernel.setArg<cl_float>(1, maxStepSizeFloat);
        selectSizeKernel.setArg<cl_float>(2, (cl_float) integrator.getErrorTolerance());
    }
    cl.executeKernel(selectSizeKernel, blockSize, blockSize);

    // Call the first integration kernel.

    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);
    integration.computeVirtualSites();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif

    // Update the time and step count.

    double dt = cl.getIntegrationUtilities().getLastStepSize();
    double time = cl.getTime()+dt;
    if (useDouble) {
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
        if (dt == maxStepSizeFloat)
            time = maxTime; // Avoid round-off error
    }
    cl.setTime(time);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    return dt;
}

double OpenCLIntegrateVariableVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void OpenCLIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLKernelSources::langevin, defines, "");
    kernel1 = cl::Kernel(program, "integrateLangevinPart1");
    kernel2 = cl::Kernel(program, "integrateLangevinPart2");
    selectSizeKernel = cl::Kernel(program, "selectLangevinStepSize");
    params.initialize(cl, 3, cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(cl_double) : sizeof(cl_float), "langevinParams");
    blockSize = min(256, system.getNumParticles());
    blockSize = max(blockSize, params.getSize());
    blockSize = min(blockSize, (int) selectSizeKernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE>(cl.getDevice()));
}

double OpenCLIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl::Buffer>(0, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(1, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(3, params.getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, integration.getRandom().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(0, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 1);
        kernel2.setArg<cl::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(3, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(4, integration.getStepSize().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(5, cl.getVelm().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(6, cl.getForce().getDeviceBuffer());
        selectSizeKernel.setArg<cl::Buffer>(7, params.getDeviceBuffer());
        int elementSize = (useDouble ? sizeof(cl_double) : sizeof(cl_float));
        selectSizeKernel.setArg(8, params.getSize()*elementSize, NULL);
        selectSizeKernel.setArg(9, blockSize*elementSize, NULL);
    }

    // Select the step size to use.

    double maxStepSize = maxTime-cl.getTime();
    float maxStepSizeFloat = (float) maxStepSize;
    if (useDouble) {
        selectSizeKernel.setArg<cl_double>(0, maxStepSize);
        selectSizeKernel.setArg<cl_double>(1, integrator.getErrorTolerance());
        selectSizeKernel.setArg<cl_double>(2, integrator.getFriction());
        selectSizeKernel.setArg<cl_double>(3, BOLTZ*integrator.getTemperature());
    }
    else {
        selectSizeKernel.setArg<cl_float>(0, maxStepSizeFloat);
        selectSizeKernel.setArg<cl_float>(1, (cl_float) integrator.getErrorTolerance());
        selectSizeKernel.setArg<cl_float>(2, (cl_float) integrator.getFriction());
        selectSizeKernel.setArg<cl_float>(3, (cl_float) (BOLTZ*integrator.getTemperature()));
    }
    cl.executeKernel(selectSizeKernel, blockSize, blockSize);

    // Call the first integration kernel.

    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cl.getPaddedNumAtoms()));
    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);
    integration.computeVirtualSites();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif

    // Update the time and step count.

    double dt = cl.getIntegrationUtilities().getLastStepSize();
    double time = cl.getTime()+dt;
    if (useDouble) {
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
        if (dt == maxStepSizeFloat)
            time = maxTime; // Avoid round-off error
    }
    cl.setTime(time);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    return dt;
}

double OpenCLIntegrateVariableLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

class OpenCLIntegrateCustomStepKernel::ReorderListener : public OpenCLContext::ReorderListener {
public:
    ReorderListener(OpenCLContext& cl, vector<OpenCLArray>& perDofValues, vector<vector<mm_float4> >& localPerDofValuesFloat, vector<vector<mm_double4> >& localPerDofValuesDouble, vector<bool>& deviceValuesAreCurrent) :
            cl(cl), perDofValues(perDofValues), localPerDofValuesFloat(localPerDofValuesFloat), localPerDofValuesDouble(localPerDofValuesDouble), deviceValuesAreCurrent(deviceValuesAreCurrent) {
        int numAtoms = cl.getNumAtoms();
        lastAtomOrder.resize(numAtoms);
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = cl.getAtomIndex()[i];
    }
    void execute() {
        // Reorder the per-DOF variables to reflect the new atom order.

        if (perDofValues.size() == 0)
            return;
        int numAtoms = cl.getNumAtoms();
        const vector<int>& order = cl.getAtomIndex();
        for (int index = 0; index < perDofValues.size(); index++) {
            if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
                if (deviceValuesAreCurrent[index])
                    perDofValues[index].download(localPerDofValuesDouble[index]);
                vector<mm_double4> swap(numAtoms);
                for (int i = 0; i < numAtoms; i++)
                    swap[lastAtomOrder[i]] = localPerDofValuesDouble[index][i];
                for (int i = 0; i < numAtoms; i++)
                    localPerDofValuesDouble[index][i] = swap[order[i]];
                perDofValues[index].upload(localPerDofValuesDouble[index]);
            }
            else {
                if (deviceValuesAreCurrent[index])
                    perDofValues[index].download(localPerDofValuesFloat[index]);
                vector<mm_float4> swap(numAtoms);
                for (int i = 0; i < numAtoms; i++)
                    swap[lastAtomOrder[i]] = localPerDofValuesFloat[index][i];
                for (int i = 0; i < numAtoms; i++)
                    localPerDofValuesFloat[index][i] = swap[order[i]];
                perDofValues[index].upload(localPerDofValuesFloat[index]);
            }
            deviceValuesAreCurrent[index] = true;
        }
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = order[i];
    }
private:
    OpenCLContext& cl;
    vector<OpenCLArray>& perDofValues;
    vector<vector<mm_float4> >& localPerDofValuesFloat;
    vector<vector<mm_double4> >& localPerDofValuesDouble;
    vector<bool>& deviceValuesAreCurrent;
    vector<int> lastAtomOrder;
};

class OpenCLIntegrateCustomStepKernel::DerivFunction : public CustomFunction {
public:
    DerivFunction(map<string, double>& energyParamDerivs, const string& param) : energyParamDerivs(energyParamDerivs), param(param) {
    }
    int getNumArguments() const {
        return 0;
    }
    double evaluate(const double* arguments) const {
        return energyParamDerivs[param];
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        return 0;
    }
    CustomFunction* clone() const {
        return new DerivFunction(energyParamDerivs, param);
    }
private:
    map<string, double>& energyParamDerivs;
    string param;
};

void OpenCLIntegrateCustomStepKernel::initialize(const System& system, const CustomIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    numGlobalVariables = integrator.getNumGlobalVariables();
    int elementSize = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    sumBuffer.initialize(cl, system.getNumParticles(), elementSize, "sumBuffer");
    summedValue.initialize(cl, 1, elementSize, "summedValue");
    perDofValues.resize(integrator.getNumPerDofVariables());
    localPerDofValuesFloat.resize(perDofValues.size());
    localPerDofValuesDouble.resize(perDofValues.size());
    for (int i = 0; i < perDofValues.size(); i++)
        perDofValues[i].initialize(cl, system.getNumParticles(), 4*elementSize, "perDofVariables");
    localValuesAreCurrent.resize(integrator.getNumPerDofVariables(), false);
    deviceValuesAreCurrent.resize(integrator.getNumPerDofVariables(), false);
    cl.addReorderListener(new ReorderListener(cl, perDofValues, localPerDofValuesFloat, localPerDofValuesDouble, deviceValuesAreCurrent));
    SimTKOpenMMUtilities::setRandomNumberSeed(integrator.getRandomNumberSeed());
}

string OpenCLIntegrateCustomStepKernel::createPerDofComputation(const string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator,
        const string& forceName, const string& energyName, vector<const TabulatedFunction*>& functions, vector<pair<string, string> >& functionNames) {
    string tempType = (cl.getSupportsDoublePrecision() ? "double3" : "float3");
    string convert = (cl.getSupportsDoublePrecision() ? "convert_double3" : "");
    map<string, Lepton::ParsedExpression> expressions;
    expressions[tempType+" tempResult = "] = expr;
    map<string, string> variables;
    variables["x"] = convert+"(position.xyz)";
    variables["v"] = convert+"(velocity.xyz)";
    variables[forceName] = convert+"(f.xyz)";
    variables["gaussian"] = convert+"(gaussian.xyz)";
    variables["uniform"] = convert+"(uniform.xyz)";
    variables["m"] = "mass";
    variables["dt"] = "stepSize";
    if (energyName != "")
        variables[energyName] = "energy";
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        variables[integrator.getGlobalVariableName(i)] = "globals["+cl.intToString(globalVariableIndex[i])+"]";
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        variables[integrator.getPerDofVariableName(i)] = convert+"(perDof"+cl.intToString(i)+")";
    for (int i = 0; i < (int) parameterNames.size(); i++)
        variables[parameterNames[i]] = "globals["+cl.intToString(parameterVariableIndex[i])+"]";
    vector<pair<ExpressionTreeNode, string> > variableNodes;
    findExpressionsForDerivs(expr.getRootNode(), variableNodes);
    for (auto& var : variables)
        variableNodes.push_back(make_pair(ExpressionTreeNode(new Operation::Variable(var.first)), var.second));
    string result = cl.getExpressionUtilities().createExpressions(expressions, variableNodes, functions, functionNames, "temp", tempType);
    if (variable == "x")
        result += "position.x = tempResult.x; position.y = tempResult.y; position.z = tempResult.z;\n";
    else if (variable == "v")
        result += "velocity.x = tempResult.x; velocity.y = tempResult.y; velocity.z = tempResult.z;\n";
    else if (variable == "")
        result += "sum[index] = tempResult.x+tempResult.y+tempResult.z;\n";
    else {
        for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
            if (variable == integrator.getPerDofVariableName(i)) {
                string varName = "perDof"+cl.intToString(i);
                result += varName+".x = tempResult.x; "+varName+".y = tempResult.y; "+varName+".z = tempResult.z;\n";
            }
    }
    return result;
}

void OpenCLIntegrateCustomStepKernel::prepareForComputation(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    bool useDouble = cl.getUseDoublePrecision() || cl.getUseMixedPrecision();
    string tempType = (cl.getSupportsDoublePrecision() ? "double3" : "float3");
    string perDofType = (useDouble ? "double4" : "float4");
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // Initialize various data structures.
        
        const map<string, double>& params = context.getParameters();
        for (auto& param : params)
            parameterNames.push_back(param.first);
        kernels.resize(integrator.getNumComputations());
        requiredGaussian.resize(integrator.getNumComputations(), 0);
        requiredUniform.resize(integrator.getNumComputations(), 0);
        needsGlobals.resize(numSteps, false);
        globalExpressions.resize(numSteps);
        stepType.resize(numSteps);
        stepTarget.resize(numSteps);
        merged.resize(numSteps, false);
        modifiesParameters = false;
        sumWorkGroupSize = cl.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
        if (sumWorkGroupSize > 512)
            sumWorkGroupSize = 512;
        map<string, string> defines;
        defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
        defines["WORK_GROUP_SIZE"] = cl.intToString(sumWorkGroupSize);

        // Record the tabulated functions.

        map<string, Lepton::CustomFunction*> functions;
        vector<pair<string, string> > functionNames;
        vector<const TabulatedFunction*> functionList;
        vector<string> tableTypes;
        tabulatedFunctions.resize(integrator.getNumTabulatedFunctions());
        for (int i = 0; i < integrator.getNumTabulatedFunctions(); i++) {
            functionList.push_back(&integrator.getTabulatedFunction(i));
            string name = integrator.getTabulatedFunctionName(i);
            string arrayName = "table"+cl.intToString(i);
            functionNames.push_back(make_pair(name, arrayName));
            functions[name] = createReferenceTabulatedFunction(integrator.getTabulatedFunction(i));
            int width;
            vector<float> f = cl.getExpressionUtilities().computeFunctionCoefficients(integrator.getTabulatedFunction(i), width);
            tabulatedFunctions[i].initialize<float>(cl, f.size(), "TabulatedFunction");
            tabulatedFunctions[i].upload(f);
            if (width == 1)
                tableTypes.push_back("float");
            else
                tableTypes.push_back("float"+cl.intToString(width));
        }

        // Record information about all the computation steps.

        vector<string> variable(numSteps);
        vector<int> forceGroup;
        vector<vector<Lepton::ParsedExpression> > expression;
        CustomIntegratorUtilities::analyzeComputations(context, integrator, expression, comparisons, blockEnd, invalidatesForces, needsForces, needsEnergy, computeBothForceAndEnergy, forceGroup, functions);
        for (int step = 0; step < numSteps; step++) {
            string expr;
            integrator.getComputationStep(step, stepType[step], variable[step], expr);
            if (stepType[step] == CustomIntegrator::WhileBlockStart)
                blockEnd[blockEnd[step]] = step; // Record where to branch back to.
            if (stepType[step] == CustomIntegrator::ComputeGlobal || stepType[step] == CustomIntegrator::IfBlockStart || stepType[step] == CustomIntegrator::WhileBlockStart)
                for (auto& expr : expression[step])
                    globalExpressions[step].push_back(ParsedExpression(replaceDerivFunctions(expr.getRootNode(), context)).createCompiledExpression());
        }
        for (int step = 0; step < numSteps; step++) {
            for (auto& expr : globalExpressions[step])
                expressionSet.registerExpression(expr);
        }
        
        // Record the indices for variables in the CompiledExpressionSet.
        
        gaussianVariableIndex = expressionSet.getVariableIndex("gaussian");
        uniformVariableIndex = expressionSet.getVariableIndex("uniform");
        dtVariableIndex = expressionSet.getVariableIndex("dt");
        for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
            globalVariableIndex.push_back(expressionSet.getVariableIndex(integrator.getGlobalVariableName(i)));
        for (auto& name : parameterNames)
            parameterVariableIndex.push_back(expressionSet.getVariableIndex(name));

        // Record the variable names and flags for the force and energy in each step.

        forceGroupFlags.resize(numSteps, -1);
        vector<string> forceGroupName;
        vector<string> energyGroupName;
        for (int i = 0; i < 32; i++) {
            stringstream fname;
            fname << "f" << i;
            forceGroupName.push_back(fname.str());
            stringstream ename;
            ename << "energy" << i;
            energyGroupName.push_back(ename.str());
        }
        vector<string> forceName(numSteps, "f");
        vector<string> energyName(numSteps, "energy");
        stepEnergyVariableIndex.resize(numSteps, expressionSet.getVariableIndex("energy"));
        for (int step = 0; step < numSteps; step++) {
            if (needsForces[step] && forceGroup[step] > -1)
                forceName[step] = forceGroupName[forceGroup[step]];
            if (needsEnergy[step] && forceGroup[step] > -1) {
                energyName[step] = energyGroupName[forceGroup[step]];
                stepEnergyVariableIndex[step] = expressionSet.getVariableIndex(energyName[step]);
            }
            if (forceGroup[step] > -1)
                forceGroupFlags[step] = 1<<forceGroup[step];
            if (forceGroupFlags[step] == -2 && step > 0)
                forceGroupFlags[step] = forceGroupFlags[step-1];
            if (forceGroupFlags[step] != -2 && savedForces.find(forceGroupFlags[step]) == savedForces.end()) {
                savedForces[forceGroupFlags[step]] = OpenCLArray();
                savedForces[forceGroupFlags[step]].initialize(cl, cl.getForce().getSize(), cl.getForce().getElementSize(), "savedForces");
            }
        }
        
        // Allocate space for storing global values, both on the host and the device.
        
        localGlobalValues.resize(expressionSet.getNumVariables());
        int elementSize = (cl.getUseDoublePrecision() || cl.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
        globalValues.initialize(cl, expressionSet.getNumVariables(), elementSize, "globalValues");
        for (int i = 0; i < integrator.getNumGlobalVariables(); i++) {
            localGlobalValues[globalVariableIndex[i]] = initialGlobalVariables[i];
            expressionSet.setVariable(globalVariableIndex[i], initialGlobalVariables[i]);
        }
        for (int i = 0; i < (int) parameterVariableIndex.size(); i++) {
            double value = context.getParameter(parameterNames[i]);
            localGlobalValues[parameterVariableIndex[i]] = value;
            expressionSet.setVariable(parameterVariableIndex[i], value);
        }
        int numContextParams = context.getParameters().size();
        localPerDofEnergyParamDerivs.resize(numContextParams);
        perDofEnergyParamDerivs.initialize(cl, max(1, numContextParams), elementSize, "perDofEnergyParamDerivs");
        
        // Record information about the targets of steps that will be stored in global variables.
        
        for (int step = 0; step < numSteps; step++) {
            if (stepType[step] == CustomIntegrator::ComputeGlobal || stepType[step] == CustomIntegrator::ComputeSum) {
                if (variable[step] == "dt")
                    stepTarget[step].type = DT;
                for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
                    if (variable[step] == integrator.getGlobalVariableName(i))
                        stepTarget[step].type = VARIABLE;
                for (auto& name : parameterNames)
                    if (variable[step] == name) {
                        stepTarget[step].type = PARAMETER;
                        modifiesParameters = true;
                    }
                stepTarget[step].variableIndex = expressionSet.getVariableIndex(variable[step]);
            }
        }

        // Identify which per-DOF steps are going to require global variables or context parameters.

        for (int step = 0; step < numSteps; step++) {
            if (stepType[step] == CustomIntegrator::ComputePerDof || stepType[step] == CustomIntegrator::ComputeSum) {
                for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
                    if (usesVariable(expression[step][0], integrator.getGlobalVariableName(i)))
                        needsGlobals[step] = true;
                for (auto& name : parameterNames)
                    if (usesVariable(expression[step][0], name))
                        needsGlobals[step] = true;
            }
        }
        
        // Determine how each step will represent the position (as just a value, or a value plus a delta).
        
        hasAnyConstraints = (context.getSystem().getNumConstraints() > 0);
        vector<bool> storePosAsDelta(numSteps, false);
        vector<bool> loadPosAsDelta(numSteps, false);
        if (hasAnyConstraints) {
            bool beforeConstrain = false;
            for (int step = numSteps-1; step >= 0; step--) {
                if (stepType[step] == CustomIntegrator::ConstrainPositions)
                    beforeConstrain = true;
                else if (stepType[step] == CustomIntegrator::ComputePerDof && variable[step] == "x" && beforeConstrain) {
                    storePosAsDelta[step] = true;
                    beforeConstrain = false;
                }
            }
            bool storedAsDelta = false;
            for (int step = 0; step < numSteps; step++) {
                loadPosAsDelta[step] = storedAsDelta;
                if (storePosAsDelta[step] == true)
                    storedAsDelta = true;
                if (stepType[step] == CustomIntegrator::ConstrainPositions)
                    storedAsDelta = false;
            }
        }
        
        // Identify steps that can be merged into a single kernel.
        
        for (int step = 1; step < numSteps; step++) {
            if (invalidatesForces[step-1] || forceGroupFlags[step] != forceGroupFlags[step-1])
                continue;
            if (stepType[step-1] == CustomIntegrator::ComputePerDof && stepType[step] == CustomIntegrator::ComputePerDof)
                merged[step] = true;
        }
        for (int step = numSteps-1; step > 0; step--)
            if (merged[step]) {
                needsForces[step-1] = (needsForces[step] || needsForces[step-1]);
                needsEnergy[step-1] = (needsEnergy[step] || needsEnergy[step-1]);
                needsGlobals[step-1] = (needsGlobals[step] || needsGlobals[step-1]);
                computeBothForceAndEnergy[step-1] = (computeBothForceAndEnergy[step] || computeBothForceAndEnergy[step-1]);
            }
        
        // Loop over all steps and create the kernels for them.
        
        for (int step = 0; step < numSteps; step++) {
            if ((stepType[step] == CustomIntegrator::ComputePerDof || stepType[step] == CustomIntegrator::ComputeSum) && !merged[step]) {
                // Compute a per-DOF value.
                
                stringstream compute;
                for (int i = 0; i < perDofValues.size(); i++)
                    compute << tempType<<" perDof"<<cl.intToString(i)<<" = convert_"<<tempType<<"(perDofValues"<<cl.intToString(i)<<"[index].xyz);\n";
                int numGaussian = 0, numUniform = 0;
                for (int j = step; j < numSteps && (j == step || merged[j]); j++) {
                    numGaussian += numAtoms*usesVariable(expression[j][0], "gaussian");
                    numUniform += numAtoms*usesVariable(expression[j][0], "uniform");
                    compute << "{\n";
                    if (numGaussian > 0)
                        compute << "float4 gaussian = gaussianValues[gaussianIndex+index];\n";
                    if (numUniform > 0)
                        compute << "float4 uniform = uniformValues[uniformIndex+index];\n";
                    compute << createPerDofComputation(stepType[j] == CustomIntegrator::ComputePerDof ? variable[j] : "", expression[j][0], integrator, forceName[j], energyName[j], functionList, functionNames);
                    if (variable[j] == "x") {
                        if (storePosAsDelta[j]) {
                            if (cl.getSupportsDoublePrecision())
                                compute << "posDelta[index] = convert_mixed4(convert_double4(position)-convert_double4(loadPos(posq, posqCorrection, index)));\n";
                            else
                                compute << "posDelta[index] = position-posq[index];\n";
                        }
                        else
                            compute << "storePos(posq, posqCorrection, index, position);\n";
                    }
                    else if (variable[j] == "v")
                        compute << "velm[index] = convert_mixed4(velocity);\n";
                    else {
                        for (int i = 0; i < perDofValues.size(); i++)
                            compute << "perDofValues"<<cl.intToString(i)<<"[index] = ("<<perDofType<<") (perDof"<<cl.intToString(i)<<".x, perDof"<<cl.intToString(i)<<".y, perDof"<<cl.intToString(i)<<".z, 0);\n";
                    }
                    if (numGaussian > 0)
                        compute << "gaussianIndex += NUM_ATOMS;\n";
                    if (numUniform > 0)
                        compute << "uniformIndex += NUM_ATOMS;\n";
                    compute << "}\n";
                }
                map<string, string> replacements;
                replacements["COMPUTE_STEP"] = compute.str();
                stringstream args;
                for (int i = 0; i < perDofValues.size(); i++) {
                    string valueName = "perDofValues"+cl.intToString(i);
                    args << ", __global " << perDofType << "* restrict " << valueName;
                }
                for (int i = 0; i < (int) tableTypes.size(); i++)
                    args << ", __global const " << tableTypes[i]<< "* restrict table" << i;
                replacements["PARAMETER_ARGUMENTS"] = args.str();
                if (loadPosAsDelta[step])
                    defines["LOAD_POS_AS_DELTA"] = "1";
                else if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
                    defines.erase("LOAD_POS_AS_DELTA");
                cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customIntegratorPerDof, replacements), defines);
                cl::Kernel kernel = cl::Kernel(program, "computePerDof");
                kernels[step].push_back(kernel);
                requiredGaussian[step] = numGaussian;
                requiredUniform[step] = numUniform;
                int index = 0;
                kernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
                setPosqCorrectionArg(cl, kernel, index++);
                kernel.setArg<cl::Buffer>(index++, integration.getPosDelta().getDeviceBuffer());
                kernel.setArg<cl::Buffer>(index++, cl.getVelm().getDeviceBuffer());
                kernel.setArg<cl::Buffer>(index++, cl.getForce().getDeviceBuffer());
                kernel.setArg<cl::Buffer>(index++, integration.getStepSize().getDeviceBuffer());
                kernel.setArg<cl::Buffer>(index++, globalValues.getDeviceBuffer());
                kernel.setArg<cl::Buffer>(index++, sumBuffer.getDeviceBuffer());
                index += 4;
                kernel.setArg<cl::Buffer>(index++, perDofEnergyParamDerivs.getDeviceBuffer());
                for (auto& array : perDofValues)
                    kernel.setArg<cl::Memory>(index++, array.getDeviceBuffer());
                for (auto& array : tabulatedFunctions)
                    kernel.setArg<cl::Buffer>(index++, array.getDeviceBuffer());
                if (stepType[step] == CustomIntegrator::ComputeSum) {
                    // Create a second kernel for this step that sums the values.

                    program = cl.createProgram(OpenCLKernelSources::customIntegrator, defines);
                    kernel = cl::Kernel(program, useDouble ? "computeDoubleSum" : "computeFloatSum");
                    kernels[step].push_back(kernel);
                    index = 0;
                    kernel.setArg<cl::Buffer>(index++, sumBuffer.getDeviceBuffer());
                    kernel.setArg<cl::Buffer>(index++, summedValue.getDeviceBuffer());
                    kernel.setArg<cl_int>(index++, numAtoms);
                }
            }
            else if (stepType[step] == CustomIntegrator::ConstrainPositions) {
                // Apply position constraints.

                cl::Program program = cl.createProgram(OpenCLKernelSources::customIntegrator, defines);
                cl::Kernel kernel = cl::Kernel(program, "applyPositionDeltas");
                kernels[step].push_back(kernel);
                int index = 0;
                kernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
                setPosqCorrectionArg(cl, kernel, index++);
                kernel.setArg<cl::Buffer>(index++, integration.getPosDelta().getDeviceBuffer());
            }
        }
        
        // Initialize the random number generator.
        
        int maxUniformRandoms = 1;
        for (int required : requiredUniform)
            maxUniformRandoms = max(maxUniformRandoms, required);
        uniformRandoms.initialize<mm_float4>(cl, maxUniformRandoms, "uniformRandoms");
        randomSeed.initialize<mm_int4>(cl, cl.getNumThreadBlocks()*OpenCLContext::ThreadBlockSize, "randomSeed");
        vector<mm_int4> seed(randomSeed.getSize());
        int rseed = integrator.getRandomNumberSeed();
        // A random seed of 0 means use a unique one
        if (rseed == 0)
            rseed = osrngseed();
        unsigned int r = (unsigned int) (rseed+1);
        for (auto& s : seed) {
            s.x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            s.y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            s.z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            s.w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        }
        randomSeed.upload(seed);
        cl::Program randomProgram = cl.createProgram(OpenCLKernelSources::customIntegrator, defines);
        randomKernel = cl::Kernel(randomProgram, "generateRandomNumbers");
        randomKernel.setArg<cl_int>(0, maxUniformRandoms);
        randomKernel.setArg<cl::Buffer>(1, uniformRandoms.getDeviceBuffer());
        randomKernel.setArg<cl::Buffer>(2, randomSeed.getDeviceBuffer());
        
        // Create the kernel for computing kinetic energy.

        stringstream computeKE;
        for (int i = 0; i < perDofValues.size(); i++)
            computeKE << tempType<<" perDof"<<cl.intToString(i)<<" = convert_"<<tempType<<"(perDofValues"<<cl.intToString(i)<<"[index].xyz);\n";
        Lepton::ParsedExpression keExpression = Lepton::Parser::parse(integrator.getKineticEnergyExpression()).optimize();
        computeKE << createPerDofComputation("", keExpression, integrator, "f", "", functionList, functionNames);
        map<string, string> replacements;
        replacements["COMPUTE_STEP"] = computeKE.str();
        stringstream args;
        for (int i = 0; i < perDofValues.size(); i++) {
            string valueName = "perDofValues"+cl.intToString(i);
            args << ", __global " << perDofType << "* restrict " << valueName;
        }
        for (int i = 0; i < (int) tableTypes.size(); i++)
            args << ", __global const " << tableTypes[i]<< "* restrict table" << i;
        replacements["PARAMETER_ARGUMENTS"] = args.str();
        if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
            defines.erase("LOAD_POS_AS_DELTA");
        cl::Program program = cl.createProgram(cl.replaceStrings(OpenCLKernelSources::customIntegratorPerDof, replacements), defines);
        kineticEnergyKernel = cl::Kernel(program, "computePerDof");
        int index = 0;
        kineticEnergyKernel.setArg<cl::Buffer>(index++, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kineticEnergyKernel, index++);
        kineticEnergyKernel.setArg<cl::Buffer>(index++, integration.getPosDelta().getDeviceBuffer());
        kineticEnergyKernel.setArg<cl::Buffer>(index++, cl.getVelm().getDeviceBuffer());
        kineticEnergyKernel.setArg<cl::Buffer>(index++, cl.getForce().getDeviceBuffer());
        kineticEnergyKernel.setArg<cl::Buffer>(index++, integration.getStepSize().getDeviceBuffer());
        kineticEnergyKernel.setArg<cl::Buffer>(index++, globalValues.getDeviceBuffer());
        kineticEnergyKernel.setArg<cl::Buffer>(index++, sumBuffer.getDeviceBuffer());
        index += 2;
        kineticEnergyKernel.setArg<cl::Buffer>(index++, uniformRandoms.getDeviceBuffer());
        if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision())
            kineticEnergyKernel.setArg<cl_double>(index++, 0.0);
        else
            kineticEnergyKernel.setArg<cl_float>(index++, 0.0f);
        kineticEnergyKernel.setArg<cl::Buffer>(index++, perDofEnergyParamDerivs.getDeviceBuffer());
        for (auto& array : perDofValues)
            kineticEnergyKernel.setArg<cl::Buffer>(index++, array.getDeviceBuffer());
        for (auto& array : tabulatedFunctions)
            kineticEnergyKernel.setArg<cl::Buffer>(index++, array.getDeviceBuffer());
        keNeedsForce = usesVariable(keExpression, "f");

        // Create a second kernel to sum the values.

        program = cl.createProgram(OpenCLKernelSources::customIntegrator, defines);
        sumKineticEnergyKernel = cl::Kernel(program, useDouble ? "computeDoubleSum" : "computeFloatSum");
        index = 0;
        sumKineticEnergyKernel.setArg<cl::Buffer>(index++, sumBuffer.getDeviceBuffer());
        sumKineticEnergyKernel.setArg<cl::Buffer>(index++, summedValue.getDeviceBuffer());
        sumKineticEnergyKernel.setArg<cl_int>(index++, numAtoms);

        // Delete the custom functions.

        for (auto& function : functions)
            delete function.second;
    }

    // Make sure all values (variables, parameters, etc.) are up to date.
    
    for (int i = 0; i < perDofValues.size(); i++) {
        if (!deviceValuesAreCurrent[i]) {
            if (useDouble)
                perDofValues[i].upload(localPerDofValuesDouble[i]);
            else
                perDofValues[i].upload(localPerDofValuesFloat[i]);
            deviceValuesAreCurrent[i] = true;
        }
        localValuesAreCurrent[i] = false;
    }
    double stepSize = integrator.getStepSize();
    recordGlobalValue(stepSize, GlobalTarget(DT, dtVariableIndex), integrator);
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        double value = context.getParameter(parameterNames[i]);
        if (value != localGlobalValues[parameterVariableIndex[i]]) {
            localGlobalValues[parameterVariableIndex[i]] = value;
            deviceGlobalsAreCurrent = false;
        }
    }
}

ExpressionTreeNode OpenCLIntegrateCustomStepKernel::replaceDerivFunctions(const ExpressionTreeNode& node, ContextImpl& context) {
    // This is called recursively to identify calls to the deriv() function inside global expressions,
    // and replace them with a custom function that returns the correct value.
    
    const Operation& op = node.getOperation();
    if (op.getId() == Operation::CUSTOM && op.getName() == "deriv") {
        string param = node.getChildren()[1].getOperation().getName();
        if (context.getParameters().find(param) == context.getParameters().end())
            throw OpenMMException("The second argument to deriv() must be a context parameter");
        needsEnergyParamDerivs = true;
        return ExpressionTreeNode(new Operation::Custom("deriv", new DerivFunction(energyParamDerivs, param)));
    }
    else {
        vector<ExpressionTreeNode> children;
        for (auto& child : node.getChildren())
            children.push_back(replaceDerivFunctions(child, context));
        return ExpressionTreeNode(op.clone(), children);
    }
}

void OpenCLIntegrateCustomStepKernel::findExpressionsForDerivs(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& variableNodes) {
    // This is called recursively to identify calls to the deriv() function inside per-DOF expressions,
    // and record the code to replace them with.
    
    const Operation& op = node.getOperation();
    if (op.getId() == Operation::CUSTOM && op.getName() == "deriv") {
        string param = node.getChildren()[1].getOperation().getName();
        int index;
        for (index = 0; index < perDofEnergyParamDerivNames.size() && param != perDofEnergyParamDerivNames[index]; index++)
            ;
        if (index == perDofEnergyParamDerivNames.size())
            perDofEnergyParamDerivNames.push_back(param);
        variableNodes.push_back(make_pair(node, "energyParamDerivs["+cl.intToString(index)+"]"));
        needsEnergyParamDerivs = true;
    }
    else {
        for (auto& child : node.getChildren())
            findExpressionsForDerivs(child, variableNodes);
    }
}

void OpenCLIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    prepareForComputation(context, integrator, forcesAreValid);
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    if (!forcesAreValid)
        savedEnergy.clear();
    
    // Loop over computation steps in the integrator and execute them.

    for (int step = 0; step < numSteps; ) {
        int nextStep = step+1;
        int forceGroups = forceGroupFlags[step];
        int lastForceGroups = context.getLastForceGroups();
        bool haveForces = (!needsForces[step] || (forcesAreValid && lastForceGroups == forceGroups));
        bool haveEnergy = (!needsEnergy[step] || savedEnergy.find(forceGroups) != savedEnergy.end());
        if (!haveForces || !haveEnergy) {
            if (forcesAreValid) {
                if (savedForces.find(lastForceGroups) != savedForces.end() && validSavedForces.find(lastForceGroups) == validSavedForces.end()) {
                    // The forces are still valid.  We just need a different force group right now.  Save the old
                    // forces in case we need them again.

                    cl.getForce().copyTo(savedForces[lastForceGroups]);
                    validSavedForces.insert(lastForceGroups);
                }
            }
            else
                validSavedForces.clear();
            
            // Recompute forces and/or energy.  Figure out what is actually needed
            // between now and the next time they get invalidated again.
            
            bool computeForce = (needsForces[step] || computeBothForceAndEnergy[step]);
            bool computeEnergy = (needsEnergy[step] || computeBothForceAndEnergy[step]);
            if (!computeEnergy && validSavedForces.find(forceGroups) != validSavedForces.end()) {
                // We can just restore the forces we saved earlier.
                
                savedForces[forceGroups].copyTo(cl.getForce());
                context.getLastForceGroups() = forceGroups;
            }
            else {
                recordChangedParameters(context);
                energy = context.calcForcesAndEnergy(computeForce, computeEnergy, forceGroups);
                savedEnergy[forceGroups] = energy;
                if (needsEnergyParamDerivs) {
                    context.getEnergyParameterDerivatives(energyParamDerivs);
                    if (perDofEnergyParamDerivNames.size() > 0) {
                        for (int i = 0; i < perDofEnergyParamDerivNames.size(); i++)
                            localPerDofEnergyParamDerivs[i] = energyParamDerivs[perDofEnergyParamDerivNames[i]];
                        perDofEnergyParamDerivs.upload(localPerDofEnergyParamDerivs, true, true);
                    }
                }
                forcesAreValid = true;
            }
        }
        if (needsEnergy[step])
            energy = savedEnergy[forceGroups];
        if (needsGlobals[step] && !deviceGlobalsAreCurrent) {
            // Upload the global values to the device.
            
            globalValues.upload(localGlobalValues, true, true);
            deviceGlobalsAreCurrent = true;
        }
        bool stepInvalidatesForces = invalidatesForces[step];
        if (stepType[step] == CustomIntegrator::ComputePerDof && !merged[step]) {
            kernels[step][0].setArg<cl_uint>(9, integration.prepareRandomNumbers(requiredGaussian[step]));
            kernels[step][0].setArg<cl::Buffer>(8, integration.getRandom().getDeviceBuffer());
            kernels[step][0].setArg<cl::Buffer>(10, uniformRandoms.getDeviceBuffer());
            if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision())
                kernels[step][0].setArg<cl_double>(11, energy);
            else
                kernels[step][0].setArg<cl_float>(11, (cl_float) energy);
            if (requiredUniform[step] > 0)
                cl.executeKernel(randomKernel, numAtoms);
            cl.executeKernel(kernels[step][0], numAtoms, 128);
        }
        else if (stepType[step] == CustomIntegrator::ComputeGlobal) {
            expressionSet.setVariable(uniformVariableIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
            expressionSet.setVariable(gaussianVariableIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
            expressionSet.setVariable(stepEnergyVariableIndex[step], energy);
            recordGlobalValue(globalExpressions[step][0].evaluate(), stepTarget[step], integrator);
        }
        else if (stepType[step] == CustomIntegrator::ComputeSum) {
            kernels[step][0].setArg<cl_uint>(9, integration.prepareRandomNumbers(requiredGaussian[step]));
            kernels[step][0].setArg<cl::Buffer>(8, integration.getRandom().getDeviceBuffer());
            kernels[step][0].setArg<cl::Buffer>(10, uniformRandoms.getDeviceBuffer());
            if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision())
                kernels[step][0].setArg<cl_double>(11, energy);
            else
                kernels[step][0].setArg<cl_float>(11, (cl_float) energy);
            if (requiredUniform[step] > 0)
                cl.executeKernel(randomKernel, numAtoms);
            cl.clearBuffer(sumBuffer);
            cl.executeKernel(kernels[step][0], numAtoms, 128);
            cl.executeKernel(kernels[step][1], sumWorkGroupSize, sumWorkGroupSize);
            if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
                double value;
                summedValue.download(&value);
                recordGlobalValue(value, stepTarget[step], integrator);
            }
            else {
                float value;
                summedValue.download(&value);
                recordGlobalValue(value, stepTarget[step], integrator);
            }
        }
        else if (stepType[step] == CustomIntegrator::UpdateContextState) {
            recordChangedParameters(context);
            stepInvalidatesForces = context.updateContextState();
        }
        else if (stepType[step] == CustomIntegrator::ConstrainPositions) {
            if (hasAnyConstraints) {
                cl.getIntegrationUtilities().applyConstraints(integrator.getConstraintTolerance());
                cl.executeKernel(kernels[step][0], numAtoms);
            }
            cl.getIntegrationUtilities().computeVirtualSites();
        }
        else if (stepType[step] == CustomIntegrator::ConstrainVelocities) {
            cl.getIntegrationUtilities().applyVelocityConstraints(integrator.getConstraintTolerance());
        }
        else if (stepType[step] == CustomIntegrator::IfBlockStart) {
            if (!evaluateCondition(step))
                nextStep = blockEnd[step]+1;
        }
        else if (stepType[step] == CustomIntegrator::WhileBlockStart) {
            if (!evaluateCondition(step))
                nextStep = blockEnd[step]+1;
        }
        else if (stepType[step] == CustomIntegrator::BlockEnd) {
            if (blockEnd[step] != -1)
                nextStep = blockEnd[step]; // Return to the start of a while block.
        }
        if (stepInvalidatesForces) {
            forcesAreValid = false;
            savedEnergy.clear();
        }
        step = nextStep;
    }
    recordChangedParameters(context);

    // Update the time and step count.

    cl.setTime(cl.getTime()+integrator.getStepSize());
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    if (cl.getAtomsWereReordered()) {
        forcesAreValid = false;
        validSavedForces.clear();
    }
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

bool OpenCLIntegrateCustomStepKernel::evaluateCondition(int step) {
    expressionSet.setVariable(uniformVariableIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
    expressionSet.setVariable(gaussianVariableIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
    expressionSet.setVariable(stepEnergyVariableIndex[step], energy);
    double lhs = globalExpressions[step][0].evaluate();
    double rhs = globalExpressions[step][1].evaluate();
    switch (comparisons[step]) {
        case CustomIntegratorUtilities::EQUAL:
            return (lhs == rhs);
        case CustomIntegratorUtilities::LESS_THAN:
            return (lhs < rhs);
        case CustomIntegratorUtilities::GREATER_THAN:
            return (lhs > rhs);
        case CustomIntegratorUtilities::NOT_EQUAL:
            return (lhs != rhs);
        case CustomIntegratorUtilities::LESS_THAN_OR_EQUAL:
            return (lhs <= rhs);
        case CustomIntegratorUtilities::GREATER_THAN_OR_EQUAL:
            return (lhs >= rhs);
    }
    throw OpenMMException("Invalid comparison operator");
}

double OpenCLIntegrateCustomStepKernel::computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    prepareForComputation(context, integrator, forcesAreValid);
    if (keNeedsForce && !forcesAreValid) {
        // Compute the force.  We want to then mark that forces are valid, which means also computing
        // potential energy if any steps will expect it to be valid too.
        
        bool willNeedEnergy = false;
        for (int i = 0; i < integrator.getNumComputations(); i++)
            willNeedEnergy |= needsEnergy[i];
        energy = context.calcForcesAndEnergy(true, willNeedEnergy, -1);
        forcesAreValid = true;
    }
    cl.clearBuffer(sumBuffer);
    kineticEnergyKernel.setArg<cl::Buffer>(8, cl.getIntegrationUtilities().getRandom().getDeviceBuffer());
    kineticEnergyKernel.setArg<cl_uint>(9, 0);
    cl.executeKernel(kineticEnergyKernel, cl.getNumAtoms());
    cl.executeKernel(sumKineticEnergyKernel, sumWorkGroupSize, sumWorkGroupSize);
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        double ke;
        summedValue.download(&ke);
        return ke;
    }
    else {
        float ke;
        summedValue.download(&ke);
        return ke;
    }
}

void OpenCLIntegrateCustomStepKernel::recordGlobalValue(double value, GlobalTarget target, CustomIntegrator& integrator) {
    switch (target.type) {
        case DT:
            if (value != localGlobalValues[dtVariableIndex])
                deviceGlobalsAreCurrent = false;
            expressionSet.setVariable(dtVariableIndex, value);
            localGlobalValues[dtVariableIndex] = value;
            cl.getIntegrationUtilities().setNextStepSize(value);
            integrator.setStepSize(value);
            break;
        case VARIABLE:
        case PARAMETER:
            expressionSet.setVariable(target.variableIndex, value);
            localGlobalValues[target.variableIndex] = value;
            deviceGlobalsAreCurrent = false;
            break;
    }
}

void OpenCLIntegrateCustomStepKernel::recordChangedParameters(ContextImpl& context) {
    if (!modifiesParameters)
        return;
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        double value = context.getParameter(parameterNames[i]);
        if (value != localGlobalValues[parameterVariableIndex[i]])
            context.setParameter(parameterNames[i], localGlobalValues[parameterVariableIndex[i]]);
    }
}

void OpenCLIntegrateCustomStepKernel::getGlobalVariables(ContextImpl& context, vector<double>& values) const {
    if (!globalValues.isInitialized()) {
        // The data structures haven't been created yet, so just return the list of values that was given earlier.
        
        values = initialGlobalVariables;
        return;
    }
    values.resize(numGlobalVariables);
    for (int i = 0; i < numGlobalVariables; i++)
        values[i] = localGlobalValues[globalVariableIndex[i]];
}

void OpenCLIntegrateCustomStepKernel::setGlobalVariables(ContextImpl& context, const vector<double>& values) {
    if (numGlobalVariables == 0)
        return;
    if (!globalValues.isInitialized()) {
        // The data structures haven't been created yet, so just store the list of values.
        
        initialGlobalVariables = values;
        return;
    }
    for (int i = 0; i < numGlobalVariables; i++) {
        localGlobalValues[globalVariableIndex[i]] = values[i];
        expressionSet.setVariable(globalVariableIndex[i], values[i]);
    }
    deviceGlobalsAreCurrent = false;
}

void OpenCLIntegrateCustomStepKernel::getPerDofVariable(ContextImpl& context, int variable, vector<Vec3>& values) const {
    values.resize(perDofValues[variable].getSize());
    const vector<int>& order = cl.getAtomIndex();
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        if (!localValuesAreCurrent[variable]) {
            perDofValues[variable].download(localPerDofValuesDouble[variable]);
            localValuesAreCurrent[variable] = true;
        }
        for (int i = 0; i < (int) values.size(); i++) {
            values[order[i]][0] = localPerDofValuesDouble[variable][i].x;
            values[order[i]][1] = localPerDofValuesDouble[variable][i].y;
            values[order[i]][2] = localPerDofValuesDouble[variable][i].z;
        }
    }
    else {
        if (!localValuesAreCurrent[variable]) {
            perDofValues[variable].download(localPerDofValuesFloat[variable]);
            localValuesAreCurrent[variable] = true;
        }
        for (int i = 0; i < (int) values.size(); i++) {
            values[order[i]][0] = localPerDofValuesFloat[variable][i].x;
            values[order[i]][1] = localPerDofValuesFloat[variable][i].y;
            values[order[i]][2] = localPerDofValuesFloat[variable][i].z;
        }
    }
}

void OpenCLIntegrateCustomStepKernel::setPerDofVariable(ContextImpl& context, int variable, const vector<Vec3>& values) {
    const vector<int>& order = cl.getAtomIndex();
    localValuesAreCurrent[variable] = true;
    deviceValuesAreCurrent[variable] = false;
    if (cl.getUseDoublePrecision() || cl.getUseMixedPrecision()) {
        localPerDofValuesDouble[variable].resize(values.size());
        for (int i = 0; i < (int) values.size(); i++)
            localPerDofValuesDouble[variable][i] = mm_double4(values[order[i]][0], values[order[i]][1], values[order[i]][2], 0);
    }
    else {
        localPerDofValuesFloat[variable].resize(values.size());
        for (int i = 0; i < (int) values.size(); i++)
            localPerDofValuesFloat[variable][i] = mm_float4(values[order[i]][0], values[order[i]][1], values[order[i]][2], 0);
    }
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

void OpenCLRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    frequency = force.getFrequency();
    int numAtoms = cl.getNumAtoms();
    cmMomentum.initialize<mm_float4>(cl, (numAtoms+OpenCLContext::ThreadBlockSize-1)/OpenCLContext::ThreadBlockSize, "cmMomentum");
    double totalMass = 0.0;
    for (int i = 0; i < numAtoms; i++)
        totalMass += system.getParticleMass(i);
    map<string, string> defines;
    defines["INVERSE_TOTAL_MASS"] = cl.doubleToString(totalMass == 0 ? 0.0 : 1.0/totalMass);
    cl::Program program = cl.createProgram(OpenCLKernelSources::removeCM, defines);
    kernel1 = cl::Kernel(program, "calcCenterOfMassMomentum");
    kernel1.setArg<cl_int>(0, numAtoms);
    kernel1.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    kernel1.setArg<cl::Buffer>(2, cmMomentum.getDeviceBuffer());
    kernel1.setArg(3, OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
    kernel2 = cl::Kernel(program, "removeCenterOfMassMomentum");
    kernel2.setArg<cl_int>(0, numAtoms);
    kernel2.setArg<cl::Buffer>(1, cl.getVelm().getDeviceBuffer());
    kernel2.setArg<cl::Buffer>(2, cmMomentum.getDeviceBuffer());
    kernel2.setArg(3, OpenCLContext::ThreadBlockSize*sizeof(mm_float4), NULL);
}

void OpenCLRemoveCMMotionKernel::execute(ContextImpl& context) {
    cl.executeKernel(kernel1, cl.getNumAtoms());
    cl.executeKernel(kernel2, cl.getNumAtoms());
}
