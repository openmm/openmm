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

#include "CudaKernels.h"
#include "CudaForceInfo.h"
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/internal/CustomHbondForceImpl.h"
#include "openmm/internal/CustomManyParticleForceImpl.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "CudaBondedUtilities.h"
#include "CudaExpressionUtilities.h"
#include "CudaIntegrationUtilities.h"
#include "CudaNonbondedUtilities.h"
#include "CudaKernelSources.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "SimTKOpenMMRealType.h"
#include "SimTKOpenMMUtilities.h"
#include <algorithm>
#include <cmath>
#include <set>

using namespace OpenMM;
using namespace std;
using Lepton::ExpressionTreeNode;
using Lepton::Operation;

#define CHECK_RESULT(result, prefix) \
    if (result != CUDA_SUCCESS) { \
        std::stringstream m; \
        m<<prefix<<": "<<CudaContext::getErrorString(result)<<" ("<<result<<")"<<" at "<<__FILE__<<":"<<__LINE__; \
        throw OpenMMException(m.str());\
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
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        if (usesVariable(node.getChildren()[i], variable))
            return true;
    return false;
}

static bool usesVariable(const Lepton::ParsedExpression& expression, const string& variable) {
    return usesVariable(expression.getRootNode(), variable);
}

static pair<ExpressionTreeNode, string> makeVariable(const string& name, const string& value) {
    return make_pair(ExpressionTreeNode(new Operation::Variable(name)), value);
}

void CudaCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void CudaCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cu.setAsCurrent();
    cu.clearAutoclearBuffers();
    for (vector<CudaContext::ForcePreComputation*>::iterator iter = cu.getPreComputations().begin(); iter != cu.getPreComputations().end(); ++iter)
        (*iter)->computeForceAndEnergy(includeForces, includeEnergy, groups);
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    bool includeNonbonded = ((groups&(1<<nb.getForceGroup())) != 0);
    cu.setComputeForceCount(cu.getComputeForceCount()+1);
    if (includeNonbonded)
        nb.prepareInteractions();
}

double CudaCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cu.getBondedUtilities().computeInteractions(groups);
    if ((groups&(1<<cu.getNonbondedUtilities().getForceGroup())) != 0)
        cu.getNonbondedUtilities().computeInteractions();
    double sum = 0.0;
    for (vector<CudaContext::ForcePostComputation*>::iterator iter = cu.getPostComputations().begin(); iter != cu.getPostComputations().end(); ++iter)
        sum += (*iter)->computeForceAndEnergy(includeForces, includeEnergy, groups);
    cu.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy) {
        CudaArray& energyArray = cu.getEnergyBuffer();
        if (cu.getUseDoublePrecision()) {
            double* energy = (double*) cu.getPinnedBuffer();
            energyArray.download(energy);
            for (int i = 0; i < energyArray.getSize(); i++)
                sum += energy[i];
        }
        else {
            float* energy = (float*) cu.getPinnedBuffer();
            energyArray.download(energy);
            for (int i = 0; i < energyArray.getSize(); i++)
                sum += energy[i];
        }
    }
    return sum;
}

void CudaUpdateStateDataKernel::initialize(const System& system) {
}

double CudaUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return cu.getTime();
}

void CudaUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    vector<CudaContext*>& contexts = cu.getPlatformData().contexts;
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->setTime(time);
}

void CudaUpdateStateDataKernel::getPositions(ContextImpl& context, vector<Vec3>& positions) {
    cu.setAsCurrent();
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    double4 periodicBoxSize = cu.getPeriodicBoxSize();
    if (cu.getUseDoublePrecision()) {
        double4* posq = (double4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            double4 pos = posq[i];
            int4 offset = cu.getPosCellOffsets()[i];
            positions[order[i]] = Vec3(pos.x-offset.x*periodicBoxSize.x, pos.y-offset.y*periodicBoxSize.y, pos.z-offset.z*periodicBoxSize.z);
        }
    }
    else if (cu.getUseMixedPrecision()) {
        float4* posq = (float4*) cu.getPinnedBuffer();
        vector<float4> posCorrection;
        cu.getPosq().download(posq);
        cu.getPosqCorrection().download(posCorrection);
        for (int i = 0; i < numParticles; ++i) {
            float4 pos1 = posq[i];
            float4 pos2 = posCorrection[i];
            int4 offset = cu.getPosCellOffsets()[i];
            positions[order[i]] = Vec3((double)pos1.x+(double)pos2.x-offset.x*periodicBoxSize.x, (double)pos1.y+(double)pos2.y-offset.y*periodicBoxSize.y, (double)pos1.z+(double)pos2.z-offset.z*periodicBoxSize.z);
        }
    }
    else {
        float4* posq = (float4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            float4 pos = posq[i];
            int4 offset = cu.getPosCellOffsets()[i];
            positions[order[i]] = Vec3(pos.x-offset.x*periodicBoxSize.x, pos.y-offset.y*periodicBoxSize.y, pos.z-offset.z*periodicBoxSize.z);
        }
    }
}

void CudaUpdateStateDataKernel::setPositions(ContextImpl& context, const vector<Vec3>& positions) {
    cu.setAsCurrent();
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cu.getUseDoublePrecision()) {
        double4* posq = (double4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            double4& pos = posq[i];
            const Vec3& p = positions[order[i]];
            pos.x = p[0];
            pos.y = p[1];
            pos.z = p[2];
        }
        for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
            posq[i] = make_double4(0.0, 0.0, 0.0, 0.0);
        cu.getPosq().upload(posq);
    }
    else {
        float4* posq = (float4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            float4& pos = posq[i];
            const Vec3& p = positions[order[i]];
            pos.x = (float) p[0];
            pos.y = (float) p[1];
            pos.z = (float) p[2];
        }
        for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
            posq[i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cu.getPosq().upload(posq);
    }
    if (cu.getUseMixedPrecision()) {
        float4* posCorrection = (float4*) cu.getPinnedBuffer();
        for (int i = 0; i < numParticles; ++i) {
            float4& c = posCorrection[i];
            const Vec3& p = positions[order[i]];
            c.x = (float) (p[0]-(float)p[0]);
            c.y = (float) (p[1]-(float)p[1]);
            c.z = (float) (p[2]-(float)p[2]);
            c.w = 0;
        }
        for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
            posCorrection[i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cu.getPosqCorrection().upload(posCorrection);
    }
    for (int i = 0; i < (int) cu.getPosCellOffsets().size(); i++)
        cu.getPosCellOffsets()[i] = make_int4(0, 0, 0, 0);
    cu.reorderAtoms();
}

void CudaUpdateStateDataKernel::getVelocities(ContextImpl& context, vector<Vec3>& velocities) {
    cu.setAsCurrent();
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        double4* velm = (double4*) cu.getPinnedBuffer();
        cu.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            double4 vel = velm[i];
            int4 offset = cu.getPosCellOffsets()[i];
            velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
        }
    }
    else {
        float4* velm = (float4*) cu.getPinnedBuffer();
        cu.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            float4 vel = velm[i];
            int4 offset = cu.getPosCellOffsets()[i];
            velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
        }
    }
}

void CudaUpdateStateDataKernel::setVelocities(ContextImpl& context, const vector<Vec3>& velocities) {
    cu.setAsCurrent();
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        double4* velm = (double4*) cu.getPinnedBuffer();
        cu.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            double4& vel = velm[i];
            const Vec3& p = velocities[order[i]];
            vel.x = p[0];
            vel.y = p[1];
            vel.z = p[2];
        }
        for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
            velm[i] = make_double4(0.0, 0.0, 0.0, 0.0);
        cu.getVelm().upload(velm);
    }
    else {
        float4* velm = (float4*) cu.getPinnedBuffer();
        cu.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            float4& vel = velm[i];
            const Vec3& p = velocities[order[i]];
            vel.x = p[0];
            vel.y = p[1];
            vel.z = p[2];
        }
        for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
            velm[i] = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cu.getVelm().upload(velm);
    }
}

void CudaUpdateStateDataKernel::getForces(ContextImpl& context, vector<Vec3>& forces) {
    cu.setAsCurrent();
    long long* force = (long long*) cu.getPinnedBuffer();
    cu.getForce().download(force);
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    int paddedNumParticles = cu.getPaddedNumAtoms();
    forces.resize(numParticles);
    double scale = 1.0/(double) 0x100000000LL;
    for (int i = 0; i < numParticles; ++i)
        forces[order[i]] = Vec3(scale*force[i], scale*force[i+paddedNumParticles], scale*force[i+paddedNumParticles*2]);
}

void CudaUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    double4 box = cu.getPeriodicBoxSize();
    a = Vec3(box.x, 0, 0);
    b = Vec3(0, box.y, 0);
    c = Vec3(0, 0, box.z);
}

void CudaUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) const {
    vector<CudaContext*>& contexts = cu.getPlatformData().contexts;
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->setPeriodicBoxSize(a[0], b[1], c[2]);
}

void CudaUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    cu.setAsCurrent();
    int version = 1;
    stream.write((char*) &version, sizeof(int));
    int precision = (cu.getUseDoublePrecision() ? 2 : cu.getUseMixedPrecision() ? 1 : 0);
    stream.write((char*) &precision, sizeof(int));
    double time = cu.getTime();
    stream.write((char*) &time, sizeof(double));
    int stepCount = cu.getStepCount();
    stream.write((char*) &stepCount, sizeof(int));
    int stepsSinceReorder = cu.getStepsSinceReorder();
    stream.write((char*) &stepsSinceReorder, sizeof(int));
    char* buffer = (char*) cu.getPinnedBuffer();
    cu.getPosq().download(buffer);
    stream.write(buffer, cu.getPosq().getSize()*cu.getPosq().getElementSize());
    if (cu.getUseMixedPrecision()) {
        cu.getPosqCorrection().download(buffer);
        stream.write(buffer, cu.getPosqCorrection().getSize()*cu.getPosqCorrection().getElementSize());
    }
    cu.getVelm().download(buffer);
    stream.write(buffer, cu.getVelm().getSize()*cu.getVelm().getElementSize());
    stream.write((char*) &cu.getAtomIndex()[0], sizeof(int)*cu.getAtomIndex().size());
    stream.write((char*) &cu.getPosCellOffsets()[0], sizeof(int4)*cu.getPosCellOffsets().size());
    double4 box = cu.getPeriodicBoxSize();
    stream.write((char*) &box, sizeof(double4));
    cu.getIntegrationUtilities().createCheckpoint(stream);
    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void CudaUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    cu.setAsCurrent();
    int version;
    stream.read((char*) &version, sizeof(int));
    if (version != 1)
        throw OpenMMException("Checkpoint was created with a different version of OpenMM");
    int precision;
    stream.read((char*) &precision, sizeof(int));
    int expectedPrecision = (cu.getUseDoublePrecision() ? 2 : cu.getUseMixedPrecision() ? 1 : 0);
    if (precision != expectedPrecision)
        throw OpenMMException("Checkpoint was created with a different numeric precision");
    double time;
    stream.read((char*) &time, sizeof(double));
    int stepCount, stepsSinceReorder;
    stream.read((char*) &stepCount, sizeof(int));
    stream.read((char*) &stepsSinceReorder, sizeof(int));
    vector<CudaContext*>& contexts = cu.getPlatformData().contexts;
    for (int i = 0; i < (int) contexts.size(); i++) {
        contexts[i]->setTime(time);
        contexts[i]->setStepCount(stepCount);
        contexts[i]->setStepsSinceReorder(stepsSinceReorder);
    }
    char* buffer = (char*) cu.getPinnedBuffer();
    stream.read(buffer, cu.getPosq().getSize()*cu.getPosq().getElementSize());
    cu.getPosq().upload(buffer);
    if (cu.getUseMixedPrecision()) {
        stream.read(buffer, cu.getPosqCorrection().getSize()*cu.getPosqCorrection().getElementSize());
        cu.getPosqCorrection().upload(buffer);
    }
    stream.read(buffer, cu.getVelm().getSize()*cu.getVelm().getElementSize());
    cu.getVelm().upload(buffer);
    stream.read((char*) &cu.getAtomIndex()[0], sizeof(int)*cu.getAtomIndex().size());
    cu.getAtomIndexArray().upload(cu.getAtomIndex());
    stream.read((char*) &cu.getPosCellOffsets()[0], sizeof(int4)*cu.getPosCellOffsets().size());
    double4 box;
    stream.read((char*) &box, sizeof(double4));
    for (int i = 0; i < (int) contexts.size(); i++)
        contexts[i]->setPeriodicBoxSize(box.x, box.y, box.z);
    cu.getIntegrationUtilities().loadCheckpoint(stream);
    SimTKOpenMMUtilities::loadCheckpoint(stream);
    for (int i = 0; i < cu.getReorderListeners().size(); i++)
        cu.getReorderListeners()[i]->execute();
}

void CudaApplyConstraintsKernel::initialize(const System& system) {
}

void CudaApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
    cu.setAsCurrent();
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        map<string, string> defines;
        CUmodule module = cu.createModule(CudaKernelSources::constraints, defines);
        applyDeltasKernel = cu.getKernel(module, "applyPositionDeltas");
    }
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    cu.clearBuffer(integration.getPosDelta());
    integration.applyConstraints(tol);
    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    int numAtoms = cu.getNumAtoms();
    void* args[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &posCorrection, &cu.getIntegrationUtilities().getPosDelta().getDevicePointer()};
    cu.executeKernel(applyDeltasKernel, args, cu.getNumAtoms());
    integration.computeVirtualSites();
}

void CudaApplyConstraintsKernel::applyToVelocities(ContextImpl& context, double tol) {
    cu.getIntegrationUtilities().applyVelocityConstraints(tol);
}

void CudaVirtualSitesKernel::initialize(const System& system) {
}

void CudaVirtualSitesKernel::computePositions(ContextImpl& context) {
    cu.getIntegrationUtilities().computeVirtualSites();
}

class CudaHarmonicBondForceInfo : public CudaForceInfo {
public:
    CudaHarmonicBondForceInfo(const HarmonicBondForce& force) : force(force) {
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

CudaCalcHarmonicBondForceKernel::~CudaCalcHarmonicBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = CudaArray::create<float2>(cu, numBonds, "bondParams");
    vector<float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = make_float2((float) length, (float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaKernelSources::harmonicBondForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
    cu.addForce(new CudaHarmonicBondForceInfo(force));
}

double CudaCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(startIndex+i, atom1, atom2, length, k);
        paramVector[i] = make_float2((float) length, (float) k);
    }
    params->upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCustomBondForceInfo : public CudaForceInfo {
public:
    CudaCustomBondForceInfo(const CustomBondForce& force) : force(force) {
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

CudaCalcCustomBondForceKernel::~CudaCalcCustomBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
}

void CudaCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = new CudaParameterSet(cu, force.getNumPerBondParameters(), numBonds, "customBondParams");
    vector<vector<float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    cu.addForce(new CudaCustomBondForceInfo(force));

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
    expressions["float dEdR = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["r"] = "r";
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customBondGlobals");
        globals->upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals->getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cu.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
}

double CudaCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    return 0.0;
}

void CudaCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
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
    
    cu.invalidateMolecules();
}

class CudaHarmonicAngleForceInfo : public CudaForceInfo {
public:
    CudaHarmonicAngleForceInfo(const HarmonicAngleForce& force) : force(force) {
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

CudaCalcHarmonicAngleForceKernel::~CudaCalcHarmonicAngleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = CudaArray::create<float2>(cu, numAngles, "angleParams");
    vector<float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], angle, k);
        paramVector[i] = make_float2((float) angle, (float) k);

    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaKernelSources::harmonicAngleForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float2");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::angleForce, replacements), force.getForceGroup());
    cu.addForce(new CudaHarmonicAngleForceInfo(force));
}

double CudaCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0)
        return;
    
    // Record the per-angle parameters.
    
    vector<float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        int atom1, atom2, atom3;
        double angle, k;
        force.getAngleParameters(startIndex+i, atom1, atom2, atom3, angle, k);
        paramVector[i] = make_float2((float) angle, (float) k);
    }
    params->upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCustomAngleForceInfo : public CudaForceInfo {
public:
    CudaCustomAngleForceInfo(const CustomAngleForce& force) : force(force) {
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

CudaCalcCustomAngleForceKernel::~CudaCalcCustomAngleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
}

void CudaCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = new CudaParameterSet(cu, force.getNumPerAngleParameters(), numAngles, "customAngleParams");
    vector<vector<float> > paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        vector<double> parameters;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    cu.addForce(new CudaCustomAngleForceInfo(force));

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
    expressions["float dEdAngle = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["theta"] = "theta";
    for (int i = 0; i < force.getNumPerAngleParameters(); i++) {
        const string& name = force.getPerAngleParameterName(i);
        variables[name] = "angleParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customAngleGlobals");
        globals->upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals->getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" angleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cu.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::angleForce, replacements), force.getForceGroup());
}

double CudaCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    return 0.0;
}

void CudaCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
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
    
    cu.invalidateMolecules();
}

class CudaPeriodicTorsionForceInfo : public CudaForceInfo {
public:
    CudaPeriodicTorsionForceInfo(const PeriodicTorsionForce& force) : force(force) {
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

CudaCalcPeriodicTorsionForceKernel::~CudaCalcPeriodicTorsionForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params = CudaArray::create<float4>(cu, numTorsions, "periodicTorsionParams");
    vector<float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], periodicity, phase, k);
        paramVector[i] = make_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaKernelSources::periodicTorsionForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params->getDevicePointer(), "float4");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
    cu.addForce(new CudaPeriodicTorsionForceInfo(force));
}

double CudaCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4, periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, periodicity, phase, k);
        paramVector[i] = make_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params->upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaRBTorsionForceInfo : public CudaForceInfo {
public:
    CudaRBTorsionForceInfo(const RBTorsionForce& force) : force(force) {
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

CudaCalcRBTorsionForceKernel::~CudaCalcRBTorsionForceKernel() {
    cu.setAsCurrent();
    if (params1 != NULL)
        delete params1;
    if (params2 != NULL)
        delete params2;
}

void CudaCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params1 = CudaArray::create<float4>(cu, numTorsions, "rbTorsionParams1");
    params2 = CudaArray::create<float2>(cu, numTorsions, "rbTorsionParams2");
    vector<float4> paramVector1(numTorsions);
    vector<float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], c0, c1, c2, c3, c4, c5);
        paramVector1[i] = make_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = make_float2((float) c4, (float) c5);

    }
    params1->upload(paramVector1);
    params2->upload(paramVector2);
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = CudaKernelSources::rbTorsionForce;
    replacements["PARAMS1"] = cu.getBondedUtilities().addArgument(params1->getDevicePointer(), "float4");
    replacements["PARAMS2"] = cu.getBondedUtilities().addArgument(params2->getDevicePointer(), "float2");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
    cu.addForce(new CudaRBTorsionForceInfo(force));
}

double CudaCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0)
        return;
    
    // Record the per-torsion parameters.
    
    vector<float4> paramVector1(numTorsions);
    vector<float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int atom1, atom2, atom3, atom4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        paramVector1[i] = make_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = make_float2((float) c4, (float) c5);
    }
    params1->upload(paramVector1);
    params2->upload(paramVector2);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCMAPTorsionForceInfo : public CudaForceInfo {
public:
    CudaCMAPTorsionForceInfo(const CMAPTorsionForce& force) : force(force) {
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

CudaCalcCMAPTorsionForceKernel::~CudaCalcCMAPTorsionForceKernel() {
    if (coefficients != NULL)
        delete coefficients;
    if (mapPositions != NULL)
        delete mapPositions;
    if (torsionMaps != NULL)
        delete torsionMaps;
}

void CudaCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    int numMaps = force.getNumMaps();
    vector<float4> coeffVec;
    vector<int2> mapPositionsVec(numMaps);
    vector<double> energy;
    vector<vector<double> > c;
    int currentPosition = 0;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        mapPositionsVec[i] = make_int2(currentPosition, size);
        currentPosition += 4*size*size;
        for (int j = 0; j < size*size; j++) {
            coeffVec.push_back(make_float4((float) c[j][0], (float) c[j][1], (float) c[j][2], (float) c[j][3]));
            coeffVec.push_back(make_float4((float) c[j][4], (float) c[j][5], (float) c[j][6], (float) c[j][7]));
            coeffVec.push_back(make_float4((float) c[j][8], (float) c[j][9], (float) c[j][10], (float) c[j][11]));
            coeffVec.push_back(make_float4((float) c[j][12], (float) c[j][13], (float) c[j][14], (float) c[j][15]));
        }
    }
    vector<vector<int> > atoms(numTorsions, vector<int>(8));
    vector<int> torsionMapsVec(numTorsions);
    for (int i = 0; i < numTorsions; i++)
        force.getTorsionParameters(startIndex+i, torsionMapsVec[i], atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], atoms[i][6], atoms[i][7]);
    coefficients = CudaArray::create<float4>(cu, coeffVec.size(), "cmapTorsionCoefficients");
    mapPositions = CudaArray::create<int2>(cu, numMaps, "cmapTorsionMapPositions");
    torsionMaps = CudaArray::create<int>(cu, numTorsions, "cmapTorsionMaps");
    coefficients->upload(coeffVec);
    mapPositions->upload(mapPositionsVec);
    torsionMaps->upload(torsionMapsVec);
    map<string, string> replacements;
    replacements["COEFF"] = cu.getBondedUtilities().addArgument(coefficients->getDevicePointer(), "float4");
    replacements["MAP_POS"] = cu.getBondedUtilities().addArgument(mapPositions->getDevicePointer(), "int2");
    replacements["MAPS"] = cu.getBondedUtilities().addArgument(torsionMaps->getDevicePointer(), "int");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::cmapTorsionForce, replacements), force.getForceGroup());
    cu.addForce(new CudaCMAPTorsionForceInfo(force));
}

double CudaCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

class CudaCustomTorsionForceInfo : public CudaForceInfo {
public:
    CudaCustomTorsionForceInfo(const CustomTorsionForce& force) : force(force) {
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

CudaCalcCustomTorsionForceKernel::~CudaCalcCustomTorsionForceKernel() {
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
}

void CudaCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params = new CudaParameterSet(cu, force.getNumPerTorsionParameters(), numTorsions, "customTorsionParams");
    vector<vector<float> > paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        vector<double> parameters;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    cu.addForce(new CudaCustomTorsionForceInfo(force));

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
    expressions["float dEdAngle = "] = forceExpression;

    // Create the kernels.

    map<string, string> variables;
    variables["theta"] = "theta";
    for (int i = 0; i < force.getNumPerTorsionParameters(); i++) {
        const string& name = force.getPerTorsionParameterName(i);
        variables[name] = "torsionParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customTorsionGlobals");
        globals->upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals->getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" torsionParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cu.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
}

double CudaCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    return 0.0;
}

void CudaCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
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
    
    cu.invalidateMolecules();
}

class CudaNonbondedForceInfo : public CudaForceInfo {
public:
    CudaNonbondedForceInfo(const NonbondedForce& force) : force(force) {
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

class CudaCalcNonbondedForceKernel::PmeIO : public CalcPmeReciprocalForceKernel::IO {
public:
    PmeIO(CudaContext& cu, CUfunction addForcesKernel) : cu(cu), addForcesKernel(addForcesKernel), forceTemp(NULL) {
        forceTemp = CudaArray::create<float4>(cu, cu.getNumAtoms(), "PmeForce");
    }
    ~PmeIO() {
        if (forceTemp != NULL)
            delete forceTemp;
    }
    float* getPosq() {
        cu.setAsCurrent();
        cu.getPosq().download(posq);
        return (float*) &posq[0];
    }
    void setForce(float* force) {
        forceTemp->upload(force);
        void* args[] = {&forceTemp->getDevicePointer(), &cu.getForce().getDevicePointer()};
        cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    }
private:
    CudaContext& cu;
    vector<float4> posq;
    CudaArray* forceTemp;
    CUfunction addForcesKernel;
};

class CudaCalcNonbondedForceKernel::PmePreComputation : public CudaContext::ForcePreComputation {
public:
    PmePreComputation(CudaContext& cu, Kernel& pme, CalcPmeReciprocalForceKernel::IO& io) : cu(cu), pme(pme), io(io) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        Vec3 boxSize(cu.getPeriodicBoxSize().x, cu.getPeriodicBoxSize().y, cu.getPeriodicBoxSize().z);
        pme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxSize, includeEnergy);
    }
private:
    CudaContext& cu;
    Kernel pme;
    CalcPmeReciprocalForceKernel::IO& io;
};

class CudaCalcNonbondedForceKernel::PmePostComputation : public CudaContext::ForcePostComputation {
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

class CudaCalcNonbondedForceKernel::SyncStreamPreComputation : public CudaContext::ForcePreComputation {
public:
    SyncStreamPreComputation(CudaContext& cu, CUstream stream, CUevent event, int forceGroup) : cu(cu), stream(stream), event(event), forceGroup(forceGroup) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0) {
            cuEventRecord(event, cu.getCurrentStream());
            cuStreamWaitEvent(stream, event, 0);
        }
    }
private:
    CudaContext& cu;
    CUstream stream;
    CUevent event;
    int forceGroup;
};

class CudaCalcNonbondedForceKernel::SyncStreamPostComputation : public CudaContext::ForcePostComputation {
public:
    SyncStreamPostComputation(CudaContext& cu, CUevent event, int forceGroup) : cu(cu), event(event), forceGroup(forceGroup) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0)
            cuStreamWaitEvent(cu.getCurrentStream(), event, 0);
        return 0.0;
    }
private:
    CudaContext& cu;
    CUevent event;
    int forceGroup;
};

CudaCalcNonbondedForceKernel::~CudaCalcNonbondedForceKernel() {
    cu.setAsCurrent();
    if (sigmaEpsilon != NULL)
        delete sigmaEpsilon;
    if (exceptionParams != NULL)
        delete exceptionParams;
    if (cosSinSums != NULL)
        delete cosSinSums;
    if (directPmeGrid != NULL)
        delete directPmeGrid;
    if (reciprocalPmeGrid != NULL)
        delete reciprocalPmeGrid;
    if (pmeBsplineModuliX != NULL)
        delete pmeBsplineModuliX;
    if (pmeBsplineModuliY != NULL)
        delete pmeBsplineModuliY;
    if (pmeBsplineModuliZ != NULL)
        delete pmeBsplineModuliZ;
    if (pmeAtomRange != NULL)
        delete pmeAtomRange;
    if (pmeAtomGridIndex != NULL)
        delete pmeAtomGridIndex;
    if (sort != NULL)
        delete sort;
    if (pmeio != NULL)
        delete pmeio;
    if (hasInitializedFFT) {
        cufftDestroy(fftForward);
        cufftDestroy(fftBackward);
        if (usePmeStream) {
            cuStreamDestroy(pmeStream);
            cuEventDestroy(pmeSyncEvent);
        }
    }
}

/**
 * Select a size for an FFT that is a multiple of 2, 3, 5, and 7.
 */
static int findFFTDimension(int minimum) {
    if (minimum < 1)
        return 1;
    while (true) {
        // Attempt to factor the current value.

        int unfactored = minimum;
        for (int factor = 2; factor < 8; factor++) {
            while (unfactored > 1 && unfactored%factor == 0)
                unfactored /= factor;
        }
        if (unfactored == 1)
            return minimum;
        minimum++;
    }
}

void CudaCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    cu.setAsCurrent();

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
    sigmaEpsilon = CudaArray::create<float2>(cu, cu.getPaddedNumAtoms(), "sigmaEpsilon");
    CudaArray& posq = cu.getPosq();
    vector<double4> temp(posq.getSize());
    float4* posqf = (float4*) &temp[0];
    double4* posqd = (double4*) &temp[0];
    vector<float2> sigmaEpsilonVector(cu.getPaddedNumAtoms(), make_float2(0, 0));
    vector<vector<int> > exclusionList(numParticles);
    double sumSquaredCharges = 0.0;
    hasCoulomb = false;
    hasLJ = false;
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        if (cu.getUseDoublePrecision())
            posqd[i] = make_double4(0, 0, 0, charge);
        else
            posqf[i] = make_float4(0, 0, 0, (float) charge);
        sigmaEpsilonVector[i] = make_float2((float) (0.5*sigma), (float) (2.0*sqrt(epsilon)));
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
    posq.upload(&temp[0]);
    sigmaEpsilon->upload(sigmaEpsilonVector);
    bool useCutoff = (force.getNonbondedMethod() != NonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != NonbondedForce::NoCutoff && force.getNonbondedMethod() != NonbondedForce::CutoffNonPeriodic);
    map<string, string> defines;
    defines["HAS_COULOMB"] = (hasCoulomb ? "1" : "0");
    defines["HAS_LENNARD_JONES"] = (hasLJ ? "1" : "0");
    defines["USE_LJ_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    if (useCutoff) {
        // Compute the reaction field constants.

        double reactionFieldK = pow(force.getCutoffDistance(), -3.0)*(force.getReactionFieldDielectric()-1.0)/(2.0*force.getReactionFieldDielectric()+1.0);
        double reactionFieldC = (1.0 / force.getCutoffDistance())*(3.0*force.getReactionFieldDielectric())/(2.0*force.getReactionFieldDielectric()+1.0);
        defines["REACTION_FIELD_K"] = cu.doubleToString(reactionFieldK);
        defines["REACTION_FIELD_C"] = cu.doubleToString(reactionFieldC);
        
        // Compute the switching coefficients.
        
        if (force.getUseSwitchingFunction()) {
            defines["LJ_SWITCH_CUTOFF"] = cu.doubleToString(force.getSwitchingDistance());
            defines["LJ_SWITCH_C3"] = cu.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
            defines["LJ_SWITCH_C4"] = cu.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
            defines["LJ_SWITCH_C5"] = cu.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
        }
    }
    if (force.getUseDispersionCorrection() && cu.getContextIndex() == 0)
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
    alpha = 0;
    ewaldSelfEnergy = 0.0;
    if (force.getNonbondedMethod() == NonbondedForce::Ewald) {
        // Compute the Ewald parameters.

        int kmaxx, kmaxy, kmaxz;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        if (cu.getContextIndex() == 0) {
            ewaldSelfEnergy = -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI);

            // Create the reciprocal space kernels.

            map<string, string> replacements;
            replacements["NUM_ATOMS"] = cu.intToString(numParticles);
            replacements["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
            replacements["KMAX_X"] = cu.intToString(kmaxx);
            replacements["KMAX_Y"] = cu.intToString(kmaxy);
            replacements["KMAX_Z"] = cu.intToString(kmaxz);
            replacements["EXP_COEFFICIENT"] = cu.doubleToString(-1.0/(4.0*alpha*alpha));
            replacements["ONE_4PI_EPS0"] = cu.doubleToString(ONE_4PI_EPS0);
            replacements["M_PI"] = cu.doubleToString(M_PI);
            CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaKernelSources::ewald, replacements);
            ewaldSumsKernel = cu.getKernel(module, "calculateEwaldCosSinSums");
            ewaldForcesKernel = cu.getKernel(module, "calculateEwaldForces");
            int elementSize = (cu.getUseDoublePrecision() ? sizeof(double2) : sizeof(float2));
            cosSinSums = new CudaArray(cu, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), elementSize, "cosSinSums");
        }
    }
    else if (force.getNonbondedMethod() == NonbondedForce::PME) {
        // Compute the PME parameters.

        int gridSizeX, gridSizeY, gridSizeZ;

        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ);
        gridSizeX = findFFTDimension(gridSizeX);
        gridSizeY = findFFTDimension(gridSizeY);
        gridSizeZ = findFFTDimension(gridSizeZ);

        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        if (cu.getContextIndex() == 0) {
            ewaldSelfEnergy = -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI);
            pmeDefines["PME_ORDER"] = cu.intToString(PmeOrder);
            pmeDefines["NUM_ATOMS"] = cu.intToString(numParticles);
            pmeDefines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
            pmeDefines["RECIP_EXP_FACTOR"] = cu.doubleToString(M_PI*M_PI/(alpha*alpha));
            pmeDefines["GRID_SIZE_X"] = cu.intToString(gridSizeX);
            pmeDefines["GRID_SIZE_Y"] = cu.intToString(gridSizeY);
            pmeDefines["GRID_SIZE_Z"] = cu.intToString(gridSizeZ);
            pmeDefines["EPSILON_FACTOR"] = cu.doubleToString(sqrt(ONE_4PI_EPS0));
            pmeDefines["M_PI"] = cu.doubleToString(M_PI);
            if (cu.getUseDoublePrecision())
                pmeDefines["USE_DOUBLE_PRECISION"] = "1";
            CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaKernelSources::pme, pmeDefines);
            if (cu.getPlatformData().useCpuPme) {
                // Create the CPU PME kernel.

                try {
                    cpuPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), *cu.getPlatformData().context);
                    cpuPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSizeX, gridSizeY, gridSizeZ, numParticles, alpha);
                    CUfunction addForcesKernel = cu.getKernel(module, "addForces");
                    pmeio = new PmeIO(cu, addForcesKernel);
                    cu.addPreComputation(new PmePreComputation(cu, cpuPme, *pmeio));
                    cu.addPostComputation(new PmePostComputation(cpuPme, *pmeio));
                }
                catch (OpenMMException& ex) {
                    // The CPU PME plugin isn't available.
                }
            }
            if (pmeio == NULL) {
                pmeGridIndexKernel = cu.getKernel(module, "findAtomGridIndex");
                pmeSpreadChargeKernel = cu.getKernel(module, "gridSpreadCharge");
                pmeConvolutionKernel = cu.getKernel(module, "reciprocalConvolution");
                pmeInterpolateForceKernel = cu.getKernel(module, "gridInterpolateForce");
                pmeEvalEnergyKernel = cu.getKernel(module, "gridEvaluateEnergy");
                pmeFinishSpreadChargeKernel = cu.getKernel(module, "finishSpreadCharge");
                cuFuncSetCacheConfig(pmeSpreadChargeKernel, CU_FUNC_CACHE_PREFER_L1);
                cuFuncSetCacheConfig(pmeInterpolateForceKernel, CU_FUNC_CACHE_PREFER_L1);

                // Create required data structures.

                int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));

                directPmeGrid = new CudaArray(cu, gridSizeX*gridSizeY*gridSizeZ, cu.getComputeCapability() >= 2.0 ? elementSize : sizeof(long long), "originalPmeGrid");
                reciprocalPmeGrid = new CudaArray(cu, gridSizeX*gridSizeY*(gridSizeZ/2+1), 2*elementSize, "reciprocalPmeGrid");

                cu.addAutoclearBuffer(*directPmeGrid);

                pmeBsplineModuliX = new CudaArray(cu, gridSizeX, elementSize, "pmeBsplineModuliX");
                pmeBsplineModuliY = new CudaArray(cu, gridSizeY, elementSize, "pmeBsplineModuliY");
                pmeBsplineModuliZ = new CudaArray(cu, gridSizeZ, elementSize, "pmeBsplineModuliZ");
                pmeAtomRange = CudaArray::create<int>(cu, gridSizeX*gridSizeY*gridSizeZ+1, "pmeAtomRange");
                pmeAtomGridIndex = CudaArray::create<int2>(cu, numParticles, "pmeAtomGridIndex");
                sort = new CudaSort(cu, new SortTrait(), cu.getNumAtoms());

                cufftResult result = cufftPlan3d(&fftForward, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_D2Z : CUFFT_R2C);
                if (result != CUFFT_SUCCESS)
                    throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
                result = cufftPlan3d(&fftBackward, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_Z2D : CUFFT_C2R);
                if (result != CUFFT_SUCCESS)
                    throw OpenMMException("Error initializing FFT: "+cu.intToString(result));

                cufftSetCompatibilityMode(fftForward, CUFFT_COMPATIBILITY_NATIVE);
                cufftSetCompatibilityMode(fftBackward, CUFFT_COMPATIBILITY_NATIVE);
                
                // Prepare for doing PME on its own stream.
                
                int cufftVersion;
                cufftGetVersion(&cufftVersion);
                usePmeStream = (cu.getComputeCapability() < 5.0 && numParticles < 130000 && cufftVersion >= 6000); // Workarounds for various CUDA bugs
                if (usePmeStream) {
                    cuStreamCreate(&pmeStream, CU_STREAM_NON_BLOCKING);
                    cufftSetStream(fftForward, pmeStream);
                    cufftSetStream(fftBackward, pmeStream);
                    CHECK_RESULT(cuEventCreate(&pmeSyncEvent, CU_EVENT_DISABLE_TIMING), "Error creating event for NonbondedForce");
                    int recipForceGroup = force.getReciprocalSpaceForceGroup();
                    if (recipForceGroup < 0)
                        recipForceGroup = force.getForceGroup();
                    cu.addPreComputation(new SyncStreamPreComputation(cu, pmeStream, pmeSyncEvent, recipForceGroup));
                    cu.addPostComputation(new SyncStreamPostComputation(cu, pmeSyncEvent, recipForceGroup));
                }
                hasInitializedFFT = true;

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
                    vector<double> moduli(ndata);
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
                        if (moduli[i] < 1.0e-7)
                            moduli[i] = (moduli[i-1]+moduli[i+1])*0.5;
                    if (cu.getUseDoublePrecision()) {
                        if (dim == 0)
                            pmeBsplineModuliX->upload(moduli);
                        else if (dim == 1)
                            pmeBsplineModuliY->upload(moduli);
                        else
                            pmeBsplineModuliZ->upload(moduli);
                    }
                    else {
                        vector<float> modulif(ndata);
                        for (int i = 0; i < ndata; i++)
                            modulif[i] = (float) moduli[i];
                        if (dim == 0)
                            pmeBsplineModuliX->upload(modulif);
                        else if (dim == 1)
                            pmeBsplineModuliY->upload(modulif);
                        else
                            pmeBsplineModuliZ->upload(modulif);
                    }
                }
            }
        }
    }

    // Add the interaction to the default nonbonded kernel.
   
    string source = cu.replaceStrings(CudaKernelSources::coulombLennardJones, defines);
    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
    if (hasLJ)
        cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo("sigmaEpsilon", "float", 2, sizeof(float2), sigmaEpsilon->getDevicePointer()));

    // Initialize the exceptions.

    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    if (numExceptions > 0) {
        exceptionAtoms.resize(numExceptions);
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        exceptionParams = CudaArray::create<float4>(cu, numExceptions, "exceptionParams");
        vector<float4> exceptionParamsVector(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            exceptionParamsVector[i] = make_float4((float) (ONE_4PI_EPS0*chargeProd), (float) sigma, (float) (4.0*epsilon), 0.0f);
            exceptionAtoms[i] = make_pair(atoms[i][0], atoms[i][1]);
        }
        exceptionParams->upload(exceptionParamsVector);
        map<string, string> replacements;
        replacements["PARAMS"] = cu.getBondedUtilities().addArgument(exceptionParams->getDevicePointer(), "float4");
        cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::nonbondedExceptions, replacements), force.getForceGroup());
    }
    cu.addForce(new CudaNonbondedForceInfo(force));
}

double CudaCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    if (cosSinSums != NULL && includeReciprocal) {
        void* sumsArgs[] = {&cu.getEnergyBuffer().getDevicePointer(), &cu.getPosq().getDevicePointer(), &cosSinSums->getDevicePointer(), cu.getPeriodicBoxSizePointer()};
        cu.executeKernel(ewaldSumsKernel, sumsArgs, cosSinSums->getSize());
        void* forcesArgs[] = {&cu.getForce().getDevicePointer(), &cu.getPosq().getDevicePointer(), &cosSinSums->getDevicePointer(), cu.getPeriodicBoxSizePointer()};
        cu.executeKernel(ewaldForcesKernel, forcesArgs, cu.getNumAtoms());
    }
    if (directPmeGrid != NULL && includeReciprocal) {
        if (usePmeStream)
            cu.setCurrentStream(pmeStream);
        void* gridIndexArgs[] = {&cu.getPosq().getDevicePointer(), &pmeAtomGridIndex->getDevicePointer(), cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeGridIndexKernel, gridIndexArgs, cu.getNumAtoms());

        sort->sort(*pmeAtomGridIndex);

        void* spreadArgs[] = {&cu.getPosq().getDevicePointer(), &directPmeGrid->getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                cu.getInvPeriodicBoxSizePointer(), &pmeAtomGridIndex->getDevicePointer()};
        cu.executeKernel(pmeSpreadChargeKernel, spreadArgs, cu.getNumAtoms(), 128);

        if (cu.getUseDoublePrecision() || cu.getComputeCapability() < 2.0) {
            void* finishSpreadArgs[] = {&directPmeGrid->getDevicePointer()};
            cu.executeKernel(pmeFinishSpreadChargeKernel, finishSpreadArgs, directPmeGrid->getSize());
        }

        if (cu.getUseDoublePrecision())
            cufftExecD2Z(fftForward, (double*) directPmeGrid->getDevicePointer(), (double2*) reciprocalPmeGrid->getDevicePointer());
        else
            cufftExecR2C(fftForward, (float*) directPmeGrid->getDevicePointer(), (float2*) reciprocalPmeGrid->getDevicePointer());

        if (includeEnergy) {
            void* computeEnergyArgs[] = {&reciprocalPmeGrid->getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(), &pmeBsplineModuliX->getDevicePointer(), &pmeBsplineModuliY->getDevicePointer(), &pmeBsplineModuliZ->getDevicePointer(), cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
            cu.executeKernel(pmeEvalEnergyKernel, computeEnergyArgs, cu.getNumAtoms());
        }

        void* convolutionArgs[] = {&reciprocalPmeGrid->getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(), &pmeBsplineModuliX->getDevicePointer(), &pmeBsplineModuliY->getDevicePointer(), &pmeBsplineModuliZ->getDevicePointer(), cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer()};
        cu.executeKernel(pmeConvolutionKernel, convolutionArgs, cu.getNumAtoms());

        if (cu.getUseDoublePrecision())
            cufftExecZ2D(fftBackward, (double2*) reciprocalPmeGrid->getDevicePointer(), (double*) directPmeGrid->getDevicePointer());
        else
            cufftExecC2R(fftBackward, (float2*) reciprocalPmeGrid->getDevicePointer(), (float*)  directPmeGrid->getDevicePointer());


        void* interpolateArgs[] = {&cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(), &directPmeGrid->getDevicePointer(),
                cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer(), &pmeAtomGridIndex->getDevicePointer()};
        cu.executeKernel(pmeInterpolateForceKernel, interpolateArgs, cu.getNumAtoms(), 128);
        if (usePmeStream) {
            cuEventRecord(pmeSyncEvent, pmeStream);
            cu.restoreDefaultStream();
        }
    }
    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
    if (dispersionCoefficient != 0.0 && includeDirect) {
        double4 boxSize = cu.getPeriodicBoxSize();
        energy += dispersionCoefficient/(boxSize.x*boxSize.y*boxSize.z);
    }
    return energy;
}

void CudaCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
    // Make sure the new parameters are acceptable.
    
    cu.setAsCurrent();
    if (force.getNumParticles() != cu.getNumAtoms())
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
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    
    // Record the per-particle parameters.
    
    CudaArray& posq = cu.getPosq();
    posq.download(cu.getPinnedBuffer());
    float4* posqf = (float4*) cu.getPinnedBuffer();
    double4* posqd = (double4*) cu.getPinnedBuffer();
    vector<float2> sigmaEpsilonVector(cu.getPaddedNumAtoms(), make_float2(0, 0));
    double sumSquaredCharges = 0.0;
    const vector<int>& order = cu.getAtomIndex();
    for (int i = 0; i < force.getNumParticles(); i++) {
        int index = order[i];
        double charge, sigma, epsilon;
        force.getParticleParameters(index, charge, sigma, epsilon);
        if (cu.getUseDoublePrecision())
            posqd[i].w = charge;
        else
            posqf[i].w = (float) charge;
        sigmaEpsilonVector[index] = make_float2((float) (0.5*sigma), (float) (2.0*sqrt(epsilon)));
        sumSquaredCharges += charge*charge;
    }
    posq.upload(cu.getPinnedBuffer());
    sigmaEpsilon->upload(sigmaEpsilonVector);
    
    // Record the exceptions.
    
    if (numExceptions > 0) {
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        vector<float4> exceptionParamsVector(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            exceptionParamsVector[i] = make_float4((float) (ONE_4PI_EPS0*chargeProd), (float) sigma, (float) (4.0*epsilon), 0.0f);
        }
        exceptionParams->upload(exceptionParamsVector);
    }
    
    // Compute other values.
    
    NonbondedForce::NonbondedMethod method = force.getNonbondedMethod();
    if (method == NonbondedForce::Ewald || method == NonbondedForce::PME)
        ewaldSelfEnergy = (cu.getContextIndex() == 0 ? -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI) : 0.0);
    if (force.getUseDispersionCorrection() && cu.getContextIndex() == 0 && (method == NonbondedForce::CutoffPeriodic || method == NonbondedForce::Ewald || method == NonbondedForce::PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
    cu.invalidateMolecules();
}

class CudaCustomNonbondedForceInfo : public CudaForceInfo {
public:
    CudaCustomNonbondedForceInfo(const CustomNonbondedForce& force) : force(force) {
        if (force.getNumInteractionGroups() > 0) {
            groupsForParticle.resize(force.getNumParticles());
            for (int i = 0; i < force.getNumInteractionGroups(); i++) {
                set<int> set1, set2;
                force.getInteractionGroupParameters(i, set1, set2);
                for (set<int>::const_iterator iter = set1.begin(); iter != set1.end(); ++iter)
                    groupsForParticle[*iter].insert(2*i);
                for (set<int>::const_iterator iter = set2.begin(); iter != set2.end(); ++iter)
                    groupsForParticle[*iter].insert(2*i+1);
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

CudaCalcCustomNonbondedForceKernel::~CudaCalcCustomNonbondedForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
    if (interactionGroupData != NULL)
        delete interactionGroupData;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
    if (forceCopy != NULL)
        delete forceCopy;
}

void CudaCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    cu.setAsCurrent();
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = (force.getNumInteractionGroups() == 0 ? "custom"+cu.intToString(forceIndex)+"_" : "");

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customNonbondedParameters");
    if (force.getNumGlobalParameters() > 0)
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customNonbondedGlobals");
    vector<vector<float> > paramVector(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
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
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions.push_back(CudaArray::create<float>(cu, f.size(), "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName, "float", width, width*sizeof(float), tabulatedFunctions[tabulatedFunctions.size()-1]->getDevicePointer()));
    }

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
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
        string value = "globals["+cu.intToString(i)+"]";
        variables.push_back(makeVariable(name, prefix+value));
    }
    stringstream compute;
    compute << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, prefix+"temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["USE_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    if (force.getUseSwitchingFunction()) {
        // Compute the switching coefficients.
        
        replacements["SWITCH_CUTOFF"] = cu.doubleToString(force.getSwitchingDistance());
        replacements["SWITCH_C3"] = cu.doubleToString(10/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 3.0));
        replacements["SWITCH_C4"] = cu.doubleToString(15/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 4.0));
        replacements["SWITCH_C5"] = cu.doubleToString(6/pow(force.getSwitchingDistance()-force.getCutoffDistance(), 5.0));
    }
    string source = cu.replaceStrings(CudaKernelSources::customNonbonded, replacements);
    if (force.getNumInteractionGroups() > 0)
        initInteractionGroups(force, source);
    else {
        cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"params"+cu.intToString(i+1), buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        if (globals != NULL) {
            globals->upload(globalParamValues);
            cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(float), globals->getDevicePointer()));
        }
    }
    cu.addForce(new CudaCustomNonbondedForceInfo(force));
    
    // Record information for the long range correction.
    
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic && force.getUseLongRangeCorrection() && cu.getContextIndex() == 0) {
        forceCopy = new CustomNonbondedForce(force);
        hasInitializedLongRangeCorrection = false;
    }
    else {
        longRangeCoefficient = 0.0;
        hasInitializedLongRangeCorrection = true;
    }
}

void CudaCalcCustomNonbondedForceKernel::initInteractionGroups(const CustomNonbondedForce& force, const string& interactionSource) {
    // Process groups to form tiles.
    
    vector<vector<int> > atomLists;
    vector<pair<int, int> > tiles;
    map<pair<int, int>, int> duplicateInteractions;
    for (int group = 0; group < force.getNumInteractionGroups(); group++) {
        // Get the list of atoms in this group and sort them.
        
        set<int> set1, set2;
        force.getInteractionGroupParameters(group, set1, set2);
        vector<int> atoms1, atoms2;
        atoms1.insert(atoms1.begin(), set1.begin(), set1.end());
        atoms2.insert(atoms2.begin(), set2.begin(), set2.end());
        sort(atoms1.begin(), atoms1.end());
        sort(atoms2.begin(), atoms2.end());
        
        // Find how many tiles we will create for this group.
        
        int tileWidth = min(min(32, (int) atoms1.size()), (int) atoms2.size());
        if (tileWidth == 0)
            continue;
        int numBlocks1 = (atoms1.size()+tileWidth-1)/tileWidth;
        int numBlocks2 = (atoms2.size()+tileWidth-1)/tileWidth;
        
        // Add the tiles.
        
        for (int i = 0; i < numBlocks1; i++)
            for (int j = 0; j < numBlocks2; j++)
                tiles.push_back(make_pair(atomLists.size()+i, atomLists.size()+numBlocks1+j));
        
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
        
        // If this group contains duplicate interactions, record that we need to skip them once.
        
        for (int i = 0; i < (int) atoms1.size(); i++) {
            int a1 = atoms1[i];
            if (set2.find(a1) == set2.end())
                continue;
            for (int j = 0; j < (int) atoms2.size() && atoms2[j] < a1; j++) {
                int a2 = atoms2[j];
                if (set1.find(a2) != set1.end()) {
                    pair<int, int> key = make_pair(a2, a1);
                    if (duplicateInteractions.find(key) == duplicateInteractions.end())
                        duplicateInteractions[key] = 0;
                    duplicateInteractions[key]++;
                }
            }
        }
    }
    
    // Build a lookup table for quickly identifying excluded interactions.
    
    set<pair<int, int> > exclusions;
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int p1, p2;
        force.getExclusionParticles(i, p1, p2);
        exclusions.insert(make_pair(min(p1, p2), max(p1, p2)));
    }
    
    // Build the exclusion flags for each tile.  While we're at it, filter out tiles
    // where all interactions are excluded, and sort the tiles by size.

    vector<vector<int> > exclusionFlags(tiles.size());
    vector<pair<int, int> > tileOrder;
    for (int tile = 0; tile < tiles.size(); tile++) {
        if (atomLists[tiles[tile].first].size() < atomLists[tiles[tile].second].size()) {
            // For efficiency, we want the first axis to be the larger one.
            
            int swap = tiles[tile].first;
            tiles[tile].first = tiles[tile].second;
            tiles[tile].second = swap;
        }
        vector<int>& atoms1 = atomLists[tiles[tile].first];
        vector<int>& atoms2 = atomLists[tiles[tile].second];
        vector<int> flags(atoms1.size(), (int) (1LL<<atoms2.size())-1);
        int numExcluded = 0;
        for (int i = 0; i < (int) atoms1.size(); i++)
            for (int j = 0; j < (int) atoms2.size(); j++) {
                int a1 = atoms1[i];
                int a2 = atoms2[j];
                bool isExcluded = false;
                pair<int, int> key = make_pair(min(a1, a2), max(a1, a2));
                if (a1 == a2 || exclusions.find(key) != exclusions.end())
                    isExcluded = true; // This is an excluded interaction.
                else if (duplicateInteractions.find(key) != duplicateInteractions.end() && duplicateInteractions[key] > 0) {
                    // Both atoms are in both sets, so skip duplicate interactions.
                    
                    isExcluded = true;
                    duplicateInteractions[key]--;
                }
                if (isExcluded) {
                    flags[i] &= -1-(1<<j);
                    numExcluded++;
                }
            }
        if (numExcluded == atoms1.size()*atoms2.size())
            continue; // All interactions are excluded.
        tileOrder.push_back(make_pair((int) -atoms2.size(), tile));
        exclusionFlags[tile] = flags;
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
    vector<int4> groupData;
    for (int tileSet = 0; tileSet < numTileSets; tileSet++) {
        int indexInTileSet = 0;
        for (int i = tileSetStart[tileSet]; i < tileSetStart[tileSet+1]; i++) {
            int tile = tileOrder[i].second;
            vector<int>& atoms1 = atomLists[tiles[tile].first];
            vector<int>& atoms2 = atomLists[tiles[tile].second];
            int range = indexInTileSet + ((indexInTileSet+atoms1.size())<<16);
            int allFlags = (1<<atoms2.size())-1;
            for (int j = 0; j < (int) atoms1.size(); j++) {
                int a1 = atoms1[j];
                int a2 = (j < atoms2.size() ? atoms2[j] : 0);
                int flags = (exclusionFlags[tile].size() > 0 ? exclusionFlags[tile][j] : allFlags);
                groupData.push_back(make_int4(a1, a2, range, flags<<indexInTileSet));
            }
            indexInTileSet += atoms1.size();
        }
        for (; indexInTileSet < 32; indexInTileSet++)
            groupData.push_back(make_int4(0, 0, 0, 0));
    }
    interactionGroupData = CudaArray::create<int4>(cu, groupData.size(), "interactionGroupData");
    interactionGroupData->upload(groupData);
    
    // Create the kernel.
    
    map<string, string> replacements;
    replacements["COMPUTE_INTERACTION"] = interactionSource;
    const string suffixes[] = {"x", "y", "z", "w"};
    stringstream localData;
    int localDataSize = 0;
    vector<CudaNonbondedUtilities::ParameterInfo>& buffers = params->getBuffers(); 
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
        args<<", const "<<buffers[i].getType()<<"* __restrict__ global_params"<<(i+1);
    if (globals != NULL)
        args<<", const float* __restrict__ globals";
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    stringstream load1;
    for (int i = 0; i < (int) buffers.size(); i++)
        load1<<buffers[i].getType()<<" params"<<(i+1)<<"1 = global_params"<<(i+1)<<"[atom1];\n";
    replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
    stringstream loadLocal2;
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getNumComponents() == 1)
            loadLocal2<<"localData[threadIdx.x].params"<<(i+1)<<" = global_params"<<(i+1)<<"[atom2];\n";
        else {
            loadLocal2<<buffers[i].getType()<<" temp_params"<<(i+1)<<" = global_params"<<(i+1)<<"[atom2];\n";
            for (int j = 0; j < buffers[i].getNumComponents(); ++j)
                loadLocal2<<"localData[threadIdx.x].params"<<(i+1)<<"_"<<suffixes[j]<<" = temp_params"<<(i+1)<<"."<<suffixes[j]<<";\n";
        }
    }
    replacements["LOAD_LOCAL_PARAMETERS"] = loadLocal2.str();
    stringstream load2;
    for (int i = 0; i < (int) buffers.size(); i++) {
        if (buffers[i].getNumComponents() == 1)
            load2<<buffers[i].getType()<<" params"<<(i+1)<<"2 = localData[localIndex].params"<<(i+1)<<";\n";
        else {
            load2<<buffers[i].getType()<<" params"<<(i+1)<<"2 = make_"<<buffers[i].getType()<<"(";
            for (int j = 0; j < buffers[i].getNumComponents(); ++j) {
                if (j > 0)
                    load2<<", ";
                load2<<"localData[localIndex].params"<<(i+1)<<"_"<<suffixes[j];
            }
            load2<<");\n";
        }
    }
    replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
    map<string, string> defines;
    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
        defines["USE_CUTOFF"] = "1";
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic)
        defines["USE_PERIODIC"] = "1";
    defines["LOCAL_MEMORY_SIZE"] = cu.intToString(max(32, cu.getNonbondedUtilities().getForceThreadBlockSize()));
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["TILE_SIZE"] = "32";
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*numTileSets/numContexts;
    int endIndex = (cu.getContextIndex()+1)*numTileSets/numContexts;
    defines["FIRST_TILE"] = cu.intToString(startIndex);
    defines["LAST_TILE"] = cu.intToString(endIndex);
    if ((localDataSize/4)%2 == 0 && !cu.getUseDoublePrecision())
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    CUmodule program = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::customNonbondedGroups, replacements), defines);
    interactionGroupKernel = cu.getKernel(program, "computeInteractionGroups");
    numGroupThreadBlocks = cu.getNonbondedUtilities().getNumForceThreadBlocks();
}

double CudaCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed) {
            globals->upload(globalParamValues);
            if (forceCopy != NULL) {
                longRangeCoefficient = CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner());
                hasInitializedLongRangeCorrection = true;
            }
        }
    }
    if (!hasInitializedLongRangeCorrection) {
        longRangeCoefficient = CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner());
        hasInitializedLongRangeCorrection = true;
    }
    if (interactionGroupData != NULL) {
        if (!hasInitializedKernel) {
            hasInitializedKernel = true;
            interactionGroupArgs.push_back(&cu.getForce().getDevicePointer());
            interactionGroupArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
            interactionGroupArgs.push_back(&cu.getPosq().getDevicePointer());
            interactionGroupArgs.push_back(&interactionGroupData->getDevicePointer());
            interactionGroupArgs.push_back(cu.getPeriodicBoxSizePointer());
            interactionGroupArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            for (int i = 0; i < (int) params->getBuffers().size(); i++)
                interactionGroupArgs.push_back(&params->getBuffers()[i].getMemory());
            if (globals != NULL)
                interactionGroupArgs.push_back(&globals->getDevicePointer());
        }
        int forceThreadBlockSize = cu.getNonbondedUtilities().getForceThreadBlockSize();
        cu.executeKernel(interactionGroupKernel, &interactionGroupArgs[0], numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    double4 boxSize = cu.getPeriodicBoxSize();
    return longRangeCoefficient/(boxSize.x*boxSize.y*boxSize.z);
}

void CudaCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
    cu.setAsCurrent();
    int numParticles = force.getNumParticles();
    if (numParticles != cu.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<vector<float> > paramVector(numParticles);
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // If necessary, recompute the long range correction.
    
    if (forceCopy != NULL) {
        longRangeCoefficient = CustomNonbondedForceImpl::calcLongRangeCorrection(force, context.getOwner());
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaGBSAOBCForceInfo : public CudaForceInfo {
public:
    CudaGBSAOBCForceInfo(const GBSAOBCForce& force) : force(force) {
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

CudaCalcGBSAOBCForceKernel::~CudaCalcGBSAOBCForceKernel() {
    cu.setAsCurrent();
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

void CudaCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    cu.setAsCurrent();
    if (cu.getPlatformData().contexts.size() > 1)
        throw OpenMMException("GBSAOBCForce does not support using multiple CUDA devices");
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    params = CudaArray::create<float2>(cu, cu.getPaddedNumAtoms(), "gbsaObcParams");
    if (cu.getUseDoublePrecision()) {
        bornRadii = CudaArray::create<double>(cu, cu.getPaddedNumAtoms(), "bornRadii");
        obcChain = CudaArray::create<double>(cu, cu.getPaddedNumAtoms(), "obcChain");
    }
    else {
        bornRadii = CudaArray::create<float>(cu, cu.getPaddedNumAtoms(), "bornRadii");
        obcChain = CudaArray::create<float>(cu, cu.getPaddedNumAtoms(), "obcChain");
    }
    bornSum = CudaArray::create<long long>(cu, cu.getPaddedNumAtoms(), "bornSum");
    bornForce = CudaArray::create<long long>(cu, cu.getPaddedNumAtoms(), "bornForce");
    cu.addAutoclearBuffer(*bornSum);
    cu.addAutoclearBuffer(*bornForce);
    CudaArray& posq = cu.getPosq();
    vector<double4> temp(posq.getSize());
    float4* posqf = (float4*) &temp[0];
    double4* posqd = (double4*) &temp[0];
    vector<float2> paramsVector(cu.getPaddedNumAtoms(), make_float2(1, 1));
    const double dielectricOffset = 0.009;
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        paramsVector[i] = make_float2((float) radius, (float) (scalingFactor*radius));
        if (cu.getUseDoublePrecision())
            posqd[i] = make_double4(0, 0, 0, charge);
        else
            posqf[i] = make_float4(0, 0, 0, (float) charge);
    }
    posq.upload(&temp[0]);
    params->upload(paramsVector);
    prefactor = -ONE_4PI_EPS0*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
    surfaceAreaFactor = -6.0*4*M_PI*force.getSurfaceAreaEnergy();
    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
    string source = CudaKernelSources::gbsaObc2;
    nb.addInteraction(useCutoff, usePeriodic, false, force.getCutoffDistance(), vector<vector<int> >(), source, force.getForceGroup());
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("obcParams", "float", 2, sizeof(float2), params->getDevicePointer()));;
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("bornForce", "long long", 1, sizeof(long long), bornForce->getDevicePointer()));;
    cu.addForce(new CudaGBSAOBCForceInfo(force));
}

double CudaCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    if (!hasCreatedKernels) {
        // These Kernels cannot be created in initialize(), because the CudaNonbondedUtilities has not been initialized yet then.

        hasCreatedKernels = true;
        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
        map<string, string> defines;
        if (nb.getUseCutoff())
            defines["USE_CUTOFF"] = "1";
        if (nb.getUsePeriodic())
            defines["USE_PERIODIC"] = "1";
        if (cu.getComputeCapability() >= 3.0 && !cu.getUseDoublePrecision())
            defines["ENABLE_SHUFFLE"] = "1";
        defines["CUTOFF_SQUARED"] = cu.doubleToString(nb.getCutoffDistance()*nb.getCutoffDistance());
        defines["CUTOFF"] = cu.doubleToString(nb.getCutoffDistance());
        defines["PREFACTOR"] = cu.doubleToString(prefactor);
        defines["SURFACE_AREA_FACTOR"] = cu.doubleToString(surfaceAreaFactor);
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        defines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
        defines["FORCE_WORK_GROUP_SIZE"] = cu.intToString(nb.getForceThreadBlockSize());
        defines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
        int numExclusionTiles = nb.getExclusionTiles().getSize();
        defines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(numExclusionTiles);
        int numContexts = cu.getPlatformData().contexts.size();
        int startExclusionIndex = cu.getContextIndex()*numExclusionTiles/numContexts;
        int endExclusionIndex = (cu.getContextIndex()+1)*numExclusionTiles/numContexts;
        defines["FIRST_EXCLUSION_TILE"] = cu.intToString(startExclusionIndex);
        defines["LAST_EXCLUSION_TILE"] = cu.intToString(endExclusionIndex);
        map<string, string> replacements;
        CUmodule module = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::gbsaObc1, replacements), defines);
        computeBornSumKernel = cu.getKernel(module, "computeBornSum");
        computeSumArgs.push_back(&bornSum->getDevicePointer());
        computeSumArgs.push_back(&cu.getPosq().getDevicePointer());
        computeSumArgs.push_back(&params->getDevicePointer());
        if (nb.getUseCutoff()) {
            computeSumArgs.push_back(&nb.getInteractingTiles().getDevicePointer());
            computeSumArgs.push_back(&nb.getInteractionCount().getDevicePointer());
            computeSumArgs.push_back(cu.getPeriodicBoxSizePointer());
            computeSumArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            computeSumArgs.push_back(&maxTiles);
            computeSumArgs.push_back(&nb.getBlockCenters().getDevicePointer());
            computeSumArgs.push_back(&nb.getBlockBoundingBoxes().getDevicePointer());
            computeSumArgs.push_back(&nb.getInteractingAtoms().getDevicePointer());
        }
        else
            computeSumArgs.push_back(&maxTiles);
        computeSumArgs.push_back(&nb.getExclusionTiles().getDevicePointer());
        force1Kernel = cu.getKernel(module, "computeGBSAForce1");
        force1Args.push_back(&cu.getForce().getDevicePointer());
        force1Args.push_back(&bornForce->getDevicePointer());
        force1Args.push_back(&cu.getEnergyBuffer().getDevicePointer());
        force1Args.push_back(&cu.getPosq().getDevicePointer());
        force1Args.push_back(&bornRadii->getDevicePointer());
        if (nb.getUseCutoff()) {
            force1Args.push_back(&nb.getInteractingTiles().getDevicePointer());
            force1Args.push_back(&nb.getInteractionCount().getDevicePointer());
            force1Args.push_back(cu.getPeriodicBoxSizePointer());
            force1Args.push_back(cu.getInvPeriodicBoxSizePointer());
            force1Args.push_back(&maxTiles);
            force1Args.push_back(&nb.getBlockCenters().getDevicePointer());
            force1Args.push_back(&nb.getBlockBoundingBoxes().getDevicePointer());
            force1Args.push_back(&nb.getInteractingAtoms().getDevicePointer());
        }
        else
            force1Args.push_back(&maxTiles);
        force1Args.push_back(&nb.getExclusionTiles().getDevicePointer());
        reduceBornSumKernel = cu.getKernel(module, "reduceBornSum");
        reduceBornForceKernel = cu.getKernel(module, "reduceBornForce");
    }
    if (nb.getUseCutoff()) {
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            computeSumArgs[3] = &nb.getInteractingTiles().getDevicePointer();
            force1Args[5] = &nb.getInteractingTiles().getDevicePointer();
            computeSumArgs[10] = &nb.getInteractingAtoms().getDevicePointer();
            force1Args[12] = &nb.getInteractingAtoms().getDevicePointer();
        }
    }
    cu.executeKernel(computeBornSumKernel, &computeSumArgs[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    float alpha = 1.0f, beta = 0.8f, gamma = 4.85f;
    void* reduceSumArgs[] = {&alpha, &beta, &gamma, &bornSum->getDevicePointer(), &params->getDevicePointer(),
            &bornRadii->getDevicePointer(), &obcChain->getDevicePointer()};
    cu.executeKernel(reduceBornSumKernel, reduceSumArgs, cu.getPaddedNumAtoms());
    cu.executeKernel(force1Kernel, &force1Args[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    void* reduceForceArgs[] = {&bornForce->getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(), &params->getDevicePointer(),
            &bornRadii->getDevicePointer(), &obcChain->getDevicePointer()};
    cu.executeKernel(reduceBornForceKernel, &reduceForceArgs[0], cu.getPaddedNumAtoms());
    return 0.0;
}

void CudaCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
    // Make sure the new parameters are acceptable.
    
    cu.setAsCurrent();
    int numParticles = force.getNumParticles();
    if (numParticles != cu.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    CudaArray& posq = cu.getPosq();
    float4* posqf = (float4*) cu.getPinnedBuffer();
    double4* posqd = (double4*) cu.getPinnedBuffer();
    posq.download(cu.getPinnedBuffer());
    vector<float2> paramsVector(cu.getPaddedNumAtoms(), make_float2(1, 1));
    const double dielectricOffset = 0.009;
    for (int i = 0; i < numParticles; i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        paramsVector[i] = make_float2((float) radius, (float) (scalingFactor*radius));
        if (cu.getUseDoublePrecision())
            posqd[i].w = charge;
        else
            posqf[i].w = (float) charge;
    }
    posq.upload(cu.getPinnedBuffer());
    params->upload(paramsVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCustomGBForceInfo : public CudaForceInfo {
public:
    CudaCustomGBForceInfo(const CustomGBForce& force) : force(force) {
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

CudaCalcCustomGBForceKernel::~CudaCalcCustomGBForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (computedValues != NULL)
        delete computedValues;
    if (energyDerivs != NULL)
        delete energyDerivs;
    if (energyDerivChain != NULL)
        delete energyDerivChain;
    if (longEnergyDerivs != NULL)
        delete longEnergyDerivs;
    if (globals != NULL)
        delete globals;
    if (valueBuffers != NULL)
        delete valueBuffers;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
}

void CudaCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    cu.setAsCurrent();
    if (cu.getPlatformData().contexts.size() > 1)
        throw OpenMMException("CustomGBForce does not support using multiple CUDA devices");
    bool useExclusionsForValue = false;
    numComputedValues = force.getNumComputedValues();
    vector<string> computedValueNames(force.getNumComputedValues());
    vector<string> computedValueExpressions(force.getNumComputedValues());
    if (force.getNumComputedValues() > 0) {
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(0, computedValueNames[0], computedValueExpressions[0], type);
        if (type == CustomGBForce::SingleParticle)
            throw OpenMMException("CudaPlatform requires that the first computed value for a CustomGBForce be of type ParticlePair or ParticlePairNoExclusions.");
        useExclusionsForValue = (type == CustomGBForce::ParticlePair);
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            force.getComputedValueParameters(i, computedValueNames[i], computedValueExpressions[i], type);
            if (type != CustomGBForce::SingleParticle)
                throw OpenMMException("CudaPlatform requires that a CustomGBForce only have one computed value of type ParticlePair or ParticlePairNoExclusions.");
        }
    }
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "custom"+cu.intToString(forceIndex)+"_";

    // Record parameters and exclusions.

    int numParticles = force.getNumParticles();
    int paddedNumParticles = cu.getPaddedNumAtoms();
    int numParams = force.getNumPerParticleParameters();
    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), paddedNumParticles, "customGBParameters", true);
    computedValues = new CudaParameterSet(cu, force.getNumComputedValues(), paddedNumParticles, "customGBComputedValues", true, cu.getUseDoublePrecision());
    if (force.getNumGlobalParameters() > 0)
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customGBGlobals");
    vector<vector<float> > paramVector(paddedNumParticles, vector<float>(numParams, 0));
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
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
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions.push_back(CudaArray::create<float>(cu, f.size(), "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName, "float", width, width*sizeof(float), tabulatedFunctions[tabulatedFunctions.size()-1]->getDevicePointer()));
        tableArgs << ", const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* __restrict__ " << arrayName;
    }

    // Record the global parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals != NULL)
        globals->upload(globalParamValues);

    // Record derivatives of expressions needed for the chain rule terms.

    vector<vector<Lepton::ParsedExpression> > valueGradientExpressions(force.getNumComputedValues());
    vector<vector<Lepton::ParsedExpression> > valueDerivExpressions(force.getNumComputedValues());
    needParameterGradient = false;
    for (int i = 1; i < force.getNumComputedValues(); i++) {
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
        valueGradientExpressions[i].push_back(ex.differentiate("x").optimize());
        valueGradientExpressions[i].push_back(ex.differentiate("y").optimize());
        valueGradientExpressions[i].push_back(ex.differentiate("z").optimize());
        if (!isZeroExpression(valueGradientExpressions[i][0]) || !isZeroExpression(valueGradientExpressions[i][1]) || !isZeroExpression(valueGradientExpressions[i][2]))
            needParameterGradient = true;
         for (int j = 0; j < i; j++)
            valueDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
    }
    vector<vector<Lepton::ParsedExpression> > energyDerivExpressions(force.getNumEnergyTerms());
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
    }
    longEnergyDerivs = CudaArray::create<long long>(cu, force.getNumComputedValues()*cu.getPaddedNumAtoms(), "customGBLongEnergyDerivatives");
    energyDerivs = new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms(), "customGBEnergyDerivatives", true);
    energyDerivChain = new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms(), "customGBEnergyDerivativeChain", true);
 
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
            string value = "globals["+cu.intToString(i)+"]";
            variables.push_back(makeVariable(name, value));
        }
        map<string, Lepton::ParsedExpression> n2ValueExpressions;
        stringstream n2ValueSource;
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[0], functions).optimize();
        n2ValueExpressions["tempValue1 = "] = ex;
        n2ValueExpressions["tempValue2 = "] = ex.renameVariables(rename);
        n2ValueSource << cu.getExpressionUtilities().createExpressions(n2ValueExpressions, variables, functionList, functionDefinitions, "temp");
        map<string, string> replacements;
        string n2ValueStr = n2ValueSource.str();
        replacements["COMPUTE_VALUE"] = n2ValueStr;
        stringstream extraArgs, atomParams, loadLocal1, loadLocal2, load1, load2;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", const float* globals";
        pairValueUsesParam.resize(params->getBuffers().size(), false);
        int atomParamSize = 6;
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cu.intToString(i+1);
            if (n2ValueStr.find(paramName+"1") != n2ValueStr.npos || n2ValueStr.find(paramName+"2") != n2ValueStr.npos) {
                extraArgs << ", const " << buffer.getType() << "* __restrict__ global_" << paramName;
                atomParams << buffer.getType() << " " << paramName << ";\n";
                loadLocal1 << "localData[localAtomIndex]." << paramName << " = " << paramName << "1;\n";
                loadLocal2 << "localData[localAtomIndex]." << paramName << " = global_" << paramName << "[j];\n";
                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
                load2 << buffer.getType() << " " << paramName << "2 = localData[atom2]." << paramName << ";\n";
                pairValueUsesParam[i] = true;
                atomParamSize += buffer.getNumComponents();
            }
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        if (useCutoff)
            pairValueDefines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            pairValueDefines["USE_PERIODIC"] = "1";
        if (useExclusionsForValue)
            pairValueDefines["USE_EXCLUSIONS"] = "1";
        if (atomParamSize%2 == 0 && !cu.getUseDoublePrecision())
            pairValueDefines["NEED_PADDING"] = "1";
        pairValueDefines["WARPS_PER_GROUP"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize()/CudaContext::TileSize);
        pairValueDefines["THREAD_BLOCK_SIZE"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize());
        pairValueDefines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
        pairValueDefines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        pairValueDefines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        pairValueDefines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
        pairValueDefines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
        pairValueSrc = cu.replaceStrings(CudaKernelSources::customGBValueN2, replacements);
        if (useExclusionsForValue)
            cu.getNonbondedUtilities().requestExclusions(exclusionList);
    }
    {
        // Create the kernel to reduce the N2 value and calculate other values.

        stringstream reductionSource, extraArgs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", const float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cu.intToString(i+1);
            extraArgs << ", const " << buffer.getType() << "* __restrict__ " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cu.intToString(i+1);
            extraArgs << ", " << buffer.getType() << "* __restrict__ global_" << valueName;
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
            variables[force.getGlobalParameterName(i)] = "globals["+cu.intToString(i)+"]";
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            variables[computedValueNames[i-1]] = "local_values"+computedValues->getParameterSuffix(i-1);
            map<string, Lepton::ParsedExpression> valueExpressions;
            valueExpressions["local_values"+computedValues->getParameterSuffix(i)+" = "] = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
            reductionSource << cu.getExpressionUtilities().createExpressions(valueExpressions, variables, functionList, functionDefinitions, "value"+cu.intToString(i)+"_temp");
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            string valueName = "values"+cu.intToString(i+1);
            reductionSource << "global_" << valueName << "[index] = local_" << valueName << ";\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["COMPUTE_VALUES"] = reductionSource.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::customGBValuePerParticle, replacements), defines);
        perParticleValueKernel = cu.getKernel(module, "computePerParticleValues");
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
            variables.push_back(makeVariable(force.getGlobalParameterName(i), "globals["+cu.intToString(i)+"]"));
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
            for (int j = 0; j < force.getNumComputedValues(); j++) {
                if (needChainForValue[j]) {
                    string index = cu.intToString(j+1);
                    n2EnergyExpressions["/*"+cu.intToString(i+1)+"*/ deriv"+index+"_1 += "] = energyDerivExpressions[i][2*j];
                    n2EnergyExpressions["/*"+cu.intToString(i+1)+"*/ deriv"+index+"_2 += "] = energyDerivExpressions[i][2*j+1];
                }
            }
            if (exclude)
                n2EnergySource << "if (!isExcluded) {\n";
            n2EnergySource << cu.getExpressionUtilities().createExpressions(n2EnergyExpressions, variables, functionList, functionDefinitions, "temp");
            if (exclude)
                n2EnergySource << "}\n";
        }
        map<string, string> replacements;
        string n2EnergyStr = n2EnergySource.str();
        replacements["COMPUTE_INTERACTION"] = n2EnergyStr;
        stringstream extraArgs, atomParams, loadLocal1, loadLocal2, clearLocal, load1, load2, declare1, recordDeriv, storeDerivs1, storeDerivs2;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", const float* globals";
        pairEnergyUsesParam.resize(params->getBuffers().size(), false);
        int atomParamSize = 7;
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cu.intToString(i+1);
            if (n2EnergyStr.find(paramName+"1") != n2EnergyStr.npos || n2EnergyStr.find(paramName+"2") != n2EnergyStr.npos) {
                extraArgs << ", const " << buffer.getType() << "* __restrict__ global_" << paramName;
                atomParams << buffer.getType() << " " << paramName << ";\n";
                loadLocal1 << "localData[localAtomIndex]." << paramName << " = " << paramName << "1;\n";
                loadLocal2 << "localData[localAtomIndex]." << paramName << " = global_" << paramName << "[j];\n";
                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
                load2 << buffer.getType() << " " << paramName << "2 = localData[atom2]." << paramName << ";\n";
                pairEnergyUsesParam[i] = true;
                atomParamSize += buffer.getNumComponents();
            }
        }
        pairEnergyUsesValue.resize(computedValues->getBuffers().size(), false);
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cu.intToString(i+1);
            if (n2EnergyStr.find(valueName+"1") != n2EnergyStr.npos || n2EnergyStr.find(valueName+"2") != n2EnergyStr.npos) {
                extraArgs << ", const " << buffer.getType() << "* __restrict__ global_" << valueName;
                atomParams << buffer.getType() << " " << valueName << ";\n";
                loadLocal1 << "localData[localAtomIndex]." << valueName << " = " << valueName << "1;\n";
                loadLocal2 << "localData[localAtomIndex]." << valueName << " = global_" << valueName << "[j];\n";
                load1 << buffer.getType() << " " << valueName << "1 = global_" << valueName << "[atom1];\n";
                load2 << buffer.getType() << " " << valueName << "2 = localData[atom2]." << valueName << ";\n";
                pairEnergyUsesValue[i] = true;
                atomParamSize += buffer.getNumComponents();
            }
        }
        extraArgs << ", unsigned long long* __restrict__ derivBuffers";
        for (int i = 0; i < force.getNumComputedValues(); i++) {
            string index = cu.intToString(i+1);
            atomParams << "real deriv" << index << ";\n";
            clearLocal << "localData[localAtomIndex].deriv" << index << " = 0;\n";
            declare1 << "real deriv" << index << "_1 = 0;\n";
            load2 << "real deriv" << index << "_2 = 0;\n";
            recordDeriv << "localData[atom2].deriv" << index << " += deriv" << index << "_2;\n";
            storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")\n";
            storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")\n";
            atomParamSize++;
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
        replacements["CLEAR_LOCAL_DERIVATIVES"] = clearLocal.str();
        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
        replacements["DECLARE_ATOM1_DERIVATIVES"] = declare1.str();
        replacements["RECORD_DERIVATIVE_2"] = recordDeriv.str();
        replacements["STORE_DERIVATIVES_1"] = storeDerivs1.str();
        replacements["STORE_DERIVATIVES_2"] = storeDerivs2.str();
        if (useCutoff)
            pairEnergyDefines["USE_CUTOFF"] = "1";
        if (usePeriodic)
            pairEnergyDefines["USE_PERIODIC"] = "1";
        if (anyExclusions)
            pairEnergyDefines["USE_EXCLUSIONS"] = "1";
        if (atomParamSize%2 != 0 && !cu.getUseDoublePrecision())
            pairEnergyDefines["NEED_PADDING"] = "1";
        pairEnergyDefines["THREAD_BLOCK_SIZE"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize());
        pairEnergyDefines["WARPS_PER_GROUP"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize()/CudaContext::TileSize);
        pairEnergyDefines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
        pairEnergyDefines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        pairEnergyDefines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        pairEnergyDefines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
        pairEnergyDefines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
        pairEnergySrc = cu.replaceStrings(CudaKernelSources::customGBEnergyN2, replacements);
    }
    {
        // Create the kernel to reduce the derivatives and calculate per-particle energy terms.

        stringstream compute, extraArgs, load;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", const float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cu.intToString(i+1);
            extraArgs << ", const " << buffer.getType() << "* __restrict__ " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cu.intToString(i+1);
            extraArgs << ", const " << buffer.getType() << "* __restrict__ " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string index = cu.intToString(i+1);
            extraArgs << ", " << buffer.getType() << "* __restrict__ derivBuffers" << index;
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        for (int i = 0; i < (int) energyDerivChain->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivChain->getBuffers()[i];
            string index = cu.intToString(i+1);
            extraArgs << ", " << buffer.getType() << "* __restrict__ derivChain" << index;
        }
        extraArgs << ", const long long* __restrict__ derivBuffersIn";
        for (int i = 0; i < energyDerivs->getNumParameters(); ++i)
            load << "derivBuffers" << energyDerivs->getParameterSuffix(i, "[index]") <<
                    " = RECIP(0x100000000)*derivBuffersIn[index+PADDED_NUM_ATOMS*" << cu.intToString(i) << "];\n";
        
        // Compute the various expressions.
        
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cu.intToString(i)+"]";
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
            expressions["/*"+cu.intToString(i+1)+"*/ energy += "] = parsed;
            for (int j = 0; j < force.getNumComputedValues(); j++)
                expressions["/*"+cu.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j)+" += "] = energyDerivExpressions[i][j];
            Lepton::ParsedExpression gradx = parsed.differentiate("x").optimize();
            Lepton::ParsedExpression grady = parsed.differentiate("y").optimize();
            Lepton::ParsedExpression gradz = parsed.differentiate("z").optimize();
            if (!isZeroExpression(gradx))
                expressions["/*"+cu.intToString(i+1)+"*/ force.x -= "] = gradx;
            if (!isZeroExpression(grady))
                expressions["/*"+cu.intToString(i+1)+"*/ force.y -= "] = grady;
            if (!isZeroExpression(gradz))
                expressions["/*"+cu.intToString(i+1)+"*/ force.z -= "] = gradz;
        }
        for (int i = 1; i < force.getNumComputedValues(); i++)
            for (int j = 0; j < i; j++)
                expressions["real dV"+cu.intToString(i)+"dV"+cu.intToString(j)+" = "] = valueDerivExpressions[i][j];
        compute << cu.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "temp");
        
        // Record values.
        
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            string index = cu.intToString(i+1);
            compute << "derivBuffers" << index << "[index] = deriv" << index << ";\n";
        }
        compute << "forceBuffers[index] += (long long) (force.x*0x100000000);\n";
        compute << "forceBuffers[index+PADDED_NUM_ATOMS] += (long long) (force.y*0x100000000);\n";
        compute << "forceBuffers[index+PADDED_NUM_ATOMS*2] += (long long) (force.z*0x100000000);\n";
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            compute << "real totalDeriv"<<i<<" = dV"<<i<<"dV0";
            for (int j = 1; j < i; j++)
                compute << " + totalDeriv"<<j<<"*dV"<<i<<"dV"<<j;
            compute << ";\n";
            compute << "deriv"<<(i+1)<<" *= totalDeriv"<<i<<";\n";
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            string index = cu.intToString(i+1);
            compute << "derivChain" << index << "[index] = deriv" << index << ";\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["LOAD_DERIVATIVES"] = load.str();
        replacements["COMPUTE_ENERGY"] = compute.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::customGBEnergyPerParticle, replacements), defines);
        perParticleEnergyKernel = cu.getKernel(module, "computePerParticleEnergy");
    }
    if (needParameterGradient) {
        // Create the kernel to compute chain rule terms for computed values that depend explicitly on particle coordinates.

        stringstream compute, extraArgs;
        if (force.getNumGlobalParameters() > 0)
            extraArgs << ", const float* globals";
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = "params"+cu.intToString(i+1);
            extraArgs << ", const " << buffer.getType() << "* __restrict__ " << paramName;
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string valueName = "values"+cu.intToString(i+1);
            extraArgs << ", const " << buffer.getType() << "* __restrict__ " << valueName;
        }
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
            string index = cu.intToString(i+1);
            extraArgs << ", " << buffer.getType() << "* __restrict__ derivBuffers" << index;
            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
        }
        map<string, string> variables;
        variables["x"] = "pos.x";
        variables["y"] = "pos.y";
        variables["z"] = "pos.z";
        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
        for (int i = 0; i < force.getNumGlobalParameters(); i++)
            variables[force.getGlobalParameterName(i)] = "globals["+cu.intToString(i)+"]";
        for (int i = 0; i < force.getNumComputedValues(); i++)
            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            string is = cu.intToString(i);
            compute << "real3 dV"<<is<<"dR = make_real3(0);\n";
            for (int j = 1; j < i; j++) {
                if (!isZeroExpression(valueDerivExpressions[i][j])) {
                    map<string, Lepton::ParsedExpression> derivExpressions;
                    string js = cu.intToString(j);
                    derivExpressions["real dV"+is+"dV"+js+" = "] = valueDerivExpressions[i][j];
                    compute << cu.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, "temp_"+is+"_"+js);
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
            compute << cu.getExpressionUtilities().createExpressions(gradientExpressions, variables, functionList, functionDefinitions, "temp");
        }
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            string is = cu.intToString(i);
            compute << "force -= deriv"<<energyDerivs->getParameterSuffix(i)<<"*dV"<<is<<"dR;\n";
        }
        map<string, string> replacements;
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["COMPUTE_FORCES"] = compute.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        CUmodule module = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::customGBGradientChainRule, replacements), defines);
        gradientChainRuleKernel = cu.getKernel(module, "computeGradientChainRuleTerms");
    }
    {
        // Create the code to calculate chain rules terms as part of the default nonbonded kernel.

        vector<pair<ExpressionTreeNode, string> > globalVariables;
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cu.intToString(i)+"]";
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
        chainSource << cu.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, prefix+"temp0_");
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
        string source = cu.replaceStrings(CudaKernelSources::customGBChainRule, replacements);
        vector<CudaNonbondedUtilities::ParameterInfo> parameters;
        vector<CudaNonbondedUtilities::ParameterInfo> arguments;
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            string paramName = prefix+"params"+cu.intToString(i+1);
            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
                parameters.push_back(CudaNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
            string paramName = prefix+"values"+cu.intToString(i+1);
            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
                parameters.push_back(CudaNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        for (int i = 0; i < (int) energyDerivChain->getBuffers().size(); i++) {
            if (needChainForValue[i]) { 
                CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivChain->getBuffers()[i];
                string paramName = prefix+"dEdV"+cu.intToString(i+1);
                parameters.push_back(CudaNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
            }
        }
        if (globals != NULL) {
            globals->upload(globalParamValues);
            arguments.push_back(CudaNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(float), globals->getDevicePointer()));
        }
        cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, force.getNumExclusions() > 0, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
        for (int i = 0; i < (int) parameters.size(); i++)
            cu.getNonbondedUtilities().addParameter(parameters[i]);
        for (int i = 0; i < (int) arguments.size(); i++)
            cu.getNonbondedUtilities().addArgument(arguments[i]);
    }
    cu.addForce(new CudaCustomGBForceInfo(force));
    cu.addAutoclearBuffer(*longEnergyDerivs);
}

double CudaCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // These two kernels can't be compiled in initialize(), because the nonbonded utilities object
        // has not yet been initialized then.

        {
            int numExclusionTiles = cu.getNonbondedUtilities().getExclusionTiles().getSize();
            pairValueDefines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(numExclusionTiles);
            int numContexts = cu.getPlatformData().contexts.size();
            int startExclusionIndex = cu.getContextIndex()*numExclusionTiles/numContexts;
            int endExclusionIndex = (cu.getContextIndex()+1)*numExclusionTiles/numContexts;
            pairValueDefines["FIRST_EXCLUSION_TILE"] = cu.intToString(startExclusionIndex);
            pairValueDefines["LAST_EXCLUSION_TILE"] = cu.intToString(endExclusionIndex);
            pairValueDefines["CUTOFF"] = cu.doubleToString(nb.getCutoffDistance());
            CUmodule module = cu.createModule(CudaKernelSources::vectorOps+pairValueSrc, pairValueDefines);
            pairValueKernel = cu.getKernel(module, "computeN2Value");
            pairValueSrc = "";
            pairValueDefines.clear();
        }
        {
            int numExclusionTiles = cu.getNonbondedUtilities().getExclusionTiles().getSize();
            pairEnergyDefines["NUM_TILES_WITH_EXCLUSIONS"] = cu.intToString(numExclusionTiles);
            int numContexts = cu.getPlatformData().contexts.size();
            int startExclusionIndex = cu.getContextIndex()*numExclusionTiles/numContexts;
            int endExclusionIndex = (cu.getContextIndex()+1)*numExclusionTiles/numContexts;
            pairEnergyDefines["FIRST_EXCLUSION_TILE"] = cu.intToString(startExclusionIndex);
            pairEnergyDefines["LAST_EXCLUSION_TILE"] = cu.intToString(endExclusionIndex);
            pairEnergyDefines["CUTOFF"] = cu.doubleToString(nb.getCutoffDistance());
            CUmodule module = cu.createModule(CudaKernelSources::vectorOps+pairEnergySrc, pairEnergyDefines);
            pairEnergyKernel = cu.getKernel(module, "computeN2Energy");
            pairEnergySrc = "";
            pairEnergyDefines.clear();
        }

        // Set arguments for kernels.
        
        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
        valueBuffers = CudaArray::create<long long>(cu, cu.getPaddedNumAtoms(), "customGBValueBuffers");
        cu.addAutoclearBuffer(*valueBuffers);
        cu.clearBuffer(valueBuffers->getDevicePointer(), sizeof(long long)*valueBuffers->getSize());
        pairValueArgs.push_back(&cu.getPosq().getDevicePointer());
        pairValueArgs.push_back(&cu.getNonbondedUtilities().getExclusions().getDevicePointer());
        pairValueArgs.push_back(&cu.getNonbondedUtilities().getExclusionTiles().getDevicePointer());
        pairValueArgs.push_back(&valueBuffers->getDevicePointer());
        if (nb.getUseCutoff()) {
            pairValueArgs.push_back(&nb.getInteractingTiles().getDevicePointer());
            pairValueArgs.push_back(&nb.getInteractionCount().getDevicePointer());
            pairValueArgs.push_back(cu.getPeriodicBoxSizePointer());
            pairValueArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            pairValueArgs.push_back(&maxTiles);
            pairValueArgs.push_back(&nb.getBlockCenters().getDevicePointer());
            pairValueArgs.push_back(&nb.getBlockBoundingBoxes().getDevicePointer());
            pairValueArgs.push_back(&nb.getInteractingAtoms().getDevicePointer());
        }
        else
            pairValueArgs.push_back(&maxTiles);
        if (globals != NULL)
            pairValueArgs.push_back(&globals->getDevicePointer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            if (pairValueUsesParam[i])
                pairValueArgs.push_back(&params->getBuffers()[i].getMemory());
        }
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            pairValueArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
        perParticleValueArgs.push_back(&cu.getPosq().getDevicePointer());
        perParticleValueArgs.push_back(&valueBuffers->getDevicePointer());
        if (globals != NULL)
            perParticleValueArgs.push_back(&globals->getDevicePointer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++)
            perParticleValueArgs.push_back(&params->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
            perParticleValueArgs.push_back(&computedValues->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            perParticleValueArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
        pairEnergyArgs.push_back(&cu.getForce().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getPosq().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getNonbondedUtilities().getExclusions().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getNonbondedUtilities().getExclusionTiles().getDevicePointer());
        if (nb.getUseCutoff()) {
            pairEnergyArgs.push_back(&nb.getInteractingTiles().getDevicePointer());
            pairEnergyArgs.push_back(&nb.getInteractionCount().getDevicePointer());
            pairEnergyArgs.push_back(cu.getPeriodicBoxSizePointer());
            pairEnergyArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            pairEnergyArgs.push_back(&maxTiles);
            pairEnergyArgs.push_back(&nb.getBlockCenters().getDevicePointer());
            pairEnergyArgs.push_back(&nb.getBlockBoundingBoxes().getDevicePointer());
            pairEnergyArgs.push_back(&nb.getInteractingAtoms().getDevicePointer());
        }
        else
            pairEnergyArgs.push_back(&maxTiles);
        if (globals != NULL)
            pairEnergyArgs.push_back(&globals->getDevicePointer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            if (pairEnergyUsesParam[i])
                pairEnergyArgs.push_back(&params->getBuffers()[i].getMemory());
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            if (pairEnergyUsesValue[i])
                pairEnergyArgs.push_back(&computedValues->getBuffers()[i].getMemory());
        }
        pairEnergyArgs.push_back(&longEnergyDerivs->getDevicePointer());
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            pairEnergyArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
        perParticleEnergyArgs.push_back(&cu.getForce().getDevicePointer());
        perParticleEnergyArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        perParticleEnergyArgs.push_back(&cu.getPosq().getDevicePointer());
        if (globals != NULL)
            perParticleEnergyArgs.push_back(&globals->getDevicePointer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++)
            perParticleEnergyArgs.push_back(&params->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
            perParticleEnergyArgs.push_back(&computedValues->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
            perParticleEnergyArgs.push_back(&energyDerivs->getBuffers()[i].getMemory());
        for (int i = 0; i < (int) energyDerivChain->getBuffers().size(); i++)
            perParticleEnergyArgs.push_back(&energyDerivChain->getBuffers()[i].getMemory());
        perParticleEnergyArgs.push_back(&longEnergyDerivs->getDevicePointer());
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            perParticleEnergyArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
        if (needParameterGradient) {
            gradientChainRuleArgs.push_back(&cu.getForce().getDevicePointer());
            gradientChainRuleArgs.push_back(&cu.getPosq().getDevicePointer());
            if (globals != NULL)
                gradientChainRuleArgs.push_back(&globals->getDevicePointer());
            for (int i = 0; i < (int) params->getBuffers().size(); i++)
                gradientChainRuleArgs.push_back(&params->getBuffers()[i].getMemory());
            for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
                gradientChainRuleArgs.push_back(&computedValues->getBuffers()[i].getMemory());
            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
                gradientChainRuleArgs.push_back(&energyDerivs->getBuffers()[i].getMemory());
        }
    }
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    if (nb.getUseCutoff()) {
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            pairValueArgs[4] = &nb.getInteractingTiles().getDevicePointer();
            pairEnergyArgs[5] = &nb.getInteractingTiles().getDevicePointer();
            pairValueArgs[11] = &nb.getInteractingAtoms().getDevicePointer();
            pairEnergyArgs[12] = &nb.getInteractingAtoms().getDevicePointer();
        }
    }
    cu.executeKernel(pairValueKernel, &pairValueArgs[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cu.executeKernel(perParticleValueKernel, &perParticleValueArgs[0], cu.getPaddedNumAtoms());
    cu.executeKernel(pairEnergyKernel, &pairEnergyArgs[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cu.executeKernel(perParticleEnergyKernel, &perParticleEnergyArgs[0], cu.getPaddedNumAtoms());
    if (needParameterGradient)
        cu.executeKernel(gradientChainRuleKernel, &gradientChainRuleArgs[0], cu.getPaddedNumAtoms());
    return 0.0;
}

void CudaCalcCustomGBForceKernel::copyParametersToContext(ContextImpl& context, const CustomGBForce& force) {
    cu.setAsCurrent();
    int numParticles = force.getNumParticles();
    if (numParticles != cu.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<vector<float> > paramVector(cu.getPaddedNumAtoms(), vector<float>(force.getNumPerParticleParameters(), 0));
    vector<double> parameters;
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCustomExternalForceInfo : public CudaForceInfo {
public:
    CudaCustomExternalForceInfo(const CustomExternalForce& force, int numParticles) : force(force), indices(numParticles, -1) {
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

CudaCalcCustomExternalForceKernel::~CudaCalcCustomExternalForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
}

void CudaCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumParticles()/numContexts;
    numParticles = endIndex-startIndex;
    if (numParticles == 0)
        return;
    vector<vector<int> > atoms(numParticles, vector<int>(1));
    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customExternalParams");
    vector<vector<float> > paramVector(numParticles);
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        force.getParticleParameters(startIndex+i, atoms[i][0], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    cu.addForce(new CudaCustomExternalForceInfo(force, system.getNumParticles()));

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
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
    variables["x"] = "pos1.x";
    variables["y"] = "pos1.y";
    variables["z"] = "pos1.z";
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        variables[name] = "particleParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customExternalGlobals");
        globals->upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals->getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    stringstream compute;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" particleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    compute << cu.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
    map<string, string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::customExternalForce, replacements), force.getForceGroup());
}

double CudaCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    return 0.0;
}

void CudaCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumParticles()/numContexts;
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
    
    cu.invalidateMolecules();
}

class CudaCustomHbondForceInfo : public CudaForceInfo {
public:
    CudaCustomHbondForceInfo(const CustomHbondForce& force) : force(force) {
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

CudaCalcCustomHbondForceKernel::~CudaCalcCustomHbondForceKernel() {
    cu.setAsCurrent();
    if (donorParams != NULL)
        delete donorParams;
    if (acceptorParams != NULL)
        delete acceptorParams;
    if (donors != NULL)
        delete donors;
    if (acceptors != NULL)
        delete acceptors;
    if (globals != NULL)
        delete globals;
    if (donorExclusions != NULL)
        delete donorExclusions;
    if (acceptorExclusions != NULL)
        delete acceptorExclusions;
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
        applyToAcceptor << forceNames[atom]<<" += trim("<<value<<");\n";
    else
        applyToDonor << forceNames[atom-3]<<" += trim("<<value<<");\n";
}

void CudaCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
    // Record the lists of donors and acceptors, and the parameters for each one.

    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumDonors()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumDonors()/numContexts;
    numDonors = endIndex-startIndex;
    numAcceptors = force.getNumAcceptors();
    if (numDonors == 0 || numAcceptors == 0)
        return;
    int numParticles = system.getNumParticles();
    donors = CudaArray::create<int4>(cu, numDonors, "customHbondDonors");
    acceptors = CudaArray::create<int4>(cu, numAcceptors, "customHbondAcceptors");
    donorParams = new CudaParameterSet(cu, force.getNumPerDonorParameters(), numDonors, "customHbondDonorParameters");
    acceptorParams = new CudaParameterSet(cu, force.getNumPerAcceptorParameters(), numAcceptors, "customHbondAcceptorParameters");
    if (force.getNumGlobalParameters() > 0)
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customHbondGlobals");
    vector<vector<float> > donorParamVector(numDonors);
    vector<int4> donorVector(numDonors);
    for (int i = 0; i < numDonors; i++) {
        vector<double> parameters;
        force.getDonorParameters(startIndex+i, donorVector[i].x, donorVector[i].y, donorVector[i].z, parameters);
        donorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            donorParamVector[i][j] = (float) parameters[j];
    }
    donors->upload(donorVector);
    donorParams->setParameterValues(donorParamVector);
    vector<vector<float> > acceptorParamVector(numAcceptors);
    vector<int4> acceptorVector(numAcceptors);
    for (int i = 0; i < numAcceptors; i++) {
        vector<double> parameters;
        force.getAcceptorParameters(i, acceptorVector[i].x, acceptorVector[i].y, acceptorVector[i].z, parameters);
        acceptorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            acceptorParamVector[i][j] = (float) parameters[j];
    }
    acceptors->upload(acceptorVector);
    acceptorParams->setParameterValues(acceptorParamVector);
    cu.addForce(new CudaCustomHbondForceInfo(force));

    // Record exclusions.

    vector<int4> donorExclusionVector(numDonors, make_int4(-1, -1, -1, -1));
    vector<int4> acceptorExclusionVector(numAcceptors, make_int4(-1, -1, -1, -1));
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
            throw OpenMMException("CustomHbondForce: CudaPlatform does not support more than four exclusions per donor");
        if (acceptorExclusionVector[acceptor].x == -1)
            acceptorExclusionVector[acceptor].x = donor;
        else if (acceptorExclusionVector[acceptor].y == -1)
            acceptorExclusionVector[acceptor].y = donor;
        else if (acceptorExclusionVector[acceptor].z == -1)
            acceptorExclusionVector[acceptor].z = donor;
        else if (acceptorExclusionVector[acceptor].w == -1)
            acceptorExclusionVector[acceptor].w = donor;
        else
            throw OpenMMException("CustomHbondForce: CudaPlatform does not support more than four exclusions per acceptor");
    }
    donorExclusions = CudaArray::create<int4>(cu, numDonors, "customHbondDonorExclusions");
    acceptorExclusions = CudaArray::create<int4>(cu, numAcceptors, "customHbondAcceptorExclusions");
    donorExclusions->upload(donorExclusionVector);
    acceptorExclusions->upload(acceptorExclusionVector);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions.push_back(CudaArray::create<float>(cu, f.size(), "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        tableArgs << ", const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* __restrict__ " << arrayName;
    }

    // Record information about parameters.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
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
        variables[name] = "globals["+cu.intToString(i)+"]";
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
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+");\n");
            computedDeltas.insert(deltaName);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real r_"+deltaName+" = SQRT(delta"+deltaName+".w);\n");
        variables[iter->first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[0]]+");\n");
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[2]]+");\n");
            computedDeltas.insert(deltaName2);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real "+angleName+" = computeAngle(delta"+deltaName1+", delta"+deltaName2+");\n");
        variables[iter->first] = angleName;
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
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
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+");\n");
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[1]]+");\n");
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName3+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[3]]+");\n");
            computedDeltas.insert(deltaName3);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 "+crossName1+" = computeCross(delta"+deltaName1+", delta"+deltaName2+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 "+crossName2+" = computeCross(delta"+deltaName2+", delta"+deltaName3+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real "+dihedralName+" = computeAngle("+crossName1+", "+crossName2+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, dihedralName+" *= (delta"+deltaName1+".x*"+crossName2+".x + delta"+deltaName1+".y*"+crossName2+".y + delta"+deltaName1+".z*"+crossName2+".z < 0 ? -1 : 1);\n");
        variables[iter->first] = dihedralName;
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }

    // Next it needs to load parameters from global memory.

    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", const float* __restrict__ globals";
    for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
        extraArgs << ", const "+buffer.getType()+"* __restrict__ donor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" donorParams"+cu.intToString(i+1)+" = donor"+buffer.getName()+"[index];\n");
    }
    for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
        extraArgs << ", const "+buffer.getType()+"* __restrict__ acceptor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" acceptorParams"+cu.intToString(i+1)+" = acceptor"+buffer.getName()+"[index];\n");
    }

    // Now evaluate the expressions.

    computeAcceptor << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");
    forceExpressions["energy += "] = energyExpression;
    computeDonor << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*delta"+deltaName;
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "-"+value);
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], value);
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 crossProd = cross(delta"+deltaName2+", delta"+deltaName1+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real lengthCross = max(SQRT(dot(crossProd,crossProd)), 1e-6f);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 deltaCross0 = -cross(trim(delta"+deltaName1+"), crossProd)*dEdAngle"+cu.intToString(index)+"/(delta"+deltaName1+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 deltaCross2 = cross(trim(delta"+deltaName2+"), crossProd)*dEdAngle"+cu.intToString(index)+"/(delta"+deltaName2+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 deltaCross1 = -(deltaCross0+deltaCross2);\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "deltaCross0");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "deltaCross1");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "deltaCross2");
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
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real r = SQRT(delta"+deltaName2+".w);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 ff;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.x = (-dEdDihedral"+cu.intToString(index)+"*r)/"+crossName1+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.y = (delta"+deltaName1+".x*delta"+deltaName2+".x + delta"+deltaName1+".y*delta"+deltaName2+".y + delta"+deltaName1+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.z = (delta"+deltaName3+".x*delta"+deltaName2+".x + delta"+deltaName3+".y*delta"+deltaName2+".y + delta"+deltaName3+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.w = (dEdDihedral"+cu.intToString(index)+"*r)/"+crossName2+".w;\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 internalF0 = ff.x*"+crossName1+";\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 internalF3 = ff.w*"+crossName2+";\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 s = ff.y*internalF0 - ff.z*internalF3;\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "internalF0");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "s-internalF0");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "-s-internalF3");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[3], "internalF3");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
    }

    // Generate the kernels.

    map<string, string> replacements;
    replacements["COMPUTE_DONOR_FORCE"] = computeDonor.str();
    replacements["COMPUTE_ACCEPTOR_FORCE"] = computeAcceptor.str();
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string, string> defines;
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["NUM_DONORS"] = cu.intToString(numDonors);
    defines["NUM_ACCEPTORS"] = cu.intToString(numAcceptors);
    defines["M_PI"] = cu.doubleToString(M_PI);
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff) {
        defines["USE_CUTOFF"] = "1";
        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff && force.getNonbondedMethod() != CustomHbondForce::CutoffNonPeriodic)
        defines["USE_PERIODIC"] = "1";
    if (force.getNumExclusions() > 0)
        defines["USE_EXCLUSIONS"] = "1";
    CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaKernelSources::customHbondForce, replacements), defines);
    donorKernel = cu.getKernel(module, "computeDonorForces");
    acceptorKernel = cu.getKernel(module, "computeAcceptorForces");
}

double CudaCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (numDonors == 0 || numAcceptors == 0)
        return 0.0;
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
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
        donorArgs.push_back(&cu.getForce().getDevicePointer());
        donorArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        donorArgs.push_back(&cu.getPosq().getDevicePointer());
        donorArgs.push_back(&donorExclusions->getDevicePointer());
        donorArgs.push_back(&donors->getDevicePointer());
        donorArgs.push_back(&acceptors->getDevicePointer());
        donorArgs.push_back(cu.getPeriodicBoxSizePointer());
        donorArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        if (globals != NULL)
            donorArgs.push_back(&globals->getDevicePointer());
        for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
            donorArgs.push_back(&buffer.getMemory());
        }
        for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
            donorArgs.push_back(&buffer.getMemory());
        }
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            donorArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
        index = 0;
        acceptorArgs.push_back(&cu.getForce().getDevicePointer());
        acceptorArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        acceptorArgs.push_back(&cu.getPosq().getDevicePointer());
        acceptorArgs.push_back(&acceptorExclusions->getDevicePointer());
        acceptorArgs.push_back(&donors->getDevicePointer());
        acceptorArgs.push_back(&acceptors->getDevicePointer());
        acceptorArgs.push_back(cu.getPeriodicBoxSizePointer());
        acceptorArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        if (globals != NULL)
            acceptorArgs.push_back(&globals->getDevicePointer());
        for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
            acceptorArgs.push_back(&buffer.getMemory());
        }
        for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
            acceptorArgs.push_back(&buffer.getMemory());
        }
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            acceptorArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
    }
    int sharedMemorySize = 3*CudaContext::ThreadBlockSize*sizeof(float4);
    cu.executeKernel(donorKernel, &donorArgs[0], max(numDonors, numAcceptors), CudaContext::ThreadBlockSize, sharedMemorySize);
    cu.executeKernel(acceptorKernel, &acceptorArgs[0], max(numDonors, numAcceptors), CudaContext::ThreadBlockSize, sharedMemorySize);
    return 0.0;
}

void CudaCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumDonors()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumDonors()/numContexts;
    if (numDonors != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of donors has changed");
    if (numAcceptors != force.getNumAcceptors())
        throw OpenMMException("updateParametersInContext: The number of acceptors has changed");
    
    // Record the per-donor parameters.
    
    if (numDonors > 0) {
        vector<vector<float> > donorParamVector(numDonors);
        vector<double> parameters;
        for (int i = 0; i < numDonors; i++) {
            int d1, d2, d3;
            force.getDonorParameters(startIndex+i, d1, d2, d3, parameters);
            donorParamVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                donorParamVector[i][j] = (float) parameters[j];
        }
        donorParams->setParameterValues(donorParamVector);
    }
    
    // Record the per-acceptor parameters.
    
    if (numAcceptors > 0) {
        vector<vector<float> > acceptorParamVector(numAcceptors);
        vector<double> parameters;
        for (int i = 0; i < numAcceptors; i++) {
            int a1, a2, a3;
            force.getAcceptorParameters(i, a1, a2, a3, parameters);
            acceptorParamVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++)
                acceptorParamVector[i][j] = (float) parameters[j];
        }
        acceptorParams->setParameterValues(acceptorParamVector);
    }
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCustomCompoundBondForceInfo : public CudaForceInfo {
public:
    CudaCustomCompoundBondForceInfo(const CustomCompoundBondForce& force) : force(force) {
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

CudaCalcCustomCompoundBondForceKernel::~CudaCalcCustomCompoundBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (globals != NULL)
        delete globals;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
}

void CudaCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    int particlesPerBond = force.getNumParticlesPerBond();
    vector<vector<int> > atoms(numBonds, vector<int>(particlesPerBond));
    params = new CudaParameterSet(cu, force.getNumPerBondParameters(), numBonds, "customCompoundBondParams");
    vector<vector<float> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<double> parameters;
        force.getBondParameters(startIndex+i, atoms[i], parameters);
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            paramVector[i][j] = (float) parameters[j];
    }
    params->setParameterValues(paramVector);
    cu.addForce(new CudaCustomCompoundBondForceInfo(force));

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        CudaArray* array = CudaArray::create<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions.push_back(array);
        array->upload(f);
        string arrayName = cu.getBondedUtilities().addArgument(array->getDevicePointer(), width == 1 ? "float" : "float"+cu.intToString(width));
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
        string index = cu.intToString(i+1);
        variables["x"+index] = "pos"+index+".x";
        variables["y"+index] = "pos"+index+".y";
        variables["z"+index] = "pos"+index+".z";
    }
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    if (force.getNumGlobalParameters() > 0) {
        globals = CudaArray::create<float>(cu, force.getNumGlobalParameters(), "customCompoundBondGlobals");
        globals->upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals->getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
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
        string index = cu.intToString(i+1);
        atomNames.push_back("P"+index);
        posNames.push_back("pos"+index);
    }
    stringstream compute;
    int index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<");\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
        variables[iter->first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<");\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<");\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = ccb_computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        variables[iter->first] = angleName;
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
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
            compute<<"real4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<");\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[1]]<<");\n";
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            compute<<"real4 delta"<<deltaName3<<" = ccb_delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[3]]<<");\n";
            computedDeltas.insert(deltaName3);
        }
        compute<<"real4 "<<crossName1<<" = ccb_computeCross(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        compute<<"real4 "<<crossName2<<" = ccb_computeCross(delta"<<deltaName2<<", delta"<<deltaName3<<");\n";
        compute<<"real "<<dihedralName<<" = ccb_computeAngle("<<crossName1<<", "<<crossName2<<");\n";
        compute<<dihedralName<<" *= (delta"<<deltaName1<<".x*"<<crossName2<<".x + delta"<<deltaName1<<".y*"<<crossName2<<".y + delta"<<deltaName1<<".z*"<<crossName2<<".z < 0 ? -1 : 1);\n";
        variables[iter->first] = dihedralName;
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    compute << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    vector<string> forceNames;
    for (int i = 0; i < particlesPerBond; i++) {
        string istr = cu.intToString(i+1);
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
            compute<<cu.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*ccb_trim(delta"+deltaName+")";
        compute<<forceNames[atoms[0]]<<" += "<<"-"<<value<<";\n";
        compute<<forceNames[atoms[1]]<<" += "<<value<<";\n";
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute<<"{\n";
        compute<<"real3 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(ccb_trim(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(ccb_trim(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real3 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[atoms[0]]<<" += deltaCross0;\n";
        compute<<forceNames[atoms[1]]<<" += deltaCross1;\n";
        compute<<forceNames[atoms[2]]<<" += deltaCross2;\n";
        compute<<"}\n";
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        compute<<"{\n";
        compute<<"real r = sqrt(delta"<<deltaName2<<".w);\n";
        compute<<"real4 ff;\n";
        compute<<"ff.x = (-dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real3 internalF0 = ff.x*ccb_trim("<<crossName1<<");\n";
        compute<<"real3 internalF3 = ff.w*ccb_trim("<<crossName2<<");\n";
        compute<<"real3 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[atoms[0]]<<" += internalF0;\n";
        compute<<forceNames[atoms[1]]<<" += s-internalF0;\n";
        compute<<forceNames[atoms[2]]<<" += -s-internalF3;\n";
        compute<<forceNames[atoms[3]]<<" += internalF3;\n";
        compute<<"}\n";
    }
    cu.getBondedUtilities().addInteraction(atoms, compute.str(), force.getForceGroup());
    map<string, string> replacements;
    replacements["M_PI"] = cu.doubleToString(M_PI);
    cu.getBondedUtilities().addPrefixCode(cu.replaceStrings(CudaKernelSources::customCompoundBond, replacements));;
}

double CudaCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (globals != NULL) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            globals->upload(globalParamValues);
    }
    return 0.0;
}

void CudaCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
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
    
    cu.invalidateMolecules();
}

class CudaCustomManyParticleForceInfo : public CudaForceInfo {
public:
    CudaCustomManyParticleForceInfo(const CustomManyParticleForce& force) : force(force) {
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

CudaCalcCustomManyParticleForceKernel::~CudaCalcCustomManyParticleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
    if (orderIndex != NULL)
        delete orderIndex;
    if (particleOrder != NULL)
        delete particleOrder;
    if (particleTypes != NULL)
        delete particleTypes;
    if (exclusions != NULL)
        delete exclusions;
    if (exclusionStartIndex != NULL)
        delete exclusionStartIndex;
    if (blockCenter != NULL)
        delete blockCenter;
    if (blockBoundingBox != NULL)
        delete blockBoundingBox;
    if (neighborPairs != NULL)
        delete neighborPairs;
    if (numNeighborPairs != NULL)
        delete numNeighborPairs;
    if (neighborStartIndex != NULL)
        delete neighborStartIndex;
    if (neighbors != NULL)
        delete neighbors;
    if (numNeighborsForAtom != NULL)
        delete numNeighborsForAtom;
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        delete tabulatedFunctions[i];
}

void CudaCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {
    cu.setAsCurrent();
    int numParticles = force.getNumParticles();
    int particlesPerSet = force.getNumParticlesPerSet();
    bool centralParticleMode = (force.getPermutationMode() == CustomManyParticleForce::UniqueCentralParticle);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    forceWorkgroupSize = 128;
    findNeighborsWorkgroupSize = 128;
    
    // Record parameter values.
    
    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customManyParticleParameters");
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
    cu.addForce(new CudaCustomManyParticleForceInfo(force));

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions.push_back(CudaArray::create<float>(cu, f.size(), "TabulatedFunction"));
        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
        tableArgs << ", const float";
        if (width > 1)
            tableArgs << width;
        tableArgs << "* __restrict__ " << arrayName;
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
        string index = cu.intToString(i+1);
        variables.push_back(makeVariable("x"+index, "pos"+index+".x"));
        variables.push_back(makeVariable("y"+index, "pos"+index+".y"));
        variables.push_back(makeVariable("z"+index, "pos"+index+".z"));
    }
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        for (int j = 0; j < particlesPerSet; j++) {
            string index = cu.intToString(j+1);
            variables.push_back(makeVariable(name+index, "params"+params->getParameterSuffix(i, index)));
        }
    }
    if (force.getNumGlobalParameters() > 0) {
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cu.intToString(i)+"]";
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
        particleTypes = CudaArray::create<int>(cu, particleTypesVec.size(), "customManyParticleTypes");
        orderIndex = CudaArray::create<int>(cu, orderIndexVec.size(), "customManyParticleOrderIndex");
        particleOrder = CudaArray::create<int>(cu, particleOrderVec.size()*particlesPerSet, "customManyParticleOrder");
        particleTypes->upload(particleTypesVec);
        orderIndex->upload(orderIndexVec);
        vector<int> flattenedOrder(particleOrder->getSize());
        for (int i = 0; i < (int) particleOrderVec.size(); i++)
            for (int j = 0; j < particlesPerSet; j++)
                flattenedOrder[i*particlesPerSet+j] = particleOrderVec[i][j];
        particleOrder->upload(flattenedOrder);
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
        exclusions = CudaArray::create<int>(cu, exclusionsVec.size(), "customManyParticleExclusions");
        exclusionStartIndex = CudaArray::create<int>(cu, exclusionStartIndexVec.size(), "customManyParticleExclusionStart");
        exclusions->upload(exclusionsVec);
        exclusionStartIndex->upload(exclusionStartIndexVec);
    }
    
    // Build data structures for the neighbor list.
    
    if (nonbondedMethod != NoCutoff) {
        int numAtomBlocks = cu.getNumAtomBlocks();
        int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter = new CudaArray(cu, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox = new CudaArray(cu, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        numNeighborPairs = CudaArray::create<int>(cu, 1, "customManyParticleNumNeighborPairs");
        neighborStartIndex = CudaArray::create<int>(cu, numParticles+1, "customManyParticleNeighborStartIndex");
        numNeighborsForAtom = CudaArray::create<int>(cu, numParticles, "customManyParticleNumNeighborsForAtom");
        CHECK_RESULT(cuEventCreate(&event, CU_EVENT_DISABLE_TIMING), "Error creating event for CustomManyParticleForce");

        // Select a size for the array that holds the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxNeighborPairs = 150*numParticles;
        neighborPairs = CudaArray::create<int2>(cu, maxNeighborPairs, "customManyParticleNeighborPairs");
        neighbors = CudaArray::create<int>(cu, maxNeighborPairs, "customManyParticleNeighbors");
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
        string index = cu.intToString(i+1);
        atomNames.push_back("P"+index);
        posNames.push_back("pos"+index);
    }
    stringstream compute;
    int index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize);\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
        variables.push_back(makeVariable(iter->first, "r_"+deltaName));
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
        if (computedDeltas.count(deltaName1) == 0) {
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<", periodicBoxSize, invPeriodicBoxSize);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<", periodicBoxSize, invPeriodicBoxSize);\n";
            computedDeltas.insert(deltaName2);
        }
        compute<<"real "<<angleName<<" = computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        variables.push_back(makeVariable(iter->first, angleName));
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
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
            compute<<"real4 delta"<<deltaName1<<" = delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize);\n";
            computedDeltas.insert(deltaName1);
        }
        if (computedDeltas.count(deltaName2) == 0) {
            compute<<"real4 delta"<<deltaName2<<" = delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize);\n";
            computedDeltas.insert(deltaName2);
        }
        if (computedDeltas.count(deltaName3) == 0) {
            compute<<"real4 delta"<<deltaName3<<" = delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[3]]<<", periodicBoxSize, invPeriodicBoxSize);\n";
            computedDeltas.insert(deltaName3);
        }
        compute<<"real4 "<<crossName1<<" = computeCross(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
        compute<<"real4 "<<crossName2<<" = computeCross(delta"<<deltaName2<<", delta"<<deltaName3<<");\n";
        compute<<"real "<<dihedralName<<" = computeAngle("<<crossName1<<", "<<crossName2<<");\n";
        compute<<dihedralName<<" *= (delta"<<deltaName1<<".x*"<<crossName2<<".x + delta"<<deltaName1<<".y*"<<crossName2<<".y + delta"<<deltaName1<<".z*"<<crossName2<<".z < 0 ? -1 : 1);\n";
        variables.push_back(makeVariable(iter->first, dihedralName));
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        compute<<buffer.getType()<<" params"<<(i+1)<<" = global_params"<<(i+1)<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    compute << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Apply forces to atoms.

    vector<string> forceNames;
    for (int i = 0; i < particlesPerSet; i++) {
        string istr = cu.intToString(i+1);
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
            compute<<cu.getExpressionUtilities().createExpressions(expressions, variables, functionList, functionDefinitions, "coordtemp");
        compute<<"}\n";
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*trim(delta"+deltaName+")";
        compute<<forceNames[atoms[0]]<<" += "<<"-"<<value<<";\n";
        compute<<forceNames[atoms[1]]<<" += "<<value<<";\n";
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
        compute<<"{\n";
        compute<<"real3 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(trim(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(trim(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
        compute<<"real3 deltaCross1 = -(deltaCross0+deltaCross2);\n";
        compute<<forceNames[atoms[0]]<<" += deltaCross0;\n";
        compute<<forceNames[atoms[1]]<<" += deltaCross1;\n";
        compute<<forceNames[atoms[2]]<<" += deltaCross2;\n";
        compute<<"}\n";
    }
    index = 0;
    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
        const vector<int>& atoms = iter->second;
        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
        compute<<"{\n";
        compute<<"real r = sqrt(delta"<<deltaName2<<".w);\n";
        compute<<"real4 ff;\n";
        compute<<"ff.x = (-dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real3 internalF0 = ff.x*trim("<<crossName1<<");\n";
        compute<<"real3 internalF3 = ff.w*trim("<<crossName2<<");\n";
        compute<<"real3 s = ff.y*internalF0 - ff.z*internalF3;\n";
        compute<<forceNames[atoms[0]]<<" += internalF0;\n";
        compute<<forceNames[atoms[1]]<<" += s-internalF0;\n";
        compute<<forceNames[atoms[2]]<<" += -s-internalF3;\n";
        compute<<forceNames[atoms[3]]<<" += internalF3;\n";
        compute<<"}\n";
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
        loadData<<"real3 pos"<<(i+1)<<" = trim(posq[atom"<<(i+1)<<"]);\n";
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
            verifyCutoff<<"real3 pos"<<(i+1)<<" = trim(posq[p"<<(i+1)<<"]);\n";
        if (!centralParticleMode) {
            for (int i = 1; i < particlesPerSet; i++) {
                for (int j = i+1; j < particlesPerSet; j++)
                    verifyCutoff<<"includeInteraction &= (delta(pos"<<(i+1)<<", pos"<<(j+1)<<", periodicBoxSize, invPeriodicBoxSize).w < CUTOFF_SQUARED);\n";
            }
        }
    }
    if (force.getNumExclusions() > 0) {
        int startCheckFrom = (nonbondedMethod == NoCutoff ? 0 : 1);
        for (int i = startCheckFrom; i < particlesPerSet; i++)
            for (int j = i+1; j < particlesPerSet; j++)
                verifyExclusions<<"includeInteraction &= !isInteractionExcluded(p"<<(i+1)<<", p"<<(j+1)<<", exclusions, exclusionStartIndex);\n";
    }
    string computeTypeIndex = "particleTypes[p"+cu.intToString(particlesPerSet)+"]";
    for (int i = particlesPerSet-2; i >= 0; i--)
        computeTypeIndex = "particleTypes[p"+cu.intToString(i+1)+"]+"+cu.intToString(numTypes)+"*("+computeTypeIndex+")";
    
    // Create replacements for extra arguments.
    
    stringstream extraArgs;
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArgs<<", const "<<buffer.getType()<<"* __restrict__ global_params"<<(i+1);
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
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["M_PI"] = cu.doubleToString(M_PI);
    defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    defines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
    defines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
    defines["NUM_GLOBALS"] = cu.intToString(max(1, force.getNumGlobalParameters()));
    defines["FIND_NEIGHBORS_WORKGROUP_SIZE"] = cu.intToString(findNeighborsWorkgroupSize);
    CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaKernelSources::customManyParticle, replacements), defines);
    forceKernel = cu.getKernel(module, "computeInteraction");
    blockBoundsKernel = cu.getKernel(module, "findBlockBounds");
    neighborsKernel = cu.getKernel(module, "findNeighbors");
    startIndicesKernel = cu.getKernel(module, "computeNeighborStartIndices");
    copyPairsKernel = cu.getKernel(module, "copyPairsToNeighborList");
    cuFuncSetCacheConfig(forceKernel, CU_FUNC_CACHE_PREFER_L1);
    cuFuncSetCacheConfig(neighborsKernel, CU_FUNC_CACHE_PREFER_L1);
    size_t bytes;
    CHECK_RESULT(cuModuleGetGlobal(&globalsPtr, &bytes, module, "globals"), "Error getting address for constant memory")
    cuMemcpyHtoD(globalsPtr, &globalParamValues[0], globalParamValues.size()*sizeof(float));
}

double CudaCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        
        // Set arguments for the force kernel.
        
        forceArgs.push_back(&cu.getForce().getDevicePointer());
        forceArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        forceArgs.push_back(&cu.getPosq().getDevicePointer());
        forceArgs.push_back(cu.getPeriodicBoxSizePointer());
        forceArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        if (nonbondedMethod != NoCutoff) {
            forceArgs.push_back(&neighbors->getDevicePointer());
            forceArgs.push_back(&neighborStartIndex->getDevicePointer());
        }
        if (particleTypes != NULL) {
            forceArgs.push_back(&particleTypes->getDevicePointer());
            forceArgs.push_back(&orderIndex->getDevicePointer());
            forceArgs.push_back(&particleOrder->getDevicePointer());
        }
        if (exclusions != NULL) {
            forceArgs.push_back(&exclusions->getDevicePointer());
            forceArgs.push_back(&exclusionStartIndex->getDevicePointer());
        }
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            forceArgs.push_back(&buffer.getMemory());
        }
        for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
            forceArgs.push_back(&tabulatedFunctions[i]->getDevicePointer());
        
        if (nonbondedMethod != NoCutoff) {
            // Set arguments for the block bounds kernel.

            blockBoundsArgs.push_back(cu.getPeriodicBoxSizePointer());
            blockBoundsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            blockBoundsArgs.push_back(&cu.getPosq().getDevicePointer());
            blockBoundsArgs.push_back(&blockCenter->getDevicePointer());
            blockBoundsArgs.push_back(&blockBoundingBox->getDevicePointer());
            blockBoundsArgs.push_back(&numNeighborPairs->getDevicePointer());

            // Set arguments for the neighbor list kernel.

            neighborsArgs.push_back(cu.getPeriodicBoxSizePointer());
            neighborsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            neighborsArgs.push_back(&cu.getPosq().getDevicePointer());
            neighborsArgs.push_back(&blockCenter->getDevicePointer());
            neighborsArgs.push_back(&blockBoundingBox->getDevicePointer());
            neighborsArgs.push_back(&neighborPairs->getDevicePointer());
            neighborsArgs.push_back(&numNeighborPairs->getDevicePointer());
            neighborsArgs.push_back(&numNeighborsForAtom->getDevicePointer());
            neighborsArgs.push_back(&maxNeighborPairs);
            if (exclusions != NULL) {
                neighborsArgs.push_back(&exclusions->getDevicePointer());
                neighborsArgs.push_back(&exclusionStartIndex->getDevicePointer());
            }
            
            // Set arguments for the kernel to find neighbor list start indices.
            
            startIndicesArgs.push_back(&numNeighborsForAtom->getDevicePointer());
            startIndicesArgs.push_back(&neighborStartIndex->getDevicePointer());
            startIndicesArgs.push_back(&numNeighborPairs->getDevicePointer());
            startIndicesArgs.push_back(&maxNeighborPairs);

            // Set arguments for the kernel to assemble the final neighbor list.
            
            copyPairsArgs.push_back(&neighborPairs->getDevicePointer());
            copyPairsArgs.push_back(&neighbors->getDevicePointer());
            copyPairsArgs.push_back(&numNeighborPairs->getDevicePointer());
            copyPairsArgs.push_back(&maxNeighborPairs);
            copyPairsArgs.push_back(&numNeighborsForAtom->getDevicePointer());
            copyPairsArgs.push_back(&neighborStartIndex->getDevicePointer());
       }
    }
    if (globalParamValues.size() > 0) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i])
                changed = true;
            globalParamValues[i] = value;
        }
        if (changed)
            cuMemcpyHtoD(globalsPtr, &globalParamValues[0], globalParamValues.size()*sizeof(float));
    }
    while (true) {
        int* numPairs = (int*) cu.getPinnedBuffer();
        if (nonbondedMethod != NoCutoff) {
            cu.executeKernel(blockBoundsKernel, &blockBoundsArgs[0], cu.getNumAtomBlocks());
            cu.executeKernel(neighborsKernel, &neighborsArgs[0], cu.getNumAtoms(), findNeighborsWorkgroupSize);

            // We need to make sure there was enough memory for the neighbor list.  Download the
            // information asynchronously so kernels can be running at the same time.

            numNeighborPairs->download(numPairs, false);
            CHECK_RESULT(cuEventRecord(event, 0), "Error recording event for CustomManyParticleForce");
            cu.executeKernel(startIndicesKernel, &startIndicesArgs[0], 256, 256, 256*sizeof(int));
            cu.executeKernel(copyPairsKernel, &copyPairsArgs[0], maxNeighborPairs);
        }
        int maxThreads = min(cu.getNumAtoms()*forceWorkgroupSize, cu.getEnergyBuffer().getSize());
        cu.executeKernel(forceKernel, &forceArgs[0], maxThreads, forceWorkgroupSize);
        if (nonbondedMethod != NoCutoff) {
            // Make sure there was enough memory for the neighbor list.

            CHECK_RESULT(cuEventSynchronize(event), "Error synchronizing on event for CustomManyParticleForce");
            if (*numPairs > maxNeighborPairs) {
                // Resize the arrays and run the calculation again.

                delete neighborPairs;
                neighborPairs = NULL;
                delete neighbors;
                neighbors = NULL;
                maxNeighborPairs = (int) (1.1*(*numPairs));
                neighborPairs = CudaArray::create<int2>(cu, maxNeighborPairs, "customManyParticleNeighborPairs");
                neighbors = CudaArray::create<int>(cu, maxNeighborPairs, "customManyParticleNeighbors");
                forceArgs[5] = &neighbors->getDevicePointer();
                neighborsArgs[5] = &neighborPairs->getDevicePointer();
                copyPairsArgs[0] = &neighborPairs->getDevicePointer();
                copyPairsArgs[1] = &neighbors->getDevicePointer();
                continue;
            }
        }
        break;
    }
    return 0.0;
}

void CudaCalcCustomManyParticleForceKernel::copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) {
    cu.setAsCurrent();
    int numParticles = force.getNumParticles();
    if (numParticles != cu.getNumAtoms())
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
    
    cu.invalidateMolecules();
}

CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
}

void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::verlet, defines, "");
    kernel1 = cu.getKernel(module, "integrateVerletPart1");
    kernel2 = cu.getKernel(module, "integrateVerletPart2");
    prevStepSize = -1.0;
}

void CudaIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    if (dt != prevStepSize) {
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            vector<double2> stepSizeVec(1);
            stepSizeVec[0] = make_double2(dt, dt);
            cu.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        else {
            vector<float2> stepSizeVec(1);
            stepSizeVec[0] = make_float2((float) dt, (float) dt);
            cu.getIntegrationUtilities().getStepSize().upload(stepSizeVec);
        }
        prevStepSize = dt;
    }

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::langevin, defines, "");
    kernel1 = cu.getKernel(module, "integrateLangevinPart1");
    kernel2 = cu.getKernel(module, "integrateLangevinPart2");
    params = new CudaArray(cu, 3, cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float), "langevinParams");
    prevStepSize = -1.0;
}

void CudaIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize/tau);
        double fscale = (1-vscale)*tau;
        double noisescale = sqrt(2*kT/tau)*sqrt(0.5*(1-vscale*vscale)*tau);
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            vector<double> p(params->getSize());
            p[0] = vscale;
            p[1] = fscale;
            p[2] = noisescale;
            params->upload(p);
            double2 ss = make_double2(0, stepSize);
            integration.getStepSize().upload(&ss);
        }
        else {
            vector<float> p(params->getSize());
            p[0] = (float) vscale;
            p[1] = (float) fscale;
            p[2] = (float) noisescale;
            params->upload(p);
            float2 ss = make_float2(0, (float) stepSize);
            integration.getStepSize().upload(&ss);
        }
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Call the first integration kernel.

    int randomIndex = integration.prepareRandomNumbers(cu.getPaddedNumAtoms());
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &params->getDevicePointer(), &integration.getStepSize().getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &posCorrection, &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getStepSize().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+stepSize);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
}

void CudaIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::brownian, defines, "");
    kernel1 = cu.getKernel(module, "integrateBrownianPart1");
    kernel2 = cu.getKernel(module, "integrateBrownianPart2");
    prevStepSize = -1.0;
}

void CudaIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
    double tauDt = tau*stepSize;
    double noise = sqrt(2.0f*BOLTZ*temperature*stepSize*tau);
    float stepSizeFloat = (float) stepSize;
    float tauDtFloat = (float) tauDt;
    float noiseFloat = (float) noise;
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();

    // Call the first integration kernel.

    int randomIndex = integration.prepareRandomNumbers(cu.getPaddedNumAtoms());
    void* args1[] = {&numAtoms, &paddedNumAtoms, useDouble ? (void*) &tauDt : (void*) &tauDtFloat,
            useDouble ? (void*) &noise : (void*) &noiseFloat,
            &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&numAtoms, useDouble ? (void*) &stepSize : (void*) &stepSizeFloat,
            &cu.getPosq().getDevicePointer(), &posCorrection, &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+stepSize);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateBrownianStepKernel::computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0);
}

CudaIntegrateVariableVerletStepKernel::~CudaIntegrateVariableVerletStepKernel() {
}

void CudaIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::verlet, defines, "");
    kernel1 = cu.getKernel(module, "integrateVerletPart1");
    kernel2 = cu.getKernel(module, "integrateVerletPart2");
    selectSizeKernel = cu.getKernel(module, "selectVerletStepSize");
    blockSize = min(256, system.getNumParticles());
}

double CudaIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();

    // Select the step size to use.

    double maxStepSize = maxTime-cu.getTime();
    float maxStepSizeFloat = (float) maxStepSize;
    double tol = integrator.getErrorTolerance();
    float tolFloat = (float) tol;
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();
    void* argsSelect[] = {&numAtoms, &paddedNumAtoms, useDouble ? (void*) &maxStepSize : (void*) &maxStepSizeFloat,
            useDouble ? (void*) &tol : (void*) &tolFloat,
            &cu.getIntegrationUtilities().getStepSize().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer()};
    int sharedSize = blockSize*(useDouble ? sizeof(double) : sizeof(float));
    cu.executeKernel(selectSizeKernel, argsSelect, blockSize, blockSize, sharedSize);

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    double dt, time;
    if (useDouble) {
        double2 stepSize;
        cu.getIntegrationUtilities().getStepSize().download(&stepSize);
        dt = stepSize.y;
        time = cu.getTime()+dt;
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
        float2 stepSize;
        cu.getIntegrationUtilities().getStepSize().download(&stepSize);
        dt = stepSize.y;
        time = cu.getTime()+dt;
        if (dt == maxStepSizeFloat)
            time = maxTime; // Avoid round-off error
    }
    cu.setTime(time);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
    return dt;
}

double CudaIntegrateVariableVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

CudaIntegrateVariableLangevinStepKernel::~CudaIntegrateVariableLangevinStepKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::langevin, defines, "");
    kernel1 = cu.getKernel(module, "integrateLangevinPart1");
    kernel2 = cu.getKernel(module, "integrateLangevinPart2");
    selectSizeKernel = cu.getKernel(module, "selectLangevinStepSize");
    params = new CudaArray(cu, 3, cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float), "langevinParams");
    blockSize = min(256, system.getNumParticles());
    blockSize = max(blockSize, params->getSize());
}

double CudaIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();

    // Select the step size to use.

    double maxStepSize = maxTime-cu.getTime();
    float maxStepSizeFloat = (float) maxStepSize;
    double tol = integrator.getErrorTolerance();
    float tolFloat = (float) tol;
    double tau = integrator.getFriction() == 0.0 ? 0.0 : 1.0/integrator.getFriction();
    float tauFloat = (float) tau;
    double kT = BOLTZ*integrator.getTemperature();
    float kTFloat = (float) kT;
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();
    void* argsSelect[] = {&numAtoms, &paddedNumAtoms, useDouble ? (void*) &maxStepSize : (void*) &maxStepSizeFloat,
            useDouble ? (void*) &tol : (void*) &tolFloat,
            useDouble ? (void*) &tau : (void*) &tauFloat,
            useDouble ? (void*) &kT : (void*) &kTFloat,
            &cu.getIntegrationUtilities().getStepSize().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &params->getDevicePointer()};
    int sharedSize = 2*blockSize*(useDouble ? sizeof(double) : sizeof(float));
    cu.executeKernel(selectSizeKernel, argsSelect, blockSize, blockSize, sharedSize);

    // Call the first integration kernel.

    int randomIndex = integration.prepareRandomNumbers(cu.getPaddedNumAtoms());
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &params->getDevicePointer(), &integration.getStepSize().getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &posCorrection, &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getStepSize().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    double dt, time;
    if (useDouble) {
        double2 stepSize;
        cu.getIntegrationUtilities().getStepSize().download(&stepSize);
        dt = stepSize.y;
        time = cu.getTime()+dt;
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
        float2 stepSize;
        cu.getIntegrationUtilities().getStepSize().download(&stepSize);
        dt = stepSize.y;
        time = cu.getTime()+dt;
        if (dt == maxStepSizeFloat)
            time = maxTime; // Avoid round-off error
    }
    cu.setTime(time);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
    return dt;
}

double CudaIntegrateVariableLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

class CudaIntegrateCustomStepKernel::ReorderListener : public CudaContext::ReorderListener {
public:
    ReorderListener(CudaContext& cu, CudaParameterSet& perDofValues, vector<vector<float> >& localPerDofValuesFloat, vector<vector<double> >& localPerDofValuesDouble, bool& deviceValuesAreCurrent) :
            cu(cu), perDofValues(perDofValues), localPerDofValuesFloat(localPerDofValuesFloat), localPerDofValuesDouble(localPerDofValuesDouble), deviceValuesAreCurrent(deviceValuesAreCurrent) {
        int numAtoms = cu.getNumAtoms();
        lastAtomOrder.resize(numAtoms);
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = cu.getAtomIndex()[i];
    }
    void execute() {
        // Reorder the per-DOF variables to reflect the new atom order.

        if (perDofValues.getNumParameters() == 0)
            return;
        int numAtoms = cu.getNumAtoms();
        const vector<int>& order = cu.getAtomIndex();
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            if (deviceValuesAreCurrent)
                perDofValues.getParameterValues(localPerDofValuesDouble);
            vector<vector<double> > swap(3*numAtoms);
            for (int i = 0; i < numAtoms; i++) {
                swap[3*lastAtomOrder[i]] = localPerDofValuesDouble[3*i];
                swap[3*lastAtomOrder[i]+1] = localPerDofValuesDouble[3*i+1];
                swap[3*lastAtomOrder[i]+2] = localPerDofValuesDouble[3*i+2];
            }
            for (int i = 0; i < numAtoms; i++) {
                localPerDofValuesDouble[3*i] = swap[3*order[i]];
                localPerDofValuesDouble[3*i+1] = swap[3*order[i]+1];
                localPerDofValuesDouble[3*i+2] = swap[3*order[i]+2];
            }
            perDofValues.setParameterValues(localPerDofValuesDouble);
        }
        else {
            if (deviceValuesAreCurrent)
                perDofValues.getParameterValues(localPerDofValuesFloat);
            vector<vector<float> > swap(3*numAtoms);
            for (int i = 0; i < numAtoms; i++) {
                swap[3*lastAtomOrder[i]] = localPerDofValuesFloat[3*i];
                swap[3*lastAtomOrder[i]+1] = localPerDofValuesFloat[3*i+1];
                swap[3*lastAtomOrder[i]+2] = localPerDofValuesFloat[3*i+2];
            }
            for (int i = 0; i < numAtoms; i++) {
                localPerDofValuesFloat[3*i] = swap[3*order[i]];
                localPerDofValuesFloat[3*i+1] = swap[3*order[i]+1];
                localPerDofValuesFloat[3*i+2] = swap[3*order[i]+2];
            }
            perDofValues.setParameterValues(localPerDofValuesFloat);
        }
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = order[i];
        deviceValuesAreCurrent = true;
    }
private:
    CudaContext& cu;
    CudaParameterSet& perDofValues;
    vector<vector<float> >& localPerDofValuesFloat;
    vector<vector<double> >& localPerDofValuesDouble;
    bool& deviceValuesAreCurrent;
    vector<int> lastAtomOrder;
};

CudaIntegrateCustomStepKernel::~CudaIntegrateCustomStepKernel() {
    cu.setAsCurrent();
    if (globalValues != NULL)
        delete globalValues;
    if (contextParameterValues != NULL)
        delete contextParameterValues;
    if (sumBuffer != NULL)
        delete sumBuffer;
    if (potentialEnergy != NULL)
        delete potentialEnergy;
    if (kineticEnergy != NULL)
        delete kineticEnergy;
    if (uniformRandoms != NULL)
        delete uniformRandoms;
    if (randomSeed != NULL)
        delete randomSeed;
    if (perDofValues != NULL)
        delete perDofValues;
    for (map<int, CudaArray*>::iterator iter = savedForces.begin(); iter != savedForces.end(); ++iter)
        delete iter->second;
}

void CudaIntegrateCustomStepKernel::initialize(const System& system, const CustomIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    numGlobalVariables = integrator.getNumGlobalVariables();
    int elementSize = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    globalValues = new CudaArray(cu, max(1, numGlobalVariables), elementSize, "globalVariables");
    sumBuffer = new CudaArray(cu, ((3*system.getNumParticles()+3)/4)*4, elementSize, "sumBuffer");
    potentialEnergy = new CudaArray(cu, 1, cu.getEnergyBuffer().getElementSize(), "potentialEnergy");
    kineticEnergy = new CudaArray(cu, 1, elementSize, "kineticEnergy");
    perDofValues = new CudaParameterSet(cu, integrator.getNumPerDofVariables(), 3*system.getNumParticles(), "perDofVariables", false, cu.getUseDoublePrecision() || cu.getUseMixedPrecision());
    cu.addReorderListener(new ReorderListener(cu, *perDofValues, localPerDofValuesFloat, localPerDofValuesDouble, deviceValuesAreCurrent));
    prevStepSize = -1.0;
    SimTKOpenMMUtilities::setRandomNumberSeed(integrator.getRandomNumberSeed());
}

string CudaIntegrateCustomStepKernel::createGlobalComputation(const string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator, const string& energyName) {
    map<string, Lepton::ParsedExpression> expressions;
    if (variable == "dt")
        expressions["dt[0].y = "] = expr;
    else {
        for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
            if (variable == integrator.getGlobalVariableName(i))
                expressions["globals["+cu.intToString(i)+"] = "] = expr;
        for (int i = 0; i < (int) parameterNames.size(); i++)
            if (variable == parameterNames[i]) {
                expressions["params["+cu.intToString(i)+"] = "] = expr;
                modifiesParameters = true;
            }
    }
    if (expressions.size() == 0)
        throw OpenMMException("Unknown global variable: "+variable);
    map<string, string> variables;
    variables["dt"] = "dt[0].y";
    variables["uniform"] = "uniform";
    variables["gaussian"] = "gaussian";
    variables[energyName] = "energy[0]";
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        variables[integrator.getGlobalVariableName(i)] = "globals["+cu.intToString(i)+"]";
    for (int i = 0; i < (int) parameterNames.size(); i++)
        variables[parameterNames[i]] = "params["+cu.intToString(i)+"]";
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    return cu.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp");
}

string CudaIntegrateCustomStepKernel::createPerDofComputation(const string& variable, const Lepton::ParsedExpression& expr, int component, CustomIntegrator& integrator, const string& forceName, const string& energyName) {
    const string suffixes[] = {".x", ".y", ".z"};
    string suffix = suffixes[component];
    map<string, Lepton::ParsedExpression> expressions;
    if (variable == "x")
        expressions["position"+suffix+" = "] = expr;
    else if (variable == "v")
        expressions["velocity"+suffix+" = "] = expr;
    else if (variable == "")
        expressions["sum[3*index+"+cu.intToString(component)+"] = "] = expr;
    else {
        for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
            if (variable == integrator.getPerDofVariableName(i))
                expressions["perDof"+suffix.substr(1)+perDofValues->getParameterSuffix(i)+" = "] = expr;
    }
    if (expressions.size() == 0)
        throw OpenMMException("Unknown per-DOF variable: "+variable);
    map<string, string> variables;
    variables["x"] = "position"+suffix;
    variables["v"] = "velocity"+suffix;
    variables[forceName] = "f"+suffix;
    variables["gaussian"] = "gaussian"+suffix;
    variables["uniform"] = "uniform"+suffix;
    variables["m"] = "mass";
    variables["dt"] = "stepSize";
    if (energyName != "")
        variables[energyName] = "energy[0]";
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        variables[integrator.getGlobalVariableName(i)] = "globals["+cu.intToString(i)+"]";
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        variables[integrator.getPerDofVariableName(i)] = "perDof"+suffix.substr(1)+perDofValues->getParameterSuffix(i);
    for (int i = 0; i < (int) parameterNames.size(); i++)
        variables[parameterNames[i]] = "params["+cu.intToString(i)+"]";
    vector<const TabulatedFunction*> functions;
    vector<pair<string, string> > functionNames;
    return cu.getExpressionUtilities().createExpressions(expressions, variables, functions, functionNames, "temp"+cu.intToString(component)+"_", "double");
}

void CudaIntegrateCustomStepKernel::prepareForComputation(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // Initialize various data structures.
        
        const map<string, double>& params = context.getParameters();
        if (useDouble) {
            contextParameterValues = CudaArray::create<double>(cu, max(1, (int) params.size()), "contextParameters");
            contextValuesDouble.resize(contextParameterValues->getSize());
            for (map<string, double>::const_iterator iter = params.begin(); iter != params.end(); ++iter) {
                contextValuesDouble[parameterNames.size()] = iter->second;
                parameterNames.push_back(iter->first);
            }
            contextParameterValues->upload(contextValuesDouble);
        }
        else {
            contextParameterValues = CudaArray::create<float>(cu, max(1, (int) params.size()), "contextParameters");
            contextValuesFloat.resize(contextParameterValues->getSize());
            for (map<string, double>::const_iterator iter = params.begin(); iter != params.end(); ++iter) {
                contextValuesFloat[parameterNames.size()] = (float) iter->second;
                parameterNames.push_back(iter->first);
            }
            contextParameterValues->upload(contextValuesFloat);
        }
        kernels.resize(integrator.getNumComputations());
        kernelArgs.resize(integrator.getNumComputations());
        requiredGaussian.resize(integrator.getNumComputations(), 0);
        requiredUniform.resize(integrator.getNumComputations(), 0);
        needsForces.resize(numSteps, false);
        needsEnergy.resize(numSteps, false);
        forceGroup.resize(numSteps, -2);
        invalidatesForces.resize(numSteps, false);
        merged.resize(numSteps, false);
        modifiesParameters = false;
        map<string, string> defines;
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        defines["WORK_GROUP_SIZE"] = cu.intToString(CudaContext::ThreadBlockSize);
        defines["SUM_BUFFER_SIZE"] = "0";
        defines["SUM_OUTPUT_INDEX"] = "0";
        
        // Build a list of all variables that affect the forces, so we can tell which
        // steps invalidate them.
        
        set<string> affectsForce;
        affectsForce.insert("x");
        for (vector<ForceImpl*>::const_iterator iter = context.getForceImpls().begin(); iter != context.getForceImpls().end(); ++iter) {
            const map<string, double> params = (*iter)->getDefaultParameters();
            for (map<string, double>::const_iterator param = params.begin(); param != params.end(); ++param)
                affectsForce.insert(param->first);
        }
        
        // Record information about all the computation steps.
        
        stepType.resize(numSteps);
        vector<string> variable(numSteps);
        vector<Lepton::ParsedExpression> expression(numSteps);
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
        for (int step = 0; step < numSteps; step++) {
            string expr;
            integrator.getComputationStep(step, stepType[step], variable[step], expr);
            if (expr.size() > 0) {
                expression[step] = Lepton::Parser::parse(expr).optimize();
                if (usesVariable(expression[step], "f")) {
                    needsForces[step] = true;
                    forceGroup[step] = -1;
                }
                if (usesVariable(expression[step], "energy")) {
                    needsEnergy[step] = true;
                    forceGroup[step] = -1;
                }
                for (int i = 0; i < 32; i++) {
                    if (usesVariable(expression[step], forceGroupName[i])) {
                        if (forceGroup[step] != -2)
                            throw OpenMMException("A single computation step cannot depend on multiple force groups");
                        needsForces[step] = true;
                        forceGroup[step] = 1<<i;
                        forceName[step] = forceGroupName[i];
                    }
                    if (usesVariable(expression[step], energyGroupName[i])) {
                        if (forceGroup[step] != -2)
                            throw OpenMMException("A single computation step cannot depend on multiple force groups");
                        needsEnergy[step] = true;
                        forceGroup[step] = 1<<i;
                        energyName[step] = energyGroupName[i];
                    }
                }
            }
            invalidatesForces[step] = (stepType[step] == CustomIntegrator::ConstrainPositions || affectsForce.find(variable[step]) != affectsForce.end());
            if (forceGroup[step] == -2 && step > 0)
                forceGroup[step] = forceGroup[step-1];
            if (forceGroup[step] != -2 && savedForces.find(forceGroup[step]) == savedForces.end())
                savedForces[forceGroup[step]] = new CudaArray(cu, cu.getForce().getSize(), cu.getForce().getElementSize(), "savedForces");
        }
        
        // Determine how each step will represent the position (as just a value, or a value plus a delta).
        
        vector<bool> storePosAsDelta(numSteps, false);
        vector<bool> loadPosAsDelta(numSteps, false);
        bool beforeConstrain = false;
        for (int step = numSteps-1; step >= 0; step--) {
            if (stepType[step] == CustomIntegrator::ConstrainPositions)
                beforeConstrain = true;
            else if (stepType[step] == CustomIntegrator::ComputePerDof && variable[step] == "x" && beforeConstrain)
                storePosAsDelta[step] = true;
        }
        bool storedAsDelta = false;
        for (int step = 0; step < numSteps; step++) {
            loadPosAsDelta[step] = storedAsDelta;
            if (storePosAsDelta[step] == true)
                storedAsDelta = true;
            if (stepType[step] == CustomIntegrator::ConstrainPositions)
                storedAsDelta = false;
        }
        
        // Identify steps that can be merged into a single kernel.
        
        for (int step = 1; step < numSteps; step++) {
            if (needsForces[step] || needsEnergy[step])
                continue;
            if (stepType[step-1] == CustomIntegrator::ComputeGlobal && stepType[step] == CustomIntegrator::ComputeGlobal &&
                    !usesVariable(expression[step], "uniform") && !usesVariable(expression[step], "gaussian"))
                merged[step] = true;
            if (stepType[step-1] == CustomIntegrator::ComputePerDof && stepType[step] == CustomIntegrator::ComputePerDof)
                merged[step] = true;
        }
        
        // Loop over all steps and create the kernels for them.
        
        for (int step = 0; step < numSteps; step++) {
            if ((stepType[step] == CustomIntegrator::ComputePerDof || stepType[step] == CustomIntegrator::ComputeSum) && !merged[step]) {
                // Compute a per-DOF value.
                
                stringstream compute;
                for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
                    CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
                    compute << buffer.getType()<<" perDofx"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index];\n";
                    compute << buffer.getType()<<" perDofy"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index+1];\n";
                    compute << buffer.getType()<<" perDofz"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index+2];\n";
                }
                int numGaussian = 0, numUniform = 0;
                for (int j = step; j < numSteps && (j == step || merged[j]); j++) {
                    numGaussian += numAtoms*usesVariable(expression[j], "gaussian");
                    numUniform += numAtoms*usesVariable(expression[j], "uniform");
                    compute << "{\n";
                    if (numGaussian > 0)
                        compute << "float4 gaussian = gaussianValues[gaussianIndex+index];\n";
                    if (numUniform > 0)
                        compute << "float4 uniform = uniformValues[uniformIndex+index];\n";
                    for (int i = 0; i < 3; i++)
                        compute << createPerDofComputation(stepType[j] == CustomIntegrator::ComputePerDof ? variable[j] : "", expression[j], i, integrator, forceName[j], energyName[j]);
                    if (variable[j] == "x") {
                        if (storePosAsDelta[j])
                            compute << "posDelta[index] = convertFromDouble4(position-convertToDouble4(loadPos(posq, posqCorrection, index)));\n";
                        else
                            compute << "storePos(posq, posqCorrection, index, convertFromDouble4(position));\n";
                    }
                    else if (variable[j] == "v")
                        compute << "velm[index] = convertFromDouble4(velocity);\n";
                    else {
                        for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
                            CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
                            compute << "perDofValues"<<cu.intToString(i+1)<<"[3*index] = perDofx"<<cu.intToString(i+1)<<";\n";
                            compute << "perDofValues"<<cu.intToString(i+1)<<"[3*index+1] = perDofy"<<cu.intToString(i+1)<<";\n";
                            compute << "perDofValues"<<cu.intToString(i+1)<<"[3*index+2] = perDofz"<<cu.intToString(i+1)<<";\n";
                        }
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
                for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
                    CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
                    string valueName = "perDofValues"+cu.intToString(i+1);
                    args << ", " << buffer.getType() << "* __restrict__ " << valueName;
                }
                replacements["PARAMETER_ARGUMENTS"] = args.str();
                if (loadPosAsDelta[step])
                    defines["LOAD_POS_AS_DELTA"] = "1";
                else if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
                    defines.erase("LOAD_POS_AS_DELTA");
                CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaKernelSources::customIntegratorPerDof, replacements), defines);
                CUfunction kernel = cu.getKernel(module, "computePerDof");
                kernels[step].push_back(kernel);
                requiredGaussian[step] = numGaussian;
                requiredUniform[step] = numUniform;
                vector<void*> args1;
                args1.push_back(&cu.getPosq().getDevicePointer());
                args1.push_back(NULL);
                args1.push_back(&integration.getPosDelta().getDevicePointer());
                args1.push_back(&cu.getVelm().getDevicePointer());
                args1.push_back(&cu.getForce().getDevicePointer());
                args1.push_back(&integration.getStepSize().getDevicePointer());
                args1.push_back(&globalValues->getDevicePointer());
                args1.push_back(&contextParameterValues->getDevicePointer());
                args1.push_back(&sumBuffer->getDevicePointer());
                args1.push_back(NULL);
                args1.push_back(NULL);
                args1.push_back(NULL);
                args1.push_back(&potentialEnergy->getDevicePointer());
                for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++)
                    args1.push_back(&perDofValues->getBuffers()[i].getMemory());
                kernelArgs[step].push_back(args1);
                if (stepType[step] == CustomIntegrator::ComputeSum) {
                    // Create a second kernel for this step that sums the values.

                    vector<void*> args2;
                    args2.push_back(&sumBuffer->getDevicePointer());
                    bool found = false;
                    for (int j = 0; j < integrator.getNumGlobalVariables() && !found; j++)
                        if (variable[step] == integrator.getGlobalVariableName(j)) {
                            args2.push_back(&globalValues->getDevicePointer());
                            defines["SUM_OUTPUT_INDEX"] = cu.intToString(j);
                            found = true;
                        }
                    for (int j = 0; j < (int) parameterNames.size() && !found; j++)
                        if (variable[step] == parameterNames[j]) {
                            args2.push_back(&contextParameterValues->getDevicePointer());
                            defines["SUM_OUTPUT_INDEX"] = cu.intToString(j);
                            found = true;
                            modifiesParameters = true;
                        }
                    if (!found)
                        throw OpenMMException("Unknown global variable: "+variable[step]);
                    defines["SUM_BUFFER_SIZE"] = cu.intToString(3*numAtoms);
                    module = cu.createModule(CudaKernelSources::customIntegrator, defines);
                    kernel = cu.getKernel(module, useDouble ? "computeDoubleSum" : "computeFloatSum");
                    kernels[step].push_back(kernel);
                    kernelArgs[step].push_back(args2);
                }
            }
            else if (stepType[step] == CustomIntegrator::ComputeGlobal && !merged[step]) {
                // Compute a global value.

                stringstream compute;
                for (int i = step; i < numSteps && (i == step || merged[i]); i++)
                    compute << "{\n" << createGlobalComputation(variable[i], expression[i], integrator, energyName[i]) << "}\n";
                map<string, string> replacements;
                replacements["COMPUTE_STEP"] = compute.str();
                CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::customIntegratorGlobal, replacements), defines);
                CUfunction kernel = cu.getKernel(module, "computeGlobal");
                kernels[step].push_back(kernel);
                vector<void*> args;
                args.push_back(&integration.getStepSize().getDevicePointer());
                args.push_back(&globalValues->getDevicePointer());
                args.push_back(&contextParameterValues->getDevicePointer());
                args.push_back(NULL);
                args.push_back(NULL);
                args.push_back(&potentialEnergy->getDevicePointer());
                kernelArgs[step].push_back(args);
            }
            else if (stepType[step] == CustomIntegrator::ConstrainPositions) {
                // Apply position constraints.

                CUmodule module = cu.createModule(CudaKernelSources::customIntegrator, defines);
                CUfunction kernel = cu.getKernel(module, "applyPositionDeltas");
                kernels[step].push_back(kernel);
                vector<void*> args;
                args.push_back(&cu.getPosq().getDevicePointer());
                args.push_back(NULL);
                args.push_back(&integration.getPosDelta().getDevicePointer());
                kernelArgs[step].push_back(args);
            }
        }
        
        // Initialize the random number generator.
        
        int maxUniformRandoms = 1;
        for (int i = 0; i < (int) requiredUniform.size(); i++)
            maxUniformRandoms = max(maxUniformRandoms, requiredUniform[i]);
        uniformRandoms = CudaArray::create<float4>(cu, maxUniformRandoms, "uniformRandoms");
        randomSeed = CudaArray::create<int4>(cu, cu.getNumThreadBlocks()*CudaContext::ThreadBlockSize, "randomSeed");
        vector<int4> seed(randomSeed->getSize());
        unsigned int r = integrator.getRandomNumberSeed()+1;
        for (int i = 0; i < randomSeed->getSize(); i++) {
            seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
            seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
        }
        randomSeed->upload(seed);
        CUmodule randomProgram = cu.createModule(CudaKernelSources::customIntegrator, defines);
        randomKernel = cu.getKernel(randomProgram, "generateRandomNumbers");
        
        // Create the kernel for summing the potential energy.

        defines["SUM_OUTPUT_INDEX"] = "0";
        defines["SUM_BUFFER_SIZE"] = cu.intToString(cu.getEnergyBuffer().getSize());
        CUmodule module = cu.createModule(CudaKernelSources::customIntegrator, defines);
        sumPotentialEnergyKernel = cu.getKernel(module, cu.getUseDoublePrecision() ? "computeDoubleSum" : "computeFloatSum");
        
        // Create the kernel for computing kinetic energy.

        stringstream computeKE;
        for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
            const CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
            computeKE << buffer.getType()<<" perDofx"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index];\n";
            computeKE << buffer.getType()<<" perDofy"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index+1];\n";
            computeKE << buffer.getType()<<" perDofz"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index+2];\n";
        }
        Lepton::ParsedExpression keExpression = Lepton::Parser::parse(integrator.getKineticEnergyExpression()).optimize();
        for (int i = 0; i < 3; i++)
            computeKE << createPerDofComputation("", keExpression, i, integrator, "f", "");
        map<string, string> replacements;
        replacements["COMPUTE_STEP"] = computeKE.str();
        stringstream args;
        for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
            const CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
            string valueName = "perDofValues"+cu.intToString(i+1);
            args << ", " << buffer.getType() << "* __restrict__ " << valueName;
        }
        replacements["PARAMETER_ARGUMENTS"] = args.str();
        defines["SUM_BUFFER_SIZE"] = cu.intToString(3*numAtoms);
        if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
            defines.erase("LOAD_POS_AS_DELTA");
        module = cu.createModule(cu.replaceStrings(CudaKernelSources::customIntegratorPerDof, replacements), defines);
        kineticEnergyKernel = cu.getKernel(module, "computePerDof");
        kineticEnergyArgs.push_back(&cu.getPosq().getDevicePointer());
        kineticEnergyArgs.push_back(NULL);
        kineticEnergyArgs.push_back(&integration.getPosDelta().getDevicePointer());
        kineticEnergyArgs.push_back(&cu.getVelm().getDevicePointer());
        kineticEnergyArgs.push_back(&cu.getForce().getDevicePointer());
        kineticEnergyArgs.push_back(&integration.getStepSize().getDevicePointer());
        kineticEnergyArgs.push_back(&globalValues->getDevicePointer());
        kineticEnergyArgs.push_back(&contextParameterValues->getDevicePointer());
        kineticEnergyArgs.push_back(&sumBuffer->getDevicePointer());
        kineticEnergyArgs.push_back(NULL);
        kineticEnergyArgs.push_back(NULL);
        kineticEnergyArgs.push_back(&uniformRandoms->getDevicePointer());
        kineticEnergyArgs.push_back(&potentialEnergy->getDevicePointer());
        for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++)
            kineticEnergyArgs.push_back(&perDofValues->getBuffers()[i].getMemory());
        keNeedsForce = usesVariable(keExpression, "f");

        // Create a second kernel to sum the values.

        defines["SUM_BUFFER_SIZE"] = cu.intToString(3*numAtoms);
        module = cu.createModule(CudaKernelSources::customIntegrator, defines);
        sumKineticEnergyKernel = cu.getKernel(module, useDouble ? "computeDoubleSum" : "computeFloatSum");
    }
    
    // Make sure all values (variables, parameters, etc.) stored on the device are up to date.
    
    if (!deviceValuesAreCurrent) {
        if (useDouble)
            perDofValues->setParameterValues(localPerDofValuesDouble);
        else
            perDofValues->setParameterValues(localPerDofValuesFloat);
        deviceValuesAreCurrent = true;
    }
    localValuesAreCurrent = false;
    double stepSize = integrator.getStepSize();
    if (stepSize != prevStepSize) {
        if (useDouble) {
            double size[] = {0, stepSize};
            integration.getStepSize().upload(size);
        }
        else {
            float size[] = {0, (float) stepSize};
            integration.getStepSize().upload(size);
        }
        prevStepSize = stepSize;
    }
    bool paramsChanged = false;
    if (useDouble) {
        for (int i = 0; i < (int) parameterNames.size(); i++) {
            double value = context.getParameter(parameterNames[i]);
            if (value != contextValuesDouble[i]) {
                contextValuesDouble[i] = value;
                paramsChanged = true;
            }
        }
        if (paramsChanged)
            contextParameterValues->upload(contextValuesDouble);
    }
    else {
        for (int i = 0; i < (int) parameterNames.size(); i++) {
            float value = (float) context.getParameter(parameterNames[i]);
            if (value != contextValuesFloat[i]) {
                contextValuesFloat[i] = value;
                paramsChanged = true;
            }
        }
        if (paramsChanged)
            contextParameterValues->upload(contextValuesFloat);
    }
}

void CudaIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    prepareForComputation(context, integrator, forcesAreValid);
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int numSteps = integrator.getNumComputations();

    // Loop over computation steps in the integrator and execute them.

    int maxUniformRandoms = uniformRandoms->getSize();
    void* randomArgs[] = {&maxUniformRandoms, &uniformRandoms->getDevicePointer(), &randomSeed->getDevicePointer()};
    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    for (int i = 0; i < numSteps; i++) {
        int lastForceGroups = context.getLastForceGroups();
        if ((needsForces[i] || needsEnergy[i]) && (!forcesAreValid || lastForceGroups != forceGroup[i])) {
            if (forcesAreValid && savedForces.find(lastForceGroups) != savedForces.end()) {
                // The forces are still valid.  We just need a different force group right now.  Save the old
                // forces in case we need them again.
                
                cu.getForce().copyTo(*savedForces[lastForceGroups]);
                validSavedForces.insert(lastForceGroups);
            }
            else
                validSavedForces.clear();
            
            // Recompute forces and/or energy.  Figure out what is actually needed
            // between now and the next time they get invalidated again.
            
            bool computeForce = false, computeEnergy = false;
            for (int j = i; ; j++) {
                if (needsForces[j])
                    computeForce = true;
                if (needsEnergy[j])
                    computeEnergy = true;
                if (invalidatesForces[j])
                    break;
                if (j == numSteps-1)
                    j = -1;
                if (j == i-1)
                    break;
            }
            if (!computeEnergy && validSavedForces.find(forceGroup[i]) != validSavedForces.end()) {
                // We can just restore the forces we saved earlier.
                
                savedForces[forceGroup[i]]->copyTo(cu.getForce());
            }
            else {
                recordChangedParameters(context);
                context.calcForcesAndEnergy(computeForce, computeEnergy, forceGroup[i]);
                if (computeEnergy) {
                    void* args[] = {&cu.getEnergyBuffer().getDevicePointer(), &potentialEnergy->getDevicePointer()};
                    cu.executeKernel(sumPotentialEnergyKernel, &args[0], CudaContext::ThreadBlockSize, CudaContext::ThreadBlockSize);
                }
            }
            forcesAreValid = true;
        }
        if (stepType[i] == CustomIntegrator::ComputePerDof && !merged[i]) {
            int randomIndex = integration.prepareRandomNumbers(requiredGaussian[i]);
            kernelArgs[i][0][1] = &posCorrection;
            kernelArgs[i][0][9] = &integration.getRandom().getDevicePointer();
            kernelArgs[i][0][10] = &randomIndex;
            kernelArgs[i][0][11] = &uniformRandoms->getDevicePointer();
            if (requiredUniform[i] > 0)
                cu.executeKernel(randomKernel, &randomArgs[0], numAtoms);
            cu.executeKernel(kernels[i][0], &kernelArgs[i][0][0], numAtoms);
        }
        else if (stepType[i] == CustomIntegrator::ComputeGlobal && !merged[i]) {
            float uniform = SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber();
            float gauss = SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
            kernelArgs[i][0][3] = &uniform;
            kernelArgs[i][0][4] = &gauss;
            cu.executeKernel(kernels[i][0], &kernelArgs[i][0][0], 1, 1);
        }
        else if (stepType[i] == CustomIntegrator::ComputeSum) {
            int randomIndex = integration.prepareRandomNumbers(requiredGaussian[i]);
            kernelArgs[i][0][1] = &posCorrection;
            kernelArgs[i][0][9] = &integration.getRandom().getDevicePointer();
            kernelArgs[i][0][10] = &randomIndex;
            kernelArgs[i][0][11] = &uniformRandoms->getDevicePointer();
            if (requiredUniform[i] > 0)
                cu.executeKernel(randomKernel, &randomArgs[0], numAtoms);
            cu.clearBuffer(*sumBuffer);
            cu.executeKernel(kernels[i][0], &kernelArgs[i][0][0], numAtoms);
            cu.executeKernel(kernels[i][1], &kernelArgs[i][1][0], CudaContext::ThreadBlockSize, CudaContext::ThreadBlockSize);
        }
        else if (stepType[i] == CustomIntegrator::UpdateContextState) {
            recordChangedParameters(context);
            context.updateContextState();
        }
        else if (stepType[i] == CustomIntegrator::ConstrainPositions) {
            cu.getIntegrationUtilities().applyConstraints(integrator.getConstraintTolerance());
            kernelArgs[i][0][1] = &posCorrection;
            cu.executeKernel(kernels[i][0], &kernelArgs[i][0][0], numAtoms);
            cu.getIntegrationUtilities().computeVirtualSites();
        }
        else if (stepType[i] == CustomIntegrator::ConstrainVelocities) {
            cu.getIntegrationUtilities().applyVelocityConstraints(integrator.getConstraintTolerance());
        }
        if (invalidatesForces[i])
            forcesAreValid = false;
    }
    recordChangedParameters(context);

    // Update the time and step count.

    cu.setTime(cu.getTime()+integrator.getStepSize());
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
    if (cu.getAtomsWereReordered()) {
        forcesAreValid = false;
        validSavedForces.clear();
    }
}

double CudaIntegrateCustomStepKernel::computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    prepareForComputation(context, integrator, forcesAreValid);
    if (keNeedsForce && !forcesAreValid) {
        // Compute the force.  We want to then mark that forces are valid, which means also computing
        // potential energy if any steps will expect it to be valid too.
        
        bool willNeedEnergy = false;
        for (int i = 0; i < integrator.getNumComputations(); i++)
            willNeedEnergy |= needsEnergy[i];
        context.calcForcesAndEnergy(true, willNeedEnergy, -1);
        if (willNeedEnergy) {
            void* args[] = {&cu.getEnergyBuffer().getDevicePointer(), &potentialEnergy->getDevicePointer()};
            cu.executeKernel(sumPotentialEnergyKernel, &args[0], CudaContext::ThreadBlockSize, CudaContext::ThreadBlockSize);
        }
        forcesAreValid = true;
    }
    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    int randomIndex = 0;
    kineticEnergyArgs[1] = &posCorrection;
    kineticEnergyArgs[9] = &cu.getIntegrationUtilities().getRandom().getDevicePointer();
    kineticEnergyArgs[10] = &randomIndex;
    cu.clearBuffer(*sumBuffer);
    cu.executeKernel(kineticEnergyKernel, &kineticEnergyArgs[0], cu.getNumAtoms());
    void* args[] = {&sumBuffer->getDevicePointer(), &kineticEnergy->getDevicePointer()};
    cu.executeKernel(sumKineticEnergyKernel, args, CudaContext::ThreadBlockSize, CudaContext::ThreadBlockSize);
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        double ke;
        kineticEnergy->download(&ke);
        return ke;
    }
    else {
        float ke;
        kineticEnergy->download(&ke);
        return ke;
    }
}

void CudaIntegrateCustomStepKernel::recordChangedParameters(ContextImpl& context) {
    if (!modifiesParameters)
        return;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        contextParameterValues->download(contextValuesDouble);
        for (int i = 0; i < (int) parameterNames.size(); i++) {
            double value = context.getParameter(parameterNames[i]);
            if (value != contextValuesDouble[i])
                context.setParameter(parameterNames[i], contextValuesDouble[i]);
        }
    }
    else {
        contextParameterValues->download(contextValuesFloat);
        for (int i = 0; i < (int) parameterNames.size(); i++) {
            float value = (float) context.getParameter(parameterNames[i]);
            if (value != contextValuesFloat[i])
                context.setParameter(parameterNames[i], contextValuesFloat[i]);
        }
    }
}

void CudaIntegrateCustomStepKernel::getGlobalVariables(ContextImpl& context, vector<double>& values) const {
    values.resize(numGlobalVariables);
    if (numGlobalVariables == 0)
        return;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        globalValues->download(values);
    else {
        vector<float> buffer;
        globalValues->download(buffer);
        for (int i = 0; i < numGlobalVariables; i++)
            values[i] = buffer[i];
    }
}

void CudaIntegrateCustomStepKernel::setGlobalVariables(ContextImpl& context, const vector<double>& values) {
    if (numGlobalVariables == 0)
        return;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        globalValues->upload(values);
    else {
        vector<float> buffer(numGlobalVariables);
        for (int i = 0; i < numGlobalVariables; i++)
            buffer[i] = (float) values[i];
        globalValues->upload(buffer);
    }
}

void CudaIntegrateCustomStepKernel::getPerDofVariable(ContextImpl& context, int variable, vector<Vec3>& values) const {
    values.resize(perDofValues->getNumObjects()/3);
    const vector<int>& order = cu.getAtomIndex();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        if (!localValuesAreCurrent) {
            perDofValues->getParameterValues(localPerDofValuesDouble);
            localValuesAreCurrent = true;
        }
        for (int i = 0; i < (int) values.size(); i++)
            for (int j = 0; j < 3; j++)
                values[order[i]][j] = localPerDofValuesDouble[3*i+j][variable];
    }
    else {
        if (!localValuesAreCurrent) {
            perDofValues->getParameterValues(localPerDofValuesFloat);
            localValuesAreCurrent = true;
        }
        for (int i = 0; i < (int) values.size(); i++)
            for (int j = 0; j < 3; j++)
                values[order[i]][j] = localPerDofValuesFloat[3*i+j][variable];
    }
}

void CudaIntegrateCustomStepKernel::setPerDofVariable(ContextImpl& context, int variable, const vector<Vec3>& values) {
    const vector<int>& order = cu.getAtomIndex();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        if (!localValuesAreCurrent) {
            perDofValues->getParameterValues(localPerDofValuesDouble);
            localValuesAreCurrent = true;
        }
        for (int i = 0; i < (int) values.size(); i++)
            for (int j = 0; j < 3; j++)
                localPerDofValuesDouble[3*i+j][variable] = values[order[i]][j];
    }
    else {
        if (!localValuesAreCurrent) {
            perDofValues->getParameterValues(localPerDofValuesFloat);
            localValuesAreCurrent = true;
        }
        for (int i = 0; i < (int) values.size(); i++)
            for (int j = 0; j < 3; j++)
                localPerDofValuesFloat[3*i+j][variable] = (float) values[order[i]][j];
    }
    deviceValuesAreCurrent = false;
}

CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
    cu.setAsCurrent();
    if (atomGroups != NULL)
        delete atomGroups;
}

void CudaApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    cu.setAsCurrent();
    randomSeed = thermostat.getRandomNumberSeed();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::andersenThermostat, defines);
    kernel = cu.getKernel(module, "applyAndersenThermostat");
    cu.getIntegrationUtilities().initRandomNumberGenerator(randomSeed);

    // Create the arrays with the group definitions.

    vector<vector<int> > groups = AndersenThermostatImpl::calcParticleGroups(system);
    atomGroups = CudaArray::create<int>(cu, cu.getNumAtoms(), "atomGroups");
    vector<int> atoms(atomGroups->getSize());
    for (int i = 0; i < (int) groups.size(); i++) {
        for (int j = 0; j < (int) groups[i].size(); j++)
            atoms[groups[i][j]] = i;
    }
    atomGroups->upload(atoms);
}

void CudaApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    cu.setAsCurrent();
    float frequency = (float) context.getParameter(AndersenThermostat::CollisionFrequency());
    float kT = (float) (BOLTZ*context.getParameter(AndersenThermostat::Temperature()));
    int randomIndex = cu.getIntegrationUtilities().prepareRandomNumbers(cu.getPaddedNumAtoms());
    int numAtoms = cu.getNumAtoms();
    void* args[] = {&numAtoms, &frequency, &kT, &cu.getVelm().getDevicePointer(), &cu.getIntegrationUtilities().getStepSize().getDevicePointer(),
            &cu.getIntegrationUtilities().getRandom().getDevicePointer(), &randomIndex, &atomGroups->getDevicePointer()};
    cu.executeKernel(kernel, args, cu.getNumAtoms());
}

CudaApplyMonteCarloBarostatKernel::~CudaApplyMonteCarloBarostatKernel() {
    cu.setAsCurrent();
    if (savedPositions != NULL)
        delete savedPositions;
    if (moleculeAtoms != NULL)
        delete moleculeAtoms;
    if (moleculeStartIndex != NULL)
        delete moleculeStartIndex;
}

void CudaApplyMonteCarloBarostatKernel::initialize(const System& system, const Force& thermostat) {
    cu.setAsCurrent();
    savedPositions = new CudaArray(cu, cu.getPaddedNumAtoms(), cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4), "savedPositions");
    CUmodule module = cu.createModule(CudaKernelSources::monteCarloBarostat);
    kernel = cu.getKernel(module, "scalePositions");
}

void CudaApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ) {
    cu.setAsCurrent();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;

        // Create the arrays with the molecule definitions.

        vector<vector<int> > molecules = context.getMolecules();
        numMolecules = molecules.size();
        moleculeAtoms = CudaArray::create<int>(cu, cu.getNumAtoms(), "moleculeAtoms");
        moleculeStartIndex = CudaArray::create<int>(cu, numMolecules+1, "moleculeStartIndex");
        vector<int> atoms(moleculeAtoms->getSize());
        vector<int> startIndex(moleculeStartIndex->getSize());
        int index = 0;
        for (int i = 0; i < numMolecules; i++) {
            startIndex[i] = index;
            for (int j = 0; j < (int) molecules[i].size(); j++)
                atoms[index++] = molecules[i][j];
        }
        startIndex[numMolecules] = index;
        moleculeAtoms->upload(atoms);
        moleculeStartIndex->upload(startIndex);

        // Initialize the kernel arguments.
        
    }
    int bytesToCopy = cu.getPosq().getSize()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4));
    CUresult result = cuMemcpyDtoD(savedPositions->getDevicePointer(), cu.getPosq().getDevicePointer(), bytesToCopy);
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error saving positions for MC barostat: "<<cu.getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }        
    float scalefX = (float) scaleX;
    float scalefY = (float) scaleY;
    float scalefZ = (float) scaleZ;
    void* args[] = {&scalefX, &scalefY, &scalefZ, &numMolecules, cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer(), 
		    &cu.getPosq().getDevicePointer(), &moleculeAtoms->getDevicePointer(), &moleculeStartIndex->getDevicePointer()};
    cu.executeKernel(kernel, args, cu.getNumAtoms());
    for (int i = 0; i < (int) cu.getPosCellOffsets().size(); i++)
        cu.getPosCellOffsets()[i] = make_int4(0, 0, 0, 0);
    lastAtomOrder = cu.getAtomIndex();
}

void CudaApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    cu.setAsCurrent();
    int bytesToCopy = cu.getPosq().getSize()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4));
    CUresult result = cuMemcpyDtoD(cu.getPosq().getDevicePointer(), savedPositions->getDevicePointer(), bytesToCopy);
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error restoring positions for MC barostat: "<<cu.getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
}

CudaRemoveCMMotionKernel::~CudaRemoveCMMotionKernel() {
    cu.setAsCurrent();
    if (cmMomentum != NULL)
        delete cmMomentum;
}

void CudaRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    cu.setAsCurrent();
    frequency = force.getFrequency();
    int numAtoms = cu.getNumAtoms();
    cmMomentum = CudaArray::create<float4>(cu, (numAtoms+CudaContext::ThreadBlockSize-1)/CudaContext::ThreadBlockSize, "cmMomentum");
    double totalMass = 0.0;
    for (int i = 0; i < numAtoms; i++)
        totalMass += system.getParticleMass(i);
    map<string, string> defines;
    defines["INVERSE_TOTAL_MASS"] = cu.doubleToString(totalMass == 0 ? 0.0 : 1.0/totalMass);
    CUmodule module = cu.createModule(CudaKernelSources::removeCM, defines);
    kernel1 = cu.getKernel(module, "calcCenterOfMassMomentum");
    kernel2 = cu.getKernel(module, "removeCenterOfMassMomentum");
}

void CudaRemoveCMMotionKernel::execute(ContextImpl& context) {
    cu.setAsCurrent();
    int numAtoms = cu.getNumAtoms();
    void* args[] = {&numAtoms, &cu.getVelm().getDevicePointer(), &cmMomentum->getDevicePointer()};
    cu.executeKernel(kernel1, args, cu.getNumAtoms(), cu.ThreadBlockSize, cu.ThreadBlockSize*sizeof(float4));
    cu.executeKernel(kernel2, args, cu.getNumAtoms(), cu.ThreadBlockSize, cu.ThreadBlockSize*sizeof(float4));
}
