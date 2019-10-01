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

#include "CudaKernels.h"
#include "CudaForceInfo.h"
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
#include "CudaBondedUtilities.h"
#include "CudaExpressionUtilities.h"
#include "CudaIntegrationUtilities.h"
#include "CudaNonbondedUtilities.h"
#include "CudaKernelSources.h"
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

void CudaCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void CudaCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cu.setForcesValid(true);
    cu.setAsCurrent();
    cu.clearAutoclearBuffers();
    for (auto computation : cu.getPreComputations())
        computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    cu.setComputeForceCount(cu.getComputeForceCount()+1);
    nb.prepareInteractions(groups);
    map<string, double>& derivs = cu.getEnergyParamDerivWorkspace();
    for (auto& param : context.getParameters())
        derivs[param.first] = 0;
}

double CudaCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    cu.setAsCurrent();
    cu.getBondedUtilities().computeInteractions(groups);
    cu.getNonbondedUtilities().computeInteractions(groups, includeForces, includeEnergy);
    double sum = 0.0;
    for (auto computation : cu.getPostComputations())
        sum += computation->computeForceAndEnergy(includeForces, includeEnergy, groups);
    cu.getIntegrationUtilities().distributeForcesFromVirtualSites();
    if (includeEnergy)
        sum += cu.reduceEnergy();
    if (!cu.getForcesValid())
        valid = false;
    return sum;
}

void CudaUpdateStateDataKernel::initialize(const System& system) {
}

double CudaUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return cu.getTime();
}

void CudaUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    vector<CudaContext*>& contexts = cu.getPlatformData().contexts;
    for (auto ctx : contexts)
        ctx->setTime(time);
}

void CudaUpdateStateDataKernel::getPositions(ContextImpl& context, vector<Vec3>& positions) {
    cu.setAsCurrent();
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    vector<float4> posCorrection;
    if (cu.getUseDoublePrecision()) {
        double4* posq = (double4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq);
    }
    else if (cu.getUseMixedPrecision()) {
        float4* posq = (float4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq, false);
        posCorrection.resize(numParticles);
        cu.getPosqCorrection().download(posCorrection);
    }
    else {
        float4* posq = (float4*) cu.getPinnedBuffer();
        cu.getPosq().download(posq);
    }
    
    // Filling in the output array is done in parallel for speed.
    
    cu.getPlatformData().threads.execute([&] (ThreadPool& threads, int threadIndex) {
        // Compute the position of each particle to return to the user.  This is done in parallel for speed.
        
        const vector<int>& order = cu.getAtomIndex();
        int numParticles = cu.getNumAtoms();
        Vec3 boxVectors[3];
        cu.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (cu.getUseDoublePrecision()) {
            double4* posq = (double4*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                double4 pos = posq[i];
                int4 offset = cu.getPosCellOffsets()[i];
                positions[order[i]] = Vec3(pos.x, pos.y, pos.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
        else if (cu.getUseMixedPrecision()) {
            float4* posq = (float4*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                float4 pos1 = posq[i];
                float4 pos2 = posCorrection[i];
                int4 offset = cu.getPosCellOffsets()[i];
                positions[order[i]] = Vec3((double)pos1.x+(double)pos2.x, (double)pos1.y+(double)pos2.y, (double)pos1.z+(double)pos2.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
        else {
            float4* posq = (float4*) cu.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                float4 pos = posq[i];
                int4 offset = cu.getPosCellOffsets()[i];
                positions[order[i]] = Vec3(pos.x, pos.y, pos.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
    });
    cu.getPlatformData().threads.waitForThreads();
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
    for (auto& offset : cu.getPosCellOffsets())
        offset = make_int4(0, 0, 0, 0);
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

void CudaUpdateStateDataKernel::getEnergyParameterDerivatives(ContextImpl& context, map<string, double>& derivs) {
    const vector<string>& paramDerivNames = cu.getEnergyParamDerivNames();
    int numDerivs = paramDerivNames.size();
    if (numDerivs == 0)
        return;
    derivs = cu.getEnergyParamDerivWorkspace();
    CudaArray& derivArray = cu.getEnergyParamDerivBuffer();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
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

void CudaUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    cu.getPeriodicBoxVectors(a, b, c);
}

void CudaUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) {
    vector<CudaContext*>& contexts = cu.getPlatformData().contexts;

    // If any particles have been wrapped to the first periodic box, we need to unwrap them
    // to avoid changing their positions.

    vector<Vec3> positions;
    for (auto& offset : cu.getPosCellOffsets()) {
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

void CudaUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    cu.setAsCurrent();
    int version = 2;
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
    Vec3 boxVectors[3];
    cu.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    stream.write((char*) boxVectors, 3*sizeof(Vec3));
    cu.getIntegrationUtilities().createCheckpoint(stream);
    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void CudaUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    cu.setAsCurrent();
    int version;
    stream.read((char*) &version, sizeof(int));
    if (version != 2)
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
    for (auto ctx : contexts) {
        ctx->setTime(time);
        ctx->setStepCount(stepCount);
        ctx->setStepsSinceReorder(stepsSinceReorder);
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
    Vec3 boxVectors[3];
    stream.read((char*) &boxVectors, 3*sizeof(Vec3));
    for (auto ctx : contexts)
        ctx->setPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    cu.getIntegrationUtilities().loadCheckpoint(stream);
    SimTKOpenMMUtilities::loadCheckpoint(stream);
    for (auto listener : cu.getReorderListeners())
        listener->execute();
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

class CudaCalcHarmonicBondForceKernel::ForceInfo : public CudaForceInfo {
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

void CudaCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params.initialize<float2>(cu, numBonds, "bondParams");
    vector<float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = make_float2((float) length, (float) k);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CudaKernelSources::harmonicBondForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params.getDevicePointer(), "float2");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cu.addForce(info);
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
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCalcCustomBondForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomBondForceKernel::~CudaCalcCustomBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
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
    info = new ForceInfo(force);
    cu.addForce(info);

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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customBondGlobals");
        globals.upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals.getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cu.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
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
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
}

double CudaCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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

class CudaCalcHarmonicAngleForceKernel::ForceInfo : public CudaForceInfo {
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

void CudaCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params.initialize<float2>(cu, numAngles, "angleParams");
    vector<float2> paramVector(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angle, k;
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], angle, k);
        paramVector[i] = make_float2((float) angle, (float) k);

    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CudaKernelSources::harmonicAngleForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params.getDevicePointer(), "float2");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::angleForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cu.addForce(info);
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
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCalcCustomAngleForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomAngleForceKernel::~CudaCalcCustomAngleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
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
    info = new ForceInfo(force);
    cu.addForce(info);

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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customAngleGlobals");
        globals.upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals.getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cu.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
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
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::angleForce, replacements), force.getForceGroup());
}

double CudaCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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

class CudaCalcPeriodicTorsionForceKernel::ForceInfo : public CudaForceInfo {
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

void CudaCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params.initialize<float4>(cu, numTorsions, "periodicTorsionParams");
    vector<float4> paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        int periodicity;
        double phase, k;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], periodicity, phase, k);
        paramVector[i] = make_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params.upload(paramVector);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CudaKernelSources::periodicTorsionForce;
    replacements["PARAMS"] = cu.getBondedUtilities().addArgument(params.getDevicePointer(), "float4");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cu.addForce(info);
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
    params.upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCalcRBTorsionForceKernel::ForceInfo : public CudaForceInfo {
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

void CudaCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    cu.setAsCurrent();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params1.initialize<float4>(cu, numTorsions, "rbTorsionParams1");
    params2.initialize<float2>(cu, numTorsions, "rbTorsionParams2");
    vector<float4> paramVector1(numTorsions);
    vector<float2> paramVector2(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], c0, c1, c2, c3, c4, c5);
        paramVector1[i] = make_float4((float) c0, (float) c1, (float) c2, (float) c3);
        paramVector2[i] = make_float2((float) c4, (float) c5);

    }
    params1.upload(paramVector1);
    params2.upload(paramVector2);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = CudaKernelSources::rbTorsionForce;
    replacements["PARAMS1"] = cu.getBondedUtilities().addArgument(params1.getDevicePointer(), "float4");
    replacements["PARAMS2"] = cu.getBondedUtilities().addArgument(params2.getDevicePointer(), "float2");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cu.addForce(info);
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
    params1.upload(paramVector1);
    params2.upload(paramVector2);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCalcCMAPTorsionForceKernel::ForceInfo : public CudaForceInfo {
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
    mapPositionsVec.resize(numMaps);
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
    coefficients.initialize<float4>(cu, coeffVec.size(), "cmapTorsionCoefficients");
    mapPositions.initialize<int2>(cu, numMaps, "cmapTorsionMapPositions");
    torsionMaps.initialize<int>(cu, numTorsions, "cmapTorsionMaps");
    coefficients.upload(coeffVec);
    mapPositions.upload(mapPositionsVec);
    torsionMaps.upload(torsionMapsVec);
    map<string, string> replacements;
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COEFF"] = cu.getBondedUtilities().addArgument(coefficients.getDevicePointer(), "float4");
    replacements["MAP_POS"] = cu.getBondedUtilities().addArgument(mapPositions.getDevicePointer(), "int2");
    replacements["MAPS"] = cu.getBondedUtilities().addArgument(torsionMaps.getDevicePointer(), "int");
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::cmapTorsionForce, replacements), force.getForceGroup());
    info = new ForceInfo(force);
    cu.addForce(info);
}

double CudaCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    int numMaps = force.getNumMaps();
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (mapPositions.getSize() != numMaps)
        throw OpenMMException("updateParametersInContext: The number of maps has changed");
    if (torsionMaps.getSize() != numTorsions)
        throw OpenMMException("updateParametersInContext: The number of CMAP torsions has changed");

    // Update the maps.

    vector<float4> coeffVec;
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
            coeffVec.push_back(make_float4((float) c[j][0], (float) c[j][1], (float) c[j][2], (float) c[j][3]));
            coeffVec.push_back(make_float4((float) c[j][4], (float) c[j][5], (float) c[j][6], (float) c[j][7]));
            coeffVec.push_back(make_float4((float) c[j][8], (float) c[j][9], (float) c[j][10], (float) c[j][11]));
            coeffVec.push_back(make_float4((float) c[j][12], (float) c[j][13], (float) c[j][14], (float) c[j][15]));
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

class CudaCalcCustomTorsionForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomTorsionForceKernel::~CudaCalcCustomTorsionForceKernel() {
    if (params != NULL)
        delete params;
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
    info = new ForceInfo(force);
    cu.addForce(info);

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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customTorsionGlobals");
        globals.upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals.getDevicePointer(), "float");
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = argName+"["+cu.intToString(i)+"]";
            variables[name] = value;
        }
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cu.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        expressions[derivVariable+" += "] = derivExpression;
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
    replacements["APPLY_PERIODIC"] = (force.usesPeriodicBoundaryConditions() ? "1" : "0");
    replacements["COMPUTE_FORCE"] = compute.str();
    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
}

double CudaCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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

class CudaCalcNonbondedForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const NonbondedForce& force) : force(force) {
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
    PmeIO(CudaContext& cu, CUfunction addForcesKernel) : cu(cu), addForcesKernel(addForcesKernel) {
        forceTemp.initialize<float4>(cu, cu.getNumAtoms(), "PmeForce");
    }
    float* getPosq() {
        cu.setAsCurrent();
        cu.getPosq().download(posq);
        return (float*) &posq[0];
    }
    void setForce(float* force) {
        forceTemp.upload(force);
        void* args[] = {&forceTemp.getDevicePointer(), &cu.getForce().getDevicePointer()};
        cu.executeKernel(addForcesKernel, args, cu.getNumAtoms());
    }
private:
    CudaContext& cu;
    vector<float4> posq;
    CudaArray forceTemp;
    CUfunction addForcesKernel;
};

class CudaCalcNonbondedForceKernel::PmePreComputation : public CudaContext::ForcePreComputation {
public:
    PmePreComputation(CudaContext& cu, Kernel& pme, CalcPmeReciprocalForceKernel::IO& io) : cu(cu), pme(pme), io(io) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        Vec3 boxVectors[3] = {Vec3(cu.getPeriodicBoxSize().x, 0, 0), Vec3(0, cu.getPeriodicBoxSize().y, 0), Vec3(0, 0, cu.getPeriodicBoxSize().z)};
        pme.getAs<CalcPmeReciprocalForceKernel>().beginComputation(io, boxVectors, includeEnergy);
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
    SyncStreamPostComputation(CudaContext& cu, CUevent event, CUfunction addEnergyKernel, CudaArray& pmeEnergyBuffer, int forceGroup) : cu(cu), event(event),
            addEnergyKernel(addEnergyKernel), pmeEnergyBuffer(pmeEnergyBuffer), forceGroup(forceGroup) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        if ((groups&(1<<forceGroup)) != 0) {
            cuStreamWaitEvent(cu.getCurrentStream(), event, 0);
            if (includeEnergy) {
                int bufferSize = pmeEnergyBuffer.getSize();
                void* args[] = {&pmeEnergyBuffer.getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(), &bufferSize};
                cu.executeKernel(addEnergyKernel, args, bufferSize);
            }
        }
        return 0.0;
    }
private:
    CudaContext& cu;
    CUevent event;
    CUfunction addEnergyKernel;
    CudaArray& pmeEnergyBuffer;
    int forceGroup;
};

CudaCalcNonbondedForceKernel::~CudaCalcNonbondedForceKernel() {
    cu.setAsCurrent();
    if (sort != NULL)
        delete sort;
    if (fft != NULL)
        delete fft;
    if (dispersionFft != NULL)
        delete dispersionFft;
    if (pmeio != NULL)
        delete pmeio;
    if (hasInitializedFFT) {
        if (useCudaFFT) {
            cufftDestroy(fftForward);
            cufftDestroy(fftBackward);
            if (doLJPME) {
                cufftDestroy(dispersionFftForward);
                cufftDestroy(dispersionFftBackward);                
            }
        }
        if (usePmeStream) {
            cuStreamDestroy(pmeStream);
            cuEventDestroy(pmeSyncEvent);
            cuEventDestroy(paramsSyncEvent);
        }
    }
}

void CudaCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    cu.setAsCurrent();
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "nonbonded"+cu.intToString(forceIndex)+"_";

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
    vector<float4> baseParticleParamVec(cu.getPaddedNumAtoms(), make_float4(0, 0, 0, 0));
    vector<vector<int> > exclusionList(numParticles);
    hasCoulomb = false;
    hasLJ = false;
    for (int i = 0; i < numParticles; i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        baseParticleParamVec[i] = make_float4(charge, sigma, epsilon, 0);
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
    usePosqCharges = hasCoulomb ? cu.requestPosqCharges() : false;

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
    if (force.getUseDispersionCorrection() && cu.getContextIndex() == 0 && !doLJPME)
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
        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        if (cu.getContextIndex() == 0) {
            paramsDefines["INCLUDE_EWALD"] = "1";
            paramsDefines["EWALD_SELF_ENERGY_SCALE"] = cu.doubleToString(ONE_4PI_EPS0*alpha/sqrt(M_PI));
            for (int i = 0; i < numParticles; i++)
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);

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
            cosSinSums.initialize(cu, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), elementSize, "cosSinSums");
        }
    }
    else if (((nonbondedMethod == PME || nonbondedMethod == LJPME) && hasCoulomb) || doLJPME) {
        // Compute the PME parameters.

        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ, false);
        gridSizeX = CudaFFT3D::findLegalDimension(gridSizeX);
        gridSizeY = CudaFFT3D::findLegalDimension(gridSizeY);
        gridSizeZ = CudaFFT3D::findLegalDimension(gridSizeZ);
        if (doLJPME) {
            NonbondedForceImpl::calcPMEParameters(system, force, dispersionAlpha, dispersionGridSizeX,
                                                  dispersionGridSizeY, dispersionGridSizeZ, true);
            dispersionGridSizeX = CudaFFT3D::findLegalDimension(dispersionGridSizeX);
            dispersionGridSizeY = CudaFFT3D::findLegalDimension(dispersionGridSizeY);
            dispersionGridSizeZ = CudaFFT3D::findLegalDimension(dispersionGridSizeZ);
        }

        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
        defines["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
        defines["USE_EWALD"] = "1";
        defines["DO_LJPME"] = doLJPME ? "1" : "0";
        if (doLJPME)
            defines["EWALD_DISPERSION_ALPHA"] = cu.doubleToString(dispersionAlpha);
        if (cu.getContextIndex() == 0) {
            paramsDefines["INCLUDE_EWALD"] = "1";
            paramsDefines["EWALD_SELF_ENERGY_SCALE"] = cu.doubleToString(ONE_4PI_EPS0*alpha/sqrt(M_PI));
            for (int i = 0; i < numParticles; i++)
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);
            if (doLJPME) {
                paramsDefines["INCLUDE_LJPME"] = "1";
                paramsDefines["LJPME_SELF_ENERGY_SCALE"] = cu.doubleToString(pow(dispersionAlpha, 6)/3.0);
                for (int i = 0; i < numParticles; i++)
                    ewaldSelfEnergy += baseParticleParamVec[i].z*pow(baseParticleParamVec[i].y*dispersionAlpha, 6)/3.0;
            }
            char deviceName[100];
            cuDeviceGetName(deviceName, 100, cu.getDevice());
            usePmeStream = (!cu.getPlatformData().disablePmeStream && string(deviceName) != "GeForce GTX 980"); // Using a separate stream is slower on GTX 980
            map<string, string> pmeDefines;
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
            if (usePmeStream)
                pmeDefines["USE_PME_STREAM"] = "1";
            if (cu.getPlatformData().deterministicForces)
                pmeDefines["USE_DETERMINISTIC_FORCES"] = "1";
            map<string, string> replacements;
            replacements["CHARGE"] = (usePosqCharges ? "pos.w" : "charges[atom]");
            CUmodule module = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::pme, replacements), pmeDefines);
            if (cu.getPlatformData().useCpuPme && !doLJPME && usePosqCharges) {
                // Create the CPU PME kernel.

                try {
                    cpuPme = getPlatform().createKernel(CalcPmeReciprocalForceKernel::Name(), *cu.getPlatformData().context);
                    cpuPme.getAs<CalcPmeReciprocalForceKernel>().initialize(gridSizeX, gridSizeY, gridSizeZ, numParticles, alpha, cu.getPlatformData().deterministicForces);
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
                if (doLJPME) {
                    pmeDefines["EWALD_ALPHA"] = cu.doubleToString(dispersionAlpha);
                    pmeDefines["GRID_SIZE_X"] = cu.intToString(dispersionGridSizeX);
                    pmeDefines["GRID_SIZE_Y"] = cu.intToString(dispersionGridSizeY);
                    pmeDefines["GRID_SIZE_Z"] = cu.intToString(dispersionGridSizeZ);
                    pmeDefines["RECIP_EXP_FACTOR"] = cu.doubleToString(M_PI*M_PI/(dispersionAlpha*dispersionAlpha));
                    pmeDefines["USE_LJPME"] = "1";
                    pmeDefines["CHARGE_FROM_SIGEPS"] = "1";
                    double invRCut6 = pow(force.getCutoffDistance(), -6);
                    double dalphaR = dispersionAlpha * force.getCutoffDistance();
                    double dar2 = dalphaR*dalphaR;
                    double dar4 = dar2*dar2;
                    double multShift6 = -invRCut6*(1.0 - exp(-dar2) * (1.0 + dar2 + 0.5*dar4));
                    defines["INVCUT6"] = cu.doubleToString(invRCut6);
                    defines["MULTSHIFT6"] = cu.doubleToString(multShift6);
                    module = cu.createModule(CudaKernelSources::vectorOps+CudaKernelSources::pme, pmeDefines);
                    pmeDispersionFinishSpreadChargeKernel = cu.getKernel(module, "finishSpreadCharge");
                    pmeDispersionGridIndexKernel = cu.getKernel(module, "findAtomGridIndex");
                    pmeDispersionSpreadChargeKernel = cu.getKernel(module, "gridSpreadCharge");
                    pmeDispersionConvolutionKernel = cu.getKernel(module, "reciprocalConvolution");
                    pmeEvalDispersionEnergyKernel = cu.getKernel(module, "gridEvaluateEnergy");
                    pmeInterpolateDispersionForceKernel = cu.getKernel(module, "gridInterpolateForce");
                    cuFuncSetCacheConfig(pmeDispersionSpreadChargeKernel, CU_FUNC_CACHE_PREFER_L1);
                }

                // Create required data structures.

                int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
                int roundedZSize = PmeOrder*(int) ceil(gridSizeZ/(double) PmeOrder);
                int gridElements = gridSizeX*gridSizeY*roundedZSize;
                if (doLJPME) {
                    roundedZSize = PmeOrder*(int) ceil(dispersionGridSizeZ/(double) PmeOrder);
                    gridElements = max(gridElements, dispersionGridSizeX*dispersionGridSizeY*roundedZSize);
                }
                pmeGrid1.initialize(cu, gridElements, 2*elementSize, "pmeGrid1");
                pmeGrid2.initialize(cu, gridElements, 2*elementSize, "pmeGrid2");
                cu.addAutoclearBuffer(pmeGrid2);
                pmeBsplineModuliX.initialize(cu, gridSizeX, elementSize, "pmeBsplineModuliX");
                pmeBsplineModuliY.initialize(cu, gridSizeY, elementSize, "pmeBsplineModuliY");
                pmeBsplineModuliZ.initialize(cu, gridSizeZ, elementSize, "pmeBsplineModuliZ");
                if (doLJPME) {
                    pmeDispersionBsplineModuliX.initialize(cu, dispersionGridSizeX, elementSize, "pmeDispersionBsplineModuliX");
                    pmeDispersionBsplineModuliY.initialize(cu, dispersionGridSizeY, elementSize, "pmeDispersionBsplineModuliY");
                    pmeDispersionBsplineModuliZ.initialize(cu, dispersionGridSizeZ, elementSize, "pmeDispersionBsplineModuliZ");
                }
                pmeAtomGridIndex.initialize<int2>(cu, numParticles, "pmeAtomGridIndex");
                int energyElementSize = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
                pmeEnergyBuffer.initialize(cu, cu.getNumThreadBlocks()*CudaContext::ThreadBlockSize, energyElementSize, "pmeEnergyBuffer");
                cu.clearBuffer(pmeEnergyBuffer);
                sort = new CudaSort(cu, new SortTrait(), cu.getNumAtoms());
                int cufftVersion;
                cufftGetVersion(&cufftVersion);
                useCudaFFT = (cufftVersion >= 7050); // There was a critical bug in version 7.0
                if (useCudaFFT) {
                    cufftResult result = cufftPlan3d(&fftForward, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_D2Z : CUFFT_R2C);
                    if (result != CUFFT_SUCCESS)
                        throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
                    result = cufftPlan3d(&fftBackward, gridSizeX, gridSizeY, gridSizeZ, cu.getUseDoublePrecision() ? CUFFT_Z2D : CUFFT_C2R);
                    if (result != CUFFT_SUCCESS)
                        throw OpenMMException("Error initializing FFT: "+cu.intToString(result));
                    if (doLJPME) {
                        result = cufftPlan3d(&dispersionFftForward, dispersionGridSizeX, dispersionGridSizeY, 
                                                dispersionGridSizeZ, cu.getUseDoublePrecision() ? CUFFT_D2Z : CUFFT_R2C);
                        if (result != CUFFT_SUCCESS)
                            throw OpenMMException("Error initializing disperison FFT: "+cu.intToString(result));
                        result = cufftPlan3d(&dispersionFftBackward, dispersionGridSizeX, dispersionGridSizeY,
                                             dispersionGridSizeZ, cu.getUseDoublePrecision() ? CUFFT_Z2D : CUFFT_C2R);
                        if (result != CUFFT_SUCCESS)
                            throw OpenMMException("Error initializing disperison FFT: "+cu.intToString(result));
                    }
                }
                else {
                    fft = new CudaFFT3D(cu, gridSizeX, gridSizeY, gridSizeZ, true);
                    if (doLJPME)
                        dispersionFft = new CudaFFT3D(cu, dispersionGridSizeX, dispersionGridSizeY, dispersionGridSizeZ, true);
                }

                // Prepare for doing PME on its own stream.

                if (usePmeStream) {
                    cuStreamCreate(&pmeStream, CU_STREAM_NON_BLOCKING);
                    if (useCudaFFT) {
                        cufftSetStream(fftForward, pmeStream);
                        cufftSetStream(fftBackward, pmeStream);
                        if (doLJPME) {
                            cufftSetStream(dispersionFftForward, pmeStream);
                            cufftSetStream(dispersionFftBackward, pmeStream);
                        }
                    }
                    CHECK_RESULT(cuEventCreate(&pmeSyncEvent, CU_EVENT_DISABLE_TIMING), "Error creating event for NonbondedForce");
                    CHECK_RESULT(cuEventCreate(&paramsSyncEvent, CU_EVENT_DISABLE_TIMING), "Error creating event for NonbondedForce");
                    int recipForceGroup = force.getReciprocalSpaceForceGroup();
                    if (recipForceGroup < 0)
                        recipForceGroup = force.getForceGroup();
                    cu.addPreComputation(new SyncStreamPreComputation(cu, pmeStream, pmeSyncEvent, recipForceGroup));
                    cu.addPostComputation(new SyncStreamPostComputation(cu, pmeSyncEvent, cu.getKernel(module, "addEnergy"), pmeEnergyBuffer, recipForceGroup));
                }
                hasInitializedFFT = true;

                // Initialize the b-spline moduli.

                for (int grid = 0; grid < 2; grid++) {
                    int xsize, ysize, zsize;
                    CudaArray *xmoduli, *ymoduli, *zmoduli;
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
        int numContexts = cu.getPlatformData().contexts.size();
        int startIndex = cu.getContextIndex()*force.getNumExceptions()/numContexts;
        int endIndex = (cu.getContextIndex()+1)*force.getNumExceptions()/numContexts;
        int numExclusions = endIndex-startIndex;
        if (numExclusions > 0) {
            paramsDefines["HAS_EXCLUSIONS"] = "1";
            vector<vector<int> > atoms(numExclusions, vector<int>(2));
            exclusionAtoms.initialize<int2>(cu, numExclusions, "exclusionAtoms");
            exclusionParams.initialize<float4>(cu, numExclusions, "exclusionParams");
            vector<int2> exclusionAtomsVec(numExclusions);
            for (int i = 0; i < numExclusions; i++) {
                int j = i+startIndex;
                exclusionAtomsVec[i] = make_int2(exclusions[j].first, exclusions[j].second);
                atoms[i][0] = exclusions[j].first;
                atoms[i][1] = exclusions[j].second;
            }
            exclusionAtoms.upload(exclusionAtomsVec);
            map<string, string> replacements;
            replacements["PARAMS"] = cu.getBondedUtilities().addArgument(exclusionParams.getDevicePointer(), "float4");
            replacements["EWALD_ALPHA"] = cu.doubleToString(alpha);
            replacements["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
            replacements["DO_LJPME"] = doLJPME ? "1" : "0";
            if (doLJPME)
                replacements["EWALD_DISPERSION_ALPHA"] = cu.doubleToString(dispersionAlpha);
            cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::pmeExclusions, replacements), force.getForceGroup());
        }
    }

    // Add the interaction to the default nonbonded kernel.

    string source = cu.replaceStrings(CudaKernelSources::coulombLennardJones, defines);
    charges.initialize(cu, cu.getPaddedNumAtoms(), cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "charges");
    baseParticleParams.initialize<float4>(cu, cu.getPaddedNumAtoms(), "baseParticleParams");
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
        cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"charge", "real", 1, charges.getElementSize(), charges.getDevicePointer()));
    sigmaEpsilon.initialize<float2>(cu, cu.getPaddedNumAtoms(), "sigmaEpsilon");
    if (hasLJ) {
        replacements["SIGMA_EPSILON1"] = prefix+"sigmaEpsilon1";
        replacements["SIGMA_EPSILON2"] = prefix+"sigmaEpsilon2";
        cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"sigmaEpsilon", "float", 2, sizeof(float2), sigmaEpsilon.getDevicePointer()));
    }
    source = cu.replaceStrings(source, replacements);
    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup(), true);

    // Initialize the exceptions.

    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*exceptions.size()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*exceptions.size()/numContexts;
    int numExceptions = endIndex-startIndex;
    if (numExceptions > 0) {
        paramsDefines["HAS_EXCEPTIONS"] = "1";
        exceptionAtoms.resize(numExceptions);
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        exceptionParams.initialize<float4>(cu, numExceptions, "exceptionParams");
        baseExceptionParams.initialize<float4>(cu, numExceptions, "baseExceptionParams");
        vector<float4> baseExceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            baseExceptionParamsVec[i] = make_float4(chargeProd, sigma, epsilon, 0);
            exceptionAtoms[i] = make_pair(atoms[i][0], atoms[i][1]);
        }
        baseExceptionParams.upload(baseExceptionParamsVec);
        map<string, string> replacements;
        replacements["PARAMS"] = cu.getBondedUtilities().addArgument(exceptionParams.getDevicePointer(), "float4");
        cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::nonbondedExceptions, replacements), force.getForceGroup());
    }
    
    // Initialize parameter offsets.

    vector<vector<float4> > particleOffsetVec(force.getNumParticles());
    vector<vector<float4> > exceptionOffsetVec(force.getNumExceptions());
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
        particleOffsetVec[particle].push_back(make_float4(charge, sigma, epsilon, paramIndex));
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
        exceptionOffsetVec[exceptionIndex[exception]].push_back(make_float4(charge, sigma, epsilon, paramIndex));
    }
    paramValues.resize(paramNames.size(), 0.0);
    particleParamOffsets.initialize<float4>(cu, max(force.getNumParticleParameterOffsets(), 1), "particleParamOffsets");
    exceptionParamOffsets.initialize<float4>(cu, max(force.getNumExceptionParameterOffsets(), 1), "exceptionParamOffsets");
    particleOffsetIndices.initialize<int>(cu, cu.getPaddedNumAtoms()+1, "particleOffsetIndices");
    exceptionOffsetIndices.initialize<int>(cu, force.getNumExceptions()+1, "exceptionOffsetIndices");
    vector<int> particleOffsetIndicesVec, exceptionOffsetIndicesVec;
    vector<float4> p, e;
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
    globalParams.initialize(cu, max((int) paramValues.size(), 1), cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "globalParams");
    recomputeParams = true;
    
    // Initialize the kernel for updating parameters.
    
    CUmodule module = cu.createModule(CudaKernelSources::nonbondedParameters, paramsDefines);
    computeParamsKernel = cu.getKernel(module, "computeParameters");
    computeExclusionParamsKernel = cu.getKernel(module, "computeExclusionParameters");
    info = new ForceInfo(force);
    cu.addForce(info);
}

double CudaCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
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
        bool computeSelfEnergy = (includeEnergy && includeReciprocal);
        int numAtoms = cu.getPaddedNumAtoms();
        vector<void*> paramsArgs = {&cu.getEnergyBuffer().getDevicePointer(), &computeSelfEnergy, &globalParams.getDevicePointer(), &numAtoms,
                &baseParticleParams.getDevicePointer(), &cu.getPosq().getDevicePointer(), &charges.getDevicePointer(), &sigmaEpsilon.getDevicePointer(),
                &particleParamOffsets.getDevicePointer(), &particleOffsetIndices.getDevicePointer()};
        int numExceptions;
        if (exceptionParams.isInitialized()) {
            numExceptions = exceptionParams.getSize();
            paramsArgs.push_back(&numExceptions);
            paramsArgs.push_back(&baseExceptionParams.getDevicePointer());
            paramsArgs.push_back(&exceptionParams.getDevicePointer());
            paramsArgs.push_back(&exceptionParamOffsets.getDevicePointer());
            paramsArgs.push_back(&exceptionOffsetIndices.getDevicePointer());
        }
        cu.executeKernel(computeParamsKernel, &paramsArgs[0], cu.getPaddedNumAtoms());
        if (exclusionParams.isInitialized()) {
            int numExclusions = exclusionParams.getSize();
            vector<void*> exclusionParamsArgs = {&cu.getPosq().getDevicePointer(), &charges.getDevicePointer(), &sigmaEpsilon.getDevicePointer(),
                    &numExclusions, &exclusionAtoms.getDevicePointer(), &exclusionParams.getDevicePointer()};
            cu.executeKernel(computeExclusionParamsKernel, &exclusionParamsArgs[0], numExclusions);
        }
        if (usePmeStream) {
            cuEventRecord(paramsSyncEvent, cu.getCurrentStream());
            cuStreamWaitEvent(pmeStream, paramsSyncEvent, 0);
        }
        if (hasOffsets)
            energy = 0.0; // The Ewald self energy was computed in the kernel.
        recomputeParams = false;
    }
    
    // Do reciprocal space calculations.
    
    if (cosSinSums.isInitialized() && includeReciprocal) {
        void* sumsArgs[] = {&cu.getEnergyBuffer().getDevicePointer(), &cu.getPosq().getDevicePointer(), &cosSinSums.getDevicePointer(), cu.getPeriodicBoxSizePointer()};
        cu.executeKernel(ewaldSumsKernel, sumsArgs, cosSinSums.getSize());
        void* forcesArgs[] = {&cu.getForce().getDevicePointer(), &cu.getPosq().getDevicePointer(), &cosSinSums.getDevicePointer(), cu.getPeriodicBoxSizePointer()};
        cu.executeKernel(ewaldForcesKernel, forcesArgs, cu.getNumAtoms());
    }
    if (pmeGrid1.isInitialized() && includeReciprocal) {
        if (usePmeStream)
            cu.setCurrentStream(pmeStream);

        // Invert the periodic box vectors.

        Vec3 boxVectors[3];
        cu.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double determinant = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
        double scale = 1.0/determinant;
        double3 recipBoxVectors[3];
        recipBoxVectors[0] = make_double3(boxVectors[1][1]*boxVectors[2][2]*scale, 0, 0);
        recipBoxVectors[1] = make_double3(-boxVectors[1][0]*boxVectors[2][2]*scale, boxVectors[0][0]*boxVectors[2][2]*scale, 0);
        recipBoxVectors[2] = make_double3((boxVectors[1][0]*boxVectors[2][1]-boxVectors[1][1]*boxVectors[2][0])*scale, -boxVectors[0][0]*boxVectors[2][1]*scale, boxVectors[0][0]*boxVectors[1][1]*scale);
        float3 recipBoxVectorsFloat[3];
        void* recipBoxVectorPointer[3];
        if (cu.getUseDoublePrecision()) {
            recipBoxVectorPointer[0] = &recipBoxVectors[0];
            recipBoxVectorPointer[1] = &recipBoxVectors[1];
            recipBoxVectorPointer[2] = &recipBoxVectors[2];
        }
        else {
            recipBoxVectorsFloat[0] = make_float3((float) recipBoxVectors[0].x, 0, 0);
            recipBoxVectorsFloat[1] = make_float3((float) recipBoxVectors[1].x, (float) recipBoxVectors[1].y, 0);
            recipBoxVectorsFloat[2] = make_float3((float) recipBoxVectors[2].x, (float) recipBoxVectors[2].y, (float) recipBoxVectors[2].z);
            recipBoxVectorPointer[0] = &recipBoxVectorsFloat[0];
            recipBoxVectorPointer[1] = &recipBoxVectorsFloat[1];
            recipBoxVectorPointer[2] = &recipBoxVectorsFloat[2];
        }

        // Execute the reciprocal space kernels.

        if (hasCoulomb) {
            void* gridIndexArgs[] = {&cu.getPosq().getDevicePointer(), &pmeAtomGridIndex.getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                    cu.getInvPeriodicBoxSizePointer(), cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
                    recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2]};
            cu.executeKernel(pmeGridIndexKernel, gridIndexArgs, cu.getNumAtoms());

            sort->sort(pmeAtomGridIndex);

            void* spreadArgs[] = {&cu.getPosq().getDevicePointer(), &pmeGrid2.getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                    cu.getInvPeriodicBoxSizePointer(), cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
                    recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2], &pmeAtomGridIndex.getDevicePointer(),
                    &charges.getDevicePointer()};
            cu.executeKernel(pmeSpreadChargeKernel, spreadArgs, cu.getNumAtoms(), 128);

            void* finishSpreadArgs[] = {&pmeGrid2.getDevicePointer(), &pmeGrid1.getDevicePointer()};
            cu.executeKernel(pmeFinishSpreadChargeKernel, finishSpreadArgs, gridSizeX*gridSizeY*gridSizeZ, 256);

            if (useCudaFFT) {
                if (cu.getUseDoublePrecision())
                    cufftExecD2Z(fftForward, (double*) pmeGrid1.getDevicePointer(), (double2*) pmeGrid2.getDevicePointer());
                else
                    cufftExecR2C(fftForward, (float*) pmeGrid1.getDevicePointer(), (float2*) pmeGrid2.getDevicePointer());
            }
            else {
                fft->execFFT(pmeGrid1, pmeGrid2, true);
            }

            if (includeEnergy) {
                void* computeEnergyArgs[] = {&pmeGrid2.getDevicePointer(), usePmeStream ? &pmeEnergyBuffer.getDevicePointer() : &cu.getEnergyBuffer().getDevicePointer(),
                        &pmeBsplineModuliX.getDevicePointer(), &pmeBsplineModuliY.getDevicePointer(), &pmeBsplineModuliZ.getDevicePointer(),
                        cu.getPeriodicBoxSizePointer(), recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2]};
                cu.executeKernel(pmeEvalEnergyKernel, computeEnergyArgs, gridSizeX*gridSizeY*gridSizeZ);
            }

            void* convolutionArgs[] = {&pmeGrid2.getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(),
                    &pmeBsplineModuliX.getDevicePointer(), &pmeBsplineModuliY.getDevicePointer(), &pmeBsplineModuliZ.getDevicePointer(),
                    cu.getPeriodicBoxSizePointer(), recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2]};
            cu.executeKernel(pmeConvolutionKernel, convolutionArgs, gridSizeX*gridSizeY*gridSizeZ, 256);

            if (useCudaFFT) {
                if (cu.getUseDoublePrecision())
                    cufftExecZ2D(fftBackward, (double2*) pmeGrid2.getDevicePointer(), (double*) pmeGrid1.getDevicePointer());
                else
                    cufftExecC2R(fftBackward, (float2*) pmeGrid2.getDevicePointer(), (float*)  pmeGrid1.getDevicePointer());
            }
            else {
                fft->execFFT(pmeGrid2, pmeGrid1, false);
            }

            void* interpolateArgs[] = {&cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(), &pmeGrid1.getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                    cu.getInvPeriodicBoxSizePointer(), cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
                    recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2], &pmeAtomGridIndex.getDevicePointer(),
                    &charges.getDevicePointer()};
            cu.executeKernel(pmeInterpolateForceKernel, interpolateArgs, cu.getNumAtoms(), 128);
        }

        if (doLJPME && hasLJ) {
            if (!hasCoulomb) {
                void* gridIndexArgs[] = {&cu.getPosq().getDevicePointer(), &pmeAtomGridIndex.getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                        cu.getInvPeriodicBoxSizePointer(), cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
                        recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2]};
                cu.executeKernel(pmeDispersionGridIndexKernel, gridIndexArgs, cu.getNumAtoms());

                sort->sort(pmeAtomGridIndex);
                cu.clearBuffer(pmeEnergyBuffer);
            }

            cu.clearBuffer(pmeGrid2);
            void* spreadArgs[] = {&cu.getPosq().getDevicePointer(), &pmeGrid2.getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                    cu.getInvPeriodicBoxSizePointer(), cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
                    recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2], &pmeAtomGridIndex.getDevicePointer(),
                    &sigmaEpsilon.getDevicePointer()};
            cu.executeKernel(pmeDispersionSpreadChargeKernel, spreadArgs, cu.getNumAtoms(), 128);

            void* finishSpreadArgs[] = {&pmeGrid2.getDevicePointer(), &pmeGrid1.getDevicePointer()};
            cu.executeKernel(pmeDispersionFinishSpreadChargeKernel, finishSpreadArgs, dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ, 256);

            if (useCudaFFT) {
                if (cu.getUseDoublePrecision())
                    cufftExecD2Z(dispersionFftForward, (double*) pmeGrid1.getDevicePointer(), (double2*) pmeGrid2.getDevicePointer());
                else
                    cufftExecR2C(dispersionFftForward, (float*) pmeGrid1.getDevicePointer(), (float2*) pmeGrid2.getDevicePointer());
            }
            else {
                dispersionFft->execFFT(pmeGrid1, pmeGrid2, true);
            }

            if (includeEnergy) {
                void* computeEnergyArgs[] = {&pmeGrid2.getDevicePointer(), usePmeStream ? &pmeEnergyBuffer.getDevicePointer() : &cu.getEnergyBuffer().getDevicePointer(),
                        &pmeDispersionBsplineModuliX.getDevicePointer(), &pmeDispersionBsplineModuliY.getDevicePointer(), &pmeDispersionBsplineModuliZ.getDevicePointer(),
                        cu.getPeriodicBoxSizePointer(), recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2]};
                cu.executeKernel(pmeEvalDispersionEnergyKernel, computeEnergyArgs, dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ);
            }

            void* convolutionArgs[] = {&pmeGrid2.getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(),
                    &pmeDispersionBsplineModuliX.getDevicePointer(), &pmeDispersionBsplineModuliY.getDevicePointer(), &pmeDispersionBsplineModuliZ.getDevicePointer(),
                    cu.getPeriodicBoxSizePointer(), recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2]};
            cu.executeKernel(pmeDispersionConvolutionKernel, convolutionArgs, dispersionGridSizeX*dispersionGridSizeY*dispersionGridSizeZ, 256);

            if (useCudaFFT) {
                if (cu.getUseDoublePrecision())
                    cufftExecZ2D(dispersionFftBackward, (double2*) pmeGrid2.getDevicePointer(), (double*) pmeGrid1.getDevicePointer());
                else
                    cufftExecC2R(dispersionFftBackward, (float2*) pmeGrid2.getDevicePointer(), (float*)  pmeGrid1.getDevicePointer());
            }
            else {
                dispersionFft->execFFT(pmeGrid2, pmeGrid1, false);
            }

            void* interpolateArgs[] = {&cu.getPosq().getDevicePointer(), &cu.getForce().getDevicePointer(), &pmeGrid1.getDevicePointer(), cu.getPeriodicBoxSizePointer(),
                    cu.getInvPeriodicBoxSizePointer(), cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
                    recipBoxVectorPointer[0], recipBoxVectorPointer[1], recipBoxVectorPointer[2], &pmeAtomGridIndex.getDevicePointer(),
                    &sigmaEpsilon.getDevicePointer()};
            cu.executeKernel(pmeInterpolateDispersionForceKernel, interpolateArgs, cu.getNumAtoms(), 128);
        }
        if (usePmeStream) {
            cuEventRecord(pmeSyncEvent, pmeStream);
            cu.restoreDefaultStream();
        }
    }

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
    
    vector<float4> baseParticleParamVec(cu.getPaddedNumAtoms(), make_float4(0, 0, 0, 0));
    const vector<int>& order = cu.getAtomIndex();
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, sigma, epsilon;
        force.getParticleParameters(i, charge, sigma, epsilon);
        baseParticleParamVec[i] = make_float4(charge, sigma, epsilon, 0);
    }
    baseParticleParams.upload(baseParticleParamVec);
    
    // Record the exceptions.
    
    if (numExceptions > 0) {
        vector<vector<int> > atoms(numExceptions, vector<int>(2));
        vector<float4> baseExceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double chargeProd, sigma, epsilon;
            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
            baseExceptionParamsVec[i] = make_float4(chargeProd, sigma, epsilon, 0);
        }
        baseExceptionParams.upload(baseExceptionParamsVec);
    }
    
    // Compute other values.
    
    ewaldSelfEnergy = 0.0;
    if (nonbondedMethod == Ewald || nonbondedMethod == PME || nonbondedMethod == LJPME) {
        if (cu.getContextIndex() == 0) {
            for (int i = 0; i < force.getNumParticles(); i++) {
                ewaldSelfEnergy -= baseParticleParamVec[i].x*baseParticleParamVec[i].x*ONE_4PI_EPS0*alpha/sqrt(M_PI);
                if (doLJPME)
                    ewaldSelfEnergy += baseParticleParamVec[i].z*pow(baseParticleParamVec[i].y*dispersionAlpha, 6)/3.0;
            }
        }
    }
    if (force.getUseDispersionCorrection() && cu.getContextIndex() == 0 && (nonbondedMethod == CutoffPeriodic || nonbondedMethod == Ewald || nonbondedMethod == PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
    cu.invalidateMolecules();
    recomputeParams = true;
}

void CudaCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (cu.getPlatformData().useCpuPme)
        cpuPme.getAs<CalcPmeReciprocalForceKernel>().getPMEParameters(alpha, nx, ny, nz);
    else {
        alpha = this->alpha;
        nx = gridSizeX;
        ny = gridSizeY;
        nz = gridSizeZ;
    }
}

void CudaCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (!doLJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    if (cu.getPlatformData().useCpuPme)
        //cpuPme.getAs<CalcPmeReciprocalForceKernel>().getLJPMEParameters(alpha, nx, ny, nz);
        throw OpenMMException("getPMEParametersInContext: CPUPME has not been implemented for LJPME yet.");
    else {
        alpha = this->dispersionAlpha;
        nx = dispersionGridSizeX;
        ny = dispersionGridSizeY;
        nz = dispersionGridSizeZ;
    }
}

class CudaCalcCustomNonbondedForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomNonbondedForce& force) : force(force) {
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

CudaCalcCustomNonbondedForceKernel::~CudaCalcCustomNonbondedForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customNonbondedGlobals");
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
    vector<string> tableTypes;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName, "float", width, width*sizeof(float), tabulatedFunctions[i].getDevicePointer()));
        if (width == 1)
            tableTypes.push_back("float");
        else
            tableTypes.push_back("float"+cu.intToString(width));
    }

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
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
        string value = "globals["+cu.intToString(i)+"]";
        variables.push_back(makeVariable(name, prefix+value));
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cu.getNonbondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        forceExpressions[derivVariable+" += interactionScale*switchValue*"] = derivExpression;
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
        initInteractionGroups(force, source, tableTypes);
    else {
        cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup(), true);
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
            cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"params"+cu.intToString(i+1), buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
        }
        if (globals.isInitialized()) {
            globals.upload(globalParamValues);
            cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(float), globals.getDevicePointer()));
        }
    }
    info = new ForceInfo(force);
    cu.addForce(info);
    
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

void CudaCalcCustomNonbondedForceKernel::initInteractionGroups(const CustomNonbondedForce& force, const string& interactionSource, const vector<string>& tableTypes) {
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
    interactionGroupData.initialize<int4>(cu, groupData.size(), "interactionGroupData");
    interactionGroupData.upload(groupData);
    numGroupTiles.initialize<int>(cu, 1, "numGroupTiles");

    // Allocate space for a neighbor list, if necessary.

    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && groupData.size() > cu.getNumThreadBlocks()) {
        filteredGroupData.initialize<int4>(cu, groupData.size(), "filteredGroupData");
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
    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
        args << ", const " << tableTypes[i]<< "* __restrict__ table" << i;
    if (globals.isInitialized())
        args<<", const float* __restrict__ globals";
    if (hasParamDerivs)
        args << ", mixed* __restrict__ energyParamDerivs";
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
    stringstream initDerivs, saveDerivs;
    const vector<string>& allParamDerivNames = cu.getEnergyParamDerivNames();
    int numDerivs = allParamDerivNames.size();
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cu.getNonbondedUtilities().addEnergyParameterDerivative(paramName);
        initDerivs<<"mixed "<<derivVariable<<" = 0;\n";
        for (int index = 0; index < numDerivs; index++)
            if (allParamDerivNames[index] == paramName)
                saveDerivs<<"energyParamDerivs[(blockIdx.x*blockDim.x+threadIdx.x)*"<<numDerivs<<"+"<<index<<"] += "<<derivVariable<<";\n";
    }
    replacements["INIT_DERIVATIVES"] = initDerivs.str();
    replacements["SAVE_DERIVATIVES"] = saveDerivs.str();
    map<string, string> defines;
    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
        defines["USE_CUTOFF"] = "1";
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic)
        defines["USE_PERIODIC"] = "1";
    int localMemorySize = max(32, cu.getNonbondedUtilities().getForceThreadBlockSize());
    defines["LOCAL_MEMORY_SIZE"] = cu.intToString(localMemorySize);
    defines["WARPS_IN_BLOCK"] = cu.intToString(localMemorySize/32);
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
    double paddedCutoff = cu.getNonbondedUtilities().padCutoff(cutoff);
    defines["PADDED_CUTOFF_SQUARED"] = cu.doubleToString(paddedCutoff*paddedCutoff);
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    defines["TILE_SIZE"] = "32";
    defines["NUM_TILES"] = cu.intToString(numTileSets);
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*numTileSets/numContexts;
    int endIndex = (cu.getContextIndex()+1)*numTileSets/numContexts;
    defines["FIRST_TILE"] = cu.intToString(startIndex);
    defines["LAST_TILE"] = cu.intToString(endIndex);
    if ((localDataSize/4)%2 == 0 && !cu.getUseDoublePrecision())
        defines["PARAMETER_SIZE_IS_EVEN"] = "1";
    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::customNonbondedGroups, replacements), defines);
    interactionGroupKernel = cu.getKernel(module, "computeInteractionGroups");
    prepareNeighborListKernel = cu.getKernel(module, "prepareToBuildNeighborList");
    buildNeighborListKernel = cu.getKernel(module, "buildNeighborList");
    numGroupThreadBlocks = cu.getNonbondedUtilities().getNumForceThreadBlocks();
}

double CudaCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    useNeighborList = (filteredGroupData.isInitialized() && cu.getNonbondedUtilities().getUseCutoff());
    if (useNeighborList && cu.getContextIndex() > 0) {
        // When using a neighbor list, run the whole calculation on a single device.
        return 0.0;
    }
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
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
            interactionGroupArgs.push_back(&cu.getForce().getDevicePointer());
            interactionGroupArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
            interactionGroupArgs.push_back(&cu.getPosq().getDevicePointer());
            interactionGroupArgs.push_back(&(useNeighborList ? filteredGroupData : interactionGroupData).getDevicePointer());
            interactionGroupArgs.push_back(&numGroupTiles.getDevicePointer());
            interactionGroupArgs.push_back(&useNeighborList);
            interactionGroupArgs.push_back(cu.getPeriodicBoxSizePointer());
            interactionGroupArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            interactionGroupArgs.push_back(cu.getPeriodicBoxVecXPointer());
            interactionGroupArgs.push_back(cu.getPeriodicBoxVecYPointer());
            interactionGroupArgs.push_back(cu.getPeriodicBoxVecZPointer());
            for (auto& buffer : params->getBuffers())
                interactionGroupArgs.push_back(&buffer.getMemory());
            for (auto& function : tabulatedFunctions)
                interactionGroupArgs.push_back(&function.getDevicePointer());
            if (globals.isInitialized())
                interactionGroupArgs.push_back(&globals.getDevicePointer());
            if (hasParamDerivs)
                interactionGroupArgs.push_back(&cu.getEnergyParamDerivBuffer().getDevicePointer());
            if (useNeighborList) {
                // Initialize kernels for building the interaction group neighbor list.

                prepareNeighborListArgs.push_back(&cu.getNonbondedUtilities().getRebuildNeighborList().getDevicePointer());
                prepareNeighborListArgs.push_back(&numGroupTiles.getDevicePointer());
                buildNeighborListArgs.push_back(&cu.getNonbondedUtilities().getRebuildNeighborList().getDevicePointer());
                buildNeighborListArgs.push_back(&numGroupTiles.getDevicePointer());
                buildNeighborListArgs.push_back(&cu.getPosq().getDevicePointer());
                buildNeighborListArgs.push_back(&interactionGroupData.getDevicePointer());
                buildNeighborListArgs.push_back(&filteredGroupData.getDevicePointer());
                buildNeighborListArgs.push_back(cu.getPeriodicBoxSizePointer());
                buildNeighborListArgs.push_back(cu.getInvPeriodicBoxSizePointer());
                buildNeighborListArgs.push_back(cu.getPeriodicBoxVecXPointer());
                buildNeighborListArgs.push_back(cu.getPeriodicBoxVecYPointer());
                buildNeighborListArgs.push_back(cu.getPeriodicBoxVecZPointer());
            }
        }
        int forceThreadBlockSize = cu.getNonbondedUtilities().getForceThreadBlockSize();
        if (useNeighborList) {
            // Rebuild the neighbor list, if necessary.

            cu.executeKernel(prepareNeighborListKernel, &prepareNeighborListArgs[0], 1, 1);
            cu.executeKernel(buildNeighborListKernel, &buildNeighborListArgs[0], numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
        }
        cu.executeKernel(interactionGroupKernel, &interactionGroupArgs[0], numGroupThreadBlocks*forceThreadBlockSize, forceThreadBlockSize);
    }
    double4 boxSize = cu.getPeriodicBoxSize();
    double volume = boxSize.x*boxSize.y*boxSize.z;
    map<string, double>& derivs = cu.getEnergyParamDerivWorkspace();
    for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
        derivs[forceCopy->getEnergyParameterDerivativeName(i)] += longRangeCoefficientDerivs[i]/volume;
    return longRangeCoefficient/volume;
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
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCalcGBSAOBCForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const GBSAOBCForce& force) : force(force) {
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

void CudaCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    cu.setAsCurrent();
    if (cu.getPlatformData().contexts.size() > 1)
        throw OpenMMException("GBSAOBCForce does not support using multiple CUDA devices");
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "obc"+cu.intToString(forceIndex)+"_";
    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
    params.initialize<float2>(cu, cu.getPaddedNumAtoms(), "gbsaObcParams");
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    charges.initialize(cu, cu.getPaddedNumAtoms(), elementSize, "gbsaObcCharges");
    bornRadii.initialize(cu, cu.getPaddedNumAtoms(), elementSize, "bornRadii");
    obcChain.initialize(cu, cu.getPaddedNumAtoms(), elementSize, "obcChain");
    bornSum.initialize<long long>(cu, cu.getPaddedNumAtoms(), "bornSum");
    bornForce.initialize<long long>(cu, cu.getPaddedNumAtoms(), "bornForce");
    cu.addAutoclearBuffer(bornSum);
    cu.addAutoclearBuffer(bornForce);
    CudaArray& posq = cu.getPosq();
    vector<double> chargeVec(cu.getPaddedNumAtoms());
    vector<float2> paramsVector(cu.getPaddedNumAtoms(), make_float2(1, 1));
    const double dielectricOffset = 0.009;
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        chargeVec[i] = charge;
        paramsVector[i] = make_float2((float) radius, (float) (scalingFactor*radius));
    }
    charges.upload(chargeVec, true);
    params.upload(paramsVector);
    prefactor = -ONE_4PI_EPS0*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
    surfaceAreaFactor = -6.0*4*M_PI*force.getSurfaceAreaEnergy();
    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
    cutoff = force.getCutoffDistance();
    string source = CudaKernelSources::gbsaObc2;
    map<string, string> replacements;
    replacements["CHARGE1"] = prefix+"charge1";
    replacements["CHARGE2"] = prefix+"charge2";
    replacements["OBC_PARAMS1"] = prefix+"obcParams1";
    replacements["OBC_PARAMS2"] = prefix+"obcParams2";
    replacements["BORN_FORCE1"] = prefix+"bornForce1";
    replacements["BORN_FORCE2"] = prefix+"bornForce2";
    source = cu.replaceStrings(source, replacements);
    nb.addInteraction(useCutoff, usePeriodic, false, cutoff, vector<vector<int> >(), source, force.getForceGroup());
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"charge", "float", 1, sizeof(float), charges.getDevicePointer()));;
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"obcParams", "float", 2, sizeof(float2), params.getDevicePointer()));
    nb.addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"bornForce", "long long", 1, sizeof(long long), bornForce.getDevicePointer()));
    info = new ForceInfo(force);
    cu.addForce(info);
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
        defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
        defines["CUTOFF"] = cu.doubleToString(cutoff);
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
        computeSumArgs.push_back(&bornSum.getDevicePointer());
        computeSumArgs.push_back(&cu.getPosq().getDevicePointer());
        computeSumArgs.push_back(&charges.getDevicePointer());
        computeSumArgs.push_back(&params.getDevicePointer());
        if (nb.getUseCutoff()) {
            computeSumArgs.push_back(&nb.getInteractingTiles().getDevicePointer());
            computeSumArgs.push_back(&nb.getInteractionCount().getDevicePointer());
            computeSumArgs.push_back(cu.getPeriodicBoxSizePointer());
            computeSumArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            computeSumArgs.push_back(cu.getPeriodicBoxVecXPointer());
            computeSumArgs.push_back(cu.getPeriodicBoxVecYPointer());
            computeSumArgs.push_back(cu.getPeriodicBoxVecZPointer());
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
        force1Args.push_back(&bornForce.getDevicePointer());
        force1Args.push_back(&cu.getEnergyBuffer().getDevicePointer());
        force1Args.push_back(&cu.getPosq().getDevicePointer());
        force1Args.push_back(&charges.getDevicePointer());
        force1Args.push_back(&bornRadii.getDevicePointer());
        force1Args.push_back(NULL);
        if (nb.getUseCutoff()) {
            force1Args.push_back(&nb.getInteractingTiles().getDevicePointer());
            force1Args.push_back(&nb.getInteractionCount().getDevicePointer());
            force1Args.push_back(cu.getPeriodicBoxSizePointer());
            force1Args.push_back(cu.getInvPeriodicBoxSizePointer());
            force1Args.push_back(cu.getPeriodicBoxVecXPointer());
            force1Args.push_back(cu.getPeriodicBoxVecYPointer());
            force1Args.push_back(cu.getPeriodicBoxVecZPointer());
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
    force1Args[6] = &includeEnergy;
    if (nb.getUseCutoff()) {
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            computeSumArgs[4] = &nb.getInteractingTiles().getDevicePointer();
            force1Args[7] = &nb.getInteractingTiles().getDevicePointer();
            computeSumArgs[14] = &nb.getInteractingAtoms().getDevicePointer();
            force1Args[17] = &nb.getInteractingAtoms().getDevicePointer();
        }
    }
    cu.executeKernel(computeBornSumKernel, &computeSumArgs[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    float alpha = 1.0f, beta = 0.8f, gamma = 4.85f;
    void* reduceSumArgs[] = {&alpha, &beta, &gamma, &bornSum.getDevicePointer(), &params.getDevicePointer(),
            &bornRadii.getDevicePointer(), &obcChain.getDevicePointer()};
    cu.executeKernel(reduceBornSumKernel, reduceSumArgs, cu.getPaddedNumAtoms());
    cu.executeKernel(force1Kernel, &force1Args[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    void* reduceForceArgs[] = {&bornForce.getDevicePointer(), &cu.getEnergyBuffer().getDevicePointer(), &params.getDevicePointer(),
            &bornRadii.getDevicePointer(), &obcChain.getDevicePointer()};
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
    
    vector<double> chargeVector(cu.getPaddedNumAtoms(), 0.0);
    vector<float2> paramsVector(cu.getPaddedNumAtoms());
    const double dielectricOffset = 0.009;
    for (int i = 0; i < numParticles; i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        chargeVector[i] = charge;
        radius -= dielectricOffset;
        paramsVector[i] = make_float2((float) radius, (float) (scalingFactor*radius));
    }
    for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
        paramsVector[i] = make_float2(1, 1);
    charges.upload(chargeVector, true);
    params.upload(paramsVector);
    
    // Mark that the current reordering may be invalid.
    
    cu.invalidateMolecules();
}

class CudaCalcCustomGBForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomGBForce& force) : force(force) {
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
    for (auto d : dValuedParam)
        delete d;
}

void CudaCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    cu.setAsCurrent();
    if (cu.getPlatformData().contexts.size() > 1)
        throw OpenMMException("CustomGBForce does not support using multiple CUDA devices");
    cutoff = force.getCutoffDistance();
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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customGBGlobals");
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
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName, "float", width, width*sizeof(float), tabulatedFunctions[i].getDevicePointer()));
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
    longEnergyDerivs.initialize<long long>(cu, force.getNumComputedValues()*cu.getPaddedNumAtoms(), "customGBLongEnergyDerivatives");
    energyDerivs = new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms(), "customGBEnergyDerivatives", true);
    energyDerivChain = new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms(), "customGBEnergyDerivativeChain", true);
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    needEnergyParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    dValue0dParam.resize(force.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        dValuedParam.push_back(new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms(), "dValuedParam", true, cu.getUseDoublePrecision()));
        dValue0dParam[i].initialize<long long>(cu, cu.getPaddedNumAtoms(), "dValue0dParam");
        cu.addAutoclearBuffer(dValue0dParam[i]);
        string name = force.getEnergyParameterDerivativeName(i);
        cu.addEnergyParameterDerivative(name);
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
            string value = "globals["+cu.intToString(i)+"]";
            variables.push_back(makeVariable(name, value));
        }
        map<string, Lepton::ParsedExpression> n2ValueExpressions;
        stringstream n2ValueSource;
        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[0], functions).optimize();
        n2ValueExpressions["tempValue1 = "] = ex;
        n2ValueExpressions["tempValue2 = "] = ex.renameVariables(rename);
        for (int i = 0; i < valueParamDerivExpressions[0].size(); i++) {
            string variableBase = "temp_dValue0dParam"+cu.intToString(i+1);
            if (!isZeroExpression(valueParamDerivExpressions[0][i])) {
                n2ValueExpressions[variableBase+"_1 = "] = valueParamDerivExpressions[0][i];
                n2ValueExpressions[variableBase+"_2 = "] = valueParamDerivExpressions[0][i].renameVariables(rename);
            }
        }
        n2ValueSource << cu.getExpressionUtilities().createExpressions(n2ValueExpressions, variables, functionList, functionDefinitions, "temp");
        map<string, string> replacements;
        string n2ValueStr = n2ValueSource.str();
        replacements["COMPUTE_VALUE"] = n2ValueStr;
        stringstream extraArgs, atomParams, loadLocal1, loadLocal2, load1, load2, tempDerivs1, tempDerivs2, storeDeriv1, storeDeriv2;
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
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string derivName = "dValue0dParam"+cu.intToString(i+1);
            extraArgs << ", unsigned long long* __restrict__ global_" << derivName;
            atomParams << "real " << derivName << ";\n";
            loadLocal2 << "localData[localAtomIndex]." << derivName << " = 0;\n";
            load1 << "real " << derivName << " = 0;\n";
            if (!isZeroExpression(valueParamDerivExpressions[0][i])) {
                load2 << "real temp_" << derivName << "_1 = 0;\n";
                load2 << "real temp_" << derivName << "_2 = 0;\n";
                tempDerivs1 << derivName << " += temp_" << derivName << "_1;\n";
                tempDerivs2 << "localData[tbx+tj]." << derivName << " += temp_" << derivName << "_2;\n";
                storeDeriv1 << "atomicAdd(&global_" << derivName << "[offset1], static_cast<unsigned long long>((long long) (" << derivName << "*0x100000000)));\n";
                storeDeriv2 << "atomicAdd(&global_" << derivName << "[offset2], static_cast<unsigned long long>((long long) (localData[threadIdx.x]." << derivName << "*0x100000000)));\n";
            }
        }
        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
        replacements["ATOM_PARAMETER_DATA"] = atomParams.str();
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
        if (atomParamSize%2 == 0 && !cu.getUseDoublePrecision())
            pairValueDefines["NEED_PADDING"] = "1";
        pairValueDefines["WARPS_PER_GROUP"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize()/CudaContext::TileSize);
        pairValueDefines["THREAD_BLOCK_SIZE"] = cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize());
        pairValueDefines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
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

        stringstream reductionSource, extraArgs, deriv0;
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
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string variableName = "dValuedParam_0_"+cu.intToString(i);
            extraArgs << ", const long long* __restrict__ dValue0dParam" << i;
            deriv0 << "real " << variableName << " = (1.0f/0x100000000)*dValue0dParam" << i << "[index];\n";
            for (int j = 0; j < dValuedParam[i]->getBuffers().size(); j++)
                extraArgs << ", real* __restrict__ global_dValuedParam_" << j << "_" << i;
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
        if (needEnergyParamDerivs) {
            map<string, Lepton::ParsedExpression> derivExpressions;
            for (int i = 1; i < force.getNumComputedValues(); i++) {
                for (int j = 0; j < valueParamDerivExpressions[i].size(); j++)
                    derivExpressions["real dValuedParam_"+cu.intToString(i)+"_"+cu.intToString(j)+" = "] = valueParamDerivExpressions[i][j];
                for (int j = 0; j < i; j++)
                    derivExpressions["real dVdV_"+cu.intToString(i)+"_"+cu.intToString(j)+" = "] = valueDerivExpressions[i][j];
            }
            reductionSource << cu.getExpressionUtilities().createExpressions(derivExpressions, variables, functionList, functionDefinitions, "derivChain_temp");
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
            for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
                n2EnergyExpressions["energyParamDeriv"+cu.intToString(j)+" += interactionScale*"] = energyParamDerivExpressions[i][j];
            if (exclude)
                n2EnergySource << "if (!isExcluded) {\n";
            n2EnergySource << cu.getExpressionUtilities().createExpressions(n2EnergyExpressions, variables, functionList, functionDefinitions, "temp");
            if (exclude)
                n2EnergySource << "}\n";
        }
        map<string, string> replacements;
        string n2EnergyStr = n2EnergySource.str();
        replacements["COMPUTE_INTERACTION"] = n2EnergyStr;
        stringstream extraArgs, atomParams, loadLocal1, loadLocal2, clearLocal, load1, load2, declare1, recordDeriv, storeDerivs1, storeDerivs2, initParamDerivs, saveParamDerivs;
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
        if (needEnergyParamDerivs) {
            extraArgs << ", mixed* __restrict__ energyParamDerivs";
            const vector<string>& allParamDerivNames = cu.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[(blockIdx.x*blockDim.x+threadIdx.x)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
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
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
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
        pairEnergyDefines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
        pairEnergyDefines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        pairEnergyDefines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        pairEnergyDefines["NUM_BLOCKS"] = cu.intToString(cu.getNumAtomBlocks());
        pairEnergyDefines["TILE_SIZE"] = cu.intToString(CudaContext::TileSize);
        pairEnergySrc = cu.replaceStrings(CudaKernelSources::customGBEnergyN2, replacements);
    }
    {
        // Create the kernel to reduce the derivatives and calculate per-particle energy terms.

        stringstream compute, extraArgs, load, initParamDerivs, saveParamDerivs;
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
        if (needEnergyParamDerivs) {
            extraArgs << ", mixed* __restrict__ energyParamDerivs";
            const vector<string>& allParamDerivNames = cu.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[(blockIdx.x*blockDim.x+threadIdx.x)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
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
            for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
                expressions["/*"+cu.intToString(i+1)+"*/ energyParamDeriv"+cu.intToString(j)+" += "] = energyParamDerivExpressions[i][j];
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
        replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
        replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
        map<string, string> defines;
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::customGBEnergyPerParticle, replacements), defines);
        perParticleEnergyKernel = cu.getKernel(module, "computePerParticleEnergy");
    }
    if (needParameterGradient || needEnergyParamDerivs) {
        // Create the kernel to compute chain rule terms for computed values that depend explicitly on particle coordinates, and for
        // derivatives with respect to global parameters.

        stringstream compute, extraArgs, initParamDerivs, saveParamDerivs;
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
        if (needEnergyParamDerivs) {
            extraArgs << ", mixed* __restrict__ energyParamDerivs";
            const vector<string>& allParamDerivNames = cu.getEnergyParamDerivNames();
            int numDerivs = allParamDerivNames.size();
            for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
                for (int j = 0; j < dValuedParam[i]->getBuffers().size(); j++)
                    extraArgs << ", real* __restrict__ dValuedParam_" << j << "_" << i;
                initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
                for (int index = 0; index < numDerivs; index++)
                    if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                        saveParamDerivs << "energyParamDerivs[(blockIdx.x*blockDim.x+threadIdx.x)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
            }
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
        if (needParameterGradient) {
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
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        CUmodule module = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::customGBGradientChainRule, replacements), defines);
        gradientChainRuleKernel = cu.getKernel(module, "computeGradientChainRuleTerms");
    }
    {
        // Create the code to calculate chain rule terms as part of the default nonbonded kernel.

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
        if (globals.isInitialized()) {
            globals.upload(globalParamValues);
            arguments.push_back(CudaNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(float), globals.getDevicePointer()));
        }
        cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, force.getNumExclusions() > 0, cutoff, exclusionList, source, force.getForceGroup());
        for (auto param : parameters)
            cu.getNonbondedUtilities().addParameter(param);
        for (auto arg : arguments)
            cu.getNonbondedUtilities().addArgument(arg);
    }
    info = new ForceInfo(force);
    cu.addForce(info);
    cu.addAutoclearBuffer(longEnergyDerivs);
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
            pairValueDefines["CUTOFF"] = cu.doubleToString(cutoff);
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
            pairEnergyDefines["CUTOFF"] = cu.doubleToString(cutoff);
            CUmodule module = cu.createModule(CudaKernelSources::vectorOps+pairEnergySrc, pairEnergyDefines);
            pairEnergyKernel = cu.getKernel(module, "computeN2Energy");
            pairEnergySrc = "";
            pairEnergyDefines.clear();
        }

        // Set arguments for kernels.
        
        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
        valueBuffers.initialize<long long>(cu, cu.getPaddedNumAtoms(), "customGBValueBuffers");
        cu.addAutoclearBuffer(valueBuffers);
        cu.clearBuffer(valueBuffers.getDevicePointer(), sizeof(long long)*valueBuffers.getSize());
        pairValueArgs.push_back(&cu.getPosq().getDevicePointer());
        pairValueArgs.push_back(&cu.getNonbondedUtilities().getExclusions().getDevicePointer());
        pairValueArgs.push_back(&cu.getNonbondedUtilities().getExclusionTiles().getDevicePointer());
        pairValueArgs.push_back(&valueBuffers.getDevicePointer());
        if (nb.getUseCutoff()) {
            pairValueArgs.push_back(&nb.getInteractingTiles().getDevicePointer());
            pairValueArgs.push_back(&nb.getInteractionCount().getDevicePointer());
            pairValueArgs.push_back(cu.getPeriodicBoxSizePointer());
            pairValueArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            pairValueArgs.push_back(cu.getPeriodicBoxVecXPointer());
            pairValueArgs.push_back(cu.getPeriodicBoxVecYPointer());
            pairValueArgs.push_back(cu.getPeriodicBoxVecZPointer());
            pairValueArgs.push_back(&maxTiles);
            pairValueArgs.push_back(&nb.getBlockCenters().getDevicePointer());
            pairValueArgs.push_back(&nb.getBlockBoundingBoxes().getDevicePointer());
            pairValueArgs.push_back(&nb.getInteractingAtoms().getDevicePointer());
        }
        else
            pairValueArgs.push_back(&maxTiles);
        if (globals.isInitialized())
            pairValueArgs.push_back(&globals.getDevicePointer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            if (pairValueUsesParam[i])
                pairValueArgs.push_back(&params->getBuffers()[i].getMemory());
        }
        for (auto& d : dValue0dParam)
            pairValueArgs.push_back(&d.getDevicePointer());
        for (auto& function : tabulatedFunctions)
            pairValueArgs.push_back(&function.getDevicePointer());
        perParticleValueArgs.push_back(&cu.getPosq().getDevicePointer());
        perParticleValueArgs.push_back(&valueBuffers.getDevicePointer());
        if (globals.isInitialized())
            perParticleValueArgs.push_back(&globals.getDevicePointer());
        for (auto& buffer : params->getBuffers())
            perParticleValueArgs.push_back(&buffer.getMemory());
        for (auto& buffer : computedValues->getBuffers())
            perParticleValueArgs.push_back(&buffer.getMemory());
        for (int i = 0; i < dValuedParam.size(); i++) {
            perParticleValueArgs.push_back(&dValue0dParam[i].getDevicePointer());
            for (int j = 0; j < dValuedParam[i]->getBuffers().size(); j++)
                perParticleValueArgs.push_back(&dValuedParam[i]->getBuffers()[j].getMemory());
        }
        for (auto& function : tabulatedFunctions)
            perParticleValueArgs.push_back(&function.getDevicePointer());
        pairEnergyArgs.push_back(&cu.getForce().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getPosq().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getNonbondedUtilities().getExclusions().getDevicePointer());
        pairEnergyArgs.push_back(&cu.getNonbondedUtilities().getExclusionTiles().getDevicePointer());
        pairEnergyArgs.push_back(NULL);
        if (nb.getUseCutoff()) {
            pairEnergyArgs.push_back(&nb.getInteractingTiles().getDevicePointer());
            pairEnergyArgs.push_back(&nb.getInteractionCount().getDevicePointer());
            pairEnergyArgs.push_back(cu.getPeriodicBoxSizePointer());
            pairEnergyArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            pairEnergyArgs.push_back(cu.getPeriodicBoxVecXPointer());
            pairEnergyArgs.push_back(cu.getPeriodicBoxVecYPointer());
            pairEnergyArgs.push_back(cu.getPeriodicBoxVecZPointer());
            pairEnergyArgs.push_back(&maxTiles);
            pairEnergyArgs.push_back(&nb.getBlockCenters().getDevicePointer());
            pairEnergyArgs.push_back(&nb.getBlockBoundingBoxes().getDevicePointer());
            pairEnergyArgs.push_back(&nb.getInteractingAtoms().getDevicePointer());
        }
        else
            pairEnergyArgs.push_back(&maxTiles);
        if (globals.isInitialized())
            pairEnergyArgs.push_back(&globals.getDevicePointer());
        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
            if (pairEnergyUsesParam[i])
                pairEnergyArgs.push_back(&params->getBuffers()[i].getMemory());
        }
        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
            if (pairEnergyUsesValue[i])
                pairEnergyArgs.push_back(&computedValues->getBuffers()[i].getMemory());
        }
        pairEnergyArgs.push_back(&longEnergyDerivs.getDevicePointer());
        if (needEnergyParamDerivs)
            pairEnergyArgs.push_back(&cu.getEnergyParamDerivBuffer().getDevicePointer());
        for (auto& function : tabulatedFunctions)
            pairEnergyArgs.push_back(&function.getDevicePointer());
        perParticleEnergyArgs.push_back(&cu.getForce().getDevicePointer());
        perParticleEnergyArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        perParticleEnergyArgs.push_back(&cu.getPosq().getDevicePointer());
        if (globals.isInitialized())
            perParticleEnergyArgs.push_back(&globals.getDevicePointer());
        for (auto& buffer : params->getBuffers())
            perParticleEnergyArgs.push_back(&buffer.getMemory());
        for (auto& buffer : computedValues->getBuffers())
            perParticleEnergyArgs.push_back(&buffer.getMemory());
        for (auto& buffer : energyDerivs->getBuffers())
            perParticleEnergyArgs.push_back(&buffer.getMemory());
        for (auto& buffer : energyDerivChain->getBuffers())
            perParticleEnergyArgs.push_back(&buffer.getMemory());
        perParticleEnergyArgs.push_back(&longEnergyDerivs.getDevicePointer());
        if (needEnergyParamDerivs)
            perParticleEnergyArgs.push_back(&cu.getEnergyParamDerivBuffer().getDevicePointer());
        for (auto& function : tabulatedFunctions)
            perParticleEnergyArgs.push_back(&function.getDevicePointer());
        if (needParameterGradient || needEnergyParamDerivs) {
            gradientChainRuleArgs.push_back(&cu.getForce().getDevicePointer());
            gradientChainRuleArgs.push_back(&cu.getPosq().getDevicePointer());
            if (globals.isInitialized())
                gradientChainRuleArgs.push_back(&globals.getDevicePointer());
            for (auto& buffer : params->getBuffers())
                gradientChainRuleArgs.push_back(&buffer.getMemory());
            for (auto& buffer : computedValues->getBuffers())
                gradientChainRuleArgs.push_back(&buffer.getMemory());
            for (auto& buffer : energyDerivs->getBuffers())
                gradientChainRuleArgs.push_back(&buffer.getMemory());
            if (needEnergyParamDerivs) {
                gradientChainRuleArgs.push_back(&cu.getEnergyParamDerivBuffer().getDevicePointer());
                for (auto d : dValuedParam)
                    for (auto& buffer : d->getBuffers())
                        gradientChainRuleArgs.push_back(&buffer.getMemory());
            }
            for (auto& function : tabulatedFunctions)
                gradientChainRuleArgs.push_back(&function.getDevicePointer());
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
    pairEnergyArgs[5] = &includeEnergy;
    if (nb.getUseCutoff()) {
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            pairValueArgs[4] = &nb.getInteractingTiles().getDevicePointer();
            pairEnergyArgs[6] = &nb.getInteractingTiles().getDevicePointer();
            pairValueArgs[14] = &nb.getInteractingAtoms().getDevicePointer();
            pairEnergyArgs[16] = &nb.getInteractingAtoms().getDevicePointer();
        }
    }
    cu.executeKernel(pairValueKernel, &pairValueArgs[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cu.executeKernel(perParticleValueKernel, &perParticleValueArgs[0], cu.getPaddedNumAtoms());
    cu.executeKernel(pairEnergyKernel, &pairEnergyArgs[0], nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    cu.executeKernel(perParticleEnergyKernel, &perParticleEnergyArgs[0], cu.getPaddedNumAtoms());
    if (needParameterGradient || needEnergyParamDerivs)
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

class CudaCalcCustomExternalForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomExternalForceKernel::~CudaCalcCustomExternalForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
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
    info = new ForceInfo(force, system.getNumParticles());
    cu.addForce(info);

    // Record information for the expressions.

    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    map<string, Lepton::CustomFunction*> customFunctions;
    customFunctions["periodicdistance"] = cu.getExpressionUtilities().getPeriodicDistancePlaceholder();
    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), customFunctions).optimize();
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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customExternalGlobals");
        globals.upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals.getDevicePointer(), "float");
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

class CudaCalcCustomHbondForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomHbondForce& force) : force(force) {
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
    donors.initialize<int4>(cu, numDonors, "customHbondDonors");
    acceptors.initialize<int4>(cu, numAcceptors, "customHbondAcceptors");
    donorParams = new CudaParameterSet(cu, force.getNumPerDonorParameters(), numDonors, "customHbondDonorParameters");
    acceptorParams = new CudaParameterSet(cu, force.getNumPerAcceptorParameters(), numAcceptors, "customHbondAcceptorParameters");
    if (force.getNumGlobalParameters() > 0)
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customHbondGlobals");
    vector<vector<float> > donorParamVector(numDonors);
    vector<int4> donorVector(numDonors);
    for (int i = 0; i < numDonors; i++) {
        vector<double> parameters;
        force.getDonorParameters(startIndex+i, donorVector[i].x, donorVector[i].y, donorVector[i].z, parameters);
        donorParamVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++)
            donorParamVector[i][j] = (float) parameters[j];
    }
    donors.upload(donorVector);
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
    acceptors.upload(acceptorVector);
    acceptorParams->setParameterValues(acceptorParamVector);
    info = new ForceInfo(force);
    cu.addForce(info);

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
    donorExclusions.initialize<int4>(cu, numDonors, "customHbondDonorExclusions");
    acceptorExclusions.initialize<int4>(cu, numAcceptors, "customHbondAcceptorExclusions");
    donorExclusions.upload(donorExclusionVector);
    acceptorExclusions.upload(acceptorExclusionVector);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
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
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real4 delta"+deltaName+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n");
            computedDeltas.insert(deltaName);
        }
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real r_"+deltaName+" = SQRT(delta"+deltaName+".w);\n");
        variables[distance.first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
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
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
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
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Next it needs to load parameters from global memory.

    if (force.getNumGlobalParameters() > 0)
        extraArgs << ", const float* __restrict__ globals";
    for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
        extraArgs << ", const "+buffer.getType()+"* __restrict__ donor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" donorParams"+cu.intToString(i+1)+" = donor"+buffer.getName()+"[donorIndex];\n");
    }
    for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
        extraArgs << ", const "+buffer.getType()+"* __restrict__ acceptor"+buffer.getName();
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" acceptorParams"+cu.intToString(i+1)+" = acceptor"+buffer.getName()+"[acceptorIndex];\n");
    }

    // Now evaluate the expressions.

    computeAcceptor << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");
    forceExpressions["energy += "] = energyExpression;
    computeDonor << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to atoms.

    index = 0;
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*delta"+deltaName;
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
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 crossProd = cross(delta"+deltaName2+", delta"+deltaName1+");\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real lengthCross = max(SQRT(dot(crossProd,crossProd)), 1e-6f);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 deltaCross0 = -cross(trim(delta"+deltaName1+"), crossProd)*dEdAngle"+cu.intToString(index)+"/(delta"+deltaName1+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 deltaCross2 = cross(trim(delta"+deltaName2+"), crossProd)*dEdAngle"+cu.intToString(index)+"/(delta"+deltaName2+".w*lengthCross);\n");
        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "real3 deltaCross1 = -(deltaCross0+deltaCross2);\n");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "deltaCross0");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "deltaCross1");
        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "deltaCross2");
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
        index++;
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
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        int index = 0;
        donorArgs.push_back(&cu.getForce().getDevicePointer());
        donorArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        donorArgs.push_back(&cu.getPosq().getDevicePointer());
        donorArgs.push_back(&donorExclusions.getDevicePointer());
        donorArgs.push_back(&donors.getDevicePointer());
        donorArgs.push_back(&acceptors.getDevicePointer());
        donorArgs.push_back(cu.getPeriodicBoxSizePointer());
        donorArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        donorArgs.push_back(cu.getPeriodicBoxVecXPointer());
        donorArgs.push_back(cu.getPeriodicBoxVecYPointer());
        donorArgs.push_back(cu.getPeriodicBoxVecZPointer());
        if (globals.isInitialized())
            donorArgs.push_back(&globals.getDevicePointer());
        for (auto& buffer : donorParams->getBuffers())
            donorArgs.push_back(&buffer.getMemory());
        for (auto& buffer : acceptorParams->getBuffers())
            donorArgs.push_back(&buffer.getMemory());
        for (auto& function : tabulatedFunctions)
            donorArgs.push_back(&function.getDevicePointer());
        index = 0;
        acceptorArgs.push_back(&cu.getForce().getDevicePointer());
        acceptorArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        acceptorArgs.push_back(&cu.getPosq().getDevicePointer());
        acceptorArgs.push_back(&acceptorExclusions.getDevicePointer());
        acceptorArgs.push_back(&donors.getDevicePointer());
        acceptorArgs.push_back(&acceptors.getDevicePointer());
        acceptorArgs.push_back(cu.getPeriodicBoxSizePointer());
        acceptorArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        acceptorArgs.push_back(cu.getPeriodicBoxVecXPointer());
        acceptorArgs.push_back(cu.getPeriodicBoxVecYPointer());
        acceptorArgs.push_back(cu.getPeriodicBoxVecZPointer());
        if (globals.isInitialized())
            acceptorArgs.push_back(&globals.getDevicePointer());
        for (auto& buffer : donorParams->getBuffers())
            acceptorArgs.push_back(&buffer.getMemory());
        for (auto& buffer : acceptorParams->getBuffers())
            acceptorArgs.push_back(&buffer.getMemory());
        for (auto& function : tabulatedFunctions)
            acceptorArgs.push_back(&function.getDevicePointer());
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

class CudaCalcCustomCentroidBondForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomCentroidBondForceKernel::~CudaCalcCustomCentroidBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
}

void CudaCalcCustomCentroidBondForceKernel::initialize(const System& system, const CustomCentroidBondForce& force) {
    cu.setAsCurrent();
    numBonds = force.getNumBonds();
    if (numBonds == 0)
        return;
    info = new ForceInfo(force);
    cu.addForce(info);
    
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
    groupParticles.initialize<int>(cu, groupParticleVec.size(), "groupParticles");
    groupParticles.upload(groupParticleVec);
    if (cu.getUseDoublePrecision()) {
        groupWeights.initialize<double>(cu, groupParticleVec.size(), "groupWeights");
        centerPositions.initialize<double4>(cu, numGroups, "centerPositions");
    }
    else {
        groupWeights.initialize<float>(cu, groupParticleVec.size(), "groupWeights");
        centerPositions.initialize<float4>(cu, numGroups, "centerPositions");
    }
    groupWeights.upload(groupWeightVec, true);
    groupOffsets.initialize<int>(cu, groupOffsetVec.size(), "groupOffsets");
    groupOffsets.upload(groupOffsetVec);
    groupForces.initialize<long long>(cu, numGroups*3, "groupForces");
    cu.addAutoclearBuffer(groupForces);
    
    // Record the bonds.
    
    int groupsPerBond = force.getNumGroupsPerBond();
    vector<int> bondGroupVec(numBonds*groupsPerBond);
    params = new CudaParameterSet(cu, force.getNumPerBondParameters(), numBonds, "customCentroidBondParams");
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
    bondGroups.initialize<int>(cu, bondGroupVec.size(), "bondGroups");
    bondGroups.upload(bondGroupVec);
    
    // Record the arguments to the force kernel.
    
    groupForcesArgs.push_back(&groupForces.getDevicePointer());
    groupForcesArgs.push_back(NULL); // Energy buffer hasn't been created yet
    groupForcesArgs.push_back(&centerPositions.getDevicePointer());
    groupForcesArgs.push_back(&bondGroups.getDevicePointer());
    groupForcesArgs.push_back(cu.getPeriodicBoxSizePointer());
    groupForcesArgs.push_back(cu.getInvPeriodicBoxSizePointer());
    groupForcesArgs.push_back(cu.getPeriodicBoxVecXPointer());
    groupForcesArgs.push_back(cu.getPeriodicBoxVecYPointer());
    groupForcesArgs.push_back(cu.getPeriodicBoxVecZPointer());
    needEnergyParamDerivs = (force.getNumEnergyParameterDerivatives() > 0);
    if (needEnergyParamDerivs)
        groupForcesArgs.push_back(NULL); // Derivatives buffer hasn't been created yet

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream extraArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        extraArgs << ", const float";
        if (width > 1)
            extraArgs << width;
        extraArgs << "* __restrict__ " << arrayName;
        groupForcesArgs.push_back(&tabulatedFunctions[i].getDevicePointer());
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
        string index = cu.intToString(i+1);
        variables["x"+index] = "pos"+index+".x";
        variables["y"+index] = "pos"+index+".y";
        variables["z"+index] = "pos"+index+".z";
    }
    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
        const string& name = force.getPerBondParameterName(i);
        variables[name] = "bondParams"+params->getParameterSuffix(i);
    }
    if (needEnergyParamDerivs)
        extraArgs << ", mixed* __restrict__ energyParamDerivs";
    if (force.getNumGlobalParameters() > 0) {
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customCentroidBondGlobals");
        globals.upload(globalParamValues);
        extraArgs << ", const float* __restrict__ globals";
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            const string& name = force.getGlobalParameterName(i);
            string value = "globals["+cu.intToString(i)+"]";
            variables[name] = value;
        }
        groupForcesArgs.push_back(&globals.getDevicePointer());
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
        string index = cu.intToString(i+1);
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
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
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
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
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
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArgs<<", const "<<buffer.getType()<<"* __restrict__ globalParams"<<i;
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = globalParams"<<i<<"[index];\n";
        groupForcesArgs.push_back(&buffer.getMemory());
    }
    forceExpressions["energy += "] = energyExpression;
    if (needEnergyParamDerivs) {
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
            string paramName = force.getEnergyParameterDerivativeName(i);
            cu.addEnergyParameterDerivative(paramName);
            Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
            forceExpressions[string("energyParamDeriv")+cu.intToString(i)+" += "] = derivExpression;
            initParamDerivs << "mixed energyParamDeriv" << i << " = 0;\n";
        }
        const vector<string>& allParamDerivNames = cu.getEnergyParamDerivNames();
        int numDerivs = allParamDerivNames.size();
        for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
            for (int index = 0; index < numDerivs; index++)
                if (allParamDerivNames[index] == force.getEnergyParameterDerivativeName(i))
                    saveParamDerivs << "energyParamDerivs[(blockIdx.x*blockDim.x+threadIdx.x)*" << numDerivs << "+" << index << "] += energyParamDeriv" << i << ";\n";
    }
    compute << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp");

    // Finally, apply forces to groups.

    vector<string> forceNames;
    for (int i = 0; i < groupsPerBond; i++) {
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
    for (auto& distance : distances) {
        const vector<int>& groups = distance.second;
        string deltaName = atomNames[groups[0]]+atomNames[groups[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*trim(delta"+deltaName+")";
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
        compute<<"real3 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(trim(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(trim(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
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
        compute<<"ff.x = (-dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
        compute<<"ff.w = (dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
        compute<<"real3 internalF0 = ff.x*trim("<<crossName1<<");\n";
        compute<<"real3 internalF3 = ff.w*trim("<<crossName2<<");\n";
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
        compute<<"atomicAdd(&groupForce[group"<<(i+1)<<"], static_cast<unsigned long long>((long long) (force"<<(i+1)<<".x*0x100000000)));\n";
        compute<<"atomicAdd(&groupForce[group"<<(i+1)<<"+NUM_GROUPS], static_cast<unsigned long long>((long long) (force"<<(i+1)<<".y*0x100000000)));\n";
        compute<<"atomicAdd(&groupForce[group"<<(i+1)<<"+NUM_GROUPS*2], static_cast<unsigned long long>((long long) (force"<<(i+1)<<".z*0x100000000)));\n";
        compute<<"__threadfence_block();\n";
    }
    map<string, string> replacements;
    replacements["M_PI"] = cu.doubleToString(M_PI);
    replacements["NUM_GROUPS"] = cu.intToString(numGroups);
    replacements["NUM_BONDS"] = cu.intToString(numBonds);
    replacements["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    replacements["EXTRA_ARGS"] = extraArgs.str();
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["INIT_PARAM_DERIVS"] = initParamDerivs.str();
    replacements["SAVE_PARAM_DERIVS"] = saveParamDerivs.str();
    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+cu.replaceStrings(CudaKernelSources::customCentroidBond, replacements));
    computeCentersKernel = cu.getKernel(module, "computeGroupCenters");
    groupForcesKernel = cu.getKernel(module, "computeGroupForces");
    applyForcesKernel = cu.getKernel(module, "applyForcesToAtoms");
}

double CudaCalcCustomCentroidBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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
    void* computeCentersArgs[] = {&cu.getPosq().getDevicePointer(), &groupParticles.getDevicePointer(), &groupWeights.getDevicePointer(),
            &groupOffsets.getDevicePointer(), &centerPositions.getDevicePointer()};
    cu.executeKernel(computeCentersKernel, computeCentersArgs, CudaContext::TileSize*numGroups);
    groupForcesArgs[1] = &cu.getEnergyBuffer().getDevicePointer();
    if (needEnergyParamDerivs)
        groupForcesArgs[9] = &cu.getEnergyParamDerivBuffer().getDevicePointer();
    cu.executeKernel(groupForcesKernel, &groupForcesArgs[0], numBonds);
    void* applyForcesArgs[] = {&groupParticles.getDevicePointer(), &groupWeights.getDevicePointer(), &groupOffsets.getDevicePointer(),
            &groupForces.getDevicePointer(), &cu.getForce().getDevicePointer()};
    cu.executeKernel(applyForcesKernel, applyForcesArgs, CudaContext::TileSize*numGroups);
    return 0.0;
}

void CudaCalcCustomCentroidBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force) {
    cu.setAsCurrent();
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
    
    cu.invalidateMolecules();
}

class CudaCalcCustomCompoundBondForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomCompoundBondForceKernel::~CudaCalcCustomCompoundBondForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
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
    info = new ForceInfo(force);
    cu.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
        string arrayName = cu.getBondedUtilities().addArgument(tabulatedFunctions[i].getDevicePointer(), width == 1 ? "float" : "float"+cu.intToString(width));
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
        globals.initialize<float>(cu, force.getNumGlobalParameters(), "customCompoundBondGlobals");
        globals.upload(globalParamValues);
        string argName = cu.getBondedUtilities().addArgument(globals.getDevicePointer(), "float");
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
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", "<<force.usesPeriodicBoundaryConditions()<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
        variables[distance.first] = "r_"+deltaName;
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
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
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
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
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
    }

    // Now evaluate the expressions.

    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
    }
    forceExpressions["energy += "] = energyExpression;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string paramName = force.getEnergyParameterDerivativeName(i);
        string derivVariable = cu.getBondedUtilities().addEnergyParameterDerivative(paramName);
        Lepton::ParsedExpression derivExpression = energyExpression.differentiate(paramName).optimize();
        forceExpressions[derivVariable+" += "] = derivExpression;
    }
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
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*ccb_trim(delta"+deltaName+")";
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
        compute<<"real3 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(ccb_trim(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(ccb_trim(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
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
        index++;
    }
    cu.getBondedUtilities().addInteraction(atoms, compute.str(), force.getForceGroup());
    map<string, string> replacements;
    replacements["M_PI"] = cu.doubleToString(M_PI);
    cu.getBondedUtilities().addPrefixCode(cu.replaceStrings(CudaKernelSources::customCompoundBond, replacements));
}

double CudaCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
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

class CudaCalcCustomManyParticleForceKernel::ForceInfo : public CudaForceInfo {
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

CudaCalcCustomManyParticleForceKernel::~CudaCalcCustomManyParticleForceKernel() {
    cu.setAsCurrent();
    if (params != NULL)
        delete params;
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
    info = new ForceInfo(force);
    cu.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = "table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
        tabulatedFunctions[i].upload(f);
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
        particleTypes.initialize<int>(cu, particleTypesVec.size(), "customManyParticleTypes");
        orderIndex.initialize<int>(cu, orderIndexVec.size(), "customManyParticleOrderIndex");
        particleOrder.initialize<int>(cu, particleOrderVec.size()*particlesPerSet, "customManyParticleOrder");
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
        exclusions.initialize<int>(cu, exclusionsVec.size(), "customManyParticleExclusions");
        exclusionStartIndex.initialize<int>(cu, exclusionStartIndexVec.size(), "customManyParticleExclusionStart");
        exclusions.upload(exclusionsVec);
        exclusionStartIndex.upload(exclusionStartIndexVec);
    }
    
    // Build data structures for the neighbor list.
    
    if (nonbondedMethod != NoCutoff) {
        int numAtomBlocks = cu.getNumAtomBlocks();
        int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
        blockCenter.initialize(cu, numAtomBlocks, 4*elementSize, "blockCenter");
        blockBoundingBox.initialize(cu, numAtomBlocks, 4*elementSize, "blockBoundingBox");
        numNeighborPairs.initialize<int>(cu, 1, "customManyParticleNumNeighborPairs");
        neighborStartIndex.initialize<int>(cu, numParticles+1, "customManyParticleNeighborStartIndex");
        numNeighborsForAtom.initialize<int>(cu, numParticles, "customManyParticleNumNeighborsForAtom");
        CHECK_RESULT(cuEventCreate(&event, CU_EVENT_DISABLE_TIMING), "Error creating event for CustomManyParticleForce");

        // Select a size for the array that holds the neighbor list.  We have to make a fairly
        // arbitrary guess, but if this turns out to be too small we'll increase it later.

        maxNeighborPairs = 150*numParticles;
        neighborPairs.initialize<int2>(cu, maxNeighborPairs, "customManyParticleNeighborPairs");
        neighbors.initialize<int>(cu, maxNeighborPairs, "customManyParticleNeighbors");
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
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        if (computedDeltas.count(deltaName) == 0) {
            compute<<"real4 delta"<<deltaName<<" = delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<", periodicBoxSize, invPeriodicBoxSize, periodicBoxVecX, periodicBoxVecY, periodicBoxVecZ);\n";
            computedDeltas.insert(deltaName);
        }
        compute<<"real r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
        variables.push_back(makeVariable(distance.first, "r_"+deltaName));
        forceExpressions["real dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(distance.first).optimize();
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
        forceExpressions["real dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(angle.first).optimize();
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
        forceExpressions["real dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(dihedral.first).optimize();
        index++;
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
    for (auto& distance : distances) {
        const vector<int>& atoms = distance.second;
        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*trim(delta"+deltaName+")";
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
        compute<<"real3 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
        compute<<"real lengthCross = max(SQRT(dot(crossProd, crossProd)), 1e-6f);\n";
        compute<<"real3 deltaCross0 = -cross(trim(delta"<<deltaName1<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
        compute<<"real3 deltaCross2 = cross(trim(delta"<<deltaName2<<"), crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
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
        forceArgs.push_back(cu.getPeriodicBoxVecXPointer());
        forceArgs.push_back(cu.getPeriodicBoxVecYPointer());
        forceArgs.push_back(cu.getPeriodicBoxVecZPointer());
        if (nonbondedMethod != NoCutoff) {
            forceArgs.push_back(&neighbors.getDevicePointer());
            forceArgs.push_back(&neighborStartIndex.getDevicePointer());
        }
        if (particleTypes.isInitialized()) {
            forceArgs.push_back(&particleTypes.getDevicePointer());
            forceArgs.push_back(&orderIndex.getDevicePointer());
            forceArgs.push_back(&particleOrder.getDevicePointer());
        }
        if (exclusions.isInitialized()) {
            forceArgs.push_back(&exclusions.getDevicePointer());
            forceArgs.push_back(&exclusionStartIndex.getDevicePointer());
        }
        for (auto& buffer : params->getBuffers())
            forceArgs.push_back(&buffer.getMemory());
        for (auto& function : tabulatedFunctions)
            forceArgs.push_back(&function.getDevicePointer());
        
        if (nonbondedMethod != NoCutoff) {
            // Set arguments for the block bounds kernel.

            blockBoundsArgs.push_back(cu.getPeriodicBoxSizePointer());
            blockBoundsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            blockBoundsArgs.push_back(cu.getPeriodicBoxVecXPointer());
            blockBoundsArgs.push_back(cu.getPeriodicBoxVecYPointer());
            blockBoundsArgs.push_back(cu.getPeriodicBoxVecZPointer());
            blockBoundsArgs.push_back(&cu.getPosq().getDevicePointer());
            blockBoundsArgs.push_back(&blockCenter.getDevicePointer());
            blockBoundsArgs.push_back(&blockBoundingBox.getDevicePointer());
            blockBoundsArgs.push_back(&numNeighborPairs.getDevicePointer());

            // Set arguments for the neighbor list kernel.

            neighborsArgs.push_back(cu.getPeriodicBoxSizePointer());
            neighborsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            neighborsArgs.push_back(cu.getPeriodicBoxVecXPointer());
            neighborsArgs.push_back(cu.getPeriodicBoxVecYPointer());
            neighborsArgs.push_back(cu.getPeriodicBoxVecZPointer());
            neighborsArgs.push_back(&cu.getPosq().getDevicePointer());
            neighborsArgs.push_back(&blockCenter.getDevicePointer());
            neighborsArgs.push_back(&blockBoundingBox.getDevicePointer());
            neighborsArgs.push_back(&neighborPairs.getDevicePointer());
            neighborsArgs.push_back(&numNeighborPairs.getDevicePointer());
            neighborsArgs.push_back(&numNeighborsForAtom.getDevicePointer());
            neighborsArgs.push_back(&maxNeighborPairs);
            if (exclusions.isInitialized()) {
                neighborsArgs.push_back(&exclusions.getDevicePointer());
                neighborsArgs.push_back(&exclusionStartIndex.getDevicePointer());
            }
            
            // Set arguments for the kernel to find neighbor list start indices.
            
            startIndicesArgs.push_back(&numNeighborsForAtom.getDevicePointer());
            startIndicesArgs.push_back(&neighborStartIndex.getDevicePointer());
            startIndicesArgs.push_back(&numNeighborPairs.getDevicePointer());
            startIndicesArgs.push_back(&maxNeighborPairs);

            // Set arguments for the kernel to assemble the final neighbor list.
            
            copyPairsArgs.push_back(&neighborPairs.getDevicePointer());
            copyPairsArgs.push_back(&neighbors.getDevicePointer());
            copyPairsArgs.push_back(&numNeighborPairs.getDevicePointer());
            copyPairsArgs.push_back(&maxNeighborPairs);
            copyPairsArgs.push_back(&numNeighborsForAtom.getDevicePointer());
            copyPairsArgs.push_back(&neighborStartIndex.getDevicePointer());
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

            numNeighborPairs.download(numPairs, false);
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

                maxNeighborPairs = (int) (1.1*(*numPairs));
                neighborPairs.resize(maxNeighborPairs);
                neighbors.resize(maxNeighborPairs);
                forceArgs[5] = &neighbors.getDevicePointer();
                neighborsArgs[5] = &neighborPairs.getDevicePointer();
                copyPairsArgs[0] = &neighborPairs.getDevicePointer();
                copyPairsArgs[1] = &neighbors.getDevicePointer();
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

class CudaCalcGayBerneForceKernel::ForceInfo : public CudaForceInfo {
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

class CudaCalcGayBerneForceKernel::ReorderListener : public CudaContext::ReorderListener {
public:
    ReorderListener(CudaCalcGayBerneForceKernel& owner) : owner(owner) {
    }
    void execute() {
        owner.sortAtoms();
    }
private:
    CudaCalcGayBerneForceKernel& owner;
};

void CudaCalcGayBerneForceKernel::initialize(const System& system, const GayBerneForce& force) {
    // Initialize interactions.

    int numParticles = force.getNumParticles();
    sigParams.initialize<float4>(cu, cu.getPaddedNumAtoms(), "sigParams");
    epsParams.initialize<float2>(cu, cu.getPaddedNumAtoms(), "epsParams");
    scale.initialize<float4>(cu, cu.getPaddedNumAtoms(), "scale");
    axisParticleIndices.initialize<int2>(cu, cu.getPaddedNumAtoms(), "axisParticleIndices");
    sortedParticles.initialize<int>(cu, cu.getPaddedNumAtoms(), "sortedParticles");
    aMatrix.initialize<float>(cu, 9*cu.getPaddedNumAtoms(), "aMatrix");
    bMatrix.initialize<float>(cu, 9*cu.getPaddedNumAtoms(), "bMatrix");
    gMatrix.initialize<float>(cu, 9*cu.getPaddedNumAtoms(), "gMatrix");
    vector<float4> sigParamsVector(cu.getPaddedNumAtoms(), make_float4(0, 0, 0, 0));
    vector<float2> epsParamsVector(cu.getPaddedNumAtoms(), make_float2(0, 0));
    vector<float4> scaleVector(cu.getPaddedNumAtoms(), make_float4(0, 0, 0, 0));
    vector<int2> axisParticleVector(cu.getPaddedNumAtoms(), make_int2(0, 0));
    isRealParticle.resize(cu.getPaddedNumAtoms());
    for (int i = 0; i < numParticles; i++) {
        int xparticle, yparticle;
        double sigma, epsilon, sx, sy, sz, ex, ey, ez;
        force.getParticleParameters(i, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
        axisParticleVector[i] = make_int2(xparticle, yparticle);
        sigParamsVector[i] = make_float4((float) (0.5*sigma), (float) (0.25*sx*sx), (float) (0.25*sy*sy), (float) (0.25*sz*sz));
        epsParamsVector[i] = make_float2((float) sqrt(epsilon), (float) (0.125*(sx*sy + sz*sz)*sqrt(sx*sy)));
        scaleVector[i] = make_float4((float) (1/sqrt(ex)), (float) (1/sqrt(ey)), (float) (1/sqrt(ez)), 0);
        isRealParticle[i] = (epsilon != 0.0);
    }
    sigParams.upload(sigParamsVector);
    epsParams.upload(epsParamsVector);
    scale.upload(scaleVector);
    axisParticleIndices.upload(axisParticleVector);
    
    // Record exceptions and exclusions.

    vector<float2> exceptionParamsVec;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, sigma, epsilon);
        if (epsilon != 0.0) {
            exceptionParamsVec.push_back(make_float2((float) sigma, (float) epsilon));
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
    numExceptions = exceptionParamsVec.size();
    exclusions.initialize<int>(cu, max(1, (int) excludedPairs.size()), "exclusions");
    exclusionStartIndex.initialize<int>(cu, numRealParticles+1, "exclusionStartIndex");
    exceptionParticles.initialize<int4>(cu, max(1, numExceptions), "exceptionParticles");
    exceptionParams.initialize<float2>(cu, max(1, numExceptions), "exceptionParams");
    if (numExceptions > 0)
        exceptionParams.upload(exceptionParamsVec);
    
    // Create data structures used for the neighbor list.

    int numAtomBlocks = (numRealParticles+31)/32;
    int elementSize = (cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    blockCenter.initialize(cu, numAtomBlocks, 4*elementSize, "blockCenter");
    blockBoundingBox.initialize(cu, numAtomBlocks, 4*elementSize, "blockBoundingBox");
    sortedPos.initialize(cu, numRealParticles, 4*elementSize, "sortedPos");
    maxNeighborBlocks = numRealParticles*2;
    neighbors.initialize<int>(cu, maxNeighborBlocks*32, "neighbors");
    neighborIndex.initialize<int>(cu, maxNeighborBlocks, "neighborIndex");
    neighborBlockCount.initialize<int>(cu, 1, "neighborBlockCount");
    if (force.getNonbondedMethod() != GayBerneForce::NoCutoff)
        CHECK_RESULT(cuEventCreate(&event, CU_EVENT_DISABLE_TIMING), "Error creating event for CustomManyParticleForce");

    // Create array for accumulating torques.
    
    torque.initialize<long long>(cu, 3*cu.getPaddedNumAtoms(), "torque");
    cu.addAutoclearBuffer(torque);

    // Create the kernels.
    
    nonbondedMethod = force.getNonbondedMethod();
    bool useCutoff = (nonbondedMethod != GayBerneForce::NoCutoff);
    bool usePeriodic = (nonbondedMethod == GayBerneForce::CutoffPeriodic);
    map<string, string> defines;
    defines["USE_SWITCH"] = (useCutoff && force.getUseSwitchingFunction() ? "1" : "0");
    double cutoff = force.getCutoffDistance();
    defines["CUTOFF_SQUARED"] = cu.doubleToString(cutoff*cutoff);
    if (useCutoff) {
        defines["USE_CUTOFF"] = 1;
        if (usePeriodic)
            defines["USE_PERIODIC"] = "1";
        
        // Compute the switching coefficients.
        
        if (force.getUseSwitchingFunction()) {
            defines["SWITCH_CUTOFF"] = cu.doubleToString(force.getSwitchingDistance());
            defines["SWITCH_C3"] = cu.doubleToString(10/pow(force.getSwitchingDistance()-cutoff, 3.0));
            defines["SWITCH_C4"] = cu.doubleToString(15/pow(force.getSwitchingDistance()-cutoff, 4.0));
            defines["SWITCH_C5"] = cu.doubleToString(6/pow(force.getSwitchingDistance()-cutoff, 5.0));
        }
    }
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaKernelSources::gayBerne, defines);
    framesKernel = cu.getKernel(module, "computeEllipsoidFrames");
    blockBoundsKernel = cu.getKernel(module, "findBlockBounds");
    neighborsKernel = cu.getKernel(module, "findNeighbors");
    forceKernel = cu.getKernel(module, "computeForce");
    torqueKernel = cu.getKernel(module, "applyTorques");
    info = new ForceInfo(force);
    cu.addForce(info);
    cu.addReorderListener(new ReorderListener(*this));
}

double CudaCalcGayBerneForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        sortAtoms();
        framesArgs.push_back(&numRealParticles);
        framesArgs.push_back(&cu.getPosq().getDevicePointer());
        framesArgs.push_back(&axisParticleIndices.getDevicePointer());
        framesArgs.push_back(&sigParams.getDevicePointer());
        framesArgs.push_back(&scale.getDevicePointer());
        framesArgs.push_back(&aMatrix.getDevicePointer());
        framesArgs.push_back(&bMatrix.getDevicePointer());
        framesArgs.push_back(&gMatrix.getDevicePointer());
        framesArgs.push_back(&sortedParticles.getDevicePointer());
        blockBoundsArgs.push_back(&numRealParticles);
        blockBoundsArgs.push_back(cu.getPeriodicBoxSizePointer());
        blockBoundsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        blockBoundsArgs.push_back(cu.getPeriodicBoxVecXPointer());
        blockBoundsArgs.push_back(cu.getPeriodicBoxVecYPointer());
        blockBoundsArgs.push_back(cu.getPeriodicBoxVecZPointer());
        blockBoundsArgs.push_back(&sortedParticles.getDevicePointer());
        blockBoundsArgs.push_back(&cu.getPosq().getDevicePointer());
        blockBoundsArgs.push_back(&sortedPos.getDevicePointer());
        blockBoundsArgs.push_back(&blockCenter.getDevicePointer());
        blockBoundsArgs.push_back(&blockBoundingBox.getDevicePointer());
        blockBoundsArgs.push_back(&neighborBlockCount.getDevicePointer());
        neighborsArgs.push_back(&numRealParticles);
        neighborsArgs.push_back(&maxNeighborBlocks);
        neighborsArgs.push_back(cu.getPeriodicBoxSizePointer());
        neighborsArgs.push_back(cu.getInvPeriodicBoxSizePointer());
        neighborsArgs.push_back(cu.getPeriodicBoxVecXPointer());
        neighborsArgs.push_back(cu.getPeriodicBoxVecYPointer());
        neighborsArgs.push_back(cu.getPeriodicBoxVecZPointer());
        neighborsArgs.push_back(&sortedPos.getDevicePointer());
        neighborsArgs.push_back(&blockCenter.getDevicePointer());
        neighborsArgs.push_back(&blockBoundingBox.getDevicePointer());
        neighborsArgs.push_back(&neighbors.getDevicePointer());
        neighborsArgs.push_back(&neighborIndex.getDevicePointer());
        neighborsArgs.push_back(&neighborBlockCount.getDevicePointer());
        neighborsArgs.push_back(&exclusions.getDevicePointer());
        neighborsArgs.push_back(&exclusionStartIndex.getDevicePointer());
        forceArgs.push_back(&cu.getForce().getDevicePointer());
        forceArgs.push_back(&torque.getDevicePointer());
        forceArgs.push_back(&numRealParticles);
        forceArgs.push_back(&numExceptions);
        forceArgs.push_back(&cu.getEnergyBuffer().getDevicePointer());
        forceArgs.push_back(&sortedPos.getDevicePointer());
        forceArgs.push_back(&sigParams.getDevicePointer());
        forceArgs.push_back(&epsParams.getDevicePointer());
        forceArgs.push_back(&sortedParticles.getDevicePointer());
        forceArgs.push_back(&aMatrix.getDevicePointer());
        forceArgs.push_back(&bMatrix.getDevicePointer());
        forceArgs.push_back(&gMatrix.getDevicePointer());
        forceArgs.push_back(&exclusions.getDevicePointer());
        forceArgs.push_back(&exclusionStartIndex.getDevicePointer());
        forceArgs.push_back(&exceptionParticles.getDevicePointer());
        forceArgs.push_back(&exceptionParams.getDevicePointer());
        if (nonbondedMethod != GayBerneForce::NoCutoff) {
            forceArgs.push_back(&maxNeighborBlocks);
            forceArgs.push_back(&neighbors.getDevicePointer());
            forceArgs.push_back(&neighborIndex.getDevicePointer());
            forceArgs.push_back(&neighborBlockCount.getDevicePointer());
            forceArgs.push_back(cu.getPeriodicBoxSizePointer());
            forceArgs.push_back(cu.getInvPeriodicBoxSizePointer());
            forceArgs.push_back(cu.getPeriodicBoxVecXPointer());
            forceArgs.push_back(cu.getPeriodicBoxVecYPointer());
            forceArgs.push_back(cu.getPeriodicBoxVecZPointer());
        }
        torqueArgs.push_back(&cu.getForce().getDevicePointer());
        torqueArgs.push_back(&torque.getDevicePointer());
        torqueArgs.push_back(&numRealParticles);
        torqueArgs.push_back(&cu.getPosq().getDevicePointer());
        torqueArgs.push_back(&axisParticleIndices.getDevicePointer());
        torqueArgs.push_back(&sortedParticles.getDevicePointer());
    }
    cu.executeKernel(framesKernel, &framesArgs[0], numRealParticles);
    cu.executeKernel(blockBoundsKernel, &blockBoundsArgs[0], (numRealParticles+31)/32);
    if (nonbondedMethod == GayBerneForce::NoCutoff) {
        cu.executeKernel(forceKernel, &forceArgs[0], cu.getNonbondedUtilities().getNumForceThreadBlocks()*cu.getNonbondedUtilities().getForceThreadBlockSize());
    }
    else {
        while (true) {
            cu.executeKernel(neighborsKernel, &neighborsArgs[0], numRealParticles);
            int* count = (int*) cu.getPinnedBuffer();
            neighborBlockCount.download(count, false);
            CHECK_RESULT(cuEventRecord(event, 0), "Error recording event for GayBerneForce");
            cu.executeKernel(forceKernel, &forceArgs[0], cu.getNonbondedUtilities().getNumForceThreadBlocks()*cu.getNonbondedUtilities().getForceThreadBlockSize());
            CHECK_RESULT(cuEventSynchronize(event), "Error synchronizing on event for GayBerneForce");
            if (*count <= maxNeighborBlocks)
                break;
            
            // There wasn't enough room for the neighbor list, so we need to recreate it.

            maxNeighborBlocks = (int) ceil((*count)*1.1);
            neighbors.resize(maxNeighborBlocks*32);
            neighborIndex.resize(maxNeighborBlocks);
            neighborsArgs[10] = &neighbors.getDevicePointer();
            neighborsArgs[11] = &neighborIndex.getDevicePointer();
            forceArgs[17] = &neighbors.getDevicePointer();
            forceArgs[18] = &neighborIndex.getDevicePointer();
        }
    }
    cu.executeKernel(torqueKernel, &torqueArgs[0], numRealParticles);
    return 0.0;
}

void CudaCalcGayBerneForceKernel::copyParametersToContext(ContextImpl& context, const GayBerneForce& force) {
    // Make sure the new parameters are acceptable.
    
    if (force.getNumParticles() != cu.getNumAtoms())
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
    
    vector<float4> sigParamsVector(cu.getPaddedNumAtoms(), make_float4(0, 0, 0, 0));
    vector<float2> epsParamsVector(cu.getPaddedNumAtoms(), make_float2(0, 0));
    vector<float4> scaleVector(cu.getPaddedNumAtoms(), make_float4(0, 0, 0, 0));
    for (int i = 0; i < force.getNumParticles(); i++) {
        int xparticle, yparticle;
        double sigma, epsilon, sx, sy, sz, ex, ey, ez;
        force.getParticleParameters(i, sigma, epsilon, xparticle, yparticle, sx, sy, sz, ex, ey, ez);
        sigParamsVector[i] = make_float4((float) (0.5*sigma), (float) (0.25*sx*sx), (float) (0.25*sy*sy), (float) (0.25*sz*sz));
        epsParamsVector[i] = make_float2((float) sqrt(epsilon), (float) (0.125*(sx*sy + sz*sz)*sqrt(sx*sy)));
        scaleVector[i] = make_float4((float) (1/sqrt(ex)), (float) (1/sqrt(ey)), (float) (1/sqrt(ez)), 0);
        if (epsilon != 0.0 && !isRealParticle[i])
            throw OpenMMException("updateParametersInContext: The set of ignored particles (ones with epsilon=0) has changed");
    }
    sigParams.upload(sigParamsVector);
    epsParams.upload(epsParamsVector);
    scale.upload(scaleVector);
    
    // Record the exceptions.
    
    if (numExceptions > 0) {
        vector<float2> exceptionParamsVec(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            int atom1, atom2;
            double sigma, epsilon;
            force.getExceptionParameters(exceptions[i], atom1, atom2, sigma, epsilon);
            exceptionParamsVec[i] = make_float2((float) sigma, (float) epsilon);
        }
        exceptionParams.upload(exceptionParamsVec);
    }
    cu.invalidateMolecules();
    sortAtoms();
}

void CudaCalcGayBerneForceKernel::sortAtoms() {
    // Sort the list of atoms by type to avoid thread divergence.  This is executed every time
    // the atoms are reordered.
    
    int nextIndex = 0;
    vector<int> particles(cu.getPaddedNumAtoms(), 0);
    const vector<int>& order = cu.getAtomIndex();
    vector<int> inverseOrder(order.size(), -1);
    for (int i = 0; i < cu.getNumAtoms(); i++) {
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
        vector<int4> exceptionParticlesVec(numExceptions);
        for (int i = 0; i < numExceptions; i++)
            exceptionParticlesVec[i] = make_int4(exceptionAtoms[i].first, exceptionAtoms[i].second, inverseOrder[exceptionAtoms[i].first], inverseOrder[exceptionAtoms[i].second]);
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

class CudaCalcCustomCVForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(CudaForceInfo& force) : force(force) {
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
    CudaForceInfo& force;
};

class CudaCalcCustomCVForceKernel::ReorderListener : public CudaContext::ReorderListener {
public:
    ReorderListener(CudaContext& cu, CudaArray& invAtomOrder) : cu(cu), invAtomOrder(invAtomOrder) {
    }
    void execute() {
        vector<int> invOrder(cu.getPaddedNumAtoms());
        const vector<int>& order = cu.getAtomIndex();
        for (int i = 0; i < order.size(); i++)
            invOrder[order[i]] = i;
        invAtomOrder.upload(invOrder);
    }
private:
    CudaContext& cu;
    CudaArray& invAtomOrder;
};

void CudaCalcCustomCVForceKernel::initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext) {
    int numCVs = force.getNumCollectiveVariables();
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < numCVs; i++)
        variableNames.push_back(force.getCollectiveVariableName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string name = force.getEnergyParameterDerivativeName(i);
        paramDerivNames.push_back(name);
        cu.addEnergyParameterDerivative(name);
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

    CudaContext& cu2 = *reinterpret_cast<CudaPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    for (auto& param : cu2.getEnergyParamDerivNames())
        cu.addEnergyParameterDerivative(param);
    
    // Create arrays for storing information.
    
    int elementSize = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    cvForces.resize(numCVs);
    for (int i = 0; i < numCVs; i++)
        cvForces[i].initialize<long long>(cu, 3*cu.getPaddedNumAtoms(), "cvForce");
    invAtomOrder.initialize<int>(cu, cu.getPaddedNumAtoms(), "invAtomOrder");
    innerInvAtomOrder.initialize<int>(cu, cu.getPaddedNumAtoms(), "innerInvAtomOrder");
    
    // Create the kernels.
    
    stringstream args, add;
    for (int i = 0; i < numCVs; i++) {
        args << ", long long* __restrict__ force" << i << ", real dEdV" << i;
        add << "forces[i] += (long long) (force" << i << "[i]*dEdV" << i << ");\n";
    }
    map<string, string> replacements;
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    replacements["ADD_FORCES"] = add.str();
    CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaKernelSources::customCVForce, replacements));
    copyStateKernel = cu.getKernel(module, "copyState");
    copyForcesKernel = cu.getKernel(module, "copyForces");
    addForcesKernel = cu.getKernel(module, "addForces");

    // This context needs to respect all forces in the inner context when reordering atoms.

    for (CudaForceInfo* info : cu2.getForceInfos())
        cu.addForce(new ForceInfo(*info));
}

double CudaCalcCustomCVForceKernel::execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy) {
    copyState(context, innerContext);
    int numCVs = variableNames.size();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    CudaContext& cu2 = *reinterpret_cast<CudaPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    vector<double> cvValues;
    vector<map<string, double> > cvDerivs(numCVs);
    void* copyForcesArgs[] = {NULL, &invAtomOrder.getDevicePointer(), &cu2.getForce().getDevicePointer(), &cu2.getAtomIndexArray().getDevicePointer(), &numAtoms, &paddedNumAtoms};
    for (int i = 0; i < numCVs; i++) {
        cvValues.push_back(innerContext.calcForcesAndEnergy(true, true, 1<<i));
        copyForcesArgs[0] = &cvForces[i].getDevicePointer();
        cu.executeKernel(copyForcesKernel, copyForcesArgs, numAtoms);
        innerContext.getEnergyParameterDerivatives(cvDerivs[i]);
    }
    
    // Compute the energy and forces.
    
    map<string, double> variables;
    for (auto& name : globalParameterNames)
        variables[name] = context.getParameter(name);
    for (int i = 0; i < numCVs; i++)
        variables[variableNames[i]] = cvValues[i];
    double energy = energyExpression.evaluate(variables);
    int bufferSize = cu.getForce().getSize();
    vector<void*> addForcesArgs;
    addForcesArgs.push_back(&cu.getForce().getDevicePointer());
    addForcesArgs.push_back(&bufferSize);
    vector<double> dEdV(numCVs);
    vector<float> dEdVFloat(numCVs);
    for (int i = 0; i < numCVs; i++) {
        dEdV[i] = variableDerivExpressions[i].evaluate(variables);
        dEdVFloat[i] = (float) dEdV[i];
        addForcesArgs.push_back(&cvForces[i].getDevicePointer());
        if (cu.getUseDoublePrecision())
            addForcesArgs.push_back(&dEdV[i]);
        else
            addForcesArgs.push_back(&dEdVFloat[i]);
    }
    cu.executeKernel(addForcesKernel, &addForcesArgs[0], numAtoms);
    
    // Compute the energy parameter derivatives.
    
    map<string, double>& energyParamDerivs = cu.getEnergyParamDerivWorkspace();
    for (int i = 0; i < paramDerivExpressions.size(); i++)
        energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate(variables);
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate(variables);
        for (auto& deriv : cvDerivs[i])
            energyParamDerivs[deriv.first] += dEdV*deriv.second;
    }
    return energy;
}

void CudaCalcCustomCVForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext) {
    int numAtoms = cu.getNumAtoms();
    CudaContext& cu2 = *reinterpret_cast<CudaPlatform::PlatformData*>(innerContext.getPlatformData())->contexts[0];
    if (!hasInitializedListeners) {
        hasInitializedListeners = true;
        
        // Initialize the listeners.
        
        ReorderListener* listener1 = new ReorderListener(cu, invAtomOrder);
        ReorderListener* listener2 = new ReorderListener(cu2, innerInvAtomOrder);
        cu.addReorderListener(listener1);
        cu2.addReorderListener(listener2);
        listener1->execute();
        listener2->execute();
    }
    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    CUdeviceptr posCorrection2 = (cu2.getUseMixedPrecision() ? cu2.getPosqCorrection().getDevicePointer() : 0);
    void* copyStateArgs[] = {&cu.getPosq().getDevicePointer(), &posCorrection, &cu.getVelm().getDevicePointer(), &cu.getAtomIndexArray().getDevicePointer(),
        &cu2.getPosq().getDevicePointer(), &posCorrection2,& cu2.getVelm().getDevicePointer(), &innerInvAtomOrder.getDevicePointer(), &numAtoms};
    cu.executeKernel(copyStateKernel, copyStateArgs, numAtoms);
    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext.setPeriodicBoxVectors(a, b, c);
    innerContext.setTime(context.getTime());
    map<string, double> innerParameters = innerContext.getParameters();
    for (auto& param : innerParameters)
        innerContext.setParameter(param.first, context.getParameter(param.first));
}

void CudaCalcCustomCVForceKernel::copyParametersToContext(ContextImpl& context, const CustomCVForce& force) {
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

class CudaCalcRMSDForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const RMSDForce& force) : force(force) {
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

void CudaCalcRMSDForceKernel::initialize(const System& system, const RMSDForce& force) {
    // Create data structures.
    
    bool useDouble = cu.getUseDoublePrecision();
    int elementSize = (useDouble ? sizeof(double) : sizeof(float));
    int numParticles = force.getParticles().size();
    if (numParticles == 0)
        numParticles = system.getNumParticles();
    referencePos.initialize(cu, system.getNumParticles(), 4*elementSize, "referencePos");
    particles.initialize<int>(cu, numParticles, "particles");
    buffer.initialize(cu, 13, elementSize, "buffer");
    recordParameters(force);
    info = new ForceInfo(force);
    cu.addForce(info);
    
    // Create the kernels.

    CUmodule module = cu.createModule(CudaKernelSources::vectorOps+CudaKernelSources::rmsd);
    kernel1 = cu.getKernel(module, "computeRMSDPart1");
    kernel2 = cu.getKernel(module, "computeRMSDForces");
}

void CudaCalcRMSDForceKernel::recordParameters(const RMSDForce& force) {
    // Record the parameters and center the reference positions.
    
    vector<int> particleVec = force.getParticles();
    if (particleVec.size() == 0)
        for (int i = 0; i < cu.getNumAtoms(); i++)
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
    vector<double4> pos;
    for (Vec3 p : centeredPositions)
        pos.push_back(make_double4(p[0], p[1], p[2], 0));
    referencePos.upload(pos, true);

    // Record the sum of the norms of the reference positions.

    sumNormRef = 0.0;
    for (int i : particleVec) {
        Vec3 p = centeredPositions[i];
        sumNormRef += p.dot(p);
    }
}

double CudaCalcRMSDForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (cu.getUseDoublePrecision())
        return executeImpl<double>(context);
    return executeImpl<float>(context);
}

template <class REAL>
double CudaCalcRMSDForceKernel::executeImpl(ContextImpl& context) {
    // Execute the first kernel.

    int numParticles = particles.getSize();
    int blockSize = 256;
    void* args1[] = {&numParticles, &cu.getPosq().getDevicePointer(), &referencePos.getDevicePointer(),
            &particles.getDevicePointer(), &buffer.getDevicePointer()};
    cu.executeKernel(kernel1, args1, blockSize, blockSize, blockSize*sizeof(REAL));
    
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
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    void* args2[] = {&numParticles, &paddedNumAtoms, &cu.getPosq().getDevicePointer(), &referencePos.getDevicePointer(),
            &particles.getDevicePointer(), &buffer.getDevicePointer(), &cu.getForce().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numParticles);
    return rmsd;
}

void CudaCalcRMSDForceKernel::copyParametersToContext(ContextImpl& context, const RMSDForce& force) {
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
    cu.invalidateMolecules(info);
}

void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::verlet, defines, "");
    kernel1 = cu.getKernel(module, "integrateVerletPart1");
    kernel2 = cu.getKernel(module, "integrateVerletPart2");
}

void CudaIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    cu.getIntegrationUtilities().setNextStepSize(dt);

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void CudaIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::langevin, defines, "");
    kernel1 = cu.getKernel(module, "integrateLangevinPart1");
    kernel2 = cu.getKernel(module, "integrateLangevinPart2");
    params.initialize(cu, 3, cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float), "langevinParams");
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
    cu.getIntegrationUtilities().setNextStepSize(stepSize);
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize*friction);
        double fscale = (friction == 0 ? stepSize : (1-vscale)/friction);
        double noisescale = sqrt(kT*(1-vscale*vscale));
        vector<double> p(params.getSize());
        p[0] = vscale;
        p[1] = fscale;
        p[2] = noisescale;
        params.upload(p, true);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Call the first integration kernel.

    int randomIndex = integration.prepareRandomNumbers(cu.getPaddedNumAtoms());
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &params.getDevicePointer(), &integration.getStepSize().getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &posCorrection, &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getStepSize().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+stepSize);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
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
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&numAtoms, useDouble ? (void*) &stepSize : (void*) &stepSizeFloat,
            &cu.getPosq().getDevicePointer(), &posCorrection, &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+stepSize);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateBrownianStepKernel::computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0);
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
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    double dt = cu.getIntegrationUtilities().getLastStepSize();
    double time = cu.getTime()+dt;
    if (useDouble) {
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
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

void CudaIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaKernelSources::langevin, defines, "");
    kernel1 = cu.getKernel(module, "integrateLangevinPart1");
    kernel2 = cu.getKernel(module, "integrateLangevinPart2");
    selectSizeKernel = cu.getKernel(module, "selectLangevinStepSize");
    params.initialize(cu, 3, cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float), "langevinParams");
    blockSize = min(256, system.getNumParticles());
    blockSize = max(blockSize, params.getSize());
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
    double friction = integrator.getFriction();
    float frictionFloat = (float) friction;
    double kT = BOLTZ*integrator.getTemperature();
    float kTFloat = (float) kT;
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();
    void* argsSelect[] = {&numAtoms, &paddedNumAtoms, useDouble ? (void*) &maxStepSize : (void*) &maxStepSizeFloat,
            useDouble ? (void*) &tol : (void*) &tolFloat,
            useDouble ? (void*) &friction : (void*) &frictionFloat,
            useDouble ? (void*) &kT : (void*) &kTFloat,
            &cu.getIntegrationUtilities().getStepSize().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &params.getDevicePointer()};
    int sharedSize = 2*blockSize*(useDouble ? sizeof(double) : sizeof(float));
    cu.executeKernel(selectSizeKernel, argsSelect, blockSize, blockSize, sharedSize);

    // Call the first integration kernel.

    int randomIndex = integration.prepareRandomNumbers(cu.getPaddedNumAtoms());
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer(),
            &params.getDevicePointer(), &integration.getStepSize().getDevicePointer(), &integration.getRandom().getDevicePointer(), &randomIndex};
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args2[] = {&numAtoms, &cu.getPosq().getDevicePointer(), &posCorrection, &integration.getPosDelta().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getStepSize().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    double dt = cu.getIntegrationUtilities().getLastStepSize();
    double time = cu.getTime()+dt;
    if (useDouble) {
        if (dt == maxStepSize)
            time = maxTime; // Avoid round-off error
    }
    else {
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
    ReorderListener(CudaContext& cu, vector<CudaArray>& perDofValues, vector<vector<float4> >& localPerDofValuesFloat, vector<vector<double4> >& localPerDofValuesDouble, vector<bool>& deviceValuesAreCurrent) :
            cu(cu), perDofValues(perDofValues), localPerDofValuesFloat(localPerDofValuesFloat), localPerDofValuesDouble(localPerDofValuesDouble), deviceValuesAreCurrent(deviceValuesAreCurrent) {
        int numAtoms = cu.getNumAtoms();
        lastAtomOrder.resize(numAtoms);
        for (int i = 0; i < numAtoms; i++)
            lastAtomOrder[i] = cu.getAtomIndex()[i];
    }
    void execute() {
        // Reorder the per-DOF variables to reflect the new atom order.

        if (perDofValues.size() == 0)
            return;
        int numAtoms = cu.getNumAtoms();
        const vector<int>& order = cu.getAtomIndex();
        for (int index = 0; index < perDofValues.size(); index++) {
            if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
                if (deviceValuesAreCurrent[index])
                    perDofValues[index].download(localPerDofValuesDouble[index]);
                vector<double4> swap(numAtoms);
                for (int i = 0; i < numAtoms; i++)
                    swap[lastAtomOrder[i]] = localPerDofValuesDouble[index][i];
                for (int i = 0; i < numAtoms; i++)
                    localPerDofValuesDouble[index][i] = swap[order[i]];
                perDofValues[index].upload(localPerDofValuesDouble[index]);
            }
            else {
                if (deviceValuesAreCurrent[index])
                    perDofValues[index].download(localPerDofValuesFloat[index]);
                vector<float4> swap(numAtoms);
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
    CudaContext& cu;
    vector<CudaArray>& perDofValues;
    vector<vector<float4> >& localPerDofValuesFloat;
    vector<vector<double4> >& localPerDofValuesDouble;
    vector<bool>& deviceValuesAreCurrent;
    vector<int> lastAtomOrder;
};

class CudaIntegrateCustomStepKernel::DerivFunction : public CustomFunction {
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

void CudaIntegrateCustomStepKernel::initialize(const System& system, const CustomIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    numGlobalVariables = integrator.getNumGlobalVariables();
    int elementSize = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
    sumBuffer.initialize(cu, system.getNumParticles(), elementSize, "sumBuffer");
    summedValue.initialize(cu, 1, elementSize, "summedValue");
    perDofValues.resize(integrator.getNumPerDofVariables());
    localPerDofValuesFloat.resize(perDofValues.size());
    localPerDofValuesDouble.resize(perDofValues.size());
    for (int i = 0; i < perDofValues.size(); i++)
        perDofValues[i].initialize(cu, system.getNumParticles(), 4*elementSize, "perDofVariables");
    localValuesAreCurrent.resize(integrator.getNumPerDofVariables(), false);
    deviceValuesAreCurrent.resize(integrator.getNumPerDofVariables(), false);
    cu.addReorderListener(new ReorderListener(cu, perDofValues, localPerDofValuesFloat, localPerDofValuesDouble, deviceValuesAreCurrent));
    SimTKOpenMMUtilities::setRandomNumberSeed(integrator.getRandomNumberSeed());
}

string CudaIntegrateCustomStepKernel::createPerDofComputation(const string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator,
        const string& forceName, const string& energyName, vector<const TabulatedFunction*>& functions, vector<pair<string, string> >& functionNames) {
    map<string, Lepton::ParsedExpression> expressions;
    expressions["double3 tempResult = "] = expr;
    map<string, string> variables;
    variables["x"] = "trimTo3(position)";
    variables["v"] = "trimTo3(velocity)";
    variables[forceName] = "trimTo3(f)";
    variables["gaussian"] = "trimTo3(gaussian)";
    variables["uniform"] = "trimTo3(uniform)";
    variables["m"] = "mass";
    variables["dt"] = "stepSize";
    if (energyName != "")
        variables[energyName] = "make_double3(energy)";
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        variables[integrator.getGlobalVariableName(i)] = "make_double3(globals["+cu.intToString(globalVariableIndex[i])+"])";
    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
        variables[integrator.getPerDofVariableName(i)] = "perDof"+cu.intToString(i);
    for (int i = 0; i < (int) parameterNames.size(); i++)
        variables[parameterNames[i]] = "make_double3(globals["+cu.intToString(parameterVariableIndex[i])+"])";
    vector<pair<ExpressionTreeNode, string> > variableNodes;
    findExpressionsForDerivs(expr.getRootNode(), variableNodes);
    for (auto& var : variables)
        variableNodes.push_back(make_pair(ExpressionTreeNode(new Operation::Variable(var.first)), var.second));
    string result = cu.getExpressionUtilities().createExpressions(expressions, variableNodes, functions, functionNames, "temp", "double3");
    if (variable == "x")
        result += "position.x = tempResult.x; position.y = tempResult.y; position.z = tempResult.z;\n";
    else if (variable == "v")
        result += "velocity.x = tempResult.x; velocity.y = tempResult.y; velocity.z = tempResult.z;\n";
    else if (variable == "")
        result += "sum[index] = tempResult.x+tempResult.y+tempResult.z;\n";
    else {
        for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
            if (variable == integrator.getPerDofVariableName(i)) {
                string varName = "perDof"+cu.intToString(i);
                result += varName+".x = tempResult.x; "+varName+".y = tempResult.y; "+varName+".z = tempResult.z;\n";
            }
    }
    return result;
}

void CudaIntegrateCustomStepKernel::prepareForComputation(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();
    string perDofType = (useDouble ? "double4" : "float4");
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        
        // Initialize various data structures.
        
        const map<string, double>& params = context.getParameters();
        for (auto& param : params)
            parameterNames.push_back(param.first);
        kernels.resize(integrator.getNumComputations());
        kernelArgs.resize(integrator.getNumComputations());
        requiredGaussian.resize(integrator.getNumComputations(), 0);
        requiredUniform.resize(integrator.getNumComputations(), 0);
        needsGlobals.resize(numSteps, false);
        globalExpressions.resize(numSteps);
        stepType.resize(numSteps);
        stepTarget.resize(numSteps);
        merged.resize(numSteps, false);
        modifiesParameters = false;
        sumWorkGroupSize = 512;
        map<string, string> defines;
        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
        defines["WORK_GROUP_SIZE"] = cu.intToString(sumWorkGroupSize);
        defines["SUM_BUFFER_SIZE"] = "0";

        // Record the tabulated functions.

        map<string, Lepton::CustomFunction*> functions;
        vector<pair<string, string> > functionNames;
        vector<const TabulatedFunction*> functionList;
        vector<string> tableTypes;
        tabulatedFunctions.resize(integrator.getNumTabulatedFunctions());
        for (int i = 0; i < integrator.getNumTabulatedFunctions(); i++) {
            functionList.push_back(&integrator.getTabulatedFunction(i));
            string name = integrator.getTabulatedFunctionName(i);
            string arrayName = "table"+cu.intToString(i);
            functionNames.push_back(make_pair(name, arrayName));
            functions[name] = createReferenceTabulatedFunction(integrator.getTabulatedFunction(i));
            int width;
            vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(integrator.getTabulatedFunction(i), width);
            tabulatedFunctions[i].initialize<float>(cu, f.size(), "TabulatedFunction");
            tabulatedFunctions[i].upload(f);
            if (width == 1)
                tableTypes.push_back("float");
            else
                tableTypes.push_back("float"+cu.intToString(width));
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
                savedForces[forceGroupFlags[step]] = CudaArray();
                savedForces[forceGroupFlags[step]].initialize(cu, cu.getForce().getSize(), cu.getForce().getElementSize(), "savedForces");
            }
        }
        
        // Allocate space for storing global values, both on the host and the device.
        
        localGlobalValues.resize(expressionSet.getNumVariables());
        int elementSize = (cu.getUseDoublePrecision() || cu.getUseMixedPrecision() ? sizeof(double) : sizeof(float));
        globalValues.initialize(cu, expressionSet.getNumVariables(), elementSize, "globalValues");
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
        perDofEnergyParamDerivs.initialize(cu, max(1, numContextParams), elementSize, "perDofEnergyParamDerivs");
        
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
                    compute << "double3 perDof"<<cu.intToString(i)<<" = trimTo3(convertToDouble4(perDofValues"<<cu.intToString(i)<<"[index]));\n";
                int numGaussian = 0, numUniform = 0;
                for (int j = step; j < numSteps && (j == step || merged[j]); j++) {
                    numGaussian += numAtoms*usesVariable(expression[j][0], "gaussian");
                    numUniform += numAtoms*usesVariable(expression[j][0], "uniform");
                    compute << "{\n";
                    if (numGaussian > 0)
                        compute << "double4 gaussian = convertToDouble4(gaussianValues[gaussianIndex+index]);\n";
                    if (numUniform > 0)
                        compute << "double4 uniform = convertToDouble4(uniformValues[uniformIndex+index]);\n";
                    compute << createPerDofComputation(stepType[j] == CustomIntegrator::ComputePerDof ? variable[j] : "", expression[j][0], integrator, forceName[j], energyName[j], functionList, functionNames);
                    if (variable[j] == "x") {
                        if (storePosAsDelta[j])
                            compute << "posDelta[index] = convertFromDouble4(position-convertToDouble4(loadPos(posq, posqCorrection, index)));\n";
                        else
                            compute << "storePos(posq, posqCorrection, index, convertFromDouble4(position));\n";
                    }
                    else if (variable[j] == "v")
                        compute << "velm[index] = convertFromDouble4(velocity);\n";
                    else {
                        for (int i = 0; i < perDofValues.size(); i++)
                            compute << "perDofValues"<<cu.intToString(i)<<"[index] = make_"<<perDofType<<"(perDof"<<cu.intToString(i)<<".x, perDof"<<cu.intToString(i)<<".y, perDof"<<cu.intToString(i)<<".z, 0);\n";
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
                    string valueName = "perDofValues"+cu.intToString(i);
                    args << ", " << perDofType << "* __restrict__ " << valueName;
                }
                for (int i = 0; i < (int) tableTypes.size(); i++)
                    args << ", const " << tableTypes[i]<< "* __restrict__ table" << i;
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
                args1.push_back(&globalValues.getDevicePointer());
                args1.push_back(&sumBuffer.getDevicePointer());
                args1.push_back(NULL);
                args1.push_back(NULL);
                args1.push_back(NULL);
                if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
                    args1.push_back(&energy);
                else
                    args1.push_back(&energyFloat);
                args1.push_back(&perDofEnergyParamDerivs.getDevicePointer());
                for (auto& array : perDofValues)
                    args1.push_back(&array.getDevicePointer());
                for (auto& array : tabulatedFunctions)
                    args1.push_back(&array.getDevicePointer());
                kernelArgs[step].push_back(args1);
                if (stepType[step] == CustomIntegrator::ComputeSum) {
                    // Create a second kernel for this step that sums the values.

                    vector<void*> args2;
                    args2.push_back(&sumBuffer.getDevicePointer());
                    args2.push_back(&summedValue.getDevicePointer());
                    defines["SUM_BUFFER_SIZE"] = cu.intToString(numAtoms);
                    module = cu.createModule(CudaKernelSources::customIntegrator, defines);
                    kernel = cu.getKernel(module, useDouble ? "computeDoubleSum" : "computeFloatSum");
                    kernels[step].push_back(kernel);
                    kernelArgs[step].push_back(args2);
                }
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
        for (int required : requiredUniform)
            maxUniformRandoms = max(maxUniformRandoms, required);
        uniformRandoms.initialize<float4>(cu, maxUniformRandoms, "uniformRandoms");
        randomSeed.initialize<int4>(cu, cu.getNumThreadBlocks()*CudaContext::ThreadBlockSize, "randomSeed");
        vector<int4> seed(randomSeed.getSize());
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
        CUmodule randomProgram = cu.createModule(CudaKernelSources::customIntegrator, defines);
        randomKernel = cu.getKernel(randomProgram, "generateRandomNumbers");
        
        // Create the kernel for computing kinetic energy.

        stringstream computeKE;
        for (int i = 0; i < perDofValues.size(); i++)
            computeKE << "double3 perDof"<<cu.intToString(i)<<" = trimTo3(convertToDouble4(perDofValues"<<cu.intToString(i)<<"[index]));\n";
        Lepton::ParsedExpression keExpression = Lepton::Parser::parse(integrator.getKineticEnergyExpression()).optimize();
        computeKE << createPerDofComputation("", keExpression, integrator, "f", "", functionList, functionNames);
        map<string, string> replacements;
        replacements["COMPUTE_STEP"] = computeKE.str();
        stringstream args;
        for (int i = 0; i < perDofValues.size(); i++) {
            string valueName = "perDofValues"+cu.intToString(i);
            args << ", " << perDofType << "* __restrict__ " << valueName;
        }
        for (int i = 0; i < (int) tableTypes.size(); i++)
            args << ", const " << tableTypes[i]<< "* __restrict__ table" << i;
        replacements["PARAMETER_ARGUMENTS"] = args.str();
        defines["SUM_BUFFER_SIZE"] = cu.intToString(numAtoms);
        if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
            defines.erase("LOAD_POS_AS_DELTA");
        CUmodule module = cu.createModule(cu.replaceStrings(CudaKernelSources::vectorOps+CudaKernelSources::customIntegratorPerDof, replacements), defines);
        kineticEnergyKernel = cu.getKernel(module, "computePerDof");
        kineticEnergyArgs.push_back(&cu.getPosq().getDevicePointer());
        kineticEnergyArgs.push_back(NULL);
        kineticEnergyArgs.push_back(&integration.getPosDelta().getDevicePointer());
        kineticEnergyArgs.push_back(&cu.getVelm().getDevicePointer());
        kineticEnergyArgs.push_back(&cu.getForce().getDevicePointer());
        kineticEnergyArgs.push_back(&integration.getStepSize().getDevicePointer());
        kineticEnergyArgs.push_back(&globalValues.getDevicePointer());
        kineticEnergyArgs.push_back(&sumBuffer.getDevicePointer());
        kineticEnergyArgs.push_back(NULL);
        kineticEnergyArgs.push_back(NULL);
        kineticEnergyArgs.push_back(&uniformRandoms.getDevicePointer());
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
            kineticEnergyArgs.push_back(&energy);
        else
            kineticEnergyArgs.push_back(&energyFloat);
        kineticEnergyArgs.push_back(&perDofEnergyParamDerivs.getDevicePointer());
        for (auto& array : perDofValues)
            kineticEnergyArgs.push_back(&array.getDevicePointer());
        for (auto& array : tabulatedFunctions)
            kineticEnergyArgs.push_back(&array.getDevicePointer());
        keNeedsForce = usesVariable(keExpression, "f");

        // Create a second kernel to sum the values.

        defines["SUM_BUFFER_SIZE"] = cu.intToString(numAtoms);
        module = cu.createModule(CudaKernelSources::customIntegrator, defines);
        sumKineticEnergyKernel = cu.getKernel(module, useDouble ? "computeDoubleSum" : "computeFloatSum");

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

ExpressionTreeNode CudaIntegrateCustomStepKernel::replaceDerivFunctions(const ExpressionTreeNode& node, ContextImpl& context) {
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

void CudaIntegrateCustomStepKernel::findExpressionsForDerivs(const ExpressionTreeNode& node, vector<pair<ExpressionTreeNode, string> >& variableNodes) {
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
        variableNodes.push_back(make_pair(node, "make_double3(energyParamDerivs["+cu.intToString(index)+"])"));
        needsEnergyParamDerivs = true;
    }
    else {
        for (auto& child : node.getChildren())
            findExpressionsForDerivs(child, variableNodes);
    }
}

void CudaIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    prepareForComputation(context, integrator, forcesAreValid);
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int numSteps = integrator.getNumComputations();
    if (!forcesAreValid)
        savedEnergy.clear();

    // Loop over computation steps in the integrator and execute them.

    int maxUniformRandoms = uniformRandoms.getSize();
    void* randomArgs[] = {&maxUniformRandoms, &uniformRandoms.getDevicePointer(), &randomSeed.getDevicePointer()};
    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
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

                    cu.getForce().copyTo(savedForces[lastForceGroups]);
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
                
                savedForces[forceGroups].copyTo(cu.getForce());
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
                        perDofEnergyParamDerivs.upload(localPerDofEnergyParamDerivs, true);
                    }
                }
            }
            forcesAreValid = true;
        }
        if (needsEnergy[step]) {
            energy = savedEnergy[forceGroups];
            energyFloat = (float) energy;
        }
        if (needsGlobals[step] && !deviceGlobalsAreCurrent) {
            // Upload the global values to the device.
            
            globalValues.upload(localGlobalValues, true);
            deviceGlobalsAreCurrent = true;
        }
        bool stepInvalidatesForces = invalidatesForces[step];
        if (stepType[step] == CustomIntegrator::ComputePerDof && !merged[step]) {
            int randomIndex = integration.prepareRandomNumbers(requiredGaussian[step]);
            kernelArgs[step][0][1] = &posCorrection;
            kernelArgs[step][0][8] = &integration.getRandom().getDevicePointer();
            kernelArgs[step][0][9] = &randomIndex;
            kernelArgs[step][0][10] = &uniformRandoms.getDevicePointer();
            if (requiredUniform[step] > 0)
                cu.executeKernel(randomKernel, &randomArgs[0], numAtoms);
            cu.executeKernel(kernels[step][0], &kernelArgs[step][0][0], numAtoms, 128);
        }
        else if (stepType[step] == CustomIntegrator::ComputeGlobal) {
            expressionSet.setVariable(uniformVariableIndex, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
            expressionSet.setVariable(gaussianVariableIndex, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
            expressionSet.setVariable(stepEnergyVariableIndex[step], energy);
            recordGlobalValue(globalExpressions[step][0].evaluate(), stepTarget[step], integrator);
        }
        else if (stepType[step] == CustomIntegrator::ComputeSum) {
            int randomIndex = integration.prepareRandomNumbers(requiredGaussian[step]);
            kernelArgs[step][0][1] = &posCorrection;
            kernelArgs[step][0][8] = &integration.getRandom().getDevicePointer();
            kernelArgs[step][0][9] = &randomIndex;
            kernelArgs[step][0][10] = &uniformRandoms.getDevicePointer();
            if (requiredUniform[step] > 0)
                cu.executeKernel(randomKernel, &randomArgs[0], numAtoms);
            cu.clearBuffer(sumBuffer);
            cu.executeKernel(kernels[step][0], &kernelArgs[step][0][0], numAtoms, 128);
            cu.executeKernel(kernels[step][1], &kernelArgs[step][1][0], sumWorkGroupSize, sumWorkGroupSize);
            if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
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
                cu.getIntegrationUtilities().applyConstraints(integrator.getConstraintTolerance());
                kernelArgs[step][0][1] = &posCorrection;
                cu.executeKernel(kernels[step][0], &kernelArgs[step][0][0], numAtoms);
            }
            cu.getIntegrationUtilities().computeVirtualSites();
        }
        else if (stepType[step] == CustomIntegrator::ConstrainVelocities) {
            cu.getIntegrationUtilities().applyVelocityConstraints(integrator.getConstraintTolerance());
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

    cu.setTime(cu.getTime()+integrator.getStepSize());
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
    if (cu.getAtomsWereReordered()) {
        forcesAreValid = false;
        validSavedForces.clear();
    }
}

bool CudaIntegrateCustomStepKernel::evaluateCondition(int step) {
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

double CudaIntegrateCustomStepKernel::computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    prepareForComputation(context, integrator, forcesAreValid);
    if (keNeedsForce && !forcesAreValid) {
        // Compute the force.  We want to then mark that forces are valid, which means also computing
        // potential energy if any steps will expect it to be valid too.
        
        bool willNeedEnergy = false;
        for (int i = 0; i < integrator.getNumComputations(); i++)
            willNeedEnergy |= needsEnergy[i];
        energy = context.calcForcesAndEnergy(true, willNeedEnergy, -1);
        energyFloat = (float) energy;
        forcesAreValid = true;
    }
    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    int randomIndex = 0;
    kineticEnergyArgs[1] = &posCorrection;
    kineticEnergyArgs[8] = &cu.getIntegrationUtilities().getRandom().getDevicePointer();
    kineticEnergyArgs[9] = &randomIndex;
    cu.clearBuffer(sumBuffer);
    cu.executeKernel(kineticEnergyKernel, &kineticEnergyArgs[0], cu.getNumAtoms());
    void* args[] = {&sumBuffer.getDevicePointer(), &summedValue.getDevicePointer()};
    cu.executeKernel(sumKineticEnergyKernel, args, sumWorkGroupSize, sumWorkGroupSize);
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
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

void CudaIntegrateCustomStepKernel::recordGlobalValue(double value, GlobalTarget target, CustomIntegrator& integrator) {
    switch (target.type) {
        case DT:
            if (value != localGlobalValues[dtVariableIndex])
                deviceGlobalsAreCurrent = false;
            expressionSet.setVariable(dtVariableIndex, value);
            localGlobalValues[dtVariableIndex] = value;
            cu.getIntegrationUtilities().setNextStepSize(value);
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

void CudaIntegrateCustomStepKernel::recordChangedParameters(ContextImpl& context) {
    if (!modifiesParameters)
        return;
    for (int i = 0; i < (int) parameterNames.size(); i++) {
        double value = context.getParameter(parameterNames[i]);
        if (value != localGlobalValues[parameterVariableIndex[i]])
            context.setParameter(parameterNames[i], localGlobalValues[parameterVariableIndex[i]]);
    }
}

void CudaIntegrateCustomStepKernel::getGlobalVariables(ContextImpl& context, vector<double>& values) const {
    if (!globalValues.isInitialized()) {
        // The data structures haven't been created yet, so just return the list of values that was given earlier.
        
        values = initialGlobalVariables;
        return;
    }
    values.resize(numGlobalVariables);
    for (int i = 0; i < numGlobalVariables; i++)
        values[i] = localGlobalValues[globalVariableIndex[i]];
}

void CudaIntegrateCustomStepKernel::setGlobalVariables(ContextImpl& context, const vector<double>& values) {
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

void CudaIntegrateCustomStepKernel::getPerDofVariable(ContextImpl& context, int variable, vector<Vec3>& values) const {
    values.resize(perDofValues[variable].getSize());
    const vector<int>& order = cu.getAtomIndex();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
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

void CudaIntegrateCustomStepKernel::setPerDofVariable(ContextImpl& context, int variable, const vector<Vec3>& values) {
    const vector<int>& order = cu.getAtomIndex();
    localValuesAreCurrent[variable] = true;
    deviceValuesAreCurrent[variable] = false;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
        localPerDofValuesDouble[variable].resize(values.size());
        for (int i = 0; i < (int) values.size(); i++)
            localPerDofValuesDouble[variable][i] = make_double4(values[order[i]][0], values[order[i]][1], values[order[i]][2], 0);
    }
    else {
        localPerDofValuesFloat[variable].resize(values.size());
        for (int i = 0; i < (int) values.size(); i++)
            localPerDofValuesFloat[variable][i] = make_float4(values[order[i]][0], values[order[i]][1], values[order[i]][2], 0);
    }
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
    atomGroups.initialize<int>(cu, cu.getNumAtoms(), "atomGroups");
    vector<int> atoms(atomGroups.getSize());
    for (int i = 0; i < (int) groups.size(); i++) {
        for (int j = 0; j < (int) groups[i].size(); j++)
            atoms[groups[i][j]] = i;
    }
    atomGroups.upload(atoms);
}

void CudaApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    cu.setAsCurrent();
    float frequency = (float) context.getParameter(AndersenThermostat::CollisionFrequency());
    float kT = (float) (BOLTZ*context.getParameter(AndersenThermostat::Temperature()));
    int randomIndex = cu.getIntegrationUtilities().prepareRandomNumbers(cu.getPaddedNumAtoms());
    int numAtoms = cu.getNumAtoms();
    void* args[] = {&numAtoms, &frequency, &kT, &cu.getVelm().getDevicePointer(), &cu.getIntegrationUtilities().getStepSize().getDevicePointer(),
            &cu.getIntegrationUtilities().getRandom().getDevicePointer(), &randomIndex, &atomGroups.getDevicePointer()};
    cu.executeKernel(kernel, args, cu.getNumAtoms());
}

void CudaApplyMonteCarloBarostatKernel::initialize(const System& system, const Force& thermostat) {
    cu.setAsCurrent();
    savedPositions.initialize(cu, cu.getPaddedNumAtoms(), cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4), "savedPositions");
    savedForces.initialize<long long>(cu, cu.getPaddedNumAtoms()*3, "savedForces");
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
        moleculeAtoms.initialize<int>(cu, cu.getNumAtoms(), "moleculeAtoms");
        moleculeStartIndex.initialize<int>(cu, numMolecules+1, "moleculeStartIndex");
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
        
    }
    int bytesToCopy = cu.getPosq().getSize()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4));
    CUresult result = cuMemcpyDtoD(savedPositions.getDevicePointer(), cu.getPosq().getDevicePointer(), bytesToCopy);
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error saving positions for MC barostat: "<<cu.getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
    result = cuMemcpyDtoD(savedForces.getDevicePointer(), cu.getForce().getDevicePointer(), savedForces.getSize()*savedForces.getElementSize());
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error saving forces for MC barostat: "<<cu.getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
    float scalefX = (float) scaleX;
    float scalefY = (float) scaleY;
    float scalefZ = (float) scaleZ;
    void* args[] = {&scalefX, &scalefY, &scalefZ, &numMolecules, cu.getPeriodicBoxSizePointer(), cu.getInvPeriodicBoxSizePointer(),
                    cu.getPeriodicBoxVecXPointer(), cu.getPeriodicBoxVecYPointer(), cu.getPeriodicBoxVecZPointer(),
		    &cu.getPosq().getDevicePointer(), &moleculeAtoms.getDevicePointer(), &moleculeStartIndex.getDevicePointer()};
    cu.executeKernel(kernel, args, cu.getNumAtoms());
    for (auto& offset : cu.getPosCellOffsets())
        offset = make_int4(0, 0, 0, 0);
    lastAtomOrder = cu.getAtomIndex();
}

void CudaApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    cu.setAsCurrent();
    int bytesToCopy = cu.getPosq().getSize()*(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4));
    CUresult result = cuMemcpyDtoD(cu.getPosq().getDevicePointer(), savedPositions.getDevicePointer(), bytesToCopy);
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error restoring positions for MC barostat: "<<cu.getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
    result = cuMemcpyDtoD(cu.getForce().getDevicePointer(), savedForces.getDevicePointer(), savedForces.getSize()*savedForces.getElementSize());
    if (result != CUDA_SUCCESS) {
        std::stringstream m;
        m<<"Error restoring forces for MC barostat: "<<cu.getErrorString(result)<<" ("<<result<<")";
        throw OpenMMException(m.str());
    }
}

void CudaRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    cu.setAsCurrent();
    frequency = force.getFrequency();
    int numAtoms = cu.getNumAtoms();
    cmMomentum.initialize<float4>(cu, (numAtoms+CudaContext::ThreadBlockSize-1)/CudaContext::ThreadBlockSize, "cmMomentum");
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
    void* args[] = {&numAtoms, &cu.getVelm().getDevicePointer(), &cmMomentum.getDevicePointer()};
    cu.executeKernel(kernel1, args, cu.getNumAtoms(), cu.ThreadBlockSize, cu.ThreadBlockSize*sizeof(float4));
    cu.executeKernel(kernel2, args, cu.getNumAtoms(), cu.ThreadBlockSize, cu.ThreadBlockSize*sizeof(float4));
}
