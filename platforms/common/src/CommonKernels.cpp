/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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
#include "openmm/common/CommonKernelUtilities.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/ExpressionUtilities.h"
#include "openmm/Context.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCentroidBondForceImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/internal/DPDIntegratorUtilities.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/timer.h"
#include "CommonKernelSources.h"
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

void CommonUpdateStateDataKernel::initialize(const System& system) {
}

double CommonUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return cc.getTime();
}

void CommonUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    for (auto ctx : cc.getAllContexts())
        ctx->setTime(time);
}

long long CommonUpdateStateDataKernel::getStepCount(const ContextImpl& context) const {
    return cc.getStepCount();
}

void CommonUpdateStateDataKernel::setStepCount(const ContextImpl& context, long long count) {
    for (auto ctx : cc.getAllContexts())
        ctx->setStepCount(count);
}

void CommonUpdateStateDataKernel::getPositions(ContextImpl& context, vector<Vec3>& positions) {
    ContextSelector selector(cc);
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    vector<mm_float4> posCorrection;
    if (cc.getUseDoublePrecision()) {
        mm_double4* posq = (mm_double4*) cc.getPinnedBuffer();
        cc.getPosq().download(posq);
    }
    else if (cc.getUseMixedPrecision()) {
        mm_float4* posq = (mm_float4*) cc.getPinnedBuffer();
        cc.getPosq().download(posq, false);
        posCorrection.resize(numParticles);
        cc.getPosqCorrection().download(posCorrection);
    }
    else {
        mm_float4* posq = (mm_float4*) cc.getPinnedBuffer();
        cc.getPosq().download(posq);
    }
    
    // Filling in the output array is done in parallel for speed.
    
    cc.getThreadPool().execute([&] (ThreadPool& threads, int threadIndex) {
        // Compute the position of each particle to return to the user.  This is done in parallel for speed.
        
        const vector<int>& order = cc.getAtomIndex();
        int numParticles = cc.getNumAtoms();
        Vec3 boxVectors[3];
        cc.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        int numThreads = threads.getNumThreads();
        int start = threadIndex*numParticles/numThreads;
        int end = (threadIndex+1)*numParticles/numThreads;
        if (cc.getUseDoublePrecision()) {
            mm_double4* posq = (mm_double4*) cc.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                mm_double4 pos = posq[i];
                mm_int4 offset = cc.getPosCellOffsets()[i];
                positions[order[i]] = Vec3(pos.x, pos.y, pos.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
        else if (cc.getUseMixedPrecision()) {
            mm_float4* posq = (mm_float4*) cc.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                mm_float4 pos1 = posq[i];
                mm_float4 pos2 = posCorrection[i];
                mm_int4 offset = cc.getPosCellOffsets()[i];
                positions[order[i]] = Vec3((double)pos1.x+(double)pos2.x, (double)pos1.y+(double)pos2.y, (double)pos1.z+(double)pos2.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
        else {
            mm_float4* posq = (mm_float4*) cc.getPinnedBuffer();
            for (int i = start; i < end; ++i) {
                mm_float4 pos = posq[i];
                mm_int4 offset = cc.getPosCellOffsets()[i];
                positions[order[i]] = Vec3(pos.x, pos.y, pos.z)-boxVectors[0]*offset.x-boxVectors[1]*offset.y-boxVectors[2]*offset.z;
            }
        }
    });
    cc.getThreadPool().waitForThreads();
}

void CommonUpdateStateDataKernel::setPositions(ContextImpl& context, const vector<Vec3>& positions) {
    ContextSelector selector(cc);
    const vector<int>& order = cc.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cc.getUseDoublePrecision()) {
        mm_double4* posq = (mm_double4*) cc.getPinnedBuffer();
        cc.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4& pos = posq[i];
            const Vec3& p = positions[order[i]];
            pos.x = p[0];
            pos.y = p[1];
            pos.z = p[2];
        }
        for (int i = numParticles; i < cc.getPaddedNumAtoms(); i++)
            posq[i] = mm_double4(0.0, 0.0, 0.0, 0.0);
        cc.getPosq().upload(posq);
    }
    else {
        mm_float4* posq = (mm_float4*) cc.getPinnedBuffer();
        cc.getPosq().download(posq);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4& pos = posq[i];
            const Vec3& p = positions[order[i]];
            pos.x = (float) p[0];
            pos.y = (float) p[1];
            pos.z = (float) p[2];
        }
        for (int i = numParticles; i < cc.getPaddedNumAtoms(); i++)
            posq[i] = mm_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cc.getPosq().upload(posq);
    }
    if (cc.getUseMixedPrecision()) {
        mm_float4* posCorrection = (mm_float4*) cc.getPinnedBuffer();
        for (int i = 0; i < numParticles; ++i) {
            mm_float4& c = posCorrection[i];
            const Vec3& p = positions[order[i]];
            c.x = (float) (p[0]-(float)p[0]);
            c.y = (float) (p[1]-(float)p[1]);
            c.z = (float) (p[2]-(float)p[2]);
            c.w = 0;
        }
        for (int i = numParticles; i < cc.getPaddedNumAtoms(); i++)
            posCorrection[i] = mm_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cc.getPosqCorrection().upload(posCorrection);
    }
    for (auto& offset : cc.getPosCellOffsets())
        offset = mm_int4(0, 0, 0, 0);
    cc.reorderAtoms();
}

void CommonUpdateStateDataKernel::getVelocities(ContextImpl& context, vector<Vec3>& velocities) {
    ContextSelector selector(cc);
    const vector<int>& order = cc.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        mm_double4* velm = (mm_double4*) cc.getPinnedBuffer();
        cc.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4 vel = velm[i];
            velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
        }
    }
    else {
        mm_float4* velm = (mm_float4*) cc.getPinnedBuffer();
        cc.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4 vel = velm[i];
            velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
        }
    }
}

void CommonUpdateStateDataKernel::setVelocities(ContextImpl& context, const vector<Vec3>& velocities) {
    ContextSelector selector(cc);
    const vector<int>& order = cc.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        mm_double4* velm = (mm_double4*) cc.getPinnedBuffer();
        cc.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_double4& vel = velm[i];
            const Vec3& p = velocities[order[i]];
            vel.x = p[0];
            vel.y = p[1];
            vel.z = p[2];
        }
        for (int i = numParticles; i < cc.getPaddedNumAtoms(); i++)
            velm[i] = mm_double4(0.0, 0.0, 0.0, 0.0);
        cc.getVelm().upload(velm);
    }
    else {
        mm_float4* velm = (mm_float4*) cc.getPinnedBuffer();
        cc.getVelm().download(velm);
        for (int i = 0; i < numParticles; ++i) {
            mm_float4& vel = velm[i];
            const Vec3& p = velocities[order[i]];
            vel.x = p[0];
            vel.y = p[1];
            vel.z = p[2];
        }
        for (int i = numParticles; i < cc.getPaddedNumAtoms(); i++)
            velm[i] = mm_float4(0.0f, 0.0f, 0.0f, 0.0f);
        cc.getVelm().upload(velm);
    }
}

void CommonUpdateStateDataKernel::computeShiftedVelocities(ContextImpl& context, double timeShift, vector<Vec3>& velocities) {
    cc.getIntegrationUtilities().computeShiftedVelocities(timeShift, velocities);
}

void CommonUpdateStateDataKernel::getForces(ContextImpl& context, vector<Vec3>& forces) {
    ContextSelector selector(cc);
    long long* force = (long long*) cc.getPinnedBuffer();
    cc.getLongForceBuffer().download(force);
    const vector<int>& order = cc.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    int paddedNumParticles = cc.getPaddedNumAtoms();
    forces.resize(numParticles);
    double scale = 1.0/(double) 0x100000000LL;
    for (int i = 0; i < numParticles; ++i)
        forces[order[i]] = Vec3(scale*force[i], scale*force[i+paddedNumParticles], scale*force[i+paddedNumParticles*2]);
}

void CommonUpdateStateDataKernel::getEnergyParameterDerivatives(ContextImpl& context, map<string, double>& derivs) {
    ContextSelector selector(cc);
    const vector<string>& paramDerivNames = cc.getEnergyParamDerivNames();
    int numDerivs = paramDerivNames.size();
    if (numDerivs == 0)
        return;
    derivs = cc.getEnergyParamDerivWorkspace();
    ArrayInterface& derivArray = cc.getEnergyParamDerivBuffer();
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
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

void CommonUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    cc.getPeriodicBoxVectors(a, b, c);
}

void CommonUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) {
    if (!cc.getBoxIsTriclinic() && (b[0] != 0 || c[0] != 0 || c[1] != 0))
        throw OpenMMException("The box shape has changed from rectangular to triclinic.  To do this, you must call setDefaultPeriodicBoxVectors() on the System to specify a triclinic default box, then reinitialize the Context.");

    // If any particles have been wrapped to the first periodic box, we need to unwrap them
    // to avoid changing their positions.

    vector<Vec3> positions;
    for (auto offset : cc.getPosCellOffsets()) {
        if (offset.x != 0 || offset.y != 0 || offset.z != 0) {
            getPositions(context, positions);
            break;
        }
    }
    
    // Update the vectors.

    for (auto ctx : cc.getAllContexts())
        ctx->setPeriodicBoxVectors(a, b, c);
    if (positions.size() > 0)
        setPositions(context, positions);
}

void CommonUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    ContextSelector selector(cc);
    int version = 3;
    stream.write((char*) &version, sizeof(int));
    int precision = (cc.getUseDoublePrecision() ? 2 : cc.getUseMixedPrecision() ? 1 : 0);
    stream.write((char*) &precision, sizeof(int));
    double time = cc.getTime();
    stream.write((char*) &time, sizeof(double));
    long long stepCount = cc.getStepCount();
    stream.write((char*) &stepCount, sizeof(long long));
    int stepsSinceReorder = cc.getStepsSinceReorder();
    stream.write((char*) &stepsSinceReorder, sizeof(int));
    char* buffer = (char*) cc.getPinnedBuffer();
    cc.getPosq().download(buffer);
    stream.write(buffer, cc.getPosq().getSize()*cc.getPosq().getElementSize());
    if (cc.getUseMixedPrecision()) {
        cc.getPosqCorrection().download(buffer);
        stream.write(buffer, cc.getPosqCorrection().getSize()*cc.getPosqCorrection().getElementSize());
    }
    cc.getVelm().download(buffer);
    stream.write(buffer, cc.getVelm().getSize()*cc.getVelm().getElementSize());
    stream.write((char*) &cc.getAtomIndex()[0], sizeof(int)*cc.getAtomIndex().size());
    stream.write((char*) &cc.getPosCellOffsets()[0], sizeof(mm_int4)*cc.getPosCellOffsets().size());
    Vec3 boxVectors[3];
    cc.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    stream.write((char*) boxVectors, 3*sizeof(Vec3));
    cc.getIntegrationUtilities().createCheckpoint(stream);
    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void CommonUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    ContextSelector selector(cc);
    int version;
    stream.read((char*) &version, sizeof(int));
    if (version != 3)
        throw OpenMMException("Checkpoint was created with a different version of OpenMM");
    int precision;
    stream.read((char*) &precision, sizeof(int));
    int expectedPrecision = (cc.getUseDoublePrecision() ? 2 : cc.getUseMixedPrecision() ? 1 : 0);
    if (precision != expectedPrecision)
        throw OpenMMException("Checkpoint was created with a different numeric precision");
    double time;
    stream.read((char*) &time, sizeof(double));
    long long stepCount;
    stream.read((char*) &stepCount, sizeof(long long));
    int stepsSinceReorder;
    stream.read((char*) &stepsSinceReorder, sizeof(int));
    vector<ComputeContext*> contexts = cc.getAllContexts();
    for (auto ctx : contexts) {
        ctx->setTime(time);
        ctx->setStepCount(stepCount);
        ctx->setStepsSinceReorder(stepsSinceReorder);
    }
    char* buffer = (char*) cc.getPinnedBuffer();
    stream.read(buffer, cc.getPosq().getSize()*cc.getPosq().getElementSize());
    cc.getPosq().upload(buffer);
    if (cc.getUseMixedPrecision()) {
        stream.read(buffer, cc.getPosqCorrection().getSize()*cc.getPosqCorrection().getElementSize());
        cc.getPosqCorrection().upload(buffer);
    }
    stream.read(buffer, cc.getVelm().getSize()*cc.getVelm().getElementSize());
    cc.getVelm().upload(buffer);
    stream.read((char*) &cc.getAtomIndex()[0], sizeof(int)*cc.getAtomIndex().size());
    cc.getAtomIndexArray().upload(cc.getAtomIndex());
    stream.read((char*) &cc.getPosCellOffsets()[0], sizeof(mm_int4)*cc.getPosCellOffsets().size());
    Vec3 boxVectors[3];
    stream.read((char*) &boxVectors, 3*sizeof(Vec3));
    for (auto ctx : contexts)
        ctx->setPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    cc.getIntegrationUtilities().loadCheckpoint(stream);
    SimTKOpenMMUtilities::loadCheckpoint(stream);
    for (auto listener : cc.getReorderListeners())
        listener->execute();
    cc.validateAtomOrder();
}

void CommonApplyConstraintsKernel::initialize(const System& system) {
}

void CommonApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
    ContextSelector selector(cc);
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;
        map<string, string> defines;
        ComputeProgram program = cc.compileProgram(CommonKernelSources::constraints, defines);
        applyDeltasKernel = program->createKernel("applyPositionDeltas");
        applyDeltasKernel->addArg(cc.getNumAtoms());
        applyDeltasKernel->addArg(cc.getPosq());
        applyDeltasKernel->addArg(cc.getIntegrationUtilities().getPosDelta());
        if (cc.getUseMixedPrecision())
            applyDeltasKernel->addArg(cc.getPosqCorrection());
    }
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    cc.clearBuffer(integration.getPosDelta());
    integration.applyConstraints(tol);
    applyDeltasKernel->execute(cc.getNumAtoms());
    integration.computeVirtualSites();
}

void CommonApplyConstraintsKernel::applyToVelocities(ContextImpl& context, double tol) {
    cc.getIntegrationUtilities().applyVelocityConstraints(tol);
}

void CommonVirtualSitesKernel::initialize(const System& system) {
}

void CommonVirtualSitesKernel::computePositions(ContextImpl& context) {
    cc.getIntegrationUtilities().computeVirtualSites();
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
    ContextSelector selector(cc);
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

void CommonCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force, int firstBond, int lastBond) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0 || firstBond >= endIndex || lastBond < startIndex || firstBond > lastBond)
        return;
    firstBond = max(firstBond, startIndex);
    lastBond = min(lastBond, endIndex-1);
    
    // Record the per-bond parameters.
    
    int numToSet = lastBond-firstBond+1;
    vector<mm_float2> paramVector(numToSet);
    for (int i = 0; i < numToSet; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(firstBond+i, atom1, atom2, length, k);
        paramVector[i] = mm_float2((float) length, (float) k);
    }
    params.uploadSubArray(paramVector.data(), firstBond-startIndex, numToSet);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, false, true);
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
        thread_local static vector<double> parameters;
        force.getBondParameters(index, particle1, particle2, parameters);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        thread_local static vector<double> parameters1, parameters2;
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
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
}

void CommonCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = new ComputeParameterSet(cc, force.getNumPerBondParameters(), numBonds, "customBondParams");
    vector<vector<double> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], paramVector[i]);
    params->setParameterValues(paramVector, true);
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
    ContextSelector selector(cc);
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

void CommonCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0 || firstBond >= endIndex || lastBond < startIndex || firstBond > lastBond)
        return;
    firstBond = max(firstBond, startIndex);
    lastBond = min(lastBond, endIndex-1);
    
    // Record the per-bond parameters.
    
    int numToSet = lastBond-firstBond+1;
    vector<vector<double> > paramVector(numToSet);
    int atom1, atom2;
    for (int i = 0; i < numToSet; i++)
        force.getBondParameters(firstBond+i, atom1, atom2, paramVector[i]);
    params->setParameterValuesSubset(firstBond-startIndex, paramVector, true);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, false, true);
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
    ContextSelector selector(cc);
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

void CommonCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0 || firstAngle >= endIndex || lastAngle < startIndex || firstAngle > lastAngle)
        return;
    firstAngle = max(firstAngle, startIndex);
    lastAngle = min(lastAngle, endIndex-1);
    
    // Record the per-angle parameters.
    
    int numToSet = lastAngle-firstAngle+1;
    vector<mm_float2> paramVector(numToSet);
    for (int i = 0; i < numToSet; i++) {
        int atom1, atom2, atom3;
        double angle, k;
        force.getAngleParameters(firstAngle+i, atom1, atom2, atom3, angle, k);
        paramVector[i] = mm_float2((float) angle, (float) k);
    }
    params.uploadSubArray(paramVector.data(), firstAngle-startIndex, numToSet);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, false, true);
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
        thread_local static vector<double> parameters;
        force.getAngleParameters(index, particle1, particle2, particle3, parameters);
        particles.resize(3);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3;
        thread_local static vector<double> parameters1, parameters2;
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
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
}

void CommonCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = new ComputeParameterSet(cc, force.getNumPerAngleParameters(), numAngles, "customAngleParams");
    vector<vector<double> > paramVector(numAngles);
    for (int i = 0; i < numAngles; i++)
        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], paramVector[i]);
    params->setParameterValues(paramVector, true);
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
    ContextSelector selector(cc);
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

void CommonCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    if (numAngles == 0 || firstAngle >= endIndex || lastAngle < startIndex || firstAngle > lastAngle)
        return;
    firstAngle = max(firstAngle, startIndex);
    lastAngle = min(lastAngle, endIndex-1);
    
    // Record the per-angle parameters.
    
    int numToSet = lastAngle-firstAngle+1;
    vector<vector<double> > paramVector(numToSet);
    int atom1, atom2, atom3;
    for (int i = 0; i < numToSet; i++)
        force.getAngleParameters(firstAngle+i, atom1, atom2, atom3, paramVector[i]);
    params->setParameterValuesSubset(firstAngle-startIndex, paramVector, true);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, false, true);
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
    ContextSelector selector(cc);
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

void CommonCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0 || firstTorsion >= endIndex || lastTorsion < startIndex || firstTorsion > lastTorsion)
        return;
    firstTorsion = max(firstTorsion, startIndex);
    lastTorsion = min(lastTorsion, endIndex-1);
    
    // Record the per-torsion parameters.
    
    int numToSet = lastTorsion-firstTorsion+1;
    vector<mm_float4> paramVector(numToSet);
    for (int i = 0; i < numToSet; i++) {
        int atom1, atom2, atom3, atom4, periodicity;
        double phase, k;
        force.getTorsionParameters(firstTorsion+i, atom1, atom2, atom3, atom4, periodicity, phase, k);
        paramVector[i] = mm_float4((float) k, (float) phase, (float) periodicity, 0.0f);
    }
    params.uploadSubArray(paramVector.data(), firstTorsion-startIndex, numToSet);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, false, true);
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
    ContextSelector selector(cc);
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
    ContextSelector selector(cc);
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
    
    cc.invalidateMolecules(info, false, true);
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
        thread_local static vector<double> parameters;
        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, parameters);
        particles.resize(4);
        particles[0] = particle1;
        particles[1] = particle2;
        particles[2] = particle3;
        particles[3] = particle4;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2, particle3, particle4;
        thread_local static vector<double> parameters1, parameters2;
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
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params = new ComputeParameterSet(cc, force.getNumPerTorsionParameters(), numTorsions, "customTorsionParams");
    vector<vector<double> > paramVector(numTorsions);
    for (int i = 0; i < numTorsions; i++)
        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], paramVector[i]);
    params->setParameterValues(paramVector, true);
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
    ContextSelector selector(cc);
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

void CommonCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force, int firstTorsion, int lastTorsion) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    if (numTorsions == 0 || firstTorsion >= endIndex || lastTorsion < startIndex || firstTorsion > lastTorsion)
        return;
    firstTorsion = max(firstTorsion, startIndex);
    lastTorsion = min(lastTorsion, endIndex-1);

    // Record the per-torsion parameters.

    int numToSet = lastTorsion-firstTorsion+1;
    vector<vector<double> > paramVector(numToSet);
    int atom1, atom2, atom3, atom4;
    for (int i = 0; i < numToSet; i++)
        force.getTorsionParameters(firstTorsion+i, atom1, atom2, atom3, atom4, paramVector[i]);
    params->setParameterValuesSubset(firstTorsion-startIndex, paramVector, true);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, false, true);
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
    ContextSelector selector(cc);
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

    ContextSelector selector(cc);
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
        thread_local static vector<double> params1, params2;
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
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
}

void CommonCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumParticles()/numContexts;
    numParticles = endIndex-startIndex;
    if (numParticles == 0)
        return;
    vector<vector<int> > atoms(numParticles, vector<int>(1));
    params = new ComputeParameterSet(cc, force.getNumPerParticleParameters(), numParticles, "customExternalParams");
    vector<vector<double> > paramVector(numParticles);
    for (int i = 0; i < numParticles; i++)
        force.getParticleParameters(startIndex+i, atoms[i][0], paramVector[i]);
    params->setParameterValues(paramVector, true);
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
    ContextSelector selector(cc);
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

void CommonCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force, int firstParticle, int lastParticle) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumParticles()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumParticles()/numContexts;
    if (numParticles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    if (numParticles == 0 || firstParticle >= endIndex || lastParticle < startIndex || firstParticle > lastParticle)
        return;
    firstParticle = max(firstParticle, startIndex);
    lastParticle = min(lastParticle, endIndex-1);
    
    // Record the per-particle parameters.
    
    int numToSet = lastParticle-firstParticle+1;
    vector<vector<double> > paramVector(numToSet);
    int particle;
    for (int i = 0; i < numToSet; i++)
        force.getParticleParameters(firstParticle+i, particle, paramVector[i]);
    params->setParameterValuesSubset(firstParticle-startIndex, paramVector, true);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, true, false);
}

class CommonCalcCustomCompoundBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomCompoundBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        thread_local static vector<double> parameters;
        force.getBondParameters(index, particles, parameters);
    }
    bool areGroupsIdentical(int group1, int group2) {
        thread_local static vector<int> particles;
        thread_local static vector<double> parameters1, parameters2;
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
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
}

void CommonCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    int particlesPerBond = force.getNumParticlesPerBond();
    vector<vector<int> > atoms(numBonds, vector<int>(particlesPerBond));
    params = new ComputeParameterSet(cc, force.getNumPerBondParameters(), numBonds, "customCompoundBondParams", false, cc.getUseDoublePrecision());
    vector<vector<double> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(startIndex+i, atoms[i], paramVector[i]);
    params->setParameterValues(paramVector, true);
    info = new ForceInfo(force);
    cc.addForce(info);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    tabulatedFunctionArrays.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctionArrays[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctionArrays[i].upload(f);
        string arrayName = cc.getBondedUtilities().addArgument(tabulatedFunctionArrays[i], width == 1 ? "float" : "float"+cc.intToString(width));
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

    // Generate the kernel.

    Lepton::ParsedExpression energyExpression = CustomCompoundBondForceImpl::prepareExpression(force, functions);
    map<string, Lepton::ParsedExpression> forceExpressions;
    stringstream compute;
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
    vector<string> forceNames;
    for (int i = 0; i < particlesPerBond; i++) {
        string istr = cc.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = make_real3(0);\n";
        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
        if (!isZeroExpression(forceExpressionX))
            forceExpressions[forceName+".x -= "] = forceExpressionX;
        if (!isZeroExpression(forceExpressionY))
            forceExpressions[forceName+".y -= "] = forceExpressionY;
        if (!isZeroExpression(forceExpressionZ))
            forceExpressions[forceName+".z -= "] = forceExpressionZ;
    }
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp", "real", force.usesPeriodicBoundaryConditions());
    cc.getBondedUtilities().addInteraction(atoms, compute.str(), force.getForceGroup());
    map<string, string> replacements;
    replacements["M_PI"] = cc.doubleToString(M_PI);
    cc.getBondedUtilities().addPrefixCode(cc.replaceStrings(CommonKernelSources::pointFunctions, replacements));
}

double CommonCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
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
    ContextSelector selector(cc);
    int numContexts = cc.getNumContexts();
    int startIndex = cc.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cc.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;

    // Record the per-bond parameters.

    vector<vector<double> > paramVector(numBonds);
    vector<int> particles;
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(startIndex+i, particles, paramVector[i]);
    params->setParameterValues(paramVector, true);

    // See if any tabulated functions have changed.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            int width;
            vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
            tabulatedFunctionArrays[i].upload(f);
        }
    }

    // Mark that the current reordering may be invalid.

    cc.invalidateMolecules(info, false, true);
}

class CommonCalcCustomCentroidBondForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(const CustomCentroidBondForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        thread_local static vector<double> parameters;
        thread_local static vector<int> groups;
        force.getBondParameters(index, groups, parameters);
        for (int group : groups) {
            vector<int> groupParticles;
            vector<double> weights;
            force.getGroupParameters(group, groupParticles, weights);
            particles.insert(particles.end(), groupParticles.begin(), groupParticles.end());
        }
    }
    bool areGroupsIdentical(int group1, int group2) {
        thread_local static vector<int> groups1, groups2;
        thread_local static vector<double> parameters1, parameters2;
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
    ContextSelector selector(cc);
    if (params != NULL)
        delete params;
}

void CommonCalcCustomCentroidBondForceKernel::initialize(const System& system, const CustomCentroidBondForce& force) {
    ContextSelector selector(cc);
    numBonds = force.getNumBonds();
    if (numBonds == 0)
        return;
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
    params = new ComputeParameterSet(cc, force.getNumPerBondParameters(), numBonds, "customCentroidBondParams", false, cc.getUseDoublePrecision());
    vector<vector<double> > paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        vector<int> groups;
        force.getBondParameters(i, groups, paramVector[i]);
        for (int j = 0; j < groups.size(); j++)
            bondGroupVec[i+j*numBonds] = groups[j];
    }
    params->setParameterValues(paramVector, true);
    bondGroups.initialize<int>(cc, bondGroupVec.size(), "bondGroups");
    bondGroups.upload(bondGroupVec);

    // Record the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string, string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream extraArgs;
    tabulatedFunctionArrays.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
        string arrayName = "table"+cc.intToString(i);
        functionDefinitions.push_back(make_pair(name, arrayName));
        functions[name] = cc.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
        tabulatedFunctionArrays[i].initialize<float>(cc, f.size(), "TabulatedFunction");
        tabulatedFunctionArrays[i].upload(f);
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

    // Generate the kernel.

    Lepton::ParsedExpression energyExpression = CustomCentroidBondForceImpl::prepareExpression(force, functions);
    map<string, Lepton::ParsedExpression> forceExpressions;
    stringstream compute, initParamDerivs, saveParamDerivs;
    for (int i = 0; i < groupsPerBond; i++) {
        compute<<"int group"<<(i+1)<<" = bondGroups[index+"<<(i*numBonds)<<"];\n";
        compute<<"real4 pos"<<(i+1)<<" = centerPositions[group"<<(i+1)<<"];\n";
    }
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
    vector<string> forceNames;
    for (int i = 0; i < groupsPerBond; i++) {
        string istr = cc.intToString(i+1);
        string forceName = "force"+istr;
        forceNames.push_back(forceName);
        compute<<"real3 "<<forceName<<" = make_real3(0);\n";
        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
        if (!isZeroExpression(forceExpressionX))
            forceExpressions[forceName+".x -= "] = forceExpressionX;
        if (!isZeroExpression(forceExpressionY))
            forceExpressions[forceName+".y -= "] = forceExpressionY;
        if (!isZeroExpression(forceExpressionZ))
            forceExpressions[forceName+".z -= "] = forceExpressionZ;
    }
    compute << cc.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions, "temp", "real", force.usesPeriodicBoundaryConditions());
    
    // Save the forces to global memory.
    
    for (int i = 0; i < groupsPerBond; i++) {
        compute<<"ATOMIC_ADD(&groupForce[group"<<(i+1)<<"], (mm_ulong) realToFixedPoint(force"<<(i+1)<<".x));\n";
        compute<<"ATOMIC_ADD(&groupForce[group"<<(i+1)<<"+numParticleGroups], (mm_ulong) realToFixedPoint(force"<<(i+1)<<".y));\n";
        compute<<"ATOMIC_ADD(&groupForce[group"<<(i+1)<<"+numParticleGroups*2], (mm_ulong) realToFixedPoint(force"<<(i+1)<<".z));\n";
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
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::pointFunctions+CommonKernelSources::customCentroidBond, replacements));
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
    for (auto& function : tabulatedFunctionArrays)
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
    ContextSelector selector(cc);
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
    ContextSelector selector(cc);
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;

    // Record the per-bond parameters.

    vector<vector<double> > paramVector(numBonds);
    vector<int> particles;
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, particles, paramVector[i]);
    params->setParameterValues(paramVector, true);

    // See if any tabulated functions have changed.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            int width;
            vector<float> f = cc.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i), width);
            tabulatedFunctionArrays[i].upload(f);
        }
    }

    // Mark that the current reordering may be invalid.

    cc.invalidateMolecules(info, false, true);
}

class CommonCalcGBSAOBCForceKernel::ForceInfo : public ComputeForceInfo {
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

void CommonCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    ContextSelector selector(cc);
    if (cc.getNumContexts() > 1)
        throw OpenMMException("GBSAOBCForce does not support using multiple devices");
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
        ;
    string prefix = "obc"+cc.intToString(forceIndex)+"_";
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    params.initialize<mm_float2>(cc, cc.getPaddedNumAtoms(), "gbsaObcParams");
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    charges.initialize(cc, cc.getPaddedNumAtoms(), elementSize, "gbsaObcCharges");
    bornRadii.initialize(cc, cc.getPaddedNumAtoms(), elementSize, "bornRadii");
    obcChain.initialize(cc, cc.getPaddedNumAtoms(), elementSize, "obcChain");
    bornSum.initialize<long long>(cc, cc.getPaddedNumAtoms(), "bornSum");
    bornForce.initialize<long long>(cc, cc.getPaddedNumAtoms(), "bornForce");
    cc.addAutoclearBuffer(bornSum);
    cc.addAutoclearBuffer(bornForce);
    vector<double> chargeVec(cc.getPaddedNumAtoms());
    vector<mm_float2> paramsVector(cc.getPaddedNumAtoms(), mm_float2(1,1));
    const double dielectricOffset = 0.009;
    for (int i = 0; i < force.getNumParticles(); i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        radius -= dielectricOffset;
        chargeVec[i] = charge;
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
    }
    charges.upload(chargeVec, true);
    params.upload(paramsVector);
    prefactor = -ONE_4PI_EPS0*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
    surfaceAreaFactor = -6.0*4*M_PI*force.getSurfaceAreaEnergy();
    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
    cutoff = force.getCutoffDistance();
    string source = CommonKernelSources::gbsaObc2;
    map<string, string> replacements;
    replacements["CHARGE1"] = prefix+"charge1";
    replacements["CHARGE2"] = prefix+"charge2";
    replacements["OBC_PARAMS1"] = prefix+"obcParams1";
    replacements["OBC_PARAMS2"] = prefix+"obcParams2";
    replacements["BORN_FORCE1"] = prefix+"bornForce1";
    replacements["BORN_FORCE2"] = prefix+"bornForce2";
    source = cc.replaceStrings(source, replacements);
    nb.addInteraction(useCutoff, usePeriodic, false, cutoff, vector<vector<int> >(), source, force.getForceGroup());
    nb.addParameter(ComputeParameterInfo(charges, prefix+"charge", "float", 1));
    nb.addParameter(ComputeParameterInfo(params, prefix+"obcParams", "float", 2));
    nb.addParameter(ComputeParameterInfo(bornForce, prefix+"bornForce", "mm_long", 1));
    info = new ForceInfo(force);
    cc.addForce(info);
}

double CommonCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    bool deviceIsCpu = cc.getIsCPU();
    if (!hasCreatedKernels) {
        // These Kernels cannot be created in initialize(), because the NonbondedUtilities has not been initialized yet then.

        hasCreatedKernels = true;
        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : 0);
        int numAtomBlocks = cc.getPaddedNumAtoms()/32;
        map<string, string> defines;
        if (nb.getUseCutoff())
            defines["USE_CUTOFF"] = "1";
        if (nb.getUsePeriodic())
            defines["USE_PERIODIC"] = "1";
        defines["CUTOFF_SQUARED"] = cc.doubleToString(cutoff*cutoff);
        defines["CUTOFF"] = cc.doubleToString(cutoff);
        defines["PREFACTOR"] = cc.doubleToString(prefactor);
        defines["SURFACE_AREA_FACTOR"] = cc.doubleToString(surfaceAreaFactor);
        defines["NUM_ATOMS"] = cc.intToString(cc.getNumAtoms());
        defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
        defines["NUM_BLOCKS"] = cc.intToString(numAtomBlocks);
        defines["FORCE_WORK_GROUP_SIZE"] = cc.intToString(nb.getForceThreadBlockSize());
        defines["TILE_SIZE"] = "32";
        int numExclusionTiles = nb.getExclusionTiles().getSize();
        defines["NUM_TILES_WITH_EXCLUSIONS"] = cc.intToString(numExclusionTiles);
        defines["FIRST_EXCLUSION_TILE"] = "0";
        defines["LAST_EXCLUSION_TILE"] = cc.intToString(numExclusionTiles);
        string file;
        if (deviceIsCpu)
            file = CommonKernelSources::gbsaObc_cpu;
        else
            file = CommonKernelSources::gbsaObc;
        ComputeProgram program = cc.compileProgram(file, defines);
        computeBornSumKernel = program->createKernel("computeBornSum");
        computeBornSumKernel->addArg(bornSum);
        computeBornSumKernel->addArg(cc.getPosq());
        computeBornSumKernel->addArg(charges);
        computeBornSumKernel->addArg(params);
        if (nb.getUseCutoff()) {
            computeBornSumKernel->addArg(nb.getInteractingTiles());
            computeBornSumKernel->addArg(nb.getInteractionCount());
            for (int i = 0; i < 5; i++)
                computeBornSumKernel->addArg(); // The periodic box size arguments are set when the kernel is executed.
            computeBornSumKernel->addArg(maxTiles);
            computeBornSumKernel->addArg(nb.getBlockCenters());
            computeBornSumKernel->addArg(nb.getBlockBoundingBoxes());
            computeBornSumKernel->addArg(nb.getInteractingAtoms());
        }
        else
            computeBornSumKernel->addArg(numAtomBlocks*(numAtomBlocks+1)/2);
        computeBornSumKernel->addArg(nb.getExclusionTiles());
        force1Kernel = program->createKernel("computeGBSAForce1");
        force1Kernel->addArg(cc.getLongForceBuffer());
        force1Kernel->addArg(bornForce);
        force1Kernel->addArg(cc.getEnergyBuffer());
        force1Kernel->addArg(cc.getPosq());
        force1Kernel->addArg(charges);
        force1Kernel->addArg(bornRadii);
        force1Kernel->addArg(); // Whether to include energy.
        if (nb.getUseCutoff()) {
            force1Kernel->addArg(nb.getInteractingTiles());
            force1Kernel->addArg(nb.getInteractionCount());
            for (int i = 0; i < 5; i++)
                force1Kernel->addArg(); // The periodic box size arguments are set when the kernel is executed.
            force1Kernel->addArg(maxTiles);
            force1Kernel->addArg(nb.getBlockCenters());
            force1Kernel->addArg(nb.getBlockBoundingBoxes());
            force1Kernel->addArg(nb.getInteractingAtoms());
        }
        else
            force1Kernel->addArg(numAtomBlocks*(numAtomBlocks+1)/2);
        force1Kernel->addArg(nb.getExclusionTiles());
        program = cc.compileProgram(CommonKernelSources::gbsaObcReductions, defines);
        reduceBornSumKernel = program->createKernel("reduceBornSum");
        reduceBornSumKernel->addArg(1.0f);
        reduceBornSumKernel->addArg(0.8f);
        reduceBornSumKernel->addArg(4.85f);
        reduceBornSumKernel->addArg(bornSum);
        reduceBornSumKernel->addArg(params);
        reduceBornSumKernel->addArg(bornRadii);
        reduceBornSumKernel->addArg(obcChain);
        reduceBornForceKernel = program->createKernel("reduceBornForce");
        reduceBornForceKernel->addArg(bornForce);
        reduceBornForceKernel->addArg(cc.getEnergyBuffer());
        reduceBornForceKernel->addArg(params);
        reduceBornForceKernel->addArg(bornRadii);
        reduceBornForceKernel->addArg(obcChain);
    }
    force1Kernel->setArg(6, (int) includeEnergy);
    if (nb.getUseCutoff()) {
        setPeriodicBoxArgs(cc, computeBornSumKernel, 6);
        setPeriodicBoxArgs(cc, force1Kernel, 9);
        if (maxTiles < nb.getInteractingTiles().getSize()) {
            maxTiles = nb.getInteractingTiles().getSize();
            computeBornSumKernel->setArg(11, maxTiles);
            force1Kernel->setArg(14, maxTiles);
        }
    }
    computeBornSumKernel->execute(nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    reduceBornSumKernel->execute(cc.getPaddedNumAtoms());
    force1Kernel->execute(nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
    reduceBornForceKernel->execute(cc.getPaddedNumAtoms());
    return 0.0;
}

void CommonCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
    // Make sure the new parameters are acceptable.
    
    ContextSelector selector(cc);
    int numParticles = force.getNumParticles();
    if (numParticles != cc.getNumAtoms())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    
    // Record the per-particle parameters.
    
    vector<double> chargeVector(cc.getPaddedNumAtoms(), 0.0);
    vector<mm_float2> paramsVector(cc.getPaddedNumAtoms());
    const double dielectricOffset = 0.009;
    for (int i = 0; i < numParticles; i++) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        chargeVector[i] = charge;
        radius -= dielectricOffset;
        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
    }
    for (int i = numParticles; i < cc.getPaddedNumAtoms(); i++)
        paramsVector[i] = mm_float2(1,1);
    charges.upload(chargeVector, true);
    params.upload(paramsVector);
    
    // Mark that the current reordering may be invalid.
    
    cc.invalidateMolecules(info, true, false);
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
    // Initialize interactions.

    ContextSelector selector(cc);
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
    ContextSelector selector(cc);
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
    
    ContextSelector selector(cc);
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

class CommonCalcCustomCVForceKernel::ForceInfo : public ComputeForceInfo {
public:
    ForceInfo(ComputeForceInfo& force) : force(force) {
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

class CommonCalcCustomCVForceKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(ComputeContext& cc, ArrayInterface& invAtomOrder) : cc(cc), invAtomOrder(invAtomOrder) {
    }
    void execute() {
        vector<int> invOrder(cc.getPaddedNumAtoms());
        const vector<int>& order = cc.getAtomIndex();
        for (int i = 0; i < order.size(); i++)
            invOrder[order[i]] = i;
        invAtomOrder.upload(invOrder);
    }
private:
    ComputeContext& cc;
    ArrayInterface& invAtomOrder;
};

// This class allows us to update tabulated functions without having to recompile expressions
// that use them.
class CommonCalcCustomCVForceKernel::TabulatedFunctionWrapper : public CustomFunction {
public:
    TabulatedFunctionWrapper(vector<Lepton::CustomFunction*>& tabulatedFunctions, int index) :
            tabulatedFunctions(tabulatedFunctions), index(index) {
    }
    int getNumArguments() const {
        return tabulatedFunctions[index]->getNumArguments();
    }
    double evaluate(const double* arguments) const {
        return tabulatedFunctions[index]->evaluate(arguments);
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        return tabulatedFunctions[index]->evaluateDerivative(arguments, derivOrder);
    }
    CustomFunction* clone() const {
        return new TabulatedFunctionWrapper(tabulatedFunctions, index);
    }
private:
    vector<Lepton::CustomFunction*>& tabulatedFunctions;    
    int index;
};

void CommonCalcCustomCVForceKernel::initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext) {
    ContextSelector selector(cc);
    int numCVs = force.getNumCollectiveVariables();
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < numCVs; i++)
        variableNames.push_back(force.getCollectiveVariableName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string name = force.getEnergyParameterDerivativeName(i);
        paramDerivNames.push_back(name);
        cc.addEnergyParameterDerivative(name);
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions(), NULL);
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        tabulatedFunctions[i] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
        functions[force.getTabulatedFunctionName(i)] = new TabulatedFunctionWrapper(tabulatedFunctions, i);
    }

    // Create the expressions.

    Lepton::ParsedExpression energyExpr = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    energyExpression = energyExpr.createCompiledExpression();
    variableDerivExpressions.clear();
    for (auto& name : variableNames)
        variableDerivExpressions.push_back(energyExpr.differentiate(name).createCompiledExpression());
    paramDerivExpressions.clear();
    for (auto& name : paramDerivNames)
        paramDerivExpressions.push_back(energyExpr.differentiate(name).createCompiledExpression());
    globalValues.resize(globalParameterNames.size());
    cvValues.resize(numCVs);
    map<string, double*> variableLocations;
    for (int i = 0; i < globalParameterNames.size(); i++)
        variableLocations[globalParameterNames[i]] = &globalValues[i];
    for (int i = 0; i < numCVs; i++)
        variableLocations[variableNames[i]] = &cvValues[i];
    energyExpression.setVariableLocations(variableLocations);
    for (CompiledExpression& expr : variableDerivExpressions)
        expr.setVariableLocations(variableLocations);
    for (CompiledExpression& expr : paramDerivExpressions)
        expr.setVariableLocations(variableLocations);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;

    // Copy parameter derivatives from the inner context.

    ComputeContext& cc2 = getInnerComputeContext(innerContext);
    for (auto& param : cc2.getEnergyParamDerivNames())
        cc.addEnergyParameterDerivative(param);
    
    // Create arrays for storing information.
    
    cvForces.resize(numCVs);
    for (int i = 0; i < numCVs; i++)
        cvForces[i].initialize<long long>(cc, 3*cc.getPaddedNumAtoms(), "cvForce");
    invAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "invAtomOrder");
    innerInvAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "innerInvAtomOrder");
    
    // Create the kernels.
    
    stringstream args, add;
    for (int i = 0; i < numCVs; i++) {
        args << ", GLOBAL mm_long * RESTRICT force" << i << ", real dEdV" << i;
        add << "forces[i] += (mm_long) (force" << i << "[i]*dEdV" << i << ");\n";
    }
    map<string, string> replacements;
    replacements["PARAMETER_ARGUMENTS"] = args.str();
    replacements["ADD_FORCES"] = add.str();
    ComputeProgram program = cc.compileProgram(cc.replaceStrings(CommonKernelSources::customCVForce, replacements));
    copyStateKernel = program->createKernel("copyState");
    copyStateKernel->addArg(cc.getPosq());
    copyStateKernel->addArg(cc2.getPosq());
    if (cc.getUseMixedPrecision()) {
        copyStateKernel->addArg(cc.getPosqCorrection());
        copyStateKernel->addArg(cc2.getPosqCorrection());
    }
    copyStateKernel->addArg(cc.getVelm());
    copyStateKernel->addArg(cc2.getVelm());
    copyStateKernel->addArg(cc.getAtomIndexArray());
    copyStateKernel->addArg(innerInvAtomOrder);
    copyStateKernel->addArg(cc.getNumAtoms());
    copyForcesKernel = program->createKernel("copyForces");
    copyForcesKernel->addArg();
    copyForcesKernel->addArg(invAtomOrder);
    copyForcesKernel->addArg(cc2.getLongForceBuffer());
    copyForcesKernel->addArg(cc2.getAtomIndexArray());
    copyForcesKernel->addArg(cc.getNumAtoms());
    copyForcesKernel->addArg(cc.getPaddedNumAtoms());
    addForcesKernel = program->createKernel("addForces");
    addForcesKernel->addArg(cc.getLongForceBuffer());
    addForcesKernel->addArg((int) cc.getLongForceBuffer().getSize());
    for (int i = 0; i < numCVs; i++) {
        addForcesKernel->addArg();
        addForcesKernel->addArg();
    }

    // This context needs to respect all forces in the inner context when reordering atoms.

    for (auto* info : cc2.getForceInfos())
        cc.addForce(new ForceInfo(*info));
}

CommonCalcCustomCVForceKernel::~CommonCalcCustomCVForceKernel() {
    for (int i = 0; i < tabulatedFunctions.size(); i++)
        if (tabulatedFunctions[i] != NULL)
            delete tabulatedFunctions[i];
}

double CommonCalcCustomCVForceKernel::execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy) {
    copyState(context, innerContext);
    int numCVs = variableNames.size();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    vector<map<string, double> > cvDerivs(numCVs);
    for (int i = 0; i < numCVs; i++) {
        cvValues[i] = innerContext.calcForcesAndEnergy(true, true, 1<<i);
        ContextSelector selector(cc);
        copyForcesKernel->setArg(0, cvForces[i]);
        copyForcesKernel->execute(numAtoms);
        innerContext.getEnergyParameterDerivatives(cvDerivs[i]);
    }

    // Compute the energy and forces.

    ContextSelector selector(cc);
    for (int i = 0; i < globalParameterNames.size(); i++)
        globalValues[i] = context.getParameter(globalParameterNames[i]);
    double energy = energyExpression.evaluate();
    for (int i = 0; i < numCVs; i++) {
        double dEdV = variableDerivExpressions[i].evaluate();
        addForcesKernel->setArg(2*i+2, cvForces[i]);
        if (cc.getUseDoublePrecision())
            addForcesKernel->setArg(2*i+3, dEdV);
        else
            addForcesKernel->setArg(2*i+3, (float) dEdV);
    }
    addForcesKernel->execute(numAtoms);

    // Compute the energy parameter derivatives.

    if (paramDerivExpressions.size() > 0) {
        map<string, double>& energyParamDerivs = cc.getEnergyParamDerivWorkspace();
        for (int i = 0; i < paramDerivExpressions.size(); i++)
            energyParamDerivs[paramDerivNames[i]] += paramDerivExpressions[i].evaluate();
        for (int i = 0; i < numCVs; i++) {
            double dEdV = variableDerivExpressions[i].evaluate();
            for (auto& deriv : cvDerivs[i])
                energyParamDerivs[deriv.first] += dEdV*deriv.second;
        }
    }
    return energy;
}

void CommonCalcCustomCVForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext) {
    ContextSelector selector(cc);
    int numAtoms = cc.getNumAtoms();
    ComputeContext& cc2 = getInnerComputeContext(innerContext);
    if (!hasInitializedListeners) {
        hasInitializedListeners = true;
        
        // Initialize the listeners.
        
        ReorderListener* listener1 = new ReorderListener(cc, invAtomOrder);
        ReorderListener* listener2 = new ReorderListener(cc2, innerInvAtomOrder);
        cc.addReorderListener(listener1);
        cc2.addReorderListener(listener2);
        listener1->execute();
        listener2->execute();
    }
    cc2.reorderAtoms();
    copyStateKernel->execute(numAtoms);
    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext.setPeriodicBoxVectors(a, b, c);
    innerContext.setTime(context.getTime());
    map<string, double> innerParameters = innerContext.getParameters();
    for (auto& param : innerParameters)
        innerContext.setParameter(param.first, context.getParameter(param.first));
}

void CommonCalcCustomCVForceKernel::copyParametersToContext(ContextImpl& context, const CustomCVForce& force) {
    // Create custom functions for the tabulated functions.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        if (tabulatedFunctions[i] != NULL) {
            delete tabulatedFunctions[i];
            tabulatedFunctions[i] = NULL;
        }
        tabulatedFunctions[i] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));
    }
}

void CommonIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::verlet);
    kernel1 = program->createKernel("integrateVerletPart1");
    kernel2 = program->createKernel("integrateVerletPart2");
}

void CommonIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    ContextSelector selector(cc);
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

    flushPeriodically(cc);
}

double CommonIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

void CommonIntegrateLangevinMiddleStepKernel::initialize(const System& system, const LangevinMiddleIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    cc.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    ComputeProgram program = cc.compileProgram(CommonKernelSources::langevinMiddle);
    kernel1 = program->createKernel("integrateLangevinMiddlePart1");
    kernel2 = program->createKernel("integrateLangevinMiddlePart2");
    kernel3 = program->createKernel("integrateLangevinMiddlePart3");
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        params.initialize<double>(cc, 2, "langevinMiddleParams");
        oldDelta.initialize<mm_double4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    }
    else {
        params.initialize<float>(cc, 2, "langevinMiddleParams");
        oldDelta.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    }
    prevStepSize = -1.0;
}

void CommonIntegrateLangevinMiddleStepKernel::execute(ContextImpl& context, const LangevinMiddleIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getStepSize());
        kernel2->addArg(numAtoms);
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        kernel2->addArg(oldDelta);
        kernel2->addArg(params);
        kernel2->addArg(integration.getStepSize());
        kernel2->addArg(integration.getRandom());
        kernel2->addArg(); // Random index will be set just before it is executed.
        kernel3->addArg(numAtoms);
        kernel3->addArg(cc.getPosq());
        kernel3->addArg(cc.getVelm());
        kernel3->addArg(integration.getPosDelta());
        kernel3->addArg(oldDelta);
        kernel3->addArg(integration.getStepSize());
        if (cc.getUseMixedPrecision())
            kernel3->addArg(cc.getPosqCorrection());
    }
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    cc.getIntegrationUtilities().setNextStepSize(stepSize);
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Calculate the integration parameters.

        double kT = BOLTZ*temperature;
        double vscale = exp(-stepSize*friction);
        double noisescale = sqrt(kT*(1-vscale*vscale));
        vector<double> p(params.getSize());
        p[0] = vscale;
        p[1] = noisescale;
        params.upload(p, true);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Perform the integration.

    kernel2->setArg(7, integration.prepareRandomNumbers(cc.getPaddedNumAtoms()));
    kernel1->execute(numAtoms);
    integration.applyVelocityConstraints(integrator.getConstraintTolerance());
    kernel2->execute(numAtoms);
    integration.applyConstraints(integrator.getConstraintTolerance());
    kernel3->execute(numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cc.setTime(cc.getTime()+stepSize);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    
    // Reduce UI lag.

    flushPeriodically(cc);
}

double CommonIntegrateLangevinMiddleStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinMiddleIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.0);
}

void CommonIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    cc.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    ComputeProgram program = cc.compileProgram(CommonKernelSources::brownian);
    kernel1 = program->createKernel("integrateBrownianPart1");
    kernel2 = program->createKernel("integrateBrownianPart2");
    prevStepSize = -1.0;
}

void CommonIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(); // tauDeltaT
        kernel1->addArg(); // noiseAmplitude
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getPosDelta());
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(integration.getRandom());
        kernel1->addArg(); // Random index will be set just before it is executed.
        kernel2->addArg(numAtoms);
        kernel2->addArg(); // oneOverDeltaT
        kernel2->addArg(cc.getPosq());
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        if (cc.getUseMixedPrecision())
            kernel2->addArg(cc.getPosqCorrection());
    }
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
            kernel1->setArg(2, tau*stepSize);
            kernel1->setArg(3, sqrt(2.0f*BOLTZ*temperature*stepSize*tau));
            kernel2->setArg(1, 1.0/stepSize);
        }
        else {
            kernel1->setArg(2, (float) (tau*stepSize));
            kernel1->setArg(3, (float) (sqrt(2.0f*BOLTZ*temperature*stepSize*tau)));
            kernel2->setArg(1, (float) (1.0/stepSize));
        }
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }

    // Call the first integration kernel.

    kernel1->setArg(8, integration.prepareRandomNumbers(cc.getPaddedNumAtoms()));
    kernel1->execute(numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    kernel2->execute(numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cc.setTime(cc.getTime()+stepSize);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();
    
    // Reduce UI lag.

    flushPeriodically(cc);
}

double CommonIntegrateBrownianStepKernel::computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0);
}

void CommonIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::verlet);
    kernel1 = program->createKernel("integrateVerletPart1");
    kernel2 = program->createKernel("integrateVerletPart2");
    selectSizeKernel = program->createKernel("selectVerletStepSize");
    blockSize = min(256, system.getNumParticles());
}

double CommonIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    ContextSelector selector(cc);
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

    flushPeriodically(cc);

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

void CommonIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    cc.initializeContexts();
    ContextSelector selector(cc);
    cc.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
    ComputeProgram program = cc.compileProgram(CommonKernelSources::langevinMiddle);
    kernel1 = program->createKernel("integrateLangevinMiddlePart1");
    kernel2 = program->createKernel("integrateLangevinMiddlePart2");
    kernel3 = program->createKernel("integrateLangevinMiddlePart3");
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision()) {
        params.initialize<double>(cc, 3, "langevinMiddleParams");
        oldDelta.initialize<mm_double4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    }
    else {
        params.initialize<float>(cc, 3, "langevinMiddleParams");
        oldDelta.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    }
    selectSizeKernel = program->createKernel("selectLangevinStepSize");
    blockSize = min(256, system.getNumParticles());
    blockSize = max(blockSize, (int) params.getSize());
}

double CommonIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    bool useDouble = cc.getUseDoublePrecision() || cc.getUseMixedPrecision();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getStepSize());
        kernel2->addArg(numAtoms);
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(integration.getPosDelta());
        kernel2->addArg(oldDelta);
        kernel2->addArg(params);
        kernel2->addArg(integration.getStepSize());
        kernel2->addArg(integration.getRandom());
        kernel2->addArg(); // Random index will be set just before it is executed.
        kernel3->addArg(numAtoms);
        kernel3->addArg(cc.getPosq());
        kernel3->addArg(cc.getVelm());
        kernel3->addArg(integration.getPosDelta());
        kernel3->addArg(oldDelta);
        kernel3->addArg(integration.getStepSize());
        if (cc.getUseMixedPrecision())
            kernel3->addArg(cc.getPosqCorrection());
        selectSizeKernel->addArg(numAtoms);
        selectSizeKernel->addArg(paddedNumAtoms);
        for (int i = 0; i < 4; i++)
            selectSizeKernel->addArg();
        selectSizeKernel->addArg(integration.getStepSize());
        selectSizeKernel->addArg(cc.getVelm());
        selectSizeKernel->addArg(cc.getLongForceBuffer());
        selectSizeKernel->addArg(params);
    }

    // Select the step size to use.

    double maxStepSize = maxTime-cc.getTime();
    if (integrator.getMaximumStepSize() > 0)
        maxStepSize = min(integrator.getMaximumStepSize(), maxStepSize);
    float maxStepSizeFloat = (float) maxStepSize;
    if (useDouble) {
        selectSizeKernel->setArg(2, maxStepSize);
        selectSizeKernel->setArg(3, integrator.getErrorTolerance());
        selectSizeKernel->setArg(4, integrator.getFriction());
        selectSizeKernel->setArg(5, BOLTZ*integrator.getTemperature());
    }
    else {
        selectSizeKernel->setArg(2, maxStepSizeFloat);
        selectSizeKernel->setArg(3, (float) integrator.getErrorTolerance());
        selectSizeKernel->setArg(4, (float) integrator.getFriction());
        selectSizeKernel->setArg(5, (float) (BOLTZ*integrator.getTemperature()));
    }
    selectSizeKernel->execute(blockSize, blockSize);

    // Perform the integration.

    kernel2->setArg(7, integration.prepareRandomNumbers(cc.getPaddedNumAtoms()));
    kernel1->execute(numAtoms);
    integration.applyVelocityConstraints(integrator.getConstraintTolerance());
    kernel2->execute(numAtoms);
    integration.applyConstraints(integrator.getConstraintTolerance());
    kernel3->execute(numAtoms);
    integration.computeVirtualSites();
    
    // Reduce UI lag.

    flushPeriodically(cc);

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

double CommonIntegrateVariableLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}

class CommonIntegrateDPDStepKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(ComputeContext& cc, vector<int>& particleTypes, ComputeArray& particleTypeArray) :
            cc(cc), particleTypes(particleTypes), particleTypeArray(particleTypeArray) {
    }
    void execute() {
        // Reorder particleTypes to reflect the new atom order.

        vector<int> sortedTypes(particleTypes.size());
        const vector<int>& order = cc.getAtomIndex();
        for (int i = 0; i < particleTypes.size(); i++)
            sortedTypes[i] = particleTypes[order[i]];
        particleTypeArray.upload(sortedTypes);
    }
private:
    ComputeContext& cc;
    ComputeArray& particleTypeArray;
    vector<int> particleTypes;
};

void CommonIntegrateDPDStepKernel::initialize(const System& system, const DPDIntegrator& integrator) {
    // Record information about the integrator.

    vector<int> particleTypeVec;
    vector<vector<double> > frictionTable, cutoffTable;
    DPDIntegratorUtilities::createTypeTables(integrator, system.getNumParticles(), numTypes, particleTypeVec, frictionTable, cutoffTable, maxCutoff);
    while (particleTypeVec.size() < cc.getPaddedNumAtoms())
        particleTypeVec.push_back(0);

    // We want the NonbondedUtilities to build a neighbor list.  Add an empty interaction
    // with a nonexistant force group to it.

    cc.getNonbondedUtilities().addInteraction(true, system.usesPeriodicBoundaryConditions(), false, maxCutoff, vector<vector<int> >(), "", 32);

    // Create the arrays.

    cc.initializeContexts();
    ContextSelector selector(cc);
    if (cc.getUseDoublePrecision() || cc.getUseMixedPrecision())
        oldDelta.initialize<mm_double4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    else
        oldDelta.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "oldDelta");
    particleType.initialize<int>(cc, cc.getPaddedNumAtoms(), "dpdParticleType");
    pairParams.initialize<mm_float2>(cc, numTypes*numTypes, "dpdPairParams");
    velDelta.initialize<long long>(cc, 3*cc.getPaddedNumAtoms(), "velDelta");
    tileCounter.initialize<int>(cc, 1, "tileCounter");
    particleType.upload(particleTypeVec);
    vector<mm_float2> pairParamsVec(numTypes*numTypes);
    for (int i = 0; i < numTypes; i++)
        for (int j = 0; j < numTypes; j++)
            pairParamsVec[i+j*numTypes] = mm_float2(frictionTable[i][j], cutoffTable[i][j]);
    pairParams.upload(pairParamsVec);
    cc.addAutoclearBuffer(velDelta);
    cc.addReorderListener(new ReorderListener(cc, particleTypeVec, particleType));
    randomSeed = integrator.getRandomNumberSeed();
    if (randomSeed == 0)
        randomSeed = osrngseed(); // A seed of 0 means use a unique one
}

void CommonIntegrateDPDStepKernel::execute(ContextImpl& context, const DPDIntegrator& integrator) {
    ContextSelector selector(cc);
    IntegrationUtilities& integration = cc.getIntegrationUtilities();
    NonbondedUtilities& nb = cc.getNonbondedUtilities();
    int numAtoms = cc.getNumAtoms();
    int paddedNumAtoms = cc.getPaddedNumAtoms();
    if (!hasInitializedKernels) {
        map<string, string> defines;
        defines["M_PI"] = cc.doubleToString(M_PI);
        defines["MAX_CUTOFF"] = cc.doubleToString(maxCutoff);
        defines["TILE_SIZE"] = cc.intToString(ComputeContext::TileSize);
        blockSize = max(32, cc.getNonbondedUtilities().getForceThreadBlockSize());
        defines["WORK_GROUP_SIZE"] = cc.intToString(blockSize);
        if (context.getSystem().usesPeriodicBoundaryConditions())
            defines["USE_PERIODIC"] = "1";
        if (cc.getIsCPU())
            defines["INTEL_WORKAROUND"] = "1";
        try {
            // Workaround for errors in NVIDIA OpenCL.

            string platformName = context.getPlatform().getPropertyValue(context.getOwner(), "OpenCLPlatformName");
            if (platformName.rfind("NVIDIA", 0) == 0)
                defines["NVIDIA_WORKAROUND"] = "1";
        }
        catch (...) {
            // This isn't the OpenCL platform.
        }
        ComputeProgram program = cc.compileProgram(CommonKernelSources::dpd, defines);
        kernel1 = program->createKernel("integrateDPDPart1");
        kernel2 = program->createKernel("integrateDPDPart2");
        kernel3 = program->createKernel("integrateDPDPart3");
        kernel4 = program->createKernel("integrateDPDPart4");
        hasInitializedKernels = true;
        kernel1->addArg(numAtoms);
        kernel1->addArg(paddedNumAtoms);
        kernel1->addArg(cc.getVelm());
        kernel1->addArg(cc.getLongForceBuffer());
        kernel1->addArg(integration.getStepSize());
        kernel1->addArg(tileCounter);
        kernel2->addArg(numAtoms);
        kernel2->addArg(paddedNumAtoms);
        kernel2->addArg(cc.getPosq());
        kernel2->addArg(cc.getVelm());
        kernel2->addArg(velDelta);
        kernel2->addArg(integration.getStepSize());
        kernel2->addArg(particleType);
        kernel2->addArg(numTypes);
        kernel2->addArg(pairParams);
        for (int i = 0; i < 7; i++)
            kernel2->addArg(); // Random seed, kT, and box vectors will be set just before it is executed.
        kernel2->addArg(nb.getExclusionTiles());
        kernel2->addArg((int) nb.getExclusionTiles().getSize());
        kernel2->addArg(nb.getInteractingTiles());
        kernel2->addArg(nb.getInteractionCount());
        kernel2->addArg(nb.getBlockCenters());
        kernel2->addArg(nb.getBlockBoundingBoxes());
        kernel2->addArg(nb.getInteractingAtoms());
        kernel2->addArg(tileCounter);
        if (cc.getUseMixedPrecision())
            kernel2->addArg(cc.getPosqCorrection());
        kernel3->addArg(numAtoms);
        kernel3->addArg(paddedNumAtoms);
        kernel3->addArg(cc.getVelm());
        kernel3->addArg(velDelta);
        kernel3->addArg(integration.getPosDelta());
        kernel3->addArg(oldDelta);
        kernel3->addArg(integration.getStepSize());
        kernel4->addArg(numAtoms);
        kernel4->addArg(cc.getPosq());
        kernel4->addArg(cc.getVelm());
        kernel4->addArg(integration.getPosDelta());
        kernel4->addArg(oldDelta);
        kernel4->addArg(integration.getStepSize());
        if (cc.getUseMixedPrecision())
            kernel4->addArg(cc.getPosqCorrection());
    }

    // Perform the integration.

    vector<int> p;
    particleType.download(p);
    double stepSize = integrator.getStepSize();
    cc.getIntegrationUtilities().setNextStepSize(stepSize);
    kernel1->execute(numAtoms);
    particleType.download(p);
    integration.applyVelocityConstraints(integrator.getConstraintTolerance());
    particleType.download(p);
    kernel2->setArg(9, randomSeed+cc.getStepCount());
    kernel2->setArg(10, (float) (BOLTZ*integrator.getTemperature()));
    setPeriodicBoxArgs(cc, kernel2, 11);
    particleType.download(p);
    kernel2->execute(2*nb.getNumForceThreadBlocks()*blockSize, blockSize);
    particleType.download(p);
    kernel3->execute(numAtoms);
    particleType.download(p);
    integration.applyConstraints(integrator.getConstraintTolerance());
    particleType.download(p);
    kernel4->execute(numAtoms);
    particleType.download(p);
    integration.computeVirtualSites();

    // Update the time and step count.

    cc.setTime(cc.getTime()+stepSize);
    cc.setStepCount(cc.getStepCount()+1);
    cc.reorderAtoms();

    // Reduce UI lag.

    flushPeriodically(cc);
    particleType.download(p);
}

double CommonIntegrateDPDStepKernel::computeKineticEnergy(ContextImpl& context, const DPDIntegrator& integrator) {
    return cc.getIntegrationUtilities().computeKineticEnergy(0.0);
}

void CommonRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    ContextSelector selector(cc);
    frequency = force.getFrequency();
    int numAtoms = cc.getNumAtoms();
    cmMomentum.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "cmMomentum");
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
    ContextSelector selector(cc);
    kernel1->execute(cc.getNumAtoms(), 64);
    kernel2->execute(cc.getNumAtoms(), 64);
}

class CommonCalcRMSDForceKernel::ForceInfo : public ComputeForceInfo {
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

void CommonCalcRMSDForceKernel::initialize(const System& system, const RMSDForce& force) {
    // Create data structures.
    
    ContextSelector selector(cc);
    bool useDouble = cc.getUseDoublePrecision();
    int elementSize = (useDouble ? sizeof(double) : sizeof(float));
    int numParticles = force.getParticles().size();
    if (numParticles == 0)
        numParticles = system.getNumParticles();
    referencePos.initialize(cc, system.getNumParticles(), 4*elementSize, "referencePos");
    particles.initialize<int>(cc, numParticles, "particles");
    buffer.initialize(cc, 13, elementSize, "buffer");
    recordParameters(force);
    info = new ForceInfo(force);
    cc.addForce(info);
    
    // Create the kernels.

    blockSize = min(256, cc.getMaxThreadBlockSize());
    map<string, string> defines;
    defines["THREAD_BLOCK_SIZE"] = cc.intToString(blockSize);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::rmsd, defines);
    kernel1 = program->createKernel("computeRMSDPart1");
    kernel2 = program->createKernel("computeRMSDForces");
    kernel1->addArg();
    kernel1->addArg(cc.getPosq());
    kernel1->addArg(referencePos);
    kernel1->addArg(particles);
    kernel1->addArg(buffer);
    kernel2->addArg();
    kernel2->addArg(cc.getPaddedNumAtoms());
    kernel2->addArg(cc.getPosq());
    kernel2->addArg(referencePos);
    kernel2->addArg(particles);
    kernel2->addArg(buffer);
    kernel2->addArg(cc.getLongForceBuffer());
}

void CommonCalcRMSDForceKernel::recordParameters(const RMSDForce& force) {
    // Record the parameters and center the reference positions.
    
    vector<int> particleVec = force.getParticles();
    if (particleVec.size() == 0)
        for (int i = 0; i < cc.getNumAtoms(); i++)
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

double CommonCalcRMSDForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    ContextSelector selector(cc);
    if (cc.getUseDoublePrecision())
        return executeImpl<double>(context);
    return executeImpl<float>(context);
}

template <class REAL>
double CommonCalcRMSDForceKernel::executeImpl(ContextImpl& context) {
    // Execute the first kernel.

    int numParticles = particles.getSize();
    kernel1->setArg(0, numParticles);
    kernel1->execute(blockSize, blockSize);
    
    // Download the results, build the F matrix, and find the maximum eigenvalue
    // and eigenvector.

    vector<REAL> b;
    buffer.download(b);

    // JAMA::Eigenvalue may run into an infinite loop if we have any NaN
    for (int i = 0; i < 9; i++) {
        if (b[i] != b[i])
            throw OpenMMException("NaN encountered during RMSD force calculation");
    }
    
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
    kernel2->setArg(0, numParticles);
    kernel2->execute(numParticles);
    return rmsd;
}

void CommonCalcRMSDForceKernel::copyParametersToContext(ContextImpl& context, const RMSDForce& force) {
    ContextSelector selector(cc);
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
    cc.invalidateMolecules(info);
}

void CommonApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    ContextSelector selector(cc);
    randomSeed = thermostat.getRandomNumberSeed();
    ComputeProgram program = cc.compileProgram(CommonKernelSources::andersenThermostat);
    kernel = program->createKernel("applyAndersenThermostat");
    cc.getIntegrationUtilities().initRandomNumberGenerator(randomSeed);

    // Create the arrays with the group definitions.

    vector<vector<int> > groups = AndersenThermostatImpl::calcParticleGroups(system);
    atomGroups.initialize<int>(cc, cc.getNumAtoms(), "atomGroups");
    vector<int> atoms(atomGroups.getSize());
    for (int i = 0; i < (int) groups.size(); i++) {
        for (int j = 0; j < (int) groups[i].size(); j++)
            atoms[groups[i][j]] = i;
    }
    atomGroups.upload(atoms);
    kernel->addArg(system.getNumParticles());
    kernel->addArg();
    kernel->addArg();
    kernel->addArg(cc.getVelm());
    kernel->addArg();
    kernel->addArg(cc.getIntegrationUtilities().getRandom());
    kernel->addArg();
    kernel->addArg(atomGroups);
}

void CommonApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    ContextSelector selector(cc);
    kernel->setArg(1, (float) context.getParameter(AndersenThermostat::CollisionFrequency()));
    kernel->setArg(2, (float) (BOLTZ*context.getParameter(AndersenThermostat::Temperature())));
    double stepSize = context.getIntegrator().getStepSize();
    if (cc.getUseDoublePrecision())
        kernel->setArg(4, stepSize);
    else
        kernel->setArg(4, (float) stepSize);
    kernel->setArg(6, cc.getIntegrationUtilities().prepareRandomNumbers(cc.getPaddedNumAtoms()));
    kernel->execute(cc.getNumAtoms());
}

void CommonApplyMonteCarloBarostatKernel::initialize(const System& system, const Force& thermostat, int components, bool rigidMolecules) {
    this->components = components;
    this->rigidMolecules = rigidMolecules;
    ContextSelector selector(cc);
    savedPositions.initialize(cc, cc.getPaddedNumAtoms(), cc.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4), "savedPositions");
    savedVelocities.initialize(cc, cc.getPaddedNumAtoms(), cc.getUseDoublePrecision() || cc.getUseMixedPrecision() ? sizeof(mm_double4) : sizeof(mm_float4), "savedVelocities");
    savedLongForces.initialize<long long>(cc, cc.getPaddedNumAtoms()*3, "savedLongForces");
    try {
        cc.getFloatForceBuffer(); // This will throw an exception on the CUDA platform.
        savedFloatForces.initialize(cc, cc.getPaddedNumAtoms(), cc.getUseDoublePrecision() ? sizeof(mm_double4) : sizeof(mm_float4), "savedForces");
    }
    catch (...) {
        // The CUDA platform doesn't have a floating point force buffer, so we don't need to copy it.
    }
    energyBuffers.resize(components);
    for (int i = 0; i < components; i++)
        energyBuffers[i].initialize(cc, cc.getNumThreadBlocks(), cc.getUseDoublePrecision() || cc.getUseMixedPrecision() ? sizeof(double) : sizeof(float), "energyBuffer");
    vector<float> zeros(energyBuffers[0].getSize(), 0.0f);
    for (int i = 0; i < components; i++)
        energyBuffers[i].upload(zeros, true);
    map<string, string> defines;
    defines["WORK_GROUP_SIZE"] = cc.intToString(cc.ThreadBlockSize);
    defines["COMPONENTS"] = cc.intToString(components);
    ComputeProgram program = cc.compileProgram(CommonKernelSources::monteCarloBarostat, defines);
    kernel = program->createKernel("scalePositions");
    kineticEnergyKernel = program->createKernel("computeMolecularKineticEnergy");
}

void CommonApplyMonteCarloBarostatKernel::saveCoordinates(ContextImpl& context) {
    ContextSelector selector(cc);
    cc.getPosq().copyTo(savedPositions);
    cc.getVelm().copyTo(savedVelocities);
    cc.getLongForceBuffer().copyTo(savedLongForces);
    if (savedFloatForces.isInitialized())
        cc.getFloatForceBuffer().copyTo(savedFloatForces);
    lastPosCellOffsets = cc.getPosCellOffsets();
    lastAtomOrder = cc.getAtomIndex();
}

void CommonApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ) {
    ContextSelector selector(cc);

    // check if atoms were reordered from energy evaluation before scaling
    atomsWereReordered = cc.getAtomsWereReordered();

    if (!hasInitializedKernels) {
        hasInitializedKernels = true;

        // Create the arrays with the molecule definitions.

        vector<vector<int> > molecules;
        if (rigidMolecules)
            molecules = context.getMolecules();
        else {
            molecules.resize(cc.getNumAtoms());
            for (int i = 0; i < molecules.size(); i++)
                molecules[i].push_back(i);
        }
        numMolecules = molecules.size();
        moleculeAtoms.initialize<int>(cc, cc.getNumAtoms(), "moleculeAtoms");
        moleculeStartIndex.initialize<int>(cc, numMolecules+1, "moleculeStartIndex");
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

        kernel->addArg();
        kernel->addArg();
        kernel->addArg();
        kernel->addArg(numMolecules);
        for (int i = 0; i < 5; i++)
            kernel->addArg();
        kernel->addArg(cc.getPosq());
        kernel->addArg(moleculeAtoms);
        kernel->addArg(moleculeStartIndex);
        kineticEnergyKernel->addArg(numMolecules);
        kineticEnergyKernel->addArg(cc.getVelm());
        kineticEnergyKernel->addArg(moleculeAtoms);
        kineticEnergyKernel->addArg(moleculeStartIndex);
        for (int i = 0; i < components; i++)
            kineticEnergyKernel->addArg(energyBuffers[i]);
    }
    kernel->setArg(0, (float) scaleX);
    kernel->setArg(1, (float) scaleY);
    kernel->setArg(2, (float) scaleZ);
    setPeriodicBoxArgs(cc, kernel, 4);
    kernel->execute(numMolecules);
}

void CommonApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    ContextSelector selector(cc);
    savedPositions.copyTo(cc.getPosq());
    savedVelocities.copyTo(cc.getVelm());
    savedLongForces.copyTo(cc.getLongForceBuffer());
    cc.setPosCellOffsets(lastPosCellOffsets);
    if (savedFloatForces.isInitialized())
        savedFloatForces.copyTo(cc.getFloatForceBuffer());

    // check if atoms were reordered from energy evaluation before or after scaling
    if (atomsWereReordered || cc.getAtomsWereReordered())
        cc.setAtomIndex(lastAtomOrder);
}

void CommonApplyMonteCarloBarostatKernel::computeKineticEnergy(ContextImpl& context, vector<double>& ke) {
    ContextSelector selector(cc);
    ke.resize(components);
    kineticEnergyKernel->execute(numMolecules);
    for (int j = 0; j < components; j++) {
        ke[j] = 0.0;
        if (energyBuffers[j].getElementSize() == sizeof(float)) {
            vector<float> buffer;
            energyBuffers[j].download(buffer);
            for (int i = 0; i < buffer.size(); i++)
                ke[j] += buffer[i];
        }
        else {
            vector<double> buffer;
            energyBuffers[j].download(buffer);
            for (int i = 0; i < buffer.size(); i++)
                ke[j] += buffer[i];
        }
    }
}

class CommonCalcATMForceKernel::ReorderListener : public ComputeContext::ReorderListener {
public:
    ReorderListener(ComputeContext& cc, ArrayInterface& invAtomOrder) : cc(cc), invAtomOrder(invAtomOrder) {
    }
    void execute() {
        vector<int> invOrder(cc.getPaddedNumAtoms());
        const vector<int>& order = cc.getAtomIndex();
        for (int i = 0; i < order.size(); i++)
            invOrder[order[i]] = i;
        invAtomOrder.upload(invOrder);
    }
private:
    ComputeContext& cc;
    ArrayInterface& invAtomOrder;
};

CommonCalcATMForceKernel::~CommonCalcATMForceKernel() {
}

void CommonCalcATMForceKernel::initialize(const System& system, const ATMForce& force) {
    ContextSelector selector(cc);
    numParticles = force.getNumParticles();
    if (numParticles == 0)
        return;
    vector<mm_float4> displVector1(cc.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    vector<mm_float4> displVector0(cc.getPaddedNumAtoms(), mm_float4(0, 0, 0, 0));
    for (int i = 0; i < numParticles; i++) {
        Vec3 displacement1, displacement0;
        force.getParticleParameters(i, displacement1, displacement0);
        displVector1[i] = mm_float4(displacement1[0], displacement1[1], displacement1[2], 0);
        displVector0[i] = mm_float4(displacement0[0], displacement0[1], displacement0[2], 0);
    }
    displ1.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "displ1");
    displ1.upload(displVector1);
    displ0.initialize<mm_float4>(cc, cc.getPaddedNumAtoms(), "displ0");
    displ0.upload(displVector0);
    invAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "invAtomOrder");
    inner0InvAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "inner0InvAtomOrder");
    inner1InvAtomOrder.initialize<int>(cc, cc.getPaddedNumAtoms(), "inner1InvAtomOrder");
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        cc.addEnergyParameterDerivative(force.getEnergyParameterDerivativeName(i));
}

void CommonCalcATMForceKernel::initKernels(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1) {
    if (!hasInitializedKernel) {
        hasInitializedKernel = true;

        //inner contexts
        ComputeContext& cc0 = getInnerComputeContext(innerContext0);
        ComputeContext& cc1 = getInnerComputeContext(innerContext1);

        // Copy positions to the inner contexts.
        vector<Vec3> positions;
        context.getPositions(positions);
        innerContext0.setPositions(positions);
        innerContext1.setPositions(positions);

        // Initialize the listeners.
        ReorderListener* listener = new ReorderListener(cc, invAtomOrder);
        ReorderListener* listener0 = new ReorderListener(cc0, inner0InvAtomOrder);
        ReorderListener* listener1 = new ReorderListener(cc1, inner1InvAtomOrder);
        cc.addReorderListener(listener);
        cc0.addReorderListener(listener0);
        cc1.addReorderListener(listener1);
        listener->execute();
        listener0->execute();
        listener1->execute();

        //create CopyState kernel
        ComputeProgram program = cc.compileProgram(CommonKernelSources::atmforce);
        copyStateKernel = program->createKernel("copyState");
        copyStateKernel->addArg(numParticles);
        copyStateKernel->addArg(cc.getPosq());
        copyStateKernel->addArg(cc0.getPosq());
        copyStateKernel->addArg(cc1.getPosq());
        copyStateKernel->addArg(displ0);
        copyStateKernel->addArg(displ1);
        copyStateKernel->addArg(cc.getAtomIndexArray());
        copyStateKernel->addArg(inner0InvAtomOrder);
        copyStateKernel->addArg(inner1InvAtomOrder);
        if (cc.getUseMixedPrecision()) {
            copyStateKernel->addArg(cc.getPosqCorrection());
            copyStateKernel->addArg(cc0.getPosqCorrection());
            copyStateKernel->addArg(cc1.getPosqCorrection());
        }

        //create the HybridForce kernel
        hybridForceKernel = program->createKernel("hybridForce");
        hybridForceKernel->addArg(numParticles);
        hybridForceKernel->addArg(cc.getPaddedNumAtoms());
        hybridForceKernel->addArg(cc.getLongForceBuffer());
        hybridForceKernel->addArg(cc0.getLongForceBuffer());
        hybridForceKernel->addArg(cc1.getLongForceBuffer());
        hybridForceKernel->addArg(invAtomOrder);
        hybridForceKernel->addArg(inner0InvAtomOrder);
        hybridForceKernel->addArg(inner1InvAtomOrder);
        hybridForceKernel->addArg();
        hybridForceKernel->addArg();

        cc0.addForce(new ComputeForceInfo());
        cc1.addForce(new ComputeForceInfo());

    }
}

void CommonCalcATMForceKernel::applyForces(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1,
        double dEdu0, double dEdu1, const map<string, double>& energyParamDerivs) {
    ContextSelector selector(cc);
    initKernels(context, innerContext0, innerContext1);
    if (cc.getUseDoublePrecision()) {
        hybridForceKernel->setArg(8, dEdu0);
        hybridForceKernel->setArg(9, dEdu1);
    }
    else {
        hybridForceKernel->setArg(8, (float) dEdu0);
        hybridForceKernel->setArg(9, (float) dEdu1);
    }
    hybridForceKernel->execute(numParticles);
    map<string, double>& derivs = cc.getEnergyParamDerivWorkspace();
    for (auto deriv : energyParamDerivs)
        derivs[deriv.first] += deriv.second;
}

void CommonCalcATMForceKernel::copyState(ContextImpl& context,
        ContextImpl& innerContext0, ContextImpl& innerContext1) {
    ContextSelector selector(cc);

    initKernels(context, innerContext0, innerContext1);

    ComputeContext& cc0 = getInnerComputeContext(innerContext0);
    ComputeContext& cc1 = getInnerComputeContext(innerContext1);
    cc0.reorderAtoms();
    cc1.reorderAtoms();
    copyStateKernel->execute(numParticles);

    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext0.setPeriodicBoxVectors(a, b, c);
    innerContext0.setTime(context.getTime());
    innerContext1.setPeriodicBoxVectors(a, b, c);
    innerContext1.setTime(context.getTime());
    map<string, double> innerParameters0 = innerContext0.getParameters();
    for (auto& param : innerParameters0)
        innerContext0.setParameter(param.first, context.getParameter(param.first));
    map<string, double> innerParameters1 = innerContext1.getParameters();
    for (auto& param : innerParameters1)
        innerContext1.setParameter(param.first, context.getParameter(param.first));
}

void CommonCalcATMForceKernel::copyParametersToContext(ContextImpl& context, const ATMForce& force) {
    ContextSelector selector(cc);
    if (force.getNumParticles() != numParticles)
        throw OpenMMException("copyParametersToContext: The number of ATMMetaForce particles has changed");
    vector<mm_float4> displVector1(cc.getPaddedNumAtoms());
    vector<mm_float4> displVector0(cc.getPaddedNumAtoms());
    for (int i = 0; i < numParticles; i++) {
        Vec3 displacement1, displacement0;
        force.getParticleParameters(i, displacement1, displacement0);
        displVector1[i] = mm_float4(displacement1[0], displacement1[1], displacement1[2], 0);
        displVector0[i] = mm_float4(displacement0[0], displacement0[1], displacement0[2], 0);
    }
    displ1.upload(displVector1);
    displ0.upload(displVector0);
}

class CommonCalcCustomCPPForceKernel::StartCalculationPreComputation : public ComputeContext::ForcePreComputation {
public:
    StartCalculationPreComputation(CommonCalcCustomCPPForceKernel& owner) : owner(owner) {
    }
    void computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        owner.beginComputation(includeForces, includeEnergy, groups);
    }
    CommonCalcCustomCPPForceKernel& owner;
};

class CommonCalcCustomCPPForceKernel::ExecuteTask : public ComputeContext::WorkTask {
public:
    ExecuteTask(CommonCalcCustomCPPForceKernel& owner, bool includeForces) : owner(owner), includeForces(includeForces) {
    }
    void execute() {
        owner.executeOnWorkerThread(includeForces);
    }
    CommonCalcCustomCPPForceKernel& owner;
    bool includeForces;
};

class CommonCalcCustomCPPForceKernel::AddForcesPostComputation : public ComputeContext::ForcePostComputation {
public:
    AddForcesPostComputation(CommonCalcCustomCPPForceKernel& owner) : owner(owner) {
    }
    double computeForceAndEnergy(bool includeForces, bool includeEnergy, int groups) {
        return owner.addForces(includeForces, includeEnergy, groups);
    }
    CommonCalcCustomCPPForceKernel& owner;
};

void CommonCalcCustomCPPForceKernel::initialize(const System& system, CustomCPPForceImpl& force) {
    ContextSelector selector(cc);
    this->force = &force;
    int numParticles = system.getNumParticles();
    forcesVec.resize(numParticles);
    positionsVec.resize(numParticles);
    floatForces.resize(3*numParticles);
    int elementSize = (cc.getUseDoublePrecision() ? sizeof(double) : sizeof(float));
    forcesArray.initialize(cc, 3*numParticles, elementSize, "forces");
    map<string, string> defines;
    defines["NUM_ATOMS"] = cc.intToString(numParticles);
    defines["PADDED_NUM_ATOMS"] = cc.intToString(cc.getPaddedNumAtoms());
    ComputeProgram program = cc.compileProgram(CommonKernelSources::customCppForce, defines);
    addForcesKernel = program->createKernel("addForces");
    addForcesKernel->addArg(forcesArray);
    addForcesKernel->addArg(cc.getLongForceBuffer());
    addForcesKernel->addArg(cc.getAtomIndexArray());
    forceGroupFlag = (1<<force.getOwner().getForceGroup());
    if (cc.getNumContexts() == 1) {
        cc.addPreComputation(new StartCalculationPreComputation(*this));
        cc.addPostComputation(new AddForcesPostComputation(*this));
    }
}

double CommonCalcCustomCPPForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (cc.getNumContexts() == 1) {
        // This method does nothing.  The actual calculation is started by the pre-computation, continued on
        // the worker thread, and finished by the post-computation.

        return 0;
    }

    // When using multiple GPUs, this method is itself called from the worker thread.
    // Submitting additional tasks and waiting for them to complete would lead to
    // a deadlock.

    if (cc.getContextIndex() != 0)
        return 0.0;
    contextImpl.getPositions(positionsVec);
    executeOnWorkerThread(includeForces);
    return addForces(includeForces, includeEnergy, -1);
}

void CommonCalcCustomCPPForceKernel::beginComputation(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return;
    contextImpl.getPositions(positionsVec);
    
    // The actual force computation will be done on a different thread.
    
    cc.getWorkThread().addTask(new ExecuteTask(*this, includeForces));
}

void CommonCalcCustomCPPForceKernel::executeOnWorkerThread(bool includeForces) {
    energy = force->computeForce(contextImpl, positionsVec, forcesVec);
    if (includeForces) {
        ContextSelector selector(cc);
        int numParticles = cc.getNumAtoms();
        if (cc.getUseDoublePrecision())
            forcesArray.upload((double*) forcesVec.data());
        else {
            for (int i = 0; i < numParticles; i++) {
                floatForces[3*i] = (float) forcesVec[i][0];
                floatForces[3*i+1] = (float) forcesVec[i][1];
                floatForces[3*i+2] = (float) forcesVec[i][2];
            }
            forcesArray.upload(floatForces);
        }
    }
}

double CommonCalcCustomCPPForceKernel::addForces(bool includeForces, bool includeEnergy, int groups) {
    if ((groups&forceGroupFlag) == 0)
        return 0;

    // Wait until executeOnWorkerThread() is finished.
    
    if (cc.getNumContexts() == 1)
        cc.getWorkThread().flush();

    // Add in the forces.
    
    if (includeForces) {
        ContextSelector selector(cc);
        addForcesKernel->execute(cc.getNumAtoms());
    }
    
    // Return the energy.
    
    return energy;
}