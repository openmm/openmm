/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/internal/NonbondedForceImpl.h"
#include "CudaBondedUtilities.h"
#include "CudaExpressionUtilities.h"
#include "CudaIntegrationUtilities.h"
//#include "CudaNonbondedUtilities.h"
#include "CudaKernelSources.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "../src/SimTKUtilities/SimTKOpenMMUtilities.h"
#include <cmath>
#include <set>

using namespace OpenMM;
using namespace std;
using Lepton::ExpressionTreeNode;
using Lepton::Operation;

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
    cuCtxSetCurrent(cu.getContext());
//    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
//    bool includeNonbonded = ((groups&(1<<nb.getForceGroup())) != 0);
//    cu.setAtomsWereReordered(false);
//    if (nb.getUseCutoff() && includeNonbonded && (cu.getMoleculesAreInvalid() || cu.getComputeForceCount()%100 == 0)) {
//        cu.reorderAtoms(!cu.getMoleculesAreInvalid());
//        nb.updateNeighborListSize();
//    }
    cu.setComputeForceCount(cu.getComputeForceCount()+1);
    cu.clearAutoclearBuffers();
//    if (includeNonbonded)
//        nb.prepareInteractions();
}

double CudaCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    cu.getBondedUtilities().computeInteractions(groups);
//    if ((groups&(1<<cu.getNonbondedUtilities().getForceGroup())) != 0)
//        cu.getNonbondedUtilities().computeInteractions();
    cu.getIntegrationUtilities().distributeForcesFromVirtualSites();
    double sum = 0.0;
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
    cuCtxSetCurrent(cu.getContext());
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
    cuCtxSetCurrent(cu.getContext());
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
            pos.x = p[0];
            pos.y = p[1];
            pos.z = p[2];
        }
        for (int i = numParticles; i < cu.getPaddedNumAtoms(); i++)
            posq[i] = make_float4(0.0, 0.0, 0.0, 0.0);
        cu.getPosq().upload(posq);
    }
    for (int i = 0; i < (int) cu.getPosCellOffsets().size(); i++)
        cu.getPosCellOffsets()[i] = make_int4(0, 0, 0, 0);
}

void CudaUpdateStateDataKernel::getVelocities(ContextImpl& context, vector<Vec3>& velocities) {
    cuCtxSetCurrent(cu.getContext());
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    if (cu.getUseDoublePrecision()) {
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
    cuCtxSetCurrent(cu.getContext());
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    if (cu.getUseDoublePrecision()) {
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
            velm[i] = make_float4(0.0, 0.0, 0.0, 0.0);
        cu.getVelm().upload(velm);
    }
}

void CudaUpdateStateDataKernel::getForces(ContextImpl& context, vector<Vec3>& forces) {
    cuCtxSetCurrent(cu.getContext());
    long long* force = (long long*) cu.getPinnedBuffer();
    cu.getForce().download(force);
    const vector<int>& order = cu.getAtomIndex();
    int numParticles = context.getSystem().getNumParticles();
    int paddedNumParticles = cu.getPaddedNumAtoms();
    forces.resize(numParticles);
    double scale = 1.0/(double) 0xFFFFFFFF;
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
    cuCtxSetCurrent(cu.getContext());
//    int version = 1;
//    stream.write((char*) &version, sizeof(int));
//    double time = cu.getTime();
//    stream.write((char*) &time, sizeof(double));
//    cu.getPosq().download();
//    stream.write((char*) &cu.getPosq()[0], sizeof(mm_float4)*cu.getPosq().getSize());
//    cu.getVelm().download();
//    stream.write((char*) &cu.getVelm()[0], sizeof(mm_float4)*cu.getVelm().getSize());
//    stream.write((char*) &cu.getAtomIndex()[0], sizeof(cl_int)*cu.getAtomIndex().getSize());
//    stream.write((char*) &cu.getPosCellOffsets()[0], sizeof(mm_int4)*cu.getPosCellOffsets().size());
//    mm_float4 box = cu.getPeriodicBoxSize();
//    stream.write((char*) &box, sizeof(mm_float4));
//    cu.getIntegrationUtilities().createCheckpoint(stream);
//    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void CudaUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    cuCtxSetCurrent(cu.getContext());
//    int version;
//    stream.read((char*) &version, sizeof(int));
//    if (version != 1)
//        throw OpenMMException("Checkpoint was created with a different version of OpenMM");
//    double time;
//    stream.read((char*) &time, sizeof(double));
//    vector<CudaContext*>& contexts = cu.getPlatformData().contexts;
//    for (int i = 0; i < (int) contexts.size(); i++)
//        contexts[i]->setTime(time);
//    stream.read((char*) &cu.getPosq()[0], sizeof(mm_float4)*cu.getPosq().getSize());
//    cu.getPosq().upload();
//    stream.read((char*) &cu.getVelm()[0], sizeof(mm_float4)*cu.getVelm().getSize());
//    cu.getVelm().upload();
//    stream.read((char*) &cu.getAtomIndex()[0], sizeof(cl_int)*cu.getAtomIndex().getSize());
//    cu.getAtomIndex().upload();
//    stream.read((char*) &cu.getPosCellOffsets()[0], sizeof(mm_int4)*cu.getPosCellOffsets().size());
//    mm_float4 box;
//    stream.read((char*) &box, sizeof(mm_float4));
//    for (int i = 0; i < (int) contexts.size(); i++)
//        contexts[i]->setPeriodicBoxSize(box.x, box.y, box.z);
//    cu.getIntegrationUtilities().loadCheckpoint(stream);
//    SimTKOpenMMUtilities::loadCheckpoint(stream);
//    for (int i = 0; i < cu.getReorderListeners().size(); i++)
//        cu.getReorderListeners()[i]->execute();
}

void CudaApplyConstraintsKernel::initialize(const System& system) {
}

void CudaApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
//    if (!hasInitializedKernel) {
//        hasInitializedKernel = true;
//        map<string, string> defines;
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        cu::Program program = cu.createProgram(CudaKernelSources::constraints, defines);
//        applyDeltasKernel = cu::Kernel(program, "applyPositionDeltas");
//        applyDeltasKernel.setArg<cu::Buffer>(0, cu.getPosq().getDeviceBuffer());
//        applyDeltasKernel.setArg<cu::Buffer>(1, cu.getIntegrationUtilities().getPosDelta().getDeviceBuffer());
//    }
//    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
//    cu.clearBuffer(integration.getPosDelta());
//    integration.applyConstraints(tol);
//    cu.executeKernel(applyDeltasKernel, cu.getNumAtoms());
//    integration.computeVirtualSites();
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
    if (params != NULL)
        delete params;
}

void CudaCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = CudaArray::create<float2>(numBonds, "bondParams");
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
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    
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

//class CudaCustomBondForceInfo : public CudaForceInfo {
//public:
//    CudaCustomBondForceInfo(const CustomBondForce& force) : CudaForceInfo(0), force(force) {
//    }
//    int getNumParticleGroups() {
//        return force.getNumBonds();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int particle1, particle2;
//        vector<double> parameters;
//        force.getBondParameters(index, particle1, particle2, parameters);
//        particles.resize(2);
//        particles[0] = particle1;
//        particles[1] = particle2;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int particle1, particle2;
//        vector<double> parameters1, parameters2;
//        force.getBondParameters(group1, particle1, particle2, parameters1);
//        force.getBondParameters(group2, particle1, particle2, parameters2);
//        for (int i = 0; i < (int) parameters1.size(); i++)
//            if (parameters1[i] != parameters2[i])
//                return false;
//        return true;
//    }
//private:
//    const CustomBondForce& force;
//};
//
//CudaCalcCustomBondForceKernel::~CudaCalcCustomBondForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (globals != NULL)
//        delete globals;
//}
//
//void CudaCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
//    numBonds = endIndex-startIndex;
//    if (numBonds == 0)
//        return;
//    vector<vector<int> > atoms(numBonds, vector<int>(2));
//    params = new CudaParameterSet(cu, force.getNumPerBondParameters(), numBonds, "customBondParams");
//    vector<vector<cl_float> > paramVector(numBonds);
//    for (int i = 0; i < numBonds; i++) {
//        vector<double> parameters;
//        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    cu.addForce(new CudaCustomBondForceInfo(force));
//
//    // Record information for the expressions.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
//    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
//    map<string, Lepton::ParsedExpression> expressions;
//    expressions["energy += "] = energyExpression;
//    expressions["float dEdR = "] = forceExpression;
//
//    // Create the kernels.
//
//    map<string, string> variables;
//    variables["r"] = "r";
//    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
//        const string& name = force.getPerBondParameterName(i);
//        variables[name] = "bondParams"+params->getParameterSuffix(i);
//    }
//    if (force.getNumGlobalParameters() > 0) {
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customBondGlobals", false, CL_MEM_READ_ONLY);
//        globals->upload(globalParamValues);
//        string argName = cu.getBondedUtilities().addArgument(globals->getDeviceBuffer(), "float");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = argName+"["+cu.intToString(i)+"]";
//            variables[name] = value;
//        }
//    }
//    stringstream compute;
//    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
//        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
//    }
//    vector<pair<string, string> > functions;
//    compute << CudaExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
//    map<string, string> replacements;
//    replacements["COMPUTE_FORCE"] = compute.str();
//    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::bondForce, replacements), force.getForceGroup());
//}
//
//double CudaCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    return 0.0;
//}
//
//void CudaCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
//    if (numBonds != endIndex-startIndex)
//        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
//    
//    // Record the per-bond parameters.
//    
//    vector<vector<cl_float> > paramVector(numBonds);
//    vector<double> parameters;
//    for (int i = 0; i < numBonds; i++) {
//        int atom1, atom2;
//        force.getBondParameters(startIndex+i, atom1, atom2, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}

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
    if (params != NULL)
        delete params;
}

void CudaCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    numAngles = endIndex-startIndex;
    if (numAngles == 0)
        return;
    vector<vector<int> > atoms(numAngles, vector<int>(3));
    params = CudaArray::create<float2>(numAngles, "angleParams");
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
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
    if (numAngles != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of angles has changed");
    
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

//class CudaCustomAngleForceInfo : public CudaForceInfo {
//public:
//    CudaCustomAngleForceInfo(const CustomAngleForce& force) : CudaForceInfo(0), force(force) {
//    }
//    int getNumParticleGroups() {
//        return force.getNumAngles();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int particle1, particle2, particle3;
//        vector<double> parameters;
//        force.getAngleParameters(index, particle1, particle2, particle3, parameters);
//        particles.resize(3);
//        particles[0] = particle1;
//        particles[1] = particle2;
//        particles[2] = particle3;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int particle1, particle2, particle3;
//        vector<double> parameters1, parameters2;
//        force.getAngleParameters(group1, particle1, particle2, particle3, parameters1);
//        force.getAngleParameters(group2, particle1, particle2, particle3, parameters2);
//        for (int i = 0; i < (int) parameters1.size(); i++)
//            if (parameters1[i] != parameters2[i])
//                return false;
//        return true;
//    }
//private:
//    const CustomAngleForce& force;
//};
//
//CudaCalcCustomAngleForceKernel::~CudaCalcCustomAngleForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (globals != NULL)
//        delete globals;
//}
//
//void CudaCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
//    numAngles = endIndex-startIndex;
//    if (numAngles == 0)
//        return;
//    vector<vector<int> > atoms(numAngles, vector<int>(3));
//    params = new CudaParameterSet(cu, force.getNumPerAngleParameters(), numAngles, "customAngleParams");
//    vector<vector<cl_float> > paramVector(numAngles);
//    for (int i = 0; i < numAngles; i++) {
//        vector<double> parameters;
//        force.getAngleParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    cu.addForce(new CudaCustomAngleForceInfo(force));
//
//    // Record information for the expressions.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
//    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("theta").optimize();
//    map<string, Lepton::ParsedExpression> expressions;
//    expressions["energy += "] = energyExpression;
//    expressions["float dEdAngle = "] = forceExpression;
//
//    // Create the kernels.
//
//    map<string, string> variables;
//    variables["theta"] = "theta";
//    for (int i = 0; i < force.getNumPerAngleParameters(); i++) {
//        const string& name = force.getPerAngleParameterName(i);
//        variables[name] = "angleParams"+params->getParameterSuffix(i);
//    }
//    if (force.getNumGlobalParameters() > 0) {
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customAngleGlobals", false, CL_MEM_READ_ONLY);
//        globals->upload(globalParamValues);
//        string argName = cu.getBondedUtilities().addArgument(globals->getDeviceBuffer(), "float");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = argName+"["+cu.intToString(i)+"]";
//            variables[name] = value;
//        }
//    }
//    stringstream compute;
//    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
//        compute<<buffer.getType()<<" angleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
//    }
//    vector<pair<string, string> > functions;
//    compute << CudaExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
//    map<string, string> replacements;
//    replacements["COMPUTE_FORCE"] = compute.str();
//    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::angleForce, replacements), force.getForceGroup());
//}
//
//double CudaCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    return 0.0;
//}
//
//void CudaCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumAngles()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumAngles()/numContexts;
//    if (numAngles != endIndex-startIndex)
//        throw OpenMMException("updateParametersInContext: The number of angles has changed");
//    
//    // Record the per-angle parameters.
//    
//    vector<vector<cl_float> > paramVector(numAngles);
//    vector<double> parameters;
//    for (int i = 0; i < numAngles; i++) {
//        int atom1, atom2, atom3;
//        force.getAngleParameters(startIndex+i, atom1, atom2, atom3, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}

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
    if (params != NULL)
        delete params;
}

void CudaCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params = CudaArray::create<float4>(numTorsions, "periodicTorsionParams");
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
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    
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
    if (params1 != NULL)
        delete params1;
    if (params2 != NULL)
        delete params2;
}

void CudaCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    numTorsions = endIndex-startIndex;
    if (numTorsions == 0)
        return;
    vector<vector<int> > atoms(numTorsions, vector<int>(4));
    params1 = CudaArray::create<float4>(numTorsions, "rbTorsionParams1");
    params2 = CudaArray::create<float2>(numTorsions, "rbTorsionParams2");
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
    int numContexts = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
    if (numTorsions != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
    
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
    coefficients = CudaArray::create<float4>(coeffVec.size(), "cmapTorsionCoefficients");
    mapPositions = CudaArray::create<int2>(numMaps, "cmapTorsionMapPositions");
    torsionMaps = CudaArray::create<int>(numTorsions, "cmapTorsionMaps");
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

//class CudaCustomTorsionForceInfo : public CudaForceInfo {
//public:
//    CudaCustomTorsionForceInfo(const CustomTorsionForce& force) : CudaForceInfo(0), force(force) {
//    }
//    int getNumParticleGroups() {
//        return force.getNumTorsions();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int particle1, particle2, particle3, particle4;
//        vector<double> parameters;
//        force.getTorsionParameters(index, particle1, particle2, particle3, particle4, parameters);
//        particles.resize(4);
//        particles[0] = particle1;
//        particles[1] = particle2;
//        particles[2] = particle3;
//        particles[3] = particle4;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int particle1, particle2, particle3, particle4;
//        vector<double> parameters1, parameters2;
//        force.getTorsionParameters(group1, particle1, particle2, particle3, particle4, parameters1);
//        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, parameters2);
//        for (int i = 0; i < (int) parameters1.size(); i++)
//            if (parameters1[i] != parameters2[i])
//                return false;
//        return true;
//    }
//private:
//    const CustomTorsionForce& force;
//};
//
//CudaCalcCustomTorsionForceKernel::~CudaCalcCustomTorsionForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (globals != NULL)
//        delete globals;
//}
//
//void CudaCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
//    numTorsions = endIndex-startIndex;
//    if (numTorsions == 0)
//        return;
//    vector<vector<int> > atoms(numTorsions, vector<int>(4));
//    params = new CudaParameterSet(cu, force.getNumPerTorsionParameters(), numTorsions, "customTorsionParams");
//    vector<vector<cl_float> > paramVector(numTorsions);
//    for (int i = 0; i < numTorsions; i++) {
//        vector<double> parameters;
//        force.getTorsionParameters(startIndex+i, atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    cu.addForce(new CudaCustomTorsionForceInfo(force));
//
//    // Record information for the expressions.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
//    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("theta").optimize();
//    map<string, Lepton::ParsedExpression> expressions;
//    expressions["energy += "] = energyExpression;
//    expressions["float dEdAngle = "] = forceExpression;
//
//    // Create the kernels.
//
//    map<string, string> variables;
//    variables["theta"] = "theta";
//    for (int i = 0; i < force.getNumPerTorsionParameters(); i++) {
//        const string& name = force.getPerTorsionParameterName(i);
//        variables[name] = "torsionParams"+params->getParameterSuffix(i);
//    }
//    if (force.getNumGlobalParameters() > 0) {
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customTorsionGlobals", false, CL_MEM_READ_ONLY);
//        globals->upload(globalParamValues);
//        string argName = cu.getBondedUtilities().addArgument(globals->getDeviceBuffer(), "float");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = argName+"["+cu.intToString(i)+"]";
//            variables[name] = value;
//        }
//    }
//    stringstream compute;
//    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
//        compute<<buffer.getType()<<" torsionParams"<<(i+1)<<" = "<<argName<<"[index];\n";
//    }
//    vector<pair<string, string> > functions;
//    compute << CudaExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
//    map<string, string> replacements;
//    replacements["COMPUTE_FORCE"] = compute.str();
//    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::torsionForce, replacements), force.getForceGroup());
//}
//
//double CudaCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    return 0.0;
//}
//
//void CudaCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumTorsions()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumTorsions()/numContexts;
//    if (numTorsions != endIndex-startIndex)
//        throw OpenMMException("updateParametersInContext: The number of torsions has changed");
//    
//    // Record the per-torsion parameters.
//    
//    vector<vector<cl_float> > paramVector(numTorsions);
//    vector<double> parameters;
//    for (int i = 0; i < numTorsions; i++) {
//        int atom1, atom2, atom3, atom4;
//        force.getTorsionParameters(startIndex+i, atom1, atom2, atom3, atom4, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}
//
//class CudaNonbondedForceInfo : public CudaForceInfo {
//public:
//    CudaNonbondedForceInfo(int requiredBuffers, const NonbondedForce& force) : CudaForceInfo(requiredBuffers), force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        double charge1, charge2, sigma1, sigma2, epsilon1, epsilon2;
//        force.getParticleParameters(particle1, charge1, sigma1, epsilon1);
//        force.getParticleParameters(particle2, charge2, sigma2, epsilon2);
//        return (charge1 == charge2 && sigma1 == sigma2 && epsilon1 == epsilon2);
//    }
//    int getNumParticleGroups() {
//        return force.getNumExceptions();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int particle1, particle2;
//        double chargeProd, sigma, epsilon;
//        force.getExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
//        particles.resize(2);
//        particles[0] = particle1;
//        particles[1] = particle2;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int particle1, particle2;
//        double chargeProd1, chargeProd2, sigma1, sigma2, epsilon1, epsilon2;
//        force.getExceptionParameters(group1, particle1, particle2, chargeProd1, sigma1, epsilon1);
//        force.getExceptionParameters(group2, particle1, particle2, chargeProd2, sigma2, epsilon2);
//        return (chargeProd1 == chargeProd2 && sigma1 == sigma2 && epsilon1 == epsilon2);
//    }
//private:
//    const NonbondedForce& force;
//};
//
//CudaCalcNonbondedForceKernel::~CudaCalcNonbondedForceKernel() {
//    if (sigmaEpsilon != NULL)
//        delete sigmaEpsilon;
//    if (exceptionParams != NULL)
//        delete exceptionParams;
//    if (cosSinSums != NULL)
//        delete cosSinSums;
//    if (pmeGrid != NULL)
//        delete pmeGrid;
//    if (pmeGrid2 != NULL)
//        delete pmeGrid2;
//    if (pmeBsplineModuliX != NULL)
//        delete pmeBsplineModuliX;
//    if (pmeBsplineModuliY != NULL)
//        delete pmeBsplineModuliY;
//    if (pmeBsplineModuliZ != NULL)
//        delete pmeBsplineModuliZ;
//    if (pmeBsplineTheta != NULL)
//        delete pmeBsplineTheta;
//    if (pmeBsplineDTheta != NULL)
//        delete pmeBsplineDTheta;
//    if (pmeAtomRange != NULL)
//        delete pmeAtomRange;
//    if (pmeAtomGridIndex != NULL)
//        delete pmeAtomGridIndex;
//    if (sort != NULL)
//        delete sort;
//    if (fft != NULL)
//        delete fft;
//}
//
//void CudaCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
//
//    // Identify which exceptions are 1-4 interactions.
//
//    vector<pair<int, int> > exclusions;
//    vector<int> exceptions;
//    for (int i = 0; i < force.getNumExceptions(); i++) {
//        int particle1, particle2;
//        double chargeProd, sigma, epsilon;
//        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
//        exclusions.push_back(pair<int, int>(particle1, particle2));
//        if (chargeProd != 0.0 || epsilon != 0.0)
//            exceptions.push_back(i);
//    }
//
//    // Initialize nonbonded interactions.
//
//    int numParticles = force.getNumParticles();
//    sigmaEpsilon = new CudaArray<mm_float2>(cu, numParticles, "sigmaEpsilon");
//    CudaArray<mm_float4>& posq = cu.getPosq();
//    vector<mm_float2> sigmaEpsilonVector(numParticles);
//    vector<vector<int> > exclusionList(numParticles);
//    double sumSquaredCharges = 0.0;
//    hasCoulomb = false;
//    hasLJ = false;
//    for (int i = 0; i < numParticles; i++) {
//        double charge, sigma, epsilon;
//        force.getParticleParameters(i, charge, sigma, epsilon);
//        posq[i].w = (float) charge;
//        sigmaEpsilonVector[i] = mm_float2((float) (0.5*sigma), (float) (2.0*sqrt(epsilon)));
//        exclusionList[i].push_back(i);
//        sumSquaredCharges += charge*charge;
//        if (charge != 0.0)
//            hasCoulomb = true;
//        if (epsilon != 0.0)
//            hasLJ = true;
//    }
//    for (int i = 0; i < (int) exclusions.size(); i++) {
//        exclusionList[exclusions[i].first].push_back(exclusions[i].second);
//        exclusionList[exclusions[i].second].push_back(exclusions[i].first);
//    }
//    posq.upload();
//    sigmaEpsilon->upload(sigmaEpsilonVector);
//    bool useCutoff = (force.getNonbondedMethod() != NonbondedForce::NoCutoff);
//    bool usePeriodic = (force.getNonbondedMethod() != NonbondedForce::NoCutoff && force.getNonbondedMethod() != NonbondedForce::CutoffNonPeriodic);
//    map<string, string> defines;
//    defines["HAS_COULOMB"] = (hasCoulomb ? "1" : "0");
//    defines["HAS_LENNARD_JONES"] = (hasLJ ? "1" : "0");
//    if (useCutoff) {
//        // Compute the reaction field constants.
//
//        double reactionFieldK = pow(force.getCutoffDistance(), -3.0)*(force.getReactionFieldDielectric()-1.0)/(2.0*force.getReactionFieldDielectric()+1.0);
//        double reactionFieldC = (1.0 / force.getCutoffDistance())*(3.0*force.getReactionFieldDielectric())/(2.0*force.getReactionFieldDielectric()+1.0);
//        defines["REACTION_FIELD_K"] = cu.doubleToString(reactionFieldK);
//        defines["REACTION_FIELD_C"] = cu.doubleToString(reactionFieldC);
//    }
//    if (force.getUseDispersionCorrection() && cu.getContextIndex() == 0)
//        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
//    else
//        dispersionCoefficient = 0.0;
//    alpha = 0;
//    if (force.getNonbondedMethod() == NonbondedForce::Ewald) {
//        // Compute the Ewald parameters.
//
//        int kmaxx, kmaxy, kmaxz;
//        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
//        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
//        defines["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
//        defines["USE_EWALD"] = "1";
//        ewaldSelfEnergy = (cu.getContextIndex() == 0 ? -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI) : 0.0);
//
//        // Create the reciprocal space kernels.
//
//        map<string, string> replacements;
//        replacements["NUM_ATOMS"] = cu.intToString(numParticles);
//        replacements["KMAX_X"] = cu.intToString(kmaxx);
//        replacements["KMAX_Y"] = cu.intToString(kmaxy);
//        replacements["KMAX_Z"] = cu.intToString(kmaxz);
//        replacements["EXP_COEFFICIENT"] = cu.doubleToString(-1.0/(4.0*alpha*alpha));
//        cu::Program program = cu.createProgram(CudaKernelSources::ewald, replacements);
//        ewaldSumsKernel = cu::Kernel(program, "calculateEwaldCosSinSums");
//        ewaldForcesKernel = cu::Kernel(program, "calculateEwaldForces");
//        cosSinSums = new CudaArray<mm_float2>(cu, (2*kmaxx-1)*(2*kmaxy-1)*(2*kmaxz-1), "cosSinSums");
//    }
//    else if (force.getNonbondedMethod() == NonbondedForce::PME) {
//        // Compute the PME parameters.
//
//        int gridSizeX, gridSizeY, gridSizeZ;
//        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ);
//        gridSizeX = CudaFFT3D::findLegalDimension(gridSizeX);
//        gridSizeY = CudaFFT3D::findLegalDimension(gridSizeY);
//        gridSizeZ = CudaFFT3D::findLegalDimension(gridSizeZ);
//        defines["EWALD_ALPHA"] = cu.doubleToString(alpha);
//        defines["TWO_OVER_SQRT_PI"] = cu.doubleToString(2.0/sqrt(M_PI));
//        defines["USE_EWALD"] = "1";
//        ewaldSelfEnergy = (cu.getContextIndex() == 0 ? -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI) : 0.0);
//        pmeDefines["PME_ORDER"] = cu.intToString(PmeOrder);
//        pmeDefines["NUM_ATOMS"] = cu.intToString(numParticles);
//        pmeDefines["RECIP_EXP_FACTOR"] = cu.doubleToString(M_PI*M_PI/(alpha*alpha));
//        pmeDefines["GRID_SIZE_X"] = cu.intToString(gridSizeX);
//        pmeDefines["GRID_SIZE_Y"] = cu.intToString(gridSizeY);
//        pmeDefines["GRID_SIZE_Z"] = cu.intToString(gridSizeZ);
//        pmeDefines["EPSILON_FACTOR"] = cu.doubleToString(sqrt(ONE_4PI_EPS0));
//
//        // Create required data structures.
//
//        pmeGrid = new CudaArray<mm_float2>(cu, gridSizeX*gridSizeY*gridSizeZ, "pmeGrid");
//        cu.addAutoclearBuffer(pmeGrid->getDeviceBuffer(), pmeGrid->getSize()*2);
//        pmeGrid2 = new CudaArray<mm_float2>(cu, gridSizeX*gridSizeY*gridSizeZ, "pmeGrid2");
//        pmeBsplineModuliX = new CudaArray<cl_float>(cu, gridSizeX, "pmeBsplineModuliX");
//        pmeBsplineModuliY = new CudaArray<cl_float>(cu, gridSizeY, "pmeBsplineModuliY");
//        pmeBsplineModuliZ = new CudaArray<cl_float>(cu, gridSizeZ, "pmeBsplineModuliZ");
//        pmeBsplineTheta = new CudaArray<mm_float4>(cu, PmeOrder*numParticles, "pmeBsplineTheta");
//        bool deviceIsCpu = (cu.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
//        if (deviceIsCpu)
//            pmeBsplineDTheta = new CudaArray<mm_float4>(cu, PmeOrder*numParticles, "pmeBsplineDTheta");
//        pmeAtomRange = new CudaArray<cl_int>(cu, gridSizeX*gridSizeY*gridSizeZ+1, "pmeAtomRange");
//        pmeAtomGridIndex = new CudaArray<mm_int2>(cu, numParticles, "pmeAtomGridIndex");
//        sort = new CudaSort<SortTrait>(cu, cu.getNumAtoms());
//        fft = new CudaFFT3D(cu, gridSizeX, gridSizeY, gridSizeZ);
//
//        // Initialize the b-spline moduli.
//
//        int maxSize = max(max(gridSizeX, gridSizeY), gridSizeZ);
//        vector<double> data(PmeOrder);
//        vector<double> ddata(PmeOrder);
//        vector<double> bsplines_data(maxSize);
//        data[PmeOrder-1] = 0.0;
//        data[1] = 0.0;
//        data[0] = 1.0;
//        for (int i = 3; i < PmeOrder; i++) {
//            double div = 1.0/(i-1.0);
//            data[i-1] = 0.0;
//            for (int j = 1; j < (i-1); j++)
//                data[i-j-1] = div*(j*data[i-j-2]+(i-j)*data[i-j-1]);
//            data[0] = div*data[0];
//        }
//
//        // Differentiate.
//
//        ddata[0] = -data[0];
//        for (int i = 1; i < PmeOrder; i++)
//            ddata[i] = data[i-1]-data[i];
//        double div = 1.0/(PmeOrder-1);
//        data[PmeOrder-1] = 0.0;
//        for (int i = 1; i < (PmeOrder-1); i++)
//            data[PmeOrder-i-1] = div*(i*data[PmeOrder-i-2]+(PmeOrder-i)*data[PmeOrder-i-1]);
//        data[0] = div*data[0];
//        for (int i = 0; i < maxSize; i++)
//            bsplines_data[i] = 0.0;
//        for (int i = 1; i <= PmeOrder; i++)
//            bsplines_data[i] = data[i-1];
//
//        // Evaluate the actual bspline moduli for X/Y/Z.
//
//        for(int dim = 0; dim < 3; dim++) {
//            int ndata = (dim == 0 ? gridSizeX : dim == 1 ? gridSizeY : gridSizeZ);
//            vector<cl_float> moduli(ndata);
//            for (int i = 0; i < ndata; i++) {
//                double sc = 0.0;
//                double ss = 0.0;
//                for (int j = 0; j < ndata; j++) {
//                    double arg = (2.0*M_PI*i*j)/ndata;
//                    sc += bsplines_data[j]*cos(arg);
//                    ss += bsplines_data[j]*sin(arg);
//                }
//                moduli[i] = (float) (sc*sc+ss*ss);
//            }
//            for (int i = 0; i < ndata; i++)
//            {
//                if (moduli[i] < 1.0e-7)
//                    moduli[i] = (moduli[i-1]+moduli[i+1])*0.5f;
//            }
//            if (dim == 0)
//                pmeBsplineModuliX->upload(moduli);
//            else if (dim == 1)
//                pmeBsplineModuliY->upload(moduli);
//            else
//                pmeBsplineModuliZ->upload(moduli);
//        }
//    }
//    else
//        ewaldSelfEnergy = 0.0;
//
//    // Add the interaction to the default nonbonded kernel.
//    
//    string source = cu.replaceStrings(CudaKernelSources::coulombLennardJones, defines);
//    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
//    if (hasLJ)
//        cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo("sigmaEpsilon", "float", 2, sizeof(cl_float2), sigmaEpsilon->getDeviceBuffer()));
//
//    // Initialize the exceptions.
//
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*exceptions.size()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*exceptions.size()/numContexts;
//    int numExceptions = endIndex-startIndex;
//    if (numExceptions > 0) {
//        vector<vector<int> > atoms(numExceptions, vector<int>(2));
//        exceptionParams = new CudaArray<mm_float4>(cu, numExceptions, "exceptionParams");
//        vector<mm_float4> exceptionParamsVector(numExceptions);
//        for (int i = 0; i < numExceptions; i++) {
//            double chargeProd, sigma, epsilon;
//            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
//            exceptionParamsVector[i] = mm_float4((float) (ONE_4PI_EPS0*chargeProd), (float) sigma, (float) (4.0*epsilon), 0.0f);
//        }
//        exceptionParams->upload(exceptionParamsVector);
//        map<string, string> replacements;
//        replacements["PARAMS"] = cu.getBondedUtilities().addArgument(exceptionParams->getDeviceBuffer(), "float4");
//        cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::nonbondedExceptions, replacements), force.getForceGroup());
//    }
//    cu.addForce(new CudaNonbondedForceInfo(cu.getNonbondedUtilities().getNumForceBuffers(), force));
//}
//
//double CudaCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
//    bool deviceIsCpu = (cu.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
//    if (!hasInitializedKernel) {
//        hasInitializedKernel = true;
//        if (cosSinSums != NULL) {
//            ewaldSumsKernel.setArg<cu::Buffer>(0, cu.getEnergyBuffer().getDeviceBuffer());
//            ewaldSumsKernel.setArg<cu::Buffer>(1, cu.getPosq().getDeviceBuffer());
//            ewaldSumsKernel.setArg<cu::Buffer>(2, cosSinSums->getDeviceBuffer());
//            ewaldForcesKernel.setArg<cu::Buffer>(0, cu.getForceBuffers().getDeviceBuffer());
//            ewaldForcesKernel.setArg<cu::Buffer>(1, cu.getPosq().getDeviceBuffer());
//            ewaldForcesKernel.setArg<cu::Buffer>(2, cosSinSums->getDeviceBuffer());
//        }
//        if (pmeGrid != NULL) {
//            string file = (deviceIsCpu ? CudaKernelSources::pme_cpu : CudaKernelSources::pme);
//            cu::Program program = cu.createProgram(file, pmeDefines);
//            pmeUpdateBsplinesKernel = cu::Kernel(program, "updateBsplines");
//            pmeAtomRangeKernel = cu::Kernel(program, "findAtomRangeForGrid");
//	    if (!deviceIsCpu)
//                pmeZIndexKernel = cu::Kernel(program, "recordZIndex");
//            pmeSpreadChargeKernel = cu::Kernel(program, "gridSpreadCharge");
//            pmeConvolutionKernel = cu::Kernel(program, "reciprocalConvolution");
//            pmeInterpolateForceKernel = cu::Kernel(program, "gridInterpolateForce");
//            pmeUpdateBsplinesKernel.setArg<cu::Buffer>(0, cu.getPosq().getDeviceBuffer());
//            pmeUpdateBsplinesKernel.setArg<cu::Buffer>(1, pmeBsplineTheta->getDeviceBuffer());
//            pmeUpdateBsplinesKernel.setArg(2, CudaContext::ThreadBlockSize*PmeOrder*sizeof(mm_float4), NULL);
//            pmeUpdateBsplinesKernel.setArg<cu::Buffer>(3, pmeAtomGridIndex->getDeviceBuffer());
//            if (deviceIsCpu)
//                pmeUpdateBsplinesKernel.setArg<cu::Buffer>(6, pmeBsplineDTheta->getDeviceBuffer());
//            pmeAtomRangeKernel.setArg<cu::Buffer>(0, pmeAtomGridIndex->getDeviceBuffer());
//            pmeAtomRangeKernel.setArg<cu::Buffer>(1, pmeAtomRange->getDeviceBuffer());
//            pmeAtomRangeKernel.setArg<cu::Buffer>(2, cu.getPosq().getDeviceBuffer());
//	    if (!deviceIsCpu) {
//                pmeZIndexKernel.setArg<cu::Buffer>(0, pmeAtomGridIndex->getDeviceBuffer());
//                pmeZIndexKernel.setArg<cu::Buffer>(1, cu.getPosq().getDeviceBuffer());
//	    }
//            pmeSpreadChargeKernel.setArg<cu::Buffer>(0, cu.getPosq().getDeviceBuffer());
//            pmeSpreadChargeKernel.setArg<cu::Buffer>(1, pmeAtomGridIndex->getDeviceBuffer());
//            pmeSpreadChargeKernel.setArg<cu::Buffer>(2, pmeAtomRange->getDeviceBuffer());
//            pmeSpreadChargeKernel.setArg<cu::Buffer>(3, pmeGrid->getDeviceBuffer());
//            pmeSpreadChargeKernel.setArg<cu::Buffer>(4, pmeBsplineTheta->getDeviceBuffer());
//            pmeConvolutionKernel.setArg<cu::Buffer>(0, pmeGrid2->getDeviceBuffer());
//            pmeConvolutionKernel.setArg<cu::Buffer>(1, cu.getEnergyBuffer().getDeviceBuffer());
//            pmeConvolutionKernel.setArg<cu::Buffer>(2, pmeBsplineModuliX->getDeviceBuffer());
//            pmeConvolutionKernel.setArg<cu::Buffer>(3, pmeBsplineModuliY->getDeviceBuffer());
//            pmeConvolutionKernel.setArg<cu::Buffer>(4, pmeBsplineModuliZ->getDeviceBuffer());
//            interpolateForceThreads = (cu.getDevice().getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() > 2*128*PmeOrder*sizeof(mm_float4) ? 128 : 64);
//            pmeInterpolateForceKernel.setArg<cu::Buffer>(0, cu.getPosq().getDeviceBuffer());
//            pmeInterpolateForceKernel.setArg<cu::Buffer>(1, cu.getForceBuffers().getDeviceBuffer());
//            pmeInterpolateForceKernel.setArg<cu::Buffer>(2, pmeGrid->getDeviceBuffer());
//            if (deviceIsCpu) {
//                pmeInterpolateForceKernel.setArg<cu::Buffer>(5, pmeBsplineTheta->getDeviceBuffer());
//                pmeInterpolateForceKernel.setArg<cu::Buffer>(6, pmeBsplineDTheta->getDeviceBuffer());
//            }
//            else
//                pmeInterpolateForceKernel.setArg(5, 2*interpolateForceThreads*PmeOrder*sizeof(mm_float4), NULL);
//            if (cu.getSupports64BitGlobalAtomics()) {
//                pmeFinishSpreadChargeKernel = cu::Kernel(program, "finishSpreadCharge");
//                pmeFinishSpreadChargeKernel.setArg<cu::Buffer>(0, pmeGrid->getDeviceBuffer());
//            }
//       }
//    }
//    if (cosSinSums != NULL && cu.getContextIndex() == 0 && includeReciprocal) {
//        mm_float4 boxSize = cu.getPeriodicBoxSize();
//        mm_float4 recipBoxSize = mm_float4((float) (2*M_PI/boxSize.x), (float) (2*M_PI/boxSize.y), (float) (2*M_PI/boxSize.z), 0);
//        float recipCoefficient = (float) (ONE_4PI_EPS0*4*M_PI/(boxSize.x*boxSize.y*boxSize.z));
//        ewaldSumsKernel.setArg<mm_float4>(3, recipBoxSize);
//        ewaldSumsKernel.setArg<cl_float>(4, recipCoefficient);
//        cu.executeKernel(ewaldSumsKernel, cosSinSums->getSize());
//        ewaldForcesKernel.setArg<mm_float4>(3, recipBoxSize);
//        ewaldForcesKernel.setArg<cl_float>(4, recipCoefficient);
//        cu.executeKernel(ewaldForcesKernel, cu.getNumAtoms());
//    }
//    if (pmeGrid != NULL && cu.getContextIndex() == 0 && includeReciprocal) {
//        mm_float4 boxSize = cu.getPeriodicBoxSize();
//        mm_float4 invBoxSize = cu.getInvPeriodicBoxSize();
//        pmeUpdateBsplinesKernel.setArg<mm_float4>(4, boxSize);
//        pmeUpdateBsplinesKernel.setArg<mm_float4>(5, invBoxSize);
//        cu.executeKernel(pmeUpdateBsplinesKernel, cu.getNumAtoms());
//        if (deviceIsCpu) {
//            pmeSpreadChargeKernel.setArg<mm_float4>(5, boxSize);
//            pmeSpreadChargeKernel.setArg<mm_float4>(6, invBoxSize);
//            cu.executeKernel(pmeSpreadChargeKernel, 2*cu.getDevice().getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>(), 1);
//        }
//        else {
//            sort->sort(*pmeAtomGridIndex);
//            pmeAtomRangeKernel.setArg<mm_float4>(3, boxSize);
//            pmeAtomRangeKernel.setArg<mm_float4>(4, invBoxSize);
//            cu.executeKernel(pmeAtomRangeKernel, cu.getNumAtoms());
//            if (cu.getSupports64BitGlobalAtomics()) {
//                pmeSpreadChargeKernel.setArg<mm_float4>(5, boxSize);
//                pmeSpreadChargeKernel.setArg<mm_float4>(6, invBoxSize);
//                cu.executeKernel(pmeSpreadChargeKernel, cu.getNumAtoms(), PmeOrder*PmeOrder*PmeOrder);
//                cu.executeKernel(pmeFinishSpreadChargeKernel, pmeGrid->getSize());
//            }
//            else {
//                pmeZIndexKernel.setArg<mm_float4>(2, boxSize);
//                pmeZIndexKernel.setArg<mm_float4>(3, invBoxSize);
//                cu.executeKernel(pmeZIndexKernel, cu.getNumAtoms());
//                cu.executeKernel(pmeSpreadChargeKernel, cu.getNumAtoms());
//            }
//        }
//        fft->execFFT(*pmeGrid, *pmeGrid2, true);
//        pmeConvolutionKernel.setArg<mm_float4>(5, invBoxSize);
//        pmeConvolutionKernel.setArg<cl_float>(6, (float) (1.0/(M_PI*boxSize.x*boxSize.y*boxSize.z)));
//        cu.executeKernel(pmeConvolutionKernel, cu.getNumAtoms());
//        fft->execFFT(*pmeGrid2, *pmeGrid, false);
//        pmeInterpolateForceKernel.setArg<mm_float4>(3, boxSize);
//        pmeInterpolateForceKernel.setArg<mm_float4>(4, invBoxSize);
//        cu.executeKernel(pmeInterpolateForceKernel, cu.getNumAtoms(), interpolateForceThreads);
//    }
//    double energy = (includeReciprocal ? ewaldSelfEnergy : 0.0);
//    if (dispersionCoefficient != 0.0 && includeDirect) {
//        mm_float4 boxSize = cu.getPeriodicBoxSize();
//        energy += dispersionCoefficient/(boxSize.x*boxSize.y*boxSize.z);
//    }
//    return energy;
//}
//
//void CudaCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
//    // Make sure the new parameters are acceptable.
//    
//    if (force.getNumParticles() != cu.getNumAtoms())
//        throw OpenMMException("updateParametersInContext: The number of particles has changed");
//    if (!hasCoulomb || !hasLJ) {
//        for (int i = 0; i < force.getNumParticles(); i++) {
//            double charge, sigma, epsilon;
//            force.getParticleParameters(i, charge, sigma, epsilon);
//            if (!hasCoulomb && charge != 0.0)
//                throw OpenMMException("updateParametersInContext: The nonbonded force kernel does not include Coulomb interactions, because all charges were originally 0");
//            if (!hasLJ && epsilon != 0.0)
//                throw OpenMMException("updateParametersInContext: The nonbonded force kernel does not include Lennard-Jones interactions, because all epsilons were originally 0");
//        }
//    }
//    vector<int> exceptions;
//    for (int i = 0; i < force.getNumExceptions(); i++) {
//        int particle1, particle2;
//        double chargeProd, sigma, epsilon;
//        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
//        if (chargeProd != 0.0 || epsilon != 0.0)
//            exceptions.push_back(i);
//    }
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*exceptions.size()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*exceptions.size()/numContexts;
//    int numExceptions = endIndex-startIndex;
//    if ((exceptionParams == NULL && numExceptions > 0) || (exceptionParams != NULL && numExceptions != exceptionParams->getSize()))
//        throw OpenMMException("updateParametersInContext: The number of non-excluded exceptions has changed");
//    
//    // Record the per-particle parameters.
//    
//    CudaArray<mm_float4>& posq = cu.getPosq();
//    posq.download();
//    vector<mm_float2> sigmaEpsilonVector(force.getNumParticles());
//    double sumSquaredCharges = 0.0;
//    CudaArray<cl_int>& order = cu.getAtomIndex();
//    for (int i = 0; i < force.getNumParticles(); i++) {
//        int index = order[i];
//        double charge, sigma, epsilon;
//        force.getParticleParameters(index, charge, sigma, epsilon);
//        posq[i].w = (float) charge;
//        sigmaEpsilonVector[index] = mm_float2((float) (0.5*sigma), (float) (2.0*sqrt(epsilon)));
//        sumSquaredCharges += charge*charge;
//    }
//    posq.upload();
//    sigmaEpsilon->upload(sigmaEpsilonVector);
//    
//    // Record the exceptions.
//    
//    if (numExceptions > 0) {
//        vector<vector<int> > atoms(numExceptions, vector<int>(2));
//        vector<mm_float4> exceptionParamsVector(numExceptions);
//        for (int i = 0; i < numExceptions; i++) {
//            double chargeProd, sigma, epsilon;
//            force.getExceptionParameters(exceptions[startIndex+i], atoms[i][0], atoms[i][1], chargeProd, sigma, epsilon);
//            exceptionParamsVector[i] = mm_float4((float) (ONE_4PI_EPS0*chargeProd), (float) sigma, (float) (4.0*epsilon), 0.0f);
//        }
//        exceptionParams->upload(exceptionParamsVector);
//    }
//    
//    // Compute other values.
//    
//    NonbondedForce::NonbondedMethod method = force.getNonbondedMethod();
//    if (method == NonbondedForce::Ewald || method == NonbondedForce::PME)
//        ewaldSelfEnergy = (cu.getContextIndex() == 0 ? -ONE_4PI_EPS0*alpha*sumSquaredCharges/sqrt(M_PI) : 0.0);
//    if (force.getUseDispersionCorrection() && cu.getContextIndex() == 0 && (method == NonbondedForce::CutoffPeriodic || method == NonbondedForce::Ewald || method == NonbondedForce::PME))
//        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
//    cu.invalidateMolecules();
//}
//
//class CudaCustomNonbondedForceInfo : public CudaForceInfo {
//public:
//    CudaCustomNonbondedForceInfo(int requiredBuffers, const CustomNonbondedForce& force) : CudaForceInfo(requiredBuffers), force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        vector<double> params1;
//        vector<double> params2;
//        force.getParticleParameters(particle1, params1);
//        force.getParticleParameters(particle2, params2);
//        for (int i = 0; i < (int) params1.size(); i++)
//            if (params1[i] != params2[i])
//                return false;
//        return true;
//    }
//    int getNumParticleGroups() {
//        return force.getNumExclusions();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int particle1, particle2;
//        force.getExclusionParticles(index, particle1, particle2);
//        particles.resize(2);
//        particles[0] = particle1;
//        particles[1] = particle2;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        return true;
//    }
//private:
//    const CustomNonbondedForce& force;
//};
//
//CudaCalcCustomNonbondedForceKernel::~CudaCalcCustomNonbondedForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (globals != NULL)
//        delete globals;
//    if (tabulatedFunctionParams != NULL)
//        delete tabulatedFunctionParams;
//    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//        delete tabulatedFunctions[i];
//}
//
//void CudaCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
//    int forceIndex;
//    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
//        ;
//    string prefix = "custom"+cu.intToString(forceIndex)+"_";
//
//    // Record parameters and exclusions.
//
//    int numParticles = force.getNumParticles();
//    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customNonbondedParameters");
//    if (force.getNumGlobalParameters() > 0)
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customNonbondedGlobals", false, CL_MEM_READ_ONLY);
//    vector<vector<cl_float> > paramVector(numParticles);
//    vector<vector<int> > exclusionList(numParticles);
//    for (int i = 0; i < numParticles; i++) {
//        vector<double> parameters;
//        force.getParticleParameters(i, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//        exclusionList[i].push_back(i);
//    }
//    for (int i = 0; i < force.getNumExclusions(); i++) {
//        int particle1, particle2;
//        force.getExclusionParticles(i, particle1, particle2);
//        exclusionList[particle1].push_back(particle2);
//        exclusionList[particle2].push_back(particle1);
//    }
//    params->setParameterValues(paramVector);
//
//    // Record the tabulated functions.
//
//    CudaExpressionUtilities::FunctionPlaceholder fp;
//    map<string, Lepton::CustomFunction*> functions;
//    vector<pair<string, string> > functionDefinitions;
//    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
//    for (int i = 0; i < force.getNumFunctions(); i++) {
//        string name;
//        vector<double> values;
//        double min, max;
//        force.getFunctionParameters(i, name, values, min, max);
//        string arrayName = prefix+"table"+cu.intToString(i);
//        functionDefinitions.push_back(make_pair(name, arrayName));
//        functions[name] = &fp;
//        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), (float) values.size()-2);
//        vector<mm_float4> f = CudaExpressionUtilities::computeFunctionCoefficients(values, min, max);
//        tabulatedFunctions.push_back(new CudaArray<mm_float4>(cu, values.size()-1, "TabulatedFunction"));
//        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
//        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName, "float", 4, sizeof(cl_float4), tabulatedFunctions[tabulatedFunctions.size()-1]->getDeviceBuffer()));
//    }
//    if (force.getNumFunctions() > 0) {
//        tabulatedFunctionParams = new CudaArray<mm_float4>(cu, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
//        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
//        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(prefix+"functionParams", "float", 4, sizeof(cl_float4), tabulatedFunctionParams->getDeviceBuffer()));
//    }
//
//    // Record information for the expressions.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    if (globals != NULL)
//        globals->upload(globalParamValues);
//    bool useCutoff = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff);
//    bool usePeriodic = (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff && force.getNonbondedMethod() != CustomNonbondedForce::CutoffNonPeriodic);
//    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
//    Lepton::ParsedExpression forceExpression = energyExpression.differentiate("r").optimize();
//    map<string, Lepton::ParsedExpression> forceExpressions;
//    forceExpressions["tempEnergy += "] = energyExpression;
//    forceExpressions["tempForce -= "] = forceExpression;
//
//    // Create the kernels.
//
//    vector<pair<ExpressionTreeNode, string> > variables;
//    ExpressionTreeNode rnode(new Operation::Variable("r"));
//    variables.push_back(make_pair(rnode, "r"));
//    variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
//    variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
//    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
//        const string& name = force.getPerParticleParameterName(i);
//        variables.push_back(makeVariable(name+"1", prefix+"params"+params->getParameterSuffix(i, "1")));
//        variables.push_back(makeVariable(name+"2", prefix+"params"+params->getParameterSuffix(i, "2")));
//    }
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        const string& name = force.getGlobalParameterName(i);
//        string value = "globals["+cu.intToString(i)+"]";
//        variables.push_back(makeVariable(name, prefix+value));
//    }
//    stringstream compute;
//    compute << CudaExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, prefix+"temp", prefix+"functionParams");
//    map<string, string> replacements;
//    replacements["COMPUTE_FORCE"] = compute.str();
//    string source = cu.replaceStrings(CudaKernelSources::customNonbonded, replacements);
//    cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, true, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
//    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//        cu.getNonbondedUtilities().addParameter(CudaNonbondedUtilities::ParameterInfo(prefix+"params"+cu.intToString(i+1), buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
//    }
//    if (globals != NULL) {
//        globals->upload(globalParamValues);
//        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(cl_float), globals->getDeviceBuffer()));
//    }
//    cu.addForce(new CudaCustomNonbondedForceInfo(cu.getNonbondedUtilities().getNumForceBuffers(), force));
//}
//
//double CudaCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    return 0.0;
//}
//
//void CudaCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
//    int numParticles = force.getNumParticles();
//    if (numParticles != cu.getNumAtoms())
//        throw OpenMMException("updateParametersInContext: The number of particles has changed");
//    
//    // Record the per-particle parameters.
//    
//    vector<vector<cl_float> > paramVector(numParticles);
//    vector<double> parameters;
//    for (int i = 0; i < numParticles; i++) {
//        force.getParticleParameters(i, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}
//
//class CudaGBSAOBCForceInfo : public CudaForceInfo {
//public:
//    CudaGBSAOBCForceInfo(int requiredBuffers, const GBSAOBCForce& force) : CudaForceInfo(requiredBuffers), force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        double charge1, charge2, radius1, radius2, scale1, scale2;
//        force.getParticleParameters(particle1, charge1, radius1, scale1);
//        force.getParticleParameters(particle2, charge2, radius2, scale2);
//        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2);
//    }
//private:
//    const GBSAOBCForce& force;
//};
//
//CudaCalcGBSAOBCForceKernel::~CudaCalcGBSAOBCForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (bornSum != NULL)
//        delete bornSum;
//    if (longBornSum != NULL)
//        delete longBornSum;
//    if (bornRadii != NULL)
//        delete bornRadii;
//    if (bornForce != NULL)
//        delete bornForce;
//    if (longBornForce != NULL)
//        delete longBornForce;
//    if (obcChain != NULL)
//        delete obcChain;
//}
//
//void CudaCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
//    if (cu.getPlatformData().contexts.size() > 1)
//        throw OpenMMException("GBSAOBCForce does not support using multiple CUDA devices");
//    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
//    params = new CudaArray<mm_float2>(cu, cu.getPaddedNumAtoms(), "gbsaObcParams");
//    bornRadii = new CudaArray<cl_float>(cu, cu.getPaddedNumAtoms(), "bornRadii");
//    obcChain = new CudaArray<cl_float>(cu, cu.getPaddedNumAtoms(), "obcChain");
//    if (cu.getSupports64BitGlobalAtomics()) {
//        longBornSum = new CudaArray<cl_long>(cu, cu.getPaddedNumAtoms(), "longBornSum");
//        longBornForce = new CudaArray<cl_long>(cu, cu.getPaddedNumAtoms(), "longBornForce");
//        bornForce = new CudaArray<cl_float>(cu, cu.getPaddedNumAtoms(), "bornForce");
//        cu.addAutoclearBuffer(longBornSum->getDeviceBuffer(), 2*longBornSum->getSize());
//        cu.addAutoclearBuffer(longBornForce->getDeviceBuffer(), 2*longBornForce->getSize());
//    }
//    else {
//        bornSum = new CudaArray<cl_float>(cu, cu.getPaddedNumAtoms()*nb.getNumForceBuffers(), "bornSum");
//        bornForce = new CudaArray<cl_float>(cu, cu.getPaddedNumAtoms()*nb.getNumForceBuffers(), "bornForce");
//        cu.addAutoclearBuffer(bornSum->getDeviceBuffer(), bornSum->getSize());
//        cu.addAutoclearBuffer(bornForce->getDeviceBuffer(), bornForce->getSize());
//    }
//    CudaArray<mm_float4>& posq = cu.getPosq();
//    int numParticles = force.getNumParticles();
//    vector<mm_float2> paramsVector(numParticles);
//    const double dielectricOffset = 0.009;
//    for (int i = 0; i < numParticles; i++) {
//        double charge, radius, scalingFactor;
//        force.getParticleParameters(i, charge, radius, scalingFactor);
//        radius -= dielectricOffset;
//        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
//        posq[i].w = (float) charge;
//    }
//    posq.upload();
//    params->upload(paramsVector);
//    prefactor = -ONE_4PI_EPS0*((1.0/force.getSoluteDielectric())-(1.0/force.getSolventDielectric()));
//    bool useCutoff = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff);
//    bool usePeriodic = (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff && force.getNonbondedMethod() != GBSAOBCForce::CutoffNonPeriodic);
//    string source = CudaKernelSources::gbsaObc2;
//    nb.addInteraction(useCutoff, usePeriodic, false, force.getCutoffDistance(), vector<vector<int> >(), source, force.getForceGroup());
//    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("obcParams", "float", 2, sizeof(cl_float2), params->getDeviceBuffer()));;
//    nb.addParameter(CudaNonbondedUtilities::ParameterInfo("bornForce", "float", 1, sizeof(cl_float), bornForce->getDeviceBuffer()));;
//    cu.addForce(new CudaGBSAOBCForceInfo(nb.getNumForceBuffers(), force));
//}
//
//double CudaCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
//    bool deviceIsCpu = (cu.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
//    if (!hasCreatedKernels) {
//        // These Kernels cannot be created in initialize(), because the CudaNonbondedUtilities has not been initialized yet then.
//
//        hasCreatedKernels = true;
//        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : 0);
//        map<string, string> defines;
//        if (nb.getForceBufferPerAtomBlock())
//            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
//        if (nb.getUseCutoff())
//            defines["USE_CUTOFF"] = "1";
//        if (nb.getUsePeriodic())
//            defines["USE_PERIODIC"] = "1";
//        defines["CUTOFF_SQUARED"] = cu.doubleToString(nb.getCutoffDistance()*nb.getCutoffDistance());
//        defines["PREFACTOR"] = cu.doubleToString(prefactor);
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//        defines["NUM_BLOCKS"] = CudaExpressionUtilities::cu.intToString(cu.getNumAtomBlocks());
//        defines["FORCE_WORK_GROUP_SIZE"] = CudaExpressionUtilities::cu.intToString(nb.getForceThreadBlockSize());
//        string platformVendor = cu::Platform(cu.getDevice().getInfo<CL_DEVICE_PLATFORM>()).getInfo<CL_PLATFORM_VENDOR>();
//        if (platformVendor == "Apple")
//            defines["USE_APPLE_WORKAROUND"] = "1";
//        string file;
//        if (deviceIsCpu)
//            file = CudaKernelSources::gbsaObc_cpu;
//        else if (cu.getSIMDWidth() == 32)
//            file = CudaKernelSources::gbsaObc_nvidia;
//        else
//            file = CudaKernelSources::gbsaObc_default;
//        cu::Program program = cu.createProgram(file, defines);
//        bool useLong = (cu.getSupports64BitGlobalAtomics() && !deviceIsCpu);
//        int index = 0;
//        computeBornSumKernel = cu::Kernel(program, "computeBornSum");
//        computeBornSumKernel.setArg<cu::Buffer>(index++, (useLong ? longBornSum->getDeviceBuffer() : bornSum->getDeviceBuffer()));
//        computeBornSumKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        computeBornSumKernel.setArg<cu::Buffer>(index++, params->getDeviceBuffer());
//        if (nb.getUseCutoff()) {
//            computeBornSumKernel.setArg<cu::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
//            computeBornSumKernel.setArg<cu::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
//            index += 2; // The periodic box size arguments are set when the kernel is executed.
//            computeBornSumKernel.setArg<cl_uint>(index++, maxTiles);
//            if (cu.getSIMDWidth() == 32 || deviceIsCpu)
//                computeBornSumKernel.setArg<cu::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
//        }
//        else
//            computeBornSumKernel.setArg<cl_uint>(index++, cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
//        if (cu.getSIMDWidth() == 32) {
//            computeBornSumKernel.setArg<cu::Buffer>(index++, nb.getExclusionIndices().getDeviceBuffer());
//            computeBornSumKernel.setArg<cu::Buffer>(index++, nb.getExclusionRowIndices().getDeviceBuffer());
//        }
//        force1Kernel = cu::Kernel(program, "computeGBSAForce1");
//        index = 0;
//        force1Kernel.setArg<cu::Buffer>(index++, (useLong ? cu.getLongForceBuffer().getDeviceBuffer() : cu.getForceBuffers().getDeviceBuffer()));
//        force1Kernel.setArg<cu::Buffer>(index++, (useLong ? longBornForce->getDeviceBuffer() : bornForce->getDeviceBuffer()));
//        force1Kernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        force1Kernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        force1Kernel.setArg<cu::Buffer>(index++, bornRadii->getDeviceBuffer());
//        if (nb.getUseCutoff()) {
//            force1Kernel.setArg<cu::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
//            force1Kernel.setArg<cu::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
//            index += 2; // The periodic box size arguments are set when the kernel is executed.
//            force1Kernel.setArg<cl_uint>(index++, maxTiles);
//            if (cu.getSIMDWidth() == 32 || deviceIsCpu)
//                force1Kernel.setArg<cu::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
//        }
//        else
//            force1Kernel.setArg<cl_uint>(index++, cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
//        if (cu.getSIMDWidth() == 32) {
//            force1Kernel.setArg<cu::Buffer>(index++, nb.getExclusionIndices().getDeviceBuffer());
//            force1Kernel.setArg<cu::Buffer>(index++, nb.getExclusionRowIndices().getDeviceBuffer());
//        }
//        program = cu.createProgram(CudaKernelSources::gbsaObcReductions, defines);
//        reduceBornSumKernel = cu::Kernel(program, "reduceBornSum");
//        reduceBornSumKernel.setArg<cl_int>(0, cu.getPaddedNumAtoms());
//        reduceBornSumKernel.setArg<cl_int>(1, nb.getNumForceBuffers());
//        reduceBornSumKernel.setArg<cl_float>(2, 1.0f);
//        reduceBornSumKernel.setArg<cl_float>(3, 0.8f);
//        reduceBornSumKernel.setArg<cl_float>(4, 4.85f);
//        reduceBornSumKernel.setArg<cu::Buffer>(5, (useLong ? longBornSum->getDeviceBuffer() : bornSum->getDeviceBuffer()));
//        reduceBornSumKernel.setArg<cu::Buffer>(6, params->getDeviceBuffer());
//        reduceBornSumKernel.setArg<cu::Buffer>(7, bornRadii->getDeviceBuffer());
//        reduceBornSumKernel.setArg<cu::Buffer>(8, obcChain->getDeviceBuffer());
//        reduceBornForceKernel = cu::Kernel(program, "reduceBornForce");
//        index = 0;
//        reduceBornForceKernel.setArg<cl_int>(index++, cu.getPaddedNumAtoms());
//        reduceBornForceKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
//        reduceBornForceKernel.setArg<cu::Buffer>(index++, bornForce->getDeviceBuffer());
//        if (useLong)
//            reduceBornForceKernel.setArg<cu::Buffer>(index++, longBornForce->getDeviceBuffer());
//        reduceBornForceKernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        reduceBornForceKernel.setArg<cu::Buffer>(index++, params->getDeviceBuffer());
//        reduceBornForceKernel.setArg<cu::Buffer>(index++, bornRadii->getDeviceBuffer());
//        reduceBornForceKernel.setArg<cu::Buffer>(index++, obcChain->getDeviceBuffer());
//    }
//    if (nb.getUseCutoff()) {
//        computeBornSumKernel.setArg<mm_float4>(5, cu.getPeriodicBoxSize());
//        computeBornSumKernel.setArg<mm_float4>(6, cu.getInvPeriodicBoxSize());
//        force1Kernel.setArg<mm_float4>(7, cu.getPeriodicBoxSize());
//        force1Kernel.setArg<mm_float4>(8, cu.getInvPeriodicBoxSize());
//        if (maxTiles < nb.getInteractingTiles().getSize()) {
//            maxTiles = nb.getInteractingTiles().getSize();
//            computeBornSumKernel.setArg<cu::Buffer>(3, nb.getInteractingTiles().getDeviceBuffer());
//            computeBornSumKernel.setArg<cl_uint>(7, maxTiles);
//            force1Kernel.setArg<cu::Buffer>(5, nb.getInteractingTiles().getDeviceBuffer());
//            force1Kernel.setArg<cl_uint>(9, maxTiles);
//            if (cu.getSIMDWidth() == 32 || deviceIsCpu) {
//                computeBornSumKernel.setArg<cu::Buffer>(8, nb.getInteractionFlags().getDeviceBuffer());
//                force1Kernel.setArg<cu::Buffer>(10, nb.getInteractionFlags().getDeviceBuffer());
//            }
//        }
//    }
//    cu.executeKernel(computeBornSumKernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
//    cu.executeKernel(reduceBornSumKernel, cu.getPaddedNumAtoms());
//    cu.executeKernel(force1Kernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
//    cu.executeKernel(reduceBornForceKernel, cu.getPaddedNumAtoms());
//    return 0.0;
//}
//
//void CudaCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
//    // Make sure the new parameters are acceptable.
//    
//    int numParticles = force.getNumParticles();
//    if (numParticles != cu.getNumAtoms())
//        throw OpenMMException("updateParametersInContext: The number of particles has changed");
//    
//    // Record the per-particle parameters.
//    
//    CudaArray<mm_float4>& posq = cu.getPosq();
//    posq.download();
//    vector<mm_float2> paramsVector(numParticles);
//    const double dielectricOffset = 0.009;
//    for (int i = 0; i < numParticles; i++) {
//        double charge, radius, scalingFactor;
//        force.getParticleParameters(i, charge, radius, scalingFactor);
//        radius -= dielectricOffset;
//        paramsVector[i] = mm_float2((float) radius, (float) (scalingFactor*radius));
//        posq[i].w = (float) charge;
//    }
//    posq.upload();
//    params->upload(paramsVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}
//
//class CudaCustomGBForceInfo : public CudaForceInfo {
//public:
//    CudaCustomGBForceInfo(int requiredBuffers, const CustomGBForce& force) : CudaForceInfo(requiredBuffers), force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        vector<double> params1;
//        vector<double> params2;
//        force.getParticleParameters(particle1, params1);
//        force.getParticleParameters(particle2, params2);
//        for (int i = 0; i < (int) params1.size(); i++)
//            if (params1[i] != params2[i])
//                return false;
//        return true;
//    }
//    int getNumParticleGroups() {
//        return force.getNumExclusions();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int particle1, particle2;
//        force.getExclusionParticles(index, particle1, particle2);
//        particles.resize(2);
//        particles[0] = particle1;
//        particles[1] = particle2;
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        return true;
//    }
//private:
//    const CustomGBForce& force;
//};
//
//CudaCalcCustomGBForceKernel::~CudaCalcCustomGBForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (computedValues != NULL)
//        delete computedValues;
//    if (energyDerivs != NULL)
//        delete energyDerivs;
//    if (longEnergyDerivs != NULL)
//        delete longEnergyDerivs;
//    if (globals != NULL)
//        delete globals;
//    if (valueBuffers != NULL)
//        delete valueBuffers;
//    if (longValueBuffers != NULL)
//        delete longValueBuffers;
//    if (tabulatedFunctionParams != NULL)
//        delete tabulatedFunctionParams;
//    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//        delete tabulatedFunctions[i];
//}
//
//void CudaCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
//    if (cu.getPlatformData().contexts.size() > 1)
//        throw OpenMMException("CustomGBForce does not support using multiple CUDA devices");
//    bool useExclusionsForValue = false;
//    numComputedValues = force.getNumComputedValues();
//    vector<string> computedValueNames(force.getNumComputedValues());
//    vector<string> computedValueExpressions(force.getNumComputedValues());
//    if (force.getNumComputedValues() > 0) {
//        CustomGBForce::ComputationType type;
//        force.getComputedValueParameters(0, computedValueNames[0], computedValueExpressions[0], type);
//        if (type == CustomGBForce::SingleParticle)
//            throw OpenMMException("CudaPlatform requires that the first computed value for a CustomGBForce be of type ParticlePair or ParticlePairNoExclusions.");
//        useExclusionsForValue = (type == CustomGBForce::ParticlePair);
//        for (int i = 1; i < force.getNumComputedValues(); i++) {
//            force.getComputedValueParameters(i, computedValueNames[i], computedValueExpressions[i], type);
//            if (type != CustomGBForce::SingleParticle)
//                throw OpenMMException("CudaPlatform requires that a CustomGBForce only have one computed value of type ParticlePair or ParticlePairNoExclusions.");
//        }
//    }
//    int forceIndex;
//    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex)
//        ;
//    string prefix = "custom"+cu.intToString(forceIndex)+"_";
//
//    // Record parameters and exclusions.
//
//    int numParticles = force.getNumParticles();
//    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customGBParameters", true);
//    computedValues = new CudaParameterSet(cu, force.getNumComputedValues(), numParticles, "customGBComputedValues", true);
//    if (force.getNumGlobalParameters() > 0)
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customGBGlobals", false, CL_MEM_READ_ONLY);
//    vector<vector<cl_float> > paramVector(numParticles);
//    vector<vector<int> > exclusionList(numParticles);
//    for (int i = 0; i < numParticles; i++) {
//        vector<double> parameters;
//        force.getParticleParameters(i, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//        exclusionList[i].push_back(i);
//    }
//    for (int i = 0; i < force.getNumExclusions(); i++) {
//        int particle1, particle2;
//        force.getExclusionParticles(i, particle1, particle2);
//        exclusionList[particle1].push_back(particle2);
//        exclusionList[particle2].push_back(particle1);
//    }
//    params->setParameterValues(paramVector);
//
//    // Record the tabulated functions.
//
//    CudaExpressionUtilities::FunctionPlaceholder fp;
//    map<string, Lepton::CustomFunction*> functions;
//    vector<pair<string, string> > functionDefinitions;
//    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
//    stringstream tableArgs;
//    for (int i = 0; i < force.getNumFunctions(); i++) {
//        string name;
//        vector<double> values;
//        double min, max;
//        force.getFunctionParameters(i, name, values, min, max);
//        string arrayName = prefix+"table"+cu.intToString(i);
//        functionDefinitions.push_back(make_pair(name, arrayName));
//        functions[name] = &fp;
//        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), (float) values.size()-2);
//        vector<mm_float4> f = CudaExpressionUtilities::computeFunctionCoefficients(values, min, max);
//        tabulatedFunctions.push_back(new CudaArray<mm_float4>(cu, values.size()-1, "TabulatedFunction"));
//        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
//        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName, "float", 4, sizeof(cl_float4), tabulatedFunctions[tabulatedFunctions.size()-1]->getDeviceBuffer()));
//        tableArgs << ", __global const float4* restrict " << arrayName;
//    }
//    if (force.getNumFunctions() > 0) {
//        tabulatedFunctionParams = new CudaArray<mm_float4>(cu, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
//        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
//        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(prefix+"functionParams", "float", 4, sizeof(cl_float4), tabulatedFunctionParams->getDeviceBuffer()));
//        tableArgs << ", __global const float4* " << prefix << "functionParams";
//    }
//
//    // Record the global parameters.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    if (globals != NULL)
//        globals->upload(globalParamValues);
//
//    // Record derivatives of expressions needed for the chain rule terms.
//
//    vector<vector<Lepton::ParsedExpression> > valueGradientExpressions(force.getNumComputedValues());
//    vector<vector<Lepton::ParsedExpression> > valueDerivExpressions(force.getNumComputedValues());
//    needParameterGradient = false;
//    for (int i = 1; i < force.getNumComputedValues(); i++) {
//        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
//        valueGradientExpressions[i].push_back(ex.differentiate("x").optimize());
//        valueGradientExpressions[i].push_back(ex.differentiate("y").optimize());
//        valueGradientExpressions[i].push_back(ex.differentiate("z").optimize());
//        if (!isZeroExpression(valueGradientExpressions[i][0]) || !isZeroExpression(valueGradientExpressions[i][1]) || !isZeroExpression(valueGradientExpressions[i][2]))
//            needParameterGradient = true;
//         for (int j = 0; j < i; j++)
//            valueDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
//    }
//    vector<vector<Lepton::ParsedExpression> > energyDerivExpressions(force.getNumEnergyTerms());
//    vector<bool> needChainForValue(force.getNumComputedValues(), false);
//    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
//        string expression;
//        CustomGBForce::ComputationType type;
//        force.getEnergyTermParameters(i, expression, type);
//        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
//        for (int j = 0; j < force.getNumComputedValues(); j++) {
//            if (type == CustomGBForce::SingleParticle) {
//                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]).optimize());
//                if (!isZeroExpression(energyDerivExpressions[i].back()))
//                    needChainForValue[j] = true;
//            }
//            else {
//                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"1").optimize());
//                if (!isZeroExpression(energyDerivExpressions[i].back()))
//                    needChainForValue[j] = true;
//                energyDerivExpressions[i].push_back(ex.differentiate(computedValueNames[j]+"2").optimize());
//                if (!isZeroExpression(energyDerivExpressions[i].back()))
//                    needChainForValue[j] = true;
//            }
//        }
//    }
//    bool deviceIsCpu = (cu.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
//    bool useLong = (cu.getSupports64BitGlobalAtomics() && !deviceIsCpu);
//    if (useLong) {
//        longEnergyDerivs = new CudaArray<cl_long>(cu, force.getNumComputedValues()*cu.getPaddedNumAtoms(), "customGBLongEnergyDerivatives");
//        energyDerivs = new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms(), "customGBEnergyDerivatives", true);
//    }
//    else
//        energyDerivs = new CudaParameterSet(cu, force.getNumComputedValues(), cu.getPaddedNumAtoms()*cu.getNonbondedUtilities().getNumForceBuffers(), "customGBEnergyDerivatives", true);
// 
//    // Create the kernels.
//
//    bool useCutoff = (force.getNonbondedMethod() != CustomGBForce::NoCutoff);
//    bool usePeriodic = (force.getNonbondedMethod() != CustomGBForce::NoCutoff && force.getNonbondedMethod() != CustomGBForce::CutoffNonPeriodic);
//    {
//        // Create the N2 value kernel.
//
//        vector<pair<ExpressionTreeNode, string> > variables;
//        map<string, string> rename;
//        ExpressionTreeNode rnode(new Operation::Variable("r"));
//        variables.push_back(make_pair(rnode, "r"));
//        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
//        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
//        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
//            const string& name = force.getPerParticleParameterName(i);
//            variables.push_back(makeVariable(name+"1", "params"+params->getParameterSuffix(i, "1")));
//            variables.push_back(makeVariable(name+"2", "params"+params->getParameterSuffix(i, "2")));
//            rename[name+"1"] = name+"2";
//            rename[name+"2"] = name+"1";
//        }
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = "globals["+cu.intToString(i)+"]";
//            variables.push_back(makeVariable(name, value));
//        }
//        map<string, Lepton::ParsedExpression> n2ValueExpressions;
//        stringstream n2ValueSource;
//        Lepton::ParsedExpression ex = Lepton::Parser::parse(computedValueExpressions[0], functions).optimize();
//        n2ValueExpressions["tempValue1 = "] = ex;
//        n2ValueExpressions["tempValue2 = "] = ex.renameVariables(rename);
//        n2ValueSource << CudaExpressionUtilities::createExpressions(n2ValueExpressions, variables, functionDefinitions, "temp", prefix+"functionParams");
//        map<string, string> replacements;
//        string n2ValueStr = n2ValueSource.str();
//        replacements["COMPUTE_VALUE"] = n2ValueStr;
//        stringstream extraArgs, loadLocal1, loadLocal2, load1, load2;
//        if (force.getNumGlobalParameters() > 0)
//            extraArgs << ", __global const float* globals";
//        pairValueUsesParam.resize(params->getBuffers().size(), false);
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//            string paramName = "params"+cu.intToString(i+1);
//            if (n2ValueStr.find(paramName+"1") != n2ValueStr.npos || n2ValueStr.find(paramName+"2") != n2ValueStr.npos) {
//                extraArgs << ", __global const " << buffer.getType() << "* restrict global_" << paramName << ", __local " << buffer.getType() << "* restrict local_" << paramName;
//                loadLocal1 << "local_" << paramName << "[localAtomIndex] = " << paramName << "1;\n";
//                loadLocal2 << "local_" << paramName << "[localAtomIndex] = global_" << paramName << "[j];\n";
//                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
//                load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
//                pairValueUsesParam[i] = true;
//            }
//        }
//        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
//        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
//        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
//        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
//        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
//        map<string, string> defines;
//        if (cu.getNonbondedUtilities().getForceBufferPerAtomBlock())
//            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
//        if (useCutoff)
//            defines["USE_CUTOFF"] = "1";
//        if (usePeriodic)
//            defines["USE_PERIODIC"] = "1";
//        if (useExclusionsForValue)
//            defines["USE_EXCLUSIONS"] = "1";
//        if (cu.getSIMDWidth() == 32)
//            defines["WARPS_PER_GROUP"] = CudaExpressionUtilities::cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize()/CudaContext::TileSize);
//        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//        defines["NUM_BLOCKS"] = CudaExpressionUtilities::cu.intToString(cu.getNumAtomBlocks());
//        string file;
//        if (deviceIsCpu)
//            file = CudaKernelSources::customGBValueN2_cpu;
//        else if (cu.getSIMDWidth() == 32)
//            file = CudaKernelSources::customGBValueN2_nvidia;
//        else
//            file = CudaKernelSources::customGBValueN2_default;
//        cu::Program program = cu.createProgram(cu.replaceStrings(file, replacements), defines);
//        pairValueKernel = cu::Kernel(program, "computeN2Value");
//        if (useExclusionsForValue)
//            cu.getNonbondedUtilities().requestExclusions(exclusionList);
//    }
//    {
//        // Create the kernel to reduce the N2 value and calculate other values.
//
//        stringstream reductionSource, extraArgs;
//        if (force.getNumGlobalParameters() > 0)
//            extraArgs << ", __global const float* globals";
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//            string paramName = "params"+cu.intToString(i+1);
//            extraArgs << ", __global const " << buffer.getType() << "* restrict " << paramName;
//        }
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
//            string valueName = "values"+cu.intToString(i+1);
//            extraArgs << ", __global " << buffer.getType() << "* restrict global_" << valueName;
//            reductionSource << buffer.getType() << " local_" << valueName << ";\n";
//        }
//        reductionSource << "local_values" << computedValues->getParameterSuffix(0) << " = sum;\n";
//        map<string, string> variables;
//        variables["x"] = "pos.x";
//        variables["y"] = "pos.y";
//        variables["z"] = "pos.z";
//        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
//            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++)
//            variables[force.getGlobalParameterName(i)] = "globals["+cu.intToString(i)+"]";
//        for (int i = 1; i < force.getNumComputedValues(); i++) {
//            variables[computedValueNames[i-1]] = "local_values"+computedValues->getParameterSuffix(i-1);
//            map<string, Lepton::ParsedExpression> valueExpressions;
//            valueExpressions["local_values"+computedValues->getParameterSuffix(i)+" = "] = Lepton::Parser::parse(computedValueExpressions[i], functions).optimize();
//            reductionSource << CudaExpressionUtilities::createExpressions(valueExpressions, variables, functionDefinitions, "value"+cu.intToString(i)+"_temp", prefix+"functionParams");
//        }
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            string valueName = "values"+cu.intToString(i+1);
//            reductionSource << "global_" << valueName << "[index] = local_" << valueName << ";\n";
//        }
//        map<string, string> replacements;
//        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
//        replacements["COMPUTE_VALUES"] = reductionSource.str();
//        map<string, string> defines;
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        cu::Program program = cu.createProgram(cu.replaceStrings(CudaKernelSources::customGBValuePerParticle, replacements), defines);
//        perParticleValueKernel = cu::Kernel(program, "computePerParticleValues");
//    }
//    {
//        // Create the N2 energy kernel.
//
//        vector<pair<ExpressionTreeNode, string> > variables;
//        ExpressionTreeNode rnode(new Operation::Variable("r"));
//        variables.push_back(make_pair(rnode, "r"));
//        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
//        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
//        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
//            const string& name = force.getPerParticleParameterName(i);
//            variables.push_back(makeVariable(name+"1", "params"+params->getParameterSuffix(i, "1")));
//            variables.push_back(makeVariable(name+"2", "params"+params->getParameterSuffix(i, "2")));
//        }
//        for (int i = 0; i < force.getNumComputedValues(); i++) {
//            variables.push_back(makeVariable(computedValueNames[i]+"1", "values"+computedValues->getParameterSuffix(i, "1")));
//            variables.push_back(makeVariable(computedValueNames[i]+"2", "values"+computedValues->getParameterSuffix(i, "2")));
//        }
//        for (int i = 0; i < force.getNumGlobalParameters(); i++)
//            variables.push_back(makeVariable(force.getGlobalParameterName(i), "globals["+cu.intToString(i)+"]"));
//        stringstream n2EnergySource;
//        bool anyExclusions = (force.getNumExclusions() > 0);
//        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
//            string expression;
//            CustomGBForce::ComputationType type;
//            force.getEnergyTermParameters(i, expression, type);
//            if (type == CustomGBForce::SingleParticle)
//                continue;
//            bool exclude = (anyExclusions && type == CustomGBForce::ParticlePair);
//            map<string, Lepton::ParsedExpression> n2EnergyExpressions;
//            n2EnergyExpressions["tempEnergy += "] = Lepton::Parser::parse(expression, functions).optimize();
//            n2EnergyExpressions["dEdR += "] = Lepton::Parser::parse(expression, functions).differentiate("r").optimize();
//            if (useLong) {
//                for (int j = 0; j < force.getNumComputedValues(); j++) {
//                    if (needChainForValue[j]) {
//                        string index = cu.intToString(j+1);
//                        n2EnergyExpressions["/*"+cu.intToString(i+1)+"*/ deriv"+index+"_1 += "] = energyDerivExpressions[i][2*j];
//                        n2EnergyExpressions["/*"+cu.intToString(i+1)+"*/ deriv"+index+"_2 += "] = energyDerivExpressions[i][2*j+1];
//                    }
//                }
//            }
//            else {
//                for (int j = 0; j < force.getNumComputedValues(); j++) {
//                    if (needChainForValue[j]) {
//                        n2EnergyExpressions["/*"+cu.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j, "_1")+" += "] = energyDerivExpressions[i][2*j];
//                        n2EnergyExpressions["/*"+cu.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j, "_2")+" += "] = energyDerivExpressions[i][2*j+1];
//                    }
//                }
//            }
//            if (exclude)
//                n2EnergySource << "if (!isExcluded) {\n";
//            n2EnergySource << CudaExpressionUtilities::createExpressions(n2EnergyExpressions, variables, functionDefinitions, "temp", prefix+"functionParams");
//            if (exclude)
//                n2EnergySource << "}\n";
//        }
//        map<string, string> replacements;
//        string n2EnergyStr = n2EnergySource.str();
//        replacements["COMPUTE_INTERACTION"] = n2EnergyStr;
//        stringstream extraArgs, loadLocal1, loadLocal2, clearLocal, load1, load2, declare1, recordDeriv, storeDerivs1, storeDerivs2, declareTemps, setTemps;
//        if (force.getNumGlobalParameters() > 0)
//            extraArgs << ", __global const float* globals";
//        pairEnergyUsesParam.resize(params->getBuffers().size(), false);
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//            string paramName = "params"+cu.intToString(i+1);
//            if (n2EnergyStr.find(paramName+"1") != n2EnergyStr.npos || n2EnergyStr.find(paramName+"2") != n2EnergyStr.npos) {
//                extraArgs << ", __global const " << buffer.getType() << "* restrict global_" << paramName << ", __local " << buffer.getType() << "* restrict local_" << paramName;
//                loadLocal1 << "local_" << paramName << "[localAtomIndex] = " << paramName << "1;\n";
//                loadLocal2 << "local_" << paramName << "[localAtomIndex] = global_" << paramName << "[j];\n";
//                load1 << buffer.getType() << " " << paramName << "1 = global_" << paramName << "[atom1];\n";
//                load2 << buffer.getType() << " " << paramName << "2 = local_" << paramName << "[atom2];\n";
//                pairEnergyUsesParam[i] = true;
//            }
//        }
//        pairEnergyUsesValue.resize(computedValues->getBuffers().size(), false);
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
//            string valueName = "values"+cu.intToString(i+1);
//            if (n2EnergyStr.find(valueName+"1") != n2EnergyStr.npos || n2EnergyStr.find(valueName+"2") != n2EnergyStr.npos) {
//                extraArgs << ", __global const " << buffer.getType() << "* restrict global_" << valueName << ", __local " << buffer.getType() << "* restrict local_" << valueName;
//                loadLocal1 << "local_" << valueName << "[localAtomIndex] = " << valueName << "1;\n";
//                loadLocal2 << "local_" << valueName << "[localAtomIndex] = global_" << valueName << "[j];\n";
//                load1 << buffer.getType() << " " << valueName << "1 = global_" << valueName << "[atom1];\n";
//                load2 << buffer.getType() << " " << valueName << "2 = local_" << valueName << "[atom2];\n";
//                pairEnergyUsesValue[i] = true;
//            }
//        }
//        if (useLong) {
//            extraArgs << ", __global long* restrict derivBuffers";
//            for (int i = 0; i < force.getNumComputedValues(); i++) {
//                string index = cu.intToString(i+1);
//                extraArgs << ", __local float* restrict local_deriv" << index;
//                clearLocal << "local_deriv" << index << "[localAtomIndex] = 0.0f;\n";
//                declare1 << "float deriv" << index << "_1 = 0.0f;\n";
//                load2 << "float deriv" << index << "_2 = 0.0f;\n";
//                recordDeriv << "local_deriv" << index << "[atom2] += deriv" << index << "_2;\n";
//                storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")\n";
//                storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")\n";
//                declareTemps << "__local float tempDerivBuffer" << index << "[64];\n";
//                setTemps << "tempDerivBuffer" << index << "[get_local_id(0)] = deriv" << index << "_1;\n";
//            }
//        }
//        else {
//            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//                const CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
//                string index = cu.intToString(i+1);
//                extraArgs << ", __global " << buffer.getType() << "* restrict derivBuffers" << index << ", __local " << buffer.getType() << "* restrict local_deriv" << index;
//                clearLocal << "local_deriv" << index << "[localAtomIndex] = 0.0f;\n";
//                declare1 << buffer.getType() << " deriv" << index << "_1 = 0.0f;\n";
//                load2 << buffer.getType() << " deriv" << index << "_2 = 0.0f;\n";
//                recordDeriv << "local_deriv" << index << "[atom2] += deriv" << index << "_2;\n";
//                storeDerivs1 << "STORE_DERIVATIVE_1(" << index << ")\n";
//                storeDerivs2 << "STORE_DERIVATIVE_2(" << index << ")\n";
//                declareTemps << "__local " << buffer.getType() << " tempDerivBuffer" << index << "[64];\n";
//                setTemps << "tempDerivBuffer" << index << "[get_local_id(0)] = deriv" << index << "_1;\n";
//            }
//        }
//        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
//        replacements["LOAD_LOCAL_PARAMETERS_FROM_1"] = loadLocal1.str();
//        replacements["LOAD_LOCAL_PARAMETERS_FROM_GLOBAL"] = loadLocal2.str();
//        replacements["CLEAR_LOCAL_DERIVATIVES"] = clearLocal.str();
//        replacements["LOAD_ATOM1_PARAMETERS"] = load1.str();
//        replacements["LOAD_ATOM2_PARAMETERS"] = load2.str();
//        replacements["DECLARE_ATOM1_DERIVATIVES"] = declare1.str();
//        replacements["RECORD_DERIVATIVE_2"] = recordDeriv.str();
//        replacements["STORE_DERIVATIVES_1"] = storeDerivs1.str();
//        replacements["STORE_DERIVATIVES_2"] = storeDerivs2.str();
//        replacements["DECLARE_TEMP_BUFFERS"] = declareTemps.str();
//        replacements["SET_TEMP_BUFFERS"] = setTemps.str();
//        map<string, string> defines;
//        if (cu.getNonbondedUtilities().getForceBufferPerAtomBlock())
//            defines["USE_OUTPUT_BUFFER_PER_BLOCK"] = "1";
//        if (useCutoff)
//            defines["USE_CUTOFF"] = "1";
//        if (usePeriodic)
//            defines["USE_PERIODIC"] = "1";
//        if (anyExclusions)
//            defines["USE_EXCLUSIONS"] = "1";
//        if (cu.getSIMDWidth() == 32)
//            defines["WARPS_PER_GROUP"] = CudaExpressionUtilities::cu.intToString(cu.getNonbondedUtilities().getForceThreadBlockSize()/CudaContext::TileSize);
//        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//        defines["NUM_BLOCKS"] = CudaExpressionUtilities::cu.intToString(cu.getNumAtomBlocks());
//        string file;
//        if (deviceIsCpu)
//            file = CudaKernelSources::customGBEnergyN2_cpu;
//        else if (cu.getSIMDWidth() == 32)
//            file = CudaKernelSources::customGBEnergyN2_nvidia;
//        else
//            file = CudaKernelSources::customGBEnergyN2_default;
//        cu::Program program = cu.createProgram(cu.replaceStrings(file, replacements), defines);
//        pairEnergyKernel = cu::Kernel(program, "computeN2Energy");
//    }
//    {
//        // Create the kernel to reduce the derivatives and calculate per-particle energy terms.
//
//        stringstream compute, extraArgs, reduce;
//        if (force.getNumGlobalParameters() > 0)
//            extraArgs << ", __global const float* globals";
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//            string paramName = "params"+cu.intToString(i+1);
//            extraArgs << ", __global const " << buffer.getType() << "* restrict " << paramName;
//        }
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
//            string valueName = "values"+cu.intToString(i+1);
//            extraArgs << ", __global const " << buffer.getType() << "* restrict " << valueName;
//        }
//        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
//            string index = cu.intToString(i+1);
//            extraArgs << ", __global " << buffer.getType() << "* restrict derivBuffers" << index;
//            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
//        }
//        if (useLong) {
//            extraArgs << ", __global const long* restrict derivBuffersIn";
//            for (int i = 0; i < energyDerivs->getNumParameters(); ++i)
//                reduce << "derivBuffers" << energyDerivs->getParameterSuffix(i, "[index]") <<
//                        " = (1.0f/0xFFFFFFFF)*derivBuffersIn[index+PADDED_NUM_ATOMS*" << cu.intToString(i) << "];\n";
//        }
//        else {
//            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
//                reduce << "REDUCE_VALUE(derivBuffers" << cu.intToString(i+1) << ", " << energyDerivs->getBuffers()[i].getType() << ")\n";
//        }
//        
//        // Compute the various expressions.
//        
//        map<string, string> variables;
//        variables["x"] = "pos.x";
//        variables["y"] = "pos.y";
//        variables["z"] = "pos.z";
//        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
//            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++)
//            variables[force.getGlobalParameterName(i)] = "globals["+cu.intToString(i)+"]";
//        for (int i = 0; i < force.getNumComputedValues(); i++)
//            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
//        map<string, Lepton::ParsedExpression> expressions;
//        for (int i = 0; i < force.getNumEnergyTerms(); i++) {
//            string expression;
//            CustomGBForce::ComputationType type;
//            force.getEnergyTermParameters(i, expression, type);
//            if (type != CustomGBForce::SingleParticle)
//                continue;
//            Lepton::ParsedExpression parsed = Lepton::Parser::parse(expression, functions).optimize();
//            expressions["/*"+cu.intToString(i+1)+"*/ energy += "] = parsed;
//            for (int j = 0; j < force.getNumComputedValues(); j++)
//                expressions["/*"+cu.intToString(i+1)+"*/ deriv"+energyDerivs->getParameterSuffix(j)+" += "] = energyDerivExpressions[i][j];
//            Lepton::ParsedExpression gradx = parsed.differentiate("x").optimize();
//            Lepton::ParsedExpression grady = parsed.differentiate("y").optimize();
//            Lepton::ParsedExpression gradz = parsed.differentiate("z").optimize();
//            if (!isZeroExpression(gradx))
//                expressions["/*"+cu.intToString(i+1)+"*/ force.x -= "] = gradx;
//            if (!isZeroExpression(grady))
//                expressions["/*"+cu.intToString(i+1)+"*/ force.y -= "] = grady;
//            if (!isZeroExpression(gradz))
//                expressions["/*"+cu.intToString(i+1)+"*/ force.z -= "] = gradz;
//        }
//        for (int i = 1; i < force.getNumComputedValues(); i++)
//            for (int j = 0; j < i; j++)
//                expressions["float dV"+cu.intToString(i)+"dV"+cu.intToString(j)+" = "] = valueDerivExpressions[i][j];
//        compute << CudaExpressionUtilities::createExpressions(expressions, variables, functionDefinitions, "temp", prefix+"functionParams");
//        
//        // Record values.
//        
//        compute << "forceBuffers[index] = forceBuffers[index]+force;\n";
//        for (int i = 1; i < force.getNumComputedValues(); i++) {
//            compute << "float totalDeriv"<<i<<" = dV"<<i<<"dV0";
//            for (int j = 1; j < i; j++)
//                compute << " + totalDeriv"<<j<<"*dV"<<i<<"dV"<<j;
//            compute << ";\n";
//            compute << "deriv"<<(i+1)<<" *= totalDeriv"<<i<<";\n";
//        }
//        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//            string index = cu.intToString(i+1);
//            compute << "derivBuffers" << index << "[index] = deriv" << index << ";\n";
//        }
//        map<string, string> replacements;
//        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
//        replacements["REDUCE_DERIVATIVES"] = reduce.str();
//        replacements["COMPUTE_ENERGY"] = compute.str();
//        map<string, string> defines;
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//        cu::Program program = cu.createProgram(cu.replaceStrings(CudaKernelSources::customGBEnergyPerParticle, replacements), defines);
//        perParticleEnergyKernel = cu::Kernel(program, "computePerParticleEnergy");
//    }
//    if (needParameterGradient) {
//        // Create the kernel to compute chain rule terms for computed values that depend explicitly on particle coordinates.
//
//        stringstream compute, extraArgs;
//        if (force.getNumGlobalParameters() > 0)
//            extraArgs << ", __global const float* globals";
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//            string paramName = "params"+cu.intToString(i+1);
//            extraArgs << ", __global const " << buffer.getType() << "* restrict " << paramName;
//        }
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
//            string valueName = "values"+cu.intToString(i+1);
//            extraArgs << ", __global const " << buffer.getType() << "* restrict " << valueName;
//        }
//        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
//            string index = cu.intToString(i+1);
//            extraArgs << ", __global " << buffer.getType() << "* restrict derivBuffers" << index;
//            compute << buffer.getType() << " deriv" << index << " = derivBuffers" << index << "[index];\n";
//        }
//        map<string, string> variables;
//        variables["x"] = "pos.x";
//        variables["y"] = "pos.y";
//        variables["z"] = "pos.z";
//        for (int i = 0; i < force.getNumPerParticleParameters(); i++)
//            variables[force.getPerParticleParameterName(i)] = "params"+params->getParameterSuffix(i, "[index]");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++)
//            variables[force.getGlobalParameterName(i)] = "globals["+cu.intToString(i)+"]";
//        for (int i = 0; i < force.getNumComputedValues(); i++)
//            variables[computedValueNames[i]] = "values"+computedValues->getParameterSuffix(i, "[index]");
//        for (int i = 1; i < force.getNumComputedValues(); i++) {
//            string is = cu.intToString(i);
//            compute << "float4 dV"<<is<<"dR = (float4) 0;\n";
//            for (int j = 1; j < i; j++) {
//                if (!isZeroExpression(valueDerivExpressions[i][j])) {
//                    map<string, Lepton::ParsedExpression> derivExpressions;
//                    string js = cu.intToString(j);
//                    derivExpressions["float dV"+is+"dV"+js+" = "] = valueDerivExpressions[i][j];
//                    compute << CudaExpressionUtilities::createExpressions(derivExpressions, variables, functionDefinitions, "temp_"+is+"_"+js, prefix+"functionParams");
//                    compute << "dV"<<is<<"dR += dV"<<is<<"dV"<<js<<"*dV"<<js<<"dR;\n";
//                }
//            }
//            map<string, Lepton::ParsedExpression> gradientExpressions;
//            if (!isZeroExpression(valueGradientExpressions[i][0]))
//                gradientExpressions["dV"+is+"dR.x += "] = valueGradientExpressions[i][0];
//            if (!isZeroExpression(valueGradientExpressions[i][1]))
//                gradientExpressions["dV"+is+"dR.y += "] = valueGradientExpressions[i][1];
//            if (!isZeroExpression(valueGradientExpressions[i][2]))
//                gradientExpressions["dV"+is+"dR.z += "] = valueGradientExpressions[i][2];
//            compute << CudaExpressionUtilities::createExpressions(gradientExpressions, variables, functionDefinitions, "temp", prefix+"functionParams");
//        }
//        for (int i = 1; i < force.getNumComputedValues(); i++) {
//            string is = cu.intToString(i);
//            compute << "force -= deriv"<<energyDerivs->getParameterSuffix(i)<<"*dV"<<is<<"dR;\n";
//        }
//        map<string, string> replacements;
//        replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
//        replacements["COMPUTE_FORCES"] = compute.str();
//        map<string, string> defines;
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        cu::Program program = cu.createProgram(cu.replaceStrings(CudaKernelSources::customGBGradientChainRule, replacements), defines);
//        gradientChainRuleKernel = cu::Kernel(program, "computeGradientChainRuleTerms");
//    }
//    {
//        // Create the code to calculate chain rules terms as part of the default nonbonded kernel.
//
//        vector<pair<ExpressionTreeNode, string> > globalVariables;
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = "globals["+cu.intToString(i)+"]";
//            globalVariables.push_back(makeVariable(name, prefix+value));
//        }
//        vector<pair<ExpressionTreeNode, string> > variables = globalVariables;
//        map<string, string> rename;
//        ExpressionTreeNode rnode(new Operation::Variable("r"));
//        variables.push_back(make_pair(rnode, "r"));
//        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Square(), rnode), "r2"));
//        variables.push_back(make_pair(ExpressionTreeNode(new Operation::Reciprocal(), rnode), "invR"));
//        for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
//            const string& name = force.getPerParticleParameterName(i);
//            variables.push_back(makeVariable(name+"1", prefix+"params"+params->getParameterSuffix(i, "1")));
//            variables.push_back(makeVariable(name+"2", prefix+"params"+params->getParameterSuffix(i, "2")));
//            rename[name+"1"] = name+"2";
//            rename[name+"2"] = name+"1";
//        }
//        map<string, Lepton::ParsedExpression> derivExpressions;
//        stringstream chainSource;
//        Lepton::ParsedExpression dVdR = Lepton::Parser::parse(computedValueExpressions[0], functions).differentiate("r").optimize();
//        derivExpressions["float dV0dR1 = "] = dVdR;
//        derivExpressions["float dV0dR2 = "] = dVdR.renameVariables(rename);
//        chainSource << CudaExpressionUtilities::createExpressions(derivExpressions, variables, functionDefinitions, prefix+"temp0_", prefix+"functionParams");
//        if (needChainForValue[0]) {
//            if (useExclusionsForValue)
//                chainSource << "if (!isExcluded) {\n";
//            chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "1") << ";\n";
//            chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(0, "2") << ";\n";
//            if (useExclusionsForValue)
//                chainSource << "}\n";
//        }
//        for (int i = 1; i < force.getNumComputedValues(); i++) {
//            if (needChainForValue[i]) {
//                chainSource << "tempForce -= dV0dR1*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "1") << ";\n";
//                chainSource << "tempForce -= dV0dR2*" << prefix << "dEdV" << energyDerivs->getParameterSuffix(i, "2") << ";\n";
//            }
//        }
//        map<string, string> replacements;
//        string chainStr = chainSource.str();
//        replacements["COMPUTE_FORCE"] = chainStr;
//        string source = cu.replaceStrings(CudaKernelSources::customGBChainRule, replacements);
//        vector<CudaNonbondedUtilities::ParameterInfo> parameters;
//        vector<CudaNonbondedUtilities::ParameterInfo> arguments;
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//            string paramName = prefix+"params"+cu.intToString(i+1);
//            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
//                parameters.push_back(CudaNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
//        }
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
//            string paramName = prefix+"values"+cu.intToString(i+1);
//            if (chainStr.find(paramName+"1") != chainStr.npos || chainStr.find(paramName+"2") != chainStr.npos)
//                parameters.push_back(CudaNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
//        }
//        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//            if (needChainForValue[i]) { 
//                const CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
//                string paramName = prefix+"dEdV"+cu.intToString(i+1);
//                parameters.push_back(CudaNonbondedUtilities::ParameterInfo(paramName, buffer.getComponentType(), buffer.getNumComponents(), buffer.getSize(), buffer.getMemory()));
//            }
//        }
//        if (globals != NULL) {
//            globals->upload(globalParamValues);
//            arguments.push_back(CudaNonbondedUtilities::ParameterInfo(prefix+"globals", "float", 1, sizeof(cl_float), globals->getDeviceBuffer()));
//        }
//        cu.getNonbondedUtilities().addInteraction(useCutoff, usePeriodic, force.getNumExclusions() > 0, force.getCutoffDistance(), exclusionList, source, force.getForceGroup());
//        for (int i = 0; i < (int) parameters.size(); i++)
//            cu.getNonbondedUtilities().addParameter(parameters[i]);
//        for (int i = 0; i < (int) arguments.size(); i++)
//            cu.getNonbondedUtilities().addArgument(arguments[i]);
//    }
//    cu.addForce(new CudaCustomGBForceInfo(cu.getNonbondedUtilities().getNumForceBuffers(), force));
//    if (useLong)
//        cu.addAutoclearBuffer(longEnergyDerivs->getDeviceBuffer(), 2*longEnergyDerivs->getSize());
//    else {
//        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
//            cu.addAutoclearBuffer(buffer.getMemory(), buffer.getSize()*energyDerivs->getNumObjects()/sizeof(cl_float));
//        }
//    }
//}
//
//double CudaCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    bool deviceIsCpu = (cu.getDevice().getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
//    CudaNonbondedUtilities& nb = cu.getNonbondedUtilities();
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        maxTiles = (nb.getUseCutoff() ? nb.getInteractingTiles().getSize() : 0);
//        bool useLong = (cu.getSupports64BitGlobalAtomics() && !deviceIsCpu);
//        if (useLong) {
//            longValueBuffers = new CudaArray<cl_long>(cu, cu.getPaddedNumAtoms(), "customGBLongValueBuffers");
//            cu.addAutoclearBuffer(longValueBuffers->getDeviceBuffer(), 2*longValueBuffers->getSize());
//            cu.clearBuffer(longValueBuffers->getDeviceBuffer(), 2*longValueBuffers->getSize());
//        }
//        else {
//            valueBuffers = new CudaArray<cl_float>(cu, cu.getPaddedNumAtoms()*nb.getNumForceBuffers(), "customGBValueBuffers");
//            cu.addAutoclearBuffer(valueBuffers->getDeviceBuffer(), valueBuffers->getSize());
//            cu.clearBuffer(*valueBuffers);
//        }
//        int index = 0;
//        pairValueKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        pairValueKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*sizeof(cl_float4), NULL);
//        pairValueKernel.setArg<cu::Buffer>(index++, cu.getNonbondedUtilities().getExclusions().getDeviceBuffer());
//        pairValueKernel.setArg<cu::Buffer>(index++, cu.getNonbondedUtilities().getExclusionIndices().getDeviceBuffer());
//        pairValueKernel.setArg<cu::Buffer>(index++, cu.getNonbondedUtilities().getExclusionRowIndices().getDeviceBuffer());
//        pairValueKernel.setArg<cu::Buffer>(index++, useLong ? longValueBuffers->getDeviceBuffer() : valueBuffers->getDeviceBuffer());
//        pairValueKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*sizeof(cl_float), NULL);
//        /// \todo Eliminate this argument and make local to the kernel. For *_default.cu kernel can actually make it TileSize rather than getForceThreadBlockSize as only half the workgroup stores to it as was done with nonbonded_default.cu.
//        /// \todo Also make the previous __local argument local as was done with nonbonded_default.cu.
//        pairValueKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*sizeof(cl_float), NULL);
//        if (nb.getUseCutoff()) {
//            pairValueKernel.setArg<cu::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
//            pairValueKernel.setArg<cu::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
//            index += 2; // Periodic box size arguments are set when the kernel is executed.
//            pairValueKernel.setArg<cl_uint>(index++, maxTiles);
//            if (cu.getSIMDWidth() == 32 || deviceIsCpu)
//                pairValueKernel.setArg<cu::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
//        }
//        else
//            pairValueKernel.setArg<cl_uint>(index++, cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
//        if (globals != NULL)
//            pairValueKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            if (pairValueUsesParam[i]) {
//                const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//                pairValueKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//                pairValueKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
//            }
//        }
//        if (tabulatedFunctionParams != NULL) {
//            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//                pairValueKernel.setArg<cu::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
//            pairValueKernel.setArg<cu::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
//        }
//        index = 0;
//        perParticleValueKernel.setArg<cl_int>(index++, cu.getPaddedNumAtoms());
//        perParticleValueKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
//        perParticleValueKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        perParticleValueKernel.setArg<cu::Buffer>(index++, useLong ? longValueBuffers->getDeviceBuffer() : valueBuffers->getDeviceBuffer());
//        if (globals != NULL)
//            perParticleValueKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//        for (int i = 0; i < (int) params->getBuffers().size(); i++)
//            perParticleValueKernel.setArg<cu::Memory>(index++, params->getBuffers()[i].getMemory());
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
//            perParticleValueKernel.setArg<cu::Memory>(index++, computedValues->getBuffers()[i].getMemory());
//        if (tabulatedFunctionParams != NULL) {
//            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//                perParticleValueKernel.setArg<cu::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
//            perParticleValueKernel.setArg<cu::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
//        }
//        index = 0;
//        pairEnergyKernel.setArg<cu::Buffer>(index++, useLong ? cu.getLongForceBuffer().getDeviceBuffer() : cu.getForceBuffers().getDeviceBuffer());
//        pairEnergyKernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        pairEnergyKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*sizeof(cl_float4), NULL);
//        pairEnergyKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        pairEnergyKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*sizeof(cl_float4), NULL);
//        pairEnergyKernel.setArg<cu::Buffer>(index++, cu.getNonbondedUtilities().getExclusions().getDeviceBuffer());
//        pairEnergyKernel.setArg<cu::Buffer>(index++, cu.getNonbondedUtilities().getExclusionIndices().getDeviceBuffer());
//        pairEnergyKernel.setArg<cu::Buffer>(index++, cu.getNonbondedUtilities().getExclusionRowIndices().getDeviceBuffer());
//        /// \todo Eliminate this argument and make local to the kernel. For *_default.cu kernel can actually make it TileSize rather than getForceThreadBlockSize as only half the workgroup stores to it as was done with nonbonded_default.cu.
//        /// \todo Also make the previous __local argument local as was done with nonbonded_default.cu.
//        pairEnergyKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*sizeof(cl_float4), NULL);
//        if (nb.getUseCutoff()) {
//            pairEnergyKernel.setArg<cu::Buffer>(index++, nb.getInteractingTiles().getDeviceBuffer());
//            pairEnergyKernel.setArg<cu::Buffer>(index++, nb.getInteractionCount().getDeviceBuffer());
//            index += 2; // Periodic box size arguments are set when the kernel is executed.
//            pairEnergyKernel.setArg<cl_uint>(index++, maxTiles);
//            if (cu.getSIMDWidth() == 32 || deviceIsCpu)
//                pairEnergyKernel.setArg<cu::Buffer>(index++, nb.getInteractionFlags().getDeviceBuffer());
//        }
//        else
//            pairEnergyKernel.setArg<cl_uint>(index++, cu.getNumAtomBlocks()*(cu.getNumAtomBlocks()+1)/2);
//        if (globals != NULL)
//            pairEnergyKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//        for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//            if (pairEnergyUsesParam[i]) {
//                const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//                pairEnergyKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//                pairEnergyKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
//            }
//        }
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++) {
//            if (pairEnergyUsesValue[i]) {
//                const CudaNonbondedUtilities::ParameterInfo& buffer = computedValues->getBuffers()[i];
//                pairEnergyKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//                pairEnergyKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
//            }
//        }
//        if (useLong) {
//            pairEnergyKernel.setArg<cu::Memory>(index++, longEnergyDerivs->getDeviceBuffer());
//            for (int i = 0; i < numComputedValues; ++i)
//                pairEnergyKernel.setArg(index++, nb.getForceThreadBlockSize()*sizeof(cl_float), NULL);
//        }
//        else {
//            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++) {
//                const CudaNonbondedUtilities::ParameterInfo& buffer = energyDerivs->getBuffers()[i];
//                pairEnergyKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//                pairEnergyKernel.setArg(index++, (deviceIsCpu ? CudaContext::TileSize : nb.getForceThreadBlockSize())*buffer.getSize(), NULL);
//            }
//        }
//        if (tabulatedFunctionParams != NULL) {
//            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//                pairEnergyKernel.setArg<cu::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
//            pairEnergyKernel.setArg<cu::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
//        }
//        index = 0;
//        perParticleEnergyKernel.setArg<cl_int>(index++, cu.getPaddedNumAtoms());
//        perParticleEnergyKernel.setArg<cl_int>(index++, nb.getNumForceBuffers());
//        perParticleEnergyKernel.setArg<cu::Buffer>(index++, cu.getForceBuffers().getDeviceBuffer());
//        perParticleEnergyKernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        perParticleEnergyKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        if (globals != NULL)
//            perParticleEnergyKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//        for (int i = 0; i < (int) params->getBuffers().size(); i++)
//            perParticleEnergyKernel.setArg<cu::Memory>(index++, params->getBuffers()[i].getMemory());
//        for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
//            perParticleEnergyKernel.setArg<cu::Memory>(index++, computedValues->getBuffers()[i].getMemory());
//        for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
//            perParticleEnergyKernel.setArg<cu::Memory>(index++, energyDerivs->getBuffers()[i].getMemory());
//        if (useLong)
//            perParticleEnergyKernel.setArg<cu::Memory>(index++, longEnergyDerivs->getDeviceBuffer());
//        if (tabulatedFunctionParams != NULL) {
//            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//                perParticleEnergyKernel.setArg<cu::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
//            perParticleEnergyKernel.setArg<cu::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
//        }
//        if (needParameterGradient) {
//            index = 0;
//            gradientChainRuleKernel.setArg<cu::Buffer>(index++, cu.getForceBuffers().getDeviceBuffer());
//            gradientChainRuleKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//            if (globals != NULL)
//                gradientChainRuleKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//            for (int i = 0; i < (int) params->getBuffers().size(); i++)
//                gradientChainRuleKernel.setArg<cu::Memory>(index++, params->getBuffers()[i].getMemory());
//            for (int i = 0; i < (int) computedValues->getBuffers().size(); i++)
//                gradientChainRuleKernel.setArg<cu::Memory>(index++, computedValues->getBuffers()[i].getMemory());
//            for (int i = 0; i < (int) energyDerivs->getBuffers().size(); i++)
//                gradientChainRuleKernel.setArg<cu::Memory>(index++, energyDerivs->getBuffers()[i].getMemory());
//        }
//    }
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    if (nb.getUseCutoff()) {
//        pairValueKernel.setArg<mm_float4>(10, cu.getPeriodicBoxSize());
//        pairValueKernel.setArg<mm_float4>(11, cu.getInvPeriodicBoxSize());
//        pairEnergyKernel.setArg<mm_float4>(11, cu.getPeriodicBoxSize());
//        pairEnergyKernel.setArg<mm_float4>(12, cu.getInvPeriodicBoxSize());
//        if (maxTiles < nb.getInteractingTiles().getSize()) {
//            maxTiles = nb.getInteractingTiles().getSize();
//            pairValueKernel.setArg<cu::Buffer>(8, nb.getInteractingTiles().getDeviceBuffer());
//            pairValueKernel.setArg<cl_uint>(12, maxTiles);
//            pairEnergyKernel.setArg<cu::Buffer>(9, nb.getInteractingTiles().getDeviceBuffer());
//            pairEnergyKernel.setArg<cl_uint>(13, maxTiles);
//            if (cu.getSIMDWidth() == 32 || deviceIsCpu) {
//                pairValueKernel.setArg<cu::Buffer>(13, nb.getInteractionFlags().getDeviceBuffer());
//                pairEnergyKernel.setArg<cu::Buffer>(14, nb.getInteractionFlags().getDeviceBuffer());
//            }
//        }
//    }
//    cu.executeKernel(pairValueKernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
//    cu.executeKernel(perParticleValueKernel, cu.getPaddedNumAtoms());
//    cu.executeKernel(pairEnergyKernel, nb.getNumForceThreadBlocks()*nb.getForceThreadBlockSize(), nb.getForceThreadBlockSize());
//    cu.executeKernel(perParticleEnergyKernel, cu.getPaddedNumAtoms());
//    if (needParameterGradient)
//        cu.executeKernel(gradientChainRuleKernel, cu.getPaddedNumAtoms());
//    return 0.0;
//}
//
//void CudaCalcCustomGBForceKernel::copyParametersToContext(ContextImpl& context, const CustomGBForce& force) {
//    int numParticles = force.getNumParticles();
//    if (numParticles != cu.getNumAtoms())
//        throw OpenMMException("updateParametersInContext: The number of particles has changed");
//    
//    // Record the per-particle parameters.
//    
//    vector<vector<cl_float> > paramVector(numParticles);
//    vector<double> parameters;
//    for (int i = 0; i < numParticles; i++) {
//        force.getParticleParameters(i, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}
//
//class CudaCustomExternalForceInfo : public CudaForceInfo {
//public:
//    CudaCustomExternalForceInfo(const CustomExternalForce& force, int numParticles) : CudaForceInfo(0), force(force), indices(numParticles, -1) {
//        vector<double> params;
//        for (int i = 0; i < force.getNumParticles(); i++) {
//            int particle;
//            force.getParticleParameters(i, particle, params);
//            indices[particle] = i;
//        }
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        particle1 = indices[particle1];
//        particle2 = indices[particle2];
//        if (particle1 == -1 && particle2 == -1)
//            return true;
//        if (particle1 == -1 || particle2 == -1)
//            return false;
//        int temp;
//        vector<double> params1;
//        vector<double> params2;
//        force.getParticleParameters(particle1, temp, params1);
//        force.getParticleParameters(particle2, temp, params2);
//        for (int i = 0; i < (int) params1.size(); i++)
//            if (params1[i] != params2[i])
//                return false;
//        return true;
//    }
//private:
//    const CustomExternalForce& force;
//    vector<int> indices;
//};
//
//CudaCalcCustomExternalForceKernel::~CudaCalcCustomExternalForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (globals != NULL)
//        delete globals;
//}
//
//void CudaCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumParticles()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumParticles()/numContexts;
//    numParticles = endIndex-startIndex;
//    if (numParticles == 0)
//        return;
//    vector<vector<int> > atoms(numParticles, vector<int>(1));
//    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customExternalParams");
//    vector<vector<cl_float> > paramVector(numParticles);
//    for (int i = 0; i < numParticles; i++) {
//        vector<double> parameters;
//        force.getParticleParameters(startIndex+i, atoms[i][0], parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    cu.addForce(new CudaCustomExternalForceInfo(force, system.getNumParticles()));
//
//    // Record information for the expressions.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    Lepton::ParsedExpression energyExpression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
//    Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x").optimize();
//    Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y").optimize();
//    Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z").optimize();
//    map<string, Lepton::ParsedExpression> expressions;
//    expressions["energy += "] = energyExpression;
//    expressions["float dEdX = "] = forceExpressionX;
//    expressions["float dEdY = "] = forceExpressionY;
//    expressions["float dEdZ = "] = forceExpressionZ;
//
//    // Create the kernels.
//
//    map<string, string> variables;
//    variables["x"] = "pos1.x";
//    variables["y"] = "pos1.y";
//    variables["z"] = "pos1.z";
//    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
//        const string& name = force.getPerParticleParameterName(i);
//        variables[name] = "particleParams"+params->getParameterSuffix(i);
//    }
//    if (force.getNumGlobalParameters() > 0) {
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customExternalGlobals", false, CL_MEM_READ_ONLY);
//        globals->upload(globalParamValues);
//        string argName = cu.getBondedUtilities().addArgument(globals->getDeviceBuffer(), "float");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = argName+"["+cu.intToString(i)+"]";
//            variables[name] = value;
//        }
//    }
//    stringstream compute;
//    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
//        compute<<buffer.getType()<<" particleParams"<<(i+1)<<" = "<<argName<<"[index];\n";
//    }
//    vector<pair<string, string> > functions;
//    compute << CudaExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
//    map<string, string> replacements;
//    replacements["COMPUTE_FORCE"] = compute.str();
//    cu.getBondedUtilities().addInteraction(atoms, cu.replaceStrings(CudaKernelSources::customExternalForce, replacements), force.getForceGroup());
//}
//
//double CudaCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    return 0.0;
//}
//
//void CudaCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumParticles()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumParticles()/numContexts;
//    if (numParticles != endIndex-startIndex)
//        throw OpenMMException("updateParametersInContext: The number of particles has changed");
//    
//    // Record the per-particle parameters.
//    
//    vector<vector<cl_float> > paramVector(numParticles);
//    vector<double> parameters;
//    for (int i = 0; i < numParticles; i++) {
//        int particle;
//        force.getParticleParameters(startIndex+i, particle, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}
//
//class CudaCustomHbondForceInfo : public CudaForceInfo {
//public:
//    CudaCustomHbondForceInfo(int requiredBuffers, const CustomHbondForce& force) : CudaForceInfo(requiredBuffers), force(force) {
//    }
//    bool areParticlesIdentical(int particle1, int particle2) {
//        return true;
//    }
//    int getNumParticleGroups() {
//        return force.getNumDonors()+force.getNumAcceptors()+force.getNumExclusions();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        int p1, p2, p3;
//        vector<double> parameters;
//        if (index < force.getNumDonors()) {
//            force.getDonorParameters(index, p1, p2, p3, parameters);
//            particles.clear();
//            particles.push_back(p1);
//            if (p2 > -1)
//                particles.push_back(p2);
//            if (p3 > -1)
//                particles.push_back(p3);
//            return;
//        }
//        index -= force.getNumDonors();
//        if (index < force.getNumAcceptors()) {
//            force.getAcceptorParameters(index, p1, p2, p3, parameters);
//            particles.clear();
//            particles.push_back(p1);
//            if (p2 > -1)
//                particles.push_back(p2);
//            if (p3 > -1)
//                particles.push_back(p3);
//            return;
//        }
//        index -= force.getNumAcceptors();
//        int donor, acceptor;
//        force.getExclusionParticles(index, donor, acceptor);
//        particles.clear();
//        force.getDonorParameters(donor, p1, p2, p3, parameters);
//        particles.push_back(p1);
//        if (p2 > -1)
//            particles.push_back(p2);
//        if (p3 > -1)
//            particles.push_back(p3);
//        force.getAcceptorParameters(acceptor, p1, p2, p3, parameters);
//        particles.push_back(p1);
//        if (p2 > -1)
//            particles.push_back(p2);
//        if (p3 > -1)
//            particles.push_back(p3);
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        int p1, p2, p3;
//        vector<double> params1, params2;
//        if (group1 < force.getNumDonors() && group2 < force.getNumDonors()) {
//            force.getDonorParameters(group1, p1, p2, p3, params1);
//            force.getDonorParameters(group2, p1, p2, p3, params2);
//            return (params1 == params2 && params1 == params2);
//        }
//        if (group1 < force.getNumDonors() || group2 < force.getNumDonors())
//            return false;
//        group1 -= force.getNumDonors();
//        group2 -= force.getNumDonors();
//        if (group1 < force.getNumAcceptors() && group2 < force.getNumAcceptors()) {
//            force.getAcceptorParameters(group1, p1, p2, p3, params1);
//            force.getAcceptorParameters(group2, p1, p2, p3, params2);
//            return (params1 == params2 && params1 == params2);
//        }
//        if (group1 < force.getNumAcceptors() || group2 < force.getNumAcceptors())
//            return false;
//        return true;
//    }
//private:
//    const CustomHbondForce& force;
//};
//
//CudaCalcCustomHbondForceKernel::~CudaCalcCustomHbondForceKernel() {
//    if (donorParams != NULL)
//        delete donorParams;
//    if (acceptorParams != NULL)
//        delete acceptorParams;
//    if (donors != NULL)
//        delete donors;
//    if (acceptors != NULL)
//        delete acceptors;
//    if (donorBufferIndices != NULL)
//        delete donorBufferIndices;
//    if (acceptorBufferIndices != NULL)
//        delete acceptorBufferIndices;
//    if (globals != NULL)
//        delete globals;
//    if (donorExclusions != NULL)
//        delete donorExclusions;
//    if (acceptorExclusions != NULL)
//        delete acceptorExclusions;
//    if (tabulatedFunctionParams != NULL)
//        delete tabulatedFunctionParams;
//    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//        delete tabulatedFunctions[i];
//}
//
//static void addDonorAndAcceptorCode(stringstream& computeDonor, stringstream& computeAcceptor, const string& value) {
//    computeDonor << value;
//    computeAcceptor << value;
//}
//
//static void applyDonorAndAcceptorForces(stringstream& applyToDonor, stringstream& applyToAcceptor, int atom, const string& value) {
//    string forceNames[] = {"f1", "f2", "f3"};
//    if (atom < 3)
//        applyToAcceptor << forceNames[atom]<<".xyz += "<<value<<";\n";
//    else
//        applyToDonor << forceNames[atom-3]<<".xyz += "<<value<<";\n";
//}
//
//void CudaCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {
//    // Record the lists of donors and acceptors, and the parameters for each one.
//
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumDonors()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumDonors()/numContexts;
//    numDonors = endIndex-startIndex;
//    numAcceptors = force.getNumAcceptors();
//    if (numDonors == 0 || numAcceptors == 0)
//        return;
//    int numParticles = system.getNumParticles();
//    donors = new CudaArray<mm_int4>(cu, numDonors, "customHbondDonors");
//    acceptors = new CudaArray<mm_int4>(cu, numAcceptors, "customHbondAcceptors");
//    donorParams = new CudaParameterSet(cu, force.getNumPerDonorParameters(), numDonors, "customHbondDonorParameters");
//    acceptorParams = new CudaParameterSet(cu, force.getNumPerAcceptorParameters(), numAcceptors, "customHbondAcceptorParameters");
//    if (force.getNumGlobalParameters() > 0)
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customHbondGlobals", false, CL_MEM_READ_ONLY);
//    vector<vector<cl_float> > donorParamVector(numDonors);
//    vector<mm_int4> donorVector(numDonors);
//    for (int i = 0; i < numDonors; i++) {
//        vector<double> parameters;
//        force.getDonorParameters(startIndex+i, donorVector[i].x, donorVector[i].y, donorVector[i].z, parameters);
//        donorParamVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            donorParamVector[i][j] = (cl_float) parameters[j];
//    }
//    donors->upload(donorVector);
//    donorParams->setParameterValues(donorParamVector);
//    vector<vector<cl_float> > acceptorParamVector(numAcceptors);
//    vector<mm_int4> acceptorVector(numAcceptors);
//    for (int i = 0; i < numAcceptors; i++) {
//        vector<double> parameters;
//        force.getAcceptorParameters(i, acceptorVector[i].x, acceptorVector[i].y, acceptorVector[i].z, parameters);
//        acceptorParamVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            acceptorParamVector[i][j] = (cl_float) parameters[j];
//    }
//    acceptors->upload(acceptorVector);
//    acceptorParams->setParameterValues(acceptorParamVector);
//
//    // Select an output buffer index for each donor and acceptor.
//
//    donorBufferIndices = new CudaArray<mm_int4>(cu, numDonors, "customHbondDonorBuffers");
//    acceptorBufferIndices = new CudaArray<mm_int4>(cu, numAcceptors, "customHbondAcceptorBuffers");
//    vector<mm_int4> donorBufferVector(numDonors);
//    vector<mm_int4> acceptorBufferVector(numAcceptors);
//    vector<int> donorBufferCounter(numParticles, 0);
//    for (int i = 0; i < numDonors; i++)
//        donorBufferVector[i] = mm_int4(donorVector[i].x > -1 ? donorBufferCounter[donorVector[i].x]++ : 0,
//                                       donorVector[i].y > -1 ? donorBufferCounter[donorVector[i].y]++ : 0,
//                                       donorVector[i].z > -1 ? donorBufferCounter[donorVector[i].z]++ : 0, 0);
//    vector<int> acceptorBufferCounter(numParticles, 0);
//    for (int i = 0; i < numAcceptors; i++)
//        acceptorBufferVector[i] = mm_int4(acceptorVector[i].x > -1 ? acceptorBufferCounter[acceptorVector[i].x]++ : 0,
//                                       acceptorVector[i].y > -1 ? acceptorBufferCounter[acceptorVector[i].y]++ : 0,
//                                       acceptorVector[i].z > -1 ? acceptorBufferCounter[acceptorVector[i].z]++ : 0, 0);
//    donorBufferIndices->upload(donorBufferVector);
//    acceptorBufferIndices->upload(acceptorBufferVector);
//    int maxBuffers = 1;
//    for (int i = 0; i < (int) donorBufferCounter.size(); i++)
//        maxBuffers = max(maxBuffers, donorBufferCounter[i]);
//    for (int i = 0; i < (int) acceptorBufferCounter.size(); i++)
//        maxBuffers = max(maxBuffers, acceptorBufferCounter[i]);
//    cu.addForce(new CudaCustomHbondForceInfo(maxBuffers, force));
//
//    // Record exclusions.
//
//    vector<mm_int4> donorExclusionVector(numDonors, mm_int4(-1, -1, -1, -1));
//    vector<mm_int4> acceptorExclusionVector(numAcceptors, mm_int4(-1, -1, -1, -1));
//    for (int i = 0; i < force.getNumExclusions(); i++) {
//        int donor, acceptor;
//        force.getExclusionParticles(i, donor, acceptor);
//        if (donor < startIndex || donor >= endIndex)
//            continue;
//        donor -= startIndex;
//        if (donorExclusionVector[donor].x == -1)
//            donorExclusionVector[donor].x = acceptor;
//        else if (donorExclusionVector[donor].y == -1)
//            donorExclusionVector[donor].y = acceptor;
//        else if (donorExclusionVector[donor].z == -1)
//            donorExclusionVector[donor].z = acceptor;
//        else if (donorExclusionVector[donor].w == -1)
//            donorExclusionVector[donor].w = acceptor;
//        else
//            throw OpenMMException("CustomHbondForce: CudaPlatform does not support more than four exclusions per donor");
//        if (acceptorExclusionVector[acceptor].x == -1)
//            acceptorExclusionVector[acceptor].x = donor;
//        else if (acceptorExclusionVector[acceptor].y == -1)
//            acceptorExclusionVector[acceptor].y = donor;
//        else if (acceptorExclusionVector[acceptor].z == -1)
//            acceptorExclusionVector[acceptor].z = donor;
//        else if (acceptorExclusionVector[acceptor].w == -1)
//            acceptorExclusionVector[acceptor].w = donor;
//        else
//            throw OpenMMException("CustomHbondForce: CudaPlatform does not support more than four exclusions per acceptor");
//    }
//    donorExclusions = new CudaArray<mm_int4>(cu, numDonors, "customHbondDonorExclusions");
//    acceptorExclusions = new CudaArray<mm_int4>(cu, numDonors, "customHbondAcceptorExclusions");
//    donorExclusions->upload(donorExclusionVector);
//    acceptorExclusions->upload(acceptorExclusionVector);
//
//    // Record the tabulated functions.
//
//    CudaExpressionUtilities::FunctionPlaceholder fp;
//    map<string, Lepton::CustomFunction*> functions;
//    vector<pair<string, string> > functionDefinitions;
//    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
//    stringstream tableArgs;
//    for (int i = 0; i < force.getNumFunctions(); i++) {
//        string name;
//        vector<double> values;
//        double min, max;
//        force.getFunctionParameters(i, name, values, min, max);
//        string arrayName = "table"+cu.intToString(i);
//        functionDefinitions.push_back(make_pair(name, arrayName));
//        functions[name] = &fp;
//        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), (float) values.size()-2);
//        vector<mm_float4> f = CudaExpressionUtilities::computeFunctionCoefficients(values, min, max);
//        tabulatedFunctions.push_back(new CudaArray<mm_float4>(cu, values.size()-1, "TabulatedFunction"));
//        tabulatedFunctions[tabulatedFunctions.size()-1]->upload(f);
//        tableArgs << ", __global const float4* restrict " << arrayName;
//    }
//    if (force.getNumFunctions() > 0) {
//        tabulatedFunctionParams = new CudaArray<mm_float4>(cu, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
//        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
//        tableArgs << ", __global const float4* restrict functionParams";
//    }
//
//    // Record information about parameters.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    if (globals != NULL)
//        globals->upload(globalParamValues);
//    map<string, string> variables;
//    for (int i = 0; i < force.getNumPerDonorParameters(); i++) {
//        const string& name = force.getPerDonorParameterName(i);
//        variables[name] = "donorParams"+donorParams->getParameterSuffix(i);
//    }
//    for (int i = 0; i < force.getNumPerAcceptorParameters(); i++) {
//        const string& name = force.getPerAcceptorParameterName(i);
//        variables[name] = "acceptorParams"+acceptorParams->getParameterSuffix(i);
//    }
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        const string& name = force.getGlobalParameterName(i);
//        variables[name] = "globals["+cu.intToString(i)+"]";
//    }
//
//    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
//    // and dihedrals the expression depends on.
//
//    map<string, vector<int> > distances;
//    map<string, vector<int> > angles;
//    map<string, vector<int> > dihedrals;
//    Lepton::ParsedExpression energyExpression = CustomHbondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
//    map<string, Lepton::ParsedExpression> forceExpressions;
//    set<string> computedDeltas;
//    computedDeltas.insert("D1A1");
//    string atomNames[] = {"A1", "A2", "A3", "D1", "D2", "D3"};
//    string atomNamesLower[] = {"a1", "a2", "a3", "d1", "d2", "d3"};
//    stringstream computeDonor, computeAcceptor, extraArgs;
//    int index = 0;
//    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
//        if (computedDeltas.count(deltaName) == 0) {
//            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+");\n");
//            computedDeltas.insert(deltaName);
//        }
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float r_"+deltaName+" = sqrt(delta"+deltaName+".w);\n");
//        variables[iter->first] = "r_"+deltaName;
//        forceExpressions["float dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
//        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
//        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
//        if (computedDeltas.count(deltaName1) == 0) {
//            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[0]]+");\n");
//            computedDeltas.insert(deltaName1);
//        }
//        if (computedDeltas.count(deltaName2) == 0) {
//            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[1]]+", "+atomNamesLower[atoms[2]]+");\n");
//            computedDeltas.insert(deltaName2);
//        }
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float "+angleName+" = computeAngle(delta"+deltaName1+", delta"+deltaName2+");\n");
//        variables[iter->first] = angleName;
//        forceExpressions["float dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
//        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
//        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
//        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
//        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
//        string dihedralName = "dihedral_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]]+atomNames[atoms[3]];
//        if (computedDeltas.count(deltaName1) == 0) {
//            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName1+" = delta("+atomNamesLower[atoms[0]]+", "+atomNamesLower[atoms[1]]+");\n");
//            computedDeltas.insert(deltaName1);
//        }
//        if (computedDeltas.count(deltaName2) == 0) {
//            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName2+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[1]]+");\n");
//            computedDeltas.insert(deltaName2);
//        }
//        if (computedDeltas.count(deltaName3) == 0) {
//            addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 delta"+deltaName3+" = delta("+atomNamesLower[atoms[2]]+", "+atomNamesLower[atoms[3]]+");\n");
//            computedDeltas.insert(deltaName3);
//        }
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 "+crossName1+" = computeCross(delta"+deltaName1+", delta"+deltaName2+");\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 "+crossName2+" = computeCross(delta"+deltaName2+", delta"+deltaName3+");\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float "+dihedralName+" = computeAngle("+crossName1+", "+crossName2+");\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, dihedralName+" *= (delta"+deltaName1+".x*"+crossName2+".x + delta"+deltaName1+".y*"+crossName2+".y + delta"+deltaName1+".z*"+crossName2+".z < 0 ? -1 : 1);\n");
//        variables[iter->first] = dihedralName;
//        forceExpressions["float dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
//    }
//
//    // Next it needs to load parameters from global memory.
//
//    if (force.getNumGlobalParameters() > 0)
//        extraArgs << ", __global const float* restrict globals";
//    for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
//        extraArgs << ", __global const "+buffer.getType()+"* restrict donor"+buffer.getName();
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" donorParams"+cu.intToString(i+1)+" = donor"+buffer.getName()+"[index];\n");
//    }
//    for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
//        extraArgs << ", __global const "+buffer.getType()+"* restrict acceptor"+buffer.getName();
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, buffer.getType()+" acceptorParams"+cu.intToString(i+1)+" = acceptor"+buffer.getName()+"[index];\n");
//    }
//
//    // Now evaluate the expressions.
//
//    computeAcceptor << CudaExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, "temp", "functionParams");
//    forceExpressions["energy += "] = energyExpression;
//    computeDonor << CudaExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, "temp", "functionParams");
//
//    // Finally, apply forces to atoms.
//
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
//        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "-"+value);
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], value);
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
//        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 crossProd = cross(delta"+deltaName2+", delta"+deltaName1+");\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float lengthCross = max(length(crossProd), 1e-6f);\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 deltaCross0 = -cross(delta"+deltaName1+", crossProd)*dEdAngle"+cu.intToString(index)+"/(delta"+deltaName1+".w*lengthCross);\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 deltaCross2 = cross(delta"+deltaName2+", crossProd)*dEdAngle"+cu.intToString(index)+"/(delta"+deltaName2+".w*lengthCross);\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 deltaCross1 = -(deltaCross0+deltaCross2);\n");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "deltaCross0.xyz");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "deltaCross1.xyz");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "deltaCross2.xyz");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
//        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
//        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
//        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
//        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "{\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float r = sqrt(delta"+deltaName2+".w);\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 ff;\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.x = (-dEdDihedral"+cu.intToString(index)+"*r)/"+crossName1+".w;\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.y = (delta"+deltaName1+".x*delta"+deltaName2+".x + delta"+deltaName1+".y*delta"+deltaName2+".y + delta"+deltaName1+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.z = (delta"+deltaName3+".x*delta"+deltaName2+".x + delta"+deltaName3+".y*delta"+deltaName2+".y + delta"+deltaName3+".z*delta"+deltaName2+".z)/delta"+deltaName2+".w;\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "ff.w = (dEdDihedral"+cu.intToString(index)+"*r)/"+crossName2+".w;\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 internalF0 = ff.x*"+crossName1+";\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 internalF3 = ff.w*"+crossName2+";\n");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "float4 s = ff.y*internalF0 - ff.z*internalF3;\n");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[0], "internalF0.xyz");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[1], "s.xyz-internalF0.xyz");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[2], "-s.xyz-internalF3.xyz");
//        applyDonorAndAcceptorForces(computeDonor, computeAcceptor, atoms[3], "internalF3.xyz");
//        addDonorAndAcceptorCode(computeDonor, computeAcceptor, "}\n");
//    }
//
//    // Generate the kernels.
//
//    map<string, string> replacements;
//    replacements["COMPUTE_DONOR_FORCE"] = computeDonor.str();
//    replacements["COMPUTE_ACCEPTOR_FORCE"] = computeAcceptor.str();
//    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
//    map<string, string> defines;
//    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//    defines["NUM_DONORS"] = cu.intToString(numDonors);
//    defines["NUM_ACCEPTORS"] = cu.intToString(numAcceptors);
//    defines["M_PI"] = cu.doubleToString(M_PI);
//    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff) {
//        defines["USE_CUTOFF"] = "1";
//        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
//    }
//    if (force.getNonbondedMethod() != CustomHbondForce::NoCutoff && force.getNonbondedMethod() != CustomHbondForce::CutoffNonPeriodic)
//        defines["USE_PERIODIC"] = "1";
//    if (force.getNumExclusions() > 0)
//        defines["USE_EXCLUSIONS"] = "1";
//    cu::Program program = cu.createProgram(cu.replaceStrings(CudaKernelSources::customHbondForce, replacements), defines);
//    donorKernel = cu::Kernel(program, "computeDonorForces");
//    acceptorKernel = cu::Kernel(program, "computeAcceptorForces");
//}
//
//double CudaCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (numDonors == 0 || numAcceptors == 0)
//        return 0.0;
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    if (!hasInitializedKernel) {
//        hasInitializedKernel = true;
//        int index = 0;
//        donorKernel.setArg<cu::Buffer>(index++, cu.getForceBuffers().getDeviceBuffer());
//        donorKernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        donorKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        donorKernel.setArg<cu::Buffer>(index++, donorExclusions->getDeviceBuffer());
//        donorKernel.setArg<cu::Buffer>(index++, donors->getDeviceBuffer());
//        donorKernel.setArg<cu::Buffer>(index++, acceptors->getDeviceBuffer());
//        donorKernel.setArg<cu::Buffer>(index++, donorBufferIndices->getDeviceBuffer());
//        donorKernel.setArg(index++, 3*CudaContext::ThreadBlockSize*sizeof(mm_float4), NULL);
//        index += 2; // Periodic box size arguments are set when the kernel is executed.
//        if (globals != NULL)
//            donorKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//        for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
//            donorKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//        }
//        for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
//            donorKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//        }
//        if (tabulatedFunctionParams != NULL) {
//            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//                donorKernel.setArg<cu::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
//            donorKernel.setArg<cu::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
//        }
//        index = 0;
//        acceptorKernel.setArg<cu::Buffer>(index++, cu.getForceBuffers().getDeviceBuffer());
//        acceptorKernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        acceptorKernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//        acceptorKernel.setArg<cu::Buffer>(index++, acceptorExclusions->getDeviceBuffer());
//        acceptorKernel.setArg<cu::Buffer>(index++, donors->getDeviceBuffer());
//        acceptorKernel.setArg<cu::Buffer>(index++, acceptors->getDeviceBuffer());
//        acceptorKernel.setArg<cu::Buffer>(index++, acceptorBufferIndices->getDeviceBuffer());
//        acceptorKernel.setArg(index++, 3*CudaContext::ThreadBlockSize*sizeof(mm_float4), NULL);
//        index += 2; // Periodic box size arguments are set when the kernel is executed.
//        if (globals != NULL)
//            acceptorKernel.setArg<cu::Buffer>(index++, globals->getDeviceBuffer());
//        for (int i = 0; i < (int) donorParams->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = donorParams->getBuffers()[i];
//            acceptorKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//        }
//        for (int i = 0; i < (int) acceptorParams->getBuffers().size(); i++) {
//            const CudaNonbondedUtilities::ParameterInfo& buffer = acceptorParams->getBuffers()[i];
//            acceptorKernel.setArg<cu::Memory>(index++, buffer.getMemory());
//        }
//        if (tabulatedFunctionParams != NULL) {
//            for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//                acceptorKernel.setArg<cu::Buffer>(index++, tabulatedFunctions[i]->getDeviceBuffer());
//            acceptorKernel.setArg<cu::Buffer>(index++, tabulatedFunctionParams->getDeviceBuffer());
//        }
//    }
//    donorKernel.setArg<mm_float4>(8, cu.getPeriodicBoxSize());
//    donorKernel.setArg<mm_float4>(9, cu.getInvPeriodicBoxSize());
//    cu.executeKernel(donorKernel, max(numDonors, numAcceptors));
//    acceptorKernel.setArg<mm_float4>(8, cu.getPeriodicBoxSize());
//    acceptorKernel.setArg<mm_float4>(9, cu.getInvPeriodicBoxSize());
//    cu.executeKernel(acceptorKernel, max(numDonors, numAcceptors));
//    return 0.0;
//}
//
//void CudaCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumDonors()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumDonors()/numContexts;
//    if (numDonors != endIndex-startIndex)
//        throw OpenMMException("updateParametersInContext: The number of donors has changed");
//    if (numAcceptors != force.getNumAcceptors())
//        throw OpenMMException("updateParametersInContext: The number of acceptors has changed");
//    
//    // Record the per-donor parameters.
//    
//    vector<vector<cl_float> > donorParamVector(numDonors);
//    vector<double> parameters;
//    for (int i = 0; i < numDonors; i++) {
//        int d1, d2, d3;
//        force.getDonorParameters(startIndex+i, d1, d2, d3, parameters);
//        donorParamVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            donorParamVector[i][j] = (cl_float) parameters[j];
//    }
//    donorParams->setParameterValues(donorParamVector);
//    
//    // Record the per-acceptor parameters.
//    
//    vector<vector<cl_float> > acceptorParamVector(numAcceptors);
//    for (int i = 0; i < numAcceptors; i++) {
//        int a1, a2, a3;
//        force.getAcceptorParameters(i, a1, a2, a3, parameters);
//        acceptorParamVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            acceptorParamVector[i][j] = (cl_float) parameters[j];
//    }
//    acceptorParams->setParameterValues(acceptorParamVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}
//
//class CudaCustomCompoundBondForceInfo : public CudaForceInfo {
//public:
//    CudaCustomCompoundBondForceInfo(const CustomCompoundBondForce& force) : CudaForceInfo(0), force(force) {
//    }
//    int getNumParticleGroups() {
//        return force.getNumBonds();
//    }
//    void getParticlesInGroup(int index, vector<int>& particles) {
//        vector<double> parameters;
//        force.getBondParameters(index, particles, parameters);
//    }
//    bool areGroupsIdentical(int group1, int group2) {
//        vector<int> particles;
//        vector<double> parameters1, parameters2;
//        force.getBondParameters(group1, particles, parameters1);
//        force.getBondParameters(group2, particles, parameters2);
//        for (int i = 0; i < (int) parameters1.size(); i++)
//            if (parameters1[i] != parameters2[i])
//                return false;
//        return true;
//    }
//private:
//    const CustomCompoundBondForce& force;
//};
//
//CudaCalcCustomCompoundBondForceKernel::~CudaCalcCustomCompoundBondForceKernel() {
//    if (params != NULL)
//        delete params;
//    if (globals != NULL)
//        delete globals;
//    if (tabulatedFunctionParams != NULL)
//        delete tabulatedFunctionParams;
//    for (int i = 0; i < (int) tabulatedFunctions.size(); i++)
//        delete tabulatedFunctions[i];
//}
//
//void CudaCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
//    numBonds = endIndex-startIndex;
//    if (numBonds == 0)
//        return;
//    int particlesPerBond = force.getNumParticlesPerBond();
//    vector<vector<int> > atoms(numBonds, vector<int>(particlesPerBond));
//    params = new CudaParameterSet(cu, force.getNumPerBondParameters(), numBonds, "customCompoundBondParams");
//    vector<vector<cl_float> > paramVector(numBonds);
//    for (int i = 0; i < numBonds; i++) {
//        vector<double> parameters;
//        force.getBondParameters(startIndex+i, atoms[i], parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    cu.addForce(new CudaCustomCompoundBondForceInfo(force));
//
//    // Record the tabulated functions.
//
//    CudaExpressionUtilities::FunctionPlaceholder fp;
//    map<string, Lepton::CustomFunction*> functions;
//    vector<pair<string, string> > functionDefinitions;
//    vector<mm_float4> tabulatedFunctionParamsVec(force.getNumFunctions());
//    stringstream tableArgs;
//    for (int i = 0; i < force.getNumFunctions(); i++) {
//        string name;
//        vector<double> values;
//        double min, max;
//        force.getFunctionParameters(i, name, values, min, max);
//        functions[name] = &fp;
//        tabulatedFunctionParamsVec[i] = mm_float4((float) min, (float) max, (float) ((values.size()-1)/(max-min)), (float) values.size()-2);
//        vector<mm_float4> f = CudaExpressionUtilities::computeFunctionCoefficients(values, min, max);
//        CudaArray<mm_float4>* array = new CudaArray<mm_float4>(cu, values.size()-1, "TabulatedFunction");
//        tabulatedFunctions.push_back(array);
//        array->upload(f);
//        string arrayName = cu.getBondedUtilities().addArgument(array->getDeviceBuffer(), "float4");
//        functionDefinitions.push_back(make_pair(name, arrayName));
//    }
//    string functionParamsName;
//    if (force.getNumFunctions() > 0) {
//        tabulatedFunctionParams = new CudaArray<mm_float4>(cu, tabulatedFunctionParamsVec.size(), "tabulatedFunctionParameters", false, CL_MEM_READ_ONLY);
//        tabulatedFunctionParams->upload(tabulatedFunctionParamsVec);
//        functionParamsName = cu.getBondedUtilities().addArgument(tabulatedFunctionParams->getDeviceBuffer(), "float4");
//    }
//    
//    // Record information about parameters.
//
//    globalParamNames.resize(force.getNumGlobalParameters());
//    globalParamValues.resize(force.getNumGlobalParameters());
//    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//        globalParamNames[i] = force.getGlobalParameterName(i);
//        globalParamValues[i] = (cl_float) force.getGlobalParameterDefaultValue(i);
//    }
//    map<string, string> variables;
//    for (int i = 0; i < particlesPerBond; i++) {
//        string index = cu.intToString(i+1);
//        variables["x"+index] = "pos"+index+".x";
//        variables["y"+index] = "pos"+index+".y";
//        variables["z"+index] = "pos"+index+".z";
//    }
//    for (int i = 0; i < force.getNumPerBondParameters(); i++) {
//        const string& name = force.getPerBondParameterName(i);
//        variables[name] = "bondParams"+params->getParameterSuffix(i);
//    }
//    if (force.getNumGlobalParameters() > 0) {
//        globals = new CudaArray<cl_float>(cu, force.getNumGlobalParameters(), "customCompoundBondGlobals", false, CL_MEM_READ_ONLY);
//        globals->upload(globalParamValues);
//        string argName = cu.getBondedUtilities().addArgument(globals->getDeviceBuffer(), "float");
//        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
//            const string& name = force.getGlobalParameterName(i);
//            string value = argName+"["+cu.intToString(i)+"]";
//            variables[name] = value;
//        }
//    }
//
//    // Now to generate the kernel.  First, it needs to calculate all distances, angles,
//    // and dihedrals the expression depends on.
//
//    map<string, vector<int> > distances;
//    map<string, vector<int> > angles;
//    map<string, vector<int> > dihedrals;
//    Lepton::ParsedExpression energyExpression = CustomCompoundBondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
//    map<string, Lepton::ParsedExpression> forceExpressions;
//    set<string> computedDeltas;
//    vector<string> atomNames, posNames;
//    for (int i = 0; i < particlesPerBond; i++) {
//        string index = cu.intToString(i+1);
//        atomNames.push_back("P"+index);
//        posNames.push_back("pos"+index);
//    }
//    stringstream compute;
//    int index = 0;
//    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
//        if (computedDeltas.count(deltaName) == 0) {
//            compute<<"float4 delta"<<deltaName<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<");\n";
//            computedDeltas.insert(deltaName);
//        }
//        compute<<"float r_"<<deltaName<<" = sqrt(delta"<<deltaName<<".w);\n";
//        variables[iter->first] = "r_"+deltaName;
//        forceExpressions["float dEdDistance"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
//        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
//        string angleName = "angle_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]];
//        if (computedDeltas.count(deltaName1) == 0) {
//            compute<<"float4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[0]]<<");\n";
//            computedDeltas.insert(deltaName1);
//        }
//        if (computedDeltas.count(deltaName2) == 0) {
//            compute<<"float4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[1]]<<", "<<posNames[atoms[2]]<<");\n";
//            computedDeltas.insert(deltaName2);
//        }
//        compute<<"float "<<angleName<<" = ccb_computeAngle(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
//        variables[iter->first] = angleName;
//        forceExpressions["float dEdAngle"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
//        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
//        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
//        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
//        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
//        string dihedralName = "dihedral_"+atomNames[atoms[0]]+atomNames[atoms[1]]+atomNames[atoms[2]]+atomNames[atoms[3]];
//        if (computedDeltas.count(deltaName1) == 0) {
//            compute<<"float4 delta"<<deltaName1<<" = ccb_delta("<<posNames[atoms[0]]<<", "<<posNames[atoms[1]]<<");\n";
//            computedDeltas.insert(deltaName1);
//        }
//        if (computedDeltas.count(deltaName2) == 0) {
//            compute<<"float4 delta"<<deltaName2<<" = ccb_delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[1]]<<");\n";
//            computedDeltas.insert(deltaName2);
//        }
//        if (computedDeltas.count(deltaName3) == 0) {
//            compute<<"float4 delta"<<deltaName3<<" = ccb_delta("<<posNames[atoms[2]]<<", "<<posNames[atoms[3]]<<");\n";
//            computedDeltas.insert(deltaName3);
//        }
//        compute<<"float4 "<<crossName1<<" = ccb_computeCross(delta"<<deltaName1<<", delta"<<deltaName2<<");\n";
//        compute<<"float4 "<<crossName2<<" = ccb_computeCross(delta"<<deltaName2<<", delta"<<deltaName3<<");\n";
//        compute<<"float "<<dihedralName<<" = ccb_computeAngle("<<crossName1<<", "<<crossName2<<");\n";
//        compute<<dihedralName<<" *= (delta"<<deltaName1<<".x*"<<crossName2<<".x + delta"<<deltaName1<<".y*"<<crossName2<<".y + delta"<<deltaName1<<".z*"<<crossName2<<".z < 0 ? -1 : 1);\n";
//        variables[iter->first] = dihedralName;
//        forceExpressions["float dEdDihedral"+cu.intToString(index)+" = "] = energyExpression.differentiate(iter->first).optimize();
//    }
//
//    // Now evaluate the expressions.
//
//    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
//        const CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
//        string argName = cu.getBondedUtilities().addArgument(buffer.getMemory(), buffer.getType());
//        compute<<buffer.getType()<<" bondParams"<<(i+1)<<" = "<<argName<<"[index];\n";
//    }
//    forceExpressions["energy += "] = energyExpression;
//    compute << CudaExpressionUtilities::createExpressions(forceExpressions, variables, functionDefinitions, "temp", functionParamsName);
//
//    // Finally, apply forces to atoms.
//
//    vector<string> forceNames;
//    for (int i = 0; i < particlesPerBond; i++) {
//        string istr = cu.intToString(i+1);
//        string forceName = "force"+istr;
//        forceNames.push_back(forceName);
//        compute<<"float4 "<<forceName<<" = (float4) (0.0f, 0.0f, 0.0f, 0.0f);\n";
//        compute<<"{\n";
//        Lepton::ParsedExpression forceExpressionX = energyExpression.differentiate("x"+istr).optimize();
//        Lepton::ParsedExpression forceExpressionY = energyExpression.differentiate("y"+istr).optimize();
//        Lepton::ParsedExpression forceExpressionZ = energyExpression.differentiate("z"+istr).optimize();
//        map<string, Lepton::ParsedExpression> expressions;
//        if (!isZeroExpression(forceExpressionX))
//            expressions[forceName+".x -= "] = forceExpressionX;
//        if (!isZeroExpression(forceExpressionY))
//            expressions[forceName+".y -= "] = forceExpressionY;
//        if (!isZeroExpression(forceExpressionZ))
//            expressions[forceName+".z -= "] = forceExpressionZ;
//        if (expressions.size() > 0)
//            compute<<CudaExpressionUtilities::createExpressions(expressions, variables, functionDefinitions, "coordtemp", functionParamsName);
//        compute<<"}\n";
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = distances.begin(); iter != distances.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName = atomNames[atoms[0]]+atomNames[atoms[1]];
//        string value = "(dEdDistance"+cu.intToString(index)+"/r_"+deltaName+")*delta"+deltaName+".xyz";
//        compute<<forceNames[atoms[0]]<<".xyz += "<<"-"<<value<<";\n";
//        compute<<forceNames[atoms[1]]<<".xyz += "<<value<<";\n";
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = angles.begin(); iter != angles.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[1]]+atomNames[atoms[0]];
//        string deltaName2 = atomNames[atoms[1]]+atomNames[atoms[2]];
//        compute<<"{\n";
//        compute<<"float4 crossProd = cross(delta"<<deltaName2<<", delta"<<deltaName1<<");\n";
//        compute<<"float lengthCross = max(length(crossProd), 1e-6f);\n";
//        compute<<"float4 deltaCross0 = -cross(delta"<<deltaName1<<", crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName1<<".w*lengthCross);\n";
//        compute<<"float4 deltaCross2 = cross(delta"<<deltaName2<<", crossProd)*dEdAngle"<<cu.intToString(index)<<"/(delta"<<deltaName2<<".w*lengthCross);\n";
//        compute<<"float4 deltaCross1 = -(deltaCross0+deltaCross2);\n";
//        compute<<forceNames[atoms[0]]<<".xyz += deltaCross0.xyz;\n";
//        compute<<forceNames[atoms[1]]<<".xyz += deltaCross1.xyz;\n";
//        compute<<forceNames[atoms[2]]<<".xyz += deltaCross2.xyz;\n";
//        compute<<"}\n";
//    }
//    index = 0;
//    for (map<string, vector<int> >::const_iterator iter = dihedrals.begin(); iter != dihedrals.end(); ++iter, ++index) {
//        const vector<int>& atoms = iter->second;
//        string deltaName1 = atomNames[atoms[0]]+atomNames[atoms[1]];
//        string deltaName2 = atomNames[atoms[2]]+atomNames[atoms[1]];
//        string deltaName3 = atomNames[atoms[2]]+atomNames[atoms[3]];
//        string crossName1 = "cross_"+deltaName1+"_"+deltaName2;
//        string crossName2 = "cross_"+deltaName2+"_"+deltaName3;
//        compute<<"{\n";
//        compute<<"float r = sqrt(delta"<<deltaName2<<".w);\n";
//        compute<<"float4 ff;\n";
//        compute<<"ff.x = (-dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName1<<".w;\n";
//        compute<<"ff.y = (delta"<<deltaName1<<".x*delta"<<deltaName2<<".x + delta"<<deltaName1<<".y*delta"<<deltaName2<<".y + delta"<<deltaName1<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
//        compute<<"ff.z = (delta"<<deltaName3<<".x*delta"<<deltaName2<<".x + delta"<<deltaName3<<".y*delta"<<deltaName2<<".y + delta"<<deltaName3<<".z*delta"<<deltaName2<<".z)/delta"<<deltaName2<<".w;\n";
//        compute<<"ff.w = (dEdDihedral"<<cu.intToString(index)<<"*r)/"<<crossName2<<".w;\n";
//        compute<<"float4 internalF0 = ff.x*"<<crossName1<<";\n";
//        compute<<"float4 internalF3 = ff.w*"<<crossName2<<";\n";
//        compute<<"float4 s = ff.y*internalF0 - ff.z*internalF3;\n";
//        compute<<forceNames[atoms[0]]<<".xyz += internalF0.xyz;\n";
//        compute<<forceNames[atoms[1]]<<".xyz += s.xyz-internalF0.xyz;\n";
//        compute<<forceNames[atoms[2]]<<".xyz += -s.xyz-internalF3.xyz;\n";
//        compute<<forceNames[atoms[3]]<<".xyz += internalF3.xyz;\n";
//        compute<<"}\n";
//    }
//    cu.getBondedUtilities().addInteraction(atoms, compute.str(), force.getForceGroup());
//    map<string, string> replacements;
//    replacements["M_PI"] = cu.doubleToString(M_PI);
//    cu.getBondedUtilities().addPrefixCode(cu.replaceStrings(CudaKernelSources::customCompoundBond, replacements));;
//}
//
//double CudaCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
//    if (globals != NULL) {
//        bool changed = false;
//        for (int i = 0; i < (int) globalParamNames.size(); i++) {
//            cl_float value = (cl_float) context.getParameter(globalParamNames[i]);
//            if (value != globalParamValues[i])
//                changed = true;
//            globalParamValues[i] = value;
//        }
//        if (changed)
//            globals->upload(globalParamValues);
//    }
//    return 0.0;
//}
//
//void CudaCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
//    int numContexts = cu.getPlatformData().contexts.size();
//    int startIndex = cu.getContextIndex()*force.getNumBonds()/numContexts;
//    int endIndex = (cu.getContextIndex()+1)*force.getNumBonds()/numContexts;
//    if (numBonds != endIndex-startIndex)
//        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
//    
//    // Record the per-bond parameters.
//    
//    vector<vector<cl_float> > paramVector(numBonds);
//    vector<int> particles;
//    vector<double> parameters;
//    for (int i = 0; i < numBonds; i++) {
//        force.getBondParameters(startIndex+i, particles, parameters);
//        paramVector[i].resize(parameters.size());
//        for (int j = 0; j < (int) parameters.size(); j++)
//            paramVector[i][j] = (cl_float) parameters[j];
//    }
//    params->setParameterValues(paramVector);
//    
//    // Mark that the current reordering may be invalid.
//    
//    cu.invalidateMolecules();
//}

CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
}

void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    map<string, string> defines;
    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(CudaKernelSources::verlet, defines, "");
    kernel1 = cu.getKernel(module, "integrateVerletPart1");
    kernel2 = cu.getKernel(module, "integrateVerletPart2");
    prevStepSize = -1.0;
}

void CudaIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    double dt = integrator.getStepSize();
    if (dt != prevStepSize) {
        if (cu.getUseDoublePrecision()) {
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

    void* args1[] = {&cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(),
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
}

//CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
//    if (params != NULL)
//        delete params;
//}
//
//void CudaIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
//    cu.getPlatformData().initializeContexts(system);
//    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
//    map<string, string> defines;
//    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//    cu::Program program = cu.createProgram(CudaKernelSources::langevin, defines, "");
//    kernel1 = cu::Kernel(program, "integrateLangevinPart1");
//    kernel2 = cu::Kernel(program, "integrateLangevinPart2");
//    params = new CudaArray<cl_float>(cu, 3, "langevinParams");
//    prevStepSize = -1.0;
//}
//
//void CudaIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
//    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
//    int numAtoms = cu.getNumAtoms();
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        kernel1.setArg<cu::Buffer>(0, cu.getVelm().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(1, cu.getForce().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(3, params->getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(4, integration.getStepSize().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(5, integration.getRandom().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(0, cu.getPosq().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(1, integration.getPosDelta().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(2, cu.getVelm().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(3, integration.getStepSize().getDeviceBuffer());
//    }
//    double temperature = integrator.getTemperature();
//    double friction = integrator.getFriction();
//    double stepSize = integrator.getStepSize();
//    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
//        // Calculate the integration parameters.
//
//        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
//        double kT = BOLTZ*temperature;
//        double vscale = exp(-stepSize/tau);
//        double fscale = (1-vscale)*tau;
//        double noisescale = sqrt(2*kT/tau)*sqrt(0.5*(1-vscale*vscale)*tau);
//        vector<cl_float> p(params->getSize());
//        p[0] = (cl_float) vscale;
//        p[1] = (cl_float) fscale;
//        p[2] = (cl_float) noisescale;
//        params->upload(p);
//        integration.getStepSize()[0].y = (cl_float) stepSize;
//        integration.getStepSize().upload();
//        prevTemp = temperature;
//        prevFriction = friction;
//        prevStepSize = stepSize;
//    }
//
//    // Call the first integration kernel.
//
//    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cu.getPaddedNumAtoms()));
//    cu.executeKernel(kernel1, numAtoms);
//
//    // Apply constraints.
//
//    integration.applyConstraints(integrator.getConstraintTolerance());
//
//    // Call the second integration kernel.
//
//    cu.executeKernel(kernel2, numAtoms);
//    integration.computeVirtualSites();
//
//    // Update the time and step count.
//
//    cu.setTime(cu.getTime()+stepSize);
//    cu.setStepCount(cu.getStepCount()+1);
//}
//
//CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
//}
//
//void CudaIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
//    cu.getPlatformData().initializeContexts(system);
//    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
//    map<string, string> defines;
//    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//    cu::Program program = cu.createProgram(CudaKernelSources::brownian, defines, "");
//    kernel1 = cu::Kernel(program, "integrateBrownianPart1");
//    kernel2 = cu::Kernel(program, "integrateBrownianPart2");
//    prevStepSize = -1.0;
//}
//
//void CudaIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
//    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
//    int numAtoms = cu.getNumAtoms();
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        kernel1.setArg<cu::Buffer>(2, cu.getForce().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(3, integration.getPosDelta().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(4, cu.getVelm().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(5, integration.getRandom().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(1, cu.getPosq().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(2, cu.getVelm().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(3, integration.getPosDelta().getDeviceBuffer());
//    }
//    double temperature = integrator.getTemperature();
//    double friction = integrator.getFriction();
//    double stepSize = integrator.getStepSize();
//    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
//        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
//        kernel1.setArg<cl_float>(0, (cl_float) (tau*stepSize));
//        kernel1.setArg<cl_float>(1, (cl_float) (sqrt(2.0f*BOLTZ*temperature*stepSize*tau)));
//        kernel2.setArg<cl_float>(0, (cl_float) (1.0/stepSize));
//        prevTemp = temperature;
//        prevFriction = friction;
//        prevStepSize = stepSize;
//    }
//
//    // Call the first integration kernel.
//
//    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cu.getPaddedNumAtoms()));
//    cu.executeKernel(kernel1, numAtoms);
//
//    // Apply constraints.
//
//    integration.applyConstraints(integrator.getConstraintTolerance());
//
//    // Call the second integration kernel.
//
//    cu.executeKernel(kernel2, numAtoms);
//    integration.computeVirtualSites();
//
//    // Update the time and step count.
//
//    cu.setTime(cu.getTime()+stepSize);
//    cu.setStepCount(cu.getStepCount()+1);
//}
//
//CudaIntegrateVariableVerletStepKernel::~CudaIntegrateVariableVerletStepKernel() {
//}
//
//void CudaIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
//    cu.getPlatformData().initializeContexts(system);
//    cu::Program program = cu.createProgram(CudaKernelSources::verlet, "");
//    kernel1 = cu::Kernel(program, "integrateVerletPart1");
//    kernel2 = cu::Kernel(program, "integrateVerletPart2");
//    selectSizeKernel = cu::Kernel(program, "selectVerletStepSize");
//    blockSize = min(min(256, system.getNumParticles()), (int) cu.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
//}
//
//double CudaIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
//    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
//    int numAtoms = cu.getNumAtoms();
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        kernel1.setArg<cl_int>(0, numAtoms);
//        kernel1.setArg<cu::Buffer>(1, cu.getIntegrationUtilities().getStepSize().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(2, cu.getPosq().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(3, cu.getVelm().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(4, cu.getForce().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
//        kernel2.setArg<cl_int>(0, numAtoms);
//        kernel2.setArg<cu::Buffer>(1, cu.getIntegrationUtilities().getStepSize().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(2, cu.getPosq().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(3, cu.getVelm().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(4, integration.getPosDelta().getDeviceBuffer());
//        selectSizeKernel.setArg<cl_int>(0, numAtoms);
//        selectSizeKernel.setArg<cu::Buffer>(3, cu.getIntegrationUtilities().getStepSize().getDeviceBuffer());
//        selectSizeKernel.setArg<cu::Buffer>(4, cu.getVelm().getDeviceBuffer());
//        selectSizeKernel.setArg<cu::Buffer>(5, cu.getForce().getDeviceBuffer());
//        selectSizeKernel.setArg(6, blockSize*sizeof(cl_float), NULL);
//    }
//
//    // Select the step size to use.
//
//    float maxStepSize = (float)(maxTime-cu.getTime());
//    selectSizeKernel.setArg<cl_float>(1, maxStepSize);
//    selectSizeKernel.setArg<cl_float>(2, (cl_float) integrator.getErrorTolerance());
//    cu.executeKernel(selectSizeKernel, blockSize, blockSize);
//
//    // Call the first integration kernel.
//
//    cu.executeKernel(kernel1, numAtoms);
//
//    // Apply constraints.
//
//    integration.applyConstraints(integrator.getConstraintTolerance());
//
//    // Call the second integration kernel.
//
//    cu.executeKernel(kernel2, numAtoms);
//    integration.computeVirtualSites();
//
//    // Update the time and step count.
//
//    cu.getIntegrationUtilities().getStepSize().download();
//    double dt = cu.getIntegrationUtilities().getStepSize()[0].y;
//    double time = cu.getTime()+dt;
//    if (dt == maxStepSize)
//        time = maxTime; // Avoid round-off error
//    cu.setTime(time);
//    cu.setStepCount(cu.getStepCount()+1);
//    return dt;
//}
//
//CudaIntegrateVariableLangevinStepKernel::~CudaIntegrateVariableLangevinStepKernel() {
//    if (params != NULL)
//        delete params;
//}
//
//void CudaIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
//    cu.getPlatformData().initializeContexts(system);
//    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
//    map<string, string> defines;
//    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
//    cu::Program program = cu.createProgram(CudaKernelSources::langevin, defines, "");
//    kernel1 = cu::Kernel(program, "integrateLangevinPart1");
//    kernel2 = cu::Kernel(program, "integrateLangevinPart2");
//    selectSizeKernel = cu::Kernel(program, "selectLangevinStepSize");
//    params = new CudaArray<cl_float>(cu, 3, "langevinParams");
//    blockSize = min(256, system.getNumParticles());
//    blockSize = max(blockSize, params->getSize());
//    blockSize = min(blockSize, (int) cu.getDevice().getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>());
//}
//
//double CudaIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
//    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
//    int numAtoms = cu.getNumAtoms();
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        kernel1.setArg<cu::Buffer>(0, cu.getVelm().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(1, cu.getForce().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(2, integration.getPosDelta().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(3, params->getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(4, integration.getStepSize().getDeviceBuffer());
//        kernel1.setArg<cu::Buffer>(5, integration.getRandom().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(0, cu.getPosq().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(1, integration.getPosDelta().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(2, cu.getVelm().getDeviceBuffer());
//        kernel2.setArg<cu::Buffer>(3, integration.getStepSize().getDeviceBuffer());
//        selectSizeKernel.setArg<cu::Buffer>(4, integration.getStepSize().getDeviceBuffer());
//        selectSizeKernel.setArg<cu::Buffer>(5, cu.getVelm().getDeviceBuffer());
//        selectSizeKernel.setArg<cu::Buffer>(6, cu.getForce().getDeviceBuffer());
//        selectSizeKernel.setArg<cu::Buffer>(7, params->getDeviceBuffer());
//        selectSizeKernel.setArg(8, params->getSize()*sizeof(cl_float), NULL);
//        selectSizeKernel.setArg(9, blockSize*sizeof(cl_float), NULL);
//    }
//
//    // Select the step size to use.
//
//    float maxStepSize = (float)(maxTime-cu.getTime());
//    selectSizeKernel.setArg<cl_float>(0, maxStepSize);
//    selectSizeKernel.setArg<cl_float>(1, (cl_float) integrator.getErrorTolerance());
//    selectSizeKernel.setArg<cl_float>(2, (cl_float) (integrator.getFriction() == 0.0 ? 0.0 : 1.0/integrator.getFriction()));
//    selectSizeKernel.setArg<cl_float>(3, (cl_float) (BOLTZ*integrator.getTemperature()));
//    cu.executeKernel(selectSizeKernel, blockSize, blockSize);
//
//    // Call the first integration kernel.
//
//    kernel1.setArg<cl_uint>(6, integration.prepareRandomNumbers(cu.getPaddedNumAtoms()));
//    cu.executeKernel(kernel1, numAtoms);
//
//    // Apply constraints.
//
//    integration.applyConstraints(integrator.getConstraintTolerance());
//
//    // Call the second integration kernel.
//
//    cu.executeKernel(kernel2, numAtoms);
//    integration.computeVirtualSites();
//
//    // Update the time and step count.
//
//    cu.getIntegrationUtilities().getStepSize().download();
//    double dt = cu.getIntegrationUtilities().getStepSize()[0].y;
//    double time = cu.getTime()+dt;
//    if (dt == maxStepSize)
//        time = maxTime; // Avoid round-off error
//    cu.setTime(time);
//    cu.setStepCount(cu.getStepCount()+1);
//    return dt;
//}
//
//class CudaIntegrateCustomStepKernel::ReorderListener : public CudaContext::ReorderListener {
//public:
//    ReorderListener(CudaContext& cu, CudaParameterSet& perDofValues, vector<vector<cl_float> >& localPerDofValues, bool& deviceValuesAreCurrent) :
//            cu(cu), perDofValues(perDofValues), localPerDofValues(localPerDofValues), deviceValuesAreCurrent(deviceValuesAreCurrent) {
//        int numAtoms = cu.getNumAtoms();
//        lastAtomOrder.resize(numAtoms);
//        for (int i = 0; i < numAtoms; i++)
//            lastAtomOrder[i] = cu.getAtomIndex()[i];
//    }
//    void execute() {
//        // Reorder the per-DOF variables to reflect the new atom order.
//
//        if (perDofValues.getNumParameters() == 0)
//            return;
//        int numAtoms = cu.getNumAtoms();
//        if (deviceValuesAreCurrent)
//            perDofValues.getParameterValues(localPerDofValues);
//        vector<vector<cl_float> > swap(3*numAtoms);
//        for (int i = 0; i < numAtoms; i++) {
//            swap[3*lastAtomOrder[i]] = localPerDofValues[3*i];
//            swap[3*lastAtomOrder[i]+1] = localPerDofValues[3*i+1];
//            swap[3*lastAtomOrder[i]+2] = localPerDofValues[3*i+2];
//        }
//        CudaArray<cl_int>& order = cu.getAtomIndex();
//        for (int i = 0; i < numAtoms; i++) {
//            localPerDofValues[3*i] = swap[3*order[i]];
//            localPerDofValues[3*i+1] = swap[3*order[i]+1];
//            localPerDofValues[3*i+2] = swap[3*order[i]+2];
//        }
//        perDofValues.setParameterValues(localPerDofValues);
//        for (int i = 0; i < numAtoms; i++)
//            lastAtomOrder[i] = order[i];
//        deviceValuesAreCurrent = true;
//    }
//private:
//    CudaContext& cu;
//    CudaParameterSet& perDofValues;
//    vector<vector<cl_float> >& localPerDofValues;
//    bool& deviceValuesAreCurrent;
//    vector<int> lastAtomOrder;
//};
//
//CudaIntegrateCustomStepKernel::~CudaIntegrateCustomStepKernel() {
//    if (globalValues != NULL)
//        delete globalValues;
//    if (contextParameterValues != NULL)
//        delete contextParameterValues;
//    if (sumBuffer != NULL)
//        delete sumBuffer;
//    if (energy != NULL)
//        delete energy;
//    if (uniformRandoms != NULL)
//        delete uniformRandoms;
//    if (randomSeed != NULL)
//        delete randomSeed;
//    if (perDofValues != NULL)
//        delete perDofValues;
//}
//
//void CudaIntegrateCustomStepKernel::initialize(const System& system, const CustomIntegrator& integrator) {
//    cu.getPlatformData().initializeContexts(system);
//    cu.getIntegrationUtilities().initRandomNumberGenerator(integrator.getRandomNumberSeed());
//    numGlobalVariables = integrator.getNumGlobalVariables();
//    globalValues = new CudaArray<cl_float>(cu, max(1, numGlobalVariables), "globalVariables", true);
//    sumBuffer = new CudaArray<cl_float>(cu, 3*system.getNumParticles(), "sumBuffer");
//    energy = new CudaArray<cl_float>(cu, 1, "energy");
//    perDofValues = new CudaParameterSet(cu, integrator.getNumPerDofVariables(), 3*system.getNumParticles(), "perDofVariables");
//    cu.addReorderListener(new ReorderListener(cu, *perDofValues, localPerDofValues, deviceValuesAreCurrent));
//    prevStepSize = -1.0;
//    SimTKOpenMMUtilities::setRandomNumberSeed(integrator.getRandomNumberSeed());
//}
//
//string CudaIntegrateCustomStepKernel::createGlobalComputation(const string& variable, const Lepton::ParsedExpression& expr, CustomIntegrator& integrator, const string& energyName) {
//    map<string, Lepton::ParsedExpression> expressions;
//    if (variable == "dt")
//        expressions["dt[0].y = "] = expr;
//    else {
//        for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
//            if (variable == integrator.getGlobalVariableName(i))
//                expressions["globals["+cu.intToString(i)+"] = "] = expr;
//        for (int i = 0; i < (int) parameterNames.size(); i++)
//            if (variable == parameterNames[i]) {
//                expressions["params["+cu.intToString(i)+"] = "] = expr;
//                modifiesParameters = true;
//            }
//    }
//    if (expressions.size() == 0)
//        throw OpenMMException("Unknown global variable: "+variable);
//    map<string, string> variables;
//    variables["dt"] = "dt[0].y";
//    variables["uniform"] = "uniform";
//    variables["gaussian"] = "gaussian";
//    variables[energyName] = "energy[0]";
//    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
//        variables[integrator.getGlobalVariableName(i)] = "globals["+cu.intToString(i)+"]";
//    for (int i = 0; i < (int) parameterNames.size(); i++)
//        variables[parameterNames[i]] = "params["+cu.intToString(i)+"]";
//    vector<pair<string, string> > functions;
//    return CudaExpressionUtilities::createExpressions(expressions, variables, functions, "temp", "");
//}
//
//string CudaIntegrateCustomStepKernel::createPerDofComputation(const string& variable, const Lepton::ParsedExpression& expr, int component, CustomIntegrator& integrator, const string& forceName, const string& energyName) {
//    const string suffixes[] = {".x", ".y", ".z"};
//    string suffix = suffixes[component];
//    map<string, Lepton::ParsedExpression> expressions;
//    if (variable == "x")
//        expressions["position"+suffix+" = "] = expr;
//    else if (variable == "v")
//        expressions["velocity"+suffix+" = "] = expr;
//    else if (variable == "")
//        expressions["sum[3*index+"+cu.intToString(component)+"] = "] = expr;
//    else {
//        for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
//            if (variable == integrator.getPerDofVariableName(i))
//                expressions["perDof"+suffix.substr(1)+perDofValues->getParameterSuffix(i)+" = "] = expr;
//    }
//    if (expressions.size() == 0)
//        throw OpenMMException("Unknown per-DOF variable: "+variable);
//    map<string, string> variables;
//    variables["x"] = "position"+suffix;
//    variables["v"] = "velocity"+suffix;
//    variables[forceName] = "f"+suffix;
//    variables["gaussian"] = "gaussian"+suffix;
//    variables["uniform"] = "uniform"+suffix;
//    variables["m"] = "mass";
//    variables["dt"] = "stepSize";
//    variables[energyName] = "energy[0]";
//    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
//        variables[integrator.getGlobalVariableName(i)] = "globals["+cu.intToString(i)+"]";
//    for (int i = 0; i < integrator.getNumPerDofVariables(); i++)
//        variables[integrator.getPerDofVariableName(i)] = "perDof"+suffix.substr(1)+perDofValues->getParameterSuffix(i);
//    for (int i = 0; i < (int) parameterNames.size(); i++)
//        variables[parameterNames[i]] = "params["+cu.intToString(i)+"]";
//    vector<pair<string, string> > functions;
//    string tempType = (cu.getSupportsDoublePrecision() ? "double" : "float");
//    return CudaExpressionUtilities::createExpressions(expressions, variables, functions, "temp"+cu.intToString(component)+"_", "", tempType);
//}
//
//void CudaIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
//    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
//    int numAtoms = cu.getNumAtoms();
//    int numSteps = integrator.getNumComputations();
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        
//        // Initialize various data structures.
//        
//        const map<string, double>& params = context.getParameters();
//        contextParameterValues = new CudaArray<cl_float>(cu, max(1, (int) params.size()), "contextParameters", true);
//        for (map<string, double>::const_iterator iter = params.begin(); iter != params.end(); ++iter) {
//            contextParameterValues->set(parameterNames.size(), (float) iter->second);
//            parameterNames.push_back(iter->first);
//        }
//        contextParameterValues->upload();
//        kernels.resize(integrator.getNumComputations());
//        requiredGaussian.resize(integrator.getNumComputations(), 0);
//        requiredUniform.resize(integrator.getNumComputations(), 0);
//        needsForces.resize(numSteps, false);
//        needsEnergy.resize(numSteps, false);
//        forceGroup.resize(numSteps, -2);
//        invalidatesForces.resize(numSteps, false);
//        merged.resize(numSteps, false);
//        modifiesParameters = false;
//        map<string, string> defines;
//        defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//        defines["WORK_GROUP_SIZE"] = cu.intToString(CudaContext::ThreadBlockSize);
//        
//        // Initialize the random number generator.
//        
//        uniformRandoms = new CudaArray<mm_float4>(cu, cu.getNumAtoms(), "uniformRandoms");
//        randomSeed = new CudaArray<mm_int4>(cu, cu.getNumThreadBlocks()*CudaContext::ThreadBlockSize, "randomSeed");
//        vector<mm_int4> seed(randomSeed->getSize());
//        unsigned int r = integrator.getRandomNumberSeed()+1;
//        for (int i = 0; i < randomSeed->getSize(); i++) {
//            seed[i].x = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
//            seed[i].y = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
//            seed[i].z = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
//            seed[i].w = r = (1664525*r + 1013904223) & 0xFFFFFFFF;
//        }
//        randomSeed->upload(seed);
//        cu::Program randomProgram = cu.createProgram(CudaKernelSources::customIntegrator, defines);
//        randomKernel = cu::Kernel(randomProgram, "generateRandomNumbers");
//        randomKernel.setArg<cu::Buffer>(0, uniformRandoms->getDeviceBuffer());
//        randomKernel.setArg<cu::Buffer>(1, randomSeed->getDeviceBuffer());
//        
//        // Build a list of all variables that affect the forces, so we can tell which
//        // steps invalidate them.
//        
//        set<string> affectsForce;
//        affectsForce.insert("x");
//        for (vector<ForceImpl*>::const_iterator iter = context.getForceImpls().begin(); iter != context.getForceImpls().end(); ++iter) {
//            const map<string, double> params = (*iter)->getDefaultParameters();
//            for (map<string, double>::const_iterator param = params.begin(); param != params.end(); ++param)
//                affectsForce.insert(param->first);
//        }
//        
//        // Record information about all the computation steps.
//        
//        stepType.resize(numSteps);
//        vector<string> variable(numSteps);
//        vector<Lepton::ParsedExpression> expression(numSteps);
//        vector<string> forceGroupName;
//        vector<string> energyGroupName;
//        for (int i = 0; i < 32; i++) {
//            stringstream fname;
//            fname << "f" << i;
//            forceGroupName.push_back(fname.str());
//            stringstream ename;
//            ename << "energy" << i;
//            energyGroupName.push_back(ename.str());
//        }
//        vector<string> forceName(numSteps, "f");
//        vector<string> energyName(numSteps, "energy");
//        for (int step = 0; step < numSteps; step++) {
//            string expr;
//            integrator.getComputationStep(step, stepType[step], variable[step], expr);
//            if (expr.size() > 0) {
//                expression[step] = Lepton::Parser::parse(expr).optimize();
//                if (usesVariable(expression[step], "f")) {
//                    needsForces[step] = true;
//                    forceGroup[step] = -1;
//                }
//                if (usesVariable(expression[step], "energy")) {
//                    needsEnergy[step] = true;
//                    forceGroup[step] = -1;
//                }
//                for (int i = 0; i < 32; i++) {
//                    if (usesVariable(expression[step], forceGroupName[i])) {
//                        if (forceGroup[step] != -2)
//                            throw OpenMMException("A single computation step cannot depend on multiple force groups");
//                        needsForces[step] = true;
//                        forceGroup[step] = 1<<i;
//                        forceName[step] = forceGroupName[i];
//                    }
//                    if (usesVariable(expression[step], energyGroupName[i])) {
//                        if (forceGroup[step] != -2)
//                            throw OpenMMException("A single computation step cannot depend on multiple force groups");
//                        needsEnergy[step] = true;
//                        forceGroup[step] = 1<<i;
//                        energyName[step] = energyGroupName[i];
//                    }
//                }
//            }
//            invalidatesForces[step] = (stepType[step] == CustomIntegrator::ConstrainPositions || affectsForce.find(variable[step]) != affectsForce.end());
//            if (forceGroup[step] == -2 && step > 0)
//                forceGroup[step] = forceGroup[step-1];
//        }
//        
//        // Determine how each step will represent the position (as just a value, or a value plus a delta).
//        
//        vector<bool> storePosAsDelta(numSteps, false);
//        vector<bool> loadPosAsDelta(numSteps, false);
//        bool beforeConstrain = false;
//        for (int step = numSteps-1; step >= 0; step--) {
//            if (stepType[step] == CustomIntegrator::ConstrainPositions)
//                beforeConstrain = true;
//            else if (stepType[step] == CustomIntegrator::ComputePerDof && variable[step] == "x" && beforeConstrain)
//                storePosAsDelta[step] = true;
//        }
//        bool storedAsDelta = false;
//        for (int step = 0; step < numSteps; step++) {
//            loadPosAsDelta[step] = storedAsDelta;
//            if (storePosAsDelta[step] == true)
//                storedAsDelta = true;
//            if (stepType[step] == CustomIntegrator::ConstrainPositions)
//                storedAsDelta = false;
//        }
//        
//        // Identify steps that can be merged into a single kernel.
//        
//        for (int step = 1; step < numSteps; step++) {
//            if (needsForces[step] || needsEnergy[step])
//                continue;
//            if (stepType[step-1] == CustomIntegrator::ComputeGlobal && stepType[step] == CustomIntegrator::ComputeGlobal)
//                merged[step] = true;
//            if (stepType[step-1] == CustomIntegrator::ComputePerDof && stepType[step] == CustomIntegrator::ComputePerDof &&
//                    !usesVariable(expression[step], "uniform"))
//                merged[step] = true;
//        }
//        
//        // Loop over all steps and create the kernels for them.
//        
//        for (int step = 0; step < numSteps; step++) {
//            if ((stepType[step] == CustomIntegrator::ComputePerDof || stepType[step] == CustomIntegrator::ComputeSum) && !merged[step]) {
//                // Compute a per-DOF value.
//                
//                stringstream compute;
//                for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
//                    const CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
//                    compute << buffer.getType()<<" perDofx"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index];\n";
//                    compute << buffer.getType()<<" perDofy"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index+1];\n";
//                    compute << buffer.getType()<<" perDofz"<<cu.intToString(i+1)<<" = perDofValues"<<cu.intToString(i+1)<<"[3*index+2];\n";
//                }
//                string convert = (cu.getSupportsDoublePrecision() ? "convert_float4(" : "(");
//                int numGaussian = 0, numUniform = 0;
//                for (int j = step; j < numSteps && (j == step || merged[j]); j++) {
//                    compute << "{\n";
//                    for (int i = 0; i < 3; i++)
//                        compute << createPerDofComputation(stepType[j] == CustomIntegrator::ComputePerDof ? variable[j] : "", expression[j], i, integrator, forceName[j], energyName[j]);
//                    if (variable[j] == "x") {
//                        if (storePosAsDelta[j]) {
//                            if (cu.getSupportsDoublePrecision())
//                                compute << "posDelta[index] = convert_float4(position-convert_double4(posq[index]));\n";
//                            else
//                                compute << "posDelta[index] = position-posq[index];\n";
//                        }
//                        else
//                            compute << "posq[index] = " << convert << "position);\n";
//                    }
//                    else if (variable[j] == "v")
//                        compute << "velm[index] = " << convert << "velocity);\n";
//                    else {
//                        for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
//                            const CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
//                            compute << "perDofValues"<<cu.intToString(i+1)<<"[3*index] = perDofx"<<cu.intToString(i+1)<<";\n";
//                            compute << "perDofValues"<<cu.intToString(i+1)<<"[3*index+1] = perDofy"<<cu.intToString(i+1)<<";\n";
//                            compute << "perDofValues"<<cu.intToString(i+1)<<"[3*index+2] = perDofz"<<cu.intToString(i+1)<<";\n";
//                        }
//                    }
//                    compute << "}\n";
//                    numGaussian += numAtoms*usesVariable(expression[j], "gaussian");
//                    numUniform += numAtoms*usesVariable(expression[j], "uniform");
//                }
//                map<string, string> replacements;
//                replacements["COMPUTE_STEP"] = compute.str();
//                stringstream args;
//                for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++) {
//                    const CudaNonbondedUtilities::ParameterInfo& buffer = perDofValues->getBuffers()[i];
//                    string valueName = "perDofValues"+cu.intToString(i+1);
//                    args << ", __global " << buffer.getType() << "* restrict " << valueName;
//                }
//                replacements["PARAMETER_ARGUMENTS"] = args.str();
//                if (loadPosAsDelta[step])
//                    defines["LOAD_POS_AS_DELTA"] = "1";
//                else if (defines.find("LOAD_POS_AS_DELTA") != defines.end())
//                    defines.erase("LOAD_POS_AS_DELTA");
//                cu::Program program = cu.createProgram(cu.replaceStrings(CudaKernelSources::customIntegratorPerDof, replacements), defines);
//                cu::Kernel kernel = cu::Kernel(program, "computePerDof");
//                kernels[step].push_back(kernel);
//                requiredGaussian[step] = numGaussian;
//                requiredUniform[step] = numUniform;
//                int index = 0;
//                kernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, integration.getPosDelta().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, cu.getVelm().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, cu.getForce().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, integration.getStepSize().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, globalValues->getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, contextParameterValues->getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, sumBuffer->getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, integration.getRandom().getDeviceBuffer());
//                index++;
//                kernel.setArg<cu::Buffer>(index++, uniformRandoms->getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, energy->getDeviceBuffer());
//                for (int i = 0; i < (int) perDofValues->getBuffers().size(); i++)
//                    kernel.setArg<cu::Memory>(index++, perDofValues->getBuffers()[i].getMemory());
//                if (stepType[step] == CustomIntegrator::ComputeSum) {
//                    // Create a second kernel for this step that sums the values.
//
//                    program = cu.createProgram(CudaKernelSources::customIntegrator, defines);
//                    kernel = cu::Kernel(program, "computeSum");
//                    kernels[step].push_back(kernel);
//                    index = 0;
//                    kernel.setArg<cu::Buffer>(index++, sumBuffer->getDeviceBuffer());
//                    bool found = false;
//                    for (int j = 0; j < integrator.getNumGlobalVariables() && !found; j++)
//                        if (variable[step] == integrator.getGlobalVariableName(j)) {
//                            kernel.setArg<cu::Buffer>(index++, globalValues->getDeviceBuffer());
//                            kernel.setArg<cl_uint>(index++, j);
//                            found = true;
//                        }
//                    for (int j = 0; j < (int) parameterNames.size() && !found; j++)
//                        if (variable[step] == parameterNames[j]) {
//                            kernel.setArg<cu::Buffer>(index++, contextParameterValues->getDeviceBuffer());
//                            kernel.setArg<cl_uint>(index++, j);
//                            found = true;
//                            modifiesParameters = true;
//                        }
//                    if (!found)
//                        throw OpenMMException("Unknown global variable: "+variable[step]);
//                    kernel.setArg<cl_int>(index++, 3*numAtoms);
//                }
//            }
//            else if (stepType[step] == CustomIntegrator::ComputeGlobal && !merged[step]) {
//                // Compute a global value.
//
//                stringstream compute;
//                for (int i = step; i < numSteps && (i == step || merged[i]); i++)
//                    compute << "{\n" << createGlobalComputation(variable[i], expression[i], integrator, energyName[i]) << "}\n";
//                map<string, string> replacements;
//                replacements["COMPUTE_STEP"] = compute.str();
//                cu::Program program = cu.createProgram(cu.replaceStrings(CudaKernelSources::customIntegratorGlobal, replacements), defines);
//                cu::Kernel kernel = cu::Kernel(program, "computeGlobal");
//                kernels[step].push_back(kernel);
//                int index = 0;
//                kernel.setArg<cu::Buffer>(index++, integration.getStepSize().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, globalValues->getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, contextParameterValues->getDeviceBuffer());
//                index += 2;
//                kernel.setArg<cu::Buffer>(index++, energy->getDeviceBuffer());
//            }
//            else if (stepType[step] == CustomIntegrator::ConstrainPositions) {
//                // Apply position constraints.
//
//                cu::Program program = cu.createProgram(CudaKernelSources::customIntegrator, defines);
//                cu::Kernel kernel = cu::Kernel(program, "applyPositionDeltas");
//                kernels[step].push_back(kernel);
//                int index = 0;
//                kernel.setArg<cu::Buffer>(index++, cu.getPosq().getDeviceBuffer());
//                kernel.setArg<cu::Buffer>(index++, integration.getPosDelta().getDeviceBuffer());
//            }
//        }
//        
//        // Create the kernel for summing energy.
//
//        cu::Program program = cu.createProgram(CudaKernelSources::customIntegrator, defines);
//        sumEnergyKernel = cu::Kernel(program, "computeSum");
//        int index = 0;
//        sumEnergyKernel.setArg<cu::Buffer>(index++, cu.getEnergyBuffer().getDeviceBuffer());
//        sumEnergyKernel.setArg<cu::Buffer>(index++, energy->getDeviceBuffer());
//        sumEnergyKernel.setArg<cl_int>(index++, 0);
//        sumEnergyKernel.setArg<cl_int>(index++, cu.getEnergyBuffer().getSize());
//    }
//    
//    // Make sure all values (variables, parameters, etc.) stored on the device are up to date.
//    
//    if (!deviceValuesAreCurrent) {
//        perDofValues->setParameterValues(localPerDofValues);
//        deviceValuesAreCurrent = true;
//    }
//    localValuesAreCurrent = false;
//    double stepSize = integrator.getStepSize();
//    if (stepSize != prevStepSize) {
//        integration.getStepSize()[0].y = (cl_float) stepSize;
//        integration.getStepSize().upload();
//        prevStepSize = stepSize;
//    }
//    bool paramsChanged = false;
//    for (int i = 0; i < (int) parameterNames.size(); i++) {
//        float value = (float) context.getParameter(parameterNames[i]);
//        if (value != contextParameterValues->get(i)) {
//            contextParameterValues->set(i, value);
//            paramsChanged = true;
//        }
//    }
//    if (paramsChanged)
//        contextParameterValues->upload();
//
//    // Loop over computation steps in the integrator and execute them.
//
//    for (int i = 0; i < numSteps; i++) {
//        if ((needsForces[i] || needsEnergy[i]) && (!forcesAreValid || context.getLastForceGroups() != forceGroup[i])) {
//            // Recompute forces and/or energy.  Figure out what is actually needed
//            // between now and the next time they get invalidated again.
//            
//            bool computeForce = false, computeEnergy = false;
//            for (int j = i; ; j++) {
//                if (needsForces[j])
//                    computeForce = true;
//                if (needsEnergy[j])
//                    computeEnergy = true;
//                if (invalidatesForces[j])
//                    break;
//                if (j == numSteps-1)
//                    j = -1;
//                if (j == i-1)
//                    break;
//            }
//            recordChangedParameters(context);
//            context.calcForcesAndEnergy(computeForce, computeEnergy, forceGroup[i]);
//            if (computeEnergy)
//                cu.executeKernel(sumEnergyKernel, CudaContext::ThreadBlockSize, CudaContext::ThreadBlockSize);
//            forcesAreValid = true;
//        }
//        if (stepType[i] == CustomIntegrator::ComputePerDof && !merged[i]) {
//            kernels[i][0].setArg<cl_uint>(9, integration.prepareRandomNumbers(requiredGaussian[i]));
//            if (requiredUniform[i] > 0)
//                cu.executeKernel(randomKernel, numAtoms);
//            cu.executeKernel(kernels[i][0], numAtoms);
//        }
//        else if (stepType[i] == CustomIntegrator::ComputeGlobal && !merged[i]) {
//            kernels[i][0].setArg<cl_float>(3, SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber());
//            kernels[i][0].setArg<cl_float>(4, SimTKOpenMMUtilities::getNormallyDistributedRandomNumber());
//            cu.executeKernel(kernels[i][0], 1, 1);
//        }
//        else if (stepType[i] == CustomIntegrator::ComputeSum) {
//            kernels[i][0].setArg<cl_uint>(9, integration.prepareRandomNumbers(requiredGaussian[i]));
//            if (requiredUniform[i] > 0)
//                cu.executeKernel(randomKernel, numAtoms);
//            cu.executeKernel(kernels[i][0], numAtoms);
//            cu.executeKernel(kernels[i][1], CudaContext::ThreadBlockSize, CudaContext::ThreadBlockSize);
//        }
//        else if (stepType[i] == CustomIntegrator::UpdateContextState) {
//            recordChangedParameters(context);
//            context.updateContextState();
//        }
//        else if (stepType[i] == CustomIntegrator::ConstrainPositions) {
//            cu.getIntegrationUtilities().applyConstraints(integrator.getConstraintTolerance());
//            cu.executeKernel(kernels[i][0], numAtoms);
//            cu.getIntegrationUtilities().computeVirtualSites();
//        }
//        else if (stepType[i] == CustomIntegrator::ConstrainVelocities) {
//            cu.getIntegrationUtilities().applyVelocityConstraints(integrator.getConstraintTolerance());
//        }
//        if (invalidatesForces[i])
//            forcesAreValid = false;
//    }
//    recordChangedParameters(context);
//
//    // Update the time and step count.
//
//    cu.setTime(cu.getTime()+stepSize);
//    cu.setStepCount(cu.getStepCount()+1);
//}
//
//void CudaIntegrateCustomStepKernel::recordChangedParameters(ContextImpl& context) {
//    if (!modifiesParameters)
//        return;
//    contextParameterValues->download();
//    for (int i = 0; i < (int) parameterNames.size(); i++) {
//        float value = (float) context.getParameter(parameterNames[i]);
//        if (value != contextParameterValues->get(i))
//            context.setParameter(parameterNames[i], contextParameterValues->get(i));
//    }
//}
//
//void CudaIntegrateCustomStepKernel::getGlobalVariables(ContextImpl& context, vector<double>& values) const {
//    globalValues->download();
//    values.resize(numGlobalVariables);
//    for (int i = 0; i < numGlobalVariables; i++)
//        values[i] = globalValues->get(i);
//}
//
//void CudaIntegrateCustomStepKernel::setGlobalVariables(ContextImpl& context, const vector<double>& values) {
//    for (int i = 0; i < numGlobalVariables; i++)
//        globalValues->set(i, (float) values[i]);
//    globalValues->upload();
//}
//
//void CudaIntegrateCustomStepKernel::getPerDofVariable(ContextImpl& context, int variable, vector<Vec3>& values) const {
//    if (!localValuesAreCurrent) {
//        perDofValues->getParameterValues(localPerDofValues);
//        localValuesAreCurrent = true;
//    }
//    values.resize(perDofValues->getNumObjects()/3);
//    CudaArray<cl_int>& order = cu.getAtomIndex();
//    for (int i = 0; i < (int) values.size(); i++)
//        for (int j = 0; j < 3; j++)
//            values[order[i]][j] = localPerDofValues[3*i+j][variable];
//}
//
//void CudaIntegrateCustomStepKernel::setPerDofVariable(ContextImpl& context, int variable, const vector<Vec3>& values) {
//    if (!localValuesAreCurrent) {
//        perDofValues->getParameterValues(localPerDofValues);
//        localValuesAreCurrent = true;
//    }
//    CudaArray<cl_int>& order = cu.getAtomIndex();
//    for (int i = 0; i < (int) values.size(); i++)
//        for (int j = 0; j < 3; j++)
//            localPerDofValues[3*i+j][variable] = (float) values[order[i]][j];
//    deviceValuesAreCurrent = false;
//}
//
//CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
//    if (atomGroups != NULL)
//        delete atomGroups;
//}
//
//void CudaApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
//    randomSeed = thermostat.getRandomNumberSeed();
//    map<string, string> defines;
//    defines["NUM_ATOMS"] = cu.intToString(cu.getNumAtoms());
//    cu::Program program = cu.createProgram(CudaKernelSources::andersenThermostat, defines);
//    kernel = cu::Kernel(program, "applyAndersenThermostat");
//    cu.getIntegrationUtilities().initRandomNumberGenerator(randomSeed);
//
//    // Create the arrays with the group definitions.
//
//    vector<vector<int> > groups = AndersenThermostatImpl::calcParticleGroups(system);
//    atomGroups = new CudaArray<int>(cu, cu.getNumAtoms(), "atomGroups");
//    vector<int> atoms(atomGroups->getSize());
//    for (int i = 0; i < (int) groups.size(); i++) {
//        for (int j = 0; j < (int) groups[i].size(); j++)
//            atoms[groups[i][j]] = i;
//    }
//    atomGroups->upload(atoms);
//}
//
//void CudaApplyAndersenThermostatKernel::execute(ContextImpl& context) {
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//        kernel.setArg<cu::Buffer>(2, cu.getVelm().getDeviceBuffer());
//        kernel.setArg<cu::Buffer>(3, cu.getIntegrationUtilities().getStepSize().getDeviceBuffer());
//        kernel.setArg<cu::Buffer>(4, cu.getIntegrationUtilities().getRandom().getDeviceBuffer());
//        kernel.setArg<cu::Buffer>(6, atomGroups->getDeviceBuffer());
//    }
//    kernel.setArg<cl_float>(0, (cl_float) context.getParameter(AndersenThermostat::CollisionFrequency()));
//    kernel.setArg<cl_float>(1, (cl_float) (BOLTZ*context.getParameter(AndersenThermostat::Temperature())));
//    kernel.setArg<cl_uint>(5, cu.getIntegrationUtilities().prepareRandomNumbers(cu.getPaddedNumAtoms()));
//    cu.executeKernel(kernel, cu.getNumAtoms());
//}
//
//CudaApplyMonteCarloBarostatKernel::~CudaApplyMonteCarloBarostatKernel() {
//    if (savedPositions != NULL)
//        delete savedPositions;
//    if (moleculeAtoms != NULL)
//        delete moleculeAtoms;
//    if (moleculeStartIndex != NULL)
//        delete moleculeStartIndex;
//}
//
//void CudaApplyMonteCarloBarostatKernel::initialize(const System& system, const MonteCarloBarostat& thermostat) {
//    savedPositions = new CudaArray<mm_float4>(cu, cu.getPaddedNumAtoms(), "savedPositions");
//    cu::Program program = cu.createProgram(CudaKernelSources::monteCarloBarostat);
//    kernel = cu::Kernel(program, "scalePositions");
//}
//
//void CudaApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scale) {
//    if (!hasInitializedKernels) {
//        hasInitializedKernels = true;
//
//        // Create the arrays with the molecule definitions.
//
//        vector<vector<int> > molecules = context.getMolecules();
//        numMolecules = molecules.size();
//        moleculeAtoms = new CudaArray<int>(cu, cu.getNumAtoms(), "moleculeAtoms");
//        moleculeStartIndex = new CudaArray<int>(cu, numMolecules+1, "moleculeStartIndex");
//        vector<int> atoms(moleculeAtoms->getSize());
//        vector<int> startIndex(moleculeStartIndex->getSize());
//        int index = 0;
//        for (int i = 0; i < numMolecules; i++) {
//            startIndex[i] = index;
//            for (int j = 0; j < (int) molecules[i].size(); j++)
//                atoms[index++] = molecules[i][j];
//        }
//        startIndex[numMolecules] = index;
//        moleculeAtoms->upload(atoms);
//        moleculeStartIndex->upload(startIndex);
//
//        // Initialize the kernel arguments.
//        
//        kernel.setArg<cl_int>(1, numMolecules);
//        kernel.setArg<cu::Buffer>(4, cu.getPosq().getDeviceBuffer());
//        kernel.setArg<cu::Buffer>(5, moleculeAtoms->getDeviceBuffer());
//        kernel.setArg<cu::Buffer>(6, moleculeStartIndex->getDeviceBuffer());
//    }
//    cu.getQueue().enqueueCopyBuffer(cu.getPosq().getDeviceBuffer(), savedPositions->getDeviceBuffer(), 0, 0, cu.getPosq().getSize()*sizeof(mm_float4));
//    kernel.setArg<cl_float>(0, (cl_float) scale);
//    kernel.setArg<mm_float4>(2, cu.getPeriodicBoxSize());
//    kernel.setArg<mm_float4>(3, cu.getInvPeriodicBoxSize());
//    cu.executeKernel(kernel, cu.getNumAtoms());
//    for (int i = 0; i < (int) cu.getPosCellOffsets().size(); i++)
//        cu.getPosCellOffsets()[i] = mm_int4(0, 0, 0, 0);
//}
//
//void CudaApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
//    cu.getQueue().enqueueCopyBuffer(savedPositions->getDeviceBuffer(), cu.getPosq().getDeviceBuffer(), 0, 0, cu.getPosq().getSize()*sizeof(mm_float4));
//}

void CudaCalcKineticEnergyKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

double CudaCalcKineticEnergyKernel::execute(ContextImpl& context) {
    // We don't currently have a GPU kernel to do this, so we retrieve the velocities and calculate the energy
    // on the CPU.

    const vector<int>& order = cu.getAtomIndex();
    double energy = 0.0;
    if (cu.getUseDoublePrecision()) {
        double4* velm = (double4*) cu.getPinnedBuffer();
        cu.getVelm().download(velm);
        for (size_t i = 0; i < masses.size(); ++i) {
            double4 v = velm[i];
            energy += masses[order[i]]*(v.x*v.x+v.y*v.y+v.z*v.z);
        }
    }
    else {
        float4* velm = (float4*) cu.getPinnedBuffer();
        cu.getVelm().download(velm);
        for (size_t i = 0; i < masses.size(); ++i) {
            float4 v = velm[i];
            energy += masses[order[i]]*(v.x*v.x+v.y*v.y+v.z*v.z);
        }
    }
    return 0.5*energy;
}

//CudaRemoveCMMotionKernel::~CudaRemoveCMMotionKernel() {
//    if (cmMomentum != NULL)
//        delete cmMomentum;
//}
//
//void CudaRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
//    frequency = force.getFrequency();
//    int numAtoms = cu.getNumAtoms();
//    cmMomentum = new CudaArray<mm_float4>(cu, (numAtoms+CudaContext::ThreadBlockSize-1)/CudaContext::ThreadBlockSize, "cmMomentum");
//    double totalMass = 0.0;
//    for (int i = 0; i < numAtoms; i++)
//        totalMass += system.getParticleMass(i);
//    map<string, string> defines;
//    defines["INVERSE_TOTAL_MASS"] = cu.doubleToString(1.0/totalMass);
//    cu::Program program = cu.createProgram(CudaKernelSources::removeCM, defines);
//    kernel1 = cu::Kernel(program, "calcCenterOfMassMomentum");
//    kernel1.setArg<cl_int>(0, numAtoms);
//    kernel1.setArg<cu::Buffer>(1, cu.getVelm().getDeviceBuffer());
//    kernel1.setArg<cu::Buffer>(2, cmMomentum->getDeviceBuffer());
//    kernel1.setArg(3, CudaContext::ThreadBlockSize*sizeof(mm_float4), NULL);
//    kernel2 = cu::Kernel(program, "removeCenterOfMassMomentum");
//    kernel2.setArg<cl_int>(0, numAtoms);
//    kernel2.setArg<cu::Buffer>(1, cu.getVelm().getDeviceBuffer());
//    kernel2.setArg<cu::Buffer>(2, cmMomentum->getDeviceBuffer());
//    kernel2.setArg(3, CudaContext::ThreadBlockSize*sizeof(mm_float4), NULL);
//}
//
//void CudaRemoveCMMotionKernel::execute(ContextImpl& context) {
//    cu.executeKernel(kernel1, cu.getNumAtoms());
//    cu.executeKernel(kernel2, cu.getNumAtoms());
//}
