/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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
#include "openmm/OpenMMException.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "kernels/gputypes.h"
#include "kernels/cudaKernels.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include <cmath>

extern "C" int OPENMMCUDA_EXPORT gpuSetConstants( gpuContext gpu );

using namespace OpenMM;
using namespace std;

void CudaCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void CudaCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    _gpuContext* gpu = data.gpu;
    if (data.nonbondedMethod != NO_CUTOFF && data.computeForceCount%100 == 0)
        gpuReorderAtoms(gpu);
    if ((data.hasNonbonded && data.nonbondedMethod != NO_CUTOFF && data.nonbondedMethod != CUTOFF) ||
        (data.hasCustomNonbonded && data.customNonbondedMethod != NO_CUTOFF && data.customNonbondedMethod != CUTOFF)) {
        double minAllowedSize = 1.999999*gpu->sim.nonbondedCutoff;
        if (gpu->sim.periodicBoxSizeX < minAllowedSize || gpu->sim.periodicBoxSizeY < minAllowedSize || gpu->sim.periodicBoxSizeZ < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
    }
    data.computeForceCount++;
    if (gpu->bIncludeGBSA || gpu->bIncludeGBVI)
        kClearBornSumAndForces(gpu);
    else if (includeForces)
        kClearForces(gpu);
    if (includeEnergy)
        kClearEnergy(gpu);
}

double CudaCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    _gpuContext* gpu = data.gpu;
    if (gpu->bIncludeGBSA || gpu->bIncludeGBVI) {
        gpu->bRecalculateBornRadii = true;
        kCalculateCDLJObcGbsaForces1(gpu);
        kReduceObcGbsaBornForces(gpu);
        if (gpu->bIncludeGBSA ) {
           kCalculateObcGbsaForces2(gpu);
        } else {
           kCalculateGBVIForces2(gpu);
        }
    }
    else if (data.hasNonbonded)
        kCalculateCDLJForces(gpu);
    if (data.hasCustomNonbonded)
        kCalculateCustomNonbondedForces(gpu, data.hasNonbonded);
    kCalculateLocalForces(gpu);
    if (includeForces)
        kReduceForces(gpu);
    double energy = 0.0;
    if (includeEnergy) {
        energy = kReduceEnergy(gpu)+data.ewaldSelfEnergy;
        if (data.dispersionCoefficient != 0.0)
            energy += data.dispersionCoefficient/(gpu->sim.periodicBoxSizeX*gpu->sim.periodicBoxSizeY*gpu->sim.periodicBoxSizeZ);
    }
    return energy;
}

void CudaUpdateStateDataKernel::initialize(const System& system) {
}

double CudaUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return data.time;
}

void CudaUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    data.time = time;
}

void CudaUpdateStateDataKernel::getPositions(ContextImpl& context, std::vector<Vec3>& positions) {
    _gpuContext* gpu = data.gpu;
    gpu->psPosq4->Download();
    int* order = gpu->psAtomIndex->_pSysData;
    int numParticles = context.getSystem().getNumParticles();
    positions.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        float4 pos = (*gpu->psPosq4)[i];
        int3 offset = gpu->posCellOffsets[i];
        positions[order[i]] = Vec3(pos.x-offset.x*gpu->sim.periodicBoxSizeX, pos.y-offset.y*gpu->sim.periodicBoxSizeY, pos.z-offset.z*gpu->sim.periodicBoxSizeZ);
    }
}

void CudaUpdateStateDataKernel::setPositions(ContextImpl& context, const std::vector<Vec3>& positions) {
    _gpuContext* gpu = data.gpu;
    int* order = gpu->psAtomIndex->_pSysData;
    int numParticles = context.getSystem().getNumParticles();
    for (int i = 0; i < numParticles; ++i) {
        float4& pos = (*gpu->psPosq4)[i];
        const Vec3& p = positions[order[i]];
        pos.x = (float) p[0];
        pos.y = (float) p[1];
        pos.z = (float) p[2];
    }
    gpu->psPosq4->Upload();
    for (int i = 0; i < (int) gpu->posCellOffsets.size(); i++)
        gpu->posCellOffsets[i] = make_int3(0, 0, 0);
}

void CudaUpdateStateDataKernel::getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) {
    _gpuContext* gpu = data.gpu;
    gpu->psVelm4->Download();
    int* order = gpu->psAtomIndex->_pSysData;
    int numParticles = context.getSystem().getNumParticles();
    velocities.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        float4 vel = (*gpu->psVelm4)[i];
        velocities[order[i]] = Vec3(vel.x, vel.y, vel.z);
    }
}

void CudaUpdateStateDataKernel::setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) {
    _gpuContext* gpu = data.gpu;
    int* order = gpu->psAtomIndex->_pSysData;
    int numParticles = context.getSystem().getNumParticles();
    for (int i = 0; i < numParticles; ++i) {
        float4& vel = (*gpu->psVelm4)[i];
        const Vec3& v = velocities[order[i]];
        vel.x = (float) v[0];
        vel.y = (float) v[1];
        vel.z = (float) v[2];
    }
    gpu->psVelm4->Upload();
}

void CudaUpdateStateDataKernel::getForces(ContextImpl& context, std::vector<Vec3>& forces) {
    _gpuContext* gpu = data.gpu;
    int* order = gpu->psAtomIndex->_pSysData;
    gpu->psForce4->Download();
    int numParticles = context.getSystem().getNumParticles();
    forces.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        float4 force = (*gpu->psForce4)[i];
        forces[order[i]] = Vec3(force.x, force.y, force.z);
    }
}

void CudaUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    _gpuContext* gpu = data.gpu;
    a = Vec3(gpu->sim.periodicBoxSizeX, 0, 0);
    b = Vec3(0, gpu->sim.periodicBoxSizeY, 0);
    c = Vec3(0, 0, gpu->sim.periodicBoxSizeZ);
}

void CudaUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) const {
    _gpuContext* gpu = data.gpu;
    gpuSetPeriodicBoxSize(gpu, a[0], b[1], c[2]);
    gpuSetConstants(gpu);
}

void CudaUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    throw OpenMMException("CudaPlatform does not support checkpointing");
}

void CudaUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    throw OpenMMException("CudaPlatform does not support checkpointing");
}

void CudaApplyConstraintsKernel::initialize(const System& system) {
}

void CudaApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
    kApplyConstraints(data.gpu);
}

void CudaVirtualSitesKernel::initialize(const System& system) {
}

void CudaVirtualSitesKernel::computePositions(ContextImpl& context) {
}

class CudaCalcHarmonicBondForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const HarmonicBondForce& force) : force(force) {
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

CudaCalcHarmonicBondForceKernel::~CudaCalcHarmonicBondForceKernel() {
}

void CudaCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    data.hasBonds = true;
    numBonds = force.getNumBonds();
    vector<int> particle1(numBonds);
    vector<int> particle2(numBonds);
    vector<float> length(numBonds);
    vector<float> k(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double lengthValue, kValue;
        force.getBondParameters(i, particle1[i], particle2[i], lengthValue, kValue);
        length[i] = (float) lengthValue;
        k[i] = (float) kValue;
    }
    gpuSetBondParameters(data.gpu, particle1, particle2, length, k);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcCustomBondForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomBondForce& force) : force(force) {
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

CudaCalcCustomBondForceKernel::~CudaCalcCustomBondForceKernel() {
}

void CudaCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    numBonds = force.getNumBonds();
    vector<int> particle1(numBonds);
    vector<int> particle2(numBonds);
    vector<vector<double> > params(numBonds);
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, particle1[i], particle2[i], params[i]);
    vector<string> paramNames;
    for (int i = 0; i < force.getNumPerBondParameters(); i++)
        paramNames.push_back(force.getPerBondParameterName(i));
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    gpuSetCustomBondParameters(data.gpu, particle1, particle2, params, force.getEnergyFunction(), paramNames, globalParamNames);
    if (globalParamValues.size() > 0)
        SetCustomBondGlobalParams(globalParamValues);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    updateGlobalParams(context);
    kCalculateCustomBondForces(data.gpu);
    return 0.0;
}

void CudaCalcCustomBondForceKernel::updateGlobalParams(ContextImpl& context) {
    bool changed = false;
    for (int i = 0; i < (int) globalParamNames.size(); i++) {
        float value = (float) context.getParameter(globalParamNames[i]);
        if (value != globalParamValues[i])
            changed = true;
        globalParamValues[i] = value;
    }
    if (changed)
        SetCustomBondGlobalParams(globalParamValues);
}

void CudaCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcHarmonicAngleForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const HarmonicAngleForce& force) : force(force) {
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

CudaCalcHarmonicAngleForceKernel::~CudaCalcHarmonicAngleForceKernel() {
}

void CudaCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    data.hasAngles = true;
    numAngles = force.getNumAngles();
    const float RadiansToDegrees = (float) (180.0/3.14159265);
    vector<int> particle1(numAngles);
    vector<int> particle2(numAngles);
    vector<int> particle3(numAngles);
    vector<float> angle(numAngles);
    vector<float> k(numAngles);
    for (int i = 0; i < numAngles; i++) {
        double angleValue, kValue;
        force.getAngleParameters(i, particle1[i], particle2[i], particle3[i], angleValue, kValue);
        angle[i] = (float) (angleValue*RadiansToDegrees);
        k[i] = (float) kValue;
    }
    gpuSetBondAngleParameters(data.gpu, particle1, particle2, particle3, angle, k);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcCustomAngleForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomAngleForce& force) : force(force) {
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

CudaCalcCustomAngleForceKernel::~CudaCalcCustomAngleForceKernel() {
}

void CudaCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    numAngles = force.getNumAngles();
    vector<int> particle1(numAngles);
    vector<int> particle2(numAngles);
    vector<int> particle3(numAngles);
    vector<vector<double> > params(numAngles);
    for (int i = 0; i < numAngles; i++)
        force.getAngleParameters(i, particle1[i], particle2[i], particle3[i], params[i]);
    vector<string> paramNames;
    for (int i = 0; i < force.getNumPerAngleParameters(); i++)
        paramNames.push_back(force.getPerAngleParameterName(i));
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    gpuSetCustomAngleParameters(data.gpu, particle1, particle2, particle3, params, force.getEnergyFunction(), paramNames, globalParamNames);
    if (globalParamValues.size() > 0)
        SetCustomAngleGlobalParams(globalParamValues);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    updateGlobalParams(context);
    kCalculateCustomAngleForces(data.gpu);
    return 0.0;
}

void CudaCalcCustomAngleForceKernel::updateGlobalParams(ContextImpl& context) {
    bool changed = false;
    for (int i = 0; i < (int) globalParamNames.size(); i++) {
        float value = (float) context.getParameter(globalParamNames[i]);
        if (value != globalParamValues[i])
            changed = true;
        globalParamValues[i] = value;
    }
    if (changed)
        SetCustomAngleGlobalParams(globalParamValues);
}

void CudaCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcPeriodicTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const PeriodicTorsionForce& force) : force(force) {
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
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, periodicity2, phase2, k2);
        return (periodicity1 == periodicity2 && phase1 == phase2 && k1 == k2);
    }
private:
    const PeriodicTorsionForce& force;
};

CudaCalcPeriodicTorsionForceKernel::~CudaCalcPeriodicTorsionForceKernel() {
}

void CudaCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    data.hasPeriodicTorsions = true;
    numTorsions = force.getNumTorsions();
    const float RadiansToDegrees = (float)(180.0/3.14159265);
    vector<int> particle1(numTorsions);
    vector<int> particle2(numTorsions);
    vector<int> particle3(numTorsions);
    vector<int> particle4(numTorsions);
    vector<float> k(numTorsions);
    vector<float> phase(numTorsions);
    vector<int> periodicity(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double kValue, phaseValue;
        force.getTorsionParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], periodicity[i], phaseValue, kValue);
        k[i] = (float) kValue;
        phase[i] = (float) (phaseValue*RadiansToDegrees);
    }
    gpuSetDihedralParameters(data.gpu, particle1, particle2, particle3, particle4, k, phase, periodicity);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcRBTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const RBTorsionForce& force) : force(force) {
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
        force.getTorsionParameters(group2, particle1, particle2, particle3, particle4, c0b, c1b, c2b, c3b, c4b, c5b);
        return (c0a == c0b && c1a == c1b && c2a == c2b && c3a == c3b && c4a == c4b && c5a == c5b);
    }
private:
    const RBTorsionForce& force;
};

CudaCalcRBTorsionForceKernel::~CudaCalcRBTorsionForceKernel() {
}

void CudaCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    data.hasRB = true;
    numTorsions = force.getNumTorsions();
    vector<int> particle1(numTorsions);
    vector<int> particle2(numTorsions);
    vector<int> particle3(numTorsions);
    vector<int> particle4(numTorsions);
    vector<float> c0(numTorsions);
    vector<float> c1(numTorsions);
    vector<float> c2(numTorsions);
    vector<float> c3(numTorsions);
    vector<float> c4(numTorsions);
    vector<float> c5(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        double c[6];
        force.getTorsionParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], c[0], c[1], c[2], c[3], c[4], c[5]);
        c0[i] = (float) c[0];
        c1[i] = (float) c[1];
        c2[i] = (float) c[2];
        c3[i] = (float) c[3];
        c4[i] = (float) c[4];
        c5[i] = (float) c[5];
    }
    gpuSetRbDihedralParameters(data.gpu, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void CudaCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcCMAPTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CMAPTorsionForce& force) : force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumTorsions();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
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
    if (torsionIndices != NULL)
        delete torsionIndices;
}

void CudaCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    if (numTorsions == 0)
        return;
    int numMaps = force.getNumMaps();
    vector<float4> coeffVec;
    vector<int2> mapPositionsVec(numMaps);
    vector<double> energy;
    vector<vector<double> > c;
    int currentPosition = 0;
    mapPositions = new CUDAStream<int2>(numMaps, 1, "cmapTorsionMapPositions");
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        (*mapPositions)[i] = make_int2(currentPosition, size);
        currentPosition += 4*size*size;
        for (int j = 0; j < size*size; j++) {
            coeffVec.push_back(make_float4(c[j][0], c[j][1], c[j][2], c[j][3]));
            coeffVec.push_back(make_float4(c[j][4], c[j][5], c[j][6], c[j][7]));
            coeffVec.push_back(make_float4(c[j][8], c[j][9], c[j][10], c[j][11]));
            coeffVec.push_back(make_float4(c[j][12], c[j][13], c[j][14], c[j][15]));
        }
    }
    coefficients = new CUDAStream<float4>((int) coeffVec.size(), 1, "cmapTorsionCoefficients");;
    for (int i = 0; i < (int) coeffVec.size(); i++)
        (*coefficients)[i] = coeffVec[i];
    torsionMaps = new CUDAStream<int>(numTorsions, 1, "cmapTorsionMaps");
    torsionIndices = new CUDAStream<int4>(4*numTorsions, 1, "cmapTorsionIndices");
    vector<int> forceBufferCounter(system.getNumParticles(), 0);
    for (int i = 0; i < numTorsions; i++) {
        int map, a1, a2, a3, a4, b1, b2, b3, b4;
        force.getTorsionParameters(i, map, a1, a2, a3, a4, b1, b2, b3, b4);
        (*torsionMaps)[i] = map;
        (*torsionIndices)[i*4] = make_int4(a1, a2, a3, a4);
        (*torsionIndices)[i*4+1] = make_int4(b1, b2, b3, b4);
        (*torsionIndices)[i*4+2] = make_int4(forceBufferCounter[a1]++, forceBufferCounter[a2]++, forceBufferCounter[a3]++, forceBufferCounter[a4]++);
        (*torsionIndices)[i*4+3] = make_int4(forceBufferCounter[b1]++, forceBufferCounter[b2]++, forceBufferCounter[b3]++, forceBufferCounter[b4]++);
    }
    coefficients->Upload();
    mapPositions->Upload();
    torsionMaps->Upload();
    torsionIndices->Upload();
    int maxBuffers = 1;
    for (int i = 0; i < (int) forceBufferCounter.size(); i++)
        maxBuffers = max(maxBuffers, forceBufferCounter[i]);
    if (maxBuffers > data.gpu->sim.outputBuffers)
        data.gpu->sim.outputBuffers = maxBuffers;
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if( numTorsions )
        kCalculateCMAPTorsionForces(data.gpu, *coefficients, *mapPositions, *torsionIndices, *torsionMaps);
    return 0.0;
}

class CudaCalcCustomTorsionForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomTorsionForce& force) : force(force) {
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

CudaCalcCustomTorsionForceKernel::~CudaCalcCustomTorsionForceKernel() {
}

void CudaCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    vector<int> particle1(numTorsions);
    vector<int> particle2(numTorsions);
    vector<int> particle3(numTorsions);
    vector<int> particle4(numTorsions);
    vector<vector<double> > params(numTorsions);
    for (int i = 0; i < numTorsions; i++)
        force.getTorsionParameters(i, particle1[i], particle2[i], particle3[i], particle4[i], params[i]);
    vector<string> paramNames;
    for (int i = 0; i < force.getNumPerTorsionParameters(); i++)
        paramNames.push_back(force.getPerTorsionParameterName(i));
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    gpuSetCustomTorsionParameters(data.gpu, particle1, particle2, particle3, particle4, params, force.getEnergyFunction(), paramNames, globalParamNames);
    if (globalParamValues.size() > 0)
        SetCustomTorsionGlobalParams(globalParamValues);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    updateGlobalParams(context);
    kCalculateCustomTorsionForces(data.gpu);
    return 0.0;
}

void CudaCalcCustomTorsionForceKernel::updateGlobalParams(ContextImpl& context) {
    bool changed = false;
    for (int i = 0; i < (int) globalParamNames.size(); i++) {
        float value = (float) context.getParameter(globalParamNames[i]);
        if (value != globalParamValues[i])
            changed = true;
        globalParamValues[i] = value;
    }
    if (changed)
        SetCustomTorsionGlobalParams(globalParamValues);
}

void CudaCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
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

CudaCalcNonbondedForceKernel::~CudaCalcNonbondedForceKernel() {
}

void CudaCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {
    data.hasNonbonded = true;
    numParticles = force.getNumParticles();
    _gpuContext* gpu = data.gpu;

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
    
    {
        vector<int> particle(numParticles);
        vector<float> c6(numParticles);
        vector<float> c12(numParticles);
        vector<float> q(numParticles);
        vector<char> symbol;
        vector<vector<int> > exclusionList(numParticles);
        for (int i = 0; i < numParticles; i++) {
            double charge, radius, depth;
            force.getParticleParameters(i, charge, radius, depth);
            particle[i] = i;
            q[i] = (float) charge;
            c6[i] = (float) (4*depth*pow(radius, 6.0));
            c12[i] = (float) (4*depth*pow(radius, 12.0));
            exclusionList[i].push_back(i);
        }
        for (int i = 0; i < (int)exclusions.size(); i++) {
            exclusionList[exclusions[i].first].push_back(exclusions[i].second);
            exclusionList[exclusions[i].second].push_back(exclusions[i].first);
        }
        CudaNonbondedMethod method = NO_CUTOFF;
        if (force.getNonbondedMethod() != NonbondedForce::NoCutoff) {
            gpuSetNonbondedCutoff(gpu, (float) force.getCutoffDistance(), (float) force.getReactionFieldDielectric());
            method = CUTOFF;
        }
        if (force.getNonbondedMethod() == NonbondedForce::CutoffPeriodic) {
            method = PERIODIC;
        }
        if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
            if (force.getReciprocalSpaceForceGroup() > 0)
                throw OpenMMException("CudaPlatform does not support force groups");
            if (force.getNonbondedMethod() == NonbondedForce::Ewald) {
                double alpha;
                int kmaxx, kmaxy, kmaxz;
                NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmaxx, kmaxy, kmaxz);
                gpuSetEwaldParameters(gpu, (float) alpha, kmaxx, kmaxy, kmaxz);
                method = EWALD;
            }
            else {
                double alpha;
                int gridSizeX, gridSizeY, gridSizeZ;
                NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSizeX, gridSizeY, gridSizeZ);
                gpuSetPMEParameters(gpu, (float) alpha, gridSizeX, gridSizeY, gridSizeZ);
                method = PARTICLE_MESH_EWALD;
            }
        }
        data.nonbondedMethod = method;
        gpuSetCoulombParameters(gpu, (float) ONE_4PI_EPS0, particle, c6, c12, q, symbol, exclusionList, method);

        // Compute the Ewald self energy.

        data.ewaldSelfEnergy = 0.0;
        if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
            double selfEnergyScale = gpu->sim.epsfac*gpu->sim.alphaEwald/std::sqrt(PI);
                for (int i = 0; i < numParticles; i++)
                    data.ewaldSelfEnergy -= selfEnergyScale*q[i]*q[i];
        }

        // Compute the long range dispersion correction.

        if (force.getUseDispersionCorrection())
            data.dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
        else
            data.dispersionCoefficient = 0.0;
    }

    // Initialize 1-4 nonbonded interactions.
    
    {
        int numExceptions = exceptions.size();
        vector<int> particle1(numExceptions);
        vector<int> particle2(numExceptions);
        vector<float> c6(numExceptions);
        vector<float> c12(numExceptions);
        vector<float> q1(numExceptions);
        vector<float> q2(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double charge, sig, eps;
            force.getExceptionParameters(exceptions[i], particle1[i], particle2[i], charge, sig, eps);
            c6[i] = (float) (4*eps*pow(sig, 6.0));
            c12[i] = (float) (4*eps*pow(sig, 12.0));
            q1[i] = (float) charge;
            q2[i] = 1.0f;
        }
        gpuSetLJ14Parameters(gpu, (float) ONE_4PI_EPS0, 1.0f, particle1, particle2, c6, c12, q1, q2);
    }
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    return 0.0;
}

void CudaCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcCustomNonbondedForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomNonbondedForce& force) : force(force) {
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

CudaCalcCustomNonbondedForceKernel::~CudaCalcCustomNonbondedForceKernel() {
}

void CudaCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    data.hasCustomNonbonded = true;
    numParticles = force.getNumParticles();
    _gpuContext* gpu = data.gpu;

    // Initialize nonbonded interactions.

    vector<int> particle(numParticles);
    vector<vector<double> > parameters(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters[i]);
        particle[i] = i;
        exclusionList[i].push_back(i);
    }
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusionList[particle1].push_back(particle2);
        exclusionList[particle2].push_back(particle1);
    }
    CudaNonbondedMethod method = NO_CUTOFF;
    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
        method = CUTOFF;
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic) {
        method = PERIODIC;
    }
    data.customNonbondedMethod = method;

    // Record the tabulated functions.

    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        force.getFunctionParameters(i, name, values, min, max);
        gpuSetTabulatedFunction(gpu, i, name, values, min, max);
    }

    // Record information for the expressions.

    vector<string> paramNames;
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        paramNames.push_back(force.getPerParticleParameterName(i));
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    gpuSetCustomNonbondedParameters(gpu, parameters, exclusionList, method, (float) force.getCutoffDistance(), force.getEnergyFunction(), paramNames, globalParamNames);
    if (globalParamValues.size() > 0)
        SetCustomNonbondedGlobalParams(globalParamValues);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    updateGlobalParams(context);
    return 0.0;
}

void CudaCalcCustomNonbondedForceKernel::updateGlobalParams(ContextImpl& context) {
    bool changed = false;
    for (int i = 0; i < (int) globalParamNames.size(); i++) {
        float value = (float) context.getParameter(globalParamNames[i]);
        if (value != globalParamValues[i])
            changed = true;
        globalParamValues[i] = value;
    }
    if (changed)
        SetCustomNonbondedGlobalParams(globalParamValues);
}

void CudaCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
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

CudaCalcGBSAOBCForceKernel::~CudaCalcGBSAOBCForceKernel() {
}

void CudaCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {

    int numParticles = system.getNumParticles();
    _gpuContext* gpu = data.gpu;
    vector<float> radius(numParticles);
    vector<float> scale(numParticles);
    vector<float> charge(numParticles);
    for (int i = 0; i < numParticles; i++) {
        double particleCharge, particleRadius, scalingFactor;
        force.getParticleParameters(i, particleCharge, particleRadius, scalingFactor);
        radius[i] = (float) particleRadius;
        scale[i] = (float) scalingFactor;
        charge[i] = (float) particleCharge;
    }
    gpuSetObcParameters(gpu, (float) force.getSoluteDielectric(), (float) force.getSolventDielectric(), radius, scale, charge);
    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
	return 0.0;
}

void CudaCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

class CudaCalcGBVIForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const GBVIForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, gamma1, gamma2;
        force.getParticleParameters(particle1, charge1, radius1, gamma1);
        force.getParticleParameters(particle2, charge2, radius2, gamma2);
        return (charge1 == charge2 && radius1 == radius2 && gamma1 == gamma2);
    }
private:
    const GBVIForce& force;
};

CudaCalcGBVIForceKernel::~CudaCalcGBVIForceKernel() {
}

void CudaCalcGBVIForceKernel::initialize(const System& system, const GBVIForce& force, const std::vector<double> & inputScaledRadii) {

    int numParticles = system.getNumParticles();
    _gpuContext* gpu = data.gpu;

    vector<int> particle(numParticles);
    vector<float> radius(numParticles);
    vector<float> scaledRadii(numParticles);
    vector<float> gammas(numParticles);

    for (int i = 0; i < numParticles; i++) {
        double charge, particleRadius, gamma;
        force.getParticleParameters(i, charge, particleRadius, gamma );
        particle[i]                  = i;
        radius[i]                    = (float) particleRadius;
        gammas[i]                    = (float) gamma;
        scaledRadii[i]               = (float) inputScaledRadii[i];
    }

    int gbviBornRadiusScalingMethod;
    if( force.getBornRadiusScalingMethod() ==  GBVIForce::QuinticSpline ){
        gbviBornRadiusScalingMethod = 1;
    } else {
        gbviBornRadiusScalingMethod = 2;
    }
    gpuSetGBVIParameters(gpu, (float) force.getSoluteDielectric(), (float) force.getSolventDielectric(), particle,
                         radius, gammas, scaledRadii, gbviBornRadiusScalingMethod,
                         static_cast<float>(force.getQuinticLowerLimitFactor()),
                         static_cast<float>(force.getQuinticUpperBornRadiusLimit()) );

    data.gpu->forces.push_back(new ForceInfo(force));
}

double CudaCalcGBVIForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
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
}

void CudaCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    numParticles = force.getNumParticles();
    vector<int> particle(numParticles);
    vector<vector<double> > params(numParticles);
    for (int i = 0; i < numParticles; i++)
        force.getParticleParameters(i, particle[i], params[i]);
    vector<string> paramNames;
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        paramNames.push_back(force.getPerParticleParameterName(i));
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    gpuSetCustomExternalParameters(data.gpu, particle, params, force.getEnergyFunction(), paramNames, globalParamNames);
    if (globalParamValues.size() > 0)
        SetCustomExternalGlobalParams(globalParamValues);
    data.gpu->forces.push_back(new ForceInfo(force, system.getNumParticles()));
}

double CudaCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    updateGlobalParams(context);
    kCalculateCustomExternalForces(data.gpu);
    return 0.0;
}

void CudaCalcCustomExternalForceKernel::updateGlobalParams(ContextImpl& context) {
    bool changed = false;
    for (int i = 0; i < (int) globalParamNames.size(); i++) {
        float value = (float) context.getParameter(globalParamNames[i]);
        if (value != globalParamValues[i])
            changed = true;
        globalParamValues[i] = value;
    }
    if (changed)
        SetCustomExternalGlobalParams(globalParamValues);
}

void CudaCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
    throw OpenMMException("CudaPlatform does not support copyParametersToContext");
}

void OPENMMCUDA_EXPORT OpenMM::cudaOpenMMInitializeIntegration(const System& system, CudaPlatform::PlatformData& data, const Integrator& integrator) {

    // Initialize any terms that haven't already been handled by a Force.

    _gpuContext* gpu = data.gpu;
    if (!data.hasBonds)
        gpuSetBondParameters(gpu, vector<int>(), vector<int>(), vector<float>(), vector<float>());
    if (!data.hasAngles)
        gpuSetBondAngleParameters(gpu, vector<int>(), vector<int>(), vector<int>(), vector<float>(), vector<float>());
    if (!data.hasPeriodicTorsions)
        gpuSetDihedralParameters(gpu, vector<int>(), vector<int>(), vector<int>(), vector<int>(), vector<float>(), vector<float>(), vector<int>());
    if (!data.hasRB)
        gpuSetRbDihedralParameters(gpu, vector<int>(), vector<int>(), vector<int>(), vector<int>(), vector<float>(), vector<float>(),
                vector<float>(), vector<float>(), vector<float>(), vector<float>());
    if (!data.hasNonbonded) {
        gpuSetCoulombParameters(gpu, (float) ONE_4PI_EPS0, vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<char>(), vector<vector<int> >(), NO_CUTOFF);
        gpuSetLJ14Parameters(gpu, (float) ONE_4PI_EPS0, 1.0f, vector<int>(), vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<float>());
        if (gpu->bIncludeGBSA || gpu->bIncludeGBVI)
            throw OpenMMException("CudaPlatform requires GBSAOBCForce and GBVIForce to be used with a NonbondedForce");
    }
    
    // Set masses.
    
    int numParticles = system.getNumParticles();
    vector<float> mass(numParticles);
    for (int i = 0; i < numParticles; i++)
        mass[i] = (float) system.getParticleMass(i);
    gpuSetMass(gpu, mass);
    
    // Set constraints.
    
    int numConstraints = system.getNumConstraints();
    vector<int> particle1(numConstraints);
    vector<int> particle2(numConstraints);
    vector<float> distance(numConstraints);
    vector<float> invMass1(numConstraints);
    vector<float> invMass2(numConstraints);
    for (int i = 0; i < numConstraints; i++) {
        int particle1Index, particle2Index;
        double constraintDistance;
        system.getConstraintParameters(i, particle1Index, particle2Index, constraintDistance);
        particle1[i] = particle1Index;
        particle2[i] = particle2Index;
        distance[i] = (float) constraintDistance;
        invMass1[i] = 1.0f/mass[particle1Index];
        invMass2[i] = 1.0f/mass[particle2Index];
    }
    gpuSetConstraintParameters(gpu, particle1, particle2, distance, invMass1, invMass2, (float)integrator.getConstraintTolerance());
    
    // Finish initialization.

    gpuBuildThreadBlockWorkList(gpu);
    gpuBuildExclusionList(gpu);
    gpuBuildOutputBuffers(gpu);
    gpuSetConstants(gpu);
    if (gpu->bIncludeGBSA || gpu->bIncludeGBVI)
        kClearBornSumAndForces(gpu);
    else
        kClearForces(gpu);
    cudaThreadSynchronize();
}

CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
}

void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    cudaOpenMMInitializeIntegration(system, data, integrator);
    prevStepSize = -1.0;
}

void CudaIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    _gpuContext* gpu = data.gpu;
    double stepSize = integrator.getStepSize();
    if (stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        gpuSetVerletIntegrationParameters(gpu, (float) stepSize, 0.0f);
        gpuSetConstants(gpu);
        prevStepSize = stepSize;
    }
    kVerletUpdatePart1(gpu);
    kApplyShake(gpu);
    kApplySettle(gpu);
    kApplyCCMA(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kVerletUpdatePart2(gpu);
    data.time += stepSize;
    data.stepCount++;
}

CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
}

void CudaIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    cudaOpenMMInitializeIntegration(system, data, integrator);
    _gpuContext* gpu = data.gpu;
    gpu->seed = (unsigned long) integrator.getRandomNumberSeed();
    gpuInitializeRandoms(gpu);
    prevTemp = -1.0;
    prevFriction = -1.0;
    prevStepSize = -1.0;
}

void CudaIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    _gpuContext* gpu = data.gpu;
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        gpuSetLangevinIntegrationParameters(gpu, (float) tau, (float) stepSize, (float) temperature, 0.0f);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    kLangevinUpdatePart1(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kLangevinUpdatePart2(gpu);
    kApplyShake(gpu);
    kApplySettle(gpu);
    kApplyCCMA(gpu);
    kSetVelocitiesFromPositions(gpu);
    data.time += stepSize;
    data.stepCount++;
}

CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
}

void CudaIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    cudaOpenMMInitializeIntegration(system, data, integrator);
    _gpuContext* gpu = data.gpu;
    gpu->seed = (unsigned long) integrator.getRandomNumberSeed();
    gpuInitializeRandoms(gpu);
    prevTemp = -1.0;
    prevFriction = -1.0;
    prevStepSize = -1.0;
}

void CudaIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    _gpuContext* gpu = data.gpu;
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        gpuSetBrownianIntegrationParameters(gpu, (float) tau, (float) stepSize, (float) temperature);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    kBrownianUpdatePart1(gpu);
    kApplyShake(gpu);
    kApplySettle(gpu);
    kApplyCCMA(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kBrownianUpdatePart2(gpu);
    data.time += stepSize;
    data.stepCount++;
}

CudaIntegrateVariableVerletStepKernel::~CudaIntegrateVariableVerletStepKernel() {
}

void CudaIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    cudaOpenMMInitializeIntegration(system, data, integrator);
    prevErrorTol = -1.0;
}

double CudaIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    _gpuContext* gpu = data.gpu;
    double errorTol = integrator.getErrorTolerance();
    if (errorTol != prevErrorTol) {
        // Initialize the GPU parameters.

        gpuSetVerletIntegrationParameters(gpu, 0.0f, (float) errorTol);
        gpuSetConstants(gpu);
        prevErrorTol = errorTol;
    }
    float maxStepSize = (float)(maxTime-data.time);
    kSelectVerletStepSize(gpu, maxStepSize);
    kVerletUpdatePart1(gpu);
    kApplyShake(gpu);
    kApplySettle(gpu);
    kApplyCCMA(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kVerletUpdatePart2(gpu);
    gpu->psStepSize->Download();
    data.time += (*gpu->psStepSize)[0].y;
    if ((*gpu->psStepSize)[0].y == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
    return (*gpu->psStepSize)[0].y;
}

CudaIntegrateVariableLangevinStepKernel::~CudaIntegrateVariableLangevinStepKernel() {
}

void CudaIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    cudaOpenMMInitializeIntegration(system, data, integrator);
    _gpuContext* gpu = data.gpu;
    gpu->seed = (unsigned long) integrator.getRandomNumberSeed();
    gpuInitializeRandoms(gpu);
    prevTemp = -1.0;
    prevFriction = -1.0;
    prevErrorTol = -1.0;
}

double CudaIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    _gpuContext* gpu = data.gpu;
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double errorTol = integrator.getErrorTolerance();
    if (temperature != prevTemp || friction != prevFriction || errorTol != prevErrorTol) {
        // Initialize the GPU parameters.

        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        gpuSetLangevinIntegrationParameters(gpu, (float) tau, 0.0f, (float) temperature, (float) errorTol);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevTemp = temperature;
        prevFriction = friction;
        prevErrorTol = errorTol;
    }
    float maxStepSize = (float)(maxTime-data.time);
    kSelectLangevinStepSize(gpu, maxStepSize);
    kLangevinUpdatePart1(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kLangevinUpdatePart2(gpu);
    kApplyShake(gpu);
    kApplySettle(gpu);
    kApplyCCMA(gpu);
    kSetVelocitiesFromPositions(gpu);
    gpu->psStepSize->Download();
    data.time += (*gpu->psStepSize)[0].y;
    if ((*gpu->psStepSize)[0].y == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
    return (*gpu->psStepSize)[0].y;
}

CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
    if (atomGroups != NULL)
        delete atomGroups;
}

void CudaApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    _gpuContext* gpu = data.gpu;
    gpu->seed = (unsigned long) thermostat.getRandomNumberSeed();
    gpuInitializeRandoms(gpu);
    prevTemp = -1.0;
    prevFrequency = -1.0;
    prevStepSize = -1.0;

    // Create the arrays with the group definitions.

    vector<vector<int> > groups = AndersenThermostatImpl::calcParticleGroups(system);
    atomGroups = new CUDAStream<int>(system.getNumParticles(), 1, "atomGroups");
    for (int i = 0; i < (int) groups.size(); i++) {
        for (int j = 0; j < (int) groups[i].size(); j++)
            (*atomGroups)[groups[i][j]] = i;
    }
    atomGroups->Upload();
}

void CudaApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    _gpuContext* gpu = data.gpu;
    double temperature = context.getParameter(AndersenThermostat::Temperature());
    double frequency = context.getParameter(AndersenThermostat::CollisionFrequency());
    double stepSize = context.getIntegrator().getStepSize();
    if (temperature != prevTemp || frequency != prevFrequency || stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        gpuSetAndersenThermostatParameters(gpu, (float) temperature, (float) frequency);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevTemp = temperature;
        prevFrequency = frequency;
        prevStepSize = stepSize;
    }
    kCalculateAndersenThermostat(gpu, *atomGroups);
}

CudaApplyMonteCarloBarostatKernel::~CudaApplyMonteCarloBarostatKernel() {
    if (moleculeAtoms != NULL)
        delete moleculeAtoms;
    if (moleculeStartIndex != NULL)
        delete moleculeStartIndex;
}

void CudaApplyMonteCarloBarostatKernel::initialize(const System& system, const MonteCarloBarostat& thermostat) {
}

void CudaApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scale) {
    if (!hasInitializedMolecules) {
        hasInitializedMolecules = true;

        // Create the arrays with the molecule definitions.

        vector<vector<int> > molecules = context.getMolecules();
        numMolecules = molecules.size();
        moleculeAtoms = new CUDAStream<int>(context.getSystem().getNumParticles(), 1, "moleculeAtoms");
        moleculeStartIndex = new CUDAStream<int>(numMolecules+1, 1, "moleculeStartIndex");
        int index = 0;
        for (int i = 0; i < numMolecules; i++) {
            (*moleculeStartIndex)[i] = index;
            for (int j = 0; j < (int) molecules[i].size(); j++)
                (*moleculeAtoms)[index++] = molecules[i][j];
        }
        (*moleculeStartIndex)[numMolecules] = index;
        moleculeAtoms->Upload();
        moleculeStartIndex->Upload();
    }
    _gpuContext* gpu = data.gpu;
    gpu->psPosqP4->CopyFrom(*gpu->psPosq4);
    kScaleAtomCoordinates(gpu, scale, *moleculeAtoms, *moleculeStartIndex);
    for (int i = 0; i < (int) gpu->posCellOffsets.size(); i++)
        gpu->posCellOffsets[i] = make_int3(0, 0, 0);
}

void CudaApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    _gpuContext* gpu = data.gpu;
    gpu->psPosq4->CopyFrom(*gpu->psPosqP4);
}

void CudaCalcKineticEnergyKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

double CudaCalcKineticEnergyKernel::execute(ContextImpl& context) {
    // We don't currently have a GPU kernel to do this, so we retrieve the velocities and calculate the energy
    // on the CPU.
    
    _gpuContext* gpu = data.gpu;
    gpu->psVelm4->Download();
    double energy = 0.0;
    for (int i = 0; i < (int) masses.size(); ++i) {
        float4 v = (*gpu->psVelm4)[i];
        energy += masses[i]*(v.x*v.x+v.y*v.y+v.z*v.z);
    }
    return 0.5*energy;
}

void CudaRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    data.removeCM = true;
    data.cmMotionFrequency = force.getFrequency();
}

void CudaRemoveCMMotionKernel::execute(ContextImpl& context) {
}
