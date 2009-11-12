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
#include "openmm/LangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "kernels/gputypes.h"
#include "kernels/cudaKernels.h"
#include <cmath>

extern "C" int gpuSetConstants( gpuContext gpu );

using namespace OpenMM;
using namespace std;

void CudaCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void CudaCalcForcesAndEnergyKernel::beginForceComputation(ContextImpl& context) {
    _gpuContext* gpu = data.gpu;
    if (data.nonbondedMethod != NO_CUTOFF && data.computeForceCount%100 == 0)
        gpuReorderAtoms(gpu);
    data.computeForceCount++;
    kClearForces(gpu);
}

void CudaCalcForcesAndEnergyKernel::finishForceComputation(ContextImpl& context) {
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
    kReduceForces(gpu);
}

void CudaCalcForcesAndEnergyKernel::beginEnergyComputation(ContextImpl& context) {
    _gpuContext* gpu = data.gpu;
    if (data.nonbondedMethod != NO_CUTOFF && data.stepCount%100 == 0)
        gpuReorderAtoms(gpu);
    data.stepCount++;
    kClearEnergy(gpu);
}

double CudaCalcForcesAndEnergyKernel::finishEnergyComputation(ContextImpl& context) {
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
    return kReduceEnergy(gpu)+data.ewaldSelfEnergy;
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
        pos.x = p[0];
        pos.y = p[1];
        pos.z = p[2];
    }
    gpu->psPosq4->Upload();
    for (int i = 0; i < gpu->posCellOffsets.size(); i++)
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
        vel.x = v[0];
        vel.y = v[1];
        vel.z = v[2];
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
}

void CudaCalcHarmonicBondForceKernel::executeForces(ContextImpl& context) {
}

double CudaCalcHarmonicBondForceKernel::executeEnergy(ContextImpl& context) {
    return 0.0;
}

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
        SetCustomBondGlobalParams(&globalParamValues[0]);
}

void CudaCalcCustomBondForceKernel::executeForces(ContextImpl& context) {
    updateGlobalParams(context);
    kCalculateCustomBondForces(data.gpu);
}

double CudaCalcCustomBondForceKernel::executeEnergy(ContextImpl& context) {
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
        SetCustomBondGlobalParams(&globalParamValues[0]);
}

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
}

void CudaCalcHarmonicAngleForceKernel::executeForces(ContextImpl& context) {
}

double CudaCalcHarmonicAngleForceKernel::executeEnergy(ContextImpl& context) {
    return 0.0;
}

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
}

void CudaCalcPeriodicTorsionForceKernel::executeForces(ContextImpl& context) {
}

double CudaCalcPeriodicTorsionForceKernel::executeEnergy(ContextImpl& context) {
    return 0.0;
}

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
}

void CudaCalcRBTorsionForceKernel::executeForces(ContextImpl& context) {
}

double CudaCalcRBTorsionForceKernel::executeEnergy(ContextImpl& context) {
    return 0.0;
}

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
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        gpuSetPeriodicBoxSize(gpu, (float)boxVectors[0][0], (float)boxVectors[1][1], (float)boxVectors[2][2]);
        CudaNonbondedMethod method = NO_CUTOFF;
        if (force.getNonbondedMethod() != NonbondedForce::NoCutoff) {
            gpuSetNonbondedCutoff(gpu, (float)force.getCutoffDistance(), force.getReactionFieldDielectric());
            method = CUTOFF;
        }
        if (force.getNonbondedMethod() == NonbondedForce::CutoffPeriodic) {
            method = PERIODIC;
        }
        if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
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
        gpuSetCoulombParameters(gpu, 138.935485f, particle, c6, c12, q, symbol, exclusionList, method);

        // Compute the Ewald self energy.

        data.ewaldSelfEnergy = 0.0;
        if (force.getNonbondedMethod() == NonbondedForce::Ewald || force.getNonbondedMethod() == NonbondedForce::PME) {
            double selfEnergyScale = gpu->sim.epsfac*gpu->sim.alphaEwald/std::sqrt(PI);
                for (int i = 0; i < numParticles; i++)
                    data.ewaldSelfEnergy -= selfEnergyScale*q[i]*q[i];
        }
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
        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, particle1, particle2, c6, c12, q1, q2);
    }
}

void CudaCalcNonbondedForceKernel::executeForces(ContextImpl& context) {
}

double CudaCalcNonbondedForceKernel::executeEnergy(ContextImpl& context) {
    return 0.0;
}

CudaCalcCustomNonbondedForceKernel::~CudaCalcCustomNonbondedForceKernel() {
}

void CudaCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {
    data.hasCustomNonbonded = true;
    numParticles = force.getNumParticles();
    _gpuContext* gpu = data.gpu;

    // Identify which exceptions are actual interactions.

    vector<pair<int, int> > exclusions;
    vector<int> exceptions;
    {
        vector<double> parameters;
        for (int i = 0; i < force.getNumExceptions(); i++) {
            int particle1, particle2;
            force.getExceptionParameters(i, particle1, particle2, parameters);
            exclusions.push_back(pair<int, int>(particle1, particle2));
            if (parameters.size() > 0)
                exceptions.push_back(i);
        }
    }

    // Initialize nonbonded interactions.

    vector<int> particle(numParticles);
    vector<vector<double> > parameters(numParticles);
    vector<vector<int> > exclusionList(numParticles);
    for (int i = 0; i < numParticles; i++) {
        force.getParticleParameters(i, parameters[i]);
        particle[i] = i;
        exclusionList[i].push_back(i);
    }
    for (int i = 0; i < (int)exclusions.size(); i++) {
        exclusionList[exclusions[i].first].push_back(exclusions[i].second);
        exclusionList[exclusions[i].second].push_back(exclusions[i].first);
    }
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    gpuSetPeriodicBoxSize(gpu, (float)boxVectors[0][0], (float)boxVectors[1][1], (float)boxVectors[2][2]);
    CudaNonbondedMethod method = NO_CUTOFF;
    if (force.getNonbondedMethod() != CustomNonbondedForce::NoCutoff)
        method = CUTOFF;
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic) {
        method = PERIODIC;
    }
    data.customNonbondedMethod = method;

    // Initialize exceptions.

    int numExceptions = exceptions.size();
    vector<int> exceptionParticle1(numExceptions);
    vector<int> exceptionParticle2(numExceptions);
    vector<vector<double> > exceptionParams(numExceptions);
    for (int i = 0; i < numExceptions; i++)
        force.getExceptionParameters(exceptions[i], exceptionParticle1[i], exceptionParticle2[i], exceptionParams[i]);

    // Record the tabulated functions.

    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        bool interpolating;
        force.getFunctionParameters(i, name, values, min, max, interpolating);
        gpuSetTabulatedFunction(gpu, i, name, values, min, max, interpolating);
    }

    // Record information for the expressions.

    vector<string> paramNames;
    vector<string> combiningRules;
    for (int i = 0; i < force.getNumParameters(); i++) {
        paramNames.push_back(force.getParameterName(i));
        combiningRules.push_back(force.getParameterCombiningRule(i));
    }
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    gpuSetCustomNonbondedParameters(gpu, parameters, exclusionList, exceptionParticle1, exceptionParticle2, exceptionParams, method,
            (float)force.getCutoffDistance(), force.getEnergyFunction(), combiningRules, paramNames, globalParamNames);
    if (globalParamValues.size() > 0)
        SetCustomNonbondedGlobalParams(&globalParamValues[0]);
}

void CudaCalcCustomNonbondedForceKernel::executeForces(ContextImpl& context) {
    updateGlobalParams(context);
}

double CudaCalcCustomNonbondedForceKernel::executeEnergy(ContextImpl& context) {
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
        SetCustomNonbondedGlobalParams(&globalParamValues[0]);
}

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
}

void CudaCalcGBSAOBCForceKernel::executeForces(ContextImpl& context) {
}

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
        double charge, particleRadius, gamma, bornRadiusScaleFactor;
        force.getParticleParameters(i, charge, particleRadius, gamma );
        particle[i]                  = i;
        radius[i]                    = (float) particleRadius;
        gammas[i]                    = (float) gamma;
        scaledRadii[i]               = (float) inputScaledRadii[i];
    }
    gpuSetGBVIParameters(gpu, (float) force.getSoluteDielectric(), (float) force.getSolventDielectric(), particle,
                         radius, gammas, scaledRadii );
}

void CudaCalcGBVIForceKernel::executeForces(ContextImpl& context) {
}

double CudaCalcGBVIForceKernel::executeEnergy(ContextImpl& context) {
    return 0.0;
}

static void initializeIntegration(const System& system, CudaPlatform::PlatformData& data, const Integrator& integrator) {

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
        gpuSetCoulombParameters(gpu, 138.935485f, vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<char>(), vector<vector<int> >(), NO_CUTOFF);
        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, vector<int>(), vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<float>());
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
    kClearBornForces(gpu);
    kClearForces(gpu);
    cudaThreadSynchronize();
}

double CudaCalcGBSAOBCForceKernel::executeEnergy(ContextImpl& context) {
	return 0.0;
}

CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
}

void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    initializeIntegration(system, data, integrator);
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
    kApplyFirstShake(gpu);
    kApplyFirstSettle(gpu);
    kApplyFirstCCMA(gpu);
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
    initializeIntegration(system, data, integrator);
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
    kApplyFirstShake(gpu);
    kApplyFirstSettle(gpu);
    kApplyFirstCCMA(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kLangevinUpdatePart2(gpu);
    kApplySecondShake(gpu);
    kApplySecondSettle(gpu);
    kApplySecondCCMA(gpu);
    data.time += stepSize;
    data.stepCount++;
}

CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
}

void CudaIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    initializeIntegration(system, data, integrator);
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
    kApplyFirstShake(gpu);
    kApplyFirstSettle(gpu);
    kApplyFirstCCMA(gpu);
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
    initializeIntegration(system, data, integrator);
    prevErrorTol = -1.0;
}

void CudaIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
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
    kApplyFirstShake(gpu);
    kApplyFirstSettle(gpu);
    kApplyFirstCCMA(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kVerletUpdatePart2(gpu);
    gpu->psStepSize->Download();
    data.time += (*gpu->psStepSize)[0].y;
    if ((*gpu->psStepSize)[0].y == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
}

CudaIntegrateVariableLangevinStepKernel::~CudaIntegrateVariableLangevinStepKernel() {
}

void CudaIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    initializeIntegration(system, data, integrator);
    _gpuContext* gpu = data.gpu;
    gpu->seed = (unsigned long) integrator.getRandomNumberSeed();
    gpuInitializeRandoms(gpu);
    prevTemp = -1.0;
    prevFriction = -1.0;
    prevErrorTol = -1.0;
}

void CudaIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
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
    kApplyFirstShake(gpu);
    kApplyFirstSettle(gpu);
    kApplyFirstCCMA(gpu);
    if (data.removeCM)
        if (data.stepCount%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    kLangevinUpdatePart2(gpu);
    kApplySecondShake(gpu);
    kApplySecondSettle(gpu);
    kApplySecondCCMA(gpu);
    gpu->psStepSize->Download();
    data.time += (*gpu->psStepSize)[0].y;
    if ((*gpu->psStepSize)[0].y == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
}

CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
}

void CudaApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    _gpuContext* gpu = data.gpu;
    gpu->seed = (unsigned long) thermostat.getRandomNumberSeed();
    gpuInitializeRandoms(gpu);
    prevTemp = -1.0;
    prevFrequency = -1.0;
    prevStepSize = -1.0;
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
    kCalculateAndersenThermostat(gpu);
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
