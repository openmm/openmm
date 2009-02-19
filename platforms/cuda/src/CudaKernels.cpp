/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "OpenMMContext.h"

#include "CudaKernels.h"
#include "CudaStreamImpl.h"
#include "LangevinIntegrator.h"
#include "ReferencePlatform.h"
#include "internal/OpenMMContextImpl.h"
#include "kernels/gputypes.h"
#include "kernels/cudaKernels.h"
#include <cmath>

extern "C" int gpuSetConstants( gpuContext gpu );

using namespace OpenMM;
using namespace std;

static void calcForces(OpenMMContextImpl& context, CudaPlatform::PlatformData& data) {
    _gpuContext* gpu = data.gpu;
    if (data.useOBC) {
        gpu->bRecalculateBornRadii = true;
        kCalculateCDLJObcGbsaForces1(gpu);
        kReduceObcGbsaBornForces(gpu);
        kCalculateObcGbsaForces2(gpu);
    }
    else {
        kClearForces(gpu);
        kCalculateCDLJForces(gpu);
    }
    kCalculateLocalForces(gpu);
    kReduceBornSumAndForces(gpu);
}

static double calcEnergy(OpenMMContextImpl& context, System& system) {
    // We don't currently have GPU kernels to calculate energy, so instead we have the reference
    // platform do it.  This is VERY slow.
    
    LangevinIntegrator integrator(0.0, 1.0, 0.0);
    ReferencePlatform platform;
    OpenMMContext refContext(system, integrator, platform);
    const Stream& positions = context.getPositions();
    double* posData = new double[positions.getSize()*3];
    positions.saveToArray(posData);
    vector<Vec3> pos(positions.getSize());
    for (int i = 0; i < pos.size(); i++)
        pos[i] = Vec3(posData[3*i], posData[3*i+1], posData[3*i+2]);
    delete[] posData;
    refContext.setPositions(pos);
    return refContext.getState(State::Energy).getPotentialEnergy();
}

void CudaInitializeForcesKernel::initialize(const System& system) {
}

void CudaInitializeForcesKernel::execute(OpenMMContextImpl& context) {
}

CudaCalcHarmonicBondForceKernel::~CudaCalcHarmonicBondForceKernel() {
}

void CudaCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    if (data.primaryKernel == NULL)
        data.primaryKernel = this;
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

void CudaCalcHarmonicBondForceKernel::executeForces(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        calcForces(context, data);
}

double CudaCalcHarmonicBondForceKernel::executeEnergy(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        return calcEnergy(context, system);
    return 0.0;
}

CudaCalcHarmonicAngleForceKernel::~CudaCalcHarmonicAngleForceKernel() {
}

void CudaCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    if (data.primaryKernel == NULL)
        data.primaryKernel = this;
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

void CudaCalcHarmonicAngleForceKernel::executeForces(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        calcForces(context, data);
}

double CudaCalcHarmonicAngleForceKernel::executeEnergy(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        return calcEnergy(context, system);
    return 0.0;
}

CudaCalcPeriodicTorsionForceKernel::~CudaCalcPeriodicTorsionForceKernel() {
}

void CudaCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    if (data.primaryKernel == NULL)
        data.primaryKernel = this;
    data.hasPeriodicTorsions = true;
    numTorsions = force.getNumTorsions();
    const float RadiansToDegrees = 180.0/3.14159265;
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

void CudaCalcPeriodicTorsionForceKernel::executeForces(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        calcForces(context, data);
}

double CudaCalcPeriodicTorsionForceKernel::executeEnergy(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        return calcEnergy(context, system);
    return 0.0;
}

CudaCalcRBTorsionForceKernel::~CudaCalcRBTorsionForceKernel() {
}

void CudaCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    if (data.primaryKernel == NULL)
        data.primaryKernel = this;
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

void CudaCalcRBTorsionForceKernel::executeForces(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        calcForces(context, data);
}

double CudaCalcRBTorsionForceKernel::executeEnergy(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        return calcEnergy(context, system);
    return 0.0;
}

CudaCalcNonbondedForceKernel::~CudaCalcNonbondedForceKernel() {
}

void CudaCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force, const std::vector<std::set<int> >& exclusions) {
    if (data.primaryKernel == NULL)
        data.primaryKernel = this;
    data.hasNonbonded = true;
    numParticles = force.getNumParticles();
    num14 = force.getNumNonbonded14();
    _gpuContext* gpu = data.gpu;
    
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
            exclusionList[i] = vector<int>(exclusions[i].begin(), exclusions[i].end());
            exclusionList[i].push_back(i);
        }
        gpuSetCoulombParameters(gpu, 138.935485f, particle, c6, c12, q, symbol, exclusionList);
    }

    // Initialize 1-4 nonbonded interactions.
    
    {
        vector<int> particle1(num14);
        vector<int> particle2(num14);
        vector<float> c6(num14);
        vector<float> c12(num14);
        vector<float> q1(num14);
        vector<float> q2(num14);
        for (int i = 0; i < num14; i++) {
            double charge, sig, eps;
            force.getNonbonded14Parameters(i, particle1[i], particle2[i], charge, sig, eps);
            c6[i] = (float) (4*eps*pow(sig, 6.0));
            c12[i] = (float) (4*eps*pow(sig, 12.0));
            q1[i] = (float) charge;
            q2[i] = 1.0f;
        }
        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, particle1, particle2, c6, c12, q1, q2);
    }
}

void CudaCalcNonbondedForceKernel::executeForces(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        calcForces(context, data);
}

double CudaCalcNonbondedForceKernel::executeEnergy(OpenMMContextImpl& context) {
    if (data.primaryKernel == this)
        return calcEnergy(context, system);
    return 0.0;
}

CudaCalcGBSAOBCForceKernel::~CudaCalcGBSAOBCForceKernel() {
}

void CudaCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {

    int numParticles = system.getNumParticles();
    _gpuContext* gpu = data.gpu;
    vector<int> particle(numParticles);
    vector<float> radius(numParticles);
    vector<float> scale(numParticles);
    for (int i = 0; i < numParticles; i++) {
        double charge, particleRadius, scalingFactor;
        force.getParticleParameters(i, charge, particleRadius, scalingFactor);
        particle[i] = i;
        radius[i] = (float) particleRadius;
        scale[i] = (float) scalingFactor;
    }
    gpuSetObcParameters(gpu, (float) force.getSoluteDielectric(), (float) force.getSolventDielectric(), particle, radius, scale);
    data.useOBC = true;
}

void CudaCalcGBSAOBCForceKernel::executeForces(OpenMMContextImpl& context) {
}

static void initializeIntegration(const System& system, CudaPlatform::PlatformData& data, const Integrator& integrator) {
    
    // Set masses.
    
    _gpuContext* gpu = data.gpu;
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
    gpuSetShakeParameters(gpu, particle1, particle2, distance, invMass1, invMass2, integrator.getConstraintTolerance());
    
    // Initialize any terms that haven't already been handled by a Force.
    
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
        gpuSetCoulombParameters(gpu, 138.935485f, vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<char>(), vector<vector<int> >());
        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, vector<int>(), vector<int>(), vector<float>(), vector<float>(), vector<float>(), vector<float>());
    }
    
    // Finish initialization.
    
    gpuBuildThreadBlockWorkList(gpu);
    gpuBuildExclusionList(gpu);
    gpuBuildOutputBuffers(gpu);
    gpuSetConstants(gpu);
    kClearBornForces(gpu);
    kClearForces(gpu);
    cudaThreadSynchronize();
}

double CudaCalcGBSAOBCForceKernel::executeEnergy(OpenMMContextImpl& context) {
	return 0.0;
}

CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
}

void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    initializeIntegration(system, data, integrator);
    prevStepSize = -1.0;
}

void CudaIntegrateVerletStepKernel::execute(OpenMMContextImpl& context, const VerletIntegrator& integrator) {
    _gpuContext* gpu = data.gpu;
    double stepSize = integrator.getStepSize();
    if (stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        gpuSetVerletIntegrationParameters(gpu, (float) stepSize);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevStepSize = stepSize;
    }
    kVerletUpdatePart1(gpu);
    kApplyFirstShake(gpu);
    if (data.removeCM) {
        int step = (int) (context.getTime()/stepSize);
        if (step%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    }
    kVerletUpdatePart2(gpu);
}

CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
}

void CudaIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    initializeIntegration(system, data, integrator);

    // if LangevinIntegrator seed does not equal default value or is less than/equal to 0, then 
    // set gpu seed and redo random values

    if( integrator.getRandomNumberSeed() <= 1 ){
       _gpuContext* gpu = data.gpu;
       gpu->seed        = static_cast<unsigned long>( integrator.getRandomNumberSeed() );
       gpuInitializeRandoms( gpu );
    }
    prevStepSize = -1.0;
}

void CudaIntegrateLangevinStepKernel::execute(OpenMMContextImpl& context, const LangevinIntegrator& integrator) {
    _gpuContext* gpu = data.gpu;
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    if (temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        double tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        gpuSetIntegrationParameters(gpu, (float) tau, (float) stepSize, (float) temperature);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    kUpdatePart1(gpu);
    kApplyFirstShake(gpu);
    if (data.removeCM) {
        int step = (int) (context.getTime()/stepSize);
        if (step%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    }
    kUpdatePart2(gpu);
    kApplySecondShake(gpu);
}

CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
}

void CudaIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    initializeIntegration(system, data, integrator);
    prevStepSize = -1.0;
}

void CudaIntegrateBrownianStepKernel::execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator) {
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
    if (data.removeCM) {
        int step = (int) (context.getTime()/stepSize);
        if (step%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    }
    kBrownianUpdatePart2(gpu);
}

CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
}

void CudaApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    prevStepSize = -1.0;
}

void CudaApplyAndersenThermostatKernel::execute(OpenMMContextImpl& context) {
    _gpuContext* gpu = data.gpu;
    double temperature = context.getParameter(AndersenThermostat::Temperature());
    double frequency = context.getParameter(AndersenThermostat::CollisionFrequency());
    double stepSize = context.getIntegrator().getStepSize();
    if (temperature != prevTemp || frequency != prevFrequency || stepSize != prevStepSize) {
        // Initialize the GPU parameters.
        
        gpuSetAndersenThermostatParameters(gpu, (float) temperature, (float) (frequency*stepSize));
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

double CudaCalcKineticEnergyKernel::execute(OpenMMContextImpl& context) {
    // We don't currently have a GPU kernel to do this, so we retrieve the velocities and calculate the energy
    // on the CPU.
    
    const Stream& velocities = context.getVelocities();
    double* v = new double[velocities.getSize()*3];
    velocities.saveToArray(v);
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i)
        energy += masses[i]*(v[i*3]*v[i*3]+v[i*3+1]*v[i*3+1]+v[i*3+2]*v[i*3+2]);
    delete v;
    return 0.5*energy;
}

void CudaRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    data.removeCM = true;
    data.cmMotionFrequency = force.getFrequency();
}

void CudaRemoveCMMotionKernel::execute(OpenMMContextImpl& context) {
}
