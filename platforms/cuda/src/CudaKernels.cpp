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

CudaCalcStandardMMForceFieldKernel::~CudaCalcStandardMMForceFieldKernel() {
}

void CudaCalcStandardMMForceFieldKernel::initialize(const System& system, const StandardMMForceField& force, const std::vector<std::set<int> >& exclusions) {
    numAtoms = force.getNumAtoms();
    numBonds = force.getNumBonds();
    numAngles = force.getNumAngles();
    numPeriodicTorsions = force.getNumPeriodicTorsions();
    numRBTorsions = force.getNumRBTorsions();
    num14 = force.getNumNonbonded14();
    const float RadiansToDegrees = 180.0/3.14159265;
    _gpuContext* gpu = data.gpu;
    
    // Initialize bonds.
    
    {
        vector<int> atom1(numBonds);
        vector<int> atom2(numBonds);
        vector<float> length(numBonds);
        vector<float> k(numBonds);
        for (int i = 0; i < numBonds; i++) {
            double lengthValue, kValue;
            force.getBondParameters(i, atom1[i], atom2[i], lengthValue, kValue);
            length[i] = (float) lengthValue;
            k[i] = (float) kValue;
        }
        gpuSetBondParameters(gpu, atom1, atom2, length, k);
    }
    
    // Initialize angles.
    
    {
        vector<int> atom1(numAngles);
        vector<int> atom2(numAngles);
        vector<int> atom3(numAngles);
        vector<float> angle(numAngles);
        vector<float> k(numAngles);
        for (int i = 0; i < numAngles; i++) {
            double angleValue, kValue;
            force.getAngleParameters(i, atom1[i], atom2[i], atom3[i], angleValue, kValue);
            angle[i] = (float) (angleValue*RadiansToDegrees);
            k[i] = (float) kValue;
        }
        gpuSetBondAngleParameters(gpu, atom1, atom2, atom3, angle, k);
    }

    // Initialize periodic torsions.
    
    {
        vector<int> atom1(numPeriodicTorsions);
        vector<int> atom2(numPeriodicTorsions);
        vector<int> atom3(numPeriodicTorsions);
        vector<int> atom4(numPeriodicTorsions);
        vector<float> k(numPeriodicTorsions);
        vector<float> phase(numPeriodicTorsions);
        vector<int> periodicity(numPeriodicTorsions);
        for (int i = 0; i < numPeriodicTorsions; i++) {
            double kValue, phaseValue;
            force.getPeriodicTorsionParameters(i, atom1[i], atom2[i], atom3[i], atom4[i], periodicity[i], phaseValue, kValue);
            k[i] = (float) kValue;
            phase[i] = (float) (phaseValue*RadiansToDegrees);
        }
        gpuSetDihedralParameters(gpu, atom1, atom2, atom3, atom4, k, phase, periodicity);
    }
    
    // Initialize Ryckaert-Bellemans torsions.
    
    {
        vector<int> atom1(numRBTorsions);
        vector<int> atom2(numRBTorsions);
        vector<int> atom3(numRBTorsions);
        vector<int> atom4(numRBTorsions);
        vector<float> c0(numRBTorsions);
        vector<float> c1(numRBTorsions);
        vector<float> c2(numRBTorsions);
        vector<float> c3(numRBTorsions);
        vector<float> c4(numRBTorsions);
        vector<float> c5(numRBTorsions);
        for (int i = 0; i < numRBTorsions; i++) {
            double c[6];
            force.getRBTorsionParameters(i, atom1[i], atom2[i], atom3[i], atom4[i], c[0], c[1], c[2], c[3], c[4], c[5]);
            c0[i] = (float) c[0];
            c1[i] = (float) c[1];
            c2[i] = (float) c[2];
            c3[i] = (float) c[3];
            c4[i] = (float) c[4];
            c5[i] = (float) c[5];
        }
        gpuSetRbDihedralParameters(gpu, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
    }
    
    // Initialize nonbonded interactions.
    
    {
        vector<int> atom(numAtoms);
        vector<float> c6(numAtoms);
        vector<float> c12(numAtoms);
        vector<float> q(numAtoms);
        vector<char> symbol;
        vector<vector<int> > exclusionList(numAtoms);
        for (int i = 0; i < numAtoms; i++) {
            double charge, radius, depth;
            force.getAtomParameters(i, charge, radius, depth);
            atom[i] = i;
            q[i] = (float) charge;
            c6[i] = (float) (4*depth*pow(radius, 6.0));
            c12[i] = (float) (4*depth*pow(radius, 12.0));
            exclusionList[i] = vector<int>(exclusions[i].begin(), exclusions[i].end());
            exclusionList[i].push_back(i);
        }
        gpuSetCoulombParameters(gpu, 138.935485f, atom, c6, c12, q, symbol, exclusionList);
    }

    // Initialize 1-4 nonbonded interactions.
    
    {
        vector<int> atom1(num14);
        vector<int> atom2(num14);
        vector<float> c6(num14);
        vector<float> c12(num14);
        vector<float> q1(num14);
        vector<float> q2(num14);
        for (int i = 0; i < num14; i++) {
            double charge, sig, eps;
            force.getNonbonded14Parameters(i, atom1[i], atom2[i], charge, sig, eps);
            c6[i] = (float) (4*eps*pow(sig, 6.0));
            c12[i] = (float) (4*eps*pow(sig, 12.0));
            float q = (float) std::sqrt(charge);
            q1[i] = q;
            q2[i] = q;
        }
        gpuSetLJ14Parameters(gpu, 138.935485f, 1.0f, atom1, atom2, c6, c12, q1, q2);
    }
}

void CudaCalcStandardMMForceFieldKernel::executeForces(OpenMMContextImpl& context) {
    _gpuContext* gpu = data.gpu;
    if (data.useOBC) {
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

double CudaCalcStandardMMForceFieldKernel::executeEnergy(OpenMMContextImpl& context) {
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

CudaCalcGBSAOBCForceFieldKernel::~CudaCalcGBSAOBCForceFieldKernel() {
}

void CudaCalcGBSAOBCForceFieldKernel::initialize(const System& system, const GBSAOBCForceField& force) {
    int numAtoms = system.getNumAtoms();
    _gpuContext* gpu = data.gpu;
    vector<int> atom(numAtoms);
    vector<float> radius(numAtoms);
    vector<float> scale(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        double charge, atomRadius, scalingFactor;
        force.getAtomParameters(i, charge, atomRadius, scalingFactor);
        atom[i] = i;
        radius[i] = (float) atomRadius;
        scale[i] = (float) scalingFactor;
    }
    gpuSetObcParameters(gpu, force.getSoluteDielectric(), force.getSolventDielectric(), atom, radius, scale);
    data.useOBC = true;
}

void CudaCalcGBSAOBCForceFieldKernel::executeForces(OpenMMContextImpl& context) {
}

double CudaCalcGBSAOBCForceFieldKernel::executeEnergy(OpenMMContextImpl& context) {
}

//CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
//}
//
//void CudaIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
//}
//
//void CudaIntegrateVerletStepKernel::execute(OpenMMContextImpl& context, const VerletIntegrator& integrator) {
//}

CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
}

void CudaIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    
    // Set masses.
    
    _gpuContext* gpu = data.gpu;
    int numAtoms = system.getNumAtoms();
    vector<float> mass(numAtoms);
    for (int i = 0; i < numAtoms; i++)
        mass[i] = (float) system.getAtomMass(i);
    gpuSetMass(gpu, mass);
    
    // Set constraints.
    
    int numConstraints = system.getNumConstraints();
    vector<int> atom1(numConstraints);
    vector<int> atom2(numConstraints);
    vector<float> distance(numConstraints);
    vector<float> invMass1(numConstraints);
    vector<float> invMass2(numConstraints);
    for (int i = 0; i < numConstraints; i++) {
        int atom1Index, atom2Index;
        double constraintDistance;
        system.getConstraintParameters(i, atom1Index, atom2Index, constraintDistance);
        atom1[i] = atom1Index;
        atom2[i] = atom2Index;
        distance[i] = (float) constraintDistance;
        invMass1[i] = 1.0f/mass[atom1Index];
        invMass2[i] = 1.0f/mass[atom2Index];
    }
    gpuSetShakeParameters(gpu, atom1, atom2, distance, invMass1, invMass2);
    gpuBuildThreadBlockWorkList(gpu);
    gpuBuildExclusionList(gpu);
    gpuBuildOutputBuffers(gpu);
    gpuSetConstants(gpu);
    kCalculateObcGbsaBornSum(gpu);
    kReduceObcGbsaBornSum(gpu);
    kClearBornForces(gpu);
    kClearForces(gpu);
    cudaThreadSynchronize();
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
        gpuSetIntegrationParameters(gpu, tau, stepSize, temperature);
        gpuSetConstants(gpu);
        kGenerateRandoms(gpu);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    kUpdatePart1(gpu);
    kApplyFirstShake(gpu);
    if (data.removeCM) {
        int step = context.getTime()/stepSize;
        if (step%data.cmMotionFrequency == 0)
            gpu->bCalculateCM = true;
    }
    kUpdatePart2(gpu);
    kApplySecondShake(gpu);
}
//
//CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
//}
//
//void CudaIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
//}
//
//void CudaIntegrateBrownianStepKernel::execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator) {
//}
//
//CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
//}
//
//void CudaApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
//}
//
//void CudaApplyAndersenThermostatKernel::execute(OpenMMContextImpl& context) {
//}

void CudaCalcKineticEnergyKernel::initialize(const System& system) {
    int numAtoms = system.getNumAtoms();
    masses.resize(numAtoms);
    for (size_t i = 0; i < numAtoms; ++i)
        masses[i] = system.getAtomMass(i);
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
