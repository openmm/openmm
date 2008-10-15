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
#include "kernels/gputypes.h"
#include "kernels/cudaKernels.h"
#include <cmath>

extern "C" int gpuSetConstants( gpuContext gpu );

using namespace OpenMM;
using namespace std;

CudaCalcStandardMMForceFieldKernel::~CudaCalcStandardMMForceFieldKernel() {
}

void CudaCalcStandardMMForceFieldKernel::initialize(const vector<vector<int> >& bondIndices, const vector<vector<double> >& bondParameters,
        const vector<vector<int> >& angleIndices, const vector<vector<double> >& angleParameters,
        const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
        const vector<vector<int> >& rbTorsionIndices, const vector<vector<double> >& rbTorsionParameters,
        const vector<vector<int> >& bonded14Indices, double lj14Scale, double coulomb14Scale,
        const vector<set<int> >& exclusions, const vector<vector<double> >& nonbondedParameters,
        NonbondedMethod nonbondedMethod, double nonbondedCutoff, double periodicBoxSize[3]) {
    numAtoms = nonbondedParameters.size();
    numBonds = bondIndices.size();
    numAngles = angleIndices.size();
    numPeriodicTorsions = periodicTorsionIndices.size();
    numRBTorsions = rbTorsionIndices.size();
    num14 = bonded14Indices.size();
    const float RadiansToDegrees = 180.0/3.14159265;
    _gpuContext* gpu = data.gpu;
    
    // Initialize bonds.
    
    {
        vector<int> atom1(numBonds);
        vector<int> atom2(numBonds);
        vector<float> length(numBonds);
        vector<float> k(numBonds);
        for (int i = 0; i < numBonds; i++) {
            atom1[i] = bondIndices[i][0];
            atom2[i] = bondIndices[i][1];
            length[i] = (float) bondParameters[i][0];
            k[i] = (float) bondParameters[i][1];
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
            atom1[i] = angleIndices[i][0];
            atom2[i] = angleIndices[i][1];
            atom3[i] = angleIndices[i][2];
            angle[i] = (float) (angleParameters[i][0]*RadiansToDegrees);
            k[i] = (float) angleParameters[i][1];
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
            atom1[i] = periodicTorsionIndices[i][0];
            atom2[i] = periodicTorsionIndices[i][1];
            atom3[i] = periodicTorsionIndices[i][2];
            atom4[i] = periodicTorsionIndices[i][3];
            k[i] = (float) periodicTorsionParameters[i][0];
            phase[i] = (float) (periodicTorsionParameters[i][1]*RadiansToDegrees);
            periodicity[i] = (int) periodicTorsionParameters[i][2];
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
            atom1[i] = rbTorsionIndices[i][0];
            atom2[i] = rbTorsionIndices[i][1];
            atom3[i] = rbTorsionIndices[i][2];
            atom4[i] = rbTorsionIndices[i][3];
            c0[i] = (float) rbTorsionParameters[i][0];
            c1[i] = (float) rbTorsionParameters[i][1];
            c2[i] = (float) rbTorsionParameters[i][2];
            c3[i] = (float) rbTorsionParameters[i][3];
            c4[i] = (float) rbTorsionParameters[i][4];
            c5[i] = (float) rbTorsionParameters[i][5];
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
            atom[i] = i;
            q[i] = (float) nonbondedParameters[i][0];
            c6[i] = (float) (4*nonbondedParameters[i][2]*pow(nonbondedParameters[i][1], 6.0));
            c12[i] = (float) (4*nonbondedParameters[i][2]*pow(nonbondedParameters[i][1], 12.0));
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
            atom1[i] = bonded14Indices[i][0];
            atom2[i] = bonded14Indices[i][1];
            double sig = 0.5*(nonbondedParameters[atom1[i]][1]+nonbondedParameters[atom2[i]][1]);
            double eps = sqrt(nonbondedParameters[atom1[i]][2]*nonbondedParameters[atom2[i]][2]);
            c6[i] = (float) (4*eps*pow(sig, 6.0)*lj14Scale);
            c12[i] = (float) (4*eps*pow(sig, 12.0)*lj14Scale);
            q1[i] = (float) nonbondedParameters[atom1[i]][0];
            q2[i] = (float) nonbondedParameters[atom2[i]][0];
        }
        gpuSetLJ14Parameters(gpu, 138.935485f, (float) coulomb14Scale, atom1, atom2, c6, c12, q1, q2);
    }
}

void CudaCalcStandardMMForceFieldKernel::executeForces(const Stream& positions, Stream& forces) {
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

double CudaCalcStandardMMForceFieldKernel::executeEnergy(const Stream& positions) {
    // We don't currently have GPU kernels to calculate energy, so instead we have the reference
    // platform do it.  This is VERY slow.
    
    LangevinIntegrator integrator(0.0, 1.0, 0.0);
    ReferencePlatform platform;
    OpenMMContext context(system, integrator, platform);
    double* posData = new double[positions.getSize()*3];
    positions.saveToArray(posData);
    vector<Vec3> pos(positions.getSize());
    for (int i = 0; i < pos.size(); i++)
        pos[i] = Vec3(posData[3*i], posData[3*i+1], posData[3*i+2]);
    delete[] posData;
    context.setPositions(pos);
    return context.getState(State::Energy).getPotentialEnergy();
}

CudaCalcGBSAOBCForceFieldKernel::~CudaCalcGBSAOBCForceFieldKernel() {
}

void CudaCalcGBSAOBCForceFieldKernel::initialize(const vector<vector<double> >& atomParameters, double solventDielectric, double soluteDielectric) {
    int numAtoms = atomParameters.size();
    _gpuContext* gpu = data.gpu;
    vector<int> atom(numAtoms);
    vector<float> radius(numAtoms);
    vector<float> scale(numAtoms);
    for (int i = 0; i < numAtoms; i++) {
        atom[i] = i;
        radius[i] = (float) atomParameters[i][1];
        scale[i] = (float) atomParameters[i][2];
    }
    gpuSetObcParameters(gpu, soluteDielectric, solventDielectric, atom, radius, scale);
    data.useOBC = true;
}

void CudaCalcGBSAOBCForceFieldKernel::executeForces(const Stream& positions, Stream& forces) {
}

double CudaCalcGBSAOBCForceFieldKernel::executeEnergy(const Stream& positions) {
}

//CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
//}
//
//void CudaIntegrateVerletStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
//}
//
//void CudaIntegrateVerletStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double stepSize) {
//}

CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
}

void CudaIntegrateLangevinStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    
    // Set masses.
    
    _gpuContext* gpu = data.gpu;
    vector<float> mass(masses.size());
    for (int i = 0; i < (int) mass.size(); i++)
        mass[i] = (float) masses[i];
    gpuSetMass(gpu, mass);
    
    // Set constraints.
    
    int numConstraints = constraintLengths.size();
    vector<int> atom1(numConstraints);
    vector<int> atom2(numConstraints);
    vector<float> distance(numConstraints);
    vector<float> invMass1(numConstraints);
    vector<float> invMass2(numConstraints);
    for (int i = 0; i < numConstraints; i++) {
        atom1[i] = constraintIndices[i][0];
        atom2[i] = constraintIndices[i][1];
        distance[i] = (float) constraintLengths[i];
        invMass1[i] = 1.0f/mass[atom1[i]];
        invMass2[i] = 1.0f/mass[atom2[i]];
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

void CudaIntegrateLangevinStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
    _gpuContext* gpu = data.gpu;
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
    if (data.removeCM)
        gpu->bCalculateCM = true;
    kUpdatePart2(gpu);
    kApplySecondShake(gpu);
}
//
//CudaIntegrateBrownianStepKernel::~CudaIntegrateBrownianStepKernel() {
//}
//
//void CudaIntegrateBrownianStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
//        const vector<double>& constraintLengths) {
//}
//
//void CudaIntegrateBrownianStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
//}
//
//CudaApplyAndersenThermostatKernel::~CudaApplyAndersenThermostatKernel() {
//}
//
//void CudaApplyAndersenThermostatKernel::initialize(const vector<double>& masses) {
//}
//
//void CudaApplyAndersenThermostatKernel::execute(Stream& velocities, double temperature, double collisionFrequency, double stepSize) {
//}

void CudaCalcKineticEnergyKernel::initialize(const vector<double>& masses) {
    this->masses = masses;
}

double CudaCalcKineticEnergyKernel::execute(const Stream& velocities) {
    // We don't currently have a GPU kernel to do this, so we retrieve the velocities and calculate the energy
    // on the CPU.
    
    double* v = new double[velocities.getSize()*3];
    velocities.saveToArray(v);
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i)
        energy += masses[i]*(v[i*3]*v[i*3]+v[i*3+1]*v[i*3+1]+v[i*3+2]*v[i*3+2]);
    delete v;
    return 0.5*energy;
}

void CudaRemoveCMMotionKernel::initialize(const vector<double>& masses) {
    data.removeCM = true;
}

void CudaRemoveCMMotionKernel::execute(Stream& velocities) {
}
