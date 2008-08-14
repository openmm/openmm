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

#include "CudaKernels.h"
#include "CudaStreamImpl.h"
#include "kernels/gputypes.h"
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
    
    // Initialize bonds.
    
    gpu->sim.bonds                      = numBonds;
    CUDAStream<int4>* psBondID          = new CUDAStream<int4>(numBonds, 1);
    gpu->psBondID                       = psBondID;
    gpu->sim.pBondID                    = psBondID->_pDevStream[0];
    CUDAStream<float2>* psBondParameter = new CUDAStream<float2>(numBonds, 1);
    gpu->psBondParameter                = psBondParameter;
    gpu->sim.pBondParameter             = psBondParameter->_pDevStream[0];
    for (int i = 0; i < numBonds; i++ ) {
        psBondID->_pSysStream[0][i].x        = bondIndices[i][0];
        psBondID->_pSysStream[0][i].y        = bondIndices[i][1];
        psBondID->_pSysStream[0][i].z        = gpu->pOutputBufferCounter[psBondID->_pSysStream[0][i].x]++;
        psBondID->_pSysStream[0][i].w        = gpu->pOutputBufferCounter[psBondID->_pSysStream[0][i].y]++;
        psBondParameter->_pSysStream[0][i].x = bondParameters[i][0];
        psBondParameter->_pSysStream[0][i].y = bondParameters[i][1];
    }
    psBondID->Upload();
    psBondParameter->Upload();
    
    // Initialize angles.
    
    gpu->sim.bond_angles                     = numAngles;
    CUDAStream<int4>* psBondAngleID1         = new CUDAStream<int4>(numAngles, 1);
    gpu->psBondAngleID1                      = psBondAngleID1;
    gpu->sim.pBondAngleID1                   = psBondAngleID1->_pDevStream[0];
    CUDAStream<int2>* psBondAngleID2         = new CUDAStream<int2>(numAngles, 1);
    gpu->psBondAngleID2                      = psBondAngleID2;
    gpu->sim.pBondAngleID2                   = psBondAngleID2->_pDevStream[0];
    CUDAStream<float2>* psBondAngleParameter = new CUDAStream<float2>(numAngles, 1);
    gpu->psBondAngleParameter                = psBondAngleParameter;
    gpu->sim.pBondAngleParameter             = psBondAngleParameter->_pDevStream[0];
    for (int i = 0; i < numAngles; i++) {
        psBondAngleID1->_pSysStream[0][i].x         = angleIndices[i][0];
        psBondAngleID1->_pSysStream[0][i].y         = angleIndices[i][1];
        psBondAngleID1->_pSysStream[0][i].z         = angleIndices[i][2];
        psBondAngleID1->_pSysStream[0][i].w         = gpu->pOutputBufferCounter[psBondAngleID1->_pSysStream[0][i].x]++;
        psBondAngleID2->_pSysStream[0][i].x         = gpu->pOutputBufferCounter[psBondAngleID1->_pSysStream[0][i].y]++;
        psBondAngleID2->_pSysStream[0][i].y         = gpu->pOutputBufferCounter[psBondAngleID1->_pSysStream[0][i].z]++;
        psBondAngleParameter->_pSysStream[0][i].x   = angleParameters[i][0]*180.0/M_PI;
        psBondAngleParameter->_pSysStream[0][i].y   = angleParameters[i][1];
    }
    psBondAngleID1->Upload();
    psBondAngleID2->Upload();
    psBondAngleParameter->Upload();

    // Initialize periodic torsions.
    
    gpu->sim.dihedrals = numPeriodicTorsions;
    CUDAStream<int4>* psDihedralID1             = new CUDAStream<int4>(numPeriodicTorsions, 1);
    gpu->psDihedralID1                          = psDihedralID1;
    gpu->sim.pDihedralID1                       = psDihedralID1->_pDevStream[0];
    CUDAStream<int4>* psDihedralID2             = new CUDAStream<int4>(numPeriodicTorsions, 1);
    gpu->psDihedralID2                          = psDihedralID2;
    gpu->sim.pDihedralID2                       = psDihedralID2->_pDevStream[0];
    CUDAStream<float4>* psDihedralParameter     = new CUDAStream<float4>(numPeriodicTorsions, 1);
    gpu->psDihedralParameter                    = psDihedralParameter;
    gpu->sim.pDihedralParameter                 = psDihedralParameter->_pDevStream[0];
    for (int i = 0; i < numPeriodicTorsions; i++) {
        psDihedralID1->_pSysStream[0][i].x              = periodicTorsionIndices[i][0];
        psDihedralID1->_pSysStream[0][i].y              = periodicTorsionIndices[i][1];
        psDihedralID1->_pSysStream[0][i].z              = periodicTorsionIndices[i][2];
        psDihedralID1->_pSysStream[0][i].w              = periodicTorsionIndices[i][3];
        psDihedralID2->_pSysStream[0][i].x              = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].x]++;
        psDihedralID2->_pSysStream[0][i].y              = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].y]++;
        psDihedralID2->_pSysStream[0][i].z              = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].z]++;
        psDihedralID2->_pSysStream[0][i].w              = gpu->pOutputBufferCounter[psDihedralID1->_pSysStream[0][i].w]++;
        psDihedralParameter->_pSysStream[0][i].x        = periodicTorsionParameters[i][0];
        psDihedralParameter->_pSysStream[0][i].y        = periodicTorsionParameters[i][1];
        psDihedralParameter->_pSysStream[0][i].z        = periodicTorsionParameters[i][2];
        psDihedralParameter->_pSysStream[0][i].w        = 0.0f;
    }
    psDihedralID1->Upload();
    psDihedralID2->Upload();
    psDihedralParameter->Upload();
    
    // Initialize Ryckaert-Bellemans torsions.
    
    gpu->sim.rb_dihedrals = numRBTorsions;
    CUDAStream<int4>* psRbDihedralID1           = new CUDAStream<int4>(numRBTorsions, 1);
    gpu->psRbDihedralID1                        = psRbDihedralID1;
    gpu->sim.pRbDihedralID1                     = psRbDihedralID1->_pDevStream[0];
    CUDAStream<int4>* psRbDihedralID2           = new CUDAStream<int4>(numRBTorsions, 1);
    gpu->psRbDihedralID2                        = psRbDihedralID2;
    gpu->sim.pRbDihedralID2                     = psRbDihedralID2->_pDevStream[0];
    CUDAStream<float4>* psRbDihedralParameter1  = new CUDAStream<float4>(numRBTorsions, 1);
    gpu->psRbDihedralParameter1                 = psRbDihedralParameter1;
    gpu->sim.pRbDihedralParameter1              = psRbDihedralParameter1->_pDevStream[0];
    CUDAStream<float2>* psRbDihedralParameter2  = new CUDAStream<float2>(numRBTorsions, 1);	
    gpu->psRbDihedralParameter2                 = psRbDihedralParameter2;
    gpu->sim.pRbDihedralParameter2              = psRbDihedralParameter2->_pDevStream[0];
    for (int i = 0; i < numRBTorsions; i++) {
        psRbDihedralID1->_pSysStream[0][i].x            = rbTorsionIndices[i][0];
        psRbDihedralID1->_pSysStream[0][i].y            = rbTorsionIndices[i][1];
        psRbDihedralID1->_pSysStream[0][i].z            = rbTorsionIndices[i][2];
        psRbDihedralID1->_pSysStream[0][i].w            = rbTorsionIndices[i][3];
        psRbDihedralID2->_pSysStream[0][i].x            = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].x]++;
        psRbDihedralID2->_pSysStream[0][i].y            = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].y]++;
        psRbDihedralID2->_pSysStream[0][i].z            = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].z]++;
        psRbDihedralID2->_pSysStream[0][i].w            = gpu->pOutputBufferCounter[psRbDihedralID1->_pSysStream[0][i].w]++;
        psRbDihedralParameter1->_pSysStream[0][i].x     = rbTorsionParameters[i][0];
        psRbDihedralParameter1->_pSysStream[0][i].y     = rbTorsionParameters[i][1];
        psRbDihedralParameter1->_pSysStream[0][i].z     = rbTorsionParameters[i][2];
        psRbDihedralParameter1->_pSysStream[0][i].w     = rbTorsionParameters[i][3];
        psRbDihedralParameter2->_pSysStream[0][i].x     = rbTorsionParameters[i][4];
        psRbDihedralParameter2->_pSysStream[0][i].y     = rbTorsionParameters[i][5];
    }
    psRbDihedralID1->Upload();
    psRbDihedralID2->Upload();
    psRbDihedralParameter1->Upload();
    psRbDihedralParameter2->Upload();
    
    // Initialize nonbonded interactions.
    
    for (int i = 0; i < numAtoms; i++) {
        gpu->psPosq4->_pSysStream[0][i].w = nonbondedParameters[i][0];
        gpu->psSigEps2->_pSysStream[0][i].x = nonbondedParameters[i][1];
        gpu->psSigEps2->_pSysStream[0][i].y = nonbondedParameters[i][2];
    }
    gpu->psPosq4->Upload();
    gpu->psSigEps2->Upload();

    // Initialize 1-4 nonbonded interactions.
    
    gpu->sim.LJ14s                              = num14;
    CUDAStream<int4>* psLJ14ID                  = new CUDAStream<int4>(num14, 1);
    gpu->psLJ14ID                               = psLJ14ID;
    gpu->sim.pLJ14ID                            = psLJ14ID->_pDevStream[0];
    CUDAStream<float4>* psLJ14Parameter         = new CUDAStream<float4>(num14, 1);
    gpu->psLJ14Parameter                        = psLJ14Parameter;
    gpu->sim.pLJ14Parameter                     = psLJ14Parameter->_pDevStream[0];
    double sqrtEps = std::sqrt(138.935485);
    for (int i = 0; i < num14; i++) {
        int atom1 = bonded14Indices[i][0];
        int atom2 = bonded14Indices[i][1];
        double atom1params[] = {0.5*nonbondedParameters[atom1][1], 2.0*sqrt(nonbondedParameters[atom1][2]), nonbondedParameters[atom1][0]*sqrtEps};
        double atom2params[] = {0.5*nonbondedParameters[atom2][1], 2.0*sqrt(nonbondedParameters[atom2][2]), nonbondedParameters[atom2][0]*sqrtEps};
        psLJ14ID->_pSysStream[0][i].x          = atom1;
        psLJ14ID->_pSysStream[0][i].y          = atom2;
        psLJ14ID->_pSysStream[0][i].z          = gpu->pOutputBufferCounter[psLJ14ID->_pSysStream[0][i].x]++;
        psLJ14ID->_pSysStream[0][i].w          = gpu->pOutputBufferCounter[psLJ14ID->_pSysStream[0][i].y]++;
        psLJ14Parameter->_pSysStream[0][i].x   = atom1params[0]+atom2params[0];
        psLJ14Parameter->_pSysStream[0][i].y   = lj14Scale*(atom1params[1]*atom2params[1]);
        psLJ14Parameter->_pSysStream[0][i].z   = coulomb14Scale*(atom1params[2]*atom2params[2]);
    }
    psLJ14ID->Upload();
    psLJ14Parameter->Upload();
    
    // Initialize exclusions.
    
    // TODO
    
    // Finish initialization.
    
    gpuSetConstants(gpu);
}

void CudaCalcStandardMMForceFieldKernel::executeForces(const Stream& positions, Stream& forces) {
}

double CudaCalcStandardMMForceFieldKernel::executeEnergy(const Stream& positions) {
}

//CudaCalcGBSAOBCForceFieldKernel::~CudaCalcGBSAOBCForceFieldKernel() {
//}
//
//void CudaCalcGBSAOBCForceFieldKernel::initialize(const vector<vector<double> >& atomParameters, double solventDielectric, double soluteDielectric) {
//}
//
//void CudaCalcGBSAOBCForceFieldKernel::executeForces(const Stream& positions, Stream& forces) {
//}
//
//double CudaCalcGBSAOBCForceFieldKernel::executeEnergy(const Stream& positions) {
//}
//
//CudaIntegrateVerletStepKernel::~CudaIntegrateVerletStepKernel() {
//}
//
//void CudaIntegrateVerletStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
//}
//
//void CudaIntegrateVerletStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double stepSize) {
//}
//
//CudaIntegrateLangevinStepKernel::~CudaIntegrateLangevinStepKernel() {
//}
//
//void CudaIntegrateLangevinStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
//        const vector<double>& constraintLengths) {
//}
//
//void CudaIntegrateLangevinStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
//}
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
//
//void CudaCalcKineticEnergyKernel::initialize(const vector<double>& masses) {
//}
//
//double CudaCalcKineticEnergyKernel::execute(const Stream& velocities) {
//}
//
//void CudaRemoveCMMotionKernel::initialize(const vector<double>& masses) {
//}
//
//void CudaRemoveCMMotionKernel::execute(Stream& velocities) {
//}
