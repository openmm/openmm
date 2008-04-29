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

#include "ReferenceKernels.h"
#include "ReferenceFloatStreamImpl.h"
#include "gbsa/CpuObc.h"
#include "SimTKReference/ReferenceAngleBondIxn.h"
#include "SimTKReference/ReferenceBondForce.h"
#include "SimTKReference/ReferenceHarmonicBondIxn.h"
#include "SimTKReference/ReferenceLJCoulomb14.h"
#include "SimTKReference/ReferenceLJCoulombIxn.h"
#include "SimTKReference/ReferenceProperDihedralBond.h"
#include "SimTKReference/ReferenceRbDihedralBond.h"
#include "SimTKReference/ReferenceStochasticDynamics.h"
#include "SimTKReference/ReferenceShakeAlgorithm.h"
#include <cmath>
#include <limits>

using namespace OpenMM;
using namespace std;

int** allocateIntArray(int length, int width) {
    int** array = new int*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new int[width];
    return array;
}

RealOpenMM** allocateRealArray(int length, int width) {
    RealOpenMM** array = new RealOpenMM*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new RealOpenMM[width];
    return array;
}

int** copyToArray(const vector<vector<int> > vec) {
    if (vec.size() == 0)
        return new int*[0];
    int** array = allocateIntArray(vec.size(), vec[0].size());
    for (int i = 0; i < vec.size(); ++i)
        for (int j = 0; j < vec[i].size(); ++j)
            array[i][j] = vec[i][j];
    return array;
}

RealOpenMM** copyToArray(const vector<vector<double> > vec) {
    if (vec.size() == 0)
        return new RealOpenMM*[0];
    RealOpenMM** array = allocateRealArray(vec.size(), vec[0].size());
    for (int i = 0; i < vec.size(); ++i)
        for (int j = 0; j < vec[i].size(); ++j)
            array[i][j] = vec[i][j];
    return array;
}

void disposeIntArray(int** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

void disposeRealArray(RealOpenMM** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

ReferenceCalcStandardMMForceFieldKernel::~ReferenceCalcStandardMMForceFieldKernel() {
    disposeIntArray(bondIndexArray, numBonds);
    disposeRealArray(bondParamArray, numBonds);
    disposeIntArray(angleIndexArray, numAngles);
    disposeRealArray(angleParamArray, numAngles);
    disposeIntArray(periodicTorsionIndexArray, numPeriodicTorsions);
    disposeRealArray(periodicTorsionParamArray, numPeriodicTorsions);
    disposeIntArray(rbTorsionIndexArray, numRBTorsions);
    disposeRealArray(rbTorsionParamArray, numRBTorsions);
    disposeRealArray(atomParamArray, numAtoms);
    disposeIntArray(exclusionArray, numAtoms);
    disposeIntArray(bonded14IndexArray, num14);
    disposeRealArray(bonded14ParamArray, num14);
}

void ReferenceCalcStandardMMForceFieldKernel::initialize(const vector<vector<int> >& bondIndices, const vector<vector<double> >& bondParameters,
        const vector<vector<int> >& angleIndices, const vector<vector<double> >& angleParameters,
        const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
        const vector<vector<int> >& rbTorsionIndices, const vector<vector<double> >& rbTorsionParameters,
        const vector<vector<int> >& bonded14Indices, double lj14Scale, double coulomb14Scale,
        const vector<set<int> >& exclusions, const vector<vector<double> >& nonbondedParameters) {
    numAtoms = nonbondedParameters.size();
    numBonds = bondIndices.size();
    numAngles = angleIndices.size();
    numPeriodicTorsions = periodicTorsionIndices.size();
    numRBTorsions = rbTorsionIndices.size();
    num14 = bonded14Indices.size();
    bondIndexArray = copyToArray(bondIndices);
    bondParamArray = copyToArray(bondParameters);
    angleIndexArray = copyToArray(angleIndices);
    angleParamArray = copyToArray(angleParameters);
    periodicTorsionIndexArray = copyToArray(periodicTorsionIndices);
    periodicTorsionParamArray = copyToArray(periodicTorsionParameters);
    rbTorsionIndexArray = copyToArray(rbTorsionIndices);
    rbTorsionParamArray = copyToArray(rbTorsionParameters);
    atomParamArray = allocateRealArray(numAtoms, 3);
    RealOpenMM sqrtEps = std::sqrt(138.935485);
    for (int i = 0; i < numAtoms; ++i) {
        atomParamArray[i][0] = 0.5*nonbondedParameters[i][1];
        atomParamArray[i][1] = 2.0*sqrt(nonbondedParameters[i][2]);
        atomParamArray[i][2] = nonbondedParameters[i][0]*sqrtEps;
    }
    exclusionArray = new int*[numAtoms];
    for (int i = 0; i < numAtoms; ++i) {
        exclusionArray[i] = new int[exclusions[i].size()+1];
        exclusionArray[i][0] = exclusions[i].size();
        int index = 0;
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            exclusionArray[i][++index] = *iter;
    }
    bonded14IndexArray = copyToArray(bonded14Indices);
    bonded14ParamArray = allocateRealArray(num14, 3);
    for (int i = 0; i < num14; ++i) {
        int atom1 = bonded14Indices[i][0];
        int atom2 = bonded14Indices[i][1];
        bonded14ParamArray[i][0] = atomParamArray[atom1][0]+atomParamArray[atom2][0];
        bonded14ParamArray[i][1] = lj14Scale*(atomParamArray[atom1][1]*atomParamArray[atom2][1]);
        bonded14ParamArray[i][2] = coulomb14Scale*(atomParamArray[atom1][2]*atomParamArray[atom2][2]);
    }
}

void ReferenceCalcStandardMMForceFieldKernel::executeForces(const Stream& positions, Stream& forces) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) positions.getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) forces.getImpl()).getData();
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, 0, 0, 0, harmonicBond);
    ReferenceAngleBondIxn angleBond;
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, 0, 0, 0, angleBond);
    ReferenceProperDihedralBond periodicTorsionBond;
    refBondForce.calculateForce(numPeriodicTorsions, periodicTorsionIndexArray, posData, periodicTorsionParamArray, forceData, 0, 0, 0, periodicTorsionBond);
    ReferenceRbDihedralBond rbTorsionBond;
    refBondForce.calculateForce(numRBTorsions, rbTorsionIndexArray, posData, rbTorsionParamArray, forceData, 0, 0, 0, rbTorsionBond);
    ReferenceLJCoulombIxn clj;
    clj.calculatePairIxn(numAtoms, posData, atomParamArray, exclusionArray, 0, forceData, 0, 0);
    ReferenceLJCoulomb14 nonbonded14;
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, 0, 0, 0, nonbonded14);
}

double ReferenceCalcStandardMMForceFieldKernel::executeEnergy(const Stream& positions) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) positions.getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = allocateRealArray(numAtoms, 3);
    int arraySize = max(max(max(max(numAtoms, numBonds), numAngles), numPeriodicTorsions), numRBTorsions);
    RealOpenMM* energyArray = new RealOpenMM[arraySize];
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    for (int i = 0; i < arraySize; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, energyArray, 0, &energy, harmonicBond);
    ReferenceAngleBondIxn angleBond;
    for (int i = 0; i < arraySize; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, energyArray, 0, &energy, angleBond);
    ReferenceProperDihedralBond periodicTorsionBond;
    for (int i = 0; i < arraySize; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numPeriodicTorsions, periodicTorsionIndexArray, posData, periodicTorsionParamArray, forceData, energyArray, 0, &energy, periodicTorsionBond);
    ReferenceRbDihedralBond rbTorsionBond;
    for (int i = 0; i < arraySize; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numRBTorsions, rbTorsionIndexArray, posData, rbTorsionParamArray, forceData, energyArray, 0, &energy, rbTorsionBond);
    ReferenceLJCoulombIxn clj;
    clj.calculatePairIxn(numAtoms, posData, atomParamArray, exclusionArray, 0, forceData, 0, &energy);
    ReferenceLJCoulomb14 nonbonded14;
    for (int i = 0; i < arraySize; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, energyArray, 0, &energy, nonbonded14);
    disposeRealArray(forceData, numAtoms);
    delete[] energyArray;
    return energy;
}

ReferenceCalcGBSAOBCForceFieldKernel::~ReferenceCalcGBSAOBCForceFieldKernel() {
    if (obc) {
        delete obc->getObcParameters();
        delete obc;
    }
}

void ReferenceCalcGBSAOBCForceFieldKernel::initialize(const vector<vector<double> >& atomParameters, double solventDielectric, double soluteDielectric) {
    int numAtoms = atomParameters.size();
    charges.resize(numAtoms);
    vector<RealOpenMM> atomicRadii(numAtoms);
    vector<RealOpenMM> scaleFactors(numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        charges[i] = atomParameters[i][0];
        atomicRadii[i] = atomParameters[i][1];
        scaleFactors[i] = atomParameters[i][2];
    }
    ObcParameters* obcParameters  = new ObcParameters(numAtoms, ObcParameters::ObcTypeII);
    obcParameters->setAtomicRadii(atomicRadii, SimTKOpenMMCommon::KcalAngUnits);
    obcParameters->setScaledRadiusFactors(scaleFactors);
    obcParameters->setSolventDielectric(solventDielectric);
    obcParameters->setSoluteDielectric(soluteDielectric);
    obc = new CpuObc(obcParameters);
    obc->setIncludeAceApproximation(true);
}

void ReferenceCalcGBSAOBCForceFieldKernel::executeForces(const Stream& positions, Stream& forces) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) positions.getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) forces.getImpl()).getData();
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 0);
}

double ReferenceCalcGBSAOBCForceFieldKernel::executeEnergy(const Stream& positions) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) positions.getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = allocateRealArray(positions.getSize(), 3);
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
    disposeRealArray(forceData, positions.getSize());
    return obc->getEnergy();
}

void ReferenceIntegrateVerletStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    
}

void ReferenceIntegrateVerletStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double stepSize) {
    
}
#include <iostream>
ReferenceIntegrateLangevinStepKernel::~ReferenceIntegrateLangevinStepKernel() {
    if (dynamics)
        delete dynamics;
    if (shake)
        delete shake;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (shakeParameters)
        disposeRealArray(shakeParameters, numConstraints);
}

void ReferenceIntegrateLangevinStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    this->masses = new RealOpenMM[masses.size()];
    for (int i = 0; i < masses.size(); ++i)
        this->masses[i] = masses[i];
    numConstraints = constraintIndices.size();
    this->constraintIndices = allocateIntArray(numConstraints, 2);
    for (int i = 0; i < numConstraints; ++i) {
        this->constraintIndices[i][0] = constraintIndices[i][0];
        this->constraintIndices[i][1] = constraintIndices[i][1];
    }
    shakeParameters = allocateRealArray(constraintLengths.size(), 1);
    for (int i = 0; i < constraintLengths.size(); ++i)
        shakeParameters[i][0] = constraintLengths[i];
}

void ReferenceIntegrateLangevinStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
    RealOpenMM** posData = ((ReferenceFloatStreamImpl&) positions.getImpl()).getData();
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) velocities.getImpl()).getData();
    RealOpenMM** forceData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) forces.getImpl()).getData()); // Reference code needs to be made const correct
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete shake;
        }
        RealOpenMM tau = (friction == 0.0 ? 0.0 : 1.0/friction);
        dynamics = new ReferenceStochasticDynamics(positions.getSize(), stepSize, tau, temperature);
        shake = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, shakeParameters);
        dynamics->setReferenceShakeAlgorithm(shake);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(positions.getSize(), posData, velData, forceData, masses);
}

void ReferenceIntegrateBrownianStepKernel::initialize(const vector<double>& masses, const vector<vector<int> >& constraintIndices,
        const vector<double>& constraintLengths) {
    
}

void ReferenceIntegrateBrownianStepKernel::execute(Stream& positions, Stream& velocities, const Stream& forces, double temperature, double friction, double stepSize) {
    
}

void ReferenceApplyAndersenThermostatKernel::initialize(const vector<double>& masses) {
    
}

void ReferenceApplyAndersenThermostatKernel::execute(Stream& velocities, double temperature, double collisionFrequency, double stepSize) {
    
}

void ReferenceCalcKineticEnergyKernel::initialize(const vector<double>& masses) {
    this->masses = masses;
}

double ReferenceCalcKineticEnergyKernel::execute(const Stream& velocities) {
    RealOpenMM** velData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) velocities.getImpl()).getData()); // Reference code needs to be made const correct
    double energy = 0.0;
    for (int i = 0; i < masses.size(); ++i)
        energy += masses[i]*(velData[i][0]*velData[i][0]+velData[i][1]*velData[i][1]+velData[i][2]*velData[i][2]);
    return 0.5*energy;
}
