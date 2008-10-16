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

#include "System.h"


#include "ReferenceKernels.h"
#include "ReferenceFloatStreamImpl.h"
#include "gbsa/CpuObc.h"
#include "SimTKReference/ReferenceAndersenThermostat.h"
#include "SimTKReference/ReferenceAngleBondIxn.h"
#include "SimTKReference/ReferenceBondForce.h"
#include "SimTKReference/ReferenceBrownianDynamics.h"
#include "SimTKReference/ReferenceHarmonicBondIxn.h"
#include "SimTKReference/ReferenceLJCoulomb14.h"
#include "SimTKReference/ReferenceLJCoulombIxn.h"
#include "SimTKReference/ReferenceProperDihedralBond.h"
#include "SimTKReference/ReferenceRbDihedralBond.h"
#include "SimTKReference/ReferenceStochasticDynamics.h"
#include "SimTKReference/ReferenceShakeAlgorithm.h"
#include "SimTKReference/ReferenceVerletDynamics.h"
#include "CMMotionRemover.h"
#include "System.h"
#include "internal/OpenMMContextImpl.h"
#include "Integrator.h"
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
    for (size_t i = 0; i < vec.size(); ++i)
        for (size_t j = 0; j < vec[i].size(); ++j)
            array[i][j] = vec[i][j];
    return array;
}

RealOpenMM** copyToArray(const vector<vector<double> > vec) {
    if (vec.size() == 0)
        return new RealOpenMM*[0];
    RealOpenMM** array = allocateRealArray(vec.size(), vec[0].size());
    for (size_t i = 0; i < vec.size(); ++i)
        for (size_t j = 0; j < vec[i].size(); ++j)
            array[i][j] = static_cast<RealOpenMM>(vec[i][j]);
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
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceCalcStandardMMForceFieldKernel::initialize(const System& system, const StandardMMForceField& force, const std::vector<std::set<int> >& exclusions) {
    numAtoms = force.getNumAtoms();
    numBonds = force.getNumBonds();
    numAngles = force.getNumAngles();
    numPeriodicTorsions = force.getNumPeriodicTorsions();
    numRBTorsions = force.getNumRBTorsions();
    num14 = force.getNumNonbonded14();
    bondIndexArray = allocateIntArray(numBonds, 2);
    bondParamArray = allocateRealArray(numBonds, 2);
    angleIndexArray = allocateIntArray(numAngles, 3);
    angleParamArray = allocateRealArray(numAngles, 2);
    periodicTorsionIndexArray = allocateIntArray(numPeriodicTorsions, 4);
    periodicTorsionParamArray = allocateRealArray(numPeriodicTorsions, 3);
    rbTorsionIndexArray = allocateIntArray(numRBTorsions, 4);
    rbTorsionParamArray = allocateRealArray(numRBTorsions, 6);
    bonded14IndexArray = allocateIntArray(num14, 2);
    bonded14ParamArray = allocateRealArray(num14, 3);
    atomParamArray = allocateRealArray(numAtoms, 3);
    for (int i = 0; i < force.getNumBonds(); ++i) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(i, atom1, atom2, length, k);
        bondIndexArray[i][0] = atom1;
        bondIndexArray[i][1] = atom2;
        bondParamArray[i][0] = length;
        bondParamArray[i][1] = k;
    }
    for (int i = 0; i < force.getNumAngles(); ++i) {
        int atom1, atom2, atom3;
        double angle, k;
        force.getAngleParameters(i, atom1, atom2, atom3, angle, k);
        angleIndexArray[i][0] = atom1;
        angleIndexArray[i][1] = atom2;
        angleIndexArray[i][2] = atom3;
        angleParamArray[i][0] = angle;
        angleParamArray[i][1] = k;
    }
    for (int i = 0; i < force.getNumPeriodicTorsions(); ++i) {
        int atom1, atom2, atom3, atom4, periodicity;
        double phase, k;
        force.getPeriodicTorsionParameters(i, atom1, atom2, atom3, atom4, periodicity, phase, k);
        periodicTorsionIndexArray[i][0] = atom1;
        periodicTorsionIndexArray[i][1] = atom2;
        periodicTorsionIndexArray[i][2] = atom3;
        periodicTorsionIndexArray[i][3] = atom4;
        periodicTorsionParamArray[i][0] = k;
        periodicTorsionParamArray[i][1] = phase;
        periodicTorsionParamArray[i][2] = periodicity;
    }
    for (int i = 0; i < force.getNumRBTorsions(); ++i) {
        int atom1, atom2, atom3, atom4;
        double c0, c1, c2, c3, c4, c5;
        force.getRBTorsionParameters(i, atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5);
        rbTorsionIndexArray[i][0] = atom1;
        rbTorsionIndexArray[i][1] = atom2;
        rbTorsionIndexArray[i][2] = atom3;
        rbTorsionIndexArray[i][3] = atom4;
        rbTorsionParamArray[i][0] = c0;
        rbTorsionParamArray[i][1] = c1;
        rbTorsionParamArray[i][2] = c2;
        rbTorsionParamArray[i][3] = c3;
        rbTorsionParamArray[i][4] = c4;
        rbTorsionParamArray[i][5] = c5;
    }
    RealOpenMM sqrtEps = static_cast<RealOpenMM>( std::sqrt(138.935485) );
    for (int i = 0; i < numAtoms; ++i) {
        double charge, radius, depth;
        force.getAtomParameters(i, charge, radius, depth);
        atomParamArray[i][0] = static_cast<RealOpenMM>(0.5*radius);
        atomParamArray[i][1] = static_cast<RealOpenMM>(2.0*sqrt(depth));
        atomParamArray[i][2] = static_cast<RealOpenMM>(charge*sqrtEps);
    }
    this->exclusions = exclusions;
    exclusionArray = new int*[numAtoms];
    for (int i = 0; i < numAtoms; ++i) {
        exclusionArray[i] = new int[exclusions[i].size()+1];
        exclusionArray[i][0] = exclusions[i].size();
        int index = 0;
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            exclusionArray[i][++index] = *iter;
    }
    for (int i = 0; i < num14; ++i) {
        int atom1, atom2;
        double charge, radius, depth;
        force.getNonbonded14Parameters(i, atom1, atom2, charge, radius, depth);
        bonded14IndexArray[i][0] = atom1;
        bonded14IndexArray[i][1] = atom2;
        bonded14ParamArray[i][0] = static_cast<RealOpenMM>(radius);
        bonded14ParamArray[i][1] = static_cast<RealOpenMM>(4.0*depth);
        bonded14ParamArray[i][2] = static_cast<RealOpenMM>(charge*sqrtEps*sqrtEps);
    }
    nonbondedMethod = CalcStandardMMForceFieldKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    double boxSize[3];
    force.getPeriodicBoxSize(boxSize[0], boxSize[1], boxSize[2]);
    periodicBoxSize[0] = (RealOpenMM) boxSize[0];
    periodicBoxSize[1] = (RealOpenMM) boxSize[1];
    periodicBoxSize[2] = (RealOpenMM) boxSize[2];
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();
        
}

void ReferenceCalcStandardMMForceFieldKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
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
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numAtoms, posData, exclusions, periodic ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, 78.3f);
    }
    if (periodic)
        clj.setPeriodic(periodicBoxSize);
    clj.calculatePairIxn(numAtoms, posData, atomParamArray, exclusionArray, 0, forceData, 0, 0);
    ReferenceLJCoulomb14 nonbonded14;
    if (nonbondedMethod != NoCutoff)
        nonbonded14.setUseCutoff(nonbondedCutoff, 78.3f);
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, 0, 0, 0, nonbonded14);
}

double ReferenceCalcStandardMMForceFieldKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
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
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numAtoms, posData, exclusions, periodic ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, 78.3f);
    }
    if (periodic)
        clj.setPeriodic(periodicBoxSize);
    clj.calculatePairIxn(numAtoms, posData, atomParamArray, exclusionArray, 0, forceData, 0, &energy);
    ReferenceLJCoulomb14 nonbonded14;
    if (nonbondedMethod != NoCutoff)
        nonbonded14.setUseCutoff(nonbondedCutoff, 78.3f);
    for (int i = 0; i < arraySize; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, energyArray, 0, &energy, nonbonded14);
    disposeRealArray(forceData, numAtoms);
    delete[] energyArray;
    return energy;
}

ReferenceCalcGBSAOBCForceFieldKernel::~ReferenceCalcGBSAOBCForceFieldKernel() {
    if (obc) {
        // delete obc->getObcParameters();
        delete obc;
    }
}

void ReferenceCalcGBSAOBCForceFieldKernel::initialize(const System& system, const GBSAOBCForceField& force) {
    int numAtoms = system.getNumAtoms();
    charges.resize(numAtoms);
    vector<RealOpenMM> atomicRadii(numAtoms);
    vector<RealOpenMM> scaleFactors(numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        double charge, radius, scalingFactor;
        force.getAtomParameters(i, charge, radius, scalingFactor);
        charges[i] = static_cast<RealOpenMM>(charge);
        atomicRadii[i] = static_cast<RealOpenMM>(radius);
        scaleFactors[i] = static_cast<RealOpenMM>(scalingFactor);
    }
    ObcParameters* obcParameters  = new ObcParameters(numAtoms, ObcParameters::ObcTypeII);
    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
    obcParameters->setSolventDielectric( static_cast<RealOpenMM>(force.getSolventDielectric()) );
    obcParameters->setSoluteDielectric( static_cast<RealOpenMM>(force.getSoluteDielectric()) );
    obc = new CpuObc(obcParameters);
    obc->setIncludeAceApproximation(true);
}

void ReferenceCalcGBSAOBCForceFieldKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 0);
}

double ReferenceCalcGBSAOBCForceFieldKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumAtoms(), 3);
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
    disposeRealArray(forceData, context.getSystem().getNumAtoms());
    return obc->getEnergy();
}

ReferenceIntegrateVerletStepKernel::~ReferenceIntegrateVerletStepKernel() {
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

void ReferenceIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    int numAtoms = system.getNumAtoms();
    masses = new RealOpenMM[numAtoms];
    for (size_t i = 0; i < numAtoms; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getAtomMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    shakeParameters = allocateRealArray(numConstraints, 1);
    for (int i = 0; i < numConstraints; ++i) {
        int atom1, atom2;
        double distance;
        system.getConstraintParameters(i, atom1, atom2, distance);
        constraintIndices[i][0] = atom1;
        constraintIndices[i][1] = atom2;
        shakeParameters[i][0] = static_cast<RealOpenMM>(distance);
    }
}

void ReferenceIntegrateVerletStepKernel::execute(OpenMMContextImpl& context, const VerletIntegrator& integrator) {
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = ((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData();
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    RealOpenMM** forceData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData()); // Reference code needs to be made const correct
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete shake;
        }
        dynamics = new ReferenceVerletDynamics(context.getSystem().getNumAtoms(), static_cast<RealOpenMM>(stepSize) );
        shake = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, shakeParameters);
        dynamics->setReferenceShakeAlgorithm(shake);
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumAtoms(), posData, velData, forceData, masses);
}

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

void ReferenceIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    int numAtoms = system.getNumAtoms();
    masses = new RealOpenMM[numAtoms];
    for (size_t i = 0; i < numAtoms; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getAtomMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    shakeParameters = allocateRealArray(numConstraints, 1);
    for (int i = 0; i < numConstraints; ++i) {
        int atom1, atom2;
        double distance;
        system.getConstraintParameters(i, atom1, atom2, distance);
        constraintIndices[i][0] = atom1;
        constraintIndices[i][1] = atom2;
        shakeParameters[i][0] = static_cast<RealOpenMM>(distance);
    }
}

void ReferenceIntegrateLangevinStepKernel::execute(OpenMMContextImpl& context, const LangevinIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = ((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData();
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    RealOpenMM** forceData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData()); // Reference code needs to be made const correct
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete shake;
        }
        RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );
        dynamics = new ReferenceStochasticDynamics(
				context.getSystem().getNumAtoms(), 
				static_cast<RealOpenMM>(stepSize), 
				static_cast<RealOpenMM>(tau), 
				static_cast<RealOpenMM>(temperature) );
        shake = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, shakeParameters);
        dynamics->setReferenceShakeAlgorithm(shake);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumAtoms(), posData, velData, forceData, masses);
}

ReferenceIntegrateBrownianStepKernel::~ReferenceIntegrateBrownianStepKernel() {
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

void ReferenceIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    int numAtoms = system.getNumAtoms();
    masses = new RealOpenMM[numAtoms];
    for (size_t i = 0; i < numAtoms; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getAtomMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    shakeParameters = allocateRealArray(numConstraints, 1);
    for (int i = 0; i < numConstraints; ++i) {
        int atom1, atom2;
        double distance;
        system.getConstraintParameters(i, atom1, atom2, distance);
        constraintIndices[i][0] = atom1;
        constraintIndices[i][1] = atom2;
        shakeParameters[i][0] = static_cast<RealOpenMM>(distance);
    }
}

void ReferenceIntegrateBrownianStepKernel::execute(OpenMMContextImpl& context, const BrownianIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = ((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData();
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    RealOpenMM** forceData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData()); // Reference code needs to be made const correct
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete shake;
        }
        dynamics = new ReferenceBrownianDynamics(
				context.getSystem().getNumAtoms(), 
				static_cast<RealOpenMM>(stepSize), 
				static_cast<RealOpenMM>(friction), 
				static_cast<RealOpenMM>(temperature) );
        shake = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, shakeParameters);
        dynamics->setReferenceShakeAlgorithm(shake);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumAtoms(), posData, velData, forceData, masses);
}

ReferenceApplyAndersenThermostatKernel::~ReferenceApplyAndersenThermostatKernel() {
    if (thermostat)
        delete thermostat;
    if (masses)
        delete[] masses;
}

void ReferenceApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    int numAtoms = system.getNumAtoms();
    masses = new RealOpenMM[numAtoms];
    for (size_t i = 0; i < numAtoms; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getAtomMass(i));
    this->thermostat = new ReferenceAndersenThermostat();
}

void ReferenceApplyAndersenThermostatKernel::execute(OpenMMContextImpl& context) {
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    thermostat->applyThermostat(
			context.getVelocities().getSize(), 
			velData, 
			masses, 
			static_cast<RealOpenMM>(context.getParameter(AndersenThermostat::Temperature)), 
			static_cast<RealOpenMM>(context.getParameter(AndersenThermostat::CollisionFrequency)), 
			static_cast<RealOpenMM>(context.getIntegrator().getStepSize()) );
}

void ReferenceCalcKineticEnergyKernel::initialize(const System& system) {
    int numAtoms = system.getNumAtoms();
    masses.resize(numAtoms);
    for (size_t i = 0; i < numAtoms; ++i)
        masses[i] = system.getAtomMass(i);
}

double ReferenceCalcKineticEnergyKernel::execute(OpenMMContextImpl& context) {
    RealOpenMM** velData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData()); // Reference code needs to be made const correct
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i)
        energy += masses[i]*(velData[i][0]*velData[i][0]+velData[i][1]*velData[i][1]+velData[i][2]*velData[i][2]);
    return 0.5*energy;
}

void ReferenceRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    frequency = force.getFrequency();
    masses.resize(system.getNumAtoms());
    for (size_t i = 0; i < masses.size(); ++i)
        masses[i] = system.getAtomMass(i);
}

void ReferenceRemoveCMMotionKernel::execute(OpenMMContextImpl& context) {
    int step = std::floor(context.getTime()/context.getIntegrator().getStepSize());
    if (step%frequency != 0)
        return;
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    
    // Calculate the center of mass momentum.
    
    RealOpenMM momentum[] = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < masses.size(); ++i) {
        momentum[0] += static_cast<RealOpenMM>( masses[i]*velData[i][0] );
        momentum[1] += static_cast<RealOpenMM>( masses[i]*velData[i][1] );
        momentum[2] += static_cast<RealOpenMM>( masses[i]*velData[i][2] );
    }
    
    // Adjust the atom velocities.
    
    momentum[0] /= static_cast<RealOpenMM>( masses.size() );
    momentum[1] /= static_cast<RealOpenMM>( masses.size() );
    momentum[2] /= static_cast<RealOpenMM>( masses.size() );
    for (size_t i = 0; i < masses.size(); ++i) {
        velData[i][0] -= static_cast<RealOpenMM>( momentum[0]/masses[i] );
        velData[i][1] -= static_cast<RealOpenMM>( momentum[1]/masses[i] );
        velData[i][2] -= static_cast<RealOpenMM>( momentum[2]/masses[i] );
    }
}
