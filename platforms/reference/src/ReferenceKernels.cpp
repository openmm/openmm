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
#include "gbsa/CpuObc.h"
#include "gbsa/CpuGBVI.h"
#include "SimTKReference/ReferenceAndersenThermostat.h"
#include "SimTKReference/ReferenceAngleBondIxn.h"
#include "SimTKReference/ReferenceBondForce.h"
#include "SimTKReference/ReferenceBrownianDynamics.h"
#include "SimTKReference/ReferenceCCMAAlgorithm.h"
#include "SimTKReference/ReferenceCustomAngleIxn.h"
#include "SimTKReference/ReferenceCustomBondIxn.h"
#include "SimTKReference/ReferenceCustomExternalIxn.h"
#include "SimTKReference/ReferenceCustomGBIxn.h"
#include "SimTKReference/ReferenceCustomNonbondedIxn.h"
#include "SimTKReference/ReferenceCustomTorsionIxn.h"
#include "SimTKReference/ReferenceHarmonicBondIxn.h"
#include "SimTKReference/ReferenceLJCoulomb14.h"
#include "SimTKReference/ReferenceLJCoulombIxn.h"
#include "SimTKReference/ReferenceProperDihedralBond.h"
#include "SimTKReference/ReferenceRbDihedralBond.h"
#include "SimTKReference/ReferenceStochasticDynamics.h"
#include "SimTKReference/ReferenceVariableStochasticDynamics.h"
#include "SimTKReference/ReferenceVariableVerletDynamics.h"
#include "SimTKReference/ReferenceVerletDynamics.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include "lepton/CustomFunction.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include <cmath>
#include <limits>

using namespace OpenMM;
using namespace std;

static int** allocateIntArray(int length, int width) {
    int** array = new int*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new int[width];
    return array;
}

static RealOpenMM** allocateRealArray(int length, int width) {
    RealOpenMM** array = new RealOpenMM*[length];
    for (int i = 0; i < length; ++i)
        array[i] = new RealOpenMM[width];
    return array;
}

static int** copyToArray(const vector<vector<int> > vec) {
    if (vec.size() == 0)
        return new int*[0];
    int** array = allocateIntArray(vec.size(), vec[0].size());
    for (size_t i = 0; i < vec.size(); ++i)
        for (size_t j = 0; j < vec[i].size(); ++j)
            array[i][j] = vec[i][j];
    return array;
}

static RealOpenMM** copyToArray(const vector<vector<double> > vec) {
    if (vec.size() == 0)
        return new RealOpenMM*[0];
    RealOpenMM** array = allocateRealArray(vec.size(), vec[0].size());
    for (size_t i = 0; i < vec.size(); ++i)
        for (size_t j = 0; j < vec[i].size(); ++j)
            array[i][j] = static_cast<RealOpenMM>(vec[i][j]);
    return array;
}

static void disposeIntArray(int** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static void disposeRealArray(RealOpenMM** array, int size) {
    if (array) {
        for (int i = 0; i < size; ++i)
            delete[] array[i];
        delete[] array;
    }
}

static RealOpenMM** extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->positions;
}

static RealOpenMM** extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->velocities;
}

static RealOpenMM** extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealOpenMM**) data->forces;
}

static void findAnglesForCCMA(const System& system, vector<ReferenceCCMAAlgorithm::AngleInfo>& angles) {
    for (int i = 0; i < system.getNumForces(); i++) {
        const HarmonicAngleForce* force = dynamic_cast<const HarmonicAngleForce*>(&system.getForce(i));
        if (force != NULL) {
            for (int j = 0; j < force->getNumAngles(); j++) {
                int atom1, atom2, atom3;
                double angle, k;
                force->getAngleParameters(j, atom1, atom2, atom3, angle, k);
                angles.push_back(ReferenceCCMAAlgorithm::AngleInfo(atom1, atom2, atom3, (RealOpenMM)angle));
            }
        }
    }
}

void ReferenceCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void ReferenceCalcForcesAndEnergyKernel::beginForceComputation(ContextImpl& context) {
    int numParticles = context.getSystem().getNumParticles();
    RealOpenMM** forceData = extractForces(context);
    for (int i = 0; i < numParticles; ++i) {
        forceData[i][0] = (RealOpenMM) 0.0;
        forceData[i][1] = (RealOpenMM) 0.0;
        forceData[i][2] = (RealOpenMM) 0.0;
    }
}

void ReferenceCalcForcesAndEnergyKernel::finishForceComputation(ContextImpl& context) {
}

void ReferenceCalcForcesAndEnergyKernel::beginEnergyComputation(ContextImpl& context) {
}

double ReferenceCalcForcesAndEnergyKernel::finishEnergyComputation(ContextImpl& context) {
    return 0.0;
}

void ReferenceUpdateStateDataKernel::initialize(const System& system) {
}

double ReferenceUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return data.time;
}

void ReferenceUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    data.time = time;
}

void ReferenceUpdateStateDataKernel::getPositions(ContextImpl& context, std::vector<Vec3>& positions) {
    int numParticles = context.getSystem().getNumParticles();
    RealOpenMM** posData = extractPositions(context);
    positions.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(posData[i][0], posData[i][1], posData[i][2]);
}

void ReferenceUpdateStateDataKernel::setPositions(ContextImpl& context, const std::vector<Vec3>& positions) {
    int numParticles = context.getSystem().getNumParticles();
    RealOpenMM** posData = extractPositions(context);
    for (int i = 0; i < numParticles; ++i) {
        posData[i][0] = (RealOpenMM) positions[i][0];
        posData[i][1] = (RealOpenMM) positions[i][1];
        posData[i][2] = (RealOpenMM) positions[i][2];
    }
}

void ReferenceUpdateStateDataKernel::getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    RealOpenMM** velData = extractVelocities(context);
    velocities.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        velocities[i] = Vec3(velData[i][0], velData[i][1], velData[i][2]);
}

void ReferenceUpdateStateDataKernel::setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    RealOpenMM** velData = extractVelocities(context);
    for (int i = 0; i < numParticles; ++i) {
        velData[i][0] = (RealOpenMM) velocities[i][0];
        velData[i][1] = (RealOpenMM) velocities[i][1];
        velData[i][2] = (RealOpenMM) velocities[i][2];
    }
}

void ReferenceUpdateStateDataKernel::getForces(ContextImpl& context, std::vector<Vec3>& forces) {
    int numParticles = context.getSystem().getNumParticles();
    RealOpenMM** forceData = extractForces(context);
    forces.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        forces[i] = Vec3(forceData[i][0], forceData[i][1], forceData[i][2]);
}

ReferenceCalcHarmonicBondForceKernel::~ReferenceCalcHarmonicBondForceKernel() {
    disposeIntArray(bondIndexArray, numBonds);
    disposeRealArray(bondParamArray, numBonds);
}

void ReferenceCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    numBonds = force.getNumBonds();
    bondIndexArray = allocateIntArray(numBonds, 2);
    bondParamArray = allocateRealArray(numBonds, 2);
    for (int i = 0; i < force.getNumBonds(); ++i) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        bondParamArray[i][0] = (RealOpenMM) length;
        bondParamArray[i][1] = (RealOpenMM) k;
    }
}

void ReferenceCalcHarmonicBondForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, 0, 0, 0, harmonicBond);
}

double ReferenceCalcHarmonicBondForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numBonds];
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    for (int i = 0; i < numBonds; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, energyArray, 0, &energy, harmonicBond);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcCustomBondForceKernel::~ReferenceCalcCustomBondForceKernel() {
    disposeIntArray(bondIndexArray, numBonds);
    disposeRealArray(bondParamArray, numBonds);
}

void ReferenceCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    numBonds = force.getNumBonds();
    int numParameters = force.getNumPerBondParameters();

    // Build the arrays.

    bondIndexArray = allocateIntArray(numBonds, 2);
    bondParamArray = allocateRealArray(numBonds, numParameters);
    vector<double> params;
    for (int i = 0; i < force.getNumBonds(); ++i) {
        int particle1, particle2;
        force.getBondParameters(i, particle1, particle2, params);
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = (RealOpenMM) params[j];
    }

    // Parse the expression used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    energyExpression = expression.createProgram();
    forceExpression = expression.differentiate("r").optimize().createProgram();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerBondParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
}

void ReferenceCalcCustomBondForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceBondForce refBondForce;
    ReferenceCustomBondIxn harmonicBond(energyExpression, forceExpression, parameterNames, globalParameters);
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, 0, 0, 0, harmonicBond);
}

double ReferenceCalcCustomBondForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numBonds];
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceBondForce refBondForce;
    ReferenceCustomBondIxn harmonicBond(energyExpression, forceExpression, parameterNames, globalParameters);
    for (int i = 0; i < numBonds; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, energyArray, 0, &energy, harmonicBond);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcHarmonicAngleForceKernel::~ReferenceCalcHarmonicAngleForceKernel() {
    disposeIntArray(angleIndexArray, numAngles);
    disposeRealArray(angleParamArray, numAngles);
}

void ReferenceCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    angleIndexArray = allocateIntArray(numAngles, 3);
    angleParamArray = allocateRealArray(numAngles, 2);
    for (int i = 0; i < force.getNumAngles(); ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        angleParamArray[i][0] = (RealOpenMM) angle;
        angleParamArray[i][1] = (RealOpenMM) k;
    }
}

void ReferenceCalcHarmonicAngleForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceBondForce refBondForce;
    ReferenceAngleBondIxn angleBond;
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, 0, 0, 0, angleBond);
}

double ReferenceCalcHarmonicAngleForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numAngles];
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceAngleBondIxn angleBond;
    for (int i = 0; i < numAngles; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, energyArray, 0, &energy, angleBond);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcCustomAngleForceKernel::~ReferenceCalcCustomAngleForceKernel() {
    disposeIntArray(angleIndexArray, numAngles);
    disposeRealArray(angleParamArray, numAngles);
}

void ReferenceCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    numAngles = force.getNumAngles();
    int numParameters = force.getNumPerAngleParameters();

    // Build the arrays.

    angleIndexArray = allocateIntArray(numAngles, 3);
    angleParamArray = allocateRealArray(numAngles, numParameters);
    vector<double> params;
    for (int i = 0; i < force.getNumAngles(); ++i) {
        int particle1, particle2, particle3;
        force.getAngleParameters(i, particle1, particle2, particle3, params);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        for (int j = 0; j < numParameters; j++)
            angleParamArray[i][j] = (RealOpenMM) params[j];
    }

    // Parse the expression used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    energyExpression = expression.createProgram();
    forceExpression = expression.differentiate("theta").optimize().createProgram();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerAngleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
}

void ReferenceCalcCustomAngleForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceBondForce refBondForce;
    ReferenceCustomAngleIxn customAngle(energyExpression, forceExpression, parameterNames, globalParameters);
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, 0, 0, 0, customAngle);
}

double ReferenceCalcCustomAngleForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numAngles];
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceBondForce refBondForce;
    ReferenceCustomAngleIxn customAngle(energyExpression, forceExpression, parameterNames, globalParameters);
    for (int i = 0; i < numAngles; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, energyArray, 0, &energy, customAngle);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcPeriodicTorsionForceKernel::~ReferenceCalcPeriodicTorsionForceKernel() {
    disposeIntArray(torsionIndexArray, numTorsions);
    disposeRealArray(torsionParamArray, numTorsions);
}

void ReferenceCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray = allocateIntArray(numTorsions, 4);
    torsionParamArray = allocateRealArray(numTorsions, 3);
    for (int i = 0; i < force.getNumTorsions(); ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = (RealOpenMM) k;
        torsionParamArray[i][1] = (RealOpenMM) phase;
        torsionParamArray[i][2] = (RealOpenMM) periodicity;
    }
}

void ReferenceCalcPeriodicTorsionForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceBondForce refBondForce;
    ReferenceProperDihedralBond periodicTorsionBond;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, 0, 0, 0, periodicTorsionBond);
}

double ReferenceCalcPeriodicTorsionForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numTorsions];
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceProperDihedralBond periodicTorsionBond;
    for (int i = 0; i < numTorsions; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, energyArray, 0, &energy, periodicTorsionBond);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcRBTorsionForceKernel::~ReferenceCalcRBTorsionForceKernel() {
    disposeIntArray(torsionIndexArray, numTorsions);
    disposeRealArray(torsionParamArray, numTorsions);
}

void ReferenceCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray = allocateIntArray(numTorsions, 4);
    torsionParamArray = allocateRealArray(numTorsions, 6);
    for (int i = 0; i < force.getNumTorsions(); ++i) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = (RealOpenMM) c0;
        torsionParamArray[i][1] = (RealOpenMM) c1;
        torsionParamArray[i][2] = (RealOpenMM) c2;
        torsionParamArray[i][3] = (RealOpenMM) c3;
        torsionParamArray[i][4] = (RealOpenMM) c4;
        torsionParamArray[i][5] = (RealOpenMM) c5;
    }
}

void ReferenceCalcRBTorsionForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceBondForce refBondForce;
    ReferenceRbDihedralBond rbTorsionBond;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, 0, 0, 0, rbTorsionBond);
}

double ReferenceCalcRBTorsionForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numTorsions];
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceRbDihedralBond rbTorsionBond;
    for (int i = 0; i < numTorsions; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, energyArray, 0, &energy, rbTorsionBond);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcCustomTorsionForceKernel::~ReferenceCalcCustomTorsionForceKernel() {
    disposeIntArray(torsionIndexArray, numTorsions);
    disposeRealArray(torsionParamArray, numTorsions);
}

void ReferenceCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    int numParameters = force.getNumPerTorsionParameters();

    // Build the arrays.

    torsionIndexArray = allocateIntArray(numTorsions, 4);
    torsionParamArray = allocateRealArray(numTorsions, numParameters);
    vector<double> params;
    for (int i = 0; i < force.getNumTorsions(); ++i) {
        int particle1, particle2, particle3, particle4;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, params);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        for (int j = 0; j < numParameters; j++)
            torsionParamArray[i][j] = (RealOpenMM) params[j];
    }

    // Parse the expression used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    energyExpression = expression.createProgram();
    forceExpression = expression.differentiate("theta").optimize().createProgram();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerTorsionParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
}

void ReferenceCalcCustomTorsionForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceBondForce refBondForce;
    ReferenceCustomTorsionIxn customTorsion(energyExpression, forceExpression, parameterNames, globalParameters);
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, 0, 0, 0, customTorsion);
}

double ReferenceCalcCustomTorsionForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM* energyArray = new RealOpenMM[numTorsions];
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceBondForce refBondForce;
    ReferenceCustomTorsionIxn customTorsion(energyExpression, forceExpression, parameterNames, globalParameters);
    for (int i = 0; i < numTorsions; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, energyArray, 0, &energy, customTorsion);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    delete[] energyArray;
    return energy;
}

ReferenceCalcNonbondedForceKernel::~ReferenceCalcNonbondedForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    disposeIntArray(exclusionArray, numParticles);
    disposeIntArray(bonded14IndexArray, num14);
    disposeRealArray(bonded14ParamArray, num14);
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {

    // Identify which exceptions are 1-4 interactions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
        if (chargeProd != 0.0 || epsilon != 0.0)
            nb14s.push_back(i);
    }

    // Build the arrays.

    num14 = nb14s.size();
    bonded14IndexArray = allocateIntArray(num14, 2);
    bonded14ParamArray = allocateRealArray(num14, 3);
    particleParamArray = allocateRealArray(numParticles, 3);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, depth;
        force.getParticleParameters(i, charge, radius, depth);
        particleParamArray[i][0] = static_cast<RealOpenMM>(0.5*radius);
        particleParamArray[i][1] = static_cast<RealOpenMM>(2.0*sqrt(depth));
        particleParamArray[i][2] = static_cast<RealOpenMM>(charge);
    }
    this->exclusions = exclusions;
    exclusionArray = new int*[numParticles];
    for (int i = 0; i < numParticles; ++i) {
        exclusionArray[i] = new int[exclusions[i].size()+1];
        exclusionArray[i][0] = exclusions[i].size();
        int index = 0;
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            exclusionArray[i][++index] = *iter;
    }
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        double charge, radius, depth;
        force.getExceptionParameters(nb14s[i], particle1, particle2, charge, radius, depth);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
        bonded14ParamArray[i][0] = static_cast<RealOpenMM>(radius);
        bonded14ParamArray[i][1] = static_cast<RealOpenMM>(4.0*depth);
        bonded14ParamArray[i][2] = static_cast<RealOpenMM>(charge);
    }
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
    periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
    periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();
    if (nonbondedMethod == Ewald) {
        double alpha;
        NonbondedForceImpl::calcEwaldParameters(system, force, alpha, kmax[0], kmax[1], kmax[2]);
        ewaldAlpha = (RealOpenMM) alpha;
    }
    else if (nonbondedMethod == PME) {
        double alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSize[0], gridSize[1], gridSize[2]);
        ewaldAlpha = (RealOpenMM) alpha;
    }
    rfDielectric = (RealOpenMM)force.getReactionFieldDielectric();
}

void ReferenceCalcNonbondedForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceLJCoulombIxn clj;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, (periodic || ewald || pme) ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, rfDielectric);
    }
    if (periodic||ewald||pme)
        clj.setPeriodic(periodicBoxSize);
    if (ewald)
        clj.setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        clj.setUsePME(ewaldAlpha, gridSize);
    clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, 0);
    ReferenceBondForce refBondForce;
    ReferenceLJCoulomb14 nonbonded14;
    if (nonbondedMethod == CutoffNonPeriodic || nonbondedMethod == CutoffPeriodic)
        nonbonded14.setUseCutoff(nonbondedCutoff, rfDielectric);
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, 0, 0, 0, nonbonded14);
}

double ReferenceCalcNonbondedForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(numParticles, 3);
    RealOpenMM energy = 0;
    ReferenceLJCoulombIxn clj;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, (periodic || ewald || pme) ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, rfDielectric);
    }
    if (periodic || ewald || pme)
        clj.setPeriodic(periodicBoxSize);
    if (ewald)
        clj.setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        clj.setUsePME(ewaldAlpha, gridSize);
    clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, &energy);
    ReferenceBondForce refBondForce;
    ReferenceLJCoulomb14 nonbonded14;
    if (nonbondedMethod == CutoffNonPeriodic || nonbondedMethod == CutoffPeriodic)
        nonbonded14.setUseCutoff(nonbondedCutoff, rfDielectric);
    RealOpenMM* energyArray = new RealOpenMM[num14];
    for (int i = 0; i < num14; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, energyArray, 0, &energy, nonbonded14);
    disposeRealArray(forceData, numParticles);
    delete[] energyArray;
    return energy;
}

class ReferenceTabulatedFunction : public Lepton::CustomFunction {
public:
    ReferenceTabulatedFunction(double min, double max, const vector<double>& values, bool interpolating) :
            min(min), max(max), values(values), interpolating(interpolating) {
    }
    int getNumArguments() const {
        return 1;
    }
    /**
     * Given the function argument, find the local spline coefficients.
     */
    void findCoefficients(double& x, double* coeff) const {
        int length = values.size();
        double scale = (length-1)/(max-min);
        int index = (int) std::floor((x-min)*scale);
        double points[4];
        points[0] = (index == 0 ? 2*values[0]-values[1] : values[index-1]);
        points[1] = values[index];
        points[2] = (index > length-2 ? values[length-1] : values[index+1]);
        points[3] = (index > length-3 ? 2*values[length-1]-values[length-2] : values[index+2]);
        if (interpolating) {
            coeff[0] = points[1];
            coeff[1] = 0.5*(-points[0]+points[2]);
            coeff[2] = 0.5*(2.0*points[0]-5.0*points[1]+4.0*points[2]-points[3]);
            coeff[3] = 0.5*(-points[0]+3.0*points[1]-3.0*points[2]+points[3]);
        }
        else {
            coeff[0] = (points[0]+4.0*points[1]+points[2])/6.0;
            coeff[1] = (-3.0*points[0]+3.0*points[2])/6.0;
            coeff[2] = (3.0*points[0]-6.0*points[1]+3.0*points[2])/6.0;
            coeff[3] = (-points[0]+3.0*points[1]-3.0*points[2]+points[3])/6.0;
        }
        x = (x-min)*scale-index;
    }
    double evaluate(const double* arguments) const {
        double x = arguments[0];
        if (x < min || x > max)
            return 0.0;
        double coeff[4];
        findCoefficients(x, coeff);
        return coeff[0]+x*(coeff[1]+x*(coeff[2]+x*coeff[3]));
    }
    double evaluateDerivative(const double* arguments, const int* derivOrder) const {
        double x = arguments[0];
        if (x < min || x > max)
            return 0.0;
        double coeff[4];
        findCoefficients(x, coeff);
        double scale = (values.size()-1)/(max-min);
        return scale*(coeff[1]+x*(2.0*coeff[2]+x*3.0*coeff[3])); // We assume a first derivative, because that's the only order ever used by CustomNonbondedForce.
    }
    CustomFunction* clone() const {
        return new ReferenceTabulatedFunction(min, max, values, interpolating);
    }
    double min, max;
    vector<double> values;
    bool interpolating;
};

ReferenceCalcCustomNonbondedForceKernel::~ReferenceCalcCustomNonbondedForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    disposeIntArray(exclusionArray, numParticles);
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceCalcCustomNonbondedForceKernel::initialize(const System& system, const CustomNonbondedForce& force) {

    // Record the exclusions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }

    // Build the arrays.

    int numParameters = force.getNumPerParticleParameters();
    particleParamArray = allocateRealArray(numParticles, numParameters);
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    exclusionArray = new int*[numParticles];
    for (int i = 0; i < numParticles; ++i) {
        exclusionArray[i] = new int[exclusions[i].size()+1];
        exclusionArray[i][0] = exclusions[i].size();
        int index = 0;
        for (set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            exclusionArray[i][++index] = *iter;
    }
    nonbondedMethod = CalcCustomNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
    periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
    periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        bool interpolating;
        force.getFunctionParameters(i, name, values, min, max, interpolating);
        functions[name] = new ReferenceTabulatedFunction(min, max, values, interpolating);
    }

    // Parse the various expressions used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    energyExpression = expression.createProgram();
    forceExpression = expression.differentiate("r").optimize().createProgram();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
}

void ReferenceCalcCustomNonbondedForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceCustomNonbondedIxn ixn(energyExpression, forceExpression, parameterNames);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, periodic ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    if (periodic)
        ixn.setPeriodic(periodicBoxSize);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ixn.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, globalParameters, forceData, 0, 0);
}

double ReferenceCalcCustomNonbondedForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(numParticles, 3);
    RealOpenMM energy = 0;
    ReferenceCustomNonbondedIxn ixn(energyExpression, forceExpression, parameterNames);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, periodic ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    if (periodic)
        ixn.setPeriodic(periodicBoxSize);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ixn.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, globalParameters, forceData, 0, &energy);
    disposeRealArray(forceData, numParticles);
    return energy;
}

ReferenceCalcGBSAOBCForceKernel::~ReferenceCalcGBSAOBCForceKernel() {
    if (obc) {
        // delete obc->getObcParameters();
        delete obc;
    }
}

void ReferenceCalcGBSAOBCForceKernel::initialize(const System& system, const GBSAOBCForce& force) {
    int numParticles = system.getNumParticles();
    charges.resize(numParticles);
    vector<RealOpenMM> atomicRadii(numParticles);
    vector<RealOpenMM> scaleFactors(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        charges[i] = static_cast<RealOpenMM>(charge);
        atomicRadii[i] = static_cast<RealOpenMM>(radius);
        scaleFactors[i] = static_cast<RealOpenMM>(scalingFactor);
    }
    ObcParameters* obcParameters  = new ObcParameters(numParticles, ObcParameters::ObcTypeII);
    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
    obcParameters->setSolventDielectric( static_cast<RealOpenMM>(force.getSolventDielectric()) );
    obcParameters->setSoluteDielectric( static_cast<RealOpenMM>(force.getSoluteDielectric()) );
    if (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff)
        obcParameters->setUseCutoff(static_cast<RealOpenMM>(force.getCutoffDistance()));
    if (force.getNonbondedMethod() == GBSAOBCForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        RealOpenMM periodicBoxSize[3];
        periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
        periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
        periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
        obcParameters->setPeriodic(periodicBoxSize);
    }
    obc = new CpuObc(obcParameters);
    obc->setIncludeAceApproximation(true);
}

void ReferenceCalcGBSAOBCForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
}

double ReferenceCalcGBSAOBCForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    return obc->getEnergy();
}

ReferenceCalcGBVIForceKernel::~ReferenceCalcGBVIForceKernel() {
    if (gbvi) {
        delete gbvi;
    }
}

void ReferenceCalcGBVIForceKernel::initialize(const System& system, const GBVIForce& force, const std::vector<double> & inputScaledRadii ) {
    int numParticles = system.getNumParticles();
    charges.resize(numParticles);
    vector<RealOpenMM> atomicRadii(numParticles);
    vector<RealOpenMM> scaledRadii(numParticles);
    vector<RealOpenMM> gammas(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, gamma;
        force.getParticleParameters(i, charge, radius, gamma);
        charges[i]       = static_cast<RealOpenMM>(charge);
        atomicRadii[i]   = static_cast<RealOpenMM>(radius);
        gammas[i]        = static_cast<RealOpenMM>(gamma);
        scaledRadii[i]   = static_cast<RealOpenMM>(inputScaledRadii[i]);
    }
    GBVIParameters * gBVIParameters = new GBVIParameters(numParticles);
    gBVIParameters->setAtomicRadii(atomicRadii);
    gBVIParameters->setGammaParameters(gammas);
    gBVIParameters->setScaledRadii(scaledRadii);
    gBVIParameters->setSolventDielectric(static_cast<RealOpenMM>(force.getSolventDielectric()));
    gBVIParameters->setSoluteDielectric(static_cast<RealOpenMM>(force.getSoluteDielectric()));
    if (force.getNonbondedMethod() != GBVIForce::NoCutoff)
        gBVIParameters->setUseCutoff(static_cast<RealOpenMM>(force.getCutoffDistance()));
    if (force.getNonbondedMethod() == GBVIForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        RealOpenMM periodicBoxSize[3];
        periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
        periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
        periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
        gBVIParameters->setPeriodic(periodicBoxSize);
    }
    gbvi = new CpuGBVI(gBVIParameters);
}

void ReferenceCalcGBVIForceKernel::executeForces(ContextImpl& context) {

    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    RealOpenMM* bornRadii  = new RealOpenMM[context.getSystem().getNumParticles()];
    gbvi->computeBornRadii(posData, bornRadii, NULL ); 
    gbvi->computeBornForces(bornRadii, posData, &charges[0], forceData);
    delete[] bornRadii;
}

double ReferenceCalcGBVIForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM* bornRadii = new RealOpenMM[context.getSystem().getNumParticles()];
    gbvi->computeBornRadii(posData, bornRadii, NULL ); 
    RealOpenMM energy     = gbvi->computeBornEnergy(bornRadii ,posData, &charges[0]);
    delete[] bornRadii;
    return static_cast<double>(energy);
}

ReferenceCalcCustomGBForceKernel::~ReferenceCalcCustomGBForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceCalcCustomGBForceKernel::initialize(const System& system, const CustomGBForce& force) {
    if (force.getNumComputedValues() > 0) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(0, name, expression, type);
        if (type == CustomGBForce::SingleParticle)
            throw OpenMMException("ReferencePlatform requires that the first computed value for a CustomGBForce be of type ParticlePair or ParticlePairNoExclusions.");
        for (int i = 1; i < force.getNumComputedValues(); i++) {
            force.getComputedValueParameters(i, name, expression, type);
            if (type != CustomGBForce::SingleParticle)
                throw OpenMMException("ReferencePlatform requires that a CustomGBForce only have one computed value of type ParticlePair or ParticlePairNoExclusions.");
        }
    }

    // Record the exclusions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }

    // Build the arrays.

    int numPerParticleParameters = force.getNumPerParticleParameters();
    particleParamArray = allocateRealArray(numParticles, numPerParticleParameters);
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numPerParticleParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    for (int i = 0; i < numPerParticleParameters; i++)
        particleParameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    nonbondedMethod = CalcCustomGBForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    Vec3 boxVectors[3];
    system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
    periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
    periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++) {
        string name;
        vector<double> values;
        double min, max;
        bool interpolating;
        force.getFunctionParameters(i, name, values, min, max, interpolating);
        functions[name] = new ReferenceTabulatedFunction(min, max, values, interpolating);
    }

    // Parse the expressions for computed values.

    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(i, name, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        valueExpressions.push_back(ex.createProgram());
        valueTypes.push_back(type);
        valueNames.push_back(name);
        if (type == CustomGBForce::SingleParticle)
            valueDerivExpressions.push_back(ex.differentiate(valueNames[i-1]).optimize().createProgram());
        else
            valueDerivExpressions.push_back(ex.differentiate("r").optimize().createProgram());
    }

    // Parse the expressions for energy terms.

    energyDerivExpressions.resize(force.getNumEnergyTerms());
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        energyExpressions.push_back(ex.createProgram());
        energyTypes.push_back(type);
        if (type != CustomGBForce::SingleParticle)
            energyDerivExpressions[i].push_back(ex.differentiate("r").optimize().createProgram());
        for (int j = 0; j < force.getNumComputedValues(); j++) {
            if (type == CustomGBForce::SingleParticle)
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]).optimize().createProgram());
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"1").optimize().createProgram());
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"2").optimize().createProgram());
            }
        }
    }

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
}

void ReferenceCalcCustomGBForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    ReferenceCustomGBIxn ixn(valueExpressions, valueDerivExpressions, valueNames, valueTypes, energyExpressions,
        energyDerivExpressions, energyTypes, particleParameterNames);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, periodic ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    if (periodic)
        ixn.setPeriodic(periodicBoxSize);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ixn.calculateIxn(numParticles, posData, particleParamArray, exclusions, globalParameters, forceData, 0);
}

double ReferenceCalcCustomGBForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(numParticles, 3);
    RealOpenMM energy = 0;
    ReferenceCustomGBIxn ixn(valueExpressions, valueDerivExpressions, valueNames, valueTypes, energyExpressions,
        energyDerivExpressions, energyTypes, particleParameterNames);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, periodic ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    if (periodic)
        ixn.setPeriodic(periodicBoxSize);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ixn.calculateIxn(numParticles, posData, particleParamArray, exclusions, globalParameters, forceData, &energy);
    disposeRealArray(forceData, numParticles);
    return energy;
}

ReferenceCalcCustomExternalForceKernel::~ReferenceCalcCustomExternalForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
}

void ReferenceCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    numParticles = force.getNumParticles();
    int numParameters = force.getNumPerParticleParameters();

    // Build the arrays.

    particles.resize(numParticles);
    particleParamArray = allocateRealArray(numParticles, numParameters);
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        force.getParticleParameters(i, particles[i], params);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = (RealOpenMM) params[j];
    }

    // Parse the expression used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    energyExpression = expression.createProgram();
    forceExpressionX = expression.differentiate("x").optimize().createProgram();
    forceExpressionY = expression.differentiate("y").optimize().createProgram();
    forceExpressionZ = expression.differentiate("z").optimize().createProgram();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
}

void ReferenceCalcCustomExternalForceKernel::executeForces(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = extractForces(context);
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceCustomExternalIxn force(energyExpression, forceExpressionX, forceExpressionY, forceExpressionZ, parameterNames, globalParameters);
    for (int i = 0; i < numParticles; ++i)
        force.calculateForce(particles[i], posData, particleParamArray[i], forceData, 0);
}

double ReferenceCalcCustomExternalForceKernel::executeEnergy(ContextImpl& context) {
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceCustomExternalIxn force(energyExpression, forceExpressionX, forceExpressionY, forceExpressionZ, parameterNames, globalParameters);
    for (int i = 0; i < numParticles; ++i)
        force.calculateForce(particles[i], posData, particleParamArray[i], forceData, &energy);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    return energy;
}

ReferenceIntegrateVerletStepKernel::~ReferenceIntegrateVerletStepKernel() {
    if (dynamics)
        delete dynamics;
    if (constraints)
        delete constraints;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (constraintDistances)
        delete[] constraintDistances;
}

void ReferenceIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    constraintDistances = new RealOpenMM[numConstraints];
    for (int i = 0; i < numConstraints; ++i) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraintIndices[i][0] = particle1;
        constraintIndices[i][1] = particle2;
        constraintDistances[i] = static_cast<RealOpenMM>(distance);
    }
}

void ReferenceIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** velData = extractVelocities(context);
    RealOpenMM** forceData = extractForces(context);
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        dynamics = new ReferenceVerletDynamics(context.getSystem().getNumParticles(), static_cast<RealOpenMM>(stepSize) );
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(context.getSystem(), angles);
        constraints = new ReferenceCCMAAlgorithm(context.getSystem().getNumParticles(), numConstraints, constraintIndices, constraintDistances, masses, angles, (RealOpenMM)integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses);
    data.time += stepSize;
    data.stepCount++;
}

ReferenceIntegrateLangevinStepKernel::~ReferenceIntegrateLangevinStepKernel() {
    if (dynamics)
        delete dynamics;
    if (constraints)
        delete constraints;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (constraintDistances)
        delete[] constraintDistances;
}

void ReferenceIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    constraintDistances = new RealOpenMM[numConstraints];
    for (int i = 0; i < numConstraints; ++i) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraintIndices[i][0] = particle1;
        constraintIndices[i][1] = particle2;
        constraintDistances[i] = static_cast<RealOpenMM>(distance);
    }
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** velData = extractVelocities(context);
    RealOpenMM** forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );
        dynamics = new ReferenceStochasticDynamics(
				context.getSystem().getNumParticles(), 
				static_cast<RealOpenMM>(stepSize), 
				static_cast<RealOpenMM>(tau), 
				static_cast<RealOpenMM>(temperature) );
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(context.getSystem(), angles);
        constraints = new ReferenceCCMAAlgorithm(context.getSystem().getNumParticles(), numConstraints, constraintIndices, constraintDistances, masses, angles, (RealOpenMM)integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses);
    data.time += stepSize;
    data.stepCount++;
}

ReferenceIntegrateBrownianStepKernel::~ReferenceIntegrateBrownianStepKernel() {
    if (dynamics)
        delete dynamics;
    if (constraints)
        delete constraints;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (constraintDistances)
        delete[] constraintDistances;
}

void ReferenceIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    constraintDistances = new RealOpenMM[numConstraints];
    for (int i = 0; i < numConstraints; ++i) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraintIndices[i][0] = particle1;
        constraintIndices[i][1] = particle2;
        constraintDistances[i] = static_cast<RealOpenMM>(distance);
    }
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** velData = extractVelocities(context);
    RealOpenMM** forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        dynamics = new ReferenceBrownianDynamics(
				context.getSystem().getNumParticles(), 
				static_cast<RealOpenMM>(stepSize), 
				static_cast<RealOpenMM>(friction), 
				static_cast<RealOpenMM>(temperature) );
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(context.getSystem(), angles);
        constraints = new ReferenceCCMAAlgorithm(context.getSystem().getNumParticles(), numConstraints, constraintIndices, constraintDistances, masses, angles, (RealOpenMM)integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses);
    data.time += stepSize;
    data.stepCount++;
}

ReferenceIntegrateVariableLangevinStepKernel::~ReferenceIntegrateVariableLangevinStepKernel() {
    if (dynamics)
        delete dynamics;
    if (constraints)
        delete constraints;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (constraintDistances)
        delete[] constraintDistances;
}

void ReferenceIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    constraintDistances = new RealOpenMM[numConstraints];
    for (int i = 0; i < numConstraints; ++i) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraintIndices[i][0] = particle1;
        constraintIndices[i][1] = particle2;
        constraintDistances[i] = static_cast<RealOpenMM>(distance);
    }
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double errorTol = integrator.getErrorTolerance();
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** velData = extractVelocities(context);
    RealOpenMM** forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || errorTol != prevErrorTol) {
        // Recreate the computation objects with the new parameters.

        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );
        dynamics = new ReferenceVariableStochasticDynamics(context.getSystem().getNumParticles(), (RealOpenMM) tau, (RealOpenMM) temperature, (RealOpenMM) errorTol);
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(context.getSystem(), angles);
        constraints = new ReferenceCCMAAlgorithm(context.getSystem().getNumParticles(), numConstraints, constraintIndices, constraintDistances, masses, angles, (RealOpenMM)integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevTemp = temperature;
        prevFriction = friction;
        prevErrorTol = errorTol;
    }
    RealOpenMM maxStepSize = (RealOpenMM) (maxTime-data.time);
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses, maxStepSize);
    data.time += dynamics->getDeltaT();
    if (dynamics->getDeltaT() == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
}

ReferenceIntegrateVariableVerletStepKernel::~ReferenceIntegrateVariableVerletStepKernel() {
    if (dynamics)
        delete dynamics;
    if (constraints)
        delete constraints;
    if (masses)
        delete[] masses;
    if (constraintIndices)
        disposeIntArray(constraintIndices, numConstraints);
    if (constraintDistances)
        delete[] constraintDistances;
}

void ReferenceIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    numConstraints = system.getNumConstraints();
    constraintIndices = allocateIntArray(numConstraints, 2);
    constraintDistances = new RealOpenMM[numConstraints];
    for (int i = 0; i < numConstraints; ++i) {
        int particle1, particle2;
        double distance;
        system.getConstraintParameters(i, particle1, particle2, distance);
        constraintIndices[i][0] = particle1;
        constraintIndices[i][1] = particle2;
        constraintDistances[i] = static_cast<RealOpenMM>(distance);
    }
}

void ReferenceIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    double errorTol = integrator.getErrorTolerance();
    RealOpenMM** posData = extractPositions(context);
    RealOpenMM** velData = extractVelocities(context);
    RealOpenMM** forceData = extractForces(context);
    if (dynamics == 0 || errorTol != prevErrorTol) {
        // Recreate the computation objects with the new parameters.

        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        dynamics = new ReferenceVariableVerletDynamics(context.getSystem().getNumParticles(), (RealOpenMM) errorTol);
        vector<ReferenceCCMAAlgorithm::AngleInfo> angles;
        findAnglesForCCMA(context.getSystem(), angles);
        constraints = new ReferenceCCMAAlgorithm(context.getSystem().getNumParticles(), numConstraints, constraintIndices, constraintDistances, masses, angles, (RealOpenMM)integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevErrorTol = errorTol;
    }
    RealOpenMM maxStepSize = (RealOpenMM) (maxTime-data.time);
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses, maxStepSize);
    data.time += dynamics->getDeltaT();
    if (dynamics->getDeltaT() == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
}

ReferenceApplyAndersenThermostatKernel::~ReferenceApplyAndersenThermostatKernel() {
    if (thermostat)
        delete thermostat;
    if (masses)
        delete[] masses;
}

void ReferenceApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    int numParticles = system.getNumParticles();
    masses = new RealOpenMM[numParticles];
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    this->thermostat = new ReferenceAndersenThermostat();
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) thermostat.getRandomNumberSeed());
}

void ReferenceApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    RealOpenMM** velData = extractVelocities(context);
    thermostat->applyThermostat(
			context.getSystem().getNumParticles(),
			velData, 
			masses, 
			static_cast<RealOpenMM>(context.getParameter(AndersenThermostat::Temperature())), 
			static_cast<RealOpenMM>(context.getParameter(AndersenThermostat::CollisionFrequency())), 
			static_cast<RealOpenMM>(context.getIntegrator().getStepSize()) );
}

void ReferenceCalcKineticEnergyKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

double ReferenceCalcKineticEnergyKernel::execute(ContextImpl& context) {
    RealOpenMM** velData = extractVelocities(context);
    double energy = 0.0;
    for (size_t i = 0; i < masses.size(); ++i)
        energy += masses[i]*(velData[i][0]*velData[i][0]+velData[i][1]*velData[i][1]+velData[i][2]*velData[i][2]);
    return 0.5*energy;
}

void ReferenceRemoveCMMotionKernel::initialize(const System& system, const CMMotionRemover& force) {
    frequency = force.getFrequency();
    masses.resize(system.getNumParticles());
    for (size_t i = 0; i < masses.size(); ++i)
        masses[i] = system.getParticleMass(i);
}

void ReferenceRemoveCMMotionKernel::execute(ContextImpl& context) {
    if (data.stepCount%frequency != 0)
        return;
    RealOpenMM** velData = extractVelocities(context);
    
    // Calculate the center of mass momentum.
    
    RealOpenMM momentum[] = {0.0, 0.0, 0.0};
    RealOpenMM mass = 0.0;
    for (size_t i = 0; i < masses.size(); ++i) {
        momentum[0] += static_cast<RealOpenMM>( masses[i]*velData[i][0] );
        momentum[1] += static_cast<RealOpenMM>( masses[i]*velData[i][1] );
        momentum[2] += static_cast<RealOpenMM>( masses[i]*velData[i][2] );
        mass += static_cast<RealOpenMM>( masses[i] );
    }
    
    // Adjust the particle velocities.
    
    momentum[0] /= mass;
    momentum[1] /= mass;
    momentum[2] /= mass;
    for (size_t i = 0; i < masses.size(); ++i) {
        velData[i][0] -= momentum[0];
        velData[i][1] -= momentum[1];
        velData[i][2] -= momentum[2];
    }
}
