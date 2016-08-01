/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "ReferenceObc.h"
#include "ReferenceAndersenThermostat.h"
#include "ReferenceAngleBondIxn.h"
#include "ReferenceBondForce.h"
#include "ReferenceBrownianDynamics.h"
#include "ReferenceCCMAAlgorithm.h"
#include "ReferenceCMAPTorsionIxn.h"
#include "ReferenceConstraints.h"
#include "ReferenceCustomAngleIxn.h"
#include "ReferenceCustomBondIxn.h"
#include "ReferenceCustomCentroidBondIxn.h"
#include "ReferenceCustomCompoundBondIxn.h"
#include "ReferenceCustomDynamics.h"
#include "ReferenceCustomExternalIxn.h"
#include "ReferenceCustomGBIxn.h"
#include "ReferenceCustomHbondIxn.h"
#include "ReferenceCustomNonbondedIxn.h"
#include "ReferenceCustomManyParticleIxn.h"
#include "ReferenceCustomTorsionIxn.h"
#include "ReferenceGayBerneForce.h"
#include "ReferenceHarmonicBondIxn.h"
#include "ReferenceLJCoulomb14.h"
#include "ReferenceLJCoulombIxn.h"
#include "ReferenceMonteCarloBarostat.h"
#include "ReferenceProperDihedralBond.h"
#include "ReferenceRbDihedralBond.h"
#include "ReferenceStochasticDynamics.h"
#include "ReferenceTabulatedFunction.h"
#include "ReferenceVariableStochasticDynamics.h"
#include "ReferenceVariableVerletDynamics.h"
#include "ReferenceVerletDynamics.h"
#include "ReferenceVirtualSites.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/internal/AndersenThermostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomCentroidBondForceImpl.h"
#include "openmm/internal/CustomCompoundBondForceImpl.h"
#include "openmm/internal/CustomHbondForceImpl.h"
#include "openmm/internal/CustomNonbondedForceImpl.h"
#include "openmm/internal/CMAPTorsionForceImpl.h"
#include "openmm/internal/NonbondedForceImpl.h"
#include "openmm/Integrator.h"
#include "openmm/OpenMMException.h"
#include "SimTKOpenMMUtilities.h"
#include "lepton/CustomFunction.h"
#include "lepton/Operation.h"
#include "lepton/Parser.h"
#include "lepton/ParsedExpression.h"
#include <cmath>
#include <iostream>
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

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static RealVec& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(RealVec*) data->periodicBoxSize;
}

static RealVec* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return (RealVec*) data->periodicBoxVectors;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
}

static map<string, double>& extractEnergyParameterDerivatives(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((map<string, double>*) data->energyParameterDerivatives);
}

/**
 * Make sure an expression doesn't use any undefined variables.
 */
static void validateVariables(const Lepton::ExpressionTreeNode& node, const set<string>& variables) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("Unknown variable in expression: "+op.getName());
    for (int i = 0; i < (int) node.getChildren().size(); i++)
        validateVariables(node.getChildren()[i], variables);
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.
    
    vector<RealVec> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (masses[i] > 0)
            shiftedVel[i] = velData[i]+forceData[i]*(timeShift/masses[i]);
        else
            shiftedVel[i] = velData[i];
    }
    
    // Apply constraints to them.
    
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; i++)
        inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
    extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}

void ReferenceCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void ReferenceCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    vector<RealVec>& forceData = extractForces(context);
    if (includeForces) {
        int numParticles = context.getSystem().getNumParticles();
        for (int i = 0; i < numParticles; ++i) {
            forceData[i][0] = (RealOpenMM) 0.0;
            forceData[i][1] = (RealOpenMM) 0.0;
            forceData[i][2] = (RealOpenMM) 0.0;
        }
    }
    else
        savedForces = forceData;
    for (map<string, double>::const_iterator iter = context.getParameters().begin(); iter != context.getParameters().end(); ++iter)
        extractEnergyParameterDerivatives(context)[iter->first] = 0;
}

double ReferenceCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    if (!includeForces)
        extractForces(context) = savedForces; // Restore the forces so computing the energy doesn't overwrite the forces with incorrect values.
    else
        ReferenceVirtualSites::distributeForces(context.getSystem(), extractPositions(context), extractForces(context));
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
    vector<RealVec>& posData = extractPositions(context);
    positions.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(posData[i][0], posData[i][1], posData[i][2]);
}

void ReferenceUpdateStateDataKernel::setPositions(ContextImpl& context, const std::vector<Vec3>& positions) {
    int numParticles = context.getSystem().getNumParticles();
    vector<RealVec>& posData = extractPositions(context);
    for (int i = 0; i < numParticles; ++i) {
        posData[i][0] = (RealOpenMM) positions[i][0];
        posData[i][1] = (RealOpenMM) positions[i][1];
        posData[i][2] = (RealOpenMM) positions[i][2];
    }
}

void ReferenceUpdateStateDataKernel::getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    vector<RealVec>& velData = extractVelocities(context);
    velocities.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        velocities[i] = Vec3(velData[i][0], velData[i][1], velData[i][2]);
}

void ReferenceUpdateStateDataKernel::setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    vector<RealVec>& velData = extractVelocities(context);
    for (int i = 0; i < numParticles; ++i) {
        velData[i][0] = (RealOpenMM) velocities[i][0];
        velData[i][1] = (RealOpenMM) velocities[i][1];
        velData[i][2] = (RealOpenMM) velocities[i][2];
    }
}

void ReferenceUpdateStateDataKernel::getForces(ContextImpl& context, std::vector<Vec3>& forces) {
    int numParticles = context.getSystem().getNumParticles();
    vector<RealVec>& forceData = extractForces(context);
    forces.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        forces[i] = Vec3(forceData[i][0], forceData[i][1], forceData[i][2]);
}

void ReferenceUpdateStateDataKernel::getEnergyParameterDerivatives(ContextImpl& context, map<string, double>& derivs) {
    derivs = extractEnergyParameterDerivatives(context);
}

void ReferenceUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    RealVec* vectors = extractBoxVectors(context);
    a = vectors[0];
    b = vectors[1];
    c = vectors[2];
}

void ReferenceUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) {
    RealVec& box = extractBoxSize(context);
    box[0] = (RealOpenMM) a[0];
    box[1] = (RealOpenMM) b[1];
    box[2] = (RealOpenMM) c[2];
    RealVec* vectors = extractBoxVectors(context);
    vectors[0] = a;
    vectors[1] = b;
    vectors[2] = c;
}

void ReferenceUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    int version = 2;
    stream.write((char*) &version, sizeof(int));
    stream.write((char*) &data.time, sizeof(data.time));
    vector<RealVec>& posData = extractPositions(context);
    stream.write((char*) &posData[0], sizeof(RealVec)*posData.size());
    vector<RealVec>& velData = extractVelocities(context);
    stream.write((char*) &velData[0], sizeof(RealVec)*velData.size());
    RealVec* vectors = extractBoxVectors(context);
    stream.write((char*) vectors, 3*sizeof(RealVec));
    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void ReferenceUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    int version;
    stream.read((char*) &version, sizeof(int));
    if (version != 2)
        throw OpenMMException("Checkpoint was created with a different version of OpenMM");
    stream.read((char*) &data.time, sizeof(data.time));
    vector<RealVec>& posData = extractPositions(context);
    stream.read((char*) &posData[0], sizeof(RealVec)*posData.size());
    vector<RealVec>& velData = extractVelocities(context);
    stream.read((char*) &velData[0], sizeof(RealVec)*velData.size());
    RealVec* vectors = extractBoxVectors(context);
    stream.read((char*) vectors, 3*sizeof(RealVec));
    SimTKOpenMMUtilities::loadCheckpoint(stream);
}

void ReferenceApplyConstraintsKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    inverseMasses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
        inverseMasses[i] = 1.0/masses[i];
    }
}

ReferenceApplyConstraintsKernel::~ReferenceApplyConstraintsKernel() {
}

void ReferenceApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
    vector<RealVec>& positions = extractPositions(context);
    extractConstraints(context).apply(positions, positions, inverseMasses, tol);
    ReferenceVirtualSites::computePositions(context.getSystem(), positions);
}

void ReferenceApplyConstraintsKernel::applyToVelocities(ContextImpl& context, double tol) {
    vector<RealVec>& positions = extractPositions(context);
    vector<RealVec>& velocities = extractVelocities(context);
    extractConstraints(context).applyToVelocities(positions, velocities, inverseMasses, tol);
}

void ReferenceVirtualSitesKernel::initialize(const System& system) {
}

void ReferenceVirtualSitesKernel::computePositions(ContextImpl& context) {
    vector<RealVec>& positions = extractPositions(context);
    ReferenceVirtualSites::computePositions(context.getSystem(), positions);
}

ReferenceCalcHarmonicBondForceKernel::~ReferenceCalcHarmonicBondForceKernel() {
    disposeIntArray(bondIndexArray, numBonds);
    disposeRealArray(bondParamArray, numBonds);
}

void ReferenceCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    numBonds = force.getNumBonds();
    bondIndexArray = allocateIntArray(numBonds, 2);
    bondParamArray = allocateRealArray(numBonds, 2);
    for (int i = 0; i < numBonds; ++i) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        bondParamArray[i][0] = (RealOpenMM) length;
        bondParamArray[i][1] = (RealOpenMM) k;
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    if (usePeriodic)
        harmonicBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, includeEnergy ? &energy : NULL, harmonicBond);
    return energy;
}

void ReferenceCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    for (int i = 0; i < numBonds; ++i) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        if (particle1 != bondIndexArray[i][0] || particle2 != bondIndexArray[i][1])
            throw OpenMMException("updateParametersInContext: The set of particles in a bond has changed");
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        bondParamArray[i][0] = (RealOpenMM) length;
        bondParamArray[i][1] = (RealOpenMM) k;
    }
}

ReferenceCalcCustomBondForceKernel::~ReferenceCalcCustomBondForceKernel() {
    disposeIntArray(bondIndexArray, numBonds);
    disposeRealArray(bondParamArray, numBonds);
}

void ReferenceCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    numBonds = force.getNumBonds();
    int numParameters = force.getNumPerBondParameters();
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    bondIndexArray = allocateIntArray(numBonds, 2);
    bondParamArray = allocateRealArray(numBonds, numParameters);
    vector<double> params;
    for (int i = 0; i < numBonds; ++i) {
        int particle1, particle2;
        force.getBondParameters(i, particle1, particle2, params);
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = (RealOpenMM) params[j];
    }

    // Parse the expression used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction()).optimize();
    energyExpression = expression.createCompiledExpression();
    forceExpression = expression.differentiate("r").createCompiledExpression();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerBondParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(expression.differentiate(param).createCompiledExpression());
    }
    set<string> variables;
    variables.insert("r");
    variables.insert(parameterNames.begin(), parameterNames.end());
    variables.insert(globalParameterNames.begin(), globalParameterNames.end());
    validateVariables(expression.getRootNode(), variables);
}

double ReferenceCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceCustomBondIxn bond(energyExpression, forceExpression, parameterNames, globalParameters, energyParamDerivExpressions);
    if (usePeriodic)
        bond.setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    for (int i = 0; i < numBonds; i++)
        bond.calculateBondIxn(bondIndexArray[i], posData, bondParamArray[i], forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    int numParameters = force.getNumPerBondParameters();
    vector<double> params;
    for (int i = 0; i < numBonds; ++i) {
        int particle1, particle2;
        force.getBondParameters(i, particle1, particle2, params);
        if (particle1 != bondIndexArray[i][0] || particle2 != bondIndexArray[i][1])
            throw OpenMMException("updateParametersInContext: The set of particles in a bond has changed");
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = (RealOpenMM) params[j];
    }
}

ReferenceCalcHarmonicAngleForceKernel::~ReferenceCalcHarmonicAngleForceKernel() {
    disposeIntArray(angleIndexArray, numAngles);
    disposeRealArray(angleParamArray, numAngles);
}

void ReferenceCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    angleIndexArray = allocateIntArray(numAngles, 3);
    angleParamArray = allocateRealArray(numAngles, 2);
    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        angleParamArray[i][0] = (RealOpenMM) angle;
        angleParamArray[i][1] = (RealOpenMM) k;
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceAngleBondIxn angleBond;
    if (usePeriodic)
        angleBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, includeEnergy ? &energy : NULL, angleBond);
    return energy;
}

void ReferenceCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        if (particle1 != angleIndexArray[i][0] || particle2 != angleIndexArray[i][1] || particle3 != angleIndexArray[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        angleParamArray[i][0] = (RealOpenMM) angle;
        angleParamArray[i][1] = (RealOpenMM) k;
    }
}

ReferenceCalcCustomAngleForceKernel::~ReferenceCalcCustomAngleForceKernel() {
    disposeIntArray(angleIndexArray, numAngles);
    disposeRealArray(angleParamArray, numAngles);
}

void ReferenceCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    numAngles = force.getNumAngles();
    int numParameters = force.getNumPerAngleParameters();
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    angleIndexArray = allocateIntArray(numAngles, 3);
    angleParamArray = allocateRealArray(numAngles, numParameters);
    vector<double> params;
    for (int i = 0; i < numAngles; ++i) {
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
    energyExpression = expression.createCompiledExpression();
    forceExpression = expression.differentiate("theta").createCompiledExpression();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerAngleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(expression.differentiate(param).createCompiledExpression());
    }
    set<string> variables;
    variables.insert("theta");
    variables.insert(parameterNames.begin(), parameterNames.end());
    variables.insert(globalParameterNames.begin(), globalParameterNames.end());
    validateVariables(expression.getRootNode(), variables);
}

double ReferenceCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceCustomAngleIxn customAngle(energyExpression, forceExpression, parameterNames, globalParameters, energyParamDerivExpressions);
    if (usePeriodic)
        customAngle.setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    for (int i = 0; i < numAngles; i++)
        customAngle.calculateBondIxn(angleIndexArray[i], posData, angleParamArray[i], forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    int numParameters = force.getNumPerAngleParameters();
    vector<double> params;
    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        force.getAngleParameters(i, particle1, particle2, particle3, params);
        if (particle1 != angleIndexArray[i][0] || particle2 != angleIndexArray[i][1] || particle3 != angleIndexArray[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        for (int j = 0; j < numParameters; j++)
            angleParamArray[i][j] = (RealOpenMM) params[j];
    }
}

ReferenceCalcPeriodicTorsionForceKernel::~ReferenceCalcPeriodicTorsionForceKernel() {
    disposeIntArray(torsionIndexArray, numTorsions);
    disposeRealArray(torsionParamArray, numTorsions);
}

void ReferenceCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray = allocateIntArray(numTorsions, 4);
    torsionParamArray = allocateRealArray(numTorsions, 3);
    for (int i = 0; i < numTorsions; ++i) {
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
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceProperDihedralBond periodicTorsionBond;
    if (usePeriodic)
        periodicTorsionBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, includeEnergy ? &energy : NULL, periodicTorsionBond);
    return energy;
}

void ReferenceCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        torsionParamArray[i][0] = (RealOpenMM) k;
        torsionParamArray[i][1] = (RealOpenMM) phase;
        torsionParamArray[i][2] = (RealOpenMM) periodicity;
    }
}

ReferenceCalcRBTorsionForceKernel::~ReferenceCalcRBTorsionForceKernel() {
    disposeIntArray(torsionIndexArray, numTorsions);
    disposeRealArray(torsionParamArray, numTorsions);
}

void ReferenceCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray = allocateIntArray(numTorsions, 4);
    torsionParamArray = allocateRealArray(numTorsions, 6);
    for (int i = 0; i < numTorsions; ++i) {
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
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceRbDihedralBond rbTorsionBond;
    if (usePeriodic)
        rbTorsionBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, includeEnergy ? &energy : NULL, rbTorsionBond);
    return energy;
}

void ReferenceCalcRBTorsionForceKernel::copyParametersToContext(ContextImpl& context, const RBTorsionForce& force) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        torsionParamArray[i][0] = (RealOpenMM) c0;
        torsionParamArray[i][1] = (RealOpenMM) c1;
        torsionParamArray[i][2] = (RealOpenMM) c2;
        torsionParamArray[i][3] = (RealOpenMM) c3;
        torsionParamArray[i][4] = (RealOpenMM) c4;
        torsionParamArray[i][5] = (RealOpenMM) c5;
    }
}

void ReferenceCalcCMAPTorsionForceKernel::initialize(const System& system, const CMAPTorsionForce& force) {
    int numMaps = force.getNumMaps();
    int numTorsions = force.getNumTorsions();
    coeff.resize(numMaps);
    vector<double> energy;
    vector<vector<double> > c;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        coeff[i].resize(size*size);
        for (int j = 0; j < size*size; j++) {
            coeff[i][j].resize(16);
            for (int k = 0; k < 16; k++)
                coeff[i][j][k] = c[j][k];
        }
    }
    torsionMaps.resize(numTorsions);
    torsionIndices.resize(numTorsions);
    for (int i = 0; i < numTorsions; i++) {
        torsionIndices[i].resize(8);
        force.getTorsionParameters(i, torsionMaps[i], torsionIndices[i][0], torsionIndices[i][1], torsionIndices[i][2],
            torsionIndices[i][3], torsionIndices[i][4], torsionIndices[i][5], torsionIndices[i][6], torsionIndices[i][7]);
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcCMAPTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM totalEnergy = 0;
    ReferenceCMAPTorsionIxn torsion(coeff, torsionMaps, torsionIndices);
    if (usePeriodic)
        torsion.setPeriodic(extractBoxVectors(context));
    torsion.calculateIxn(posData, forceData, &totalEnergy);
    return totalEnergy;
}

void ReferenceCalcCMAPTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CMAPTorsionForce& force) {
    int numMaps = force.getNumMaps();
    int numTorsions = force.getNumTorsions();
    if (coeff.size() != numMaps)
        throw OpenMMException("updateParametersInContext: The number of maps has changed");
    if (torsionMaps.size() != numTorsions)
        throw OpenMMException("updateParametersInContext: The number of CMAP torsions has changed");

    // Update the maps.

    vector<double> energy;
    vector<vector<double> > c;
    for (int i = 0; i < numMaps; i++) {
        int size;
        force.getMapParameters(i, size, energy);
        if (coeff[i].size() != size*size)
            throw OpenMMException("updateParametersInContext: The size of a map has changed");
        CMAPTorsionForceImpl::calcMapDerivatives(size, energy, c);
        for (int j = 0; j < size*size; j++)
            for (int k = 0; k < 16; k++)
                coeff[i][j][k] = c[j][k];
    }

    // Update the indices.

    for (int i = 0; i < numTorsions; i++) {
        int index[8];
        force.getTorsionParameters(i, torsionMaps[i], index[0], index[1], index[2], index[3], index[4], index[5], index[6], index[7]);
        for (int j = 0; j < 8; j++)
            if (index[j] != torsionIndices[i][j])
                throw OpenMMException("updateParametersInContext: The set of particles in a CMAP torsion has changed");
    }
}

ReferenceCalcCustomTorsionForceKernel::~ReferenceCalcCustomTorsionForceKernel() {
    disposeIntArray(torsionIndexArray, numTorsions);
    disposeRealArray(torsionParamArray, numTorsions);
}

void ReferenceCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    int numParameters = force.getNumPerTorsionParameters();
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    torsionIndexArray = allocateIntArray(numTorsions, 4);
    torsionParamArray = allocateRealArray(numTorsions, numParameters);
    vector<double> params;
    for (int i = 0; i < numTorsions; ++i) {
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
    energyExpression = expression.createCompiledExpression();
    forceExpression = expression.differentiate("theta").createCompiledExpression();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerTorsionParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(expression.differentiate(param).createCompiledExpression());
    }
    set<string> variables;
    variables.insert("theta");
    variables.insert(parameterNames.begin(), parameterNames.end());
    variables.insert(globalParameterNames.begin(), globalParameterNames.end());
    validateVariables(expression.getRootNode(), variables);
}

double ReferenceCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceCustomTorsionIxn customTorsion(energyExpression, forceExpression, parameterNames, globalParameters, energyParamDerivExpressions);
    if (usePeriodic)
        customTorsion.setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    for (int i = 0; i < numTorsions; i++)
        customTorsion.calculateBondIxn(torsionIndexArray[i], posData, torsionParamArray[i], forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    int numParameters = force.getNumPerTorsionParameters();
    vector<double> params;
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, params);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        for (int j = 0; j < numParameters; j++)
            torsionParamArray[i][j] = (RealOpenMM) params[j];
    }
}

ReferenceCalcNonbondedForceKernel::~ReferenceCalcNonbondedForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
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
    if (nonbondedMethod == NoCutoff) {
        neighborList = NULL;
        useSwitchingFunction = false;
    }
    else {
        neighborList = new NeighborList();
        useSwitchingFunction = force.getUseSwitchingFunction();
        switchingDistance = force.getSwitchingDistance();
    }
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
    if (force.getUseDispersionCorrection())
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
}

double ReferenceCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceLJCoulombIxn clj;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, extractBoxVectors(context), periodic || ewald || pme, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, rfDielectric);
    }
    if (periodic || ewald || pme) {
        RealVec* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 1.999999*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        clj.setPeriodic(boxVectors);
    }
    if (ewald)
        clj.setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        clj.setUsePME(ewaldAlpha, gridSize);
    if (useSwitchingFunction)
        clj.setUseSwitchingFunction(switchingDistance);
    clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusions, 0, forceData, 0, includeEnergy ? &energy : NULL, includeDirect, includeReciprocal);
    if (includeDirect) {
        ReferenceBondForce refBondForce;
        ReferenceLJCoulomb14 nonbonded14;
        refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, includeEnergy ? &energy : NULL, nonbonded14);
        if (periodic || ewald || pme) {
            RealVec* boxVectors = extractBoxVectors(context);
            energy += dispersionCoefficient/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]);
        }
    }
    return energy;
}

void ReferenceCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force) {
    if (force.getNumParticles() != numParticles)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");
    vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (chargeProd != 0.0 || epsilon != 0.0)
            nb14s.push_back(i);
    }
    if (nb14s.size() != num14)
        throw OpenMMException("updateParametersInContext: The number of non-excluded exceptions has changed");

    // Record the values.

    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, depth;
        force.getParticleParameters(i, charge, radius, depth);
        particleParamArray[i][0] = static_cast<RealOpenMM>(0.5*radius);
        particleParamArray[i][1] = static_cast<RealOpenMM>(2.0*sqrt(depth));
        particleParamArray[i][2] = static_cast<RealOpenMM>(charge);
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
    
    // Recompute the coefficient for the dispersion correction.

    NonbondedForce::NonbondedMethod method = force.getNonbondedMethod();
    if (force.getUseDispersionCorrection() && (method == NonbondedForce::CutoffPeriodic || method == NonbondedForce::Ewald || method == NonbondedForce::PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
}

void ReferenceCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME");
    alpha = ewaldAlpha;
    nx = gridSize[0];
    ny = gridSize[1];
    nz = gridSize[2];
}

ReferenceCalcCustomNonbondedForceKernel::~ReferenceCalcCustomNonbondedForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    if (neighborList != NULL)
        delete neighborList;
    if (forceCopy != NULL)
        delete forceCopy;
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
    nonbondedMethod = CalcCustomNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff) {
        neighborList = NULL;
        useSwitchingFunction = false;
    }
    else {
        neighborList = new NeighborList();
        useSwitchingFunction = force.getUseSwitchingFunction();
        switchingDistance = force.getSwitchingDistance();
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the various expressions used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    energyExpression = expression.createCompiledExpression();
    forceExpression = expression.differentiate("r").createCompiledExpression();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParameterNames.push_back(force.getGlobalParameterName(i));
        globalParamValues[force.getGlobalParameterName(i)] = force.getGlobalParameterDefaultValue(i);
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(expression.differentiate(param).createCompiledExpression());
    }
    set<string> variables;
    variables.insert("r");
    for (int i = 0; i < numParameters; i++) {
        variables.insert(parameterNames[i]+"1");
        variables.insert(parameterNames[i]+"2");
    }
    variables.insert(globalParameterNames.begin(), globalParameterNames.end());
    validateVariables(expression.getRootNode(), variables);

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
    
    // Record information for the long range correction.
    
    if (force.getNonbondedMethod() == CustomNonbondedForce::CutoffPeriodic && force.getUseLongRangeCorrection()) {
        forceCopy = new CustomNonbondedForce(force);
        hasInitializedLongRangeCorrection = false;
    }
    else {
        longRangeCoefficient = 0.0;
        hasInitializedLongRangeCorrection = true;
    }
    
    // Record the interaction groups.
    
    for (int i = 0; i < force.getNumInteractionGroups(); i++) {
        set<int> set1, set2;
        force.getInteractionGroupParameters(i, set1, set2);
        interactionGroups.push_back(make_pair(set1, set2));
    }
}

double ReferenceCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealVec* boxVectors = extractBoxVectors(context);
    RealOpenMM energy = 0;
    ReferenceCustomNonbondedIxn ixn(energyExpression, forceExpression, parameterNames, energyParamDerivExpressions);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, extractBoxVectors(context), periodic, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    if (periodic) {
        double minAllowedSize = 2*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        ixn.setPeriodic(boxVectors);
    }
    if (interactionGroups.size() > 0)
        ixn.setInteractionGroups(interactionGroups);
    bool globalParamsChanged = false;
    for (int i = 0; i < (int) globalParameterNames.size(); i++) {
        double value = context.getParameter(globalParameterNames[i]);
        if (globalParamValues[globalParameterNames[i]] != value)
            globalParamsChanged = true;
        globalParamValues[globalParameterNames[i]] = value;
    }
    if (useSwitchingFunction)
        ixn.setUseSwitchingFunction(switchingDistance);
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    ixn.calculatePairIxn(numParticles, posData, particleParamArray, exclusions, 0, globalParamValues, forceData, 0, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    
    // Add in the long range correction.
    
    if (!hasInitializedLongRangeCorrection || (globalParamsChanged && forceCopy != NULL)) {
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
    }
    double volume = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
    energy += longRangeCoefficient/volume;
    for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += longRangeCoefficientDerivs[i]/volume;
    return energy;
}

void ReferenceCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    
    // If necessary, recompute the long range correction.
    
    if (forceCopy != NULL) {
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs);
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }
}

ReferenceCalcGBSAOBCForceKernel::~ReferenceCalcGBSAOBCForceKernel() {
    if (obc) {
        delete obc->getObcParameters();
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
    ObcParameters* obcParameters = new ObcParameters(numParticles, ObcParameters::ObcTypeII);
    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
    obcParameters->setSolventDielectric(static_cast<RealOpenMM>(force.getSolventDielectric()));
    obcParameters->setSoluteDielectric(static_cast<RealOpenMM>(force.getSoluteDielectric()));
    obcParameters->setPi4Asolv(4*M_PI*force.getSurfaceAreaEnergy());
    if (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff)
        obcParameters->setUseCutoff(static_cast<RealOpenMM>(force.getCutoffDistance()));
    isPeriodic = (force.getNonbondedMethod() == GBSAOBCForce::CutoffPeriodic);
    obc = new ReferenceObc(obcParameters);
    obc->setIncludeAceApproximation(true);
}

double ReferenceCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    if (isPeriodic)
        obc->getObcParameters()->setPeriodic(extractBoxVectors(context));
    return obc->computeBornEnergyForces(posData, charges, forceData);
}

void ReferenceCalcGBSAOBCForceKernel::copyParametersToContext(ContextImpl& context, const GBSAOBCForce& force) {
    int numParticles = force.getNumParticles();
    ObcParameters* obcParameters = obc->getObcParameters();
    if (numParticles != obcParameters->getAtomicRadii().size())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    vector<RealOpenMM> atomicRadii(numParticles);
    vector<RealOpenMM> scaleFactors(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        charges[i] = (RealOpenMM) charge;
        atomicRadii[i] = (RealOpenMM) radius;
        scaleFactors[i] = (RealOpenMM) scalingFactor;
    }
    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
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
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expressions for computed values.

    valueDerivExpressions.resize(force.getNumComputedValues());
    valueGradientExpressions.resize(force.getNumComputedValues());
    valueParamDerivExpressions.resize(force.getNumComputedValues());
    set<string> particleVariables, pairVariables;
    pairVariables.insert("r");
    particleVariables.insert("x");
    particleVariables.insert("y");
    particleVariables.insert("z");
    for (int i = 0; i < numPerParticleParameters; i++) {
        particleVariables.insert(particleParameterNames[i]);
        pairVariables.insert(particleParameterNames[i]+"1");
        pairVariables.insert(particleParameterNames[i]+"2");
    }
    particleVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    pairVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, expression;
        CustomGBForce::ComputationType type;
        force.getComputedValueParameters(i, name, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        valueExpressions.push_back(ex.createCompiledExpression());
        valueTypes.push_back(type);
        valueNames.push_back(name);
        if (i == 0) {
            valueDerivExpressions[i].push_back(ex.differentiate("r").createCompiledExpression());
            validateVariables(ex.getRootNode(), pairVariables);
        }
        else {
            valueGradientExpressions[i].push_back(ex.differentiate("x").createCompiledExpression());
            valueGradientExpressions[i].push_back(ex.differentiate("y").createCompiledExpression());
            valueGradientExpressions[i].push_back(ex.differentiate("z").createCompiledExpression());
            for (int j = 0; j < i; j++)
                valueDerivExpressions[i].push_back(ex.differentiate(valueNames[j]).createCompiledExpression());
            validateVariables(ex.getRootNode(), particleVariables);
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++) {
            string param = force.getEnergyParameterDerivativeName(j);
            energyParamDerivNames.push_back(param);
            valueParamDerivExpressions[i].push_back(ex.differentiate(param).createCompiledExpression());
        }
        particleVariables.insert(name);
        pairVariables.insert(name+"1");
        pairVariables.insert(name+"2");
    }

    // Parse the expressions for energy terms.

    energyDerivExpressions.resize(force.getNumEnergyTerms());
    energyGradientExpressions.resize(force.getNumEnergyTerms());
    energyParamDerivExpressions.resize(force.getNumEnergyTerms());
    for (int i = 0; i < force.getNumEnergyTerms(); i++) {
        string expression;
        CustomGBForce::ComputationType type;
        force.getEnergyTermParameters(i, expression, type);
        Lepton::ParsedExpression ex = Lepton::Parser::parse(expression, functions).optimize();
        energyExpressions.push_back(ex.createCompiledExpression());
        energyTypes.push_back(type);
        if (type != CustomGBForce::SingleParticle)
            energyDerivExpressions[i].push_back(ex.differentiate("r").createCompiledExpression());
        for (int j = 0; j < force.getNumComputedValues(); j++) {
            if (type == CustomGBForce::SingleParticle) {
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]).createCompiledExpression());
                energyGradientExpressions[i].push_back(ex.differentiate("x").createCompiledExpression());
                energyGradientExpressions[i].push_back(ex.differentiate("y").createCompiledExpression());
                energyGradientExpressions[i].push_back(ex.differentiate("z").createCompiledExpression());
                validateVariables(ex.getRootNode(), particleVariables);
            }
            else {
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"1").createCompiledExpression());
                energyDerivExpressions[i].push_back(ex.differentiate(valueNames[j]+"2").createCompiledExpression());
                validateVariables(ex.getRootNode(), pairVariables);
            }
        }
        for (int j = 0; j < force.getNumEnergyParameterDerivatives(); j++)
            energyParamDerivExpressions[i].push_back(ex.differentiate(force.getEnergyParameterDerivativeName(j)).createCompiledExpression());
    }

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
}

double ReferenceCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    ReferenceCustomGBIxn ixn(valueExpressions, valueDerivExpressions, valueGradientExpressions, valueParamDerivExpressions, valueNames, valueTypes,
        energyExpressions, energyDerivExpressions, energyGradientExpressions, energyParamDerivExpressions, energyTypes, particleParameterNames);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (periodic)
        ixn.setPeriodic(extractBoxVectors(context));
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, extractBoxVectors(context), periodic, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    ixn.calculateIxn(numParticles, posData, particleParamArray, exclusions, globalParameters, forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomGBForceKernel::copyParametersToContext(ContextImpl& context, const CustomGBForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
}

ReferenceCalcCustomExternalForceKernel::PeriodicDistanceFunction::PeriodicDistanceFunction(RealVec** boxVectorHandle) : boxVectorHandle(boxVectorHandle) {
}

int ReferenceCalcCustomExternalForceKernel::PeriodicDistanceFunction::getNumArguments() const {
    return 6;
}

double ReferenceCalcCustomExternalForceKernel::PeriodicDistanceFunction::evaluate(const double* arguments) const {
    RealVec* boxVectors = *boxVectorHandle;
    RealVec delta = RealVec(arguments[0], arguments[1], arguments[2])-RealVec(arguments[3], arguments[4], arguments[5]);
    delta -= boxVectors[2]*floor(delta[2]/boxVectors[2][2]+0.5);
    delta -= boxVectors[1]*floor(delta[1]/boxVectors[1][1]+0.5);
    delta -= boxVectors[0]*floor(delta[0]/boxVectors[0][0]+0.5);
    return sqrt(delta.dot(delta));
}

double ReferenceCalcCustomExternalForceKernel::PeriodicDistanceFunction::evaluateDerivative(const double* arguments, const int* derivOrder) const {
    int argIndex = -1;
    for (int i = 0; i < 6; i++) {
        if (derivOrder[i] > 0) {
            if (derivOrder[i] > 1 || argIndex != -1)
                throw OpenMMException("Unsupported derivative of periodicdistance"); // Should be impossible for this to happen.
            argIndex = i;
        }
    }
    RealVec* boxVectors = *boxVectorHandle;
    RealVec delta = RealVec(arguments[0], arguments[1], arguments[2])-RealVec(arguments[3], arguments[4], arguments[5]);
    delta -= boxVectors[2]*floor(delta[2]/boxVectors[2][2]+0.5);
    delta -= boxVectors[1]*floor(delta[1]/boxVectors[1][1]+0.5);
    delta -= boxVectors[0]*floor(delta[0]/boxVectors[0][0]+0.5);
    double r = sqrt(delta.dot(delta));
    if (r == 0)
        return 0.0;    
    if (argIndex < 3)
        return delta[argIndex]/r;
    return -delta[argIndex-3]/r;
}

Lepton::CustomFunction* ReferenceCalcCustomExternalForceKernel::PeriodicDistanceFunction::clone() const {
    return new PeriodicDistanceFunction(boxVectorHandle);
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

    map<string, Lepton::CustomFunction*> functions;
    PeriodicDistanceFunction periodicDistance(&boxVectors);
    functions["periodicdistance"] = &periodicDistance;
    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    energyExpression = expression.createCompiledExpression();
    forceExpressionX = expression.differentiate("x").createCompiledExpression();
    forceExpressionY = expression.differentiate("y").createCompiledExpression();
    forceExpressionZ = expression.differentiate("z").createCompiledExpression();
    for (int i = 0; i < numParameters; i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    set<string> variables;
    variables.insert("x");
    variables.insert("y");
    variables.insert("z");
    variables.insert(parameterNames.begin(), parameterNames.end());
    variables.insert(globalParameterNames.begin(), globalParameterNames.end());
    validateVariables(expression.getRootNode(), variables);
}

double ReferenceCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    boxVectors = extractBoxVectors(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ReferenceCustomExternalIxn force(energyExpression, forceExpressionX, forceExpressionY, forceExpressionZ, parameterNames, globalParameters);
    for (int i = 0; i < numParticles; ++i)
        force.calculateForce(particles[i], posData, particleParamArray[i], forceData, includeEnergy ? &energy : NULL);
    return energy;
}

void ReferenceCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        int particle;
        vector<double> parameters;
        force.getParticleParameters(i, particle, parameters);
        if (particle != particles[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
}

ReferenceCalcCustomHbondForceKernel::~ReferenceCalcCustomHbondForceKernel() {
    disposeRealArray(donorParamArray, numDonors);
    disposeRealArray(acceptorParamArray, numAcceptors);
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomHbondForceKernel::initialize(const System& system, const CustomHbondForce& force) {

    // Record the exclusions.

    numDonors = force.getNumDonors();
    numAcceptors = force.getNumAcceptors();
    numParticles = system.getNumParticles();
    exclusions.resize(numDonors);
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int donor, acceptor;
        force.getExclusionParticles(i, donor, acceptor);
        exclusions[donor].insert(acceptor);
    }

    // Build the arrays.

    vector<vector<int> > donorParticles(numDonors);
    int numDonorParameters = force.getNumPerDonorParameters();
    donorParamArray = allocateRealArray(numDonors, numDonorParameters);
    for (int i = 0; i < numDonors; ++i) {
        vector<double> parameters;
        int d1, d2, d3;
        force.getDonorParameters(i, d1, d2, d3, parameters);
        donorParticles[i].push_back(d1);
        donorParticles[i].push_back(d2);
        donorParticles[i].push_back(d3);
        for (int j = 0; j < numDonorParameters; j++)
            donorParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    vector<vector<int> > acceptorParticles(numAcceptors);
    int numAcceptorParameters = force.getNumPerAcceptorParameters();
    acceptorParamArray = allocateRealArray(numAcceptors, numAcceptorParameters);
    for (int i = 0; i < numAcceptors; ++i) {
        vector<double> parameters;
        int a1, a2, a3;
        force.getAcceptorParameters(i, a1, a2, a3, parameters);
        acceptorParticles[i].push_back(a1);
        acceptorParticles[i].push_back(a2);
        acceptorParticles[i].push_back(a3);
        for (int j = 0; j < numAcceptorParameters; j++)
            acceptorParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    NonbondedMethod nonbondedMethod = CalcCustomHbondForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the object used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomHbondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    vector<string> donorParameterNames;
    vector<string> acceptorParameterNames;
    for (int i = 0; i < numDonorParameters; i++)
        donorParameterNames.push_back(force.getPerDonorParameterName(i));
    for (int i = 0; i < numAcceptorParameters; i++)
        acceptorParameterNames.push_back(force.getPerAcceptorParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    ixn = new ReferenceCustomHbondIxn(donorParticles, acceptorParticles, energyExpression, donorParameterNames, acceptorParameterNames, distances, angles, dihedrals);
    isPeriodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff)
        ixn->setUseCutoff(nonbondedCutoff);

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
}

double ReferenceCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    if (isPeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    ixn->calculatePairIxn(posData, donorParamArray, acceptorParamArray, exclusions, globalParameters, forceData, includeEnergy ? &energy : NULL);
    return energy;
}

void ReferenceCalcCustomHbondForceKernel::copyParametersToContext(ContextImpl& context, const CustomHbondForce& force) {
    if (numDonors != force.getNumDonors())
        throw OpenMMException("updateParametersInContext: The number of donors has changed");
    if (numAcceptors != force.getNumAcceptors())
        throw OpenMMException("updateParametersInContext: The number of acceptors has changed");

    // Record the values.

    vector<double> parameters;
    int numDonorParameters = force.getNumPerDonorParameters();
    const vector<vector<int> >& donorAtoms = ixn->getDonorAtoms();
    for (int i = 0; i < numDonors; ++i) {
        int d1, d2, d3;
        force.getDonorParameters(i, d1, d2, d3, parameters);
        if (d1 != donorAtoms[i][0] || d2 != donorAtoms[i][1] || d3 != donorAtoms[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in a donor group has changed");
        for (int j = 0; j < numDonorParameters; j++)
            donorParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
    int numAcceptorParameters = force.getNumPerAcceptorParameters();
    const vector<vector<int> >& acceptorAtoms = ixn->getAcceptorAtoms();
    for (int i = 0; i < numAcceptors; ++i) {
        int a1, a2, a3;
        force.getAcceptorParameters(i, a1, a2, a3, parameters);
        if (a1 != acceptorAtoms[i][0] || a2 != acceptorAtoms[i][1] || a3 != acceptorAtoms[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an acceptor group has changed");
        for (int j = 0; j < numAcceptorParameters; j++)
            acceptorParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
}

ReferenceCalcCustomCentroidBondForceKernel::~ReferenceCalcCustomCentroidBondForceKernel() {
    disposeRealArray(bondParamArray, numBonds);
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomCentroidBondForceKernel::initialize(const System& system, const CustomCentroidBondForce& force) {
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    int numGroups = force.getNumGroups();
    vector<vector<int> > groupAtoms(numGroups);
    vector<double> ignored;
    for (int i = 0; i < numGroups; i++)
        force.getGroupParameters(i, groupAtoms[i], ignored);
    vector<vector<double> > normalizedWeights;
    CustomCentroidBondForceImpl::computeNormalizedWeights(force, system, normalizedWeights);
    numBonds = force.getNumBonds();
    vector<vector<int> > bondGroups(numBonds);
    int numBondParameters = force.getNumPerBondParameters();
    bondParamArray = allocateRealArray(numBonds, numBondParameters);
    for (int i = 0; i < numBonds; ++i) {
        vector<double> parameters;
        force.getBondParameters(i, bondGroups[i], parameters);
        for (int j = 0; j < numBondParameters; j++)
            bondParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the object used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomCentroidBondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    vector<string> bondParameterNames;
    for (int i = 0; i < numBondParameters; i++)
        bondParameterNames.push_back(force.getPerBondParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(energyExpression.differentiate(param).createCompiledExpression());
    }
    ixn = new ReferenceCustomCentroidBondIxn(force.getNumGroupsPerBond(), groupAtoms, normalizedWeights, bondGroups, energyExpression, bondParameterNames, distances, angles, dihedrals, energyParamDerivExpressions);

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
}

double ReferenceCalcCustomCentroidBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    if (usePeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    ixn->calculatePairIxn(posData, bondParamArray, globalParameters, forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomCentroidBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCentroidBondForce& force) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    int numParameters = force.getNumPerBondParameters();
    const vector<vector<int> >& bondGroups = ixn->getBondGroups();
    vector<int> groups;
    vector<double> params;
    for (int i = 0; i < numBonds; ++i) {
        force.getBondParameters(i, groups, params);
        for (int j = 0; j < groups.size(); j++)
            if (groups[j] != bondGroups[i][j])
                throw OpenMMException("updateParametersInContext: The set of groups in a bond has changed");
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = (RealOpenMM) params[j];
    }
}

ReferenceCalcCustomCompoundBondForceKernel::~ReferenceCalcCustomCompoundBondForceKernel() {
    disposeRealArray(bondParamArray, numBonds);
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    numBonds = force.getNumBonds();
    vector<vector<int> > bondParticles(numBonds);
    int numBondParameters = force.getNumPerBondParameters();
    bondParamArray = allocateRealArray(numBonds, numBondParameters);
    for (int i = 0; i < numBonds; ++i) {
        vector<double> parameters;
        force.getBondParameters(i, bondParticles[i], parameters);
        for (int j = 0; j < numBondParameters; j++)
            bondParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }

    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the object used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomCompoundBondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    vector<string> bondParameterNames;
    for (int i = 0; i < numBondParameters; i++)
        bondParameterNames.push_back(force.getPerBondParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(energyExpression.differentiate(param).createCompiledExpression());
    }
    ixn = new ReferenceCustomCompoundBondIxn(force.getNumParticlesPerBond(), bondParticles, energyExpression, bondParameterNames, distances, angles, dihedrals, energyParamDerivExpressions);

    // Delete the custom functions.

    for (map<string, Lepton::CustomFunction*>::iterator iter = functions.begin(); iter != functions.end(); iter++)
        delete iter->second;
}

double ReferenceCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    if (usePeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    ixn->calculatePairIxn(posData, bondParamArray, globalParameters, forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomCompoundBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomCompoundBondForce& force) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    int numParameters = force.getNumPerBondParameters();
    const vector<vector<int> >& bondAtoms = ixn->getBondAtoms();
    vector<int> particles;
    vector<double> params;
    for (int i = 0; i < numBonds; ++i) {
        force.getBondParameters(i, particles, params);
        for (int j = 0; j < particles.size(); j++)
            if (particles[j] != bondAtoms[i][j])
                throw OpenMMException("updateParametersInContext: The set of particles in a bond has changed");
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = (RealOpenMM) params[j];
    }
}

ReferenceCalcCustomManyParticleForceKernel::~ReferenceCalcCustomManyParticleForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {

    // Build the arrays.

    numParticles = system.getNumParticles();
    int numParticleParameters = force.getNumPerParticleParameters();
    particleParamArray = allocateRealArray(numParticles, numParticleParameters);
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        int type;
        force.getParticleParameters(i, parameters, type);
        for (int j = 0; j < numParticleParameters; j++)
            particleParamArray[i][j] = parameters[j];
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    ixn = new ReferenceCustomManyParticleIxn(force);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    cutoffDistance = force.getCutoffDistance();
}

double ReferenceCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    RealOpenMM energy = 0;
    map<string, double> globalParameters;
    for (int i = 0; i < (int) globalParameterNames.size(); i++)
        globalParameters[globalParameterNames[i]] = context.getParameter(globalParameterNames[i]);
    if (nonbondedMethod == CutoffPeriodic) {
        RealVec* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 2*cutoffDistance;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        ixn->setPeriodic(boxVectors);
    }
    ixn->calculateIxn(posData, particleParamArray, globalParameters, forceData, includeEnergy ? &energy : NULL);
    return energy;
}

void ReferenceCalcCustomManyParticleForceKernel::copyParametersToContext(ContextImpl& context, const CustomManyParticleForce& force) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = 0; i < numParticles; ++i) {
        vector<double> parameters;
        int type;
        force.getParticleParameters(i, parameters, type);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = static_cast<RealOpenMM>(parameters[j]);
    }
}

ReferenceCalcGayBerneForceKernel::~ReferenceCalcGayBerneForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcGayBerneForceKernel::initialize(const System& system, const GayBerneForce& force) {
    ixn = new ReferenceGayBerneForce(force);
}

double ReferenceCalcGayBerneForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return ixn->calculateForce(extractPositions(context), extractForces(context), extractBoxVectors(context));
}

void ReferenceCalcGayBerneForceKernel::copyParametersToContext(ContextImpl& context, const GayBerneForce& force) {
    delete ixn;
    ixn = NULL;
    ixn = new ReferenceGayBerneForce(force);
}

ReferenceIntegrateVerletStepKernel::~ReferenceIntegrateVerletStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
}

void ReferenceIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    double stepSize = integrator.getStepSize();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceVerletDynamics(context.getSystem().getNumParticles(), static_cast<RealOpenMM>(stepSize));
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

ReferenceIntegrateLangevinStepKernel::~ReferenceIntegrateLangevinStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateLangevinStepKernel::initialize(const System& system, const LangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateLangevinStepKernel::execute(ContextImpl& context, const LangevinIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        RealOpenMM tau = static_cast<RealOpenMM>(friction == 0.0 ? 0.0 : 1.0/friction);
        dynamics = new ReferenceStochasticDynamics(
                context.getSystem().getNumParticles(), 
                static_cast<RealOpenMM>(stepSize), 
                static_cast<RealOpenMM>(tau), 
                static_cast<RealOpenMM>(temperature));
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

ReferenceIntegrateBrownianStepKernel::~ReferenceIntegrateBrownianStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceBrownianDynamics(
                context.getSystem().getNumParticles(), 
                static_cast<RealOpenMM>(stepSize), 
                static_cast<RealOpenMM>(friction), 
                static_cast<RealOpenMM>(temperature));
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateBrownianStepKernel::computeKineticEnergy(ContextImpl& context, const BrownianIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0);
}

ReferenceIntegrateVariableLangevinStepKernel::~ReferenceIntegrateVariableLangevinStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateVariableLangevinStepKernel::initialize(const System& system, const VariableLangevinIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

double ReferenceIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double errorTol = integrator.getErrorTolerance();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || errorTol != prevErrorTol) {
        // Recreate the computation objects with the new parameters.

        if (dynamics)
            delete dynamics;
        RealOpenMM tau = static_cast<RealOpenMM>(friction == 0.0 ? 0.0 : 1.0/friction);
        dynamics = new ReferenceVariableStochasticDynamics(context.getSystem().getNumParticles(), (RealOpenMM) tau, (RealOpenMM) temperature, (RealOpenMM) errorTol);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevErrorTol = errorTol;
    }
    RealOpenMM maxStepSize = (RealOpenMM) (maxTime-data.time);
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, maxStepSize, integrator.getConstraintTolerance());
    data.time += dynamics->getDeltaT();
    if (dynamics->getDeltaT() == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
    return dynamics->getDeltaT();
}

double ReferenceIntegrateVariableLangevinStepKernel::computeKineticEnergy(ContextImpl& context, const VariableLangevinIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

ReferenceIntegrateVariableVerletStepKernel::~ReferenceIntegrateVariableVerletStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateVariableVerletStepKernel::initialize(const System& system, const VariableVerletIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
}

double ReferenceIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    double errorTol = integrator.getErrorTolerance();
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    if (dynamics == 0 || errorTol != prevErrorTol) {
        // Recreate the computation objects with the new parameters.

        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceVariableVerletDynamics(context.getSystem().getNumParticles(), (RealOpenMM) errorTol);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevErrorTol = errorTol;
    }
    RealOpenMM maxStepSize = (RealOpenMM) (maxTime-data.time);
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, maxStepSize, integrator.getConstraintTolerance());
    data.time += dynamics->getDeltaT();
    if (dynamics->getDeltaT() == maxStepSize)
        data.time = maxTime; // Avoid round-off error
    data.stepCount++;
    return dynamics->getDeltaT();
}

double ReferenceIntegrateVariableVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VariableVerletIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

ReferenceIntegrateCustomStepKernel::~ReferenceIntegrateCustomStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateCustomStepKernel::initialize(const System& system, const CustomIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    perDofValues.resize(integrator.getNumPerDofVariables());
    for (int i = 0; i < (int) perDofValues.size(); i++)
        perDofValues[i].resize(numParticles);

    // Create the computation objects.

    dynamics = new ReferenceCustomDynamics(system.getNumParticles(), integrator);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    
    // Record global variables.
    
    map<string, double> globals;
    globals["dt"] = integrator.getStepSize();
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        globals[integrator.getGlobalVariableName(i)] = globalValues[i];
    
    // Execute the step.
    
    dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
    dynamics->update(context, context.getSystem().getNumParticles(), posData, velData, forceData, masses, globals, perDofValues, forcesAreValid, integrator.getConstraintTolerance());
    
    // Record changed global variables.
    
    integrator.setStepSize(globals["dt"]);
    for (int i = 0; i < (int) globalValues.size(); i++)
        globalValues[i] = globals[integrator.getGlobalVariableName(i)];
    data.time += dynamics->getDeltaT();
    data.stepCount++;
}

double ReferenceIntegrateCustomStepKernel::computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    vector<RealVec>& posData = extractPositions(context);
    vector<RealVec>& velData = extractVelocities(context);
    vector<RealVec>& forceData = extractForces(context);
    
    // Record global variables.
    
    map<string, double> globals;
    globals["dt"] = integrator.getStepSize();
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        globals[integrator.getGlobalVariableName(i)] = globalValues[i];
    
    // Compute the kinetic energy.
    
    return dynamics->computeKineticEnergy(context, context.getSystem().getNumParticles(), posData, velData, forceData, masses, globals, perDofValues, forcesAreValid);
}

void ReferenceIntegrateCustomStepKernel::getGlobalVariables(ContextImpl& context, vector<double>& values) const {
    values = globalValues;
}

void ReferenceIntegrateCustomStepKernel::setGlobalVariables(ContextImpl& context, const vector<double>& values) {
    globalValues = values;
}

void ReferenceIntegrateCustomStepKernel::getPerDofVariable(ContextImpl& context, int variable, vector<Vec3>& values) const {
    values.resize(perDofValues[variable].size());
    for (int i = 0; i < (int) values.size(); i++)
        values[i] = perDofValues[variable][i];
}

void ReferenceIntegrateCustomStepKernel::setPerDofVariable(ContextImpl& context, int variable, const vector<Vec3>& values) {
    perDofValues[variable].resize(values.size());
    for (int i = 0; i < (int) values.size(); i++)
        perDofValues[variable][i] = values[i];
}

ReferenceApplyAndersenThermostatKernel::~ReferenceApplyAndersenThermostatKernel() {
    if (thermostat)
        delete thermostat;
}

void ReferenceApplyAndersenThermostatKernel::initialize(const System& system, const AndersenThermostat& thermostat) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = static_cast<RealOpenMM>(system.getParticleMass(i));
    this->thermostat = new ReferenceAndersenThermostat();
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) thermostat.getRandomNumberSeed());
    particleGroups = AndersenThermostatImpl::calcParticleGroups(system);
}

void ReferenceApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    vector<RealVec>& velData = extractVelocities(context);
    thermostat->applyThermostat(particleGroups, velData, masses,
        static_cast<RealOpenMM>(context.getParameter(AndersenThermostat::Temperature())),
        static_cast<RealOpenMM>(context.getParameter(AndersenThermostat::CollisionFrequency())),
        static_cast<RealOpenMM>(context.getIntegrator().getStepSize()));
}

ReferenceApplyMonteCarloBarostatKernel::~ReferenceApplyMonteCarloBarostatKernel() {
    if (barostat)
        delete barostat;
}

void ReferenceApplyMonteCarloBarostatKernel::initialize(const System& system, const Force& barostat) {
}

void ReferenceApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ) {
    if (barostat == NULL)
        barostat = new ReferenceMonteCarloBarostat(context.getSystem().getNumParticles(), context.getMolecules());
    vector<RealVec>& posData = extractPositions(context);
    RealVec* boxVectors = extractBoxVectors(context);
    barostat->applyBarostat(posData, boxVectors, scaleX, scaleY, scaleZ);
}

void ReferenceApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    vector<RealVec>& posData = extractPositions(context);
    barostat->restorePositions(posData);
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
    vector<RealVec>& velData = extractVelocities(context);
    
    // Calculate the center of mass momentum.
    
    RealOpenMM momentum[] = {0.0, 0.0, 0.0};
    RealOpenMM mass = 0.0;
    for (size_t i = 0; i < masses.size(); ++i) {
        momentum[0] += static_cast<RealOpenMM>(masses[i]*velData[i][0]);
        momentum[1] += static_cast<RealOpenMM>(masses[i]*velData[i][1]);
        momentum[2] += static_cast<RealOpenMM>(masses[i]*velData[i][2]);
        mass += static_cast<RealOpenMM>(masses[i]);
    }
    
    // Adjust the particle velocities.
    
    momentum[0] /= mass;
    momentum[1] /= mass;
    momentum[2] /= mass;
    for (size_t i = 0; i < masses.size(); ++i) {
        if (masses[i] != 0.0) {
            velData[i][0] -= momentum[0];
            velData[i][1] -= momentum[1];
            velData[i][2] -= momentum[2];
        }
    }
}
