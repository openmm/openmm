/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2024 Stanford University and the Authors.      *
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
#include "ReferenceCMAPTorsionIxn.h"
#include "ReferenceConstraints.h"
#include "ReferenceCustomAngleIxn.h"
#include "ReferenceCustomBondIxn.h"
#include "ReferenceCustomCentroidBondIxn.h"
#include "ReferenceCustomCompoundBondIxn.h"
#include "ReferenceCustomCVForce.h"
#include "ReferenceCustomDynamics.h"
#include "ReferenceCustomExternalIxn.h"
#include "ReferenceCustomGBIxn.h"
#include "ReferenceCustomHbondIxn.h"
#include "ReferenceCustomNonbondedIxn.h"
#include "ReferenceCustomManyParticleIxn.h"
#include "ReferenceCustomTorsionIxn.h"
#include "ReferenceGayBerneForce.h"
#include "ReferenceHarmonicBondIxn.h"
#include "ReferenceLangevinMiddleDynamics.h"
#include "ReferenceLJCoulomb14.h"
#include "ReferenceLJCoulombIxn.h"
#include "ReferenceMonteCarloBarostat.h"
#include "ReferenceNoseHooverChain.h"
#include "ReferenceNoseHooverDynamics.h"
#include "ReferencePointFunctions.h"
#include "ReferenceProperDihedralBond.h"
#include "ReferenceRbDihedralBond.h"
#include "ReferenceRMSDForce.h"
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

static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->positions;
}

static vector<Vec3>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->velocities;
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->forces;
}

static Vec3& extractBoxSize(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->periodicBoxSize;
}

static Vec3* extractBoxVectors(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return data->periodicBoxVectors;
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->constraints;
}

static const ReferenceVirtualSites& extractVirtualSites(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->virtualSites;
}

static map<string, double>& extractEnergyParameterDerivatives(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *data->energyParameterDerivatives;
}

/**
 * Make sure an expression doesn't use any undefined variables.
 */
static void validateVariables(const Lepton::ExpressionTreeNode& node, const set<string>& variables) {
    const Lepton::Operation& op = node.getOperation();
    if (op.getId() == Lepton::Operation::VARIABLE && variables.find(op.getName()) == variables.end())
        throw OpenMMException("Unknown variable in expression: "+op.getName());
    for (auto& child : node.getChildren())
        validateVariables(child, variables);
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3> shiftedVel(numParticles);
    context.computeShiftedVelocities(timeShift, shiftedVel);
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}

void ReferenceCalcForcesAndEnergyKernel::initialize(const System& system) {
}

void ReferenceCalcForcesAndEnergyKernel::beginComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    vector<Vec3>& forceData = extractForces(context);
    if (includeForces) {
        int numParticles = context.getSystem().getNumParticles();
        for (int i = 0; i < numParticles; ++i) {
            forceData[i][0] = 0.0;
            forceData[i][1] = 0.0;
            forceData[i][2] = 0.0;
        }
    }
    else
        savedForces = forceData;
    for (auto& param : context.getParameters())
        extractEnergyParameterDerivatives(context)[param.first] = 0;
}

double ReferenceCalcForcesAndEnergyKernel::finishComputation(ContextImpl& context, bool includeForces, bool includeEnergy, int groups, bool& valid) {
    if (!includeForces)
        extractForces(context) = savedForces; // Restore the forces so computing the energy doesn't overwrite the forces with incorrect values.
    else
        extractVirtualSites(context).distributeForces(context.getSystem(), extractPositions(context), extractForces(context));
    return 0.0;
}

void ReferenceUpdateStateDataKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

double ReferenceUpdateStateDataKernel::getTime(const ContextImpl& context) const {
    return data.time;
}

void ReferenceUpdateStateDataKernel::setTime(ContextImpl& context, double time) {
    data.time = time;
}

long long ReferenceUpdateStateDataKernel::getStepCount(const ContextImpl& context) const {
    return data.stepCount;
}

void ReferenceUpdateStateDataKernel::setStepCount(const ContextImpl& context, long long count) {
    data.stepCount = count;
}

void ReferenceUpdateStateDataKernel::getPositions(ContextImpl& context, std::vector<Vec3>& positions) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3>& posData = extractPositions(context);
    positions.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        positions[i] = Vec3(posData[i][0], posData[i][1], posData[i][2]);
}

void ReferenceUpdateStateDataKernel::setPositions(ContextImpl& context, const std::vector<Vec3>& positions) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3>& posData = extractPositions(context);
    for (int i = 0; i < numParticles; ++i) {
        posData[i][0] = positions[i][0];
        posData[i][1] = positions[i][1];
        posData[i][2] = positions[i][2];
    }
}

void ReferenceUpdateStateDataKernel::getVelocities(ContextImpl& context, std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3>& velData = extractVelocities(context);
    velocities.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        velocities[i] = Vec3(velData[i][0], velData[i][1], velData[i][2]);
}

void ReferenceUpdateStateDataKernel::setVelocities(ContextImpl& context, const std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3>& velData = extractVelocities(context);
    for (int i = 0; i < numParticles; ++i) {
        velData[i][0] = velocities[i][0];
        velData[i][1] = velocities[i][1];
        velData[i][2] = velocities[i][2];
    }
}

void ReferenceUpdateStateDataKernel::computeShiftedVelocities(ContextImpl& context, double timeShift, std::vector<Vec3>& velocities) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    
    // Compute the shifted velocities.
    
    velocities.resize(numParticles);
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (masses[i] == 0) {
            velocities[i] = velData[i];
            inverseMasses[i] = 0;
        }
        else {
            velocities[i] = velData[i]+forceData[i]*(timeShift/masses[i]);
            inverseMasses[i] = 1/masses[i];
        }
    }
    
    // Apply constraints to them.

    if (timeShift != 0)
        extractConstraints(context).applyToVelocities(posData, velocities, inverseMasses, 1e-4);
}

void ReferenceUpdateStateDataKernel::getForces(ContextImpl& context, std::vector<Vec3>& forces) {
    int numParticles = context.getSystem().getNumParticles();
    vector<Vec3>& forceData = extractForces(context);
    forces.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        forces[i] = Vec3(forceData[i][0], forceData[i][1], forceData[i][2]);
}

void ReferenceUpdateStateDataKernel::getEnergyParameterDerivatives(ContextImpl& context, map<string, double>& derivs) {
    derivs = extractEnergyParameterDerivatives(context);
}

void ReferenceUpdateStateDataKernel::getPeriodicBoxVectors(ContextImpl& context, Vec3& a, Vec3& b, Vec3& c) const {
    Vec3* vectors = extractBoxVectors(context);
    a = vectors[0];
    b = vectors[1];
    c = vectors[2];
}

void ReferenceUpdateStateDataKernel::setPeriodicBoxVectors(ContextImpl& context, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3& box = extractBoxSize(context);
    box[0] = a[0];
    box[1] = b[1];
    box[2] = c[2];
    Vec3* vectors = extractBoxVectors(context);
    vectors[0] = a;
    vectors[1] = b;
    vectors[2] = c;
}

void ReferenceUpdateStateDataKernel::createCheckpoint(ContextImpl& context, ostream& stream) {
    int version = 3;
    stream.write((char*) &version, sizeof(int));
    stream.write((char*) &data.time, sizeof(data.time));
    stream.write((char*) &data.stepCount, sizeof(long long));
    vector<Vec3>& posData = extractPositions(context);
    stream.write((char*) &posData[0], sizeof(Vec3)*posData.size());
    vector<Vec3>& velData = extractVelocities(context);
    stream.write((char*) &velData[0], sizeof(Vec3)*velData.size());
    Vec3* vectors = extractBoxVectors(context);
    stream.write((char*) vectors, 3*sizeof(Vec3));
    SimTKOpenMMUtilities::createCheckpoint(stream);
}

void ReferenceUpdateStateDataKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    int version;
    stream.read((char*) &version, sizeof(int));
    if (version != 3)
        throw OpenMMException("Checkpoint was created with a different version of OpenMM");
    stream.read((char*) &data.time, sizeof(data.time));
    stream.read((char*) &data.stepCount, sizeof(long long));
    vector<Vec3>& posData = extractPositions(context);
    stream.read((char*) &posData[0], sizeof(Vec3)*posData.size());
    vector<Vec3>& velData = extractVelocities(context);
    stream.read((char*) &velData[0], sizeof(Vec3)*velData.size());
    Vec3* vectors = extractBoxVectors(context);
    stream.read((char*) vectors, 3*sizeof(Vec3));
    SimTKOpenMMUtilities::loadCheckpoint(stream);
}

void ReferenceApplyConstraintsKernel::initialize(const System& system) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    inverseMasses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        masses[i] = system.getParticleMass(i);
        inverseMasses[i] = 1.0/masses[i];
    }
}

ReferenceApplyConstraintsKernel::~ReferenceApplyConstraintsKernel() {
}

void ReferenceApplyConstraintsKernel::apply(ContextImpl& context, double tol) {
    vector<Vec3>& positions = extractPositions(context);
    extractConstraints(context).apply(positions, positions, inverseMasses, tol);
    extractVirtualSites(context).computePositions(context.getSystem(), positions);
}

void ReferenceApplyConstraintsKernel::applyToVelocities(ContextImpl& context, double tol) {
    vector<Vec3>& positions = extractPositions(context);
    vector<Vec3>& velocities = extractVelocities(context);
    extractConstraints(context).applyToVelocities(positions, velocities, inverseMasses, tol);
}

void ReferenceVirtualSitesKernel::initialize(const System& system) {
}

void ReferenceVirtualSitesKernel::computePositions(ContextImpl& context) {
    vector<Vec3>& positions = extractPositions(context);
    extractVirtualSites(context).computePositions(context.getSystem(), positions);
}

void ReferenceCalcHarmonicBondForceKernel::initialize(const System& system, const HarmonicBondForce& force) {
    numBonds = force.getNumBonds();
    bondIndexArray.resize(numBonds, vector<int>(2));
    bondParamArray.resize(numBonds, vector<double>(2));
    for (int i = 0; i < numBonds; ++i) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        bondParamArray[i][0] = length;
        bondParamArray[i][1] = k;
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcHarmonicBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    if (usePeriodic)
        harmonicBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, includeEnergy ? &energy : NULL, harmonicBond);
    return energy;
}

void ReferenceCalcHarmonicBondForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicBondForce& force, int firstBond, int lastBond) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    for (int i = firstBond; i <= lastBond; ++i) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(i, particle1, particle2, length, k);
        if (particle1 != bondIndexArray[i][0] || particle2 != bondIndexArray[i][1])
            throw OpenMMException("updateParametersInContext: The set of particles in a bond has changed");
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        bondParamArray[i][0] = length;
        bondParamArray[i][1] = k;
    }
}

ReferenceCalcCustomBondForceKernel::~ReferenceCalcCustomBondForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomBondForceKernel::initialize(const System& system, const CustomBondForce& force) {
    numBonds = force.getNumBonds();
    int numParameters = force.getNumPerBondParameters();
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    bondIndexArray.resize(numBonds, vector<int>(2));
    bondParamArray.resize(numBonds, vector<double>(numParameters));
    vector<double> params;
    for (int i = 0; i < numBonds; ++i) {
        int particle1, particle2;
        force.getBondParameters(i, particle1, particle2, params);
        bondIndexArray[i][0] = particle1;
        bondIndexArray[i][1] = particle2;
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = params[j];
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
    ixn = new ReferenceCustomBondIxn(energyExpression, forceExpression, parameterNames, energyParamDerivExpressions);
}

double ReferenceCalcCustomBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    ixn->setGlobalParameters(globalParameters);
    if (usePeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    for (int i = 0; i < numBonds; i++)
        ixn->calculateBondIxn(bondIndexArray[i], posData, bondParamArray[i], forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomBondForceKernel::copyParametersToContext(ContextImpl& context, const CustomBondForce& force, int firstBond, int lastBond) {
    if (numBonds != force.getNumBonds())
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");

    // Record the values.

    int numParameters = force.getNumPerBondParameters();
    vector<double> params;
    for (int i = firstBond; i <= lastBond; ++i) {
        int particle1, particle2;
        force.getBondParameters(i, particle1, particle2, params);
        if (particle1 != bondIndexArray[i][0] || particle2 != bondIndexArray[i][1])
            throw OpenMMException("updateParametersInContext: The set of particles in a bond has changed");
        for (int j = 0; j < numParameters; j++)
            bondParamArray[i][j] = params[j];
    }
}

void ReferenceCalcHarmonicAngleForceKernel::initialize(const System& system, const HarmonicAngleForce& force) {
    numAngles = force.getNumAngles();
    angleIndexArray.resize(numAngles, vector<int>(3));
    angleParamArray.resize(numAngles, vector<double>(2));
    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        angleParamArray[i][0] = angle;
        angleParamArray[i][1] = k;
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcHarmonicAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceAngleBondIxn angleBond;
    if (usePeriodic)
        angleBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, includeEnergy ? &energy : NULL, angleBond);
    return energy;
}

void ReferenceCalcHarmonicAngleForceKernel::copyParametersToContext(ContextImpl& context, const HarmonicAngleForce& force, int firstAngle, int lastAngle) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    for (int i = firstAngle; i <= lastAngle; ++i) {
        int particle1, particle2, particle3;
        double angle, k;
        force.getAngleParameters(i, particle1, particle2, particle3, angle, k);
        if (particle1 != angleIndexArray[i][0] || particle2 != angleIndexArray[i][1] || particle3 != angleIndexArray[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        angleParamArray[i][0] = angle;
        angleParamArray[i][1] = k;
    }
}

ReferenceCalcCustomAngleForceKernel::~ReferenceCalcCustomAngleForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomAngleForceKernel::initialize(const System& system, const CustomAngleForce& force) {
    numAngles = force.getNumAngles();
    int numParameters = force.getNumPerAngleParameters();
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    angleIndexArray.resize(numAngles, vector<int>(3));
    angleParamArray.resize(numAngles, vector<double>(numParameters));
    vector<double> params;
    for (int i = 0; i < numAngles; ++i) {
        int particle1, particle2, particle3;
        force.getAngleParameters(i, particle1, particle2, particle3, params);
        angleIndexArray[i][0] = particle1;
        angleIndexArray[i][1] = particle2;
        angleIndexArray[i][2] = particle3;
        for (int j = 0; j < numParameters; j++)
            angleParamArray[i][j] = params[j];
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
    ixn = new ReferenceCustomAngleIxn(energyExpression, forceExpression, parameterNames, energyParamDerivExpressions);
}

double ReferenceCalcCustomAngleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    ixn->setGlobalParameters(globalParameters);
    if (usePeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    for (int i = 0; i < numAngles; i++)
        ixn->calculateBondIxn(angleIndexArray[i], posData, angleParamArray[i], forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomAngleForceKernel::copyParametersToContext(ContextImpl& context, const CustomAngleForce& force, int firstAngle, int lastAngle) {
    if (numAngles != force.getNumAngles())
        throw OpenMMException("updateParametersInContext: The number of angles has changed");

    // Record the values.

    int numParameters = force.getNumPerAngleParameters();
    vector<double> params;
    for (int i = firstAngle; i <= lastAngle; ++i) {
        int particle1, particle2, particle3;
        force.getAngleParameters(i, particle1, particle2, particle3, params);
        if (particle1 != angleIndexArray[i][0] || particle2 != angleIndexArray[i][1] || particle3 != angleIndexArray[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an angle has changed");
        for (int j = 0; j < numParameters; j++)
            angleParamArray[i][j] = params[j];
    }
}

void ReferenceCalcPeriodicTorsionForceKernel::initialize(const System& system, const PeriodicTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray.resize(numTorsions, vector<int>(4));
    torsionParamArray.resize(numTorsions, vector<double>(3));
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = k;
        torsionParamArray[i][1] = phase;
        torsionParamArray[i][2] = periodicity;
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcPeriodicTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceBondForce refBondForce;
    ReferenceProperDihedralBond periodicTorsionBond;
    if (usePeriodic)
        periodicTorsionBond.setPeriodic(extractBoxVectors(context));
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, includeEnergy ? &energy : NULL, periodicTorsionBond);
    return energy;
}

void ReferenceCalcPeriodicTorsionForceKernel::copyParametersToContext(ContextImpl& context, const PeriodicTorsionForce& force, int firstTorsion, int lastTorsion) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    for (int i = firstTorsion; i <= lastTorsion; ++i) {
        int particle1, particle2, particle3, particle4, periodicity;
        double phase, k;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, periodicity, phase, k);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        torsionParamArray[i][0] = k;
        torsionParamArray[i][1] = phase;
        torsionParamArray[i][2] = periodicity;
    }
}

void ReferenceCalcRBTorsionForceKernel::initialize(const System& system, const RBTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    torsionIndexArray.resize(numTorsions, vector<int>(4));
    torsionParamArray.resize(numTorsions, vector<double>(6));
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        double c0, c1, c2, c3, c4, c5;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        torsionParamArray[i][0] = c0;
        torsionParamArray[i][1] = c1;
        torsionParamArray[i][2] = c2;
        torsionParamArray[i][3] = c3;
        torsionParamArray[i][4] = c4;
        torsionParamArray[i][5] = c5;
    }
    usePeriodic = force.usesPeriodicBoundaryConditions();
}

double ReferenceCalcRBTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
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
        torsionParamArray[i][0] = c0;
        torsionParamArray[i][1] = c1;
        torsionParamArray[i][2] = c2;
        torsionParamArray[i][3] = c3;
        torsionParamArray[i][4] = c4;
        torsionParamArray[i][5] = c5;
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
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double totalEnergy = 0;
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
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomTorsionForceKernel::initialize(const System& system, const CustomTorsionForce& force) {
    numTorsions = force.getNumTorsions();
    int numParameters = force.getNumPerTorsionParameters();
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    torsionIndexArray.resize(numTorsions, vector<int>(4));
    torsionParamArray.resize(numTorsions, vector<double>(numParameters));
    vector<double> params;
    for (int i = 0; i < numTorsions; ++i) {
        int particle1, particle2, particle3, particle4;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, params);
        torsionIndexArray[i][0] = particle1;
        torsionIndexArray[i][1] = particle2;
        torsionIndexArray[i][2] = particle3;
        torsionIndexArray[i][3] = particle4;
        for (int j = 0; j < numParameters; j++)
            torsionParamArray[i][j] = params[j];
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
    ixn = new ReferenceCustomTorsionIxn(energyExpression, forceExpression, parameterNames, energyParamDerivExpressions);
}

double ReferenceCalcCustomTorsionForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    ixn->setGlobalParameters(globalParameters);
    if (usePeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    for (int i = 0; i < numTorsions; i++)
        ixn->calculateBondIxn(torsionIndexArray[i], posData, torsionParamArray[i], forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    return energy;
}

void ReferenceCalcCustomTorsionForceKernel::copyParametersToContext(ContextImpl& context, const CustomTorsionForce& force, int firstTorsion, int lastTorsion) {
    if (numTorsions != force.getNumTorsions())
        throw OpenMMException("updateParametersInContext: The number of torsions has changed");

    // Record the values.

    int numParameters = force.getNumPerTorsionParameters();
    vector<double> params;
    for (int i = firstTorsion; i <= lastTorsion; ++i) {
        int particle1, particle2, particle3, particle4;
        force.getTorsionParameters(i, particle1, particle2, particle3, particle4, params);
        if (particle1 != torsionIndexArray[i][0] || particle2 != torsionIndexArray[i][1] || particle3 != torsionIndexArray[i][2] || particle4 != torsionIndexArray[i][3])
            throw OpenMMException("updateParametersInContext: The set of particles in a torsion has changed");
        for (int j = 0; j < numParameters; j++)
            torsionParamArray[i][j] = params[j];
    }
}

ReferenceCalcNonbondedForceKernel::~ReferenceCalcNonbondedForceKernel() {
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force) {

    // Identify which exceptions are 1-4 interactions.

    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    vector<int> nb14s;
    map<int, int> nb14Index;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end()) {
            nb14Index[i] = nb14s.size();
            nb14s.push_back(i);
        }
    }

    // Build the arrays.

    num14 = nb14s.size();
    bonded14IndexArray.resize(num14, vector<int>(2));
    bonded14ParamArray.resize(num14, vector<double>(3));
    particleParamArray.resize(numParticles, vector<double>(3));
    baseParticleParams.resize(numParticles);
    baseExceptionParams.resize(num14);
    for (int i = 0; i < numParticles; ++i)
       force.getParticleParameters(i, baseParticleParams[i][0], baseParticleParams[i][1], baseParticleParams[i][2]);
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        force.getExceptionParameters(nb14s[i], particle1, particle2, baseExceptionParams[i][0], baseExceptionParams[i][1], baseExceptionParams[i][2]);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
    }
    for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
        string param;
        int particle;
        double charge, sigma, epsilon;
        force.getParticleParameterOffset(i, param, particle, charge, sigma, epsilon);
        particleParamOffsets[make_pair(param, particle)] = {charge, sigma, epsilon};
    }
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionParamOffsets[make_pair(param, nb14Index[exception])] = {charge, sigma, epsilon};
    }
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
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
        ewaldAlpha = alpha;
    }
    else if (nonbondedMethod == PME) {
        double alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSize[0], gridSize[1], gridSize[2], false);
        ewaldAlpha = alpha;
    }
    else if (nonbondedMethod == LJPME) {
        double alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, gridSize[0], gridSize[1], gridSize[2], false);
        ewaldAlpha = alpha;
        NonbondedForceImpl::calcPMEParameters(system, force, alpha, dispersionGridSize[0], dispersionGridSize[1], dispersionGridSize[2], true);
        ewaldDispersionAlpha = alpha;
        useSwitchingFunction = false;
    }
    if (nonbondedMethod == NoCutoff || nonbondedMethod == CutoffNonPeriodic)
        exceptionsArePeriodic = false;
    else
        exceptionsArePeriodic = force.getExceptionsUsePeriodicBoundaryConditions();
    rfDielectric = force.getReactionFieldDielectric();
    if (force.getUseDispersionCorrection())
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(system, force);
    else
        dispersionCoefficient = 0.0;
}

double ReferenceCalcNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy, bool includeDirect, bool includeReciprocal) {
    computeParameters(context);
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceLJCoulombIxn clj;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    bool ljpme = (nonbondedMethod == LJPME);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, extractBoxVectors(context), periodic || ewald || pme || ljpme, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, rfDielectric);
    }
    if (periodic || ewald || pme || ljpme) {
        Vec3* boxVectors = extractBoxVectors(context);
        double minAllowedSize = 1.999999*nonbondedCutoff;
        if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
            throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
        clj.setPeriodic(boxVectors);
        clj.setPeriodicExceptions(exceptionsArePeriodic);
    }
    if (ewald)
        clj.setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        clj.setUsePME(ewaldAlpha, gridSize);
    if (ljpme){
        clj.setUsePME(ewaldAlpha, gridSize);
        clj.setUseLJPME(ewaldDispersionAlpha, dispersionGridSize);
    }
    if (useSwitchingFunction)
        clj.setUseSwitchingFunction(switchingDistance);
    clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusions, forceData, includeEnergy ? &energy : NULL, includeDirect, includeReciprocal);
    if (includeDirect) {
        ReferenceBondForce refBondForce;
        ReferenceLJCoulomb14 nonbonded14;
        if (exceptionsArePeriodic) {
            Vec3* boxVectors = extractBoxVectors(context);
            nonbonded14.setPeriodic(boxVectors);
        }
        refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, includeEnergy ? &energy : NULL, nonbonded14);
        if (periodic || ewald || pme) {
            Vec3* boxVectors = extractBoxVectors(context);
            energy += dispersionCoefficient/(boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2]);
        }
    }
    return energy;
}

void ReferenceCalcNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const NonbondedForce& force, int firstParticle, int lastParticle, int firstException, int lastException) {
    if (force.getNumParticles() != numParticles)
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Identify which exceptions are 1-4 interactions.

    set<int> exceptionsWithOffsets;
    for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
        string param;
        int exception;
        double charge, sigma, epsilon;
        force.getExceptionParameterOffset(i, param, exception, charge, sigma, epsilon);
        exceptionsWithOffsets.insert(exception);
    }
    vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon);
        if (chargeProd != 0.0 || epsilon != 0.0 || exceptionsWithOffsets.find(i) != exceptionsWithOffsets.end())
            nb14s.push_back(i);
    }
    if (nb14s.size() != num14)
        throw OpenMMException("updateParametersInContext: The number of non-excluded exceptions has changed");

    // Record the values.

    for (int i = firstParticle; i <= lastParticle; ++i)
        force.getParticleParameters(i, baseParticleParams[i][0], baseParticleParams[i][1], baseParticleParams[i][2]);
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        force.getExceptionParameters(nb14s[i], particle1, particle2, baseExceptionParams[i][0], baseExceptionParams[i][1], baseExceptionParams[i][2]);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
    }
    
    // Recompute the coefficient for the dispersion correction.

    NonbondedForce::NonbondedMethod method = force.getNonbondedMethod();
    if (force.getUseDispersionCorrection() && (method == NonbondedForce::CutoffPeriodic || method == NonbondedForce::Ewald || method == NonbondedForce::PME))
        dispersionCoefficient = NonbondedForceImpl::calcDispersionCorrection(context.getSystem(), force);
}

void ReferenceCalcNonbondedForceKernel::getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != PME && nonbondedMethod != LJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using PME or LJPME");
    alpha = ewaldAlpha;
    nx = gridSize[0];
    ny = gridSize[1];
    nz = gridSize[2];
}

void ReferenceCalcNonbondedForceKernel::getLJPMEParameters(double& alpha, int& nx, int& ny, int& nz) const {
    if (nonbondedMethod != LJPME)
        throw OpenMMException("getPMEParametersInContext: This Context is not using LJPME");
    alpha = ewaldDispersionAlpha;
    nx = dispersionGridSize[0];
    ny = dispersionGridSize[1];
    nz = dispersionGridSize[2];
}

void ReferenceCalcNonbondedForceKernel::computeParameters(ContextImpl& context) {
    // Compute particle parameters.

    vector<double> charges(numParticles), sigmas(numParticles), epsilons(numParticles);
    for (int i = 0; i < numParticles; i++) {
        charges[i] = baseParticleParams[i][0];
        sigmas[i] = baseParticleParams[i][1];
        epsilons[i] = baseParticleParams[i][2];
    }
    for (auto& offset : particleParamOffsets) {
        double value = context.getParameter(offset.first.first);
        int index = offset.first.second;
        charges[index] += value*offset.second[0];
        sigmas[index] += value*offset.second[1];
        epsilons[index] += value*offset.second[2];
    }
    for (int i = 0; i < numParticles; i++) {
        particleParamArray[i][0] = 0.5*sigmas[i];
        particleParamArray[i][1] = 2.0*sqrt(epsilons[i]);
        particleParamArray[i][2] = charges[i];
    }

    // Compute exception parameters.

    charges.resize(num14);
    sigmas.resize(num14);
    epsilons.resize(num14);
    for (int i = 0; i < num14; i++) {
        charges[i] = baseExceptionParams[i][0];
        sigmas[i] = baseExceptionParams[i][1];
        epsilons[i] = baseExceptionParams[i][2];
    }
    for (auto& offset : exceptionParamOffsets) {
        double value = context.getParameter(offset.first.first);
        int index = offset.first.second;
        charges[index] += value*offset.second[0];
        sigmas[index] += value*offset.second[1];
        epsilons[index] += value*offset.second[2];
    }
    for (int i = 0; i < num14; i++) {
        bonded14ParamArray[i][0] = sigmas[i];
        bonded14ParamArray[i][1] = 4.0*epsilons[i];
        bonded14ParamArray[i][2] = charges[i];
    }
}

ReferenceCalcCustomNonbondedForceKernel::~ReferenceCalcCustomNonbondedForceKernel() {
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

    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        force.getParticleParameters(i, particleParamArray[i]);
    nonbondedMethod = CalcCustomNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff) {
        neighborList = NULL;
        useSwitchingFunction = false;
    }
    else {
        neighborList = new NeighborList();
        useSwitchingFunction = force.getUseSwitchingFunction();
        switchingDistance = force.getSwitchingDistance();
    }

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the expressions.
    
    createExpressions(force);
    
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

void ReferenceCalcCustomNonbondedForceKernel::createExpressions(const CustomNonbondedForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the various expressions used to calculate the force.

    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(), functions).optimize();
    energyExpression = expression.createCompiledExpression();
    forceExpression = expression.differentiate("r").createCompiledExpression();
    parameterNames.clear();
    globalParameterNames.clear();
    globalParamValues.clear();
    computedValueNames.clear();
    computedValueExpressions.clear();
    energyParamDerivNames.clear();
    energyParamDerivExpressions.clear();
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        parameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParameterNames.push_back(force.getGlobalParameterName(i));
        globalParamValues[force.getGlobalParameterName(i)] = force.getGlobalParameterDefaultValue(i);
    }
    set<string> particleVariables, pairVariables;
    particleVariables.insert("r");
    pairVariables.insert("r");
    for (auto& name : parameterNames) {
        particleVariables.insert(name);
        pairVariables.insert(name+"1");
        pairVariables.insert(name+"2");
    }
    particleVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    pairVariables.insert(globalParameterNames.begin(), globalParameterNames.end());
    for (int i = 0; i < force.getNumComputedValues(); i++) {
        string name, exp;
        force.getComputedValueParameters(i, name, exp);
        Lepton::ParsedExpression parsed = Lepton::Parser::parse(exp, functions);
        validateVariables(parsed.getRootNode(), particleVariables);
        computedValueNames.push_back(name);
        computedValueExpressions.push_back(parsed.createCompiledExpression());
    }
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(expression.differentiate(param).createCompiledExpression());
    }
    for (auto& name : computedValueNames) {
        pairVariables.insert(name+"1");
        pairVariables.insert(name+"2");
    }
    validateVariables(expression.getRootNode(), pairVariables);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

double ReferenceCalcCustomNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    Vec3* boxVectors = extractBoxVectors(context);
    double energy = 0;
    ReferenceCustomNonbondedIxn ixn(energyExpression, forceExpression, parameterNames, energyParamDerivExpressions, computedValueNames, computedValueExpressions);
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
    for (auto& name : globalParameterNames) {
        double value = context.getParameter(name);
        if (globalParamValues[name] != value)
            globalParamsChanged = true;
        globalParamValues[name] = value;
    }
    if (useSwitchingFunction)
        ixn.setUseSwitchingFunction(switchingDistance);
    vector<double> energyParamDerivValues(energyParamDerivNames.size()+1, 0.0);
    ixn.calculatePairIxn(numParticles, posData, particleParamArray, exclusions, globalParamValues, forceData, includeEnergy ? &energy : NULL, &energyParamDerivValues[0]);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    for (int i = 0; i < energyParamDerivNames.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += energyParamDerivValues[i];
    
    // Add in the long range correction.
    
    if (!hasInitializedLongRangeCorrection) {
        ThreadPool threads;
        longRangeCorrectionData = CustomNonbondedForceImpl::prepareLongRangeCorrection(*forceCopy, threads.getNumThreads());
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, threads);
        hasInitializedLongRangeCorrection = true;
    }
    else if (globalParamsChanged && forceCopy != NULL) {
        ThreadPool threads;
        CustomNonbondedForceImpl::calcLongRangeCorrection(*forceCopy, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, threads);
    }
    double volume = boxVectors[0][0]*boxVectors[1][1]*boxVectors[2][2];
    energy += longRangeCoefficient/volume;
    for (int i = 0; i < longRangeCoefficientDerivs.size(); i++)
        energyParamDerivs[energyParamDerivNames[i]] += longRangeCoefficientDerivs[i]/volume;
    return energy;
}

void ReferenceCalcCustomNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomNonbondedForce& force, int firstParticle, int lastParticle) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> parameters;
    for (int i = firstParticle; i <= lastParticle; ++i) {
        force.getParticleParameters(i, parameters);
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = parameters[j];
    }
    
    // If necessary, recompute the long range correction.
    
    if (forceCopy != NULL) {
        ThreadPool threads;
        longRangeCorrectionData = CustomNonbondedForceImpl::prepareLongRangeCorrection(force, threads.getNumThreads());
        CustomNonbondedForceImpl::calcLongRangeCorrection(force, longRangeCorrectionData, context.getOwner(), longRangeCoefficient, longRangeCoefficientDerivs, threads);
        hasInitializedLongRangeCorrection = true;
        *forceCopy = force;
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed)
        createExpressions(force);
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
    vector<double> atomicRadii(numParticles);
    vector<double> scaleFactors(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        charges[i] = charge;
        atomicRadii[i] = radius;
        scaleFactors[i] = scalingFactor;
    }
    ObcParameters* obcParameters = new ObcParameters(numParticles, ObcParameters::ObcTypeII);
    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
    obcParameters->setSolventDielectric(force.getSolventDielectric());
    obcParameters->setSoluteDielectric(force.getSoluteDielectric());
    obcParameters->setPi4Asolv(4*M_PI*force.getSurfaceAreaEnergy());
    if (force.getNonbondedMethod() != GBSAOBCForce::NoCutoff)
        obcParameters->setUseCutoff(force.getCutoffDistance());
    isPeriodic = (force.getNonbondedMethod() == GBSAOBCForce::CutoffPeriodic);
    obc = new ReferenceObc(obcParameters);
    obc->setIncludeAceApproximation(true);
}

double ReferenceCalcGBSAOBCForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
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

    vector<double> atomicRadii(numParticles);
    vector<double> scaleFactors(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, scalingFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor);
        charges[i] = charge;
        atomicRadii[i] = radius;
        scaleFactors[i] = scalingFactor;
    }
    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
}

ReferenceCalcCustomGBForceKernel::~ReferenceCalcCustomGBForceKernel() {
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

    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        force.getParticleParameters(i, particleParamArray[i]);
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        particleParameterNames.push_back(force.getPerParticleParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    nonbondedMethod = CalcCustomGBForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = force.getCutoffDistance();
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the expressions.
    
    createExpressions(force);
}

void ReferenceCalcCustomGBForceKernel::createExpressions(const CustomGBForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expressions for computed values.

    valueExpressions.clear();
    valueTypes.clear();
    valueNames.clear();
    energyParamDerivNames.clear();
    valueDerivExpressions.clear();
    valueGradientExpressions.clear();
    valueParamDerivExpressions.clear();
    valueDerivExpressions.resize(force.getNumComputedValues());
    valueGradientExpressions.resize(force.getNumComputedValues());
    valueParamDerivExpressions.resize(force.getNumComputedValues());
    set<string> particleVariables, pairVariables;
    pairVariables.insert("r");
    particleVariables.insert("x");
    particleVariables.insert("y");
    particleVariables.insert("z");
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
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

    energyExpressions.clear();
    energyTypes.clear();
    energyDerivExpressions.clear();
    energyGradientExpressions.clear();
    energyParamDerivExpressions.clear();
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

    for (auto& function : functions)
        delete function.second;
}

double ReferenceCalcCustomGBForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    ReferenceCustomGBIxn ixn(valueExpressions, valueDerivExpressions, valueGradientExpressions, valueParamDerivExpressions, valueNames, valueTypes,
        energyExpressions, energyDerivExpressions, energyGradientExpressions, energyParamDerivExpressions, energyTypes, particleParameterNames);
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    if (periodic)
        ixn.setPeriodic(extractBoxVectors(context));
    if (nonbondedMethod != NoCutoff) {
        vector<set<int> > empty(context.getSystem().getNumParticles()); // Don't omit exclusions from the neighbor list
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, empty, extractBoxVectors(context), periodic, nonbondedCutoff, 0.0);
        ixn.setUseCutoff(nonbondedCutoff, *neighborList);
    }
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
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
            particleParamArray[i][j] = parameters[j];
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed)
        createExpressions(force);
}

ReferenceCalcCustomExternalForceKernel::~ReferenceCalcCustomExternalForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomExternalForceKernel::initialize(const System& system, const CustomExternalForce& force) {
    numParticles = force.getNumParticles();
    int numParameters = force.getNumPerParticleParameters();

    // Build the arrays.

    particles.resize(numParticles);
    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        force.getParticleParameters(i, particles[i], particleParamArray[i]);

    // Parse the expression used to calculate the force.

    map<string, Lepton::CustomFunction*> functions;
    ReferencePointDistanceFunction periodicDistance(true, &boxVectors);
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
    ixn = new ReferenceCustomExternalIxn(energyExpression, forceExpressionX, forceExpressionY, forceExpressionZ, parameterNames);

}

double ReferenceCalcCustomExternalForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    boxVectors = extractBoxVectors(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    ixn->setGlobalParameters(globalParameters);
    for (int i = 0; i < numParticles; ++i)
        ixn->calculateForce(particles[i], posData, particleParamArray[i], forceData, includeEnergy ? &energy : NULL);
    return energy;
}

void ReferenceCalcCustomExternalForceKernel::copyParametersToContext(ContextImpl& context, const CustomExternalForce& force, int firstParticle, int lastParticle) {
    if (numParticles != force.getNumParticles())
        throw OpenMMException("updateParametersInContext: The number of particles has changed");

    // Record the values.

    int numParameters = force.getNumPerParticleParameters();
    vector<double> params;
    for (int i = firstParticle; i <= lastParticle; ++i) {
        int particle;
        force.getParticleParameters(i, particle, params);
        if (particle != particles[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
        for (int j = 0; j < numParameters; j++)
            particleParamArray[i][j] = params[j];
    }
}

ReferenceCalcCustomHbondForceKernel::~ReferenceCalcCustomHbondForceKernel() {
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

    donorParticles.resize(numDonors);
    donorParamArray.resize(numDonors);
    for (int i = 0; i < numDonors; ++i) {
        int d1, d2, d3;
        force.getDonorParameters(i, d1, d2, d3, donorParamArray[i]);
        donorParticles[i].push_back(d1);
        donorParticles[i].push_back(d2);
        donorParticles[i].push_back(d3);
    }
    acceptorParticles.resize(numAcceptors);
    acceptorParamArray.resize(numAcceptors);
    for (int i = 0; i < numAcceptors; ++i) {
        int a1, a2, a3;
        force.getAcceptorParameters(i, a1, a2, a3, acceptorParamArray[i]);
        acceptorParticles[i].push_back(a1);
        acceptorParticles[i].push_back(a2);
        acceptorParticles[i].push_back(a3);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    nonbondedCutoff = force.getCutoffDistance();

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the interaction.
    
    createInteraction(force);
}

void ReferenceCalcCustomHbondForceKernel::createInteraction(const CustomHbondForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Parse the expression and create the object used to calculate the interaction.

    map<string, vector<int> > distances;
    map<string, vector<int> > angles;
    map<string, vector<int> > dihedrals;
    Lepton::ParsedExpression energyExpression = CustomHbondForceImpl::prepareExpression(force, functions, distances, angles, dihedrals);
    vector<string> donorParameterNames;
    vector<string> acceptorParameterNames;
    for (int i = 0; i < force.getNumPerDonorParameters(); i++)
        donorParameterNames.push_back(force.getPerDonorParameterName(i));
    for (int i = 0; i < force.getNumPerAcceptorParameters(); i++)
        acceptorParameterNames.push_back(force.getPerAcceptorParameterName(i));
    ixn = new ReferenceCustomHbondIxn(donorParticles, acceptorParticles, energyExpression, donorParameterNames, acceptorParameterNames, distances, angles, dihedrals);
    NonbondedMethod nonbondedMethod = CalcCustomHbondForceKernel::NonbondedMethod(force.getNonbondedMethod());
    isPeriodic = (nonbondedMethod == CutoffPeriodic);
    if (nonbondedMethod != NoCutoff)
        ixn->setUseCutoff(nonbondedCutoff);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

double ReferenceCalcCustomHbondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    if (isPeriodic)
        ixn->setPeriodic(extractBoxVectors(context));
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
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
            donorParamArray[i][j] = parameters[j];
    }
    int numAcceptorParameters = force.getNumPerAcceptorParameters();
    const vector<vector<int> >& acceptorAtoms = ixn->getAcceptorAtoms();
    for (int i = 0; i < numAcceptors; ++i) {
        int a1, a2, a3;
        force.getAcceptorParameters(i, a1, a2, a3, parameters);
        if (a1 != acceptorAtoms[i][0] || a2 != acceptorAtoms[i][1] || a3 != acceptorAtoms[i][2])
            throw OpenMMException("updateParametersInContext: The set of particles in an acceptor group has changed");
        for (int j = 0; j < numAcceptorParameters; j++)
            acceptorParamArray[i][j] = parameters[j];
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete ixn;
        ixn = NULL;
        createInteraction(force);
    }
}

ReferenceCalcCustomCentroidBondForceKernel::~ReferenceCalcCustomCentroidBondForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomCentroidBondForceKernel::initialize(const System& system, const CustomCentroidBondForce& force) {
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    int numGroups = force.getNumGroups();
    groupAtoms.resize(numGroups);
    vector<double> ignored;
    for (int i = 0; i < numGroups; i++)
        force.getGroupParameters(i, groupAtoms[i], ignored);
    CustomCentroidBondForceImpl::computeNormalizedWeights(force, system, normalizedWeights);
    numBonds = force.getNumBonds();
    bondGroups.resize(numBonds);
    bondParamArray.resize(numBonds);
    for (int i = 0; i < numBonds; ++i)
        force.getBondParameters(i, bondGroups[i], bondParamArray[i]);

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the interaction.
    
    createInteraction(force);
}

void ReferenceCalcCustomCentroidBondForceKernel::createInteraction(const CustomCentroidBondForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Create implementations of point functions.

    functions["pointdistance"] = new ReferencePointDistanceFunction(usePeriodic, &boxVectors);
    functions["pointangle"] = new ReferencePointAngleFunction(usePeriodic, &boxVectors);
    functions["pointdihedral"] = new ReferencePointDihedralFunction(usePeriodic, &boxVectors);

    // Parse the expression and create the object used to calculate the interaction.

    int numGroups = force.getNumGroups();
    Lepton::ParsedExpression energyExpression = CustomCentroidBondForceImpl::prepareExpression(force, functions);
    vector<string> bondParameterNames;
    for (int i = 0; i < force.getNumPerBondParameters(); i++)
        bondParameterNames.push_back(force.getPerBondParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(energyExpression.differentiate(param).createCompiledExpression());
    }
    ixn = new ReferenceCustomCentroidBondIxn(force.getNumGroupsPerBond(), groupAtoms, normalizedWeights, bondGroups, energyExpression, bondParameterNames, energyParamDerivExpressions);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

double ReferenceCalcCustomCentroidBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    if (usePeriodic) {
        boxVectors = extractBoxVectors(context);
        ixn->setPeriodic(boxVectors);
    }
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
            bondParamArray[i][j] = params[j];
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete ixn;
        ixn = NULL;
        createInteraction(force);
    }
}

ReferenceCalcCustomCompoundBondForceKernel::~ReferenceCalcCustomCompoundBondForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomCompoundBondForceKernel::initialize(const System& system, const CustomCompoundBondForce& force) {
    usePeriodic = force.usesPeriodicBoundaryConditions();

    // Build the arrays.

    numBonds = force.getNumBonds();
    bondParticles.resize(numBonds);
    bondParamArray.resize(numBonds);
    for (int i = 0; i < numBonds; ++i)
        force.getBondParameters(i, bondParticles[i], bondParamArray[i]);

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the interaction.

    createInteraction(force);
}

void ReferenceCalcCustomCompoundBondForceKernel::createInteraction(const CustomCompoundBondForce& force) {
    // Create custom functions for the tabulated functions.

    map<string, Lepton::CustomFunction*> functions;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

    // Create implementations of point functions.

    functions["pointdistance"] = new ReferencePointDistanceFunction(usePeriodic, &boxVectors);
    functions["pointangle"] = new ReferencePointAngleFunction(usePeriodic, &boxVectors);
    functions["pointdihedral"] = new ReferencePointDihedralFunction(usePeriodic, &boxVectors);

    // Parse the expression and create the object used to calculate the interaction.

    Lepton::ParsedExpression energyExpression = CustomCompoundBondForceImpl::prepareExpression(force, functions);
    vector<string> bondParameterNames;
    for (int i = 0; i < force.getNumPerBondParameters(); i++)
        bondParameterNames.push_back(force.getPerBondParameterName(i));
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        string param = force.getEnergyParameterDerivativeName(i);
        energyParamDerivNames.push_back(param);
        energyParamDerivExpressions.push_back(energyExpression.differentiate(param).createCompiledExpression());
    }
    ixn = new ReferenceCustomCompoundBondIxn(force.getNumParticlesPerBond(), bondParticles, energyExpression, bondParameterNames, energyParamDerivExpressions);

    // Delete the custom functions.

    for (auto& function : functions)
        delete function.second;
}

double ReferenceCalcCustomCompoundBondForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    if (usePeriodic) {
        boxVectors = extractBoxVectors(context);
        ixn->setPeriodic(boxVectors);
    }
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
            bondParamArray[i][j] = params[j];
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete ixn;
        ixn = NULL;
        createInteraction(force);
    }
}

ReferenceCalcCustomManyParticleForceKernel::~ReferenceCalcCustomManyParticleForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomManyParticleForceKernel::initialize(const System& system, const CustomManyParticleForce& force) {
    // Build the arrays.

    numParticles = system.getNumParticles();
    particleParamArray.resize(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        int type;
        force.getParticleParameters(i, particleParamArray[i], type);
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));

    // Record the tabulated function update counts for future reference.

    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        tabulatedFunctionUpdateCount[force.getTabulatedFunctionName(i)] = force.getTabulatedFunction(i).getUpdateCount();

    // Create the interaction.

    ixn = new ReferenceCustomManyParticleIxn(force);
    nonbondedMethod = CalcCustomManyParticleForceKernel::NonbondedMethod(force.getNonbondedMethod());
    cutoffDistance = force.getCutoffDistance();
}

double ReferenceCalcCustomManyParticleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    if (nonbondedMethod == CutoffPeriodic) {
        Vec3* boxVectors = extractBoxVectors(context);
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
            particleParamArray[i][j] = parameters[j];
    }

    // See if any tabulated functions have changed.

    bool changed = false;
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        string name = force.getTabulatedFunctionName(i);
        if (force.getTabulatedFunction(i).getUpdateCount() != tabulatedFunctionUpdateCount[name]) {
            tabulatedFunctionUpdateCount[name] = force.getTabulatedFunction(i).getUpdateCount();
            changed = true;
        }
    }
    if (changed) {
        delete ixn;
        ixn = NULL;
        ixn = new ReferenceCustomManyParticleIxn(force);
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

ReferenceCalcCustomCVForceKernel::~ReferenceCalcCustomCVForceKernel() {
    if (ixn != NULL)
        delete ixn;
}

void ReferenceCalcCustomCVForceKernel::initialize(const System& system, const CustomCVForce& force, ContextImpl& innerContext) {
    for (int i = 0; i < force.getNumGlobalParameters(); i++)
        globalParameterNames.push_back(force.getGlobalParameterName(i));
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        energyParamDerivNames.push_back(force.getEnergyParameterDerivativeName(i));
    ixn = new ReferenceCustomCVForce(force);
}

double ReferenceCalcCustomCVForceKernel::execute(ContextImpl& context, ContextImpl& innerContext, bool includeForces, bool includeEnergy) {
    copyState(context, innerContext);
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = 0;
    map<string, double> globalParameters;
    for (auto& name : globalParameterNames)
        globalParameters[name] = context.getParameter(name);
    map<string, double>& energyParamDerivs = extractEnergyParameterDerivatives(context);
    ixn->calculateIxn(innerContext, posData, globalParameters, forceData, includeEnergy ? &energy : NULL, energyParamDerivs);
    return energy;
}

void ReferenceCalcCustomCVForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext) {
    extractPositions(innerContext) = extractPositions(context);
    extractVelocities(innerContext) = extractVelocities(context);
    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext.setPeriodicBoxVectors(a, b, c);
    innerContext.setTime(context.getTime());
    map<string, double> innerParameters = innerContext.getParameters();
    for (auto& param : innerParameters)
        innerContext.setParameter(param.first, context.getParameter(param.first));
}

void ReferenceCalcCustomCVForceKernel::copyParametersToContext(ContextImpl& context, const CustomCVForce& force) {
    ixn->updateTabulatedFunctions(force);
}

void ReferenceCalcRMSDForceKernel::initialize(const System& system, const RMSDForce& force) {
    particles = force.getParticles();
    if (particles.size() == 0)
        for (int i = 0; i < system.getNumParticles(); i++)
            particles.push_back(i);
    referencePos = force.getReferencePositions();
    Vec3 center;
    for (int i : particles)
        center += referencePos[i];
    center /= particles.size();
    for (Vec3& p : referencePos)
        p -= center;
}

double ReferenceCalcRMSDForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    ReferenceRMSDForce rmsd(referencePos, particles);
    return rmsd.calculateIxn(posData, forceData);
}

void ReferenceCalcRMSDForceKernel::copyParametersToContext(ContextImpl& context, const RMSDForce& force) {
    if (referencePos.size() != force.getReferencePositions().size())
        throw OpenMMException("updateParametersInContext: The number of reference positions has changed");
    particles = force.getParticles();
    if (particles.size() == 0)
        for (int i = 0; i < referencePos.size(); i++)
            particles.push_back(i);
    referencePos = force.getReferencePositions();
    Vec3 center;
    for (int i : particles)
        center += referencePos[i];
    center /= particles.size();
    for (Vec3& p : referencePos)
        p -= center;
}

ReferenceIntegrateVerletStepKernel::~ReferenceIntegrateVerletStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateVerletStepKernel::initialize(const System& system, const VerletIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

void ReferenceIntegrateVerletStepKernel::execute(ContextImpl& context, const VerletIntegrator& integrator) {
    double stepSize = integrator.getStepSize();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceVerletDynamics(context.getSystem().getNumParticles(), stepSize);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateVerletStepKernel::computeKineticEnergy(ContextImpl& context, const VerletIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}

ReferenceIntegrateNoseHooverStepKernel::~ReferenceIntegrateNoseHooverStepKernel() {
    if (chainPropagator)
        delete chainPropagator;
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateNoseHooverStepKernel::initialize(const System& system, const NoseHooverIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    this->chainPropagator = new ReferenceNoseHooverChain();
}

void ReferenceIntegrateNoseHooverStepKernel::execute(ContextImpl& context, const NoseHooverIntegrator& integrator, bool &forcesAreValid) {
    double stepSize = integrator.getStepSize();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceNoseHooverDynamics(context.getSystem().getNumParticles(), stepSize);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
        prevStepSize = stepSize;
    }
    dynamics->step1(context, context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance(), forcesAreValid,
                     integrator.getAllThermostatedIndividualParticles(), integrator.getAllThermostatedPairs(), integrator.getMaximumPairDistance());
    int numChains = integrator.getNumThermostats();
    for(int chain = 0; chain < numChains; ++chain) {
        const auto &thermostatChain = integrator.getThermostat(chain);
        std::pair<double, double> KEs = computeMaskedKineticEnergy(context, thermostatChain, true);
        std::pair<double, double> scaleFactors = propagateChain(context, thermostatChain, KEs, stepSize);
        scaleVelocities(context, thermostatChain, scaleFactors);
    }
    dynamics->step2(context, context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance(), forcesAreValid,
                     integrator.getAllThermostatedIndividualParticles(), integrator.getAllThermostatedPairs(), integrator.getMaximumPairDistance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateNoseHooverStepKernel::computeKineticEnergy(ContextImpl& context, const NoseHooverIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0);
}

std::pair<double, double> ReferenceIntegrateNoseHooverStepKernel::propagateChain(ContextImpl& context, const NoseHooverChain &nhc,
                                                                                     std::pair<double, double> kineticEnergy, double timeStep) {
    double absKE = kineticEnergy.first;
    double relKE = kineticEnergy.second;
    if (absKE < 1e-8) return {1.0, 1.0};  // (catches the problem of zero velocities in the first dynamics step, where we have nothing to scale)
    // Get the variables describing the NHC
    int chainLength = nhc.getChainLength();
    int chainID = nhc.getChainID();
    int numDOFs = nhc.getNumDegreesOfFreedom();
    int numMTS = nhc.getNumMultiTimeSteps();

    int nAtoms = nhc.getThermostatedAtoms().size();
    double absScale = 0;
    if (nAtoms) {
        if (chainPositions.size() < 2*chainID+1){
            chainPositions.resize(2*chainID+1);
        }
        if (chainVelocities.size() < 2*chainID+1){
            chainVelocities.resize(2*chainID+1);
        }
        auto& positions = chainPositions.at(2*chainID);
        auto& velocities = chainVelocities.at(2*chainID);
        if (positions.size() < chainLength){
            positions.resize(chainLength, 0);
        }
        if (velocities.size() < chainLength){
            velocities.resize(chainLength, 0);
        }
        double temperature = nhc.getTemperature();
        double collisionFrequency = nhc.getCollisionFrequency();
        absScale = chainPropagator->propagate(absKE, velocities, positions, numDOFs,
                                              temperature, collisionFrequency, timeStep,
                                              numMTS, nhc.getYoshidaSuzukiWeights());
    }
    double relScale = 0;
    int nPairs = nhc.getThermostatedPairs().size();
    if (nPairs) {
        if (chainPositions.size() < 2*chainID+2){
            chainPositions.resize(2*chainID+2);
        }
        if (chainVelocities.size() < 2*chainID+2){
            chainVelocities.resize(2*chainID+2);
        }
        auto& positions = chainPositions.at(2*chainID+1);
        auto& velocities = chainVelocities.at(2*chainID+1);
        if (positions.size() < chainLength){
            positions.resize(chainLength, 0);
        }
        if (velocities.size() < chainLength){
            velocities.resize(chainLength, 0);
        }
        double temperature = nhc.getRelativeTemperature();
        double collisionFrequency = nhc.getRelativeCollisionFrequency();
        relScale = chainPropagator->propagate(relKE, velocities, positions, 3*nPairs,
                                              temperature, collisionFrequency, timeStep,
                                              numMTS, nhc.getYoshidaSuzukiWeights());
    }
    return {absScale, relScale};
}

double ReferenceIntegrateNoseHooverStepKernel::computeHeatBathEnergy(ContextImpl& context, const NoseHooverChain &nhc) {
    double potentialEnergy = 0;
    double kineticEnergy = 0;
    int chainLength = nhc.getChainLength();
    int chainID = nhc.getChainID();
    int nAtoms = nhc.getThermostatedAtoms().size();
    int nPairs = nhc.getThermostatedPairs().size();
    if (nAtoms) {
        double temperature = nhc.getTemperature();
        double collisionFrequency = nhc.getCollisionFrequency();
        double kT = temperature * BOLTZ;
        int numDOFs = nhc.getNumDegreesOfFreedom();
        for(int i = 0; i < chainLength; ++i) {
            double prefac = i ? 1 : numDOFs;
            double mass = prefac * kT / (collisionFrequency * collisionFrequency);
            double velocity = chainVelocities[2*chainID][i];
            // The kinetic energy of this bead
            kineticEnergy += 0.5 * mass * velocity * velocity;
            // The potential energy of this bead
            double position = chainPositions[2*chainID][i];
            potentialEnergy += prefac * kT * position;
        }
    }
    if (nPairs) {
        double temperature = nhc.getRelativeTemperature();
        double collisionFrequency = nhc.getRelativeCollisionFrequency();
        double kT = temperature * BOLTZ;
        int numDOFs = 3 * nPairs;
        for(int i = 0; i < chainLength; ++i) {
            double prefac = i ? 1 : numDOFs;
            double mass = prefac * kT / (collisionFrequency * collisionFrequency);
            double velocity = chainVelocities[2*chainID+1][i];
            // The kinetic energy of this bead
            kineticEnergy += 0.5 * mass * velocity * velocity;
            // The potential energy of this bead
            double position = chainPositions[2*chainID+1][i];
            potentialEnergy += prefac * kT * position;
        }
    }
    return kineticEnergy + potentialEnergy;
}

std::pair<double, double> ReferenceIntegrateNoseHooverStepKernel::computeMaskedKineticEnergy(ContextImpl& context,
                                                                const NoseHooverChain &noseHooverChain, bool downloadValue) {
    const std::vector<int>& atomsList = noseHooverChain.getThermostatedAtoms();
    const std::vector<std::pair<int,int>>& pairsList = noseHooverChain.getThermostatedPairs();
    std::vector<Vec3>& velocities = extractVelocities(context);
    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    std::vector<double> masses(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);

    double comKE = 0;
    double relKE = 0;
    // kinetic energy of individual atoms
    for (const auto &m: atomsList){
        comKE += 0.5 * masses[m] * velocities[m].dot(velocities[m]);
    }
    // center of mass kinetic energy of pairs
    for (const auto &p: pairsList){
        double m1 = masses[p.first];
        double m2 = masses[p.second];
        Vec3 v1 = velocities[p.first];
        Vec3 v2 = velocities[p.second];
        double invMass = 1.0 / (m1 + m2);
        double redMass = m1 * m2 * invMass;
        double fracM1 = m1 * invMass;
        double fracM2 = m2 * invMass;
        Vec3 comVelocity = fracM1 * v1 + fracM2 * v2;
        Vec3 relVelocity = v2 - v1;

        comKE += 0.5 * (m1 + m2) * comVelocity.dot(comVelocity);
        relKE += 0.5 * redMass * relVelocity.dot(relVelocity);
    }
    // We ignore the downloadValue argument here and always return the correct value
    return {comKE, relKE};
}


void ReferenceIntegrateNoseHooverStepKernel::scaleVelocities(ContextImpl& context, const NoseHooverChain &noseHooverChain, std::pair<double, double> scaleFactors) {
    const auto& atoms = noseHooverChain.getThermostatedAtoms();
    const auto& pairs = noseHooverChain.getThermostatedPairs();
    std::vector<Vec3>& velocities = extractVelocities(context);
    double absScale = scaleFactors.first;
    double relScale = scaleFactors.second;

    const System& system = context.getSystem();
    int numParticles = system.getNumParticles();
    std::vector<double> masses(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    // scale absolute velocities
    for (const auto &a: atoms){
        velocities[a] *= absScale;
    }
    // scale relative velocities and absolute center of mass velocities for each pair
    for (const auto &p: pairs){
        int p1 = p.first;
        int p2 = p.second;
        double m1 = masses[p.first];
        double m2 = masses[p.second];
        Vec3 v1 = velocities[p.first];
        Vec3 v2 = velocities[p.second];
        double invMass = 1.0 / (m1 + m2);
        double fracM1 = m1 * invMass;
        double fracM2 = m2 * invMass;
        Vec3 comVelocity = fracM1 * v1 + fracM2 * v2;
        Vec3 relVelocity = v2 - v1;
        velocities[p1] = absScale * comVelocity - relScale * relVelocity * fracM2;
        velocities[p2] = absScale * comVelocity + relScale * relVelocity * fracM1;
    }
}

void ReferenceIntegrateNoseHooverStepKernel::createCheckpoint(ContextImpl& context, ostream& stream) const {
    size_t numChains = chainPositions.size();
    assert(numChains == chainVelocities.size());
    stream.write((char*) &numChains, sizeof(size_t));
    for (size_t i=0; i<numChains; i++){
        auto & noseHooverPositions = chainPositions.at(i);
        auto & noseHooverVelocities = chainVelocities.at(i);
        size_t numBeads = noseHooverPositions.size();
        assert(numBeads == noseHooverVelocities.size());
        stream.write((char*) &numBeads, sizeof(size_t));
        stream.write((char*) noseHooverPositions.data(), sizeof(double)*numBeads);
        stream.write((char*) noseHooverVelocities.data(), sizeof(double)*numBeads);
    }
}

void ReferenceIntegrateNoseHooverStepKernel::loadCheckpoint(ContextImpl& context, istream& stream) {
    size_t numChains, numBeads;
    stream.read((char*) &numChains, sizeof(size_t));
    chainPositions.clear();
    chainVelocities.clear();
    for (size_t i=0; i<numChains; i++){
        stream.read((char*) &numBeads, sizeof(size_t));
        std::vector<double> noseHooverPositions(numBeads);
        std::vector<double> noseHooverVelocities(numBeads);
        stream.read((char*) &noseHooverPositions[0], sizeof(double)*numBeads);
        stream.read((char*) &noseHooverVelocities[0], sizeof(double)*numBeads);
        chainPositions.push_back(noseHooverPositions);
        chainVelocities.push_back(noseHooverVelocities);
    }
}

void ReferenceIntegrateNoseHooverStepKernel::getChainStates(ContextImpl& context, vector<vector<double> >& positions, vector<vector<double> >& velocities) const {
    positions = chainPositions;
    velocities = chainVelocities;
}

void ReferenceIntegrateNoseHooverStepKernel::setChainStates(ContextImpl& context, const vector<vector<double> >& positions, const vector<vector<double> >& velocities) {
    chainPositions = positions;
    chainVelocities = velocities;
}

ReferenceIntegrateLangevinMiddleStepKernel::~ReferenceIntegrateLangevinMiddleStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateLangevinMiddleStepKernel::initialize(const System& system, const LangevinMiddleIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateLangevinMiddleStepKernel::execute(ContextImpl& context, const LangevinMiddleIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceLangevinMiddleDynamics(
                context.getSystem().getNumParticles(), 
                stepSize, 
                friction, 
                temperature);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context, posData, velData, masses, integrator.getConstraintTolerance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateLangevinMiddleStepKernel::computeKineticEnergy(ContextImpl& context, const LangevinMiddleIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.0);
}

ReferenceIntegrateBrownianStepKernel::~ReferenceIntegrateBrownianStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateBrownianStepKernel::initialize(const System& system, const BrownianIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateBrownianStepKernel::execute(ContextImpl& context, const BrownianIntegrator& integrator) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double stepSize = integrator.getStepSize();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceBrownianDynamics(
                context.getSystem().getNumParticles(), 
                stepSize, 
                friction, 
                temperature);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
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
        masses[i] = system.getParticleMass(i);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

double ReferenceIntegrateVariableLangevinStepKernel::execute(ContextImpl& context, const VariableLangevinIntegrator& integrator, double maxTime) {
    double temperature = integrator.getTemperature();
    double friction = integrator.getFriction();
    double errorTol = integrator.getErrorTolerance();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    if (dynamics == 0 || temperature != prevTemp || friction != prevFriction || errorTol != prevErrorTol) {
        // Recreate the computation objects with the new parameters.

        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceVariableStochasticDynamics(context.getSystem().getNumParticles(), friction, temperature, errorTol);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
        prevTemp = temperature;
        prevFriction = friction;
        prevErrorTol = errorTol;
    }
    double maxStepSize = maxTime-data.time;
    if (integrator.getMaximumStepSize() > 0)
        maxStepSize = min(integrator.getMaximumStepSize(), maxStepSize);
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
        masses[i] = system.getParticleMass(i);
}

double ReferenceIntegrateVariableVerletStepKernel::execute(ContextImpl& context, const VariableVerletIntegrator& integrator, double maxTime) {
    double errorTol = integrator.getErrorTolerance();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    if (dynamics == 0 || errorTol != prevErrorTol) {
        // Recreate the computation objects with the new parameters.

        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceVariableVerletDynamics(context.getSystem().getNumParticles(), errorTol);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        dynamics->setVirtualSites(extractVirtualSites(context));
        prevErrorTol = errorTol;
    }
    double maxStepSize = maxTime-data.time;
    if (integrator.getMaximumStepSize() > 0)
        maxStepSize = min(integrator.getMaximumStepSize(), maxStepSize);
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
        masses[i] = system.getParticleMass(i);
    perDofValues.resize(integrator.getNumPerDofVariables());
    for (auto& values : perDofValues)
        values.resize(numParticles);

    // Create the computation objects.

    dynamics = new ReferenceCustomDynamics(system.getNumParticles(), integrator);
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) integrator.getRandomNumberSeed());
}

void ReferenceIntegrateCustomStepKernel::execute(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    
    // Record global variables.
    
    map<string, double> globals;
    globals["dt"] = integrator.getStepSize();
    for (int i = 0; i < integrator.getNumGlobalVariables(); i++)
        globals[integrator.getGlobalVariableName(i)] = globalValues[i];
    
    // Execute the step.
    
    dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
    dynamics->setVirtualSites(extractVirtualSites(context));
    dynamics->update(context, context.getSystem().getNumParticles(), posData, velData, forceData, masses, globals, perDofValues, forcesAreValid, integrator.getConstraintTolerance());
    
    // Record changed global variables.
    
    integrator.setStepSize(globals["dt"]);
    for (int i = 0; i < (int) globalValues.size(); i++)
        globalValues[i] = globals[integrator.getGlobalVariableName(i)];
    data.time += dynamics->getDeltaT();
    data.stepCount++;
}

double ReferenceIntegrateCustomStepKernel::computeKineticEnergy(ContextImpl& context, CustomIntegrator& integrator, bool& forcesAreValid) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    
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
        masses[i] = system.getParticleMass(i);
    this->thermostat = new ReferenceAndersenThermostat();
    SimTKOpenMMUtilities::setRandomNumberSeed((unsigned int) thermostat.getRandomNumberSeed());
    particleGroups = AndersenThermostatImpl::calcParticleGroups(system);
}

void ReferenceApplyAndersenThermostatKernel::execute(ContextImpl& context) {
    vector<Vec3>& velData = extractVelocities(context);
    thermostat->applyThermostat(particleGroups, velData, masses,
        context.getParameter(AndersenThermostat::Temperature()),
        context.getParameter(AndersenThermostat::CollisionFrequency()),
        context.getIntegrator().getStepSize());
}

ReferenceApplyMonteCarloBarostatKernel::~ReferenceApplyMonteCarloBarostatKernel() {
    if (barostat)
        delete barostat;
}

void ReferenceApplyMonteCarloBarostatKernel::initialize(const System& system, const Force& barostat, bool rigidMolecules) {
    this->rigidMolecules = rigidMolecules;
}

void ReferenceApplyMonteCarloBarostatKernel::saveCoordinates(ContextImpl& context) {
    if (barostat == NULL) {
        if (rigidMolecules)
            barostat = new ReferenceMonteCarloBarostat(context.getSystem().getNumParticles(), context.getMolecules());
        else {
            vector<vector<int> > molecules(context.getSystem().getNumParticles());
            for (int i = 0; i < molecules.size(); i++)
                molecules[i].push_back(i);
            barostat = new ReferenceMonteCarloBarostat(context.getSystem().getNumParticles(), molecules);
        }
    }
    vector<Vec3>& posData = extractPositions(context);
    barostat->savePositions(posData);
}

void ReferenceApplyMonteCarloBarostatKernel::scaleCoordinates(ContextImpl& context, double scaleX, double scaleY, double scaleZ) {
    vector<Vec3>& posData = extractPositions(context);
    Vec3* boxVectors = extractBoxVectors(context);
    barostat->applyBarostat(posData, boxVectors, scaleX, scaleY, scaleZ);
}

void ReferenceApplyMonteCarloBarostatKernel::restoreCoordinates(ContextImpl& context) {
    vector<Vec3>& posData = extractPositions(context);
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
    vector<Vec3>& velData = extractVelocities(context);
    
    // Calculate the center of mass momentum.
    
    double momentum[] = {0.0, 0.0, 0.0};
    double mass = 0.0;
    for (size_t i = 0; i < masses.size(); ++i) {
        momentum[0] += masses[i]*velData[i][0];
        momentum[1] += masses[i]*velData[i][1];
        momentum[2] += masses[i]*velData[i][2];
        mass += masses[i];
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

void ReferenceCalcATMForceKernel::initialize(const System& system, const ATMForce& force) {
    numParticles = force.getNumParticles();

    //displacement map
    displ1.resize(numParticles);
    displ0.resize(numParticles);
    for (int i = 0; i < numParticles; i++) {
        Vec3 displacement1, displacement0;
        force.getParticleParameters(i, displacement1, displacement0 );
        displ1[i] = displacement1;
        displ0[i] = displacement0;
    }
}

void ReferenceCalcATMForceKernel::applyForces(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1,
        double dEdu0, double dEdu1, const map<string, double>& energyParamDerivs) {
    vector<Vec3>& force = extractForces(context);
    vector<Vec3>& force0 = extractForces(innerContext0);
    vector<Vec3>& force1 = extractForces(innerContext1);
    for (int i = 0; i < force.size(); i++)
        force[i] += dEdu0*force0[i] + dEdu1*force1[i];
    map<string, double>& derivs = extractEnergyParameterDerivatives(context);
    for (auto deriv : energyParamDerivs)
        derivs[deriv.first] += deriv.second;
}

void ReferenceCalcATMForceKernel::copyState(ContextImpl& context, ContextImpl& innerContext0, ContextImpl& innerContext1) {
    vector<Vec3>& pos = extractPositions(context);

    //in the initial state, particles are displaced by displ0
    vector<Vec3> pos0(pos);
    for (int i = 0; i < pos0.size(); i++)
        pos0[i] += displ0[i];
    innerContext0.setPositions(pos0);

    //in the target state, particles are displaced by displ1
    vector<Vec3> pos1(pos);
    for (int i = 0; i < pos1.size(); i++)
        pos1[i] += displ1[i];
    innerContext1.setPositions(pos1);

    Vec3 a, b, c;
    context.getPeriodicBoxVectors(a, b, c);
    innerContext0.setPeriodicBoxVectors(a, b, c);
    innerContext1.setPeriodicBoxVectors(a, b, c);

    innerContext0.setTime(context.getTime());
    innerContext1.setTime(context.getTime());

    map<string, double> innerParameters;

    innerParameters = innerContext0.getParameters();
    for (auto& param : innerParameters)
        innerContext0.setParameter(param.first, context.getParameter(param.first));

    innerParameters = innerContext1.getParameters();
    for (auto& param : innerParameters)
        innerContext1.setParameter(param.first, context.getParameter(param.first));

}

void ReferenceCalcATMForceKernel::copyParametersToContext(ContextImpl& context, const ATMForce& force) {
    if (force.getNumParticles() != numParticles)
          throw OpenMMException("copyParametersToContext: The number of ATMForce particles has changed");
    displ1.resize(numParticles);
    displ0.resize(numParticles);
    for (int i = 0; i < numParticles; i++) {
        Vec3 displacement1, displacement0;
        force.getParticleParameters(i, displacement1, displacement0 );
        displ1[i] = displacement1;
        displ0[i] = displacement0;
    }
}

void ReferenceCalcCustomCPPForceKernel::initialize(const System& system, CustomCPPForceImpl& force) {
    this->force = &force;
    forces.resize(system.getNumParticles());
}

double ReferenceCalcCustomCPPForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& forceData = extractForces(context);
    double energy = force->computeForce(context, posData, forces);
    if (includeForces)
        for (int i = 0; i < forces.size(); i++)
            forceData[i] += forces[i];
    return energy;
}
