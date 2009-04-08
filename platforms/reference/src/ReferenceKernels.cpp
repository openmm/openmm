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
#include "ReferenceFloatStreamImpl.h"
#include "gbsa/CpuObc.h"
#include "SimTKReference/ReferenceAndersenThermostat.h"
#include "SimTKReference/ReferenceAngleBondIxn.h"
#include "SimTKReference/ReferenceBondForce.h"
#include "SimTKReference/ReferenceBrownianDynamics.h"
#include "SimTKReference/ReferenceHarmonicBondIxn.h"
#include "SimTKReference/ReferenceLincsAlgorithm.h"
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
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
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

void ReferenceInitializeForcesKernel::initialize(const System& system) {
}

void ReferenceInitializeForcesKernel::execute(OpenMMContextImpl& context) {
    double zero[] = {0.0, 0.0, 0.0};
    context.getForces().fillWithValue(zero);
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

void ReferenceCalcHarmonicBondForceKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    ReferenceBondForce refBondForce;
    ReferenceHarmonicBondIxn harmonicBond;
    refBondForce.calculateForce(numBonds, bondIndexArray, posData, bondParamArray, forceData, 0, 0, 0, harmonicBond);
}

double ReferenceCalcHarmonicBondForceKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
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

void ReferenceCalcHarmonicAngleForceKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    ReferenceBondForce refBondForce;
    ReferenceAngleBondIxn angleBond;
    refBondForce.calculateForce(numAngles, angleIndexArray, posData, angleParamArray, forceData, 0, 0, 0, angleBond);
}

double ReferenceCalcHarmonicAngleForceKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
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

void ReferenceCalcPeriodicTorsionForceKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    ReferenceBondForce refBondForce;
    ReferenceProperDihedralBond periodicTorsionBond;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, 0, 0, 0, periodicTorsionBond);
}

double ReferenceCalcPeriodicTorsionForceKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
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

void ReferenceCalcRBTorsionForceKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    ReferenceBondForce refBondForce;
    ReferenceRbDihedralBond rbTorsionBond;
    refBondForce.calculateForce(numTorsions, torsionIndexArray, posData, torsionParamArray, forceData, 0, 0, 0, rbTorsionBond);
}

double ReferenceCalcRBTorsionForceKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
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

ReferenceCalcNonbondedForceKernel::~ReferenceCalcNonbondedForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    disposeIntArray(exclusionArray, numParticles);
    disposeIntArray(bonded14IndexArray, num14);
    disposeRealArray(bonded14ParamArray, num14);
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceCalcNonbondedForceKernel::initialize(const System& system, const NonbondedForce& force, const std::vector<std::set<int> >& exclusions) {
    numParticles = force.getNumParticles();
    num14 = force.getNumNonbonded14();
    bonded14IndexArray = allocateIntArray(num14, 2);
    bonded14ParamArray = allocateRealArray(num14, 3);
    particleParamArray = allocateRealArray(numParticles, 3);
    RealOpenMM sqrtEps = static_cast<RealOpenMM>( std::sqrt(138.935485) );
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, depth;
        force.getParticleParameters(i, charge, radius, depth);
        particleParamArray[i][0] = static_cast<RealOpenMM>(0.5*radius);
        particleParamArray[i][1] = static_cast<RealOpenMM>(2.0*sqrt(depth));
        particleParamArray[i][2] = static_cast<RealOpenMM>(charge*sqrtEps);
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
        force.getNonbonded14Parameters(i, particle1, particle2, charge, radius, depth);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
        bonded14ParamArray[i][0] = static_cast<RealOpenMM>(radius);
        bonded14ParamArray[i][1] = static_cast<RealOpenMM>(4.0*depth);
        bonded14ParamArray[i][2] = static_cast<RealOpenMM>(charge*sqrtEps*sqrtEps);
    }
    nonbondedMethod = CalcNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
    nonbondedCutoff = (RealOpenMM) force.getCutoffDistance();
    Vec3 boxVectors[3];
    force.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
    periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
    periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();
        
}

void ReferenceCalcNonbondedForceKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    ReferenceLJCoulombIxn clj;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, (periodic || ewald) ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, 78.3f);
    }
    if (periodic||ewald)
        clj.setPeriodic(periodicBoxSize);
    if (ewald) {
        clj.setRecipVectors();
        clj.calculateEwaldIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, 0);
    }
    else {
        clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, 0);
    }
    ReferenceBondForce refBondForce;
    ReferenceLJCoulomb14 nonbonded14;
    if (nonbondedMethod != NoCutoff)
        nonbonded14.setUseCutoff(nonbondedCutoff, 78.3f);
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, 0, 0, 0, nonbonded14);
}

double ReferenceCalcNonbondedForceKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = allocateRealArray(numParticles, 3);
    RealOpenMM energy = 0;
    ReferenceLJCoulombIxn clj;
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, (periodic || ewald) ? periodicBoxSize : NULL, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, 78.3f);
    }
    if (periodic || ewald)
        clj.setPeriodic(periodicBoxSize);
    if (ewald) {
        clj.setRecipVectors();
        clj.calculateEwaldIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, &energy);
    }
    else {
        clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, &energy);
    }
    ReferenceBondForce refBondForce;
    ReferenceLJCoulomb14 nonbonded14;
    if (nonbondedMethod != NoCutoff)
        nonbonded14.setUseCutoff(nonbondedCutoff, 78.3f);
    RealOpenMM* energyArray = new RealOpenMM[num14];
    for (int i = 0; i < num14; ++i)
        energyArray[i] = 0;
    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, energyArray, 0, &energy, nonbonded14);
    disposeRealArray(forceData, numParticles);
    delete[] energyArray;
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

    // If there is a NonbondedForce in this system, use it to initialize cutoffs and periodic boundary conditions.

    for (int i = 0; i < system.getNumForces(); i++) {
        const NonbondedForce* nonbonded = dynamic_cast<const NonbondedForce*>(&system.getForce(i));
        if (nonbonded != NULL) {
            if (nonbonded->getNonbondedMethod() != NonbondedForce::NoCutoff)
                obcParameters->setUseCutoff(nonbonded->getCutoffDistance());
            if (nonbonded->getNonbondedMethod() == NonbondedForce::CutoffPeriodic) {
                Vec3 boxVectors[3];
                nonbonded->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
                RealOpenMM periodicBoxSize[3];
                periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
                periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
                periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
                obcParameters->setPeriodic(periodicBoxSize);
            }
            break;
        }
    }
    obc = new CpuObc(obcParameters);
    obc->setIncludeAceApproximation(true);
}

void ReferenceCalcGBSAOBCForceKernel::executeForces(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = ((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData();
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
}

double ReferenceCalcGBSAOBCForceKernel::executeEnergy(OpenMMContextImpl& context) {
    RealOpenMM** posData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData()); // Reference code needs to be made const correct
    RealOpenMM** forceData = allocateRealArray(context.getSystem().getNumParticles(), 3);
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
    disposeRealArray(forceData, context.getSystem().getNumParticles());
    return obc->getEnergy();
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

void ReferenceIntegrateVerletStepKernel::execute(OpenMMContextImpl& context, const VerletIntegrator& integrator) {
    double stepSize = integrator.getStepSize();
    RealOpenMM** posData = ((ReferenceFloatStreamImpl&) context.getPositions().getImpl()).getData();
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    RealOpenMM** forceData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getForces().getImpl()).getData()); // Reference code needs to be made const correct
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics) {
            delete dynamics;
            delete constraints;
        }
        dynamics = new ReferenceVerletDynamics(context.getSystem().getNumParticles(), static_cast<RealOpenMM>(stepSize) );
        constraints = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, constraintDistances, integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses);
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
            delete constraints;
        }
        RealOpenMM tau = static_cast<RealOpenMM>( friction == 0.0 ? 0.0 : 1.0/friction );
        dynamics = new ReferenceStochasticDynamics(
				context.getSystem().getNumParticles(), 
				static_cast<RealOpenMM>(stepSize), 
				static_cast<RealOpenMM>(tau), 
				static_cast<RealOpenMM>(temperature) );
        constraints = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, constraintDistances, integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses);
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
            delete constraints;
        }
        dynamics = new ReferenceBrownianDynamics(
				context.getSystem().getNumParticles(), 
				static_cast<RealOpenMM>(stepSize), 
				static_cast<RealOpenMM>(friction), 
				static_cast<RealOpenMM>(temperature) );
        constraints = new ReferenceShakeAlgorithm(numConstraints, constraintIndices, constraintDistances, integrator.getConstraintTolerance());
        dynamics->setReferenceConstraintAlgorithm(constraints);
        prevTemp = temperature;
        prevFriction = friction;
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem().getNumParticles(), posData, velData, forceData, masses);
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

void ReferenceApplyAndersenThermostatKernel::execute(OpenMMContextImpl& context) {
    RealOpenMM** velData = ((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData();
    thermostat->applyThermostat(
			context.getVelocities().getSize(), 
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

double ReferenceCalcKineticEnergyKernel::execute(OpenMMContextImpl& context) {
    RealOpenMM** velData = const_cast<RealOpenMM**>(((ReferenceFloatStreamImpl&) context.getVelocities().getImpl()).getData()); // Reference code needs to be made const correct
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

void ReferenceRemoveCMMotionKernel::execute(OpenMMContextImpl& context) {
    int step = (int)std::floor(context.getTime()/context.getIntegrator().getStepSize());
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
    
    // Adjust the particle velocities.
    
    momentum[0] /= static_cast<RealOpenMM>( masses.size() );
    momentum[1] /= static_cast<RealOpenMM>( masses.size() );
    momentum[2] /= static_cast<RealOpenMM>( masses.size() );
    for (size_t i = 0; i < masses.size(); ++i) {
        velData[i][0] -= static_cast<RealOpenMM>( momentum[0]/masses[i] );
        velData[i][1] -= static_cast<RealOpenMM>( momentum[1]/masses[i] );
        velData[i][2] -= static_cast<RealOpenMM>( momentum[2]/masses[i] );
    }
}
