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

#include "ReferenceFreeEnergyKernels.h"
#include "gbsa/CpuGBVISoftcore.h"
#include "gbsa/CpuObcSoftcore.h"
#include "SimTKReference/ReferenceFreeEnergyLJCoulomb14Softcore.h"
#include "SimTKReference/ReferenceFreeEnergyLJCoulombSoftcoreIxn.h"
#include "ReferenceBondForce.h"
#include "openmm/System.h"
#include "openmm/internal/ContextImpl.h"
#include "SimTKUtilities/SimTKOpenMMUtilities.h"
#include <cmath>
#include <limits>

using namespace std;
using namespace OpenMM;

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

ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel::~ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel() {
    disposeRealArray(particleParamArray, numParticles);
    disposeIntArray(exclusionArray, numParticles);
    disposeIntArray(bonded14IndexArray, num14);
    disposeRealArray(bonded14ParamArray, num14);
    if (neighborList != NULL)
        delete neighborList;
}

void ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel::initialize(const System& system, const NonbondedSoftcoreForce& force) {

    // Identify which exceptions are 1-4 interactions.

    numParticles = force.getNumParticles();
    exclusions.resize(numParticles);
    std::vector<int> nb14s;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon, softcoreLJLambda;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda);
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
        if (chargeProd != 0.0 || epsilon != 0.0)
            nb14s.push_back(i);
    }

    // Build the arrays.

    particleParamArray = allocateRealArray(numParticles, 4);
    RealOpenMM sqrtEps = static_cast<RealOpenMM>( std::sqrt(138.935485) );
    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, depth, softcoreLJLambda;
        force.getParticleParameters(i, charge, radius, depth,softcoreLJLambda);
        particleParamArray[i][0] = static_cast<RealOpenMM>(0.5*radius);
        particleParamArray[i][1] = static_cast<RealOpenMM>(2.0*sqrt(depth));
        particleParamArray[i][2] = static_cast<RealOpenMM>(charge*sqrtEps);
        particleParamArray[i][3] = static_cast<RealOpenMM>(softcoreLJLambda);
    }

    this->exclusions = exclusions;
    exclusionArray   = new int*[numParticles];
    for (int i = 0; i < numParticles; ++i) {
        exclusionArray[i]    = new int[exclusions[i].size()+1];
        exclusionArray[i][0] = exclusions[i].size();
        int index = 0;
        for (std::set<int>::const_iterator iter = exclusions[i].begin(); iter != exclusions[i].end(); ++iter)
            exclusionArray[i][++index] = *iter;
    }

    num14              = nb14s.size();
    bonded14IndexArray = allocateIntArray(num14, 2);
    bonded14ParamArray = allocateRealArray(num14, 4);
    for (int i = 0; i < num14; ++i) {
        int particle1, particle2;
        double charge, radius, depth, softcoreLJLambda;
        force.getExceptionParameters(nb14s[i], particle1, particle2, charge, radius, depth, softcoreLJLambda);
        bonded14IndexArray[i][0] = particle1;
        bonded14IndexArray[i][1] = particle2;
        bonded14ParamArray[i][0] = static_cast<RealOpenMM>(radius);
        bonded14ParamArray[i][1] = static_cast<RealOpenMM>(4.0*depth);
        bonded14ParamArray[i][2] = static_cast<RealOpenMM>(charge*sqrtEps*sqrtEps);
        bonded14ParamArray[i][3] = static_cast<RealOpenMM>(softcoreLJLambda);
    }

    nonbondedMethod  = CalcNonbondedSoftcoreForceKernel::NonbondedSoftcoreMethod(force.getNonbondedMethod());
    nonbondedCutoff  = (RealOpenMM) force.getCutoffDistance();
    //softCoreLJLambda = (RealOpenMM) force.getSoftCoreLJLambda();

    Vec3 boxVectors[3];
    system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
    periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
    periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
    periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];

    if (nonbondedMethod == NoCutoff)
        neighborList = NULL;
    else
        neighborList = new NeighborList();

#if 0
    if (nonbondedMethod == Ewald || nonbondedMethod == PME) {
        RealOpenMM ewaldErrorTol = (RealOpenMM) force.getEwaldErrorTolerance();
        ewaldAlpha = (RealOpenMM) (std::sqrt(-std::log(ewaldErrorTol))/nonbondedCutoff);
        RealOpenMM mx = periodicBoxSize[0]/nonbondedCutoff;
        RealOpenMM my = periodicBoxSize[1]/nonbondedCutoff;
        RealOpenMM mz = periodicBoxSize[2]/nonbondedCutoff;
        RealOpenMM pi = (RealOpenMM) 3.1415926535897932385;
        kmax[0] = (int)std::ceil(-(mx/pi)*std::log(ewaldErrorTol));
        kmax[1] = (int)std::ceil(-(my/pi)*std::log(ewaldErrorTol));
        kmax[2] = (int)std::ceil(-(mz/pi)*std::log(ewaldErrorTol));
        if (kmax[0]%2 == 0)
            kmax[0]++;
        if (kmax[1]%2 == 0)
            kmax[1]++;
        if (kmax[2]%2 == 0)
            kmax[2]++;
    }
    if (nonbondedMethod == Ewald || nonbondedMethod == PME) {
        RealOpenMM ewaldErrorTol = (RealOpenMM) force.getEwaldErrorTolerance();
        ewaldAlpha = (RealOpenMM) (std::sqrt(-std::log(ewaldErrorTol))/nonbondedCutoff);
        RealOpenMM mx = periodicBoxSize[0]/nonbondedCutoff;
        RealOpenMM my = periodicBoxSize[1]/nonbondedCutoff;
        RealOpenMM mz = periodicBoxSize[2]/nonbondedCutoff;
        RealOpenMM pi = (RealOpenMM) 3.1415926535897932385;
        kmax[0] = (int)std::ceil(-(mx/pi)*std::log(ewaldErrorTol));
        kmax[1] = (int)std::ceil(-(my/pi)*std::log(ewaldErrorTol));
        kmax[2] = (int)std::ceil(-(mz/pi)*std::log(ewaldErrorTol));
        if (kmax[0]%2 == 0)
            kmax[0]++;
        if (kmax[1]%2 == 0)
            kmax[1]++;
        if (kmax[2]%2 == 0)
            kmax[2]++;
    }
#endif
    rfDielectric = (RealOpenMM)force.getReactionFieldDielectric();
}

double ReferenceFreeEnergyCalcNonbondedSoftcoreForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);

    RealOpenMM energy = 0;
    ReferenceFreeEnergyLJCoulombSoftcoreIxn clj;
   // clj.setSoftCoreLJLambda( softCoreLJLambda );
    bool periodic = (nonbondedMethod == CutoffPeriodic);
    bool ewald  = (nonbondedMethod == Ewald);
    bool pme  = (nonbondedMethod == PME);
    if (nonbondedMethod != NoCutoff) {
        computeNeighborListVoxelHash(*neighborList, numParticles, posData, exclusions, periodicBoxSize, periodic || ewald || pme, nonbondedCutoff, 0.0);
        clj.setUseCutoff(nonbondedCutoff, *neighborList, rfDielectric);
    }
    if (periodic || ewald || pme)
        clj.setPeriodic(periodicBoxSize);
    if (ewald)
        clj.setUseEwald(ewaldAlpha, kmax[0], kmax[1], kmax[2]);
    if (pme)
        clj.setUsePME(ewaldAlpha);
    clj.calculatePairIxn(numParticles, posData, particleParamArray, exclusionArray, 0, forceData, 0, &energy);
    ReferenceBondForce refBondForce;
    ReferenceFreeEnergyLJCoulomb14Softcore nonbonded14;
    if (nonbondedMethod == CutoffNonPeriodic || nonbondedMethod == CutoffPeriodic)
        nonbonded14.setUseCutoff(nonbondedCutoff, rfDielectric);

    refBondForce.calculateForce(num14, bonded14IndexArray, posData, bonded14ParamArray, forceData, &energy, nonbonded14);

    return energy;
}

ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel::~ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel() {
    if (obc) {
        delete obc;
    }
}

void ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel::initialize(const System& system, const GBSAOBCSoftcoreForce& force) {

    int numParticles = system.getNumParticles();

    charges.resize(numParticles);

    std::vector<RealOpenMM> atomicRadii(numParticles);
    std::vector<RealOpenMM> scaleFactors(numParticles);
    std::vector<RealOpenMM> nonPolarScaleFactors(numParticles);

    for (int i = 0; i < numParticles; ++i) {

        double charge, radius, scalingFactor, nonPolarScaleFactor;
        force.getParticleParameters(i, charge, radius, scalingFactor, nonPolarScaleFactor);

        charges[i]              = static_cast<RealOpenMM>(charge);
        atomicRadii[i]          = static_cast<RealOpenMM>(radius);
        scaleFactors[i]         = static_cast<RealOpenMM>(scalingFactor);
        nonPolarScaleFactors[i] = static_cast<RealOpenMM>(nonPolarScaleFactor);
    }

    ObcSoftcoreParameters* obcParameters  = new ObcSoftcoreParameters(numParticles, ObcSoftcoreParameters::ObcTypeII);

    obcParameters->setAtomicRadii(atomicRadii);
    obcParameters->setScaledRadiusFactors(scaleFactors);
    obcParameters->setNonPolarScaleFactors(nonPolarScaleFactors);

    obcParameters->setSolventDielectric( static_cast<RealOpenMM>(force.getSolventDielectric()) );
    obcParameters->setSoluteDielectric(  static_cast<RealOpenMM>(force.getSoluteDielectric()) );

    // nonPolarPrefactor is in units of kJ/mol/nm^2 convert to kcal/mol/A^2 
    // to be consistent w/ polar part of OBC calculation

    RealOpenMM nonPolarPrefactor =  static_cast<RealOpenMM>(force.getNonPolarPrefactor());
    nonPolarPrefactor           /= static_cast<RealOpenMM>(418.4);
    obcParameters->setNonPolarPrefactor( nonPolarPrefactor );

    // If there is a NonbondedForce in this system, use it to initialize cutoffs and periodic boundary conditions.

    for (int i = 0; i < system.getNumForces(); i++) {
        const NonbondedForce* nonbonded = dynamic_cast<const NonbondedForce*>(&system.getForce(i));
        if (nonbonded != NULL) {
            if (nonbonded->getNonbondedMethod() != NonbondedForce::NoCutoff)
                obcParameters->setUseCutoff(static_cast<RealOpenMM>(nonbonded->getCutoffDistance()));
            if (nonbonded->getNonbondedMethod() == NonbondedForce::CutoffPeriodic) {
                Vec3 boxVectors[3];
                system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
                RealOpenMM periodicBoxSize[3];
                periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
                periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
                periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
                obcParameters->setPeriodic(periodicBoxSize);
            }
            break;
        }
    }
    obc = new CpuObcSoftcore(obcParameters);
    obc->setIncludeAceApproximation(true);
}

double ReferenceFreeEnergyCalcGBSAOBCSoftcoreForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData   = extractPositions(context);
    vector<RealVec>& forceData = extractForces(context);
    obc->computeImplicitSolventForces(posData, &charges[0], forceData, 1);
    return obc->getEnergy();
}

ReferenceFreeEnergyCalcGBVISoftcoreForceKernel::~ReferenceFreeEnergyCalcGBVISoftcoreForceKernel() {
    if (gbviSoftcore) {
        delete gbviSoftcore;
    }
}

void ReferenceFreeEnergyCalcGBVISoftcoreForceKernel::initialize(const System& system, const GBVISoftcoreForce& force, const std::vector<double> & inputScaledRadii ) {

    int numParticles = system.getNumParticles();

    charges.resize(numParticles);
    std::vector<RealOpenMM> atomicRadii(numParticles);
    std::vector<RealOpenMM> scaledRadii(numParticles);
    std::vector<RealOpenMM> gammas(numParticles);
    std::vector<RealOpenMM> bornRadiusScaleFactors(numParticles);

    for (int i = 0; i < numParticles; ++i) {
        double charge, radius, gamma, bornRadiusScaleFactor;
        force.getParticleParameters(i, charge, radius, gamma, bornRadiusScaleFactor);
        charges[i]                = static_cast<RealOpenMM>(charge);
        atomicRadii[i]            = static_cast<RealOpenMM>(radius);
        gammas[i]                 = static_cast<RealOpenMM>(gamma);
        scaledRadii[i]            = static_cast<RealOpenMM>(inputScaledRadii[i]);
        bornRadiusScaleFactors[i] = static_cast<RealOpenMM>(bornRadiusScaleFactor);
    }

    GBVISoftcoreParameters* gBVIParameters = new GBVISoftcoreParameters(numParticles);
    gBVIParameters->setAtomicRadii(atomicRadii);
    gBVIParameters->setGammaParameters(gammas);
    gBVIParameters->setBornRadiusScaleFactors(bornRadiusScaleFactors);
    gBVIParameters->setScaledRadii(scaledRadii);

    // switching function/scaling

    // quintic spline

    if( force.getBornRadiusScalingMethod() == GBVISoftcoreForce::QuinticSpline ){
        gBVIParameters->setBornRadiusScalingSoftcoreMethod( GBVISoftcoreParameters::QuinticSpline );
        gBVIParameters->setQuinticLowerLimitFactor(         static_cast<RealOpenMM>(force.getQuinticLowerLimitFactor()) );
        gBVIParameters->setQuinticUpperBornRadiusLimit(     static_cast<RealOpenMM>(force.getQuinticUpperBornRadiusLimit()) );
    }

    gBVIParameters->setSolventDielectric( static_cast<RealOpenMM>(force.getSolventDielectric()) );
    gBVIParameters->setSoluteDielectric( static_cast<RealOpenMM>(force.getSoluteDielectric()) );

    if (force.getNonbondedMethod() != GBVISoftcoreForce::NoCutoff)
        gBVIParameters->setUseCutoff(static_cast<RealOpenMM>(force.getCutoffDistance()));
    if (force.getNonbondedMethod() == GBVISoftcoreForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        RealOpenMM periodicBoxSize[3];
        periodicBoxSize[0] = (RealOpenMM) boxVectors[0][0];
        periodicBoxSize[1] = (RealOpenMM) boxVectors[1][1];
        periodicBoxSize[2] = (RealOpenMM) boxVectors[2][2];
        gBVIParameters->setPeriodic(periodicBoxSize);
    }    
    gbviSoftcore = new CpuGBVISoftcore(gBVIParameters);
}

double ReferenceFreeEnergyCalcGBVISoftcoreForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& posData = extractPositions(context);

    vector<RealOpenMM> bornRadii(context.getSystem().getNumParticles());
    gbviSoftcore->computeBornRadii(posData, bornRadii );
    if (includeForces) {
        vector<RealVec>& forceData = extractForces(context);
        gbviSoftcore->computeBornForces(bornRadii, posData, &charges[0], forceData);
    }
    RealOpenMM energy = 0.0;
    if (includeEnergy)
        energy = gbviSoftcore->computeBornEnergy(bornRadii, posData, &charges[0]);
    return static_cast<double>(energy);
}
