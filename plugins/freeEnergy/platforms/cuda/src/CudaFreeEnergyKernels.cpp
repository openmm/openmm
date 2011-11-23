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
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "CudaFreeEnergyKernels.h"
#include "CudaForceInfo.h"

#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "kernels/gputypes.h"
#include "kernels/cudaKernels.h"
#include "kernels/GpuFreeEnergyCudaKernels.h" 

#include <cmath>
#include <map>
#include <cstring>
#include <cstdlib>
#include <typeinfo>

using namespace OpenMM;

typedef std::map< std::string, int > MapStringInt;
typedef MapStringInt::iterator MapStringIntI;
typedef MapStringInt::const_iterator MapStringIntCI;

// force names

const std::string HARMONIC_BOND_FORCE             = "HarmonicBond";
const std::string HARMONIC_ANGLE_FORCE            = "HarmonicBond";
const std::string PERIODIC_TORSION_FORCE          = "PeriodicTorsion";
const std::string RB_TORSION_FORCE                = "RbTorsion";

const std::string NB_FORCE                        = "Nb";
const std::string NB_SOFTCORE_FORCE               = "NbSoftcore";

const std::string NB_EXCEPTION_FORCE              = "NbException";
const std::string NB_EXCEPTION_SOFTCORE_FORCE     = "NbSoftcoreException";

const std::string GBSA_OBC_FORCE                  = "Obc";
const std::string GBSA_OBC_SOFTCORE_FORCE         = "ObcSoftcore";

const std::string GBVI_FORCE                      = "GBVI";
const std::string GBVI_SOFTCORE_FORCE             = "GBVISoftcore";

static void getForceMap(const System& system, MapStringInt& forceMap, FILE* log) { 

    // check forces and relevant parameters

    for(int i = 0; i < system.getNumForces(); ++i) {

        int hit                 = 0;
        const Force& force      = system.getForce(i);
         std::string forceName  = "NA";

        // bond

        if( !hit ){

            try {
               const HarmonicBondForce& harmonicBondForce = dynamic_cast<const HarmonicBondForce&>(force);
               forceMap[HARMONIC_BOND_FORCE]              = 1;
               forceName                                  = HARMONIC_BOND_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // angle

        if( !hit ){

            try {
               const HarmonicAngleForce& harmonicAngleForce = dynamic_cast<const HarmonicAngleForce&>(force);
               forceMap[HARMONIC_ANGLE_FORCE]               = 1;
               forceName                                    = HARMONIC_ANGLE_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // PeriodicTorsionForce

        if( !hit ){

            try {
               const PeriodicTorsionForce & periodicTorsionForce = dynamic_cast<const PeriodicTorsionForce&>(force);
               forceMap[PERIODIC_TORSION_FORCE]                  = 1;
               forceName                                         = PERIODIC_TORSION_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // RBTorsionForce

        if( !hit ){
            try {
               const RBTorsionForce& rBTorsionForce = dynamic_cast<const RBTorsionForce&>(force);
               forceMap[RB_TORSION_FORCE]           = 1;
               forceName                            = RB_TORSION_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // nonbonded

        if( !hit ){
            try {
               const NonbondedForce& nbForce = dynamic_cast<const NonbondedForce&>(force);
               forceMap[NB_FORCE]            = 1;
               forceName                     = NB_FORCE;
            } catch( std::bad_cast ){
            }
        }

        // nonbonded softcore

        if( !hit ){
            try {
               const NonbondedSoftcoreForce& nbForce = dynamic_cast<const NonbondedSoftcoreForce&>(force);
               forceMap[NB_SOFTCORE_FORCE]           = 1;
               forceName                             = NB_SOFTCORE_FORCE;
            } catch( std::bad_cast ){
            }
        }

        // GBSA OBC

        if( !hit ){
            try {
               const GBSAOBCForce& obcForce       = dynamic_cast<const GBSAOBCForce&>(force);
               forceMap[GBSA_OBC_FORCE]           = 1;
               forceName                          = GBSA_OBC_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // GBSA OBC softcore

        if( !hit ){
            try {
               const GBSAOBCSoftcoreForce& obcForce = dynamic_cast<const GBSAOBCSoftcoreForce&>(force);
               forceMap[GBSA_OBC_SOFTCORE_FORCE]    = 1;
               forceName                            = GBSA_OBC_SOFTCORE_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // GB/VI

        if( !hit ){
            try {
               const GBVIForce& obcForce  = dynamic_cast<const GBVIForce&>(force);
               forceMap[GBVI_FORCE]       = 1;
               forceName                  = GBVI_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // GB/VI softcore

        if( !hit ){
            try {
               const GBVISoftcoreForce& gbviForce = dynamic_cast<const GBVISoftcoreForce&>(force);
               forceMap[GBVI_SOFTCORE_FORCE]      = 1;
               forceName                          = GBVI_SOFTCORE_FORCE;
               hit++;
            } catch( std::bad_cast ){
            }
        }

        if( log ){
            (void) fprintf( log, "Map: Force %d %s\n", i, forceName.c_str() );
        }
     }
}

class CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const NonbondedSoftcoreForce& force) : force(force) {
    }    
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, sigma1, sigma2, epsilon1, epsilon2, softcoreLJLambda1, softcoreLJLambda2;
        force.getParticleParameters(particle1, charge1, sigma1, epsilon1, softcoreLJLambda1);
        force.getParticleParameters(particle2, charge2, sigma2, epsilon2, softcoreLJLambda2);
        return (charge1 == charge2 && sigma1 == sigma2 && epsilon1 == epsilon2 && softcoreLJLambda1 == softcoreLJLambda2);
    }    
    int getNumParticleGroups() {
        return force.getNumExceptions();
    }    
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon, softcoreLJLambda;
        force.getExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }    
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double chargeProd1, chargeProd2, sigma1, sigma2, epsilon1, epsilon2, softcoreLJLambda1, softcoreLJLambda2;
        force.getExceptionParameters(group1, particle1, particle2, chargeProd1, sigma1, epsilon1, softcoreLJLambda1);
        force.getExceptionParameters(group2, particle1, particle2, chargeProd2, sigma2, epsilon2, softcoreLJLambda2);
        return (chargeProd1 == chargeProd2 && sigma1 == sigma2 && epsilon1 == epsilon2 && softcoreLJLambda1 == softcoreLJLambda2);
    }    
private:
    const NonbondedSoftcoreForce& force;
};

CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::~CudaFreeEnergyCalcNonbondedSoftcoreForceKernel() {
    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "~CudaFreeEnergyCalcNonbondedSoftcoreForceKernel called.\n" );
        (void) fflush( data.getLog() );
    }
    data.decrementKernelCount();
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::initialize(const System& system, const NonbondedSoftcoreForce& force) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::initialize";

// ---------------------------------------------------------------------------------------

    if( data.getLog() ){
        (void) fprintf( data.getLog(), "%s called.\n", methodName.c_str() );
        (void) fflush( data.getLog() );
    }

    // check forces and relevant parameters

    MapStringInt forceMap;
    getForceMap( system, forceMap, data.getLog() );

    int softcore        = 0;
    if( forceMap.find( GBSA_OBC_FORCE ) != forceMap.end() ){
       setIncludeGBSA( true );
    }

    if( forceMap.find( GBSA_OBC_SOFTCORE_FORCE ) != forceMap.end() ){
       setIncludeGBSA( true );
       softcore++;
    }

    if( forceMap.find( GBVI_FORCE ) != forceMap.end() ){
       setIncludeGBVI( true );
    }
    if( forceMap.find( GBVI_SOFTCORE_FORCE ) != forceMap.end() ){
       setIncludeGBVI( true );
       softcore++;
    }

    if( forceMap.find( NB_SOFTCORE_FORCE ) != forceMap.end() ){
       softcore++;
    }

    setIncludeSoftcore( softcore );

    numParticles      = force.getNumParticles();

    // Identify which exceptions are 1-4 interactions.

    std::vector<pair<int, int> > exclusions;
    std::vector<int> exceptions;
    for (int i = 0; i < force.getNumExceptions(); i++) {
        int particle1, particle2;
        double chargeProd, sigma, epsilon, softcoreLJLambda;
        force.getExceptionParameters(i, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda);
        exclusions.push_back(pair<int, int>(particle1, particle2));
        if (chargeProd != 0.0 || epsilon != 0.0)
            exceptions.push_back(i);
    }

    // Initialize nonbonded interactions.
    
    if( numParticles > 0 ){
        std::vector<int> particle(numParticles);
        std::vector<float> c6(numParticles);
        std::vector<float> c12(numParticles);
        std::vector<float> q(numParticles);
        std::vector<float> softcoreLJLambdaArray(numParticles);
        std::vector<char> symbol;
        std::vector<std::vector<int> > exclusionList(numParticles);
        for (int i = 0; i < numParticles; i++) {
            double charge, radius, depth, softcoreLJLambda;
            force.getParticleParameters(i, charge, radius, depth, softcoreLJLambda);
            particle[i]              = i;
            q[i]                     = static_cast<float>( charge );
            c6[i]                    = static_cast<float>( (4*depth*pow(radius, 6.0)) );
            c12[i]                   = static_cast<float>( (4*depth*pow(radius, 12.0)) );
            softcoreLJLambdaArray[i] = static_cast<float>( softcoreLJLambda );
            exclusionList[i].push_back(i);
        }

        for (int i = 0; i < (int)exclusions.size(); i++) {
            exclusionList[exclusions[i].first].push_back(exclusions[i].second);
            exclusionList[exclusions[i].second].push_back(exclusions[i].first);
        }
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        freeEnergyGpuSetPeriodicBoxSize( data.getFreeEnergyGpu(), static_cast<float>(boxVectors[0][0] ), static_cast<float>(boxVectors[1][1] ), static_cast<float>(boxVectors[2][2] ));

        CudaFreeEnergyNonbondedMethod method = FREE_ENERGY_NO_CUTOFF;
        if (force.getNonbondedMethod() != NonbondedSoftcoreForce::NoCutoff) {
            method = FREE_ENERGY_CUTOFF;
        }
        if (force.getNonbondedMethod() == NonbondedSoftcoreForce::CutoffPeriodic) {
            method = FREE_ENERGY_PERIODIC;
        }

        // setup parameters

        gpuSetNonbondedSoftcoreParameters( data.getFreeEnergyGpu(), 138.935485f, particle, c6, c12, q,
                                           softcoreLJLambdaArray, symbol, exclusionList, method,
                                           static_cast<float>(force.getCutoffDistance() ), static_cast<float>(force.getReactionFieldDielectric()));
    }

    // Initialize 1-4 nonbonded interactions.

    numExceptions = exceptions.size();
    if( numExceptions > 0 ){
        std::vector<int> particle1(numExceptions);
        std::vector<int> particle2(numExceptions);
        std::vector<float> c6(numExceptions);
        std::vector<float> c12(numExceptions);
        std::vector<float> qProd(numExceptions);
        std::vector<float> softcoreLJLambdaArray(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double charge, sig, eps, softcoreLJLambda;
            force.getExceptionParameters(exceptions[i], particle1[i], particle2[i], charge, sig, eps, softcoreLJLambda);
            c6[i]                    = static_cast<float>( (4.0*eps*pow(sig, 6.0)) );
            c12[i]                   = static_cast<float>( (4.0*eps*pow(sig, 12.0)) );
            qProd[i]                 = static_cast<float>( charge );
            softcoreLJLambdaArray[i] = static_cast<float>( softcoreLJLambda );
        }
        gpuSetLJ14SoftcoreParameters( data.getFreeEnergyGpu(), 138.935485f, particle1, particle2, c6, c12, qProd, softcoreLJLambdaArray);
    } else if( data.getLog() ){
        (void) fprintf( data.getLog(), "Mo nonbonded softcore exceptions.\n" );
        (void) fflush( data.getLog() );
    }

    data.getFreeEnergyGpu()->gpuContext->forces.push_back(new ForceInfo(force));
}

double CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::execute( ContextImpl& context, bool includeForces, bool includeEnergy ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

    freeEnergyGpuContext gpu = data.getFreeEnergyGpu();

    data.initializeGpu( );
 
    // calculate nonbonded ixns here, only if implicit solvent is inactive

    if ( !getIncludeGBSA() && !getIncludeGBVI() ) {
        kCalculateCDLJSoftcoreForces(gpu);
    }

    // local LJ-14 forces

    if( getNumExceptions() > 0 ){
        kCalculateLocalSoftcoreForces(gpu);
    }

    return 0.0;
}

bool CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getIncludeGBSA( void ) const {
    return bIncludeGBSA;
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::setIncludeGBSA( bool inputIncludeGBSA ){
    bIncludeGBSA = inputIncludeGBSA;
}

bool CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getIncludeGBVI( void ) const {
    return bIncludeGBVI;
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::setIncludeGBVI( bool inputIncludeGBVI ){
    bIncludeGBVI = inputIncludeGBVI;
}

int CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getIncludeSoftcore( void ) const {
    return includeSoftcore;
}

int CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getNumExceptions( void ) const {
    return numExceptions;
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::setIncludeSoftcore( int inputIncludeSoftcore ){
    includeSoftcore = inputIncludeSoftcore;
}

class CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const GBSAOBCSoftcoreForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, scale1, scale2, particleNonPolarScalingFactor1, particleNonPolarScalingFactor2;
        force.getParticleParameters(particle1, charge1, radius1, scale1, particleNonPolarScalingFactor1);
        force.getParticleParameters(particle2, charge2, radius2, scale2, particleNonPolarScalingFactor2);
        return (charge1 == charge2 && radius1 == radius2 && scale1 == scale2 && particleNonPolarScalingFactor1 == particleNonPolarScalingFactor2);
    }
private:
    const GBSAOBCSoftcoreForce& force;
};

CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::~CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel() {
    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "~CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel called.\n" );
        (void) fflush( data.getLog() );
    }
    data.decrementKernelCount();
}

void CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::initialize(const System& system, const GBSAOBCSoftcoreForce& force) {

// ---------------------------------------------------------------------------------------

    freeEnergyGpuContext gpu  = data.getFreeEnergyGpu();

    MapStringInt forceMap;
    getForceMap( system, forceMap, log);

    // check that nonbonded (non-softcore is not active)

    if( forceMap.find( NB_FORCE ) != forceMap.end() ){ 
        throw OpenMMException( "Mixing NonbondedForce and GBSAOBCSoftoreForce is not allowed -- use NonbondedSoftcoreForce " );
    }
    if( forceMap.find( NB_SOFTCORE_FORCE ) == forceMap.end() ){ 
        throw OpenMMException( "NonbondedSoftcore force must be included w/ GBSAOBCSoftcore force." );
    }

    int numParticles = system.getNumParticles();

    std::vector<float> radius(numParticles);
    std::vector<float> scale(numParticles);
    std::vector<float> charge(numParticles);
    std::vector<float> nonPolarScalingFactors(numParticles);

    for( int ii = 0; ii < numParticles; ii++ ){
        double particleCharge, particleRadius, scalingFactor, particleNonPolarScalingFactor;
        force.getParticleParameters( ii, particleCharge, particleRadius, scalingFactor, particleNonPolarScalingFactor);
        radius[ii]                 = static_cast<float>( particleRadius);
        scale[ii]                  = static_cast<float>( scalingFactor);
        charge[ii]                 = static_cast<float>( particleCharge);
        nonPolarScalingFactors[ii] = static_cast<float>( particleNonPolarScalingFactor);
    }

    gpuSetObcSoftcoreParameters( gpu, static_cast<float>( force.getSoluteDielectric()),
                                 static_cast<float>( force.getSolventDielectric()),
                                 static_cast<float>( force.getNonPolarPrefactor()),
                                 radius, scale, charge, nonPolarScalingFactors );

    data.getFreeEnergyGpu()->gpuContext->forces.push_back(new ForceInfo(force));
    return;

}

double CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

// ---------------------------------------------------------------------------------------

    freeEnergyGpuContext freeEnergyGpu = data.getFreeEnergyGpu();
    gpuContext gpu                     = freeEnergyGpu->gpuContext;
    int call = 0;
    
    // send address's of arrays, ... to device on first call
    // required since force/energy buffers not set when CudaFreeEnergyCalcGBVISoftcoreForceKernel::initialize() was called

    data.initializeGpu( );

    // (1) clear Born force array
    // (2) calculate Born radii and sum
    // (3) loop 1
    // (4) sum/calculate Born forces
    // (5) loop 2

    kClearSoftcoreBornForces(gpu);
    kCalculateObcGbsaSoftcoreBornSum( freeEnergyGpu );
    kReduceObcGbsaSoftcoreBornSum(gpu);
    kCalculateCDLJObcGbsaSoftcoreForces1( freeEnergyGpu );

    // sum Born forces and  execute second OBC loop

    kReduceObcGbsaSoftcoreBornForces(gpu);
    kCalculateObcGbsaSoftcoreForces2( freeEnergyGpu );

    if( data.getLog() ){
        kPrintObcGbsaSoftcore( freeEnergyGpu, "Post kCalculateObcGbsaSoftcoreForces2", call, data.getLog() );
    }

    return 0.0;
}

class CudaFreeEnergyCalcGBVISoftcoreForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const GBVISoftcoreForce& force) : force(force) {
    }    
    bool areParticlesIdentical(int particle1, int particle2) {
        double charge1, charge2, radius1, radius2, gamma1, gamma2, bornRadiusScaleFactor1, bornRadiusScaleFactor2;
        force.getParticleParameters(particle1, charge1, radius1, gamma1, bornRadiusScaleFactor1);
        force.getParticleParameters(particle2, charge2, radius2, gamma2, bornRadiusScaleFactor2);
        return (charge1 == charge2 && radius1 == radius2 && gamma1 == gamma2 && bornRadiusScaleFactor1 == bornRadiusScaleFactor2);
    }    
private:
    const GBVISoftcoreForce& force;
};

CudaFreeEnergyCalcGBVISoftcoreForceKernel::~CudaFreeEnergyCalcGBVISoftcoreForceKernel() {
    if( 0 && data.getLog() ){
        (void) fprintf( data.getLog(), "~CudaFreeEnergyCalcGBVISoftcoreForceKernel called.\n" );
        (void) fflush( data.getLog() );
    }
    data.decrementKernelCount();
}

void CudaFreeEnergyCalcGBVISoftcoreForceKernel::initialize(const System& system, const GBVISoftcoreForce& force, const std::vector<double> & inputScaledRadii) {

// ---------------------------------------------------------------------------------------

    int numParticles          = system.getNumParticles();
    freeEnergyGpuContext gpu  = data.getFreeEnergyGpu();

    // check forces and relevant parameters

    MapStringInt forceMap;
    getForceMap( system, forceMap, log);

    // check that nonbonded (non-softcore is not active)

    if( forceMap.find( NB_FORCE ) != forceMap.end() ){ 
        throw OpenMMException( "Mixing NonbondedForce and GBVISoftoreForce not allowed -- use NonbondedSoftcoreForce " );
    }

    std::vector<int>   particle(numParticles);
    std::vector<float> radius(numParticles);
    std::vector<float> scaledRadii(numParticles);
    std::vector<float> gammas(numParticles);
    std::vector<float> bornRadiusScaleFactors(numParticles);

    for( int i = 0; i < numParticles; i++ ){
        double charge, particleRadius, gamma, bornRadiusScaleFactor;
        force.getParticleParameters(i, charge, particleRadius, gamma, bornRadiusScaleFactor);
        particle[i]                  = i;
        radius[i]                    = static_cast<float>( particleRadius );
        gammas[i]                    = static_cast<float>( gamma );
        scaledRadii[i]               = static_cast<float>( inputScaledRadii[i] );
        bornRadiusScaleFactors[i]    = static_cast<float>( bornRadiusScaleFactor );
    }

    std::vector<float> quinticSplineParameters;
    if( force.getBornRadiusScalingMethod() == GBVISoftcoreForce::QuinticSpline ){

        // quintic spline

        quinticSplineParameters.resize(2);
        quinticSplineParameters[0] = static_cast<float>(force.getQuinticLowerLimitFactor());
        quinticSplineParameters[1] = static_cast<float>(force.getQuinticUpperBornRadiusLimit());
        quinticSplineParameters[1] = powf( quinticSplineParameters[1], -3.0f ); 
        quinticScaling =  1;
    }

    // load parameters onto board
    // defined in kCalculateGBVISoftcore.cu

    gpuSetGBVISoftcoreParameters( gpu, static_cast<float>( force.getSoluteDielectric() ), static_cast<float>( force.getSolventDielectric() ),
                                  particle, radius, gammas, scaledRadii, bornRadiusScaleFactors, quinticSplineParameters);

    data.getFreeEnergyGpu()->gpuContext->forces.push_back(new ForceInfo(force));

    return;
}

double CudaFreeEnergyCalcGBVISoftcoreForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {

    freeEnergyGpuContext freeEnergyGpu = data.getFreeEnergyGpu();
    gpuContext gpu                     = freeEnergyGpu->gpuContext;
    
    // send address's of arrays, ... to device on first call
    // required since force/energy buffers not set when CudaFreeEnergyCalcGBVISoftcoreForceKernel::initialize() was called

    data.initializeGpu( );

    // (1) clear Born force array
    // (2) calculate Born radii and sum
    // (3) loop 1
    // (4) sum/calculate Born forces
    // (5) loop 2

    // calculate Born radii and first loop of GB/VI forces

    kClearSoftcoreBornForces(gpu);

    kCalculateGBVISoftcoreBornSum( freeEnergyGpu );

    if( quinticScaling ){
        kReduceGBVIBornSumQuinticScaling( freeEnergyGpu );
    } else {
        kReduceGBVISoftcoreBornSum( freeEnergyGpu );
    }

    kCalculateCDLJObcGbsaSoftcoreForces1( freeEnergyGpu );

    if( quinticScaling ){
        kReduceGBVIBornForcesQuinticScaling(freeEnergyGpu);
    } else {
        kReduceGBVISoftcoreBornForces( freeEnergyGpu );
    }

    // second loop of GB/VI forces

    kCalculateGBVISoftcoreForces2( freeEnergyGpu );
    if( data.getLog() ){
        kPrintGBVISoftcore( freeEnergyGpu, "Post kCalculateGBVISoftcoreForces2", 0, data.getLog() );
    }

    return 0.0;
}
