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
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "kernels/gputypes.h"
#include "kernels/cudaKernels.h"
#include "kernels/GpuFreeEnergyCudaKernels.h" 
#include "kernels/GpuLJ14Softcore.h" 

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

CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::~CudaFreeEnergyCalcNonbondedSoftcoreForceKernel() {
    if( log ){
        (void) fprintf( log, "CudaFreeEnergyCalcNonbondedSoftcoreForceKernel destructor called.\n" );
        (void) fflush( log );
    }
    delete gpuNonbondedSoftcore;
    delete gpuLJ14Softcore;
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::initialize(const System& system, const NonbondedSoftcoreForce& force) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::initialize";

// ---------------------------------------------------------------------------------------

    if( log ){
        (void) fprintf( log, "%s called.\n", methodName.c_str() );
        (void) fflush( log );
    }

    // check forces and relevant parameters

    MapStringInt forceMap;
    getForceMap( system, forceMap, log);

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
    _gpuContext* gpu  = data.gpu;

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
    
    {
        std::vector<int> particle(numParticles);
        std::vector<float> c6(numParticles);
        std::vector<float> c12(numParticles);
        std::vector<float> q(numParticles);
        std::vector<float> softcoreLJLambdaArray(numParticles);
        std::vector<char> symbol;
        std::vector<std::vector<int> > exclusionList(numParticles);
        float minSoftcoreLJLambda = 1.0e+20f;
        for (int i = 0; i < numParticles; i++) {
            double charge, radius, depth, softcoreLJLambda;
            force.getParticleParameters(i, charge, radius, depth, softcoreLJLambda);
            particle[i]              = i;
            q[i]                     = static_cast<float>( charge );
            c6[i]                    = static_cast<float>( (4*depth*pow(radius, 6.0)) );
            c12[i]                   = static_cast<float>( (4*depth*pow(radius, 12.0)) );
            softcoreLJLambdaArray[i] = static_cast<float>( softcoreLJLambda );
            if( minSoftcoreLJLambda > softcoreLJLambda ){
               minSoftcoreLJLambda = softcoreLJLambda;
            }
            exclusionList[i].push_back(i);
        }

        for (int i = 0; i < (int)exclusions.size(); i++) {
            exclusionList[exclusions[i].first].push_back(exclusions[i].second);
            exclusionList[exclusions[i].second].push_back(exclusions[i].first);
        }
        Vec3 boxVectors[3];
        system.getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        gpuSetPeriodicBoxSize(gpu, static_cast<float>(boxVectors[0][0] ), static_cast<float>(boxVectors[1][1] ), static_cast<float>(boxVectors[2][2] ));
        CudaNonbondedMethod method = NO_CUTOFF;
        if (force.getNonbondedMethod() != NonbondedSoftcoreForce::NoCutoff) {
            gpuSetNonbondedCutoff(gpu, static_cast<float>(force.getCutoffDistance() ), force.getReactionFieldDielectric());
            method = CUTOFF;
        }
        if (force.getNonbondedMethod() == NonbondedSoftcoreForce::CutoffPeriodic) {
            method = PERIODIC;
        }
        if (force.getNonbondedMethod() == NonbondedSoftcoreForce::Ewald || force.getNonbondedMethod() == NonbondedSoftcoreForce::PME) {
            double ewaldErrorTol = force.getEwaldErrorTolerance();
            double alpha = (1.0/force.getCutoffDistance())*std::sqrt(-std::log(ewaldErrorTol));
            double mx = boxVectors[0][0]/force.getCutoffDistance();
            double my = boxVectors[1][1]/force.getCutoffDistance();
            double mz = boxVectors[2][2]/force.getCutoffDistance();
            double pi = 3.1415926535897932385;
            int kmaxx = (int)std::ceil(-(mx/pi)*std::log(ewaldErrorTol));
            int kmaxy = (int)std::ceil(-(my/pi)*std::log(ewaldErrorTol));
            int kmaxz = (int)std::ceil(-(mz/pi)*std::log(ewaldErrorTol));
            if (force.getNonbondedMethod() == NonbondedSoftcoreForce::Ewald) {
                if (kmaxx%2 == 0)
                    kmaxx++;
                if (kmaxy%2 == 0)
                    kmaxy++;
                if (kmaxz%2 == 0)
                    kmaxz++;
                gpuSetEwaldParameters(gpu, static_cast<float>( alpha ), kmaxx, kmaxy, kmaxz);
                method = EWALD;
            }
            else {
                int gridSizeX = kmaxx*3;
                int gridSizeY = kmaxy*3;
                int gridSizeZ = kmaxz*3;
                gridSizeX = ((gridSizeX+3)/4)*4;
                gridSizeY = ((gridSizeY+3)/4)*4;
                gridSizeZ = ((gridSizeZ+3)/4)*4;
                gpuSetPMEParameters(gpu, static_cast<float>( alpha ), gridSizeX, gridSizeY, gridSizeZ);
                method = PARTICLE_MESH_EWALD;
            }
        }
        data.nonbondedMethod = method;

        // setup parameters

        gpuNonbondedSoftcore = gpuSetNonbondedSoftcoreParameters(gpu, 138.935485f, particle, c6, c12, q,
                                                                 softcoreLJLambdaArray, symbol, exclusionList, method);

        // Compute the Ewald self energy.

        data.ewaldSelfEnergy = 0.0;
        if (force.getNonbondedMethod() == NonbondedSoftcoreForce::Ewald || force.getNonbondedMethod() == NonbondedSoftcoreForce::PME) {
            double selfEnergyScale = gpu->sim.epsfac*gpu->sim.alphaEwald/std::sqrt(PI);
                for (int i = 0; i < numParticles; i++)
                    data.ewaldSelfEnergy -= selfEnergyScale*q[i]*q[i];
        }
    }

    // Initialize 1-4 nonbonded interactions.
    
    {
        numExceptions = exceptions.size();
        std::vector<int> particle1(numExceptions);
        std::vector<int> particle2(numExceptions);
        std::vector<float> c6(numExceptions);
        std::vector<float> c12(numExceptions);
        std::vector<float> q1(numExceptions);
        std::vector<float> q2(numExceptions);
        std::vector<float> softcoreLJLambdaArray(numExceptions);
        for (int i = 0; i < numExceptions; i++) {
            double charge, sig, eps, softcoreLJLambda;
            force.getExceptionParameters(exceptions[i], particle1[i], particle2[i], charge, sig, eps, softcoreLJLambda);
            c6[i]                    = static_cast<float>( (4*eps*pow(sig, 6.0)) );
            c12[i]                   = static_cast<float>( (4*eps*pow(sig, 12.0)) );
            q1[i]                    = static_cast<float>( charge );
            q2[i]                    = 1.0f;
            softcoreLJLambdaArray[i] = static_cast<float>( softcoreLJLambda );
        }
        gpuLJ14Softcore = gpuSetLJ14SoftcoreParameters(gpu, 138.935485f, 1.0f, particle1, particle2, c6, c12, q1, q2, softcoreLJLambdaArray);
    }
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::executeForces(ContextImpl& context) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

    _gpuContext* gpu = data.gpu;
//    static int call  = 0;
//    call++;

    // write array, ... address's to board

    if( setSim == 0 ){
        setSim++;
        if( log ){
            (void) fprintf( log, "%s Obc=%d GB/VI=%d exceptions=%d\n",
                            methodName.c_str(), getIncludeGBSA(), getIncludeGBVI(), getNumExceptions() );
            (void) fflush( log );
        }
        SetCalculateLocalSoftcoreGpuSim( gpu );
        SetCalculateCDLJSoftcoreGpuSim( gpu );

        // flip strides (unsure if this is needed)

#if 0
        (void) fprintf( stderr, "flipping gpuLJ14Softcore\n" ); fflush( stderr );
        GpuLJ14Softcore* gpuLJ14Softcore = getGpuLJ14Softcore( );
        if( gpuLJ14Softcore ){
            gpuLJ14Softcore->flipStrides( gpu );
            if( log ){
                (void) fprintf( log, "CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::executeForces flipping LJ14\n" );
                (void) fflush( log );
            }
        }
#endif

    }
 
    // calculate nonbonded ixns here, only if implicit solvent is inactive

    if ( !getIncludeGBSA() && !getIncludeGBVI() ) {
        kCalculateCDLJSoftcoreForces(gpu);
    }
  
    // local LJ-14 forces

    kCalculateLocalSoftcoreForces(gpu);
//kPrintForces(gpu, "Post kCalculateLocalSoftcoreForces", call );
//kReduceForces(gpu);
}

double CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::executeEnergy(ContextImpl& context) {
    executeForces(context);
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

bool CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getIncludeSoftcore( void ) const {
    return bIncludeSoftcore;
}

int CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getNumExceptions( void ) const {
    return numExceptions;
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::setIncludeSoftcore( bool inputIncludeSoftcore ){
    bIncludeSoftcore = inputIncludeSoftcore;
}

GpuLJ14Softcore* CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::getGpuLJ14Softcore( void ) const {
    return gpuLJ14Softcore;
}

void CudaFreeEnergyCalcNonbondedSoftcoreForceKernel::setGpuLJ14Softcore( GpuLJ14Softcore* inputGpuLJ14Softcore ){
    gpuLJ14Softcore = inputGpuLJ14Softcore;
}

CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::~CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel() {
    delete gpuObcGbsaSoftcore;
}

void CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::initialize(const System& system, const GBSAOBCSoftcoreForce& force) {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::initialize";

// ---------------------------------------------------------------------------------------

    _gpuContext* gpu = data.gpu;

    MapStringInt forceMap;
    getForceMap( system, forceMap, log);

    // check that nonbonded (non-softcore is not active)

    if( forceMap.find( NB_FORCE ) != forceMap.end() ){ 
        throw OpenMMException( "Mixing NonbondedForce and GBSAOBCSoftoreForce not allowed -- use NonbondedSoftcoreForce " );
    }
    if( forceMap.find( NB_SOFTCORE_FORCE ) == forceMap.end() ){ 
        throw OpenMMException( "NonbondedSoftcore force must be included w/ GBSAOBCSoftcore force." );
    }

    int numParticles = system.getNumParticles();

    std::vector<float> radius(numParticles);
    std::vector<float> scale(numParticles);
    std::vector<float> charge(numParticles);
    std::vector<float> nonPolarScalingFactors(numParticles);

    for (int i = 0; i < numParticles; i++) {
        double particleCharge, particleRadius, scalingFactor, particleNonPolarScalingFactor;
        force.getParticleParameters(i, particleCharge, particleRadius, scalingFactor, particleNonPolarScalingFactor);
        radius[i]                 = static_cast<float>( particleRadius);
        scale[i]                  = static_cast<float>( scalingFactor);
        charge[i]                 = static_cast<float>( particleCharge);
        nonPolarScalingFactors[i] = static_cast<float>( particleNonPolarScalingFactor);
    }

    gpuObcGbsaSoftcore = gpuSetObcSoftcoreParameters(gpu, static_cast<float>( force.getSoluteDielectric()),
                                                     static_cast<float>( force.getSolventDielectric()),
                                                     static_cast<float>( force.getNonPolarPrefactor()),
                                                     radius, scale, charge, nonPolarScalingFactors );
}

void CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::executeForces(ContextImpl& context) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

    _gpuContext* gpu = data.gpu;

    int debug        = 1;
    static int call  = 0;

    // send address's of arrays, ... to device on first call
    // required since force/energy buffers not set when CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::initialize() was called

    if( setSim == 0 ){
       setSim++;
       SetCalculateObcGbsaSoftcoreBornSumSim( gpu );
       SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( gpu );
       SetCalculateObcGbsaSoftcoreForces2Sim( gpu );
    }

    // required!!
    gpu->bRecalculateBornRadii = true;

    // calculate Born radii and first loop of Obc forces

    if( debug && log ){
        call++;
        if( log ){
            (void) fprintf( log, "\n%s: calling kCalculateCDLJObcGbsaSoftcoreForces1\n", methodName.c_str() );
            (void) fflush( log );
        }
    }

    kClearBornForces(gpu);
    kCalculateObcGbsaSoftcoreBornSum(gpu);
    kReduceObcGbsaBornSum(gpu);
    kCalculateCDLJObcGbsaSoftcoreForces1(gpu);

//kPrintForces(gpu, "Post kCalculateCDLJObcGbsaSoftcoreForces1", call );
    if( debug && log ){
        (void) fprintf( log, "\n%s: calling kReduceObcGbsaBornForces\n", methodName.c_str()  );
        (void) fflush( log );
    }

    // compute Born forces

    kReduceObcGbsaSoftcoreBornForces(gpu);

    if( debug && log ){
        (void) fprintf( log, "\n%s calling kCalculateObcGbsaForces2\n", methodName.c_str() );
        (void) fflush( log );
    }

    // second loop of Obc GBSA forces

    kCalculateObcGbsaSoftcoreForces2(gpu);
}

double CudaFreeEnergyCalcGBSAOBCSoftcoreForceKernel::executeEnergy(ContextImpl& context) {
    executeForces( context );
	 return 0.0;
}

CudaFreeEnergyCalcGBVISoftcoreForceKernel::~CudaFreeEnergyCalcGBVISoftcoreForceKernel() {
    if( log ){
        (void) fprintf( log, "CudaFreeEnergyCalcGBVISoftcoreForceKernel destructor called -- freeing gpuGBVISoftcore.\n" );
        (void) fflush( log );
    }
    delete gpuGBVISoftcore;
}

void CudaFreeEnergyCalcGBVISoftcoreForceKernel::initialize(const System& system, const GBVISoftcoreForce& force, const std::vector<double> & inputScaledRadii) {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "CudaFreeEnergyCalcGBVISoftcoreForceKernel::initialize";

// ---------------------------------------------------------------------------------------

    int numParticles = system.getNumParticles();
    _gpuContext* gpu = data.gpu;

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

    for (int i = 0; i < numParticles; i++) {
        double charge, particleRadius, gamma, bornRadiusScaleFactor;
        force.getParticleParameters(i, charge, particleRadius, gamma, bornRadiusScaleFactor);
        particle[i]                  = i;
        radius[i]                    = static_cast<float>( particleRadius );
        gammas[i]                    = static_cast<float>( gamma );
        scaledRadii[i]               = static_cast<float>( inputScaledRadii[i] );
        bornRadiusScaleFactors[i]    = static_cast<float>( bornRadiusScaleFactor );
    }

    // tanh not implemented

//    std::vector<float> tanhScaleFactors;
    std::vector<float> quinticSplineParameters;
    if( force.getBornRadiusScalingMethod() == GBVISoftcoreForce::Tanh ){
/*
        double alpha, beta, gamma;
        force.getTanhParameters( alpha, beta, gamma );
        tanhScaleFactors.resize( 3 );
        tanhScaleFactors[0] = static_cast<float>(alpha);
        tanhScaleFactors[1] = static_cast<float>(beta);
        tanhScaleFactors[2] = static_cast<float>(gamma);
*/
    } else if( force.getBornRadiusScalingMethod() == GBVISoftcoreForce::QuinticSpline ){

        // quintic spline

        quinticSplineParameters.resize(2);
        quinticSplineParameters[0] = static_cast<float>(force.getQuinticLowerLimitFactor());
        quinticSplineParameters[1] = static_cast<float>(force.getQuinticUpperBornRadiusLimit());
        quinticSplineParameters[1] = powf( quinticSplineParameters[1], -3.0f ); 
        setQuinticScaling( 1 );
    }

    // load parameters onto board
    // defined in kCalculateGBVISoftcore.cu

    gpuGBVISoftcore = gpuSetGBVISoftcoreParameters(gpu, static_cast<float>( force.getSoluteDielectric() ), static_cast<float>( force.getSolventDielectric() ),
                                                    particle, radius, gammas, scaledRadii, bornRadiusScaleFactors, quinticSplineParameters);

}

void CudaFreeEnergyCalcGBVISoftcoreForceKernel::executeForces(ContextImpl& context) {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "CudaFreeEnergyCalcGBVISoftcoreForceKernel::executeForces";

// ---------------------------------------------------------------------------------------

    _gpuContext* gpu = data.gpu;
    int debug        = 1;
    static int call  = 0;

    // send address's of arrays, ... to device on first call
    // required since force/energy buffers not set when CudaFreeEnergyCalcGBVISoftcoreForceKernel::initialize() was called

    if( setSim == 0 ){
       setSim++;
       SetCalculateGBVISoftcoreBornSumGpuSim( gpu );
       SetCalculateCDLJObcGbsaSoftcoreGpu1Sim( gpu );
       SetCalculateGBVISoftcoreForces2Sim( gpu );
    }

    // required?

    gpu->bRecalculateBornRadii = true; // fixed

    // calculate Born radii and first loop of GB/VI forces

    if( debug && log ){
        call++;
        if( log ){
            (void) fprintf( log, "\n%s: calling kCalculateCDLJObcGbsaSoftcoreForces1 & %s\n", methodName.c_str(),
                            getQuinticScaling() ? "kReduceGBVIBornSumQuinticScaling" : "kReduceGBVIBornSum" );
            (void) fflush( log );
        }
    }

    kClearBornForces(gpu);
    kCalculateGBVISoftcoreBornSum(gpu);

    if( getQuinticScaling() ){
        kReduceGBVIBornSumQuinticScaling(gpu, gpuGBVISoftcore );
    } else {
        kReduceGBVIBornSum(gpu);
    }
    kCalculateCDLJObcGbsaSoftcoreForces1(gpu);

//kPrintForces(gpu, "Post kCalculateCDLJObcGbsaSoftcoreForces1", call );

    if( debug && log ){
        (void) fprintf( log, "\n%s: calling %s\n", methodName.c_str(),
                        getQuinticScaling() ? "kReduceGBVIBornForcesQuinticScaling" : "kReduceObcGbsaBornForces" );
        (void) fflush( log );
    }

    // compute Born forces

    if( getQuinticScaling() ){
        kReduceGBVIBornForcesQuinticScaling(gpu);
    } else {
        gpu->bIncludeGBVI = true;
        kReduceObcGbsaBornForces(gpu);
        gpu->bIncludeGBVI = false;
    }

    if( debug && log ){
        (void) fprintf( log, "\n%s: calling kCalculateGBVIForces2\n", methodName.c_str() );
        (void) fflush( log );
    }

    // second loop of GB/VI forces

    kCalculateGBVISoftcoreForces2(gpu);
}

double CudaFreeEnergyCalcGBVISoftcoreForceKernel::executeEnergy(ContextImpl& context) {
    executeForces( context );
    return 0.0;
}

int CudaFreeEnergyCalcGBVISoftcoreForceKernel::getQuinticScaling( void ) const {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "CudaFreeEnergyCalcGBVISoftcoreForceKernel::getQuinticScaling";

// ---------------------------------------------------------------------------------------

    return quinticScaling;
}

void CudaFreeEnergyCalcGBVISoftcoreForceKernel::setQuinticScaling( int inputQuinticScaling) {

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "CudaFreeEnergyCalcGBVISoftcoreForceKernel::setQuinticScaling";

// ---------------------------------------------------------------------------------------

    quinticScaling = inputQuinticScaling;
}

