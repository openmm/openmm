/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs                                                   *
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

#include "ValidateOpenMM.h" 

using namespace OpenMM;

// fixed force names
    
const std::string ValidateOpenMM::HARMONIC_BOND_FORCE             = "HarmonicBond";
const std::string ValidateOpenMM::HARMONIC_ANGLE_FORCE            = "HarmonicAngle";
const std::string ValidateOpenMM::PERIODIC_TORSION_FORCE          = "PeriodicTorsion";
const std::string ValidateOpenMM::RB_TORSION_FORCE                = "RbTorsion";

const std::string ValidateOpenMM::NB_FORCE                        = "Nb";
const std::string ValidateOpenMM::NB_SOFTCORE_FORCE               = "NbSoftcore";

const std::string ValidateOpenMM::NB_EXCEPTION_FORCE              = "NbException";
const std::string ValidateOpenMM::NB_EXCEPTION_SOFTCORE_FORCE     = "NbSoftcoreException";

const std::string ValidateOpenMM::GBSA_OBC_FORCE                  = "GBSAOBC";
const std::string ValidateOpenMM::GBSA_OBC_SOFTCORE_FORCE         = "GBSAOBCSoftcore";

const std::string ValidateOpenMM::GBVI_FORCE                      = "GBVI";
const std::string ValidateOpenMM::GBVI_SOFTCORE_FORCE             = "GBVISoftcore";

const std::string ValidateOpenMM::CM_MOTION_REMOVER               = "CMMotionRemover";
const std::string ValidateOpenMM::ANDERSEN_THERMOSTAT             = "AndersenThermostat";
const std::string ValidateOpenMM::CUSTOM_BOND_FORCE               = "CustomBond";
const std::string ValidateOpenMM::CUSTOM_EXTERNAL_FORCE           = "CustomExternal";
const std::string ValidateOpenMM::CUSTOM_NONBONDED_FORCE          = "CustomNonBonded";


ValidateOpenMM::ValidateOpenMM( void ) {

    _log                                          = NULL;

    // force dependencies 
    // these may need to specialized depending on platform since
    // currently apply to Cuda platform, but may not apply to OpenCL platform

    std::vector< std::string > nbVector;
    nbVector.push_back( NB_FORCE );
    _forceDependencies[GBSA_OBC_FORCE]            = nbVector;
    _forceDependencies[GBVI_FORCE]                = nbVector;

    std::vector< std::string > nbSoftcoreVector;
    nbSoftcoreVector.push_back( NB_SOFTCORE_FORCE );
    _forceDependencies[GBSA_OBC_SOFTCORE_FORCE]   = nbSoftcoreVector;
    _forceDependencies[GBVI_SOFTCORE_FORCE]       = nbSoftcoreVector;
    
}

ValidateOpenMM::~ValidateOpenMM() {
}


int ValidateOpenMM::isNanOrInfinity( double number ){
    return (number != number || number == std::numeric_limits<double>::infinity() || number == -std::numeric_limits<double>::infinity()) ? 1 : 0;
}

FILE* ValidateOpenMM::getLog( void ) const {
    return _log;
}

void ValidateOpenMM::setLog(  FILE* log ){
    _log = log;
}

std::string ValidateOpenMM::getForceName(const Force& force) const { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ValidateForce::getForceName";

// ---------------------------------------------------------------------------------------

    std::string forceName  = "NA";

    // bond force

    try {
        const HarmonicBondForce& harmonicBondForce = dynamic_cast<const HarmonicBondForce&>(force);
        return HARMONIC_BOND_FORCE;
    } catch( std::bad_cast ){
    }

    // angle force

    try {
        const HarmonicAngleForce& harmonicAngleForce = dynamic_cast<const HarmonicAngleForce&>(force);
        return HARMONIC_ANGLE_FORCE;
    } catch( std::bad_cast ){
    }

    // periodic torsion force

    try {
        const PeriodicTorsionForce & periodicTorsionForce = dynamic_cast<const PeriodicTorsionForce&>(force);
       return PERIODIC_TORSION_FORCE;
    } catch( std::bad_cast ){
    }

    // RB torsion force

    try {
        const RBTorsionForce& rBTorsionForce = dynamic_cast<const RBTorsionForce&>(force);
       return RB_TORSION_FORCE;
    } catch( std::bad_cast ){
    }

    // nonbonded force

    try {
        const NonbondedForce& nbForce = dynamic_cast<const NonbondedForce&>(force);
        return  NB_FORCE;
    } catch( std::bad_cast ){
    }

    // GBSA OBC

    try {
       const GBSAOBCForce& obcForce       = dynamic_cast<const GBSAOBCForce&>(force);
       return GBSA_OBC_FORCE;
    } catch( std::bad_cast ){
    }

    // GB/VI

    try {
        const GBVIForce& obcForce  = dynamic_cast<const GBVIForce&>(force);
        return GBVI_FORCE;
    } catch( std::bad_cast ){
    }

#ifdef INCLUDE_FREE_ENERGY_PLUGIN

    // free energy plugin forces

    // nonbonded softcore

    try {
        const NonbondedSoftcoreForce& nbForce = dynamic_cast<const NonbondedSoftcoreForce&>(force);
        return NB_SOFTCORE_FORCE;
    } catch( std::bad_cast ){
    }

    // GBSA OBC softcore

    try {
        const GBSAOBCSoftcoreForce& obcForce = dynamic_cast<const GBSAOBCSoftcoreForce&>(force);
        return GBSA_OBC_SOFTCORE_FORCE;
    } catch( std::bad_cast ){
    }

    // GB/VI softcore

    try {
        const GBVISoftcoreForce& gbviForce = dynamic_cast<const GBVISoftcoreForce&>(force);
        return  GBVI_SOFTCORE_FORCE;
    } catch( std::bad_cast ){
    }

#endif

    // CMMotionRemover

    try {
        const CMMotionRemover& cMMotionRemover  = dynamic_cast<const CMMotionRemover&>(force);
        return CM_MOTION_REMOVER;
    } catch( std::bad_cast ){
    }

    // AndersenThermostat

    try {
        const AndersenThermostat & andersenThermostat = dynamic_cast<const AndersenThermostat&>(force);
        return ANDERSEN_THERMOSTAT;
    } catch( std::bad_cast ){
    }

    // CustomBondForce

    try {
        const CustomBondForce & customBondForce = dynamic_cast<const CustomBondForce&>(force);
        return CUSTOM_BOND_FORCE;
    } catch( std::bad_cast ){
    }

    // CustomExternalForce

    try {
        const CustomExternalForce& customExternalForce = dynamic_cast<const CustomExternalForce&>(force);
        return CUSTOM_EXTERNAL_FORCE;
    } catch( std::bad_cast ){
    }

    // CustomNonbondedForce

    try {
        const CustomNonbondedForce& customNonbondedForce = dynamic_cast<const CustomNonbondedForce&>(force);
        return CUSTOM_NONBONDED_FORCE;
    } catch( std::bad_cast ){
    }

    return forceName;
}

Force* ValidateOpenMM::copyForce(const Force& force) const { 

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ValidateForce::copyForce";

// ---------------------------------------------------------------------------------------

    // bond force

    try {
        const HarmonicBondForce& harmonicBondForce = dynamic_cast<const HarmonicBondForce&>(force);
        return new HarmonicBondForce( harmonicBondForce );
    } catch( std::bad_cast ){
    }

    // angle force

    try {
        const HarmonicAngleForce& harmonicAngleForce = dynamic_cast<const HarmonicAngleForce&>(force);
        return new HarmonicAngleForce( harmonicAngleForce );
    } catch( std::bad_cast ){
    }

    // periodic torsion force

    try {
        const PeriodicTorsionForce & periodicTorsionForce = dynamic_cast<const PeriodicTorsionForce&>(force);
        return new PeriodicTorsionForce( periodicTorsionForce );
    } catch( std::bad_cast ){
    }

    // RB torsion force

    try {
        const RBTorsionForce& rBTorsionForce = dynamic_cast<const RBTorsionForce&>(force);
        return new RBTorsionForce( rBTorsionForce );
    } catch( std::bad_cast ){
    }

    // nonbonded force

    try {
        const NonbondedForce& nbForce = dynamic_cast<const NonbondedForce&>(force);
        return new NonbondedForce( nbForce ); 
    } catch( std::bad_cast ){
    }

    // GBSA OBC

    try {
        const GBSAOBCForce& obcForce       = dynamic_cast<const GBSAOBCForce&>(force);
        return new GBSAOBCForce( obcForce );
    } catch( std::bad_cast ){
    }

    // GB/VI

    try {
        const GBVIForce& gbviForce  = dynamic_cast<const GBVIForce&>(force);
        return new GBVIForce( gbviForce );
    } catch( std::bad_cast ){
    }

#ifdef INCLUDE_FREE_ENERGY_PLUGIN

    // free energy plugin forces

    // nonbonded softcore
    try {
         const NonbondedSoftcoreForce& nbForce = dynamic_cast<const NonbondedSoftcoreForce&>(force);
         return new NonbondedSoftcoreForce( nbForce );
    } catch( std::bad_cast ){
    }

    // GBSA OBC softcore

    try {
        const GBSAOBCSoftcoreForce& obcForce = dynamic_cast<const GBSAOBCSoftcoreForce&>(force);
        return new GBSAOBCSoftcoreForce( obcForce );
    } catch( std::bad_cast ){
    }

    // GB/VI softcore

    try {
        const GBVISoftcoreForce& gbviForce = dynamic_cast<const GBVISoftcoreForce&>(force);
        return new GBVISoftcoreForce( gbviForce );
    } catch( std::bad_cast ){
    }
#endif

    // CMMotionRemover

    try {
        const CMMotionRemover& cMMotionRemover  = dynamic_cast<const CMMotionRemover&>(force);
        return new CMMotionRemover( cMMotionRemover );
    } catch( std::bad_cast ){
    }

    // AndersenThermostat

    try {
        const AndersenThermostat & andersenThermostat = dynamic_cast<const AndersenThermostat&>(force);
        return new AndersenThermostat( andersenThermostat );
    } catch( std::bad_cast ){
    }

    // CustomBondForce

    try {
        const CustomBondForce & customBondForce = dynamic_cast<const CustomBondForce&>(force);
        return new CustomBondForce( customBondForce );
    } catch( std::bad_cast ){
    }

    // CustomExternalForce

    try {
        const CustomExternalForce& customExternalForce = dynamic_cast<const CustomExternalForce&>(force);
        return new CustomExternalForce( customExternalForce );
    } catch( std::bad_cast ){
    }

    // CustomNonbondedForce

    try {
        const CustomNonbondedForce& customNonbondedForce = dynamic_cast<const CustomNonbondedForce&>(force);
        return new CustomNonbondedForce( customNonbondedForce );
    } catch( std::bad_cast ){
    }

    return NULL;
}

/** 
 * Get copy of input system, but leave out forces
 * 
 * @param systemToCopy   system to copy
 * @return system  w/o forces
 */
System* ValidateOpenMM::copySystemExcludingForces( const System& systemToCopy ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ValidateOpenMM::copySystemExcludingForces";

// ---------------------------------------------------------------------------------------

    // create new system and add masses and box dimensions

    System* system = new System();
    for(int ii = 0; ii < systemToCopy.getNumParticles(); ii++ ){
        double mass = systemToCopy.getParticleMass( ii );
        system->addParticle( mass );
    }   

    Vec3 a, b, c;
    systemToCopy.getDefaultPeriodicBoxVectors( a, b, c );
    system->setDefaultPeriodicBoxVectors( a, b, c );

    copyConstraints( systemToCopy, system );

    return system;
}

/**---------------------------------------------------------------------------------------

   Set the velocities/positions of context2 to those of context1

   @param context1                 context1 
   @param context2                 context2 

   @return 0

   --------------------------------------------------------------------------------------- */

void ValidateOpenMM::synchContexts( const Context& context1, Context& context2 ) const {

// ---------------------------------------------------------------------------------------

    //static const char* methodName  = "ValidateOpenMM::synchContexts: ";

// ---------------------------------------------------------------------------------------

    const State state                       = context1.getState(State::Positions | State::Velocities);
    const std::vector<Vec3>& positions      = state.getPositions();
    const std::vector<Vec3>& velocities     = state.getVelocities();
 
    context2.setPositions( positions );
    context2.setVelocities( velocities );
 
    return;
}

/**---------------------------------------------------------------------------------------

   Copy constraints

   @param systemToCopy         system whose constraints are to be copied
   @param system               system to add constraints to
   @param log                  log file pointer -- may be NULL

   --------------------------------------------------------------------------------------- */

void ValidateOpenMM::copyConstraints( const System& systemToCopy, System* system, FILE* log ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "copyConstraints";
   
// ---------------------------------------------------------------------------------------

    for( int ii = 0; ii < systemToCopy.getNumConstraints(); ii++ ){
        int particle1, particle2;
        double distance;
        systemToCopy.getConstraintParameters( ii, particle1, particle2, distance ); 
        system->addConstraint( particle1, particle2, distance );
    }

}

/**---------------------------------------------------------------------------------------
    
   Get force dependencies
    
   @param forceName            force to check if there exist any dependencies
   @param returnVector         vector of forces the input force is dependent on (example: GBSAOBC force requires Nonbonded force)
    
   --------------------------------------------------------------------------------------- */
    
void ValidateOpenMM::getForceDependencies( std::string forceName, StringVector& returnVector ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "getForceDependencies";
   
// ---------------------------------------------------------------------------------------

    StringStringVectorMapCI dependencyVector = _forceDependencies.find( forceName );
    if( dependencyVector != _forceDependencies.end() ){
        StringVector dependices = (*dependencyVector).second;
        for( unsigned int ii = 0; ii < dependices.size(); ii++ ){
            returnVector.push_back( dependices[ii] );
        }
    }
}

/**
 * Write masses and box dimensions to parameter file
 */

void ValidateOpenMM::writeMasses( FILE* filePtr, const System& system ) const {

    (void) fprintf( filePtr, "Masses %d\n", system.getNumParticles() );
    for(int ii = 0; ii < system.getNumParticles(); ii++ ){
       double mass = system.getParticleMass( ii );
       (void) fprintf( filePtr, "%8d %14.7e\n", ii, mass );
    }

    Vec3 a, b, c;
    system.getDefaultPeriodicBoxVectors( a, b, c);
    (void) fprintf( filePtr, "Box %14.6f %14.6f %14.6f  %14.6f %14.6f %14.6f   %14.6f %14.6f %14.6f\n",
                    a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2] );
}

/**
 * Write constraints to parameter file
 */

void ValidateOpenMM::writeConstraints( FILE* filePtr, const System& system ) const {

    (void) fprintf( filePtr, "Constraints %d\n", system.getNumConstraints() );
    for(int ii = 0; ii < system.getNumConstraints(); ii++ ){
       int particle1, particle2;
       double distance;
       system.getConstraintParameters( ii, particle1, particle2, distance );
       (void) fprintf( filePtr, "%8d %8d %8d %14.7e\n", ii, particle1, particle2, distance );
    }
}

/**
 * Write harmonicBondForce parameters to file
 */

void ValidateOpenMM::writeHarmonicBondForce( FILE* filePtr, const HarmonicBondForce& harmonicBondForce ) const {

    (void) fprintf( filePtr, "HarmonicBondForce %d\n", harmonicBondForce.getNumBonds() );
    for(int ii = 0; ii < harmonicBondForce.getNumBonds(); ii++ ){
       int particle1, particle2;
       double length, k;
       harmonicBondForce.getBondParameters( ii, particle1, particle2, length, k );
       (void) fprintf( filePtr, "%8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, length, k);
    }
}

/**
 * Write harmonicAngleForce parameters to file
 */

void ValidateOpenMM::writeHarmonicAngleForce( FILE* filePtr, const HarmonicAngleForce& harmonicAngleForce ) const {

    (void) fprintf( filePtr, "HarmonicAngleForce %d\n", harmonicAngleForce.getNumAngles() );
    for(int ii = 0; ii < harmonicAngleForce.getNumAngles(); ii++ ){
       int particle1, particle2, particle3;
       double angle, k;
       harmonicAngleForce.getAngleParameters( ii, particle1, particle2, particle3, angle, k );
       (void) fprintf( filePtr, "%8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, angle, k);
    }
}

/**
 * Write rbTorsionForce parameters to file
 */

void ValidateOpenMM::writeRbTorsionForce( FILE* filePtr, const RBTorsionForce& rbTorsionForce ) const {

    (void) fprintf( filePtr, "RBTorsionForce %d\n", rbTorsionForce.getNumTorsions() );
    for(int ii = 0; ii < rbTorsionForce.getNumTorsions(); ii++ ){
       int particle1, particle2, particle3, particle4;
       double c0, c1, c2, c3, c4, c5;
       rbTorsionForce.getTorsionParameters( ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
       (void) fprintf( filePtr, "%8d %8d %8d %8d %8d %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e\n", ii, particle1, particle2, particle3, particle4, c0, c1, c2, c3, c4, c5 );
    }
}

/**
 * Write periodicTorsionForce parameters to file
 */

void ValidateOpenMM::writePeriodicTorsionForce( FILE* filePtr, const PeriodicTorsionForce& periodicTorsionForce ) const {

    (void) fprintf( filePtr, "PeriodicTorsionForce %d\n", periodicTorsionForce.getNumTorsions() );
    for(int ii = 0; ii < periodicTorsionForce.getNumTorsions(); ii++ ){
       int particle1, particle2, particle3, particle4, periodicity;
       double phase, k;
       periodicTorsionForce.getTorsionParameters( ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
       (void) fprintf( filePtr, "%8d %8d %8d %8d %8d %8d %14.7e %14.7e\n", ii, particle1, particle2, particle3, particle4, periodicity, phase, k );
    }
}

/**
 * Write GBSAOBCForce parameters to file
 */

void ValidateOpenMM::writeGbsaObcForce( FILE* filePtr, const GBSAOBCForce& gbsaObcForce ) const {

    (void) fprintf( filePtr, "GBSAOBCForce %d\n", gbsaObcForce.getNumParticles() );
    for(int ii = 0; ii <  gbsaObcForce.getNumParticles(); ii++ ){
       double charge, radius, scalingFactor;
       gbsaObcForce.getParticleParameters( ii, charge, radius, scalingFactor );
       (void) fprintf( filePtr, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor );
    }
    (void) fprintf( filePtr, "SoluteDielectric %14.7e\n", gbsaObcForce.getSoluteDielectric() );
    (void) fprintf( filePtr, "SolventDielectric %14.7e\n", gbsaObcForce.getSolventDielectric() );
}

#ifdef INCLUDE_FREE_ENERGY_PLUGIN

/**
 * Write GBSAOBCSoftcoreForce parameters to file
 */

void ValidateOpenMM::writeGbsaObcSoftcoreForce( FILE* filePtr, const GBSAOBCSoftcoreForce& gbsaObcForce ) const {

    (void) fprintf( filePtr, "GBSAOBCSoftcoreForce %d\n", gbsaObcForce.getNumParticles() );
    for(int ii = 0; ii <  gbsaObcForce.getNumParticles(); ii++ ){
       double charge, radius, scalingFactor;
       double nonPolarScalingFactor;
       gbsaObcForce.getParticleParameters( ii, charge, radius, scalingFactor, nonPolarScalingFactor );
       (void) fprintf( filePtr, "%8d  %14.7e %14.7e %14.7e %14.7e\n", ii, charge, radius, scalingFactor, nonPolarScalingFactor );
    }
    (void) fprintf( filePtr, "SoluteDielectric %14.7e\n", gbsaObcForce.getSoluteDielectric() );
    (void) fprintf( filePtr, "SolventDielectric %14.7e\n", gbsaObcForce.getSolventDielectric() );
}

#endif

/**
 * Write GBSA GB/VI Force parameters to file
 */

void ValidateOpenMM::writeGBVIForce( FILE* filePtr, const GBVIForce& gbviForce ) const {

    (void) fprintf( filePtr, "GBVIForce %d\n", gbviForce.getNumParticles() );
    for(int ii = 0; ii <  gbviForce.getNumParticles(); ii++ ){
       double charge, radius, gamma;
       gbviForce.getParticleParameters( ii, charge, radius, gamma );
       (void) fprintf( filePtr, "%8d  %14.7e %14.7e %14.7e\n", ii, charge, radius, gamma );
    }

    (void) fprintf( filePtr, "GBVIBonds  %d\n", gbviForce.getNumBonds() );
    for(int ii = 0; ii <  gbviForce.getNumBonds(); ii++ ){
       int atomI, atomJ;
       double bondLength;
       gbviForce.getBondParameters( ii, atomI, atomJ, bondLength );
       (void) fprintf( filePtr, "%8d %8d %8d %14.7e\n", ii, atomI, atomJ, bondLength );
    }
    (void) fprintf( filePtr, "SoluteDielectric %14.7e\n", gbviForce.getSoluteDielectric() );
    (void) fprintf( filePtr, "SolventDielectric %14.7e\n", gbviForce.getSolventDielectric() );

}

#ifdef INCLUDE_FREE_ENERGY_PLUGIN

/**
 * Write GBSA GB/VI softcore force parameters to file
 */

void ValidateOpenMM::writeGBVISoftcoreForce( FILE* filePtr, const GBVISoftcoreForce& gbviSoftcoreForce ) const {

    (void) fprintf( filePtr, "GBVISoftcoreForce %d\n", gbviSoftcoreForce.getNumParticles() );
    for(int ii = 0; ii <  gbviSoftcoreForce.getNumParticles(); ii++ ){
       double charge, radius, gamma;
       double nonPolarScalingFactor;
       gbviSoftcoreForce.getParticleParameters( ii, charge, radius, gamma, nonPolarScalingFactor );
       (void) fprintf( filePtr, "%8d  %14.7e %14.7e %14.7e %14.7e\n", ii, charge, radius, gamma, nonPolarScalingFactor );
    }

    (void) fprintf( filePtr, "GBVISoftcoreBonds  %d\n", gbviSoftcoreForce.getNumBonds() );
    for(int ii = 0; ii <  gbviSoftcoreForce.getNumBonds(); ii++ ){
       int atomI, atomJ;
       double bondLength;
       gbviSoftcoreForce.getBondParameters( ii, atomI, atomJ, bondLength );
       (void) fprintf( filePtr, "%8d %8d %8d %14.7e\n", ii, atomI, atomJ, bondLength );
    }
    (void) fprintf( filePtr, "SoluteDielectric %14.7e\n", gbviSoftcoreForce.getSoluteDielectric() );
    (void) fprintf( filePtr, "SolventDielectric %14.7e\n", gbviSoftcoreForce.getSolventDielectric() );

    (void) fprintf( filePtr, "BornRadiusScalingMethod %d\n", gbviSoftcoreForce.getBornRadiusScalingMethod() );
    (void) fprintf( filePtr, "QuinticLowerLimitFactor %14.7e\n", gbviSoftcoreForce.getQuinticLowerLimitFactor() );
    (void) fprintf( filePtr, "QuinticUpperBornRadiusLimit %14.7e\n", gbviSoftcoreForce.getQuinticUpperBornRadiusLimit() );
}

#endif

/**
 * Write coordinates to file
 */

void ValidateOpenMM::writeVec3( FILE* filePtr, const std::vector<Vec3>& vect3Array ) const {

    for( unsigned int ii = 0; ii < vect3Array.size(); ii++ ){
       (void) fprintf( filePtr, "%8d  %14.7e %14.7e %14.7e\n", ii, 
                       vect3Array[ii][0], vect3Array[ii][1], vect3Array[ii][2] );
    }
}

/**
 * Write context info to file
 */

void ValidateOpenMM::writeContext( FILE* filePtr, const Context& context ) const {

    State state                  = context.getState( State::Positions | State::Velocities | State::Forces | State::Energy );
    std::vector<Vec3> positions  = state.getPositions();
    std::vector<Vec3> velocities = state.getVelocities();
    std::vector<Vec3> forces     = state.getForces();

    (void) fprintf( filePtr, "Positions %zu\n", positions.size() );
    writeVec3( filePtr, positions );

    (void) fprintf( filePtr, "Velocities %zu\n", velocities.size() );
    writeVec3( filePtr, velocities );

    (void) fprintf( filePtr, "Forces %zu\n", forces.size() );
    writeVec3( filePtr, forces );

    (void) fprintf( filePtr, "KineticEnergy %14.7e\n", state.getKineticEnergy() );
    (void) fprintf( filePtr, "PotentialEnergy %14.7e\n", state.getPotentialEnergy() );
}

/**
 * Write nonbonded parameters to file
 */

void ValidateOpenMM::writeNonbondedForce( FILE* filePtr, const NonbondedForce & nonbondedForce ) const {

    // charge and vdw parameters

    (void) fprintf( filePtr, "NonbondedForce %d\n", nonbondedForce.getNumParticles() );
    for(int ii = 0; ii < nonbondedForce.getNumParticles(); ii++ ){
       double charge, sigma, epsilon;
       nonbondedForce.getParticleParameters( ii, charge, sigma, epsilon );
       (void) fprintf( filePtr, "%8d %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon );
    }

    // cutoff, dielectric, Ewald tolerance

    (void) fprintf( filePtr, "CutoffDistance %14.7e\n", nonbondedForce.getCutoffDistance() );
    (void) fprintf( filePtr, "RFDielectric %14.7e\n", nonbondedForce.getReactionFieldDielectric() );
    (void) fprintf( filePtr, "EwaldRTolerance %14.7e\n", nonbondedForce.getEwaldErrorTolerance() );

    // cutoff mode

    std::string nonbondedForceMethod;
    switch( nonbondedForce.getNonbondedMethod() ){
        case NonbondedForce::NoCutoff:
            nonbondedForceMethod = "NoCutoff";
            break;
        case NonbondedForce::CutoffNonPeriodic:
            nonbondedForceMethod = "CutoffNonPeriodic";
            break;
        case NonbondedForce::CutoffPeriodic:
            nonbondedForceMethod = "CutoffPeriodic";
            break;
        case NonbondedForce::Ewald:
            nonbondedForceMethod = "Ewald";
            break;
        case NonbondedForce::PME:
            nonbondedForceMethod = "PME";
            break;
        default:
            nonbondedForceMethod = "Unknown";
    }
    (void) fprintf( filePtr, "NonbondedForceMethod %s\n", nonbondedForceMethod.c_str() );

    (void) fprintf( filePtr, "NonbondedForceExceptions %d\n", nonbondedForce.getNumExceptions() );
    for(int ii = 0; ii < nonbondedForce.getNumExceptions(); ii++ ){
       int particle1, particle2;
       double chargeProd, sigma, epsilon;
       nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon );
       (void) fprintf( filePtr, "%8d %8d %8d %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon );
    }

}

#ifdef INCLUDE_FREE_ENERGY_PLUGIN

/**
 * Write nonbonded softcore parameters to file
 */

void ValidateOpenMM::writeNonbondedSoftcoreForce( FILE* filePtr, const NonbondedSoftcoreForce & nonbondedSoftcoreForce ) const {

    (void) fprintf( filePtr, "NonbondedSoftcoreForce %d\n", nonbondedSoftcoreForce.getNumParticles() );
    for(int ii = 0; ii < nonbondedSoftcoreForce.getNumParticles(); ii++ ){
       double charge, sigma, epsilon, softCoreLJ;
       nonbondedSoftcoreForce.getParticleParameters( ii, charge, sigma, epsilon, softCoreLJ );
       (void) fprintf( filePtr, "%8d %14.7e %14.7e %14.7e %14.7e\n", ii, charge, sigma, epsilon, softCoreLJ );
    }

    (void) fprintf( filePtr, "CutoffDistance %14.7e\n", nonbondedSoftcoreForce.getCutoffDistance() );
    (void) fprintf( filePtr, "RFDielectric %14.7e\n", nonbondedSoftcoreForce.getReactionFieldDielectric() );
    (void) fprintf( filePtr, "EwaldRTolerance %14.7e\n", nonbondedSoftcoreForce.getEwaldErrorTolerance() );

    std::string nonbondedSoftcoreForceMethod;
    switch( nonbondedSoftcoreForce.getNonbondedMethod() ){
        case NonbondedSoftcoreForce::NoCutoff:
            nonbondedSoftcoreForceMethod = "NoCutoff";
            break;
/*
        case NonbondedSoftcoreForce::CutoffNonPeriodic:
            nonbondedSoftcoreForceMethod = "CutoffNonPeriodic";
            break;
        case NonbondedSoftcoreForce::CutoffPeriodic:
            nonbondedSoftcoreForceMethod = "CutoffPeriodic";
            break;
        case NonbondedSoftcoreForce::Ewald:
            nonbondedSoftcoreForceMethod = "Ewald";
            break;
        case NonbondedSoftcoreForce::PME:
            nonbondedSoftcoreForceMethod = "PME";
            break;
*/
        default:
            nonbondedSoftcoreForceMethod = "Unknown";
    }
    (void) fprintf( filePtr, "NonbondedSoftcoreForceMethod %s\n", nonbondedSoftcoreForceMethod.c_str() );

    (void) fprintf( filePtr, "NonbondedSoftcoreForceExceptions %d\n", nonbondedSoftcoreForce.getNumExceptions() );
    for(int ii = 0; ii < nonbondedSoftcoreForce.getNumExceptions(); ii++ ){
       int particle1, particle2;
       double chargeProd, sigma, epsilon, softCoreLJ;
       nonbondedSoftcoreForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softCoreLJ );
       (void) fprintf( filePtr, "%8d %8d %8d %14.7e %14.7e %14.7e %14.7e\n", ii, particle1, particle2, chargeProd, sigma, epsilon, softCoreLJ );
    }

    return;
}
#endif

void ValidateOpenMM::writeIntegrator( FILE* filePtr, const Integrator& integrator ) const {

    // LangevinIntegrator
 
    try {
        const LangevinIntegrator langevinIntegrator = dynamic_cast<const LangevinIntegrator&>(integrator);
        (void) fprintf( filePtr, "Integrator LangevinIntegrator\n" );
        (void) fprintf( filePtr, "StepSize %14.7e\n", langevinIntegrator.getStepSize() );
        (void) fprintf( filePtr, "ConstraintTolerance %14.7e\n", langevinIntegrator.getConstraintTolerance() );
        (void) fprintf( filePtr, "Temperature %14.7e\n", langevinIntegrator.getTemperature() );
        (void) fprintf( filePtr, "Friction %14.7e\n", langevinIntegrator.getFriction() );
        (void) fprintf( filePtr, "RandomNumberSeed %d\n", langevinIntegrator.getRandomNumberSeed() );
        return;
    } catch( std::bad_cast ){
    }
 
    // VariableLangevinIntegrator
 
    try {
        const VariableLangevinIntegrator& langevinIntegrator = dynamic_cast<const VariableLangevinIntegrator&>(integrator);
        (void) fprintf( filePtr, "Integrator VariableLangevinIntegrator\n" );
        (void) fprintf( filePtr, "StepSize %14.7e\n", langevinIntegrator.getStepSize() );
        (void) fprintf( filePtr, "ConstraintTolerance %14.7e\n", langevinIntegrator.getConstraintTolerance() );
        (void) fprintf( filePtr, "Temperature %14.7e\n", langevinIntegrator.getTemperature() );
        (void) fprintf( filePtr, "Friction %14.7e\n", langevinIntegrator.getFriction() );
        (void) fprintf( filePtr, "RandomNumberSeed %d\n", langevinIntegrator.getRandomNumberSeed() );
        (void) fprintf( filePtr, "ErrorTolerance %14.7e\n", langevinIntegrator.getErrorTolerance() );
        return;
    } catch( std::bad_cast ){
    }
 
    // VerletIntegrator
 
    try {
        const VerletIntegrator& verletIntegrator = dynamic_cast<const VerletIntegrator&>(integrator);
        (void) fprintf( filePtr, "Integrator VerletIntegrator\n" );
        (void) fprintf( filePtr, "StepSize %14.7e\n", verletIntegrator.getStepSize() );
        (void) fprintf( filePtr, "ConstraintTolerance %14.7e\n", verletIntegrator.getConstraintTolerance() );
        return;
    } catch( std::bad_cast ){
    }
     
    // VariableVerletIntegrator
 
    try {
        const VariableVerletIntegrator & variableVerletIntegrator = dynamic_cast<const VariableVerletIntegrator&>(integrator);
        (void) fprintf( filePtr, "Integrator VariableVerletIntegrator\n" );
        (void) fprintf( filePtr, "StepSize %14.7e\n", variableVerletIntegrator.getStepSize() );
        (void) fprintf( filePtr, "ConstraintTolerance %14.7e\n", variableVerletIntegrator.getConstraintTolerance() );
        (void) fprintf( filePtr, "ErrorTolerance %14.7e\n", variableVerletIntegrator.getErrorTolerance() );
        return;
    } catch( std::bad_cast ){
    }
     
    // BrownianIntegrator
 
    try {
        const BrownianIntegrator& brownianIntegrator = dynamic_cast<const BrownianIntegrator&>(integrator);
        (void) fprintf( filePtr, "Integrator BrownianIntegrator\n" );
        (void) fprintf( filePtr, "StepSize %14.7e\n", brownianIntegrator.getStepSize() );
        (void) fprintf( filePtr, "ConstraintTolerance %14.7e\n", brownianIntegrator.getConstraintTolerance() );
        (void) fprintf( filePtr, "Temperature %14.7e\n", brownianIntegrator.getTemperature() );
        (void) fprintf( filePtr, "Friction %14.7e\n", brownianIntegrator.getFriction() );
        (void) fprintf( filePtr, "RandomNumberSeed %d\n", brownianIntegrator.getRandomNumberSeed() );
        return;
    } catch( std::bad_cast ){
    }
     
    if( getLog() ){
       (void) fprintf( getLog(), "Integrator not recognized." );
       (void) fflush( getLog() );
    }

    return;
}

/**
 * Write parameter file: Custom forces not implemented
 * Mesage is sent to stderr if a force is not recognized
 *
 */

void ValidateOpenMM::writeParameterFile( const Context& context, const std::string& parameterFileName ) const {

    // open file
 
    FILE* filePtr = fopen( parameterFileName.c_str(), "w" );

    const System& system          = context.getSystem();
    const Integrator& integrator  = context.getIntegrator();

    // (void) fprintf( filePtr, "Version %s\n", versionString.c_str() );
    (void) fprintf( filePtr, "Particles %8d\n", system.getNumParticles() );
    writeMasses( filePtr, system );
    (void) fprintf( filePtr, "NumberOfForces    %8d\n", system.getNumForces() );

    // print active forces and relevant parameters

    for(int i = 0; i < system.getNumForces(); ++i) {

        int hit                       = 0;
        const Force& force            = system.getForce(i);

        // bond

        if( !hit ){

            try {
               const HarmonicBondForce& harmonicBondForce = dynamic_cast<const HarmonicBondForce&>(force);
               writeHarmonicBondForce( filePtr, harmonicBondForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // angle

        if( !hit ){
    
            try {
               const HarmonicAngleForce& harmonicAngleForce = dynamic_cast<const HarmonicAngleForce&>(force);
               writeHarmonicAngleForce( filePtr, harmonicAngleForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }

        // PeriodicTorsionForce
    
        if( !hit ){
    
            try {
               const PeriodicTorsionForce & periodicTorsionForce = dynamic_cast<const PeriodicTorsionForce&>(force);
               writePeriodicTorsionForce( filePtr, periodicTorsionForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // RBTorsionForce
    
        if( !hit ){
            try {
               const RBTorsionForce& rBTorsionForce = dynamic_cast<const RBTorsionForce&>(force);
               writeRbTorsionForce(  filePtr, rBTorsionForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // nonbonded
    
        if( !hit ){
            try {
               const NonbondedForce& nbForce = dynamic_cast<const NonbondedForce&>(force);
               writeNonbondedForce( filePtr, nbForce );
               hit++;
            } catch( std::bad_cast ){
            }
        } 

        // nonbonded softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( !hit ){
            try {
               const NonbondedSoftcoreForce& nbForce = dynamic_cast<const NonbondedSoftcoreForce&>(force);
               writeNonbondedSoftcoreForce( filePtr, nbForce );
               hit++;
            } catch( std::bad_cast ){
            }
        } 
#endif

        // GBSA OBC
    
        if( !hit ){
            try {
               const GBSAOBCForce& obcForce = dynamic_cast<const GBSAOBCForce&>(force);
               writeGbsaObcForce(  filePtr, obcForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    

        // GBSA OBC softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( !hit ){
            try {
               const GBSAOBCSoftcoreForce& obcForce = dynamic_cast<const GBSAOBCSoftcoreForce&>(force);
               writeGbsaObcSoftcoreForce(  filePtr, obcForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
#endif
    
        // GB/VI
    
        if( !hit ){
            try {
               const GBVIForce& obcForce = dynamic_cast<const GBVIForce&>(force);
               writeGBVIForce(  filePtr, obcForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
    
        // GB/VI softcore
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
        if( !hit ){
            try {
               const GBVISoftcoreForce& obcForce = dynamic_cast<const GBVISoftcoreForce&>(force);
               writeGBVISoftcoreForce(  filePtr, obcForce );
               hit++;
            } catch( std::bad_cast ){
            }
        }
#endif
    
        // COM

        if( !hit ){
    
            try {
               const CMMotionRemover& cMMotionRemover = dynamic_cast<const CMMotionRemover&>(force);
               (void) fprintf( filePtr, "CMMotionRemover %d\n", cMMotionRemover.getFrequency() );
               hit++;
            } catch( std::bad_cast ){
            }
        }

        if( !hit ){
           char buffer[1024];
           (void) sprintf( buffer, "   %2d force not recognized.\n", i );
           (void) fprintf( stderr, "%s\n", buffer );
//           throw OpenMMException( buffer );
        }

    }

    // constraints

    writeConstraints( filePtr, system );

    // context

    writeContext( filePtr, context );

    // integrator 

    writeIntegrator( filePtr, integrator );

    // close file

    (void) fclose( filePtr );
}
