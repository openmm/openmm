#ifndef VALIDATE_OPENMM_H_
#define VALIDATE_OPENMM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "OpenMM.h"
#include "../../../platforms/reference/include/ReferencePlatform.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "ValidateWindowsIncludes.h"

// free-energy plugin includes

//#define	INCLUDE_FREE_ENERGY_PLUGIN
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
#include "../../../plugins/freeEnergy/openmmapi/include/OpenMMFreeEnergy.h"
#include "../../../plugins/freeEnergy/openmmapi/include/openmm/freeEnergyKernels.h"
//#include "../../../plugins/freeEnergy/platforms/reference/include/ReferenceFreeEnergyPlatform.h"
#include "../../../plugins/freeEnergy/platforms/reference/include/ReferenceFreeEnergyKernelFactory.h"
#endif

#include <sstream>
#include <typeinfo>

#include <limits>
#include <cstdlib>
#include <cstring>
#include <cstdio>

namespace OpenMM {

typedef std::map< std::string, int > StringIntMap;
typedef StringIntMap::iterator StringIntMapI;
typedef StringIntMap::const_iterator StringIntMapCI;

typedef std::map< std::string, double > StringDoubleMap;
typedef StringDoubleMap::iterator StringDoubleMapI;
typedef StringDoubleMap::const_iterator StringDoubleMapCI;

typedef std::map< std::string, unsigned int > StringUIntMap;
typedef StringUIntMap::iterator StringUIntMapI;
typedef StringUIntMap::const_iterator StringUIntMapCI;

typedef std::vector< std::string > StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

typedef std::vector< int > IntVector;
typedef IntVector::iterator IntVectorI;
typedef IntVector::const_iterator IntVectorCI;

typedef std::map< std::string, std::vector<std::string> > StringStringVectorMap;
typedef StringStringVectorMap::iterator StringStringVectorMapI;
typedef StringStringVectorMap::const_iterator StringStringVectorMapCI;

/**
 * Base class w/ common functionality 
 */
class ValidateOpenMM {
public:

    ValidateOpenMM( void );
    ~ValidateOpenMM();

    // force names
    
    static const std::string HARMONIC_BOND_FORCE;
    static const std::string HARMONIC_ANGLE_FORCE;
    static const std::string PERIODIC_TORSION_FORCE;
    static const std::string RB_TORSION_FORCE;
    
    static const std::string NB_FORCE;
    static const std::string NB_SOFTCORE_FORCE;
    
    static const std::string NB_EXCEPTION_FORCE;
    static const std::string NB_EXCEPTION_SOFTCORE_FORCE;
    
    static const std::string GBSA_OBC_FORCE;
    static const std::string GBSA_OBC_SOFTCORE_FORCE;
    
    static const std::string GBVI_FORCE;
    static const std::string GBVI_SOFTCORE_FORCE;

    static const std::string CM_MOTION_REMOVER;
    static const std::string ANDERSEN_THERMOSTAT;
    static const std::string CUSTOM_BOND_FORCE;
    static const std::string CUSTOM_EXTERNAL_FORCE;
    static const std::string CUSTOM_NONBONDED_FORCE;

    /**
     * Return true if input number is nan or infinity
     *
     * @param number   number to test
     *
     * @return true if number is nan or infinity
     */
    
    static int isNanOrInfinity( double number );
    
    /**
     * Get force name
     * 
     * @param force      OpenMM system force 
     *
     * @return force name or "NA" if force is not recognized
     */
    std::string getForceName(const Force& force ) const;

    /**
     * Copy force 
     * 
     * @param force      OpenMM system force to copy 
     *
     * @return force or NULL if not recognized
     */
    Force* copyForce(const Force& force) const;

    /**
     * Get copy of input system, but omit forces
     * 
     * @param systemToCopy   system to copy
     *
     * @return copy of system but w/o forces
     */
    System* copySystemExcludingForces( const System& systemToCopy ) const;

    /**
     *
     * Set the velocities/positions of context2 to those of context1
     * 
     * @param context1                 context1 
     * @param context2                 context2 
     *
     * @return 0
     */
    void synchContexts( const Context& context1, Context& context2 ) const;
    
    /**
     *
     * Get log FILE* reference
     * 
     * @return log
     *
     */
    FILE* getLog( ) const;
    
    /**
     *
     * Set log FILE* reference
     * 
     * @param log  log
     *
     */
    void OPENMM_VALIDATE_EXPORT setLog( FILE* log );
    
    /**---------------------------------------------------------------------------------------
    
       Copy constraints
    
       @param systemToCopy         system whose constraints are to be copied
       @param system               system to add constraints to
       @param log                  log file pointer -- may be NULL
    
       --------------------------------------------------------------------------------------- */
    
    void copyConstraints( const System& systemToCopy, System* system, FILE* log = NULL ) const;
    
    /**---------------------------------------------------------------------------------------
    
       Get force dependencies
    
       @param forceName            force to check if there exist any dependencies
       @param returnVector         vector of forces the input force is dependent on (example: GBSAOBC force requires Nonbonded force since
                                   on Cuda platofrm they are computed in same loop and hence are inseparable)
    
       --------------------------------------------------------------------------------------- */
    
    void getForceDependencies( std::string forceName, StringVector& returnVector ) const;
    
    /**
     * Write masses to parameter file
     *
     * @param filePtr    file to write masses to
     * @param system     write masses in system
     */
    
    void writeMasses( FILE* filePtr, const System& system ) const;
    
    /**
     * Write constraints to parameter file
     *
     * @param filePtr    file to write constraints to
     * @param system     write constraints in system
     *
     */
    
    void writeConstraints( FILE* filePtr, const System& system ) const;
    
    /**
     * Write harmonicBondForce parameters to file
     *
     * @param filePtr                file to write forces to
     * @param harmonicBondForce     write harmonicBondForce parameters
     *
     */
    
    void writeHarmonicBondForce( FILE* filePtr, const HarmonicBondForce& harmonicBondForce ) const;

    /**
     * Write harmonicAngleForce parameters to file
     *
     * @param filePtr                file to write forces to
     * @param harmonicAngleForce     write harmonicAngleForce parameters
     *
     */
    
    void writeHarmonicAngleForce( FILE* filePtr, const HarmonicAngleForce& harmonicAngleForce ) const;
    
    /**
     * Write rbTorsionForce parameters to file
     *
     * @param filePtr                file to write forces to
     * @param rbTorsionForce         write rbTorsionForce parameters
     *
     */
    
    void writeRbTorsionForce( FILE* filePtr, const RBTorsionForce& rbTorsionForce ) const;
    
    /**
     * Write periodicTorsionForce parameters to file
     *
     * @param filePtr                file to write forces to
     * @param periodicTorsionForce   write periodicTorsionForce parameters
     *
     */
    
    void writePeriodicTorsionForce( FILE* filePtr, const PeriodicTorsionForce& periodicTorsionForce ) const;
    
    /**
     * Write nonbonded parameters to file
     *
     * @param filePtr                file to write forces to
     * @param nonbondedForce         write nonbondedForce parameters
     *
     */
    
    void writeNonbondedForce( FILE* filePtr, const NonbondedForce & nonbondedForce ) const;
    
    /**
     * Write GBSAOBCForce parameters to file
     *
     * @param filePtr                file to write forces to
     * @param gbsaObcForce           write gbsaObcForce parameters
     *
     */
    
    void writeGbsaObcForce( FILE* filePtr, const GBSAOBCForce& gbsaObcForce ) const;
    
    /**
     * Write GBSA GB/VI Force parameters to file
     *
     * @param filePtr                file to write forces to
     * @param gbviObcForce           write gbviObcForce parameters
     *
     */
    
    void writeGBVIForce( FILE* filePtr, const GBVIForce& gbviForce ) const;
    
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
    
    /**
     * Write nonbonded softcore parameters to file
     *
     * @param filePtr                file to write forces to
     * @param nonbondedForce         write nonbondedForce parameters
     */
    
    void writeNonbondedSoftcoreForce( FILE* filePtr, const NonbondedSoftcoreForce & nonbondedSoftcoreForce ) const;

    /**
     * Write GBSAOBCSoftcoreForce parameters to file
     *
     * @param filePtr                file to write forces to
     * @param gbsaObcForce           write gbsaObcForce parameters
     *
     */
    
    void writeGbsaObcSoftcoreForce( FILE* filePtr, const GBSAOBCSoftcoreForce& gbsaObcForce ) const;
    
    /**
     * Write GBSA GB/VI softcore force parameters to file
     *
     * @param filePtr                file to write forces to
     * @param gbviObcForce           write gbviObcForce parameters
     *
     */
    
    void writeGBVISoftcoreForce( FILE* filePtr, const GBVISoftcoreForce& gbviSoftcoreForce ) const;
    
#endif
    
    /**
     * Write coordinates, velocities, ... to file
     *
     * @param filePtr                file to write Vec3 entries to
     * @param vect3Array             write array of Vec3
     *
     */
    
    void writeVec3( FILE* filePtr, const std::vector<Vec3>& vect3Array ) const;
    
    /**
     * Write context info to file (positions, velocities, forces, energies)
     *
     * @param filePtr                file to write entries to
     * @param context                write context positions, velocities, forces, energies to file
     *
     */
    
    void writeContext( FILE* filePtr, const Context& context ) const;
    
    /**
     * Write integrator
     *
     * @param filePtr                file to write integrator info to
     * @param integrator             write integrator info (time step, seed, ... as applicable)
     *
     */
    
    void writeIntegrator( FILE* filePtr, const Integrator& integrator ) const;
    
    /**
     * Write parameter file
     * @param context                context whose entries are to be written to file
     * @param parameterFileName      file name
     *
     */
    
    void writeParameterFile( const Context& context, const std::string& parameterFileName ) const;
    
private:

    FILE* _log;

    // map of force dependencies (e.g., GBSAObc requires NB force on CudaPlatform)

    StringStringVectorMap _forceDependencies;
};


} // namespace OpenMM

#endif /*VALIDATE_OPENMM_H_*/
