#ifndef VALIDATE_OPENMM_FORCES_H_
#define VALIDATE_OPENMM_FORCES_H_

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

#include "ValidateOpenMM.h"

namespace OpenMM {

typedef std::map< int, int > MapIntInt;
typedef MapIntInt::iterator MapIntIntI;
typedef MapIntInt::const_iterator MapIntIntCI;

// Helper class for ValidateOpenMMForces class used to store force results and facitilitate comparisons of
// resulting forces

class ForceValidationResult {
public:

    ForceValidationResult( const Context& context1, const Context& context2, StringUIntMap& forceNamesMap );
    ~ForceValidationResult();

    /**
     * Get potential energy at specified platform index (0 || 1)
     * 
     * @return  potential energy for spercifed platform
     *
     * @throws OpenMMException if energyIndex is not 0 or 1
     */
    double getPotentialEnergy( int energyIndex  ) const; 

    /**
     * Get array of forces at specified platform index (0 || 1)
     * 
     * @return array of force norms
     *
     * @throws OpenMMException if forceIndex is not 0 or 1
     */
    std::vector<double> getForceNorms( int forceIndex  ) const; 

    /**
     * Get array of forces at platform index (0 || 1)
     * 
     * @return force array
     *
     * @throws OpenMMException if forceIndex is not 0 or 1
     */
    std::vector<Vec3> getForces( int forceIndex  ) const; 

    /**
     * Get maximum delta in force norm
     * 
     * @param maxIndex  return atom index of entry with maximum delta norm (optional)
     *
     * @return max delta in norm of forces
     */
    double getMaxDeltaForceNorm( int* maxIndex = NULL  ) const;

    /**
     * Get maximum relative delta in force norm
     * 
     * @param maxIndex  return atom index of entry w/ maximum relative delta norm (optional)
     *
     * @return max relative delta in norm of forces
     */
    double getMaxRelativeDeltaForceNorm( int* maxIndex = NULL  ) const;

    /**
     * Get maximum dot product between forces
     * 
     * @param maxIndex  return atom index of entry w/ maximum dot product between forces (optional)
     *
     * @return max dot product between forces
     */
    double getMaxDotProduct( int* maxIndex = NULL  ) const;

    /**
     * Get name of force associated w/ computed results
     * 
     * @return force name(s); if more than one force active in computation, 
     * then names are concatenated and separated by '::' (e.g., 'NB_FORCE::GBSA_OBC_FORCE')
     */
    std::string getForceName( void ) const;

    /**
     * Get platform name
     * 
     * @param index index of platform (0 or 1)
     *
     * @return platform name
     *
     * @throws OpenMMException if index is not 0 or 1
     */
    std::string getPlatformName( int index ) const;

    /**
     * Register index of two entries that differ by a specified tolerance
     * 
     * @param index inconsistent index
     *
     */
    void registerInconsistentForceIndex( int index, int value = 1 );

    /**
     * Clear list of entries that differ by a specified tolerance
     * 
     */
    void clearInconsistentForceIndexList( void );

    /**
     * Get list of entries that differ by a specified tolerance
     * 
     */
    void getInconsistentForceIndexList( std::vector<int>& inconsistentIndices ) const;

    /**
     * Get number of entries in inconsistent index list
     * 
     */
    int getNumberOfInconsistentForceEntries( void ) const;

    /**
     * Return true if nans were detected
     * 
     * @return true if nans were detected
     */
    int nansDetected( void ) const;

    /**
     * Determine if force norms are valid
     * 
     * @param tolerance              tolerance
     */
    void compareForceNorms( double tolerance );

    /**
     * Determine if forces are valid
     * 
     * @param tolerance   tolerance
     */
    void compareForces( double tolerance );

private:

    // computed potential energies and forces fror two platforms

    double _potentialEnergies[2];
    std::vector<Vec3> _forces[2];

    // platform and force names

    std::string _platforms[2];
    std::vector<std::string> _forceNames;

    // force norms and stat entries

    std::vector<double> _norms[2];
    std::vector<double> _normStatVectors[2];

    // map of indicies w/ inconsistent force entries

    std::map<int, int> _inconsistentForceIndicies;

    // if set, then nans detected

    int _nansDetected;

    /**
     * Calculate norms of vectors
     * 
     */
    void _calculateNorms( void );

    /**
     * Calculate norms of specified vector
     * 
     */
    void _calculateNormOfForceVector( int forceIndex );

    // stat indices

    static const int STAT_AVG = 0;
    static const int STAT_STD = 1;
    static const int STAT_MIN = 2;
    static const int STAT_ID1 = 3;
    static const int STAT_MAX = 4;
    static const int STAT_ID2 = 5;
    static const int STAT_CNT = 6;
    static const int STAT_END = 7;

    /**
     * Find vector stats
     * 
     */
    void _findStatsForDouble( const std::vector<double>& array, std::vector<double>& statVector ) const;
};

// Class used to compare forces/potential energies on two platforms

class ValidateOpenMMForces : public ValidateOpenMM {
public:

    OPENMM_VALIDATE_EXPORT ValidateOpenMMForces( void );
    OPENMM_VALIDATE_EXPORT ~ValidateOpenMMForces();

    /**
     * Validate force/energy by comparing the results between the forces/energies computed on user-provided context platform
     * with Reference platform
     * 
     * @param context          context reference
     * @param summaryString    output summary string of results of comparison (optional)
     *
     * @return number of inconsistent entries
     */
     int OPENMM_VALIDATE_EXPORT compareWithReferencePlatform(Context& context, std::string* summaryString = NULL );

    /**
     * Validate force/energy by comparing the results between the forces/energies computed on two different platforms
     * 
     * @param context          context reference
     * @param compareForces    indices of force to be tested
     * @param platform1        first platform to compute forces
     * @param platform2        second platform to compute forces
     *
     * @return ForceValidationResult reference containing results of force/energy computations
     *         on the two input platforms
     */
     ForceValidationResult* compareForce(Context& context, std::vector<int>& compareForces,
                                          Platform& platform1, Platform& platform2 ) const;

    /**
     * Compare individual forces by comparing calculations across two platforms (platform associated w/ input context and
     * comparisonPlatform)
     * 
     * @param context                 context reference
     * @param platform                comparsion platform reference 
     * @param forceValidationResults  output vector of ForceValidationResult ptrs (user is responsible for deleting
     *                                individual ForceValidationResult objects)
     */
    void compareOpenMMForces(Context& context, Platform& comparisonPlatform, std::vector<ForceValidationResult*>& forceValidationResults ) const;

    /**
     * Determine if results are consistent
     * 
     * @param forceValidationResults  vector of ForceValidationResult ptrs to check if forces are consistent
     */
    void checkForInconsistentForceEntries( std::vector<ForceValidationResult*>& forceValidationResults ) const;

    /**
     * Get total number of force entries that are inconsistent
     * 
     * @param forceValidationResults  vector of ForceValidationResult ptrs to check if forces are consistent
     */
    int getTotalNumberOfInconsistentForceEntries( std::vector<ForceValidationResult*>& forceValidationResults ) const;

    /**
     * Get summary string of results
     * 
     * @param forceValidationResults  vector of ForceValidationResult ptrs
     */
    std::string getSummary( std::vector<ForceValidationResult*>& forceValidationResults ) const;

    /**
     * Set force tolerance
     * 
     * @param tolerance     force tolerance
     */
    void setForceTolerance( double tolerance );

    /**
     * Get force tolerance
     * 
     * @return force tolerance
     */
    double getForceTolerance( void ) const;

    /* 
     * Get force tolerance for specified force
     *
     * @param forceName   name of force
     *
     * @return force tolerance
     *
     * */
      
    double getForceTolerance( const std::string& forceName ) const;
       
    /* 
     * Get max errors to print in summary string
     *
     * @return max errors to print
     *
     * */
      
    int getMaxErrorsToPrint( void ) const;
       
    /* 
     * Set max errors to print in summary string
     *
     * @param maxErrorsToPrint max errors to print
     *
     * */
      
    void setMaxErrorsToPrint( int maxErrorsToPrint );
       
    /* 
     * Return true if force is not to be validated (Andersen thermostat, CM motion remover, ...)
     *
     * @param forceName   force name
     *
     * @return true if force is not currently validated
     **/
      
    int isExcludedForce( std::string forceName ) const;
       
private:

     // initialize class entries

     void _initialize( void );

      /* 
       * Format output line
       *
       * @param tab         tab
       * @param description description
       * @param value       value
       *
       * @return string containing contents
       *
       * */
      
      std::string _getLine( const std::string& tab,
                            const std::string& description,
                            const std::string& value ) const;
      
     std::vector<ForceValidationResult*> _forceValidationResults;

     // max errors to print

     int _maxErrorsToPrint;

     // tolerence

     double _forceTolerance;

     // map of force tolerances to type (name)

     StringDoubleMap _forceTolerances;
     
     // forces to be excluded from validation

     StringIntMap _forcesToBeExcluded;
};

} // namespace OpenMM

#endif /*VALIDATE_OPENMM_FORCES_H_*/
