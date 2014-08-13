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

#include <iostream>
#include <sstream>
#include "openmm/OpenMMException.h"
#include "ValidateOpenMMForces.h"

using namespace OpenMM;

ForceValidationResult::ForceValidationResult( const Context& context1, const Context& context2, StringUIntMap& forceNamesMap  ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::ForceValidationResult";

// ---------------------------------------------------------------------------------------

    _nansDetected                             = 0;

    // calculate forces/energies

    int types                                 = State::Forces | State::Energy;
    State state1                              = context1.getState( types );
    State state2                              = context2.getState( types );

    _potentialEnergies[0]                     = state1.getPotentialEnergy();
    _potentialEnergies[1]                     = state2.getPotentialEnergy();
    if( ValidateOpenMM::isNanOrInfinity( _potentialEnergies[0] ) || ValidateOpenMM::isNanOrInfinity( _potentialEnergies[1] ) ){
       _nansDetected++;
    }

    _forces[0]                                = state1.getForces();
    _forces[1]                                = state2.getForces();

    _platforms[0]                             = context1.getPlatform().getName();
    _platforms[1]                             = context2.getPlatform().getName();

    for( StringUIntMapI ii = forceNamesMap.begin(); ii != forceNamesMap.end(); ii++ ){
       _forceNames.push_back( (*ii).first );
    }
    _calculateNorms();
}

ForceValidationResult::~ForceValidationResult(  ){
}

void ForceValidationResult::_calculateNorms( void ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::_calculateNorms";

// ---------------------------------------------------------------------------------------

    for( int jj = 0; jj < 2; jj++ ){
        if( _norms[jj].size() < 1 ){
            _calculateNormOfForceVector( jj );
            _findStatsForDouble( _norms[jj], _normStatVectors[jj] );
        }
    }
}

/**---------------------------------------------------------------------------------------

   Find stats for double array

   @param array                   array 
   @param statVector              vector of stats 

   @return 0

   --------------------------------------------------------------------------------------- */

void ForceValidationResult::_findStatsForDouble( const std::vector<double>& array, std::vector<double>& statVector ) const {

// ---------------------------------------------------------------------------------------

    //static const char* methodName  = "ForceValidationResult::_findStatsForDouble";

// ---------------------------------------------------------------------------------------

    statVector.resize( STAT_END );
 
    double avgValue   =  0.0;
    double stdValue   =  0.0;
    double minValue   =  1.0e+30;
    double maxValue   = -1.0e+30;
    int minValueIndex = 0;
    int maxValueIndex = 0;
 
    for( unsigned int ii = 0; ii < array.size(); ii++ ){
 
       double norm  =  array[ii];
 
       avgValue    += norm;
       stdValue    += norm*norm;
 
       if( norm > maxValue ){
          maxValue       = norm;
          maxValueIndex  = ii;
       }
       if( norm < minValue ){
          minValue       = norm;
          minValueIndex  = ii;
       }
    }
 
    double count  = static_cast<double>(array.size());
    double iCount = count > 0.0 ? 1.0/count : 0.0;
 
    statVector[STAT_AVG] = avgValue*iCount;
    statVector[STAT_STD] = stdValue - avgValue*avgValue*count;
    if( count > 1.0 ){
       statVector[STAT_STD] = std::sqrt( stdValue/( count - 1.0 ) );
    }
    statVector[STAT_MIN] = minValue;
    statVector[STAT_ID1] = static_cast<double>(minValueIndex);
    statVector[STAT_MAX] = maxValue;
    statVector[STAT_ID2] = static_cast<double>(maxValueIndex);
    statVector[STAT_CNT] = count;
 
    return;
}

void ForceValidationResult::_calculateNormOfForceVector( int forceIndex ){

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::_calculateNormOfForceVector";

// ---------------------------------------------------------------------------------------

    for( unsigned int ii = 0; ii < _forces[forceIndex].size(); ii++ ){
        Vec3 f1 = _forces[forceIndex][ii];
        _norms[forceIndex].push_back( std::sqrt( (f1[0]*f1[0]) + (f1[1]*f1[1]) + (f1[2]*f1[2]) ) );
        if( ValidateOpenMM::isNanOrInfinity( f1[0] ) || ValidateOpenMM::isNanOrInfinity( f1[1] ) || ValidateOpenMM::isNanOrInfinity( f1[2] ) ){
            _nansDetected++;
        }
    }

    return;
}

double ForceValidationResult::getPotentialEnergy( int energyIndex ) const {
    if( energyIndex == 0 || energyIndex == 1 ){
        return _potentialEnergies[energyIndex];
    } else {
       char buffer[1024];
       (void) sprintf( buffer, "getPotentialEnergy: energyIndex=%d is inconsistent.", energyIndex );
       throw OpenMMException( buffer );
   }
}

std::vector<Vec3> ForceValidationResult::getForces( int forceIndex  ) const {
    if( forceIndex == 0 || forceIndex == 1 ){
        return _forces[forceIndex];
    } else {
       char buffer[1024];
       (void) sprintf( buffer, "getForces: forceIndex=%d is inconsistent.", forceIndex );
       throw OpenMMException( buffer );
   }
}

std::vector<double> ForceValidationResult::getForceNorms( int forceIndex  ) const {
    if( forceIndex == 0 || forceIndex == 1 ){
        return _norms[forceIndex];
    } else {
       char buffer[1024];
       (void) sprintf( buffer, "getForceNorms: forceIndex=%d is inconsistent.", forceIndex );
       throw OpenMMException( buffer );
   }
}

int ForceValidationResult::nansDetected( void ) const {
    return _nansDetected;
}

void ForceValidationResult::registerInconsistentForceIndex( int index, int value ){
    _inconsistentForceIndicies[index] = value;
}

void ForceValidationResult::clearInconsistentForceIndexList( void ){
    _inconsistentForceIndicies.clear();
}

void ForceValidationResult::getInconsistentForceIndexList( std::vector<int>& inconsistentIndices ) const {
    for( MapIntIntCI ii = _inconsistentForceIndicies.begin(); ii != _inconsistentForceIndicies.end(); ii++ ){
        inconsistentIndices.push_back( (*ii).first );
    }
}

int ForceValidationResult::getNumberOfInconsistentForceEntries( void ) const {
    return  static_cast<int>(_inconsistentForceIndicies.size() );
}

double ForceValidationResult::getMaxDeltaForceNorm( int* maxIndex  ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::calculateNorm";

// ---------------------------------------------------------------------------------------

    double maxDelta                           = -1.0e+90;
    int maxDeltaIndex                         = -1;

    for( unsigned int ii = 0; ii < _forces[0].size(); ii++ ){
 
        Vec3 f1                = _forces[0][ii];
        double normF1          = _norms[0][ii];
     
        Vec3 f2                = _forces[1][ii];
        double normF2          = _norms[1][ii];
  
        double delta           = std::sqrt( (f1[0]-f2[0])*(f1[0]-f2[0]) + (f1[1]-f2[1])*(f1[1]-f2[1]) + (f1[2]-f2[2])*(f1[2]-f2[2]) );
        if( maxDelta < delta ){
           maxDelta      = delta;
           maxDeltaIndex = static_cast<int>(ii);
        }
    }

    if( maxIndex ){
        *maxIndex = maxDeltaIndex;
    }
         
    return maxDelta;
}

double ForceValidationResult::getMaxRelativeDeltaForceNorm( int* maxIndex  ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::getMaxRelativeDeltaForceNorm";

// ---------------------------------------------------------------------------------------

    double maxRelativeDelta    = -1.0e+90;
    int maxRelativeDeltaIndex  = -1;

    for( unsigned int ii = 0; ii < _forces[0].size(); ii++ ){
 
        Vec3 f1                = _forces[0][ii];
        double normF1          = _norms[0][ii];

        Vec3 f2                = _forces[1][ii];
        double normF2          = _norms[1][ii];

        double delta           = std::sqrt( (f1[0]-f2[0])*(f1[0]-f2[0]) + (f1[1]-f2[1])*(f1[1]-f2[1]) + (f1[2]-f2[2])*(f1[2]-f2[2]) );
        double forceSum        = 0.5*(normF1 + normF2);
        if( forceSum > 0.0 ){
            delta /= forceSum;
        }
  
        if( maxRelativeDelta < delta ){
           maxRelativeDelta      = delta;
           maxRelativeDeltaIndex = static_cast<int>(ii);
        }
    }

    if( maxIndex ){
        *maxIndex = maxRelativeDeltaIndex;
    }
         
    return maxRelativeDelta;
}

double ForceValidationResult::getMaxDotProduct( int* maxIndex  ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::getMaxDotProduct";

// ---------------------------------------------------------------------------------------

    double maxDotProduct       = -1.0e+90;
    int maxDotProductIndex     = -1;

    for( unsigned int ii = 0; ii < _forces[0].size(); ii++ ){
 
        Vec3 f1                = _forces[0][ii];
        double normF1          = _norms[0][ii];

        Vec3 f2                = _forces[1][ii];
        double normF2          = _norms[1][ii];

        double dotProduct      = f1[0]*f2[0] + f1[1]*f2[1] + f1[2]*f2[2];

        if(  (normF1*normF2) > 0.0 ){
            dotProduct     /= (normF1*normF2);
            dotProduct      = 1.0 - dotProduct;
        } else {
            dotProduct      = 0.0;
        }
 
        if( maxDotProduct < dotProduct ){
           maxDotProduct      = dotProduct;
           maxDotProductIndex = static_cast<int>(ii);
        }
    }

    if( maxIndex ){
        *maxIndex = maxDotProductIndex;
    }
         
    return maxDotProduct;
}

std::string ForceValidationResult::getForceName( void ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::getForceName";

// ---------------------------------------------------------------------------------------

    std::stringstream forceName;
    for( unsigned int ii = 0; ii < _forceNames.size(); ii++ ){
       forceName << _forceNames[ii];
       if( ii < (_forceNames.size()-1) )forceName << "::";
    }
    return forceName.str();
}

std::string ForceValidationResult::getPlatformName( int index ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ForceValidationResult::getPlatformName";

// ---------------------------------------------------------------------------------------

    if( index == 0 || index == 1 ){
        return _platforms[index];
    } else {
       char buffer[1024];
       (void) sprintf( buffer, "getPlatformName: index=%d is inconsistent.", index );
       throw OpenMMException( buffer );
   }
}

void ForceValidationResult::compareForces( double tolerance ){

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = " ForceValidationResult::compareForces";

// ---------------------------------------------------------------------------------------

    std::vector<Vec3> forces1          = getForces( 0 );
    std::vector<double> forceNorms1    = getForceNorms( 0 );

    std::vector<Vec3> forces2          = getForces( 1 );
    std::vector<double> forceNorms2    = getForceNorms( 1 );

    clearInconsistentForceIndexList( );
    for( unsigned int jj = 0; jj < forces1.size(); jj++ ){
        if( ValidateOpenMM::isNanOrInfinity( forceNorms1[jj] ) ||  ValidateOpenMM::isNanOrInfinity( forceNorms2[jj] ) ){
            registerInconsistentForceIndex( jj );
        } else {
            double delta   = std::sqrt( (forces1[jj][0]-forces2[jj][0])*(forces1[jj][0]-forces2[jj][0]) +
                                        (forces1[jj][1]-forces2[jj][1])*(forces1[jj][1]-forces2[jj][1]) +
                                        (forces1[jj][2]-forces2[jj][2])*(forces1[jj][2]-forces2[jj][2]) );

            double sum     = 0.5*( fabs( forceNorms1[jj] ) + fabs( forceNorms2[jj] ) );
            double diff    = delta > 0.0 ? delta/sum : 0.0;
            if( diff > tolerance && sum > tolerance ){
                registerInconsistentForceIndex( jj );
            }
        }
    }

}

void ForceValidationResult::compareForceNorms( double tolerance ){

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "ForceValidationResult::compareForceNorms";

// ---------------------------------------------------------------------------------------

    std::vector<double> forceNorms0    = getForceNorms( 0 );
    std::vector<double> forceNorms1    = getForceNorms( 1 );

    clearInconsistentForceIndexList( );
    for( unsigned int jj = 0; jj < forceNorms0.size(); jj++ ){
        if( ValidateOpenMM::isNanOrInfinity( forceNorms0[jj] ) || ValidateOpenMM::isNanOrInfinity( forceNorms1[jj] ) ){
            registerInconsistentForceIndex( jj );
        } else {
            double sum  = 0.5*( fabs( forceNorms0[jj] ) + fabs( forceNorms1[jj] ) );
            double diff = sum > 0.0 ? fabs( forceNorms0[jj] - forceNorms1[jj] )/sum : 0.0;
            if( diff > tolerance && sum > tolerance ){
                registerInconsistentForceIndex( jj );
            }
        }
    }
}

ValidateOpenMMForces::ValidateOpenMMForces( void ) {
    _initialize();
}

void ValidateOpenMMForces::_initialize( void ){

    _forceTolerance                               = 1.0e-02;
    _maxErrorsToPrint                             = 25;

    // force tolerances by name

    _forceTolerances[HARMONIC_BOND_FORCE]         = _forceTolerance;
    _forceTolerances[HARMONIC_ANGLE_FORCE]        = _forceTolerance;
    //_forceTolerances[PERIODIC_TORSION_FORCE]      = _forceTolerance;
    _forceTolerances[PERIODIC_TORSION_FORCE]      = 0.3;
    //_forceTolerances[RB_TORSION_FORCE]            = _forceTolerance;
    _forceTolerances[RB_TORSION_FORCE]            = 0.3;
    _forceTolerances[NB_FORCE]                    = _forceTolerance;
    _forceTolerances[NB_SOFTCORE_FORCE]           = _forceTolerance;
    _forceTolerances[GBSA_OBC_FORCE]              = _forceTolerance;
    _forceTolerances[GBSA_OBC_SOFTCORE_FORCE]     = _forceTolerance;
    _forceTolerances[GBVI_FORCE]                  = _forceTolerance;
    _forceTolerances[GBVI_SOFTCORE_FORCE]         = _forceTolerance;

    // forces to be excluded from validation

    _forcesToBeExcluded[CM_MOTION_REMOVER]        = 1;
    _forcesToBeExcluded[ANDERSEN_THERMOSTAT]      = 1;
}

ValidateOpenMMForces::~ValidateOpenMMForces( ){

    for( unsigned int ii = 0; ii < _forceValidationResults.size(); ii++ ){
         delete _forceValidationResults[ii];
    }
    _forceValidationResults.resize( 0 );
    
}

double ValidateOpenMMForces::getForceTolerance( void ) const {
    return _forceTolerance;
}

void ValidateOpenMMForces::setForceTolerance( double forceTolerance ){
    _forceTolerance = forceTolerance;
}

int ValidateOpenMMForces::compareWithReferencePlatform( Context& context, std::string* summaryString ){

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "ValidateOpenMMForces::compareWithReferencePlatform";

// ---------------------------------------------------------------------------------------

    ReferencePlatform referencePlatform;
#ifdef INCLUDE_FREE_ENERGY_PLUGIN
    ReferenceFreeEnergyKernelFactory* factory  = new ReferenceFreeEnergyKernelFactory();
    referencePlatform.registerKernelFactory(CalcNonbondedSoftcoreForceKernel::Name(), factory);
    referencePlatform.registerKernelFactory(CalcGBSAOBCSoftcoreForceKernel::Name(), factory);
    referencePlatform.registerKernelFactory(CalcGBVISoftcoreForceKernel::Name(), factory);
#endif

    compareOpenMMForces( context, referencePlatform, _forceValidationResults );

    checkForInconsistentForceEntries( _forceValidationResults );

    if( summaryString ){
        *summaryString = getSummary( _forceValidationResults );
    }

    return getTotalNumberOfInconsistentForceEntries( _forceValidationResults );
}

void ValidateOpenMMForces::compareOpenMMForces( Context& context, Platform& platform,
                                                std::vector<ForceValidationResult*>& forceValidationResults ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "ValidateOpenMMForces::compareOpenMMForces";

// ---------------------------------------------------------------------------------------

    // loop over forces
 
    const System& system         = context.getSystem();
    std::vector<int> compareAllForces;
    for( int ii = 0; ii < system.getNumForces(); ii++ ){
        std::vector<int> compareForces;
        std::string forceName  = getForceName( system.getForce( ii ) );
        if( isExcludedForce( forceName ) == 0 ){
            compareForces.push_back( ii );
            compareAllForces.push_back( ii );
            ForceValidationResult* forceValidationResult = compareForce(context, compareForces, platform, context.getPlatform() );
            forceValidationResults.push_back( forceValidationResult );
        } else if( getLog() ){
           (void) fprintf( getLog(), "Force %s is not being validated.\n", forceName.c_str() );
        }
    }
    ForceValidationResult* forceValidationResult = compareForce(context, compareAllForces, platform, context.getPlatform() );
    forceValidationResults.push_back( forceValidationResult );
}

ForceValidationResult* ValidateOpenMMForces::compareForce(Context& context, std::vector<int>& compareForces,
                                                          Platform& platform1, Platform& platform2 ) const {

// ---------------------------------------------------------------------------------------

    //static const std::string methodName      = "ValidateOpenMMForces::compareForce";

// ---------------------------------------------------------------------------------------

    // note if platforms are identical

    if( getLog() && platform1.getName().compare( platform2.getName() ) == 0 ){
        (void) fprintf( getLog(), "Note: Platforms to compares %s are identical.\n", platform1.getName().c_str() );
        (void) fflush( getLog() );
    }

    const System& system         = context.getSystem();

    // collect systemForceNameMap[forceName] = index in system
    //         systemForceNameIndex[index]   = force name

    StringIntMap systemForceNameMap;
    StringVector systemForceNameIndex;
    systemForceNameIndex.resize( system.getNumForces() );
    for( int ii = 0; ii < system.getNumForces(); ii++ ){
        std::string forceName         = getForceName( system.getForce( ii ) );
        if( forceName.compare( "NA" ) == 0 ){
            std::stringstream message;
            message << "Force at index=" << ii << " not found -- aborting!";
            std::cerr << message.str() << std::endl;
            throw OpenMM::OpenMMException(message.str());
        }
        systemForceNameMap[forceName] = ii;
        systemForceNameIndex[ii]      = forceName;
    }

    // diagnostics

    if( 0 && getLog() ){
        for( StringIntMapI ii = systemForceNameMap.begin(); ii != systemForceNameMap.end(); ii++ ){
            int index = (*ii).second;
            (void) fprintf( getLog(), "  System force map %s index=%d reverse map=%s\n", (*ii).first.c_str(), index, systemForceNameIndex[index].c_str() );
        }
        for( unsigned int ii = 0; ii < compareForces.size(); ii++ ){
           (void) fprintf( getLog(), "   ValidateOpenMMForces %u %s\n", ii, systemForceNameIndex[compareForces[ii]].c_str() );
        }
        (void) fflush( getLog() );
    }

    // get system copy and add forces to system

    System* validationSystem     = copySystemExcludingForces( system );
    StringUIntMap forceNamesMap;
    for( unsigned int ii = 0; ii < compareForces.size(); ii++ ){
        const Force& forceToCopy = system.getForce( compareForces[ii] );
        Force* force             = copyForce( forceToCopy );
        validationSystem->addForce( force );
        forceNamesMap[systemForceNameIndex[compareForces[ii]]] = ii;
    }

    // include any missing dependencies (e.g, OBC force requires NB force for Cuda platform)

    for( StringUIntMapI ii = forceNamesMap.begin(); ii != forceNamesMap.end(); ii++ ){
       std::string forceName = (*ii).first;
       StringVector dependencyVector;
       getForceDependencies( forceName, dependencyVector ); 
       for( unsigned int jj = 0; jj < dependencyVector.size(); jj++ ){
           std::string dependentForceName = dependencyVector[jj];
           StringUIntMapCI dependent      = forceNamesMap.find( dependentForceName );
           if( dependent == forceNamesMap.end() ){
              forceNamesMap[dependentForceName] = 1;
              int forceIndex                    = systemForceNameMap[dependentForceName];
              const Force& forceToCopy          = system.getForce( forceIndex );
              validationSystem->addForce( copyForce( forceToCopy ) ); 
           }
       }
    }

    // create contexts

    VerletIntegrator verletIntegrator( 0.001 );
    Context* validationContext1  = new Context( *validationSystem, verletIntegrator, platform1);
    Context* validationContext2  = new Context( *validationSystem, verletIntegrator, platform2);

    // set positions

    synchContexts( context, *validationContext1 );
    synchContexts( context, *validationContext2 );

    // diagnostics

    if( 0 && getLog() ){
        std::stringstream forceNames;
        (void) fprintf( getLog(), "    Validating system forces=%d\n", validationSystem->getNumForces() );
        for( int ii = 0; ii < validationSystem->getNumForces(); ii++ ){
            std::string forceName         = getForceName( validationSystem->getForce( ii ) );
            forceNames << forceName;
            if( ii < (validationSystem->getNumForces()-1) ){
                forceNames << "_";
            } else {
                forceNames << "Parameters.txt";
            }
            (void) fprintf( getLog(), "       force %d %s\n", ii, forceName.c_str() );
        }
        writeParameterFile( *validationContext1, forceNames.str() ); 
        (void) fflush( getLog() );
    }

    // calculate forces & build return result

    ForceValidationResult* forceValidationResult = new ForceValidationResult( *validationContext1, *validationContext2, forceNamesMap );

    delete validationContext1;
    delete validationContext2;
    delete validationSystem;
        
    return forceValidationResult;
}


void ValidateOpenMMForces::checkForInconsistentForceEntries( std::vector<ForceValidationResult*>& forceValidationResults  ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "ValidateOpenMMForces::checkForInconsistentForceEntries";

// ---------------------------------------------------------------------------------------

    for( unsigned int ii = 0; ii < forceValidationResults.size(); ii++ ){
        std::string forceName              = forceValidationResults[ii]->getForceName();
        double tolerance                   = getForceTolerance( forceName );
        forceValidationResults[ii]->compareForces( tolerance ); 
    }

}

int ValidateOpenMMForces::getTotalNumberOfInconsistentForceEntries( std::vector<ForceValidationResult*>& forceValidationResults  ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "ValidateOpenMMForces::getTotalNumberOfInconsistentForceEntries";

// ---------------------------------------------------------------------------------------

    int inconsistentEntries = 0;
    for( unsigned int ii = 0; ii < forceValidationResults.size(); ii++ ){
        inconsistentEntries  += forceValidationResults[ii]->getNumberOfInconsistentForceEntries();
    }

    return inconsistentEntries;
}

int ValidateOpenMMForces::isExcludedForce( std::string forceName ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "ValidateOpenMMForces::isExcludedForce";

// ---------------------------------------------------------------------------------------

    StringIntMapCI isPresent = _forcesToBeExcluded.find( forceName );
    return isPresent != _forcesToBeExcluded.end() ? 1 : 0;
}


/* 
 * Get force tolerance for specified force
 *
 * @param forceName   name of force
 *
 * @return force tolerance
 *
 * */

double ValidateOpenMMForces::getForceTolerance( const std::string& forceName ) const {
    
    StringDoubleMapCI forceIsPresent = _forceTolerances.find( forceName );
    if( forceIsPresent != _forceTolerances.end() ){
        return (*forceIsPresent).second;
    }
    return _forceTolerance;
}

int ValidateOpenMMForces::getMaxErrorsToPrint( void ) const {
    return _maxErrorsToPrint;
}

void ValidateOpenMMForces::setMaxErrorsToPrint( int maxErrorsToPrint ){
    _maxErrorsToPrint = maxErrorsToPrint;
}

/* 
 * Format line
 *
 * @param tab         tab
 * @param description description
 * @param value       value
 *
 * @return string containing contents
 *
 * */

std::string ValidateOpenMMForces::_getLine( const std::string& tab, 
                                            const std::string& description,
                                            const std::string& value ) const {

// ---------------------------------------------------------------------------------------

    static const std::string methodName      = "GromacsOpenMMTest::_getLine";
    static const unsigned int MAX_LINE_CHARS = 256; 
    char line[MAX_LINE_CHARS];
 
 // ---------------------------------------------------------------------------------------
 
    std::stringstream message;
    memset( line, ' ', MAX_LINE_CHARS );   
#ifdef WIN32
    (void) sprintf_s( line, MAX_LINE_CHARS, "%s %-40s %s", tab.c_str(), description.c_str(), value.c_str() );
#else
    (void) sprintf( line, "%s %-40s %s", tab.c_str(), description.c_str(), value.c_str() );
#endif
    message << std::string( line ) << std::endl;
 
    return message.str();

}

std::string ValidateOpenMMForces::getSummary( std::vector<ForceValidationResult*>& forceValidationResults ) const {

// ---------------------------------------------------------------------------------------

    // static const std::string methodName      = "ValidateOpenMMForces::getSummary";
   static const unsigned int MAX_LINE_CHARS = 256; 
   char value[MAX_LINE_CHARS];

// ---------------------------------------------------------------------------------------

#ifdef WIN32
#define LOCAL_SPRINTF0(a,b) sprintf_s( (a), MAX_LINE_CHARS, (b) );   
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#define LOCAL_SPRINTF4(a,b,c,d) sprintf_s( (a), MAX_LINE_CHARS, (b), (c), (d) );   
#define LOCAL_SPRINTF5(a,b,c,d,e) sprintf_s( (a), MAX_LINE_CHARS, (b), (c), (d), (e) );   
#define LOCAL_SPRINTF6(a,b,c,d,e,f) sprintf_s( (a), MAX_LINE_CHARS, (b), (c), (d), (e), (f) );   
#define LOCAL_SPRINTF12(a,b,c,d,e,f,g,h,i,j,k,l) sprintf_s( (a), MAX_LINE_CHARS, (b), (c), (d), (e), (f), (g), (h), (i), (j), (k), (l) );   
#else
#define LOCAL_SPRINTF0(a,b) sprintf( (a), (b) );   
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#define LOCAL_SPRINTF4(a,b,c,d) sprintf( (a), (b), (c), (d) );   
#define LOCAL_SPRINTF5(a,b,c,d,e) sprintf( (a), (b), (c), (d), (e) );   
#define LOCAL_SPRINTF6(a,b,c,d,e,f) sprintf( (a), (b), (c), (d), (e), (f) );   
#define LOCAL_SPRINTF12(a,b,c,d,e,f,g,h,i,j,k,l) sprintf( (a), (b), (c), (d), (e), (f), (g), (h), (i), (j), (k), (l) );   
#endif

    unsigned int maxMissesToPrint = static_cast<unsigned int>(getMaxErrorsToPrint());
    std::stringstream summary;
    std::string tab      = "   ";
    int useForces        = 1;
    for( unsigned int ii = 0; ii < forceValidationResults.size(); ii++ ){

        // platforms

        if( ii == 0 ){
            std::string platform1 = forceValidationResults[ii]->getPlatformName( 0 );
            std::string platform2 = forceValidationResults[ii]->getPlatformName( 1 );
            (void) LOCAL_SPRINTF4( value, "%10s  %10s\n", platform1.c_str(), platform2.c_str());
            summary << _getLine( tab, "Platforms", value );
        }

        // force name
        std::string forceName              = forceValidationResults[ii]->getForceName();
        double tolerance                   = getForceTolerance( forceName );
        summary << _getLine( tab, "Force", forceName );
        (void) LOCAL_SPRINTF( value, "%.3e ", tolerance );
        summary << _getLine( tab, "Tolerance", value );
        int index;

        // delta

        double delta = forceValidationResults[ii]->getMaxDeltaForceNorm( &index );
        (void) LOCAL_SPRINTF4( value, "%.3e at index %6d", delta, index );
        summary << _getLine( tab, "Max Delta", value );

        // relative delta

        delta        = forceValidationResults[ii]->getMaxRelativeDeltaForceNorm( &index );
        (void) LOCAL_SPRINTF4( value, "%.3e at index %6d", delta, index );
        summary << _getLine( tab, "Max Relative Delta", value );

        // PE delta

        double pe1   = forceValidationResults[ii]->getPotentialEnergy( 0 );
        double pe2   = forceValidationResults[ii]->getPotentialEnergy( 1 );
        double sum   = 0.5*(fabs( pe1 ) + fabs( pe2 ) );
        delta        = sum > 0.0 ? fabs( pe1 - pe2 )/sum : 0.0;
        (void) LOCAL_SPRINTF5( value, "%10.4e PE[%14.6e %14.6e]", delta, pe1, pe2 );
        summary << _getLine( tab, "Potential energies relative delta", value );

        // misses

        std::vector<int> inconsistentIndices;
        forceValidationResults[ii]->getInconsistentForceIndexList( inconsistentIndices );
        if( inconsistentIndices.size() > 0 ){ 
    
            std::vector<Vec3> forces1          = forceValidationResults[ii]->getForces( 0 );
            std::vector<Vec3> forces2          = forceValidationResults[ii]->getForces( 1 );
    
            std::vector<double> forceNorms1    = forceValidationResults[ii]->getForceNorms( 0 );
            std::vector<double> forceNorms2    = forceValidationResults[ii]->getForceNorms( 1 );
            for( unsigned int kk = 0; kk < inconsistentIndices.size() && kk < maxMissesToPrint; kk++ ){
                int jj = inconsistentIndices[kk];
                if( isNanOrInfinity( forceNorms1[jj] ) || isNanOrInfinity( forceNorms2[jj] ) ){
                     (void) LOCAL_SPRINTF5( value, "         nan at index %6d  norms: [%12.5e  %12.5e]",  jj, forceNorms1[jj], forceNorms2[jj] );
                     summary << _getLine( tab, "Error", value );
                } else {
    
                    double diff, sum;
                    if( useForces ){
                        double delta   = std::sqrt( (forces1[jj][0]-forces2[jj][0])*(forces1[jj][0]-forces2[jj][0]) +
                                                    (forces1[jj][1]-forces2[jj][1])*(forces1[jj][1]-forces2[jj][1]) +
                                                    (forces1[jj][2]-forces2[jj][2])*(forces1[jj][2]-forces2[jj][2]) );
    
                        sum            = 0.5*( fabs( forceNorms1[jj] ) + fabs( forceNorms2[jj] ) );
                        diff           = delta > 0.0 ? delta/sum : 0.0;
                    } else {
                        sum            = 0.5*( fabs( forceNorms1[jj] ) + fabs( forceNorms2[jj] ) );
                        diff           = sum > 0.0 ? fabs( forceNorms1[jj] - forceNorms2[jj] )/sum : 0.0;
                    }
        
                    (void) LOCAL_SPRINTF12( value, "%11.5e  at index %6d  norms: [%11.5e  %11.5e] forces: [%12.5e  %12.5e %12.5e] [%12.5e  %12.5e %12.5e]",
                                            diff, jj, forceNorms1[jj], forceNorms2[jj],
                                            forces1[jj][0], forces1[jj][1], forces1[jj][2],
                                            forces2[jj][0], forces2[jj][1], forces2[jj][2] );
                    summary << _getLine( tab, "Error", value );
                    if( kk == (maxMissesToPrint-1) ){
                        value[0] = '\0';
                        summary << _getLine( tab, "No more errors will be printed: call ValidateOpenMMForces::setMaxErrorsToPrint() to modifiy this limit.", value );
                    }
                }
            }
        }

        (void) LOCAL_SPRINTF0( value, "\n" );
        summary << _getLine( tab, " ", value );
        if( inconsistentIndices.size() > 0 ){
            (void) LOCAL_SPRINTF( value, "%u\n", static_cast<unsigned int>(inconsistentIndices.size()) );
            summary << _getLine( tab, "Total errors", value );
        }
        if( forceValidationResults[ii]->nansDetected() ){
            value[0] = '\0';
            summary << _getLine( tab, "Nans detected.", value );
        }
    }

    return summary.str();
}
