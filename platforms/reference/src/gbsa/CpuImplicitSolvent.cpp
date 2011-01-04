
/* Portions copyright (c) 2006 Stanford University and Simbios.
 * Contributors: Pande Group
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "../SimTKUtilities/SimTKOpenMMCommon.h"
#include "../SimTKUtilities/SimTKOpenMMLog.h"
#include "../SimTKUtilities/SimTKOpenMMUtilities.h"
#include "CpuImplicitSolvent.h"

#include <cstdio>

using namespace OpenMM;
using namespace std;

//#define UseGromacsMalloc 1

// Replacement new/delete w/ Gromac's smalloc() and sfree()

// static data member created by call to cpuSetImplicitSolventParameters() 
// stores parameter settings, ...
// used to calculate GBSA forces/energy

CpuImplicitSolvent* CpuImplicitSolvent::_cpuImplicitSolvent                     = NULL;

// info file related-stuff

std::string CpuImplicitSolvent::_defaultInfoFileName                            = std::string( "CpuImplicitSolventInfo" );

// key for info file: base file name 

const std::string CpuImplicitSolvent::CpuImplicitSolventBaseFileName            = std::string( "CpuImplicitSolventBaseFileName" );

// key for info file: file generation frequency 

const std::string CpuImplicitSolvent::CpuImplicitSolventFileGenerationFrequency = std::string( "CpuImplicitSolventFileGenerationFrequency" );


/**---------------------------------------------------------------------------------------

   CpuImplicitSolvent constructor

   @param implicitSolventParameters      implicitSolventParameters object
   
   --------------------------------------------------------------------------------------- */

CpuImplicitSolvent::CpuImplicitSolvent( ImplicitSolventParameters* implicitSolventParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::CpuImplicitSolvent(2)";

   // ---------------------------------------------------------------------------------------

   _initializeDataMembers( );
   _implicitSolventParameters = implicitSolventParameters;

}

/**---------------------------------------------------------------------------------------

   CpuImplicitSolvent destructor

   --------------------------------------------------------------------------------------- */

CpuImplicitSolvent::~CpuImplicitSolvent( ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::~CpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   delete _implicitSolventParameters;
}

/**---------------------------------------------------------------------------------------

   Initialize data members

   --------------------------------------------------------------------------------------- */

void CpuImplicitSolvent::_initializeDataMembers( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::_initializeDataMembers";

   // ---------------------------------------------------------------------------------------

   _includeAceApproximation             = 0;

   _forceConversionFactor               = (RealOpenMM) 1.0;
   _energyConversionFactor              = (RealOpenMM) 1.0;

   _forceCallIndex                      = 0;

   _implicitSolventEnergy               = (RealOpenMM) 0.0;

   _baseFileName                        = SimTKOpenMMCommon::NotSet;
   _outputFileFrequency                 = 1;
}

/**---------------------------------------------------------------------------------------

   Return number of atoms

   @return number of atoms

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::getNumberOfAtoms( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getNumberOfAtoms";

   // ---------------------------------------------------------------------------------------

   return _implicitSolventParameters->getNumberOfAtoms();

}

/**---------------------------------------------------------------------------------------

   Return energy 

   @return energy

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuImplicitSolvent::getEnergy( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getEnergy";

   // ---------------------------------------------------------------------------------------

   return _implicitSolventEnergy;

}

/**---------------------------------------------------------------------------------------

   Delete static _cpuImplicitSolvent object if set

   @return SimTKOpenMMCommon::DefaultReturn if _cpuImplicitSolvent was set; 
           otherwise return SimTKOpenMMCommon::ErrorReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::deleteCpuImplicitSolvent( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::deleteCpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   if( _cpuImplicitSolvent != NULL ){
      delete _cpuImplicitSolvent;
      _cpuImplicitSolvent = NULL;
      return SimTKOpenMMCommon::DefaultReturn;
   } else {
      return SimTKOpenMMCommon::ErrorReturn;
   }
}

/**---------------------------------------------------------------------------------------

   Set static member _cpuImplicitSolvent

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setCpuImplicitSolvent( CpuImplicitSolvent* cpuImplicitSolvent ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setCpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   _cpuImplicitSolvent = cpuImplicitSolvent;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get static member cpuImplicitSolvent

   @return static member cpuImplicitSolvent

   --------------------------------------------------------------------------------------- */

CpuImplicitSolvent* CpuImplicitSolvent::getCpuImplicitSolvent( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getCpuImplicitSolvent";

   // ---------------------------------------------------------------------------------------

   return _cpuImplicitSolvent;

}
/**---------------------------------------------------------------------------------------

   Set energy 

   @param energy

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setEnergy( RealOpenMM energy ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setEnergy";

   // ---------------------------------------------------------------------------------------

   _implicitSolventEnergy = energy;
   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Return ImplicitSolventParameters

   @return ImplicitSolventParameters

   --------------------------------------------------------------------------------------- */

ImplicitSolventParameters* CpuImplicitSolvent::getImplicitSolventParameters( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getImplicitSolventParameters";

   // ---------------------------------------------------------------------------------------

   return _implicitSolventParameters;

}

/**---------------------------------------------------------------------------------------

   Set ImplicitSolventParameters

   @param ImplicitSolventParameters

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setImplicitSolventParameters( ImplicitSolventParameters* implicitSolventParameters ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setImplicitSolventParameters";

   // ---------------------------------------------------------------------------------------

   _implicitSolventParameters = implicitSolventParameters;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return flag signalling whether AceApproximation for nonpolar term is to be included

   @return flag

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::includeAceApproximation( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::includeAceApproximation";

   // ---------------------------------------------------------------------------------------

   return _includeAceApproximation;

}

/**---------------------------------------------------------------------------------------

   Set flag indicating whether AceApproximation is to be included

   @param includeAceApproximation new includeAceApproximation value

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setIncludeAceApproximation( int includeAceApproximation ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setImplicitSolventParameters";

   // ---------------------------------------------------------------------------------------

   _includeAceApproximation = includeAceApproximation;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return ForceConversionFactor for units

   @return ForceConversionFactor

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuImplicitSolvent::getForceConversionFactor( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getForceConversionFactor";

   // ---------------------------------------------------------------------------------------

   return _forceConversionFactor;

}

/**---------------------------------------------------------------------------------------

   Set ForceConversionFactor

   @param ForceConversionFactor (units conversion)

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setForceConversionFactor( RealOpenMM forceConversionFactor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setForceConversionFactor";

   // ---------------------------------------------------------------------------------------

   _forceConversionFactor = forceConversionFactor;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return EnergyConversionFactor for units

   @return EnergyConversionFactor

   --------------------------------------------------------------------------------------- */

RealOpenMM CpuImplicitSolvent::getEnergyConversionFactor( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getEnergyConversionFactor";

   // ---------------------------------------------------------------------------------------

   return _energyConversionFactor;

}

/**---------------------------------------------------------------------------------------

   Set EnergyConversionFactor

   @param EnergyConversionFactor (units conversion)

   @return SimTKOpenMMCommon::DefaultReturn 

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setEnergyConversionFactor( RealOpenMM energyConversionFactor ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::setEnergyConversionFactor";

   // ---------------------------------------------------------------------------------------

   _energyConversionFactor = energyConversionFactor;

   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Return ForceCallIndex -- number of times forces have been calculated

   @return ForceCallIndex

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::getForceCallIndex( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getForceCallIndex";

   // ---------------------------------------------------------------------------------------

   return _forceCallIndex;

}

/**---------------------------------------------------------------------------------------

   Increment ForceCallIndex

   @return incremented forceCallIndex

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::incrementForceCallIndex( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::incrementForceCallIndex";

   // ---------------------------------------------------------------------------------------

   _forceCallIndex++;

   return _forceCallIndex;

}

/**---------------------------------------------------------------------------------------

   Return bornForce, a work array of size _implicitSolventParameters->getNumberOfAtoms()*sizeof( RealOpenMM )
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuImplicitSolvent::getBornForce( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getImplicitSolventBornForce";

   // ---------------------------------------------------------------------------------------

   if( _bornForce.size() == 0 ){
      _bornForce.resize(_implicitSolventParameters->getNumberOfAtoms());
   }
   return _bornForce;

}

/**---------------------------------------------------------------------------------------

   Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if it is not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuImplicitSolvent::getBornRadii( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getBornRadii";

   // ---------------------------------------------------------------------------------------

   if( _bornRadii.size() == 0 ){
      _bornRadii.resize(_implicitSolventParameters->getNumberOfAtoms(), 0.0);
   }
   return _bornRadii;
}

/**---------------------------------------------------------------------------------------

   Return Born radii: size = _implicitSolventParameters->getNumberOfAtoms()

   @return array

   --------------------------------------------------------------------------------------- */

const vector<RealOpenMM>& CpuImplicitSolvent::getBornRadiiConst( void ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getBornRadii";

   // ---------------------------------------------------------------------------------------

   return _bornRadii;
}

/**---------------------------------------------------------------------------------------

   Return Born radii temp work array of size=_implicitSolventParameters->getNumberOfAtoms()
   On first call, memory for array is allocated if not set

   @return array

   --------------------------------------------------------------------------------------- */

vector<RealOpenMM>& CpuImplicitSolvent::getBornRadiiTemp( void ){

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getImplicitSolventBornRadiiTemp";

   // ---------------------------------------------------------------------------------------

   if( _tempBornRadii.size() == 0 ){
      _tempBornRadii.resize(_implicitSolventParameters->getNumberOfAtoms(), 0.0);
   }
   return _tempBornRadii;
}

/**---------------------------------------------------------------------------------------

   Compute Born radii

   @param atomCoordinates     atomic coordinates
   @param bornRadii           output array of Born radii
   @param obcChain            output array of Obc chain derivatives

   @return SimTKOpenMMCommon::DefaultReturn or SimTKOpenMMCommon::ErrorReturn 
           if problems encountered

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeBornRadii( vector<RealVec>& atomCoordinates, vector<RealOpenMM>& bornRadii ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuImplicitSolvent::computeBornRadii";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << methodName;
   message << " Error: calling from base class.";
   SimTKOpenMMLog::printError( message );
   return SimTKOpenMMCommon::ErrorReturn;

}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces

   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces (output)
   @param updateBornRadii     if set, then Born radii are updated for current configuration; 
                              otherwise radii correspond to configuration from previous iteration


   @return SimTKOpenMMCommon::DefaultReturn; abort if cpuImplicitSolvent is not set

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeImplicitSolventForces( vector<RealVec>& atomCoordinates,
                                                      const RealOpenMM* partialCharges,
                                                      vector<RealVec>& forces, int updateBornRadii ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuImplicitSolvent::computeImplicitSolventForces";

   // ---------------------------------------------------------------------------------------

   int callId = incrementForceCallIndex();

   // check parameters have been initialized

   ImplicitSolventParameters* implicitSolventParameters =  getImplicitSolventParameters();
   if( !implicitSolventParameters ){
      std::stringstream message;
      message << methodName;
      message << " implicitSolventParameters has not been initialized!";
      SimTKOpenMMLog::printError( message );
      return SimTKOpenMMCommon::ErrorReturn; 
   }

   // check to see if Born radii have been previously calculated
   // if not, then calculate;
   // logic here assumes that the radii are intitialized to zero 
   // and then once computed, always greater than zero.

   // after first iteration Born radii are updated in force calculation (computeBornEnergyForces())
   // unless updateBornRadii is set

   vector<RealOpenMM>& bornRadii = getBornRadii();
   if( updateBornRadii || bornRadii[0] < (RealOpenMM) 0.0001 || callId == 1 ){
      computeBornRadii( atomCoordinates, bornRadii );
   }

   // compute forces

   computeBornEnergyForces( getBornRadii(), atomCoordinates, partialCharges, forces );

   return SimTKOpenMMCommon::DefaultReturn; 
}

/**---------------------------------------------------------------------------------------

   Get Born energy and forces based on J. Phys. Chem. A V101 No 16, p. 3005 (Simbios)

   @param bornRadii           Born radii -- optional; if NULL, then ImplicitSolventParameters 
                              entry is used
   @param atomCoordinates     atomic coordinates
   @param partialCharges      partial charges
   @param forces              forces

   @return SimTKOpenMMCommon::ErrorReturn since the call should be implemented 
           in a derived class

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeBornEnergyForces( vector<RealOpenMM>& bornRadii,
                                                 vector<RealVec>& atomCoordinates,
                                                 const RealOpenMM* partialCharges,
                                                 vector<RealVec>& forces ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName = "\nCpuImplicitSolvent::computeBornEnergyForces";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   message << methodName;
   message << " Error: calling from base class.";
   SimTKOpenMMLog::printError( message );
   return SimTKOpenMMCommon::ErrorReturn; 

}

/**---------------------------------------------------------------------------------------

   Get nonpolar solvation force constribution via ACE approximation

   @param implicitSolventParameters parameters
   @param vdwRadii                  Vdw radii
   @param bornRadii                 Born radii
   @param energy                    energy (output): value is incremented from input value 
   @param forces                    forces: values are incremented from input values

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::computeAceNonPolarForce( const ImplicitSolventParameters* implicitSolventParameters,
                                                 const vector<RealOpenMM>& bornRadii, RealOpenMM* energy,
                                                 vector<RealOpenMM>& forces ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::computeAceNonPolarForce";

   static const RealOpenMM minusSix = -6.0;

   // ---------------------------------------------------------------------------------------

   // compute the nonpolar solvation via ACE approximation

   const RealOpenMM probeRadius          = implicitSolventParameters->getProbeRadius();
   const RealOpenMM surfaceAreaFactor    = implicitSolventParameters->getPi4Asolv();

   const RealOpenMM* atomicRadii         = implicitSolventParameters->getAtomicRadii();
   int numberOfAtoms                     = implicitSolventParameters->getNumberOfAtoms();

   // 1 + 1 + pow + 3 + 1 + 2 FLOP

   // the original ACE equation is based on Eq.2 of

   // M. Schaefer, C. Bartels and M. Karplus, "Solution Conformations
   // and Thermodynamics of Structured Peptides: Molecular Dynamics
   // Simulation with an Implicit Solvation Model", J. Mol. Biol.,
   // 284, 835-848 (1998)  (ACE Method)

   // The original equation includes the factor (atomicRadii[atomI]/bornRadii[atomI]) to the first power,
   // whereas here the ratio is raised to the sixth power: (atomicRadii[atomI]/bornRadii[atomI])**6

   // This modification was made by Jay Ponder who observed it gave better correlations w/
   // observed values. He did not think it was important enough to write up, so there is
   // no paper to cite.

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){
      if( bornRadii[atomI] > 0.0 ){
         RealOpenMM r            = atomicRadii[atomI] + probeRadius;
         RealOpenMM ratio6       = POW( atomicRadii[atomI]/bornRadii[atomI], (RealOpenMM) 6.0 );
         RealOpenMM saTerm       = surfaceAreaFactor*r*r*ratio6;
         *energy                += saTerm;
         forces[atomI]          += minusSix*saTerm/bornRadii[atomI]; 
      }
   }

   return SimTKOpenMMCommon::DefaultReturn; 

}

/**---------------------------------------------------------------------------------------

   Get BaseFileName

   @return    baseFileName

   --------------------------------------------------------------------------------------- */

const std::string& CpuImplicitSolvent::getBaseFileName( void ) const {

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::getBaseFileName";

// ---------------------------------------------------------------------------------------

   return _baseFileName;

}

/**---------------------------------------------------------------------------------------

   Set BaseFileName

   @param    input baseFileName

   @return   SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setBaseFileName( const std::string& baseFileName ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::setBaseFileName";

// ---------------------------------------------------------------------------------------

   _baseFileName = baseFileName;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Set OutputFileFrequency

   @param    input outputFileFrequency

   @return   SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::setOutputFileFrequency( int outputFileFrequency ){

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::setOutputFileFrequency";

// ---------------------------------------------------------------------------------------

   _outputFileFrequency = outputFileFrequency;
   return SimTKOpenMMCommon::DefaultReturn;

}

/**---------------------------------------------------------------------------------------

   Get OutputFileFrequency

   @return   output file frequency

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::getOutputFileFrequency( void ) const {

// ---------------------------------------------------------------------------------------

// static const std::string methodName = "\nCpuImplicitSolvent::getOutputFileFrequency";

// ---------------------------------------------------------------------------------------

   return _outputFileFrequency;

}

/**---------------------------------------------------------------------------------------

   Read info file

   @param infoFileName       file name to read

   @return SimTKOpenMMCommon::DefaultReturn if ok; else SimTKOpenMMCommon::ErrorReturn; if major 
           problem, program will exit

   --------------------------------------------------------------------------------------- */

int CpuImplicitSolvent::readInfoFile( const std::string infoFileName ){

   // ---------------------------------------------------------------------------------------

   const int bufferSize = 2048;
   char buffer[bufferSize];

   static const std::string methodName = "\nCpuImplicitSolvent::readInfoFile";

   // ---------------------------------------------------------------------------------------

   FILE* infoFile = NULL;
#ifdef _MSC_VER
   fopen_s( &infoFile, infoFileName.c_str(), "r" );
#else
   infoFile = fopen( infoFileName.c_str(), "r" );
#endif

   if( infoFile == NULL ){

#ifdef _MSC_VER
      fopen_s( &infoFile, _defaultInfoFileName.c_str(), "r" );
#else
      infoFile = fopen( _defaultInfoFileName.c_str(), "r" );
#endif

      if( infoFile == NULL ){
         std::stringstream message;
         message << methodName << " fopen failed on info file=<" << _defaultInfoFileName << "> and file=<" << infoFileName << ">";
         SimTKOpenMMLog::printMessage( message );
         return SimTKOpenMMCommon::ErrorReturn;
      } else {
         std::stringstream message;
         message << methodName << " opened info file=<" << _defaultInfoFileName << ">.";
         SimTKOpenMMLog::printMessage( message );
      }
   } else {
      std::stringstream message;
      message << methodName << " opened info file=<" << infoFileName << ">.";
      SimTKOpenMMLog::printMessage( message );
   }

   /* Sample file:

   */

   int errors = 0;
   while( fgets( buffer, bufferSize, infoFile ) ){

      std::stringstream message;
      std::string fileName   = std::string( );
      size_t bufferLen       = strlen( buffer );
      if( bufferLen > 0 ){
         buffer[bufferLen-1] = '\0'; 
      }
      
      message << "\n<" << buffer << ">";

//fprintf( amoebaLog->getLogFile(), "\nreadInfoFile <%s>", buffer );
//fflush( amoebaLog->getLogFile() );

      StringVector tokens;
      SimTKOpenMMUtilities::tokenizeString( buffer, tokens );

      bool done = false;
      for( StringVectorI ii = tokens.begin(); ii != tokens.end() && !done; ii++ ){

         std::string token        = *ii;
         size_t tokenLength       = token.length();
         bool recognized          = false;

         // skip comments and blank lines

         if( tokenLength < 2 || !token.compare( "#" ) ){
            done       = true;
            recognized = true;

         // base file name

         } else if( !token.compare( CpuImplicitSolvent::CpuImplicitSolventBaseFileName ) ){
            recognized = true;
            setBaseFileName( *(++ii) );

         // file frequency

         } else if( !token.compare( CpuImplicitSolvent::CpuImplicitSolventFileGenerationFrequency ) ){
            recognized        = true;
            std::string value = *(++ii);
            if( SimTKOpenMMUtilities::isValidInteger( value ) ){ 
               setOutputFileFrequency( atoi( value.c_str() ) );
            } else {
               message << "\nToken=<" << token << "> is not a valid integer for key=" << CpuImplicitSolvent::CpuImplicitSolventFileGenerationFrequency;
            }

         // keys used by objects other than 'CpuImplicitSolvent'

         } else {
/*
            recognized      = true;
            std::string key = *ii;
            ii++;
            if( ii != tokens.end() ){
               _inputArguments[key] = *ii;
            } else {
               _inputArguments[key] = CpuImplicitSolventCommon::Comment;
            }
*/
         }
         if( !recognized ){
            message << "\nToken=<" << token << "> not recognized.";
         }
      }
      SimTKOpenMMLog::printMessage( message );
   }

   (void) fclose( infoFile );

   // report if errors

   if( errors ){
      std::stringstream message;
      message << "\nErrors -- aborting.\n";
      SimTKOpenMMLog::printMessage( message );
      exit(-1);
   } else {
      std::stringstream message;
      message << "\nNo errors detected parsing info file " << infoFileName << "\n";
/*
      message << "\nNumber of hash arguments=" << _amoebaInputArguments.size() << " <key,value> pairs listed below:\n";
      for( StringStringMapCI ii = _amoebaInputArguments.begin(); ii != _amoebaInputArguments.end(); ii++ ){
         message << "\n   <" << (*ii).first <<  ">=<" << (*ii).second << ">";
      }
*/
      message << "\n\n";
      SimTKOpenMMLog::printMessage( message );
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------
      
   Get string w/ state 
   
   @param title               title (optional)
      @return string containing state
      
   --------------------------------------------------------------------------------------- */

std::string CpuImplicitSolvent::getStateString( const char* title ) const {

   // ---------------------------------------------------------------------------------------

   // static const char* methodName = "\nCpuImplicitSolvent::getStateString";

   // ---------------------------------------------------------------------------------------

   std::stringstream message;
   if( title ){
      message << title;
   }

   message << "\nImplicitSolvent info:";
   message << "\nForce conversion=" << getForceConversionFactor() << " Energy conversion=" << getEnergyConversionFactor();
   message << "\nInclude ACE approximation=";
   if( includeAceApproximation() ){
      message << "Y";
   } else {
      message << "N";
   }

   // get parameters

   message << getImplicitSolventParameters()->getStateString( NULL );

   return message.str();
}

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */

/*
void* CpuImplicitSolvent::operator new( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

//   (void) fprintf( stdout, "\nCpuImplicitSolvent new called -- size=%u", size );
//   (void) fflush( stdout );

   return ptr;
} */

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */

/*
void CpuImplicitSolvent::operator delete( void *ptr ){

   // (void) fprintf( stdout, "\nCpuImplicitSolvent delete called." );
   // (void) fflush( stdout );

   sfree( ptr ); 
} */

/* ---------------------------------------------------------------------------------------

   Override C++ new w/ Gromac's smalloc/sfree (Simbios)

   @param size						bytes to allocate

   @return ptr to allocated memory

   --------------------------------------------------------------------------------------- */
/*
void* CpuImplicitSolvent::operator new[]( size_t size ){

   void *ptr;
   smalloc(ptr, (int) size); 

   // (void) fprintf( stdout, "\nCpuImplicitSolvent new[] called -- size=%u", size );
   // (void) fflush( stdout );

   return ptr;
} */

/* ---------------------------------------------------------------------------------------

   Override C++ delete w/ Gromac's sfree (Simbios)

   @param ptr						ptr to block to free

   --------------------------------------------------------------------------------------- */
/*
void CpuImplicitSolvent::operator delete[]( void *ptr ){

   // (void) fprintf( stdout, "\nCpuImplicitSolvent delete[] called." );
   // (void) fflush( stdout );

   sfree( ptr ); 
}
*/
