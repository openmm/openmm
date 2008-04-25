
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

#include "cpuGrycukInterface.h"

#include "SimTKOpenMMLog.h"
#include "SimTKOpenMMUtilities.h"
// #include "SimTKOpenMMGromacsUtilities.h"
#include "GrycukParameters.h"
#include "CpuGrycuk.h"

/**---------------------------------------------------------------------------------------

	Setup for Grycuk calculations from Gromacs

   @param numberOfAtoms            number of atoms

   @param grycukScaleFactors          array of Grycuk scale factors (one entry each atom)

   @param atomicRadii              atomic radii in Angstrom (one entry each atom)

   @param includeAceApproximation  if true, then include nonpolar 
                                   ACE term in calculations

   @param soluteDielectric         solute dielectric

   @param solventDielectric        solvent dielectric

   @param log                      log reference -- if NULL, then errors/warnings
                                   output to stderr

   The method creates a CpuGrycuk instance -- currently the Grycuk type II model is the
   default (see paper). If the Grycuk type I model is desired change

      GrycukParameters* grycukParameters  = new GrycukParameters( numberOfAtoms, GrycukParameters::GrycukTypeII );
   to
      GrycukParameters* grycukParameters  = new GrycukParameters( numberOfAtoms, GrycukParameters::GrycukTypeI  );

   The created object is a static member of the class CpuGrycuk; 
   when the force routine, cpuCalculateGrycukForces(), is called, 
   the static object is used to compute the forces and energy 

   @return 0

   --------------------------------------------------------------------------------------- */

extern "C" int 
cpuSetGrycukParameters( int numberOfAtoms, RealOpenMM* atomicRadii, RealOpenMM* grycukScaleFactors,
                        int includeAceApproximation,
                        RealOpenMM soluteDielectric, RealOpenMM solventDielectric, FILE* log ){

   // ---------------------------------------------------------------------------------------

   static const char* methodName   = "\ncpuSetGrycukParameters: ";

   // ---------------------------------------------------------------------------------------
   
   // set log file if not NULL

   if( log ){
      SimTKOpenMMLog::setSimTKOpenMMLog( log );
   }

   // set Grycuk parameters

   GrycukParameters* grycukParameters  = new GrycukParameters( numberOfAtoms );
   grycukParameters->setScaledRadiusFactors( grycukScaleFactors );
   grycukParameters->setAtomicRadii( atomicRadii, SimTKOpenMMCommon::KcalAngUnits );

   // dielectric constants

   grycukParameters->setSolventDielectric( solventDielectric );
   grycukParameters->setSoluteDielectric( soluteDielectric );

   // ---------------------------------------------------------------------------------------

   // create CpuGrycuk instance that will calculate forces
  
   CpuGrycuk* cpuGrycuk = new CpuGrycuk( grycukParameters );

   // set static member for subsequent calls to calculate forces/energy 

   CpuImplicitSolvent::setCpuImplicitSolvent( cpuGrycuk );

   // set base file name, ...

   //cpuGrycuk->readInfoFile( "CpuImplicitSolventInfo" );

   // include/do not include ACE approximation (nonpolar solvation)

   cpuGrycuk->setIncludeAceApproximation( includeAceApproximation );

   // ---------------------------------------------------------------------------------------

   // diagnostics
 
   if( log ){
      std::string state = cpuGrycuk->getStateString( methodName );
      (void) fprintf( log, "\n%s\nDone w/ setup\n", state.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   return 0;
}

/**---------------------------------------------------------------------------------------

   Get Grycuk scale factors given masses

   @param numberOfAtoms number of atoms
   @param masses        input masses 
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getGrycukScaleFactorsGivenAtomMasses( int numberOfAtoms, const RealOpenMM* masses,
                                                      RealOpenMM* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetGrycukScaleFactorsGivenAtomMasses";

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      double scaleFactor = 0.8;
      RealOpenMM mass    = masses[atomI];

      if ( mass < 1.2 && mass >= 1.0 ){        // hydrogen
         scaleFactor  = 0.85; 
      } else if( mass > 11.8 && mass < 12.2 ){ // carbon
         scaleFactor  = 0.72; 
      } else if( mass > 14.0 && mass < 15.0 ){ // nitrogen
         scaleFactor  = 0.79;
      } else if( mass > 15.5 && mass < 16.5 ){ // oxygen
         scaleFactor  = 0.85; 
      } else if( mass > 31.5 && mass < 32.5 ){ // sulphur
         scaleFactor  = 0.96;
      } else if( mass > 29.5 && mass < 30.5 ){ // phosphorus
         scaleFactor  = 0.86;
      } else {
         std::stringstream message;
         message << methodName;
         message << " Warning: mass for atom " << atomI << " mass=" << mass << "> not recognized.";
         SimTKOpenMMLog::printMessage( message );
      }

      scaleFactors[atomI] = (RealOpenMM) scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

/**---------------------------------------------------------------------------------------

   Get Grycuk scale factors given atomic numbers

   @param numberOfAtoms number of atoms
   @param atomicNumber  input atomic number for each atom
   @param scaleFactors  output atomic numbers

   @return SimTKOpenMMCommon::DefaultReturn

   --------------------------------------------------------------------------------------- */

extern "C" int getGrycukScaleFactors( int numberOfAtoms, const int* atomicNumber, RealOpenMM* scaleFactors ){

   // ---------------------------------------------------------------------------------------

   static const std::string methodName = "\ngetGrycukScaleFactors";

   // ---------------------------------------------------------------------------------------

   for( int atomI = 0; atomI < numberOfAtoms; atomI++ ){

      double scaleFactor;
      switch( atomicNumber[atomI] ){

         case 1: // hydrogen

            scaleFactor  = 0.85; 
            break;

         case 6: // carbon

            scaleFactor  = 0.72; 
            break;

         case 7: // nitrogen

            scaleFactor  = 0.79;
            break;

         case 8: // oxygen

            scaleFactor  = 0.85;
            break;

         case 15: // phosphorus

            scaleFactor  = 0.86;
            break;

         case 16: // sulphur

            scaleFactor  = 0.85;
            break;

         default:

            scaleFactor = 0.8;

            std::stringstream message;
            message << methodName;
            message << " Warning: atom number=" << atomicNumber[atomI] << " for atom " << atomI << "> not handled -- useing default value.";
            SimTKOpenMMLog::printMessage( message );
            break;
      }

      scaleFactors[atomI] = (RealOpenMM) scaleFactor;
   }

   return SimTKOpenMMCommon::DefaultReturn;
}

