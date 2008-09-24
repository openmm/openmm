/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include <cmath>
#include <limits>
#include "OpenMMException.h"
#include <sstream>

#include "BrookStreamImpl.h"
#include "BrookCalcGBSAOBCForceFieldKernel.h"
#include "force.h"
#include "kgbsa.h"
#include "kforce.h"
#include "math.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookCalcGBSAOBCForceFieldKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 *
 */

BrookCalcGBSAOBCForceFieldKernel::BrookCalcGBSAOBCForceFieldKernel( std::string name, const Platform& platform ) :
                     CalcGBSAOBCForceFieldKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcGBSAOBCForceFieldKernel::BrookCalcGBSAOBCForceFieldKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   _numberOfAtoms                    = 0;
   _brookGbsa                        = NULL;
   _log                              = NULL;

   const BrookPlatform brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
}   

/** 
 * BrookCalcGBSAOBCForceFieldKernel destructor
 * 
 */

BrookCalcGBSAOBCForceFieldKernel::~BrookCalcGBSAOBCForceFieldKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcGBSAOBCForceFieldKernel::BrookCalcGBSAOBCForceFieldKernel";
   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   delete _brookGbsa;
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCalcGBSAOBCForceFieldKernel::getLog( void ) const {
   return _log;
}

/** 
 * Set log file reference
 * 
 * @param  log file reference
 *
 * @return  DefaultReturnValue
 *
 */

int BrookCalcGBSAOBCForceFieldKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param atomParameters            vector containing atom index, charge, radius, scalingFactor
 * @param solventDielectric         solvent dielectric
 * @param soluteDielectric          solute dielectric
 *
 */

void BrookCalcGBSAOBCForceFieldKernel::initialize( const std::vector<std::vector<double> >& atomParameters, 
                                                   double solventDielectric, double soluteDielectric ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcGBSAOBCForceFieldKernel::initialize";

// ---------------------------------------------------------------------------------------

   FILE* log                 = getLog();
   _numberOfAtoms            = (int) atomParameters.size();

   // ---------------------------------------------------------------------------------------

   // bonded

   if( _brookGbsa ){
      delete _brookGbsa;
   }
   _brookGbsa              = new BrookGbsa();
   _brookGbsa->setLog( log );
    
   _brookGbsa->setup( atomParameters, solventDielectric, soluteDielectric, getPlatform() );

   if( log ){
      std::string contents = _brookGbsa->getContentsString( ); 
      (void) fprintf( log, "%s brookGbsa::contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------
    
}

/** 
 * Compute forces given atom coordinates
 * 
 * @param positions                 atom coordinates
 * @param forces                    output forces
 *
 */

void BrookCalcGBSAOBCForceFieldKernel::executeForces( const Stream& positions, Stream& forces ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName   = "BrookCalcGBSAOBCForceFieldKernel::executeForces";
   static const int PrintOn              = 0; 

// ---------------------------------------------------------------------------------------

   // OBC

   const BrookStreamImpl& positionStreamC              = dynamic_cast<const BrookStreamImpl&> (positions.getImpl());
   BrookStreamImpl& positionStream                     = const_cast<BrookStreamImpl&>         (positionStreamC);
   BrookStreamImpl& forceStream                        = dynamic_cast<BrookStreamImpl&>       (forces.getImpl());

   float includeAce                                    = (float) (_brookGbsa->includeAce());
   BrookFloatStreamInternal**  gbsaForceStreams        = _brookGbsa->getForceStreams();

   // calculate Born radii first time thru and initialize on board

   if( !_brookGbsa->haveBornRadiiBeenInitialized() ){
      _brookGbsa->calculateBornRadii( positions );  
   }

   // first major loop

   kObcLoop1( (float) _brookGbsa->getNumberOfAtoms(),
              (float) _brookGbsa->getAtomSizeCeiling(),
              (float) _brookGbsa->getDuplicationFactor(),
              (float) _brookGbsa->getAtomStreamWidth( ),
              (float) _brookGbsa->getPartialForceStreamWidth( ),
              _brookGbsa->getSoluteDielectric(),
              _brookGbsa->getSolventDielectric(),
              includeAce,

              positionStream.getBrookStream(),

              _brookGbsa->getObcBornRadii()->getBrookStream(),
              _brookGbsa->getObcAtomicRadii()->getBrookStream(),

              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream()
            );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( 1 && PrintOn && getLog() ){

      (void) fprintf( getLog(), "\nPost kObcLoop1: atms=%d ceil=%d dup=%d atomStrW=%3d prtlF=%3d diel=%.3f %.3f ACE=%.1f\n",
                      _brookGbsa->getNumberOfAtoms(),
                      _brookGbsa->getAtomSizeCeiling(),
                      _brookGbsa->getDuplicationFactor(),
                      _brookGbsa->getAtomStreamWidth( ),
                      _brookGbsa->getPartialForceStreamWidth( ),
                      _brookGbsa->getSoluteDielectric(),
                      _brookGbsa->getSolventDielectric(), includeAce );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      _brookGbsa->getObcBornRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nAtomR\n" );
      _brookGbsa->getObcAtomicRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStreams output\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

   }

// ---------------------------------------------------------------------------------------

   // gather for first loop

   kPostObcLoop1_nobranch(
              (float) _brookGbsa->getDuplicationFactor(),
              (float) _brookGbsa->getAtomStreamWidth( ),
              (float) _brookGbsa->getPartialForceStreamWidth( ),
              (float) _brookGbsa->getNumberOfAtoms(),
              (float) _brookGbsa->getAtomSizeCeiling(),
              (float) _brookGbsa->getInnerLoopUnroll(),
              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream(),
              _brookGbsa->getObcChain()->getBrookStream(),
              _brookGbsa->getObcBornRadii()->getBrookStream(),
              _brookGbsa->getObcIntermediateForce()->getBrookStream(),
              _brookGbsa->getObcBornRadii2()->getBrookStream() );
 
// ---------------------------------------------------------------------------------------

   // diagnostics

   if( PrintOn && getLog()){

      (void) fprintf( getLog(), "\nPost kPostObcLoop1_nobranch: dup=%d aStrW=%d pStrW=%d no.atms=%3d ceil=%3d Unroll=%1d\n",
                      _brookGbsa->getDuplicationFactor(),
                      _brookGbsa->getAtomStreamWidth( ),
                      _brookGbsa->getPartialForceStreamWidth( ),
                      _brookGbsa->getNumberOfAtoms(),
                      _brookGbsa->getAtomSizeCeiling(),
                      _brookGbsa->getInnerLoopUnroll() );

      (void) fprintf( getLog(), "\nForceStreams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

      (void) fprintf( getLog(), "\nObcChain\n" );
      _brookGbsa->getObcChain()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      _brookGbsa->getObcBornRadii()->printToFile( getLog() );

      // output

      (void) fprintf( getLog(), "\nObcIntermediateForce output\n" );
      _brookGbsa->getObcIntermediateForce()->printToFile( getLog() );

      // output

      (void) fprintf( getLog(), "\nObcBornRadii2 output\n" );
      _brookGbsa->getObcBornRadii2()->printToFile( getLog() );

   }

// ---------------------------------------------------------------------------------------

   // second major loop

   kObcLoop2( (float) _brookGbsa->getNumberOfAtoms(),
              (float) _brookGbsa->getAtomSizeCeiling(),
              (float) _brookGbsa->getDuplicationFactor(),
              (float) _brookGbsa->getAtomStreamWidth( ),
              (float) _brookGbsa->getPartialForceStreamWidth( ),
              positionStream.getBrookStream(),
              _brookGbsa->getObcScaledAtomicRadii()->getBrookStream(),
              _brookGbsa->getObcBornRadii2()->getBrookStream(),
              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream()
            );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "\nPost kObcLoop2: no.atms=%5d ceil=%3d dup=%3d strW=%3d pStrW=%3d\n",
                      _brookGbsa->getNumberOfAtoms(),
                      _brookGbsa->getAtomSizeCeiling(),
                      _brookGbsa->getDuplicationFactor(),
                      _brookGbsa->getAtomStreamWidth( ),
                      _brookGbsa->getPartialForceStreamWidth( ) );

      BrookStreamInternal* brookStreamInternalPos  = positionStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nPositionStream\n" );
      brookStreamInternalPos->printToFile( getLog() );

      (void) fprintf( getLog(), "\nObcScaledAtomicRadii\n" );
      _brookGbsa->getObcScaledAtomicRadii()->printToFile( getLog() );

      (void) fprintf( getLog(), "\ngetObcBornRadii2\n" );
      _brookGbsa->getObcBornRadii2()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nForceStreams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

   }

// ---------------------------------------------------------------------------------------

   // gather for second loop

    float mergeNonObcForces = 1.0f;
    float kcalMolTokJNM     = -0.4184f;
    kPostObcLoop2_nobranch(
              (float) _brookGbsa->getDuplicationFactor(),
              (float) _brookGbsa->getAtomStreamWidth( ),
              (float) _brookGbsa->getPartialForceStreamWidth( ),
              (float) _brookGbsa->getNumberOfAtoms(),
              (float) _brookGbsa->getAtomSizeCeiling(),
              (float) _brookGbsa->getInnerLoopUnroll(),
              kcalMolTokJNM,
              mergeNonObcForces,
              _brookGbsa->getObcIntermediateForce()->getBrookStream(),
              forceStream.getBrookStream(),
              gbsaForceStreams[0]->getBrookStream(),
              gbsaForceStreams[1]->getBrookStream(),
              gbsaForceStreams[2]->getBrookStream(),
              gbsaForceStreams[3]->getBrookStream(),
              _brookGbsa->getObcAtomicRadii()->getBrookStream(),
              _brookGbsa->getObcBornRadii()->getBrookStream(),
              _brookGbsa->getObcChain()->getBrookStream(),
              forceStream.getBrookStream()
           );

// ---------------------------------------------------------------------------------------

   // diagnostics

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "\nPost kPostObcLoop2_nobranch: atms=%d ceil=%d dup=%d atomStrW=%3d prtlF=%3d diel=%.3f %.3f ACE=%.1f\n",
                      _brookGbsa->getNumberOfAtoms(),
                      _brookGbsa->getAtomSizeCeiling(),
                      _brookGbsa->getDuplicationFactor(),
                      _brookGbsa->getAtomStreamWidth( ),
                      _brookGbsa->getPartialForceStreamWidth( ),
                      _brookGbsa->getSoluteDielectric(),
                      _brookGbsa->getSolventDielectric(), includeAce );

      (void) fprintf( getLog(), "\nPartialForceStreams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         gbsaForceStreams[ii]->printToFile( getLog() );
      }

      BrookStreamInternal* brookStreamInternalF  = forceStream.getBrookStreamImpl();
      (void) fprintf( getLog(), "\nForceStream\n" );
      brookStreamInternalF->printToFile( getLog() );

      (void) fprintf( getLog(), "\nChain\n" );
      _brookGbsa->getObcChain()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nBornR\n" );
      _brookGbsa->getObcBornRadii()->printToFile( getLog() );

   }

   // ---------------------------------------------------------------------------------------
}

/**
 * Execute the kernel to calculate the OBC energy
 * 
 * @param positions   atom positions
 *
 * @return  potential energy due to the OBC forces
 *
 */

double BrookCalcGBSAOBCForceFieldKernel::executeEnergy( const Stream& positions ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcGBSAOBCForceFieldKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   return (double) _brookGbsa->getEnergy( positions );

}
