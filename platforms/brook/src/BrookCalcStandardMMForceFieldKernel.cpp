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
#include "BrookCalcStandardMMForceFieldKernel.h"
#include "gpu/kforce.h"
#include "gpu/kinvmap_gather.h"
#include "StandardMMForceField.h"

using namespace OpenMM;
using namespace std;

/** 
 * BrookCalcStandardMMForceFieldKernel constructor
 * 
 * @param name                      kernel name
 * @param platform                  platform
 *
 */

BrookCalcStandardMMForceFieldKernel::BrookCalcStandardMMForceFieldKernel( std::string name, const Platform& platform ) :
                     CalcStandardMMForceFieldKernel( name, platform ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::BrookCalcStandardMMForceFieldKernel";

// ---------------------------------------------------------------------------------------

   _numberOfAtoms                           = 0;
   _brookBonded                             = NULL;
   _brookNonBonded                          = NULL;
   _refForceField                           = NULL;
   _log                                     = NULL;

   _refForceField                           = NULL;
   _refSystem                               = NULL;
   _refOpenMMContext                        = NULL;
   _referencePlatform                       = NULL;
   _refVerletIntegrator                     = NULL;

   const BrookPlatform brookPlatform = dynamic_cast<const BrookPlatform&> (platform);
   if( brookPlatform.getLog() != NULL ){
      setLog( brookPlatform.getLog() );
   }
      
}   

/** 
 * BrookCalcStandardMMForceFieldKernel destructor
 * 
 */

BrookCalcStandardMMForceFieldKernel::~BrookCalcStandardMMForceFieldKernel( ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::BrookCalcStandardMMForceFieldKernel";

// ---------------------------------------------------------------------------------------

   delete _brookBonded;
   delete _brookNonBonded;

   // deleted w/ kernel delete? If activated, program crashes

   //delete _refForceField;
/*
   delete _refSystem;
   delete _refOpenMMContext;
   delete _referencePlatform;
   delete _refVerletIntegrator;
*/
}

/** 
 * Get log file reference
 * 
 * @return  log file reference
 *
 */

FILE* BrookCalcStandardMMForceFieldKernel::getLog( void ) const {
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

int BrookCalcStandardMMForceFieldKernel::setLog( FILE* log ){
   _log = log;
   return BrookCommon::DefaultReturnValue;
}

/** 
 * Initialize the kernel, setting up the values of all the force field parameters.
 * 
 * @param bondIndices               the two atoms connected by each bond term
 * @param bondParameters            the force parameters (length, k) for each bond term
 * @param angleIndices              the three atoms connected by each angle term
 * @param angleParameters           the force parameters (angle, k) for each angle term
 * @param periodicTorsionIndices    the four atoms connected by each periodic torsion term
 * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
 * @param rbTorsionIndices          the four atoms connected by each Ryckaert-Bellemans torsion term
 * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
 * @param bonded14Indices           each element contains the indices of two atoms whose nonbonded interactions should be reduced since
 *                                  they form a bonded 1-4 pair
 * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
 * @param coulomb14Scale            the factor by which Coulomb interactions should be reduced for bonded 1-4 pairs
 * @param exclusions                the i'th element lists the indices of all atoms with which the i'th atom should not interact through
 *                                  nonbonded forces.  Bonded 1-4 pairs are also included in this list, since they should be omitted from
 *                                  the standard nonbonded calculation.
 * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each atom
 * @param nonbondedMethod           the method to use for handling long range nonbonded interactions
 * @param nonbondedCutoff           the cutoff distance for nonbonded interactions (if nonbondedMethod involves a cutoff)
 * @param periodicBoxSize           the size of the periodic box (if nonbondedMethod involves a periodic boundary conditions)
 */
 
void BrookCalcStandardMMForceFieldKernel::initialize( 
        const vector<vector<int> >& bondIndices,            const vector<vector<double> >& bondParameters,
        const vector<vector<int> >& angleIndices,           const vector<vector<double> >& angleParameters,
        const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
        const vector<vector<int> >& rbTorsionIndices,       const vector<vector<double> >& rbTorsionParameters,
        const vector<vector<int> >& bonded14Indices,        double lj14Scale, double coulomb14Scale,
        const vector<set<int> >& exclusions, const vector<vector<double> >& nonbondedParameters,
         NonbondedMethod nonbondedMethod, double nonbondedCutoff, double periodicBoxSize[3] ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::initialize";

// ---------------------------------------------------------------------------------------

   FILE* log                 = getLog();
   _numberOfAtoms            = nonbondedParameters.size();

   // ---------------------------------------------------------------------------------------

   // bonded

   if( _brookBonded ){
      delete _brookBonded;
   }
   _brookBonded              = new BrookBonded();
   _brookBonded->setLog( log );
    
   _brookBonded->setup( _numberOfAtoms, 
                        bondIndices,            bondParameters, 
                        angleIndices,           angleParameters,
                        periodicTorsionIndices, periodicTorsionParameters,
                        rbTorsionIndices,       rbTorsionParameters,
                        bonded14Indices,        nonbondedParameters,
                        lj14Scale, coulomb14Scale, getPlatform() );

   // echo contents

   if( log ){
      std::string contents = _brookBonded->getContentsString( ); 
      (void) fprintf( log, "%s brookBonded::contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   // nonbonded

   if( _brookNonBonded ){
      delete _brookNonBonded;
   }
   _brookNonBonded           = new BrookNonBonded();
   _brookNonBonded->setLog( log );

   _brookNonBonded->setup( _numberOfAtoms, nonbondedParameters, exclusions, getPlatform() );

   // echo contents

   if( log ){
      std::string contents = _brookNonBonded->getContentsString( ); 
      (void) fprintf( log, "%s brookNonBonded::contents\n%s", methodName.c_str(), contents.c_str() );
      (void) fflush( log );
   }

   // ---------------------------------------------------------------------------------------

   // used for calculating energy

/*
   _referenceCalcStandardMMForceFieldKernel = new ReferenceCalcStandardMMForceFieldKernel( getName(), getPlatform() );
   _referenceCalcStandardMMForceFieldKernel->initialize( bondIndices,            bondParameters, 
                                                         angleIndices,           angleParameters,
                                                         periodicTorsionIndices, periodicTorsionParameters,
                                                         rbTorsionIndices,       rbTorsionParameters,
                                                         bonded14Indices,        lj14Scale, coulomb14Scale,
                                                         exclusions, nonbondedParameters,
                                                         nonbondedMethod, nonbondedCutoff, periodicBoxSize );
*/
   _refForceField = new StandardMMForceField( _numberOfAtoms, bondIndices.size(), angleIndices.size(), 
                                              periodicTorsionIndices.size(),  rbTorsionIndices.size() );

   // bonds

   vector<vector<int> >::const_iterator bondI_it    = bondIndices.begin();
   vector<vector<double> >::const_iterator bondP_it = bondParameters.begin();
   for( unsigned int ii = 0; ii < bondIndices.size(); ii++, bondI_it++, bondP_it++ ){
      const vector<int>    bondInd = *bondI_it;
      const vector<double> bondPrm = *bondP_it;
      _refForceField->setBondParameters( ii, bondInd[0], bondInd[1], bondPrm[0], bondPrm[1] );
   }

   // angles

   vector<vector<int> >::const_iterator angleI_it    = angleIndices.begin();
   vector<vector<double> >::const_iterator angleP_it = angleParameters.begin();
   for( unsigned int ii = 0; ii < angleIndices.size(); ii++, angleI_it++, angleP_it++ ){
      const vector<int>    angleInd = *angleI_it;
      const vector<double> anglePrm = *angleP_it;
      _refForceField->setAngleParameters( ii, angleInd[0], angleInd[1], angleInd[2], anglePrm[0], anglePrm[1] );
   }

   // periodicTorsions

   vector<vector<int> >::const_iterator periodicTorsionI_it    = periodicTorsionIndices.begin();
   vector<vector<double> >::const_iterator periodicTorsionP_it = periodicTorsionParameters.begin();
   for( unsigned int ii = 0; ii < periodicTorsionIndices.size(); ii++, periodicTorsionI_it++, periodicTorsionP_it++ ){
      const vector<int>    periodicTorsionInd = *periodicTorsionI_it;
      const vector<double> periodicTorsionPrm = *periodicTorsionP_it;
      _refForceField->setPeriodicTorsionParameters( ii, periodicTorsionInd[0], periodicTorsionInd[1], periodicTorsionInd[2], periodicTorsionInd[3],
                                                    (int) (periodicTorsionPrm[2] + 0.001), periodicTorsionPrm[1], periodicTorsionPrm[0] );
/*
printf( "PeriodicTor: [%d %d %d %d] [%.5e %.5e %.5e]\n", periodicTorsionInd[0], periodicTorsionInd[1], periodicTorsionInd[2], periodicTorsionInd[3],
                                                         periodicTorsionPrm[2], periodicTorsionPrm[1], periodicTorsionPrm[0] ); fflush( stdout );
*/
   }

   // rbTorsions

   vector<vector<int> >::const_iterator rbTorsionI_it    = rbTorsionIndices.begin();
   vector<vector<double> >::const_iterator rbTorsionP_it = rbTorsionParameters.begin();
   for( unsigned int ii = 0; ii < rbTorsionIndices.size(); ii++, rbTorsionI_it++, rbTorsionP_it++ ){
      const vector<int>    rbTorsionInd = *rbTorsionI_it;
      const vector<double> rbTorsionPrm = *rbTorsionP_it;
      _refForceField->setRBTorsionParameters( ii, rbTorsionInd[0], rbTorsionInd[1], rbTorsionInd[2], rbTorsionInd[3],
                                              rbTorsionPrm[0], rbTorsionPrm[1], rbTorsionPrm[2], rbTorsionPrm[3], 
                                              rbTorsionPrm[4], rbTorsionPrm[5] );

/*
printf( "RbTor: [%d %d %d %d] [%.5e %.5e %.5e %.5e %.5e %.5e]\n", rbTorsionInd[0], rbTorsionInd[1], rbTorsionInd[2], rbTorsionInd[3],
                                                         rbTorsionPrm[0], rbTorsionPrm[1], rbTorsionPrm[2], rbTorsionPrm[3], rbTorsionPrm[4], rbTorsionPrm[5] ); fflush( stdout );
*/
   }

   // nonbonded

   for( unsigned int ii = 0; ii < nonbondedParameters.size(); ii++ ){

      vector<double> nonbondedParameterVector 
                                     = nonbondedParameters[ii];

      double c6                      = nonbondedParameterVector[1];
      double c12                     = nonbondedParameterVector[2];
      double charge                  = nonbondedParameterVector[0];

      // int index, double charge, double radius, double depth

      _refForceField->setAtomParameters( ii, charge, c6, c12 );

   }

   // ---------------------------------------------------------------------------------------

   return;
}

/** 
 * Execute the kernel to calculate the bonded & nonbonded forces
 * 
 * @param positions   stream of type Double3 containing the position (x, y, z) of each atom
 * @param forces      stream of type Double3 containing the force (x, y, z) on each atom.  On entry, this contains the forces that
 *                    have been calculated so far.  The kernel should add its own forces to the values already in the stream.
 */

void BrookCalcStandardMMForceFieldKernel::executeForces( const Stream& positions, Stream& forces ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::executeForces";

   static const int I_Stream                = 0;
   static const int J_Stream                = 1;
   static const int K_Stream                = 2;
   static const int L_Stream                = 3;

   static const int PrintOn                 = 0;
   static const int MaxErrorMessages        = 2;
   static       int ErrorMessages           = 0;

   static const float4 dummyParameters( 0.0, 0.0, 0.0, 0.0 );

   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------


   const BrookStreamImpl& positionStreamC              = dynamic_cast<const BrookStreamImpl&> (positions.getImpl());
   BrookStreamImpl& positionStream                     = const_cast<BrookStreamImpl&>         (positionStreamC);
   BrookStreamImpl& forceStream                        = dynamic_cast<BrookStreamImpl&>       (forces.getImpl());
   
   // nonbonded forces

   // added charge stream to knbforce_CDLJ4
   // libs generated from ~/src/gmxgpu-nsqOpenMM

   BrookFloatStreamInternal**  nonbondedForceStreams = _brookNonBonded->getForceStreams();

   float epsfac                                      = 138.935485f;

   knbforce_CDLJ4(
             (float) _brookNonBonded->getNumberOfAtoms(),
             (float) _brookNonBonded->getAtomSizeCeiling(),
             (float) _brookNonBonded->getDuplicationFactor(),
             (float) _brookNonBonded->getAtomStreamHeight( ),
             (float) _brookNonBonded->getAtomStreamWidth( ),
             (float) _brookNonBonded->getJStreamWidth( ),
             (float) _brookNonBonded->getPartialForceStreamWidth( ),
             epsfac,
             dummyParameters,
             positionStream.getBrookStream(),
             _brookNonBonded->getChargeStream()->getBrookStream(),
             _brookNonBonded->getOuterVdwStream()->getBrookStream(),
             _brookNonBonded->getInnerSigmaStream()->getBrookStream(),
             _brookNonBonded->getInnerEpsilonStream()->getBrookStream(),
             _brookNonBonded->getExclusionStream()->getBrookStream(),
             nonbondedForceStreams[0]->getBrookStream(),
             nonbondedForceStreams[1]->getBrookStream(),
             nonbondedForceStreams[2]->getBrookStream(),
             nonbondedForceStreams[3]->getBrookStream()
           );

/*
float zerof = 0.0f;
nonbondedForceStreams[0]->fillWithValue( &zerof );
nonbondedForceStreams[1]->fillWithValue( &zerof );
nonbondedForceStreams[2]->fillWithValue( &zerof );
nonbondedForceStreams[3]->fillWithValue( &zerof );
*/

   // diagnostics

   if( 1 && PrintOn ){
      (void) fprintf( getLog(), "\nPost knbforce_CDLJ4: atoms=%6d ceiling=%3d dupFac=%3d", _brookNonBonded->getNumberOfAtoms(),  
                                                                                           _brookNonBonded->getAtomSizeCeiling(),
                                                                                           _brookNonBonded->getDuplicationFactor()  );

      (void) fprintf( getLog(), "\n                      hght=%6d   width=%3d   jWid=%3d", _brookNonBonded->getAtomStreamHeight( ),
                                                                                           _brookNonBonded->getAtomStreamWidth( ),
                                                                                           _brookNonBonded->getJStreamWidth( ) );
      (void) fprintf( getLog(), "\n                      pFrc=%6d     eps=%12.5e\n",       _brookNonBonded->getPartialForceStreamWidth( ), epsfac );

      (void) fprintf( getLog(), "\nOuterVdwStreamd\n" );
      _brookNonBonded->getOuterVdwStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nInnerSigmaStream\n" );
      _brookNonBonded->getInnerSigmaStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nInnerEpsilonStream\n" );
      _brookNonBonded->getInnerEpsilonStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nExclusionStream\n" );
      _brookNonBonded->getExclusionStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nChargeStream\n" );
      _brookNonBonded->getChargeStream()->printToFile( getLog() );

      for( int ii = 0; ii < 4; ii++ ){
         (void) fprintf( getLog(), "\nForce stream %d\n", ii );
         nonbondedForceStreams[ii]->printToFile( getLog() );
      }
   }

// ---------------------------------------------------------------------------------------

   // gather forces

   kMergeFloat3_4_nobranch( (float) _brookNonBonded->getDuplicationFactor(),
                            (float) _brookNonBonded->getAtomStreamWidth(),
                            (float) _brookNonBonded->getPartialForceStreamWidth(),
                            (float) _brookNonBonded->getNumberOfAtoms(),
                            (float) _brookNonBonded->getAtomSizeCeiling(),
                            (float) _brookNonBonded->getOuterLoopUnroll(),
                            nonbondedForceStreams[0]->getBrookStream(),
                            nonbondedForceStreams[1]->getBrookStream(),
                            nonbondedForceStreams[2]->getBrookStream(),
                            nonbondedForceStreams[3]->getBrookStream(),
                            forceStream.getBrookStream() );

   // diagnostics

   if( 0 && PrintOn ){

      (void) fprintf( getLog(), "\nNB forces" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );

   }

// ---------------------------------------------------------------------------------------

   // bonded

         epsfac                                        = (float) (_brookBonded->getLJ_14Scale()*_brookBonded->getCoulombFactor());
   float width                                         = (float) (_brookBonded->getInverseMapStreamWidth());

   // bonded forces

   BrookFloatStreamInternal**  bondedParameters        = _brookBonded->getBondedParameterStreams();
   BrookFloatStreamInternal**  bondedForceStreams      = _brookBonded->getBondedForceStreams();

   BrookFloatStreamInternal**  inverseStreamMaps[4];
   inverseStreamMaps[0]                                = _brookBonded->getInverseStreamMapsStreams( 0 );
   inverseStreamMaps[1]                                = _brookBonded->getInverseStreamMapsStreams( 1 );
   inverseStreamMaps[2]                                = _brookBonded->getInverseStreamMapsStreams( 2 );
   inverseStreamMaps[3]                                = _brookBonded->getInverseStreamMapsStreams( 3 );

   kbonded_CDLJ( epsfac, 
                 (float) bondedForceStreams[0]->getStreamWidth(),
                 dummyParameters,
                 positionStream.getBrookStream(),
                 _brookBonded->getChargeStream()->getBrookStream(),
                 _brookBonded->getAtomIndicesStream()->getBrookStream(),
                 bondedParameters[0]->getBrookStream(),
                 bondedParameters[1]->getBrookStream(),
                 bondedParameters[2]->getBrookStream(),
                 bondedParameters[3]->getBrookStream(),
                 bondedParameters[4]->getBrookStream(),
                 bondedForceStreams[0]->getBrookStream(),
                 bondedForceStreams[1]->getBrookStream(),
                 bondedForceStreams[2]->getBrookStream(),
                 bondedForceStreams[3]->getBrookStream() );


   // diagnostics

   if( 1 && PrintOn ){

      int countPrintInvMap[4] = { 3, 5, 2, 4 }; 

      (void) fprintf( getLog(), "\nPost kbonded_CDLJ: epsFac=%.6f %.6f %.6f", epsfac, _brookBonded->getLJ_14Scale(), _brookBonded->getCoulombFactor());
      (void) fprintf( getLog(), "\nAtom indices stream\n" );
      _brookBonded->getAtomIndicesStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nCharge stream\n" );
      _brookBonded->getChargeStream()->printToFile( getLog() );

      for( int ii = 0; ii < 5; ii++ ){
         (void) fprintf( getLog(), "\nParam stream %d\n", ii );
         bondedParameters[ii]->printToFile( getLog() );
      }
      for( int ii = 0; ii < 4; ii++ ){
         (void) fprintf( getLog(), "\nForce stream %d\n", ii );
         bondedForceStreams[ii]->printToFile( getLog() );
      }

/*
      (void) fprintf( getLog(), "\nNB1 forces" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );
*/

      (void) fprintf( getLog(), "\nInverse map streams -- K_Stream cnt=%d\n", _brookBonded->getInverseMapStreamCount( K_Stream ) );
      for( int ii = 0; ii < 4; ii++ ){
         for( int jj = 0; jj < countPrintInvMap[ii]; jj++ ){
            (void) fprintf( getLog(), "\n   Inverse map streams index=%d %d\n", ii, jj );
            inverseStreamMaps[ii][jj]->printToFile( getLog() );
         }
      }
   }

   // gather forces

   if( _brookBonded->getInverseMapStreamCount( I_Stream ) == 3 && _brookBonded->getInverseMapStreamCount( K_Stream ) == 3 ){

      kinvmap_gather3_3( width,

                         inverseStreamMaps[I_Stream][0]->getBrookStream(),
                         inverseStreamMaps[I_Stream][1]->getBrookStream(),
                         inverseStreamMaps[I_Stream][2]->getBrookStream(),
                         bondedForceStreams[I_Stream]->getBrookStream(),

                         inverseStreamMaps[K_Stream][0]->getBrookStream(),
                         inverseStreamMaps[K_Stream][1]->getBrookStream(),
                         inverseStreamMaps[K_Stream][2]->getBrookStream(),
                         bondedForceStreams[K_Stream]->getBrookStream(),

                         forceStream.getBrookStream(), forceStream.getBrookStream() );

   } else if( _brookBonded->getInverseMapStreamCount( I_Stream ) == 3 && _brookBonded->getInverseMapStreamCount( K_Stream ) == 4 ){

      kinvmap_gather3_4( width,

                         inverseStreamMaps[I_Stream][0]->getBrookStream(),
                         inverseStreamMaps[I_Stream][1]->getBrookStream(),
                         inverseStreamMaps[I_Stream][2]->getBrookStream(),
                         bondedForceStreams[I_Stream]->getBrookStream(),

                         inverseStreamMaps[K_Stream][0]->getBrookStream(),
                         inverseStreamMaps[K_Stream][1]->getBrookStream(),
                         inverseStreamMaps[K_Stream][2]->getBrookStream(),
                         inverseStreamMaps[K_Stream][3]->getBrookStream(),
                         bondedForceStreams[K_Stream]->getBrookStream(),

                         forceStream.getBrookStream(), forceStream.getBrookStream() );

   } else if( _brookBonded->getInverseMapStreamCount( I_Stream ) == 3 && _brookBonded->getInverseMapStreamCount( K_Stream ) == 5 ){

      kinvmap_gather3_5( width,
                         inverseStreamMaps[I_Stream][0]->getBrookStream(),
                         inverseStreamMaps[I_Stream][1]->getBrookStream(),
                         inverseStreamMaps[I_Stream][2]->getBrookStream(),
                         bondedForceStreams[I_Stream]->getBrookStream(),
                         inverseStreamMaps[K_Stream][0]->getBrookStream(),
                         inverseStreamMaps[K_Stream][1]->getBrookStream(),
                         inverseStreamMaps[K_Stream][2]->getBrookStream(),
                         inverseStreamMaps[K_Stream][3]->getBrookStream(),
                         inverseStreamMaps[K_Stream][4]->getBrookStream(),
                         bondedForceStreams[K_Stream]->getBrookStream(),
                         forceStream.getBrookStream(), forceStream.getBrookStream() );

   } else {

      // case not handled -- throw an exception

      if( _brookBonded->getLog() && ErrorMessages++ < MaxErrorMessages ){
         (void) fprintf( _brookBonded->getLog(), "%s case: I-map=%d K-map=%d -- not handled.\n",
                          methodName.c_str(), _brookBonded->getInverseMapStreamCount( I_Stream ),
                                              _brookBonded->getInverseMapStreamCount( K_Stream ) );
         (void) fflush(  _brookBonded->getLog() );
      }

      kinvmap_gather3_3( width,

                         inverseStreamMaps[I_Stream][0]->getBrookStream(),
                         inverseStreamMaps[I_Stream][1]->getBrookStream(),
                         inverseStreamMaps[I_Stream][2]->getBrookStream(),
                         bondedForceStreams[I_Stream]->getBrookStream(),

                         inverseStreamMaps[K_Stream][0]->getBrookStream(),
                         inverseStreamMaps[K_Stream][1]->getBrookStream(),
                         inverseStreamMaps[K_Stream][2]->getBrookStream(),
                         bondedForceStreams[K_Stream]->getBrookStream(),

                         forceStream.getBrookStream(), forceStream.getBrookStream() );
/*
      std::stringstream message;
      message << methodName << "I-maps=" << _brookBonded->getInverseMapStreamCount( I_Stream ) << " and " << 
                               "K-maps=" << _brookBonded->getInverseMapStreamCount( K_Stream ) << " not handled.";
      throw OpenMMException( message.str() );
*/

   }

   // diagnostics

   if( 0 && PrintOn ){

      (void) fprintf( getLog(), "\nPost 3_4/3_5 && NB forces" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );

   }

   if( _brookBonded->getInverseMapStreamCount( J_Stream ) == 5 && _brookBonded->getInverseMapStreamCount( L_Stream ) == 2 ){
      kinvmap_gather5_2( width,
                         inverseStreamMaps[J_Stream][0]->getBrookStream(),
                         inverseStreamMaps[J_Stream][1]->getBrookStream(),
                         inverseStreamMaps[J_Stream][2]->getBrookStream(),
                         inverseStreamMaps[J_Stream][3]->getBrookStream(),
                         inverseStreamMaps[J_Stream][4]->getBrookStream(),
                         bondedForceStreams[J_Stream]->getBrookStream(),
                         inverseStreamMaps[L_Stream][0]->getBrookStream(),
                         inverseStreamMaps[L_Stream][1]->getBrookStream(),
                         bondedForceStreams[L_Stream]->getBrookStream(),
                         forceStream.getBrookStream(), forceStream.getBrookStream() );
   
   } else {

      // case not handled -- throw an exception

      if( _brookBonded->getLog() && ErrorMessages++ < MaxErrorMessages ){
         (void) fprintf( _brookBonded->getLog(), "%s case: J-map=%d L-map=%d -- not handled.\n",
                          methodName.c_str(), _brookBonded->getInverseMapStreamCount( J_Stream ),
                                              _brookBonded->getInverseMapStreamCount( L_Stream ) );
         (void) fflush(  _brookBonded->getLog() );
      }

      kinvmap_gather5_2( width,
                         inverseStreamMaps[J_Stream][0]->getBrookStream(),
                         inverseStreamMaps[J_Stream][1]->getBrookStream(),
                         inverseStreamMaps[J_Stream][2]->getBrookStream(),
                         inverseStreamMaps[J_Stream][3]->getBrookStream(),
                         inverseStreamMaps[J_Stream][4]->getBrookStream(),
                         bondedForceStreams[J_Stream]->getBrookStream(),
                         inverseStreamMaps[L_Stream][0]->getBrookStream(),
                         inverseStreamMaps[L_Stream][1]->getBrookStream(),
                         bondedForceStreams[L_Stream]->getBrookStream(),
                         forceStream.getBrookStream(), forceStream.getBrookStream() );
/*
      std::stringstream message;
      message << methodName << "J-maps=" << _brookBonded->getInverseMapStreamCount( J_Stream ) << " and " << 
                               "L-maps=" << _brookBonded->getInverseMapStreamCount( L_Stream ) << " not handled.";
      throw OpenMMException( message.str() );
*/

   }

   // diagnostics

   if( 1 && PrintOn ){

      (void) fprintf( getLog(), "\nFinal NB & bonded forces" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );
/*
      void* dataV = brookStreamInternalF->getData(1);
      float* data = (float*) dataV;
      (void) fprintf( getLog(), "\nFinal NB & bonded forces RAW\n" );
      for( int ii = 0; ii < _brookNonBonded->getNumberOfAtoms()*3; ii += 3 ){
         (void) fprintf( getLog(), "%d [%.6e %.6e %.6e]\n", ii, data[ii], data[ii+1], data[ii+2] );
      }
*/

   }

   // ---------------------------------------------------------------------------------------
}

/**
 * Execute the kernel to calculate the energy.
 * 
 * @param positions   atom positions
 *
 * @return  potential energy due to the StandardMMForceField
 * Currently always return 0.0 since energies not calculated on gpu
 *
 */

double BrookCalcStandardMMForceFieldKernel::executeEnergy( const Stream& atomPositions ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   const BrookStreamImpl& positionStreamC              = dynamic_cast<const BrookStreamImpl&> (atomPositions.getImpl());
   BrookStreamImpl& positionStream                     = const_cast<BrookStreamImpl&>         (positionStreamC);
   BrookOpenMMFloat* positionsF                        = (BrookOpenMMFloat*) positionStream.getData();

   OpenMMContext* context                              = getReferenceOpenMMContext( atomPositions.getSize() );

   vector<Vec3> positions( positionStream.getSize() );
   int index = 0;
   for( int ii = 0; ii < positionStream.getSize(); ii++, index += 3 ){
      positions[ii] = Vec3( positionsF[index], positionsF[index+1], positionsF[index+2] );
   }
   context->setPositions(positions);

   State state     = context->getState( State::Energy );
   double energy   = state.getPotentialEnergy();

// (void) fprintf( stdout, "BrookCalcStandardMMForceFieldKernel::executeEnergy E=%.5e\n", energy ); fflush( stdout );

   return energy;

}

/**
 * Get reference Context
 * 
 * @param numberOfAtoms  number of atoms
 *
 * @return  OpenMMContext
 *
 */

OpenMMContext* BrookCalcStandardMMForceFieldKernel::getReferenceOpenMMContext( int numberOfAtoms ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::getReferenceOpenMMContext";

// ---------------------------------------------------------------------------------------

   if( _refOpenMMContext == NULL ){
      
      _referencePlatform       = new ReferencePlatform();
      _refSystem               = new System( numberOfAtoms, 0 ); 
      _refVerletIntegrator     = new VerletIntegrator( 0.01 );

      _refSystem->addForce( _refForceField );

      _refOpenMMContext        = new OpenMMContext( *_refSystem, *_refVerletIntegrator, *_referencePlatform );
   }

   return _refOpenMMContext;

}

/**
 * Execute the kernel to calculate the energy.
 * 
 * @param positions   atom positions
 *
 * @return  potential energy due to the StandardMMForceField
 * Currently always return 0.0 since energies not calculated on gpu
 *
 */

double BrookCalcStandardMMForceFieldKernel::executeEnergyOld( const Stream& atomPositions ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName      = "BrookCalcStandardMMForceFieldKernel::executeEnergy";

// ---------------------------------------------------------------------------------------

   const BrookStreamImpl& positionStreamC              = dynamic_cast<const BrookStreamImpl&> (atomPositions.getImpl());
   BrookStreamImpl& positionStream                     = const_cast<BrookStreamImpl&>         (positionStreamC);
   BrookOpenMMFloat* positionsF                        = (BrookOpenMMFloat*) positionStream.getData();

   ReferencePlatform refPlatform;
   System system( atomPositions.getSize(), 0); 
   VerletIntegrator integrator( 0.01 );

   system.addForce( _refForceField );

   OpenMMContext context( system, integrator, refPlatform );

   vector<Vec3> positions( positionStream.getSize() );
   int index = 0;
   for( int ii = 0; ii < positionStream.getSize(); ii++, index += 3 ){
      positions[ii] = Vec3( positionsF[index], positionsF[index+1], positionsF[index+2] );
   }
   context.setPositions(positions);

   State state     = context.getState( State::Energy );
   double energy   = state.getPotentialEnergy();

   //delete forceField;

// (void) fprintf( stdout, "BrookCalcStandardMMForceFieldKernel::executeEnergy E=%.5e\n", energy ); fflush( stdout );

   return energy;

}
