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

#include <math.h>
#include <sstream>
#include "BrookBonded.h"
#include "BrookPlatform.h"
#include "BrookStreamFactory.h"
#include "OpenMMException.h"
#include "gpu/kinvmap_gather.h"
#include "gpu/invmap.h"
#include "gpu/kforce.h"

using namespace OpenMM;
using namespace std;

#define ATOMS(X,Y) (particles[ 5*(X) + (Y) + 1 ])
#define PARAMS(X,Y,Z) (params[(Y)][4*(X) + Z])

/** 
 *
 * BrookBonded constructor
 * 
 */

BrookBonded::BrookBonded( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::BrookBonded";

// ---------------------------------------------------------------------------------------

   _setupCompleted            = 0;
   _numberOfParticles         = 0;

   _coulombFactor             = (BrookOpenMMFloat) 138.935485;

   _particleIndicesStream     = NULL;
   _chargeStream              = NULL;

   // parameter streams

   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      _bondedParameters[ii] = NULL;
   }

   // inverse maps & force streams

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      _bondedForceStreams[ii]       = NULL;
      _inverseMapStreamCount[ii]    = 0;
      for( int jj = 0; jj < MaxNumberOfInverseMaps; jj++ ){
         _inverseStreamMaps[ii][jj] = NULL;
      }
   }

   _maxInverseMapStreamCount[StreamI] = 9;
   _maxInverseMapStreamCount[StreamJ] = 6;
   _maxInverseMapStreamCount[StreamK] = 6;
   _maxInverseMapStreamCount[StreamL] = 9;

   _maxNumberOfInverseMaps = _maxInverseMapStreamCount[0];
   for( int ii = 1; ii < getNumberOfForceStreams(); ii++ ){
      if( _maxNumberOfInverseMaps < _maxInverseMapStreamCount[ii] ){
         _maxNumberOfInverseMaps  = _maxInverseMapStreamCount[ii];
      }
   }

   // check that MaxNumberOfInverseMaps is big enough

   if( _maxNumberOfInverseMaps > MaxNumberOfInverseMaps ){
      std::stringstream message;
      message << methodName << " max number of inverse maps=" << _maxNumberOfInverseMaps << " is greater than hardwired value=" << MaxNumberOfInverseMaps;
      throw OpenMMException( message.str() );
   }

   _inverseMapStreamWidth = -1;
}   
 
/** 
 * BrookBonded destructor
 * 
 */

BrookBonded::~BrookBonded( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName = "BrookBonded::~BrookBonded";

// ---------------------------------------------------------------------------------------

   delete _particleIndicesStream;
   delete _chargeStream;

   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      delete _bondedParameters[ii];
   }
   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      delete _bondedForceStreams[ii];
      for( int jj = 0; jj < MaxNumberOfInverseMaps; jj++ ){
         delete _inverseStreamMaps[ii][jj];
      }
   }

}

/** 
 * Get number of parameter streams
 * 
 * @return  number of parameter streams
 *
 */

int BrookBonded::getNumberOfParameterStreams( void ) const {
   return NumberOfParameterStreams;
}

/** 
 * Get number of force streams
 * 
 * @return  number of force streams
 *
 */

int BrookBonded::getNumberOfForceStreams( void ) const {
   return NumberOfForceStreams;
}

/** 
 * Get max number of inverse maps
 * 
 * @return  max number of inverse maps
 *
 */

int BrookBonded::getMaxInverseMapStreamCount( void ) const {
   return _maxNumberOfInverseMaps;
}

/** 
 * Get max number of inverse maps for specified force stream index
 * 
 * @param index index of force stream
 *
 * @return  max number of inverse maps or -1 if index is out of range
 *
 */

int BrookBonded::getMaxInverseMapStreamCount( int index ) const {
   return (index >= 0 && index < getNumberOfForceStreams())  ? _maxInverseMapStreamCount[index] : -1;
}

/** 
 * Get width of inverse map streams
 * 
 * @return  width of inverse map streams
 *
 */

int BrookBonded::getInverseMapStreamWidth( void ) const {
   return _inverseMapStreamWidth;
}

/** 
 * Get Coulomb factor
 * 
 * @return Coulomb factor
 *
 */

BrookOpenMMFloat BrookBonded::getCoulombFactor( void ) const {
   return _coulombFactor;
}

/** 
 * Return true if force[index] stream is set
 *
 * @param    index into force stream
 * @return   true if index is valid && force[index] stream is set; else false
 *
 */

int BrookBonded::isForceStreamSet( int index ) const {
   return (index >= 0 && index < getNumberOfForceStreams() && _bondedForceStreams[index]) ? 1 : 0;
}

/** 
 * Return SetupCompleted flag
 *
 * @return SetupCompleted flag
 */
          
int BrookBonded::isSetupCompleted( void ) const {
   return _setupCompleted;
}
   
/** 
 * Set SetupCompleted flag
 *
 * @param  setupCompleted flag
 *
 * @return SetupCompleted flag
 */
          
int BrookBonded::setupCompleted( int setupCompleted ){
   _setupCompleted = setupCompleted;
   return _setupCompleted;
}
   
/** 
 * Return string showing if all inverse map streams are set
 *
 * @param    index into inverse map stream array
 *
 * @return   informative string
 *
 */

std::string BrookBonded::checkInverseMapStream( int index ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBonded::getContents";
   int ok                                   = 1;

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   if( index < 0 || index >= getNumberOfForceStreams() ){
      message << "Inverse map index=" << index << " is out of range [0, "  << getNumberOfForceStreams() << ")";
      return message.str();
   }

   message << "[ ";
   for( int ii = 0; ii < getInverseMapStreamCount( index ); ii++ ){
      message << (_inverseStreamMaps[index][ii] ? 1 : 0) << " ";
      if( !_inverseStreamMaps[index][ii] ){
         ok = 0;
      }
   }
   message << " ]";
   if( !ok ){
      message << " ERROR!";
   }

   return message.str();
   
}

/** 
 * Return true if paramsterSet[index] stream is set
 *
 * @param    index into parameter stream
 *
 * @return   true if index is valid && paramsterSet[index] stream is set; else false
 *
 */

int BrookBonded::isParameterStreamSet( int index ) const {
   return (index >= 0 && index < getNumberOfParameterStreams() && _bondedParameters[index]) ? 1 : 0;
}

/** 
 * Validate stream count
 * 
 * @return  DefaultReturnValue if count valid; else return ErrorReturnValue
 *
 * @exception OpenMMException is thrown if count is invalid
 *
 */

int BrookBonded::_validateInverseMapStreamCount( int index, int count ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_validateInverseMapStreamCount";

// ---------------------------------------------------------------------------------------

   int valid = DefaultReturnValue;

   if( index == 2 ){
      if( count > 5 ){
         valid = ErrorReturnValue;
      }
   }

   if( valid == ErrorReturnValue ){
      std::stringstream message;
      message << methodName << " input index=" << index << " has invalid count=" << count;
      throw OpenMMException( message.str() );
   }

   return DefaultReturnValue;
}

/** 
 * Get inverse map stream count
 * 
 * @return  count
 *
 * @exception OpenMMException is thrown if index is out of range [0, getNumberOfForceStreams() ]
 *
 */

int BrookBonded::getInverseMapStreamCount( int index ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::getInverseMapStreamCount";

// ---------------------------------------------------------------------------------------

   if( index < 0 || index >= getNumberOfForceStreams() ){
      std::stringstream message;
      message << methodName << " input index=" << index << " is out of range [0, " << getNumberOfForceStreams() << " )";
      throw OpenMMException( message.str() );
   }

   return _inverseMapStreamCount[index];
}

/** 
 * Get bonded particle indices stream
 * 
 * @return  particle indices stream
 *
 */

BrookFloatStreamInternal* BrookBonded::getParticleIndicesStream( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::getParticleIndicesStream";

// ---------------------------------------------------------------------------------------

   return _particleIndicesStream;
}

/** 
 * Get bonded charge stream
 * 
 * @return  charge stream
 *
 */

BrookFloatStreamInternal* BrookBonded::getChargeStream( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::getChargeStream";

// ---------------------------------------------------------------------------------------

   return _chargeStream;
}

/** 
 * Get array of bonded parameter streams
 * 
 * @return  array of bonded parameter streams
 *
 */

BrookFloatStreamInternal** BrookBonded::getBondedParameterStreams( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::getBondedParameterStreams";

// ---------------------------------------------------------------------------------------

   return _bondedParameters;
}

/** 
 * Get array of force streams
 * 
 * @return  array  force streams
 *
 */

BrookFloatStreamInternal** BrookBonded::getBondedForceStreams( void ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::getBondedForceStreams";

// ---------------------------------------------------------------------------------------

   return _bondedForceStreams;
}

/** 
 * Get array of inverse map streams
 * 
 * @param index  array index 
 *
 * @return  array inverse map streams
 *
 */

BrookFloatStreamInternal** BrookBonded::getInverseStreamMapsStreams( int index ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::getInverseStreamMapsStreams";

// ---------------------------------------------------------------------------------------

   // no checking on index -- assume ok -- speed

   return _inverseStreamMaps[index];
}

/* Gromacs sorts most bonded interactions
 * We sort those that gromacs forgets
 * The parameters are all symmetric with respect 
 * to inversion of order.
 * 
 * Also we assume that the same set of indices are not
 * used with different sets of parameters in gromacs
 * This is just an assumption for convenience here,
 * nothing stops us from evaluating such terms in the
 * kernel.
 * 
 * To begin with all interactions have -1 -1 -1 -1
 * After we fit in all the interactions, we have
 * to convert the -1's to zeros to avoid indexing
 * errors on the GPU.
 * */

/*Flips i,j,k,l to l,k,j,i while correctly shuffling the params */

void BrookBonded::_flipQuartet( int ibonded, int *particles ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::_flipQuartet";

// ---------------------------------------------------------------------------------------

   int tmp;

   //For now, simply flip the indices
   //we're just studying the packing

   tmp = ATOMS( ibonded, 0 );
   ATOMS( ibonded, 0 ) = ATOMS( ibonded, 3 );
   ATOMS( ibonded, 3 ) = ATOMS( ibonded, 0 );

   tmp = ATOMS( ibonded, 1 );
   ATOMS( ibonded, 1 ) = ATOMS( ibonded, 2 );
   ATOMS( ibonded, 2 ) = ATOMS( ibonded, 1 );

}

int BrookBonded::_matchTorsion( int i, int j, int k, int l, int nbondeds, int *particles ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::_matchTorsion";

// ---------------------------------------------------------------------------------------

   int tmp;
   if( j > k ){
      //swap into order
      tmp = i;
      i = l;
      l = tmp;

      tmp = j;
      j = k;
      k = tmp;
   }
   
   for( int n = 0; n < nbondeds; n++ ){
      if( i == ATOMS(n, 0) && j == ATOMS(n, 1) && 
          k == ATOMS(n, 2) && l == ATOMS(n, 3) ){
         return n;
      }
   }

   return ErrorReturnValue;
}

/* 
 * We try to match i,j,k against the first 3 in a quartet, then
 * against the next three. Then we check if there's a free slot
 * in the fourth component and match the middle two components
 *
 * The last bit will not be needed as long as gromacs generates
 * all torsions. But I think there are force fields that don't 
 * use all torsions.
 *
 * @return ErrorReturnValue if error; else particle index
 *
 **/

int BrookBonded::_matchAngle( int i, int j, int k, int nbondeds, 
                              int *particles, int *flag ){

// ---------------------------------------------------------------------------------------

   int n;
   static const std::string methodName = "BrookBonded::_matchAngle";

// ---------------------------------------------------------------------------------------

   // validate i > k

   if( i > k ){
      if( getLog() ){
         (void) fprintf( getLog(), "%s Invalid triplet %d-%d-%d\n", methodName.c_str(), i, j, k );
         (void) fflush( getLog() );
      }
      std::stringstream message;
      message << methodName << " invalid triplet: " << i << "-" << j << "-" << k  << std::endl;
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }
   
   for( n = 0; n < nbondeds; n++ ){
      if( j == ATOMS(n, 1) ){
         if(   (i == ATOMS(n, 0) && k == ATOMS(n, 2)) 
             || (i == ATOMS(n, 2) && k == ATOMS(n, 0)) )  {
            *flag = 0;
            return n;
         }
      }
      else if( j == ATOMS(n, 2) ){
         if( i == ATOMS(n, 1) ){
            if( ATOMS(n, 3) == -1 || k == ATOMS(n, 3) ){
               ATOMS(n, 3) = k;
               *flag = 1;
               return n;
            }
         }
         else if( k == ATOMS(n, 1) ){
            if( ATOMS(n, 3) == -1 || i == ATOMS(n, 3) ){
               ATOMS(n, 3) = i;
               *flag = 1;
               return n;
            }
         }
      }
   }
   return ErrorReturnValue;
}

/*flag = 0 means match i-j, 1 means j-k and 2 means k-l
 * Again, like above, if we have uninitialized slots, we take em
 *
 * */

int BrookBonded::_matchBond( int i, int j, int nbondeds, int *particles, int *flag ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_matchBond";

// ---------------------------------------------------------------------------------------

   if( i > j ){
      if( getLog() ){
         (void) fprintf( getLog(), "%s invalid bond %d-%d\n", methodName.c_str(), i, j );
         (void) fflush( getLog() );
      }
      std::stringstream message;
      message << methodName << " invalid bond " << i << "-" << j << std::endl;
      throw OpenMMException( message.str() );
      return ErrorReturnValue;
   }

   for( int n = 0; n < nbondeds; n++ ){

      if( ( i == ATOMS(n, 0) && j == ATOMS(n, 1) ) || ( i == ATOMS(n, 1) && j == ATOMS(n, 0) ) ){

         *flag = 0;
         return n;
      }
      
      //Try second spot

      if(  ATOMS(n, 2) == -1 ){ //One available spot
         if( i == ATOMS(n, 1) ){
            ATOMS(n, 2) = j;
            *flag = 1;
            return n;
         }
         if( j == ATOMS(n, 1) ){
            ATOMS(n, 2) = i;
            *flag = 1;
            return n;
         }
      } else if(  ( i == ATOMS(n, 1) && j == ATOMS(n, 2) ) ||( i == ATOMS(n, 2) && j == ATOMS(n, 1) ) ){
         *flag = 1;
         return n;
      }
      
      //Try third spot
      if( ATOMS(n, 3) == -1 ){
         if( i == ATOMS(n, 2) ){
            ATOMS(n, 3) = j;
            *flag = 2;
            return n;
         }
         if( j == ATOMS(n, 2) ){
            ATOMS(n, 3) = i;
            *flag = 2;
            return n;
         }
      } else if(  ( i == ATOMS(n, 2) && j == ATOMS(n, 3) ) ||( i == ATOMS(n, 3) && j == ATOMS(n, 2) ) ){
         *flag = 2;
         return n;
      }
   }

   return ErrorReturnValue;
}

int BrookBonded::_matchPair( int i, int j, int nbondeds, int *particles ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName = "BrookBonded::_matchPair";

// ---------------------------------------------------------------------------------------

   for( int n = 0; n < nbondeds; n++ ){

      //If there is an empty slot available:

      if( ATOMS(n, 0) == -1 ){

         //If one of i,j matches the l particle

         if( ATOMS(n, 3) == i ){
            ATOMS(n, 0) = j;
            return n;
         }
         if( ATOMS(n, 3) == j ){
            ATOMS(n, 0) = i;
            return n;
         }
      }

      //If the l-particle is available
      if( ATOMS(n, 3) == -1 ){
         if( ATOMS(n, 0) == i ){
            ATOMS(n, 3) = j;
            return n;
         }
         if( ATOMS(n, 0) == j ){
            ATOMS(n, 3) = i;
            return n;
         }
      }

      //Both are unavailable, both much match
      if(   ( i == ATOMS(n, 0) && j == ATOMS(n, 3) ) 
          || ( i == ATOMS(n, 3) && j == ATOMS(n, 0) ) ){
         return n;
      }
   }
   return ErrorReturnValue;
}

/** 
 * Setup Ryckaert-Bellemans parameters/particle indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param particles                 array of particle indices
 * @param params                    arrays of bond parameters
 * @param rbTorsionIndices          the four particles connected by each Ryckaert-Bellemans torsion term
 * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::_addRBTorsions( int *nbondeds, int *particles, float *params[], 
                                 const vector<vector<int> >& rbTorsionIndices, 
                                 const vector<vector<double> >& rbTorsionParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_addRBTorsions";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s number of bonds=%d\n", methodName.c_str(), rbTorsionIndices.size() );
   }

   for( unsigned int ii = 0; ii < rbTorsionIndices.size(); ii++ ){

      vector<int> particlesIndices  = rbTorsionIndices[ii];
      vector<double> rbParameters   = rbTorsionParameters[ii];

      int index = 0;
      int i     = particlesIndices[index++];
      int j     = particlesIndices[index++];
      int k     = particlesIndices[index++];
      int l     = particlesIndices[index++];
      
      int ibonded = _matchTorsion( i, j, k, l, *nbondeds, particles );
      
      if( ibonded < 0 ){
         ibonded             = *nbondeds;
         ATOMS( ibonded, 0 ) = i;
         ATOMS( ibonded, 1 ) = j;
         ATOMS( ibonded, 2 ) = k;
         ATOMS( ibonded, 3 ) = l;
         (*nbondeds)++;
      }

      // note -- we are starting w/ index c1 not c0 -- c0 is not used
      // to calculate the forces!

      index                   = 1;
      PARAMS( ibonded, 0, 0 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 0, 1 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 0, 2 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 0, 3 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 1, 0 ) = (BrookOpenMMFloat) rbParameters[index++];

      if( debug && getLog() ){
         (void) fprintf( getLog(), "   %d [%d %d %d %d] %.3e %.3e %.3e %.3e\n", ibonded, i, j, k, l,
                         rbParameters[0], rbParameters[1],
                         rbParameters[2], rbParameters[3], rbParameters[4], rbParameters[5] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup periodic torsion parameters/particle indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param particles                 array of particle indices
 * @param params                    arrays of bond parameters
 * @param periodicTorsionIndices    the four particles connected by each periodic torsion term
 * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::_addPTorsions( int *nbondeds, int *particles, BrookOpenMMFloat* params[],
                                const vector<vector<int> >& periodicTorsionIndices, 
                                const vector<vector<double> >& periodicTorsionParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_addPTorsion";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s npdih=%d\n", methodName.c_str(), periodicTorsionIndices.size() );
   }

   for( unsigned int ii = 0; ii < periodicTorsionIndices.size(); ii++ ){

      vector<int> particlesIndices      = periodicTorsionIndices[ii];
      vector<double> pTParameters   = periodicTorsionParameters[ii];

      int index = 0;
      int i     = particlesIndices[index++];
      int j     = particlesIndices[index++];
      int k     = particlesIndices[index++];
      int l     = particlesIndices[index++];

      int ibonded = _matchTorsion( i, j, k, l, *nbondeds, particles );
      
      if( ibonded < 0 ){
         ibonded = *nbondeds;
         ATOMS( ibonded, 0 ) = i;
         ATOMS( ibonded, 1 ) = j;
         ATOMS( ibonded, 2 ) = k;
         ATOMS( ibonded, 3 ) = l;
         (*nbondeds)++;
      }

      // note: parameters 0 & 2 switched

      PARAMS( ibonded, 1, 1 ) = (BrookOpenMMFloat) pTParameters[2];
      PARAMS( ibonded, 1, 2 ) = (BrookOpenMMFloat) pTParameters[1];
      PARAMS( ibonded, 1, 3 ) = (BrookOpenMMFloat) pTParameters[0];

      if( debug && getLog() ){
         (void) fprintf( getLog(), "   %d [%d %d %d %d] %.3e %.3e %.3e\n", ibonded, i, j, k, l,
                         pTParameters[0], pTParameters[1], pTParameters[2] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup angle bond parameters/particle indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param particles                 array of particle indices
 * @param params                    arrays of bond parameters
 * @param angleIndices              the angle bond particle indices
 * @param angleParameters           the angle parameters (angle in radians, force constant)
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::_addAngles( int *nbondeds, int *particles, float *params[], const std::vector<std::vector<int> >& angleIndices,
                             const std::vector<std::vector<double> >& angleParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_addAngles";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s nang=%d\n", methodName.c_str(), angleIndices.size() );
   }

   // loop over bonds

   for( unsigned int ii = 0; ii < angleIndices.size(); ii++ ){

      vector<int> particlesIndices   = angleIndices[ii];
      vector<double> angParameters   = angleParameters[ii];

      int index = 0;
      int i     = particlesIndices[index++];
      int j     = particlesIndices[index++];
      int k     = particlesIndices[index++];
 
      int flag;
      int ibonded = _matchAngle( i, j, k, *nbondeds, particles, &flag );
      if( ibonded < 0 ){
         ibonded = *nbondeds;
         ATOMS( ibonded, 0 ) = i;
         ATOMS( ibonded, 1 ) = j;
         ATOMS( ibonded, 2 ) = k;
         flag = 0;
         (*nbondeds)++;
      }
      PARAMS( ibonded, 2, flag*2 )     = (BrookOpenMMFloat) angParameters[0];
      PARAMS( ibonded, 2, flag*2 + 1 ) = (BrookOpenMMFloat) angParameters[1];
      
      if( debug && getLog() ){
         (void) fprintf( getLog(), "   %d [%d %d %d ] %.6e %.6e\n", ibonded, i, j, k,
                         angParameters[0], angParameters[1] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup harmonic bond parameters/particle indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param particles                 array of particle indices
 * @param params                    arrays of bond parameters
 * @param bondIndices               two harmonic bond particle indices
 * @param bondParameters            the force parameters (distance, k)
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::_addBonds( int *nbondeds, int *particles, float *params[], const vector<vector<int> >& bondIndices,
                            const vector<vector<double> >& bondParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_addBonds";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s nbonds=%d\n", methodName.c_str(), bondIndices.size() );
   }

   // loop over bonds

   for( unsigned int ii = 0; ii < bondIndices.size(); ii++ ){

      vector<int> particlesIndices   = bondIndices[ii];
      vector<double> bndParameters   = bondParameters[ii];

      int index = 0;
      int i     = particlesIndices[index++];
      int j     = particlesIndices[index++];

      // insure i < j

      if( i > j ){
         int k = i;
         i     = j;
         j     = k;
      }

      int flag;
      int ibonded = _matchBond( i, j, *nbondeds, particles, &flag );
int saveIbond = ibonded;
      if( ibonded < 0 ){
         ibonded = *nbondeds;
         ATOMS( ibonded, 0 ) = i;
         ATOMS( ibonded, 1 ) = j;
         flag = 0;
         (*nbondeds)++;
      }

      if( flag < 2 ){
         PARAMS( ibonded, 3, flag*2 )     = (BrookOpenMMFloat) bndParameters[0];
         PARAMS( ibonded, 3, flag*2 + 1 ) = (BrookOpenMMFloat) bndParameters[1];
      } else {
         PARAMS( ibonded, 4, 0 )          = (BrookOpenMMFloat) bndParameters[0];
         PARAMS( ibonded, 4, 1 )          = (BrookOpenMMFloat) bndParameters[1];
      }

      if( debug && getLog() ){
         (void) fprintf( getLog(), "   %d (%5d) [%6d %6d ] flag=%2d %.3e %.3e\n", ibonded, saveIbond, i, j, flag,
                         bndParameters[0], bndParameters[1] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup LJ/Coulomb 1-4 parameters/particle indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param particles                 array of particle indices
 * @param params                    arrays of bond parameters
 * @param charges                   array of charges
 * @param bonded14Indices           each element contains the indices of two particles whose nonbonded interactions should be reduced since
 *                                  they form a bonded 1-4 pair
 * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each particle
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::_addPairs( int *nbondeds, int *particles, BrookOpenMMFloat* params[],
                            BrookOpenMMFloat* charges,
                            const std::vector<std::vector<int> >& bonded14Indices,
                            const std::vector<std::vector<double> >& nonbondedParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_addPairs";
   static const double oneSixth        = 1.0/6.0;
   static const int debug              = 1;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s npairs=%d sz=%u %u\n", methodName.c_str(), bonded14Indices.size(), bonded14Indices.size(), nonbondedParameters.size() ); fflush( getLog() );
   }

   for( unsigned int ii = 0; ii < bonded14Indices.size(); ii++ ){

      int i     = bonded14Indices[ii][0];
      int j     = bonded14Indices[ii][1];

      int ibonded = _matchPair( i, j, *nbondeds, particles );
      if( ibonded < 0 ){
         ibonded = *nbondeds;
         ATOMS(ibonded, 0) = i;
         ATOMS(ibonded, 3) = j;
         (*nbondeds)++;
      }

//(void) fprintf( getLog(), "%s   %d %ibonded=%d done\n", methodName.c_str(), ii, ibonded ); fflush( getLog() );

      vector<double> iParameters  = nonbondedParameters[ii];
/*
      double c6                   = iParameters[1];
      double c12                  = 4.0*iParameters[2];
      double sig, eps;
      if( c12 != 0.0 ){
         //eps = c6*c6/c12;
         //sig = pow( c12/c6, oneSixth );
         eps = c12;
         sig = c6;
      } else {
         eps = 0.0;
         sig = 1.0;
      }
*/

      PARAMS( ibonded, 4, 2 ) = (BrookOpenMMFloat) iParameters[1];
      PARAMS( ibonded, 4, 3 ) = (BrookOpenMMFloat) (4.0*iParameters[2]);

      // a little wasteful, but ...

      charges[ibonded] = (BrookOpenMMFloat) iParameters[0];

      if( debug ){
         (void) fprintf( getLog(), "   %d [%d %d ] %.3e %.3e q=%.4f\n", ibonded, i, j, iParameters[1], iParameters[2], charges[ibonded] );
      }
   }

   return DefaultReturnValue;
}

/** 
 * Create and load inverse maps for bonded ixns
 * 
 * @param nbondeds                  number of bonded entries
 * @param nparticles                number of particles
 * @param particles                 arrays of particle indices (particles[numberOfBonds][4])
 * @param platform                  BrookPlatform reference
 * @param log                       log file reference (optional)
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::_loadInvMaps( int nbondeds, int nparticles, int *particles, int particleStreamWidth, int particleStreamSize ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::_loadInvMaps";
   static const int PrintOn            = 0;
   double dangleValue                  = 0.0;

// ---------------------------------------------------------------------------------------

   // get particle stream size

/*
   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (brookPlatform.getDefaultStreamFactory() );
   int particleStreamWidth                      = brookStreamFactory.getDefaultParticleStreamWidth();
   int particleStreamSize                       = brookPlatform.getStreamSize( getNumberOfParticles(), particleStreamWidth, NULL );
*/
   _inverseMapStreamWidth                       = particleStreamWidth;
   
// ---------------------------------------------------------------------------------------

   // allocate temp memory

   float4** invmaps = new float4*[getMaxInverseMapStreamCount()];
   float* block     = new float[4*getMaxInverseMapStreamCount()*particleStreamSize];

   //memset( block, 0, 4*getMaxInverseMapStreamCount()*particleStreamSize*sizeof( float ) );

   float* blockPtr = block;
   for( int ii = 0; ii < getMaxInverseMapStreamCount(); ii++ ){
      invmaps[ii]  = (float4*) blockPtr;
      blockPtr    += 4*particleStreamSize;
   }
   int* counts = new int[particleStreamSize];

// ---------------------------------------------------------------------------------------

   // get inverse maps and load into streams

   // create streams
   // done independently from loading since for test cases some stream counts == 0, but kernels expect stream to
   // have been created even though unused

   // load streams -- initialize all streams even if unused since gather methods will still pick up

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      for( int jj = 0; jj < getMaxInverseMapStreamCount(ii); jj++ ){
         _inverseStreamMaps[ii][jj] = new BrookFloatStreamInternal( BrookCommon::BondedInverseMapStreams, particleStreamSize,
                                                                    particleStreamWidth, BrookStreamInternal::Float4, dangleValue );
      }
   }

   if( PrintOn && getLog() ){
      (void) fprintf( getLog(), "%s force stream strms=%d nbondeds=%d max counts=[%d %d %d %d] strSz&Wd=%d %d\n", methodName.c_str(), getNumberOfForceStreams(),
                      nbondeds, getMaxInverseMapStreamCount(0), getMaxInverseMapStreamCount(1), getMaxInverseMapStreamCount(2), getMaxInverseMapStreamCount(3),
                      particleStreamSize, particleStreamWidth );
      (void) fflush( getLog() );
   }

   // load data

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      for( int jj = 0; jj < 4*getMaxInverseMapStreamCount()*particleStreamSize; jj++ ){
         block[jj] = -1.0f;
      }
      _gpuCalcInvMap( ii, 4, nbondeds, nparticles, particles, getInverseMapStreamCount( ii ), counts, invmaps, &(_inverseMapStreamCount[ii]) );
//gpuPrintInvMaps( _inverseMapStreamCount[ii], nparticles, counts, invmaps, getLog() );
      _validateInverseMapStreamCount( ii, _inverseMapStreamCount[ii] ); 
      for( int jj = 0; jj < _inverseMapStreamCount[ii]; jj++ ){
         _inverseStreamMaps[ii][jj]->loadFromArray( invmaps[jj] );

         if( PrintOn && getLog() ){
            (void) fprintf( getLog(), "%s inverseMap stream strms=%d count=%d index=%d %d InverseMapStreamCount[ii]=%d max=%d\n",
                            methodName.c_str(), getNumberOfForceStreams(), _inverseMapStreamCount[ii], ii, jj,
                            getInverseMapStreamCount( ii ), getMaxInverseMapStreamCount( ii ) );

            for( int kk = 0; kk < particleStreamSize; kk++ ){
               (void) fprintf( getLog(), "%8d [ %.1f %.1f %.1f %.1f]\n", kk, invmaps[jj][kk].x, invmaps[jj][kk].y, invmaps[jj][kk].z, invmaps[jj][kk].w  );
            }
         }

      }    

      // for small systems, must all initialize inverse maps to negative values in order to 
      // keep invalid entries from being included in forces

      if( _inverseMapStreamCount[ii] < getMaxInverseMapStreamCount( ii ) ){
         for( int jj = 0; jj < 4*particleStreamSize; jj++ ){
            block[jj] = -1.0f;
         }
         for( int jj = _inverseMapStreamCount[ii]; jj < getMaxInverseMapStreamCount( ii ); jj++ ){
             _inverseStreamMaps[ii][jj]->loadFromArray( invmaps[0] );
         }
      }    

   }

   // diagnostics 

   if( PrintOn && getLog() ){
      (void) fprintf( getLog(), "%s done\n", methodName.c_str() );
      (void) fflush( getLog() );
   }

   // free memory

   delete[] counts;
   delete[] invmaps[0];
   delete[] invmaps;

   return DefaultReturnValue;
}

/* 
 * Setup for bonded ixns
 *
 * @param numberOfParticles            number of particles
 * @param bondIndices                  vector of vector of harmonic                   bond indices    -- one entry each bond (2 particles     )
 * @param bondParameters               vector of vector of harmonic                   bond parameters -- one entry each bond (2 parameters)
 * @param angleIndices                 vector of vector of angle                      bond indices    -- one entry each bond (3 particles     )
 * @param angleParameters              vector of vector of angle                      bond parameters -- one entry each bond (2 parameters)
 * @param periodicTorsionIndices       vector of vector of periodicTorsionIndices     bond indices    -- one entry each bond (4 particles     )
 * @param periodicTorsionParameters    vector of vector of periodicTorsionParameters  bond parameters -- one entry each bond (3 parameters)
 * @param rbTorsionIndices             vector of vector of rb torsion                 bond indices    -- one entry each bond (4 particles     )
 * @param rbTorsionParameters          vector of vector of rb torsion                 bond parameters -- one entry each bond (5 parameters)
 * @param bonded14Indices              vector of vector of Lennard-Jones 14           particle indices    -- one entry each bond (2 particles     )
 * @param nonbondedParameters          vector of vector of Lennard-Jones 14           parameters      -- one entry each bond (3 parameters)
 * @param platform                     Brook platform reference
 *
 * @return always 1
 *
 * We go through the interactions in the following order
 * torsions
 * angles
 * 14 interactions
 * bonds
 * For each new interaction, we try to fit it in with 
 * those already in. This may not necessarily result in
 * the optimal fit, but should not be too bad. 
 * */

int BrookBonded::setup( int numberOfParticles,
                        BrookBondParameters* harmonicBondBrookBondParameters,
                        BrookBondParameters* harmonicAngleBrookBondParameters,
                        BrookBondParameters* periodicTorsionBrookBondParameters,
                        BrookBondParameters* rbTorsionBrookBondParameters,
                        BrookBondParameters* nonBonded14ForceParameters,  
                        int particleStreamWidth, int particleStreamSize ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::setup";
   static const int PrintOn            = 0;
   double dangleValue                  = 0.0;

// ---------------------------------------------------------------------------------------

   if( PrintOn && getLog() ){
      (void) fprintf( getLog(), "%s particles=%d\n   [%p %p %p %p %p] (bond, angle, pd, rb, 14)\n"
                      "StreamW=%d StreamSz=%d\n", methodName.c_str(), numberOfParticles, harmonicBondBrookBondParameters,
                      harmonicAngleBrookBondParameters,
                      periodicTorsionBrookBondParameters, rbTorsionBrookBondParameters, nonBonded14ForceParameters,
                      particleStreamWidth, particleStreamSize ); fflush( getLog() );
   }

   _numberOfParticles = numberOfParticles;

   // check that particle indices & parameters agree

   // allocate temp memory

   int maxBonds     = 10*numberOfParticles;
   int* particles   = new int[5*maxBonds];
   float* charges   = new BrookOpenMMFloat[maxBonds];

   BrookOpenMMFloat** params = new BrookOpenMMFloat*[getNumberOfParameterStreams()];
   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      params[ii] = new BrookOpenMMFloat[4*maxBonds];
   }

// ---------------------------------------------------------------------------------------

   // build streams

//   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (brookPlatform.getDefaultStreamFactory() );
//   int particleStreamWidth                      = brookStreamFactory.getDefaultParticleStreamWidth();

   // Initialize all particle indices to -1 to indicate empty slots
   // All parameters must be initialized to values that will 
   // produce zero for the corresponding force. 

   memset( charges, 0, maxBonds*sizeof( BrookOpenMMFloat ) );
   for( int ii = 0; ii < maxBonds; ii++ ){

      ATOMS( ii, 0 ) = -1;
      ATOMS( ii, 1 ) = -1;
      ATOMS( ii, 2 ) = -1;
      ATOMS( ii, 3 ) = -1;
   
      for( int jj = 0; jj < getNumberOfParameterStreams(); jj++ ){
         PARAMS( ii, jj, 0 ) = 0.0; 
         PARAMS( ii, jj, 1 ) = 0.0; 
         PARAMS( ii, jj, 2 ) = 0.0; 
         PARAMS( ii, jj, 3 ) = 0.0; 
      }

      //Set sigma negative to indicate no pair
      PARAMS( ii, 4, 2 ) = -1.0;
   }

   // nbondeds tracks number of ixn

   int nbondeds = 0;

   if( rbTorsionBrookBondParameters ){
      _addRBTorsions( &nbondeds, particles, params, rbTorsionBrookBondParameters->getParticleIndices(),         rbTorsionBrookBondParameters->getBondParameters() );
   }
   if( periodicTorsionBrookBondParameters ){
      _addPTorsions(  &nbondeds, particles, params, periodicTorsionBrookBondParameters->getParticleIndices(),   periodicTorsionBrookBondParameters->getBondParameters() );
   }
   if( harmonicAngleBrookBondParameters ){
      _addAngles(     &nbondeds, particles, params, harmonicAngleBrookBondParameters->getParticleIndices(),     harmonicAngleBrookBondParameters->getBondParameters() );
   }
   if( harmonicBondBrookBondParameters ){
      _addBonds(      &nbondeds, particles, params, harmonicBondBrookBondParameters->getParticleIndices(),      harmonicBondBrookBondParameters->getBondParameters() );
   }
   if( nonBonded14ForceParameters ){
      _addPairs(      &nbondeds, particles, params, charges, nonBonded14ForceParameters->getParticleIndices(), nonBonded14ForceParameters->getBondParameters() );
   }

// ---------------------------------------------------------------------------------------

   // check that number of bonds not too large for memory allocated

   if( nbondeds >= maxBonds ){
      std::stringstream message;
      message << methodName << " number of bonds=" << nbondeds << " is greater than maxBonds=" << maxBonds << " numberOfParticles=" << numberOfParticles;
      throw OpenMMException( message.str() );

   } else if( nbondeds < 1 ){

      // return if no bonds

      (void) fprintf( getLog(), "%s WARNING: particles=%d number of bonds=%d maxBonds=%d\n", methodName.c_str(), numberOfParticles, nbondeds, maxBonds );
      _setupCompleted = 1;
      return DefaultReturnValue;

   } else if( PrintOn && getLog() ){
      (void) fprintf( getLog(), "%s particles=%d number of bonds=%d maxBonds=%d\n", methodName.c_str(), numberOfParticles, nbondeds, maxBonds );
      (void) fflush( getLog() );
   }

// ---------------------------------------------------------------------------------------

   // charge stream

   _chargeStream             = new BrookFloatStreamInternal( BrookCommon::BondedChargeStream, nbondeds, particleStreamWidth,
                                                             BrookStreamInternal::Float, dangleValue );

// ---------------------------------------------------------------------------------------

   // particle indices stream

   _particleIndicesStream    = new BrookFloatStreamInternal( BrookCommon::BondedParticleIndicesStream, nbondeds, particleStreamWidth,
                                                             BrookStreamInternal::Float4, dangleValue );

   int* buffer               = new int[4*_particleIndicesStream->getStreamSize()];
   memset( buffer, 0, sizeof( int )*4*_particleIndicesStream->getStreamSize() );

   int index                 = 0;
   for( int ii = 0; ii < nbondeds; ii++ ){
      for( int jj = 0; jj < 4; jj++ ){
         buffer[index++] = ATOMS( ii, jj );
//(void) fprintf( getLog(), "%s particleIndices %d %d  %d buffer=%d particles=%d\n", methodName.c_str(), ii, jj, index, buffer[index-1], ATOMS( ii, jj ) );
      }
   }
   _particleIndicesStream->loadFromArray( buffer, BrookStreamInternal::Integer ); 
   delete[] buffer;

// ---------------------------------------------------------------------------------------

   _chargeStream->loadFromArray( charges ); 

// ---------------------------------------------------------------------------------------

   // bonded parameters

   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      _bondedParameters[ii]  = new BrookFloatStreamInternal( BrookCommon::BondedParametersStream, nbondeds, particleStreamWidth,
                                                             BrookStreamInternal::Float4, dangleValue );
      _bondedParameters[ii]->loadFromArray( params[ii] );
   }

// ---------------------------------------------------------------------------------------

   // debug stuff

   if( PrintOn && getLog() ){

      (void) fprintf( getLog(), "%s nbondeds=%d strDim [%d %d ] sz=%d\n", methodName.c_str(), nbondeds,
                      _particleIndicesStream->getStreamWidth(),
                      _particleIndicesStream->getStreamHeight(), 
                      _particleIndicesStream->getStreamSize() );

      int kIndex = 0;
      int jIndex = 1;
      int iIndex = 2;
      int lIndex = 3;
      int mIndex = 4;
      
      /*
       * float4 parm0   (rbc[1], rbc[2], rbc[3], rbc[4])
       * float4 parm1   (rbc[5], cp, phi, mult) latter three are for pdih
       * float4 parm2   (theta, k, theta, k ) angles i-j-k and j-k-l
       * float4 parm3   (x0, k, x0, k) bonds i-j, j-k
       * float4 parm4   (x0, k, sig14, eps14) bond k-l, i-l
       */
      
      (void) fprintf( getLog(), "\nParams\n" );
      int index = 0; 
      for( int ii = 0; ii < 4*nbondeds; ii += 4, index++ ){
         //  #define PARAMS(X,Y,Z) (params[(Y)][4*(X) + Z])
         (void) fprintf( getLog(), "\n%4d   [%4d %4d %4d %4d]\n"
                              "     rb[%10.6f %10.6f %10.6f %10.6f %10.6f]\n"
                              "    phi[%10.6f %6.1f %5.1f]\n"
                              "    ang[%6.3f %8.2f] [%6.3f %8.2f]\n"
                              "      h[%8.5f %14.6e] [%8.5f %14.6e] [%8.5f %14.6e]\n"
                              "     14[%15.6e %15.6e] q14=%15.6e\n", index,
                         ATOMS(index, 0), ATOMS(index, 1), ATOMS(index, 2), ATOMS(index, 3),
                         params[kIndex][ii], params[kIndex][ii+1], params[kIndex][ii+2], params[kIndex][ii+3], params[jIndex][ii],
                         params[jIndex][ii+1], params[jIndex][ii+2], params[jIndex][ii+3],
                         params[iIndex][ii], params[iIndex][ii+1], params[iIndex][ii+2], params[iIndex][ii+3], 
                         params[lIndex][ii], params[lIndex][ii+1], params[lIndex][ii+2], params[lIndex][ii+3],
                         params[mIndex][ii], params[mIndex][ii+1], params[mIndex][ii+2], params[mIndex][ii+3], charges[index] );
      }
   }

   // load inverse maps to streams

   _loadInvMaps( nbondeds, getNumberOfParticles(), particles, particleStreamWidth, particleStreamSize );
   
// ---------------------------------------------------------------------------------------

   // free memory

   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      delete[] params[ii];
   }

   delete[] params;
   delete[] particles;
   delete[] charges;

   // initialize output streams

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      _bondedForceStreams[ii] = new BrookFloatStreamInternal( BrookCommon::UnrolledForceStream, nbondeds, particleStreamWidth, 
                                                              BrookStreamInternal::Float3, dangleValue );
   }

   _setupCompleted = 1;

   return DefaultReturnValue;
}

/* 
 * Get contents of object
 *
 *
 * @param level   level of dump
 *
 * @return string containing contents
 *
 * */

std::string BrookBonded::getContentsString( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBonded::getContentsString";

   static const unsigned int MAX_LINE_CHARS = 256;
   char value[MAX_LINE_CHARS];
   static const char* Set                   = "Set";
   static const char* NotSet                = "Not set";

// ---------------------------------------------------------------------------------------

   std::stringstream message;
   std::string tab   = "   ";

#ifdef WIN32
#define LOCAL_SPRINTF(a,b,c) sprintf_s( (a), MAX_LINE_CHARS, (b), (c) );   
#define LOCAL_2_SPRINTF(a,b,c,d) sprintf_s( (a), MAX_LINE_CHARS, (b), (c), (d) );   
#else
#define LOCAL_SPRINTF(a,b,c) sprintf( (a), (b), (c) );   
#define LOCAL_2_SPRINTF(a,b,c,d) sprintf( (a), (b), (c), (d) );   
#endif

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfParticles() );
   message << _getLine( tab, "Number of particles:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getInverseMapStreamWidth() );
   message << _getLine( tab, "Inverse map stream width:", value ); 

/*
   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamWidth() );
   message << _getLine( tab, "Particle stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamHeight() );
   message << _getLine( tab, "Particle stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getParticleStreamSize() );
   message << _getLine( tab, "Particle stream size:", value ); 
*/

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfParameterStreams() );
   message << _getLine( tab, "Number of parameter streams:", value ); 
   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      char description[256];
      (void) LOCAL_SPRINTF( description, "Parameter stream %d", ii );
      message << _getLine( tab, description,  (isParameterStreamSet(ii) ? Set : NotSet) ); 
   }
 
   (void) LOCAL_SPRINTF( value, "%d", getMaxInverseMapStreamCount() )
   message << _getLine( tab, "Max inverseMap count:", value ); 

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      char description[256];
      (void) LOCAL_SPRINTF( description, "Inverse map count %d", ii );
      
      std::string checkInverseMap = checkInverseMapStream( ii );
      (void) LOCAL_2_SPRINTF( value, "%d %s", getInverseMapStreamCount( ii ), checkInverseMap.c_str() );
      message << _getLine( tab, description, value ); 
   }

   message << _getLine( tab, "Log:",                 (getLog()               ? Set : NotSet) ); 
   message << _getLine( tab, "Particle indices stream:", (getParticleIndicesStream() ? Set : NotSet) ); 
   //message << _getLine( tab, "Charge stream:",        (getChargeStream()      ? Set : NotSet) ); 
 
   (void) LOCAL_SPRINTF( value, "%d", getNumberOfForceStreams() );
   message << _getLine( tab, "Number of force streams:", value ); 

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      char description[256];
      (void) LOCAL_SPRINTF( description, "Force stream %d", ii );
      message << _getLine( tab, description,  (isForceStreamSet(ii) ? Set : NotSet) ); 
   }
 
#undef LOCAL_SPRINTF
#undef LOCAL_2_SPRINTF

   return message.str();
}

/*
 * Helper functions for building inverse maps for 
 * torsions, impropers and angles.
 * 
 * For each particle, calculates the positions at which it's
 * forces are to be picked up from and stores the position
 * in the appropriate index.
 *
 * Input: number of torsions, the particle indices, and a flag indicating
 *        whether we're doing i(0), j(1), k(2) or l(3)
 * Output: an array of counts per particle
 *         arrays of inversemaps
 *         nimaps - the number of invmaps actually used.
 *
 * @param posflag       0-niparticles-1
 * @param niparticles       3 for angles, 4 for torsions, impropers
 * @param nints         number of interactions
 * @param nparticles        number of particles
 * @param *particles        gromacs interaction list
 * @param nmaps         maximum number of inverse maps
 * @param   counts[]    output counts of how many places each particle occurs
 * @param *invmaps[]    output array of nmaps inverse maps
 * @param *nimaps,      output max number of inverse maps actually used
 *
 * @return DefaultReturnValue, unless error in which case exits w/ OpenMM exception
 *
 **/

int BrookBonded::_gpuCalcInvMap( int posflag, int niparticles, int nints, int nparticles,
                                 int *particles, int nmaps, int counts[], float4 *invmaps[],
                                 int *nimaps ){

// ---------------------------------------------------------------------------------------

   int i, j;
   int particle;
   int mapnum, mapcomp;

   static const std::string methodName      = "BrookBonded::_gpuCalcInvMap";

   static const unsigned int MAX_LINE_CHARS = 256;
   //char value[MAX_LINE_CHARS];
   static const char* Set                   = "Set";
   static const char* NotSet                = "Not set";
   static const int PrintOn                 = 0;

// ---------------------------------------------------------------------------------------

   memset( counts, 0, sizeof( int )*nparticles );

   for( i = 0; i < nmaps; i++ ){
      for( j = 0; j < nparticles; j++ ){
         invmaps[i][j] = float4( -1.0, -1.0, -1.0, -1.0 );
      }
   }
   
   //This will hold the number of imaps actually used

   *nimaps = -1;

   //Now note down the positions where each particle occurs

   if( PrintOn && getLog() ){
      (void) fprintf( getLog(), "%s: pos=%d ni=%d nints=%d nparticles=%d nmaps=<%d>\n", methodName.c_str(), posflag, niparticles, nints, nparticles, nmaps ); 
      (void) fflush( getLog() );
   }

int particleRange[2]   = { 90000000, -90000000 };
int mapnumRange[2] = { 90000000, -90000000 };

   for(  i = 0; i < nints; i++ ){
      //This is our particle
      particle = particles[ (niparticles + 1) * i + posflag + 1 ];

      //Special for merged bondeds
      if ( particle == -1 ){
         continue;
      }

if( particle < particleRange[0] ){
   particleRange[0] = particle;
}
if( particle > particleRange[1] ){
   particleRange[1] = particle;
}
      //Check to make sure we're inside the limits
      if ( counts[particle] > nmaps * 4 ){
         if( PrintOn && getLog() ){
            (void) fprintf( getLog(), "%s Particle %d has too many proper torsions (%d, max %d)\n",
                            methodName.c_str(), particle, counts[particle], nmaps*4 );
            (void) fflush( getLog() );
         }
         std::stringstream message;
         message << methodName << " Particle " << particle << " has too many proper torsions; valid range:(" << counts[particle] << ", " << nmaps*4 << ")";
         throw OpenMMException( message.str() );
      }
      
      //Which invmap will this go into

      mapnum = counts[particle] / 4;

      if ( mapnum > *nimaps )
         *nimaps = mapnum;

      //Which component will it be
      mapcomp = counts[particle] % 4;

      //Set it
      //This is silly, but otherwise I have to declare it as float*
      //and things get even more confusing. :)
      switch (mapcomp){
         case 0: invmaps[mapnum][particle].x = (float) i; break;
         case 1: invmaps[mapnum][particle].y = (float) i; break;
         case 2: invmaps[mapnum][particle].z = (float) i; break;
         case 3: invmaps[mapnum][particle].w = (float) i; break;
         default:
            if( PrintOn && getLog() ){
               (void) fprintf( getLog(), "mapcomp %d invalid -- impossible!\n", mapcomp );
               (void) fflush( getLog() );
            }
            std::stringstream message;
            message << methodName << " mapcomp " << mapcomp << " invalid -- actually impossible!";
            throw OpenMMException( message.str() );
            break;
      }
      
      counts[particle]++;

if( mapnum < mapnumRange[0] ){
   mapnumRange[0] = mapnum;
}
if( mapnum > mapnumRange[1] ){
   mapnumRange[1] = mapnum;
}

//fprintf( gpu->log, "%d particle=%d  mapcomp=%d counts[]=%d mapnum=%d\n", i, particle, mapcomp, counts[particle], mapnum );

   }

   (*nimaps)++;

if( PrintOn && getLog() ){
   (void) fprintf( getLog(), "%s mnmaps=%d Ranges: particle [%d %d] mapnum [%d %d]\n",
                   methodName.c_str(), *nimaps, particleRange[0], particleRange[1], mapnumRange[0], mapnumRange[1] );
   (void) fflush( getLog() );
}

   return DefaultReturnValue;	
}


void BrookBonded::_gpuPrintInvMaps( int nmaps, int nparticles, int counts[], float4 *invmap[], FILE* logFile ){
   int i;
   int j;
   for(  i = 0; i < nparticles; i++ ){
      fprintf( logFile, "%d %d ", i, counts[i] );
      for(  j = 0; j < nmaps; j++ ){
         fprintf( logFile, "%6.0f %6.0f %6.0f %6.0f", invmap[j][i].x, invmap[j][i].y, 
                  invmap[j][i].z, invmap[j][i].w );
      }
      fprintf( logFile, "\n");
   }
}

/* We are still plagued by kernel call overheads. This is for a big fat
 * merged inverse gather kernel:
 * Since we have 32 bit floats, we have 23 bits of mantissa or the largest
 * integer we can represent is 2^23. So it should be quite safe to add 
 * 100000 * n to the index where n is the stream in which we should do the
 * lookup. This assumes that nints < 100000, preferably nints << 100000
 * which should always be true
 * */

int BrookBonded::_gpuCalcInvMap_merged( 
      int nints,    //number of interactions
      int nparticles,   //number of particles
      int *particles,   //ijkl,ijkl,ijkl...
      int nmaps,      //maximum number of inverse maps
        int counts[],   //output counts of how many places each particle occurs
      float4 *invmaps[], //output array of nmaps inverse maps
      int *nimaps        //output max number of inverse maps actually used
      ){
   int i, j;
   int particle;
   int mapnum, mapcomp;
   int pos;
   
   for(  i = 0; i < nparticles; i++ )
      counts[i] = 0;

   for(  i = 0; i < nmaps; i++ ){
      for(  j = 0; j < nparticles; j++ ){
         invmaps[i][j] = float4( -1.0, -1.0, -1.0, -1.0 );
      }
   }

   //This will hold the number of imaps actually used
   *nimaps = -1;

   //For each particle
   for(  i = 0; i < nints; i++ ){
      for(  j = 0; j < 4; j++ ){
         
         particle = particles[ i * 4 + j ];
         
         if ( particle == -1 ){
         	//Nothing to be done for this particle, go to next
         	continue;
         }
         
         //Which map
         mapnum = counts[ particle ] / 4;
         
         //Make sure we have space
         if ( mapnum >= nmaps ){
         	printf( "Particle %d has too many bondeds(%d, max %d)\n",
         			 particle, counts[particle], nmaps * 4 );
         	return 0;
         }
         	
         if ( mapnum > *nimaps ){
         	*nimaps = mapnum;
         }

         //Which component
         mapcomp = counts[ particle ] % 4;
         
         //Encode target stream and position
         pos = 100000 * j + i;

         switch ( mapcomp ){
         	case 0: invmaps[mapnum][particle].x = (float) pos; break;
         	case 1: invmaps[mapnum][particle].y = (float) pos; break;
         	case 2: invmaps[mapnum][particle].z = (float) pos; break;
         	case 3: invmaps[mapnum][particle].w = (float) pos; break;
         }

         counts[ particle ]++;

      }
   }
   
   (*nimaps)++;
   return 1;
}

/* Repacks the invmap streams for more efficient access in the
 * merged inverse gather kernel
 *
 * buf should be nimaps * nparticles large.
 * */
int BrookBonded::_gpuRepackInvMap_merged( int nparticles, int nmaps, int *counts, 
                                          float4 *invmaps[], float4 *buf ){
   int i, j;
   int nmaps_i;

   for(  i = 0; i < nparticles; i++ ){
      for(  j = 0; j < nmaps; j++ ){
         buf[ i + j*nparticles ] = float4( -1.0f, -1.0f, -1.0f, -1.0f );
      }
   }
   
   for(  i = 0; i < nparticles; i++ ){
      
      nmaps_i = counts[i] / 4;
      if ( counts[i] % 4 ) 
         nmaps_i++;
      
      for(  j = 0; j < nmaps_i; j++ ){
         buf[ i + j * nparticles ] = invmaps[j][i];
      }
   }
   return 1;
}

/** 
 * Compute forces
 * 
 */

void BrookBonded::computeForces( BrookStreamImpl& positionStream, BrookStreamImpl& forceStream ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBonded::computeForces";

   static const int PrintOn                 = 0;

   static const int I_Stream                = 0;
   static const int J_Stream                = 1;
   static const int K_Stream                = 2;
   static const int L_Stream                = 3;

   static const int MaxErrorMessages        = 2;
   static       int ErrorMessages           = 0;

   static const float4 dummyParameters( 0.0, 0.0, 0.0, 0.0 );

   // static const int debug                   = 1;

// ---------------------------------------------------------------------------------------

   // bonded

   float epsfac                                        = (float) (getCoulombFactor());
   float width                                         = (float) (getInverseMapStreamWidth());

   // bonded forces

   BrookFloatStreamInternal**  bondedParameters        = getBondedParameterStreams();
   BrookFloatStreamInternal**  bondedForceStreams      = getBondedForceStreams();

   BrookFloatStreamInternal**  inverseStreamMaps[4];
   inverseStreamMaps[0]                                = getInverseStreamMapsStreams( 0 );
   inverseStreamMaps[1]                                = getInverseStreamMapsStreams( 1 );
   inverseStreamMaps[2]                                = getInverseStreamMapsStreams( 2 );
   inverseStreamMaps[3]                                = getInverseStreamMapsStreams( 3 );

   kbonded_CDLJ( epsfac, 
                 (float) bondedForceStreams[0]->getStreamWidth(),
                 dummyParameters,
                 positionStream.getBrookStream(),
                 getChargeStream()->getBrookStream(),
                 getParticleIndicesStream()->getBrookStream(),
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

      (void) fprintf( getLog(), "\nPost kbonded_CDLJ: epsFac=%.6f %.6f", epsfac, getCoulombFactor());
      (void) fprintf( getLog(), "\nParticle indices stream\n" );
      getParticleIndicesStream()->printToFile( getLog() );

      (void) fprintf( getLog(), "\nCharge stream\n" );
      getChargeStream()->printToFile( getLog() );

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

      (void) fprintf( getLog(), "\nInverse map streams\n" );
      for( int ii = 0; ii < 4; ii++ ){
         (void) fprintf( getLog(), "\nInverse map streams -- StreamIndex=%d cnt=%d\n", ii, getInverseMapStreamCount( ii ) );
         for( int jj = 0; jj < countPrintInvMap[ii]; jj++ ){
            (void) fprintf( getLog(), "\n   Inverse map streams index=%d %d\n", ii, jj );
            inverseStreamMaps[ii][jj]->printToFile( getLog() );
         }
      }
   }

   // gather forces

   if( getInverseMapStreamCount( I_Stream ) == 3 && getInverseMapStreamCount( K_Stream ) == 3 ){

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

   } else if( getInverseMapStreamCount( I_Stream ) == 2 && getInverseMapStreamCount( K_Stream ) == 2 ){

      kinvmap_gather2_2( width,

                         inverseStreamMaps[I_Stream][0]->getBrookStream(),
                         inverseStreamMaps[I_Stream][1]->getBrookStream(),
                         bondedForceStreams[I_Stream]->getBrookStream(),

                         inverseStreamMaps[K_Stream][0]->getBrookStream(),
                         inverseStreamMaps[K_Stream][1]->getBrookStream(),
                         bondedForceStreams[K_Stream]->getBrookStream(),

                         forceStream.getBrookStream(), forceStream.getBrookStream() );

   } else if( getInverseMapStreamCount( I_Stream ) == 3 && getInverseMapStreamCount( K_Stream ) == 4 ){

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

   } else if( getInverseMapStreamCount( I_Stream ) == 3 && getInverseMapStreamCount( K_Stream ) == 5 ){

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

   } else if( getInverseMapStreamCount( I_Stream ) == 1 && getInverseMapStreamCount( K_Stream ) == 1 ){

      kinvmap_gather1_1( width,
                         inverseStreamMaps[I_Stream][0]->getBrookStream(),
                         bondedForceStreams[I_Stream]->getBrookStream(),
                         inverseStreamMaps[K_Stream][0]->getBrookStream(),
                         bondedForceStreams[K_Stream]->getBrookStream(),
                         forceStream.getBrookStream(), forceStream.getBrookStream() );

   } else {

      // case not handled -- throw an exception

      if( getLog() && ErrorMessages++ < MaxErrorMessages && getInverseMapStreamCount( I_Stream ) > 0 && getInverseMapStreamCount( K_Stream ) > 0 ){
         (void) fprintf( getLog(), "%s case: I-map=%d K-map=%d -- not handled.\n",
                          methodName.c_str(), getInverseMapStreamCount( I_Stream ),
                                              getInverseMapStreamCount( K_Stream ) );
         (void) fflush(  getLog() );
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
      message << methodName << "I-maps=" << getInverseMapStreamCount( I_Stream ) << " and " << 
                               "K-maps=" << getInverseMapStreamCount( K_Stream ) << " not handled.";
      throw OpenMMException( message.str() );
*/

   }

   // diagnostics

   if( 0 && PrintOn ){

      (void) fprintf( getLog(), "\nPost 3_4/3_5 && NB forces" );
      BrookStreamInternal* brookStreamInternalF   = forceStream.getBrookStreamImpl();
      brookStreamInternalF->printToFile( getLog() );

   }

   if( getInverseMapStreamCount( J_Stream ) == 1 && getInverseMapStreamCount( L_Stream ) == 1 ){

      kinvmap_gather1_1( width,
                         inverseStreamMaps[J_Stream][0]->getBrookStream(),
                         bondedForceStreams[J_Stream]->getBrookStream(),
                         inverseStreamMaps[L_Stream][0]->getBrookStream(),
                         bondedForceStreams[L_Stream]->getBrookStream(),
                         forceStream.getBrookStream(), forceStream.getBrookStream() );
   
   } else if( getInverseMapStreamCount( J_Stream ) == 5 && getInverseMapStreamCount( L_Stream ) == 2 ){

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

      if( getLog() && ErrorMessages++ < MaxErrorMessages && getInverseMapStreamCount( J_Stream ) > 0 && getInverseMapStreamCount( L_Stream ) > 0 ){
         (void) fprintf( getLog(), "%s case: J-map=%d L-map=%d -- not handled.\n",
                          methodName.c_str(), getInverseMapStreamCount( J_Stream ),
                                              getInverseMapStreamCount( L_Stream ) );
         (void) fflush(  getLog() );
      }

      // this is for testing purposes a-- may need to be cleaned 
      // or add new gather functions

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
      message << methodName << "J-maps=" << getInverseMapStreamCount( J_Stream ) << " and " << 
                               "L-maps=" << getInverseMapStreamCount( L_Stream ) << " not handled.";
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
      for( int ii = 0; ii < _brookNonBonded->getNumberOfParticles()*3; ii += 3 ){
         (void) fprintf( getLog(), "%d [%.6e %.6e %.6e]\n", ii, data[ii], data[ii+1], data[ii+2] );
      }
*/

   }

   // ---------------------------------------------------------------------------------------
}
