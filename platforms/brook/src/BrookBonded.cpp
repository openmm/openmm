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
#include "gpu/invmap.h"
#include "gpu/kforce.h"

using namespace OpenMM;
using namespace std;

#define ATOMS(X,Y) (atoms[ 5*(X) + (Y) + 1 ])
#define PARAMS(X,Y,Z) (params[(Y)][4*(X) + Z])

/** 
 * BrookBonded constructor
 * 
 */

BrookBonded::BrookBonded( ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::BrookBonded";

// ---------------------------------------------------------------------------------------

   _numberOfAtoms             = 0;
   _ljScale                   = 1.0;
   _coulombFactor             = 332.0;

   _atomIndicesStream         = NULL;
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

   _invMapStreamWidth = -1;
}   
 
/** 
 * BrookBonded destructor
 * 
 */

BrookBonded::~BrookBonded( ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName = "BrookBonded::~BrookBonded";

// ---------------------------------------------------------------------------------------

   delete _atomIndicesStream;
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
 * Get width of inverse map streams
 * 
 * @return  width of inverse map streams
 *
 */

int BrookBonded::getInvMapStreamWidth( void ) const {
   return _invMapStreamWidth;
}

/** 
 * Get LJ 14 scaling parameter
 * 
 * @return LJ 14 scaling parameter
 *
 */

BrookOpenMMFloat BrookBonded::getLJ_14Scale( void ) const {
   return _ljScale;
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

int BrookBonded::validateInverseMapStreamCount( int index, int count ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::validateInverseMapStreamCount";

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
 * Get bonded atom indices stream
 * 
 * @return  atom indices stream
 *
 */

BrookFloatStreamImpl* BrookBonded::getAtomIndicesStream( void ) const {

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::getAtomIndicesStream";

// ---------------------------------------------------------------------------------------

   return _atomIndicesStream;
}

/** 
 * Get array of bonded parameter streams
 * 
 * @return  array of bonded parameter streams
 *
 */

BrookFloatStreamImpl** BrookBonded::getBondedParameterStreams( void ){

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

BrookFloatStreamImpl** BrookBonded::getBondedForceStreams( void ){

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

BrookFloatStreamImpl** BrookBonded::getInverseStreamMapsStreams( int index ){

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

void BrookBonded::flipQuartet( int ibonded, int *atoms ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::flipQuartet";

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

int BrookBonded::matchTorsion( int i, int j, int k, int l, int nbondeds, int *atoms ){

// ---------------------------------------------------------------------------------------

   // static const std::string methodName = "BrookBonded::matchTorsion";

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
 * all dihedrals. But I think there are force fields that don't 
 * use all dihedrals.
 *
 * @return ErrorReturnValue if error; else atom index
 *
 **/

int BrookBonded::matchAngle( int i, int j, int k, int nbondeds, 
                             int *atoms, int *flag ){

// ---------------------------------------------------------------------------------------

   int n;
   static const std::string methodName = "BrookBonded::matchAngle";

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

int BrookBonded::matchBond( int i, int j, int nbondeds, int *atoms, int *flag ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::matchBond";

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
      }
      else if(  ( i == ATOMS(n, 1) && j == ATOMS(n, 2) ) ||( i == ATOMS(n, 2) && j == ATOMS(n, 1) ) ){
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

int BrookBonded::matchPair( int i, int j, int nbondeds, int *atoms ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName = "BrookBonded::matchPair";

// ---------------------------------------------------------------------------------------

   for( int n = 0; n < nbondeds; n++ ){

      //If there is an empty slot available:

      if( ATOMS(n, 0) == -1 ){

         //If one of i,j matches the l atom

         if( ATOMS(n, 3) == i ){
            ATOMS(n, 0) = j;
            return n;
         }
         if( ATOMS(n, 3) == j ){
            ATOMS(n, 0) = i;
            return n;
         }
      }

      //If the l-atom is available
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
 * Setup Ryckaert-Bellemans parameters/atom indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param atoms                     array of atom indices
 * @param params                    arrays of bond parameters
 * @param rbTorsionIndices          the four atoms connected by each Ryckaert-Bellemans torsion term
 * @param rbTorsionParameters       the coefficients (in order of increasing powers) for each Ryckaert-Bellemans torsion term
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::addRBDihedrals( int *nbondeds, int *atoms, float *params[], 
                                 const vector<vector<int> >& rbTorsionIndices, 
                                 const vector<vector<double> >& rbTorsionParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::BrookBonded";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s number of bonds=%d\n", methodName.c_str(), rbTorsionIndices.size() );
   }

   for( unsigned int ii = 0; ii < rbTorsionIndices.size(); ii++ ){

      vector<int> atomsIndices      = rbTorsionIndices[ii];
      vector<double> rbParameters   = rbTorsionParameters[ii];

      int index = 0;
      int i     = atomsIndices[index++];
      int j     = atomsIndices[index++];
      int k     = atomsIndices[index++];
      int l     = atomsIndices[index++];
      
      int ibonded = matchTorsion( i, j, k, l, *nbondeds, atoms );
      
      if( ibonded < 0 ){
         ibonded             = *nbondeds;
         ATOMS( ibonded, 0 ) = i;
         ATOMS( ibonded, 1 ) = j;
         ATOMS( ibonded, 2 ) = k;
         ATOMS( ibonded, 3 ) = l;
         (*nbondeds)++;
      }

      index                   = 0;
      PARAMS( ibonded, 0, 0 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 0, 1 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 0, 2 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 0, 3 ) = (BrookOpenMMFloat) rbParameters[index++];
      PARAMS( ibonded, 1, 0 ) = (BrookOpenMMFloat) rbParameters[index++];

      if( debug && getLog() ){
         (void) fprintf( getLog(), "   %d [%d %d %d %d] %.3e %.3e %.3e %.3e\n", ibonded, i, j, k, l,
                         rbParameters[0], rbParameters[1],
                         rbParameters[2], rbParameters[3], rbParameters[4] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup periodic torsion parameters/atom indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param atoms                     array of atom indices
 * @param params                    arrays of bond parameters
 * @param periodicTorsionIndices    the four atoms connected by each periodic torsion term
 * @param periodicTorsionParameters the force parameters (k, phase, periodicity) for each periodic torsion term
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::addPDihedrals( int *nbondeds, int *atoms, BrookOpenMMFloat* params[],
                                const vector<vector<int> >& periodicTorsionIndices, 
                                const vector<vector<double> >& periodicTorsionParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::addPDihedrals";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s npdih=%d\n", methodName.c_str(), periodicTorsionIndices.size() );
   }

   for( unsigned int ii = 0; ii < periodicTorsionIndices.size(); ii++ ){

      vector<int> atomsIndices      = periodicTorsionIndices[ii];
      vector<double> pTParameters   = periodicTorsionParameters[ii];

      int index = 0;
      int i     = atomsIndices[index++];
      int j     = atomsIndices[index++];
      int k     = atomsIndices[index++];
      int l     = atomsIndices[index++];

      int ibonded = matchTorsion( i, j, k, l, *nbondeds, atoms );
      
      if( ibonded < 0 ){
         ibonded = *nbondeds;
         ATOMS( ibonded, 0 ) = i;
         ATOMS( ibonded, 1 ) = j;
         ATOMS( ibonded, 2 ) = k;
         ATOMS( ibonded, 3 ) = l;
         (*nbondeds)++;
      }
      PARAMS( ibonded, 1, 1 ) = (BrookOpenMMFloat) pTParameters[0];
      PARAMS( ibonded, 1, 2 ) = (BrookOpenMMFloat) pTParameters[1];
      PARAMS( ibonded, 1, 3 ) = (BrookOpenMMFloat) pTParameters[2];

      if( debug && getLog() ){
         (void) fprintf( getLog(), "   %d [%d %d %d %d] %.3e %.3e %.3e\n", ibonded, i, j, k, l,
                         pTParameters[0], pTParameters[1], pTParameters[2] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup angle bond parameters/atom indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param atoms                     array of atom indices
 * @param params                    arrays of bond parameters
 * @param angleIndices              the angle bond atom indices
 * @param angleParameters           the angle parameters (angle in radians, force constant)
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::addAngles( int *nbondeds, int *atoms, float *params[], const std::vector<std::vector<int> >& angleIndices,
                            const std::vector<std::vector<double> >& angleParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::addAngles";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s nang=%d\n", methodName.c_str(), angleIndices.size() );
   }

   // loop over bonds

   for( unsigned int ii = 0; ii < angleIndices.size(); ii++ ){

      vector<int> atomsIndices       = angleIndices[ii];
      vector<double> angParameters   = angleParameters[ii];

      int index = 0;
      int i     = atomsIndices[index++];
      int j     = atomsIndices[index++];
      int k     = atomsIndices[index++];
 
      int flag;
      int ibonded = matchAngle( i, j, k, *nbondeds, atoms, &flag );
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
 * Setup harmonic bond parameters/atom indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param atoms                     array of atom indices
 * @param params                    arrays of bond parameters
 * @param bondIndices               two harmonic bond atom indices
 * @param bondParameters            the force parameters (distance, k)
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::addBonds( int *nbondeds, int *atoms, float *params[], const vector<vector<int> >& bondIndices,
                           const vector<vector<double> >& bondParameters ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::addBonds";
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s nbonds=%d\n", methodName.c_str(), bondIndices.size() );
   }

   // loop over bonds

   for( unsigned int ii = 0; ii < bondIndices.size(); ii++ ){

      vector<int> atomsIndices       = bondIndices[ii];
      vector<double> bndParameters   = bondParameters[ii];

      int index = 0;
      int i     = atomsIndices[index++];
      int j     = atomsIndices[index++];

      int flag;
      int ibonded = matchBond( i, j, *nbondeds, atoms, &flag );
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
         (void) fprintf( getLog(), "   %d [%d %d ] flag=%d %.3e %.3e\n", ibonded, i, j, flag,
                         bndParameters[0], bndParameters[1] );
      }
   }

   return DefaultReturnValue;
}

/**
 * Setup LJ/Coulomb 1-4 parameters/atom indices
 * 
 * @param nbondeds                  number of bonded entries
 * @param atoms                     array of atom indices
 * @param params                    arrays of bond parameters
 * @param charges                   array of charges
 * @param bonded14Indices           each element contains the indices of two atoms whose nonbonded interactions should be reduced since
 *                                  they form a bonded 1-4 pair
 * @param nonbondedParameters       the nonbonded force parameters (charge, sigma, epsilon) for each atom
 * @param lj14Scale                 the factor by which van der Waals interactions should be reduced for bonded 1-4 pairs
 * @param log                       log reference
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::addPairs( int *nbondeds, int *atoms, BrookOpenMMFloat* params[],
                           BrookOpenMMFloat* charges,
                           const std::vector<std::vector<int> >& bonded14Indices,
                           const std::vector<std::vector<double> >& nonbondedParameters,
                           double lj14Scale, double coulombScale ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::BrookBonded";
   static const double oneSixth        = 1.0/6.0;
   static const int debug              = 0;

// ---------------------------------------------------------------------------------------

   if( debug && getLog() ){
      (void) fprintf( getLog(), "%s npairs=%d\n", methodName.c_str(), bonded14Indices.size() );
   }

   for( unsigned int ii = 0; ii < bonded14Indices.size(); ii++ ){

      std::vector<int> atomsIndices        = bonded14Indices[ii];

      int index = 0;
      int i     = atomsIndices[index++];
      int j     = atomsIndices[index++];

      int ibonded = matchPair( i, j, *nbondeds, atoms );
      if( ibonded < 0 ){
         ibonded = *nbondeds;
         ATOMS(ibonded, 0) = i;
         ATOMS(ibonded, 3) = j;
         (*nbondeds)++;
      }

      vector<double> iParameters  = nonbondedParameters[i];
      vector<double> jParameters  = nonbondedParameters[j];

      double c6                   = iParameters[0] + jParameters[0];
      double c12                  = lj14Scale*(iParameters[1] * jParameters[1]);

      double sig, eps;
      if( c12 != 0.0 ){
         eps = c6*c6/c12;
         sig = pow( c12/c6, oneSixth );
      } else {
         eps = 0.0;
         sig = 1.0;
      }

      PARAMS( ibonded, 4, 2 ) = (BrookOpenMMFloat) sig;
      PARAMS( ibonded, 4, 3 ) = (BrookOpenMMFloat) eps;

      // a little wasteful, but ...

      charges[i] = (BrookOpenMMFloat) iParameters[2];
      charges[j] = (BrookOpenMMFloat) jParameters[2];

      if( debug ){
         (void) fprintf( getLog(), "\n %d [%d %d ] %.3e %.3e", ibonded, i, j, sig, eps );
      }
   }

   return DefaultReturnValue;
}

/** 
 * Create and load inverse maps for bonded ixns
 * 
 * @param nbondeds                  number of bonded entries
 * @param natoms                    number of atoms
 * @param atoms                     arrays of atom indices (atoms[numberOfBonds][4])
 * @param platform                  BrookPlatform reference
 * @param log                       log file reference (optional)
 *
 * @return nonzero value if error
 *
 */

int BrookBonded::loadInvMaps( int nbondeds, int natoms, int *atoms, const BrookPlatform& brookPlatform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::loadInvMaps";

// ---------------------------------------------------------------------------------------

   // get atom stream size

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (brookPlatform.getDefaultStreamFactory() );
   BrookStreamInfo* brookStreamInfo             = brookStreamFactory.getBrookStreamInfo( BrookStreamFactory::AtomPositions );
   int atomStreamWidth                          = brookStreamInfo->getStreamWidth();
   int atomStreamSize                           = brookPlatform.getStreamSize( getNumberOfAtoms(), atomStreamWidth, NULL );
   _invMapStreamWidth                           = atomStreamWidth;
   
   // allocate temp memory

   float4** invmaps = new float4*[getMaxInverseMapStreamCount()];
   float* block     = new float[4*getMaxInverseMapStreamCount()*atomStreamSize];
   for( int ii = 0; ii < getMaxInverseMapStreamCount(); ii++ ){
      invmaps[ii]  = (float4*) block;
      block       += 4*atomStreamSize;
   }
   int* counts = new int[atomStreamSize];

   // get inverse maps and load into streams

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      gpuCalcInvMap( ii, 4, nbondeds, natoms, atoms, getInverseMapStreamCount( ii ), counts, invmaps, &(_inverseMapStreamCount[ii]) );
      validateInverseMapStreamCount( ii, _inverseMapStreamCount[ii] ); 
      for( int jj = 0; jj < _inverseMapStreamCount[ii]; jj++ ){
         StreamImpl* inverseStreamMaps = brookStreamFactory.createStreamImpl( BrookStreamFactory::BondedInverseMapStreams, atomStreamSize, Stream::Float4, brookPlatform );
         _inverseStreamMaps[ii][jj]    = dynamic_cast<BrookFloatStreamImpl*> ( inverseStreamMaps );
         _inverseStreamMaps[ii][jj]->loadFromArray( invmaps[jj] );
      }
   }

   if( getLog() ){
      (void) fprintf( getLog(), "%s done\n", methodName.c_str() );
      (void) fflush( getLog() );
      //gpuPrintInvMaps( bp->nimaps, natoms, counts, invmaps, gpu->log );
   }

   // free memory

   delete counts;
   delete invmaps[0];
   delete invmaps;

   return DefaultReturnValue;
}

/* 
 * Setup for bonded ixns
 *
 * @param numberOfAtoms                number of atoms
 * @param bondIndices                  vector of vector of harmonic                   bond indices    -- one entry each bond (2 atoms     )
 * @param bondParameters               vector of vector of harmonic                   bond parameters -- one entry each bond (2 parameters)
 * @param angleIndices                 vector of vector of angle                      bond indices    -- one entry each bond (3 atoms     )
 * @param angleParameters              vector of vector of angle                      bond parameters -- one entry each bond (2 parameters)
 * @param periodicTorsionIndices       vector of vector of periodicTorsionIndices     bond indices    -- one entry each bond (4 atoms     )
 * @param periodicTorsionParameters    vector of vector of periodicTorsionParameters  bond parameters -- one entry each bond (3 parameters)
 * @param rbTorsionIndices             vector of vector of rb torsion                 bond indices    -- one entry each bond (4 atoms     )
 * @param rbTorsionParameters          vector of vector of rb torsion                 bond parameters -- one entry each bond (5 parameters)
 * @param bonded14Indices              vector of vector of Lennard-Jones 14           atom indices    -- one entry each bond (2 atoms     )
 * @param nonbondedParameters          vector of vector of Lennard-Jones 14           parameters      -- one entry each bond (3 parameters)
 * @param lj14Scale                    scaling factor for 1-4 ixns
 * @param coulombScale                 Coulomb scaling factor for 1-4 ixns
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

int BrookBonded::setup( int numberOfAtoms,
                        const vector<vector<int> >& bondIndices,            const vector<vector<double> >& bondParameters,
                        const vector<vector<int> >& angleIndices,           const vector<vector<double> >& angleParameters,
                        const vector<vector<int> >& periodicTorsionIndices, const vector<vector<double> >& periodicTorsionParameters,
                        const vector<vector<int> >& rbTorsionIndices,       const vector<vector<double> >& rbTorsionParameters,
                        const vector<vector<int> >& bonded14Indices,        const vector<vector<double> >& nonbondedParameters,
                        double lj14Scale, double coulombScale,  const BrookPlatform& platform ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName = "BrookBonded::setup";

// ---------------------------------------------------------------------------------------

   _numberOfAtoms = numberOfAtoms;

   // check that atom indices & parameters agree

   if( bondIndices.size() != bondParameters.size() ){
      std::stringstream message;
      message << methodName << " number of harmonic bond atom indices=" << bondIndices.size() << " does not equal number of harmonic bond parameter entries=" << bondParameters.size();
      throw OpenMMException( message.str() );
   } else if( getLog() ){
      (void) fprintf( getLog(), "%s harmonic bonds=%d\n", methodName.c_str(), bondIndices.size() );
      (void) fflush( getLog() );
   }

   if( angleIndices.size() != angleParameters.size() ){
      std::stringstream message;
      message << methodName << " number of angle atom indices=" << angleIndices.size() << " does not equal number of angle parameter entries=" << angleParameters.size();
      throw OpenMMException( message.str() );
   } else if( getLog() ){
      (void) fprintf( getLog(), "%s angle bonds=%d\n", methodName.c_str(), angleIndices.size() );
      (void) fflush( getLog() );
   }

   if( periodicTorsionIndices.size() != periodicTorsionParameters.size() ){
      std::stringstream message;
      message << methodName << " number of periodicTorsion atom indices=" << periodicTorsionIndices.size() << " does not equal number of periodicTorsion parameter entries=" << periodicTorsionParameters.size();
      throw OpenMMException( message.str() );
   } else if( getLog() ){
      (void) fprintf( getLog(), "%s periodicTorsion bonds=%d\n", methodName.c_str(), periodicTorsionIndices.size() );
      (void) fflush( getLog() );
   }

   if( rbTorsionIndices.size() != rbTorsionParameters.size() ){
      std::stringstream message;
      message << methodName << " number of rbTorsion atom indices=" << rbTorsionIndices.size() << " does not equal number of rbTorsion parameter entries=" << rbTorsionParameters.size();
      throw OpenMMException( message.str() );
   } else if( getLog() ){
      (void) fprintf( getLog(), "%s rbTorsion bonds=%d\n", methodName.c_str(), rbTorsionIndices.size() );
      (void) fflush( getLog() );
   }

   if( (numberOfAtoms != (int) nonbondedParameters.size()) && bonded14Indices.size() > 0 ){
      std::stringstream message;
      message << methodName << " number atoms=" << numberOfAtoms << " does not equal number of nb parameter entries=" << nonbondedParameters.size();
      throw OpenMMException( message.str() );
   } else if( getLog() ){
      (void) fprintf( getLog(), "%s LJ 14 ixns=%d\n", methodName.c_str(), bonded14Indices.size() );
      (void) fflush( getLog() );
   }

   // allocate temp memory

   int maxBonds = 10*numberOfAtoms;
   int* atoms   = new int[5*maxBonds];

   BrookOpenMMFloat** params = new BrookOpenMMFloat*[getNumberOfParameterStreams()];
   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      params[ii] = new BrookOpenMMFloat[4*maxBonds];
   }
   BrookOpenMMFloat* charges = new BrookOpenMMFloat[numberOfAtoms];

   // Initialize all atom indices to -1 to indicate empty slots
   // All parameters must be initialized to values that will 
   // produce zero for the corresponding force. 

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
   addRBDihedrals( &nbondeds, atoms, params, rbTorsionIndices,       rbTorsionParameters       );
   addPDihedrals ( &nbondeds, atoms, params, periodicTorsionIndices, periodicTorsionParameters );
   addAngles(      &nbondeds, atoms, params, angleIndices,           angleParameters           );
   addBonds(       &nbondeds, atoms, params, bondIndices,            bondParameters            );

   addPairs( &nbondeds, atoms, params, charges, bonded14Indices, nonbondedParameters, lj14Scale, coulombScale );

   // check that number of bonds not too large for memory allocated

   if( nbondeds >= maxBonds ){
      std::stringstream message;
      message << methodName << " number of bonds=" << nbondeds << " is greater than maxBonds=" << maxBonds;
      throw OpenMMException( message.str() );
   } else if( getLog() ){
      (void) fprintf( getLog(), "%s atoms=%d number of bonds=%d maxBonds=%d\n", methodName.c_str(), numberOfAtoms, nbondeds, maxBonds );
      (void) fflush( getLog() );
   }

   // get factory

   const BrookStreamFactory& brookStreamFactory = dynamic_cast<const BrookStreamFactory&> (platform.getDefaultStreamFactory());

   // build streams

   StreamImpl* atomIndicesStream    = brookStreamFactory.createStreamImpl( BrookStreamFactory::BondedAtomIndicesStream, nbondeds, Stream::Float4, platform );
   _atomIndicesStream               = dynamic_cast<BrookFloatStreamImpl*> (atomIndicesStream); 
   _atomIndicesStream->loadFromArray( atoms, Stream::Integer ); 

   StreamImpl* chargeStream         = brookStreamFactory.createStreamImpl( BrookStreamFactory::BondedChargeStream, numberOfAtoms, Stream::Float, platform );
   _chargeStream                    = dynamic_cast<BrookFloatStreamImpl*> (chargeStream); 
   _chargeStream->loadFromArray( charges ); 

   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      StreamImpl* bondedParameters  = brookStreamFactory.createStreamImpl( BrookStreamFactory::BondedParametersStream, nbondeds, Stream::Float4, platform );
      _bondedParameters[ii]         = dynamic_cast<BrookFloatStreamImpl*> (bondedParameters); 
      _bondedParameters[ii]->loadFromArray( params[ii] );
   }

   // debug stuff

   if( 1 && getLog() ){

      BrookFloatStreamImpl* atomIndicesStream = dynamic_cast<BrookFloatStreamImpl*> (_atomIndicesStream); 
      (void) fprintf( getLog(), "%s nbondeds=%d strDim [%d %d ] sz=%d\n", methodName.c_str(), nbondeds,
                      atomIndicesStream->getStreamWidth(),
                      atomIndicesStream->getStreamHeight(), 
                      atomIndicesStream->getStreamSize() );

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
                              "     14[%15.6e %15.6e]\n", index,
                         ATOMS(index, 0), ATOMS(index, 1), ATOMS(index, 2), ATOMS(index, 3),
                         params[kIndex][ii], params[kIndex][ii+1], params[kIndex][ii+2], params[kIndex][ii+3], params[jIndex][ii],
                         params[jIndex][ii+1], params[jIndex][ii+2], params[jIndex][ii+3],
                         params[iIndex][ii], params[iIndex][ii+1], params[iIndex][ii+2], params[iIndex][ii+3], 
                         params[lIndex][ii], params[lIndex][ii+1], params[lIndex][ii+2], params[lIndex][ii+3],
                         params[mIndex][ii], params[mIndex][ii+1], params[mIndex][ii+2], params[mIndex][ii+3] );
      }
   }

   // load inverse maps to streams

   loadInvMaps( nbondeds, getNumberOfAtoms(), atoms, platform );
   
   // free memory

   for( int ii = 0; ii < getNumberOfParameterStreams(); ii++ ){
      delete[] params[ii];
   }
   delete[] params;
   delete[] atoms;
   delete[] charges;

   // set the fudge factors

   _ljScale        = (BrookOpenMMFloat) lj14Scale;
   _coulombFactor  = (BrookOpenMMFloat) coulombScale;

   // initialize output streams

   for( int ii = 0; ii < getNumberOfForceStreams(); ii++ ){
      StreamImpl* bondedForceStreams = brookStreamFactory.createStreamImpl( BrookStreamFactory::UnrolledForceStream, nbondeds, Stream::Float3, platform );
      _bondedForceStreams[ii]        = dynamic_cast<BrookFloatStreamImpl*> (bondedForceStreams); 
   }

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

std::string BrookBonded::getContents( int level ) const {

// ---------------------------------------------------------------------------------------

   static const std::string methodName      = "BrookBonded::getContents";

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

   (void) LOCAL_SPRINTF( value, "%d", getNumberOfAtoms() );
   message << _getLine( tab, "Number of atoms:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getLJ_14Scale() );
   message << _getLine( tab, "LJ 14 scaling:", value ); 

   (void) LOCAL_SPRINTF( value, "%.5f", getCoulombFactor() );
   message << _getLine( tab, "Coulomb factor:", value ); 

/*
   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamWidth() );
   message << _getLine( tab, "Atom stream width:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamHeight() );
   message << _getLine( tab, "Atom stream height:", value ); 

   (void) LOCAL_SPRINTF( value, "%d", getAtomStreamSize() );
   message << _getLine( tab, "Atom stream size:", value ); 
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
   message << _getLine( tab, "Atom indices stream:", (getAtomIndicesStream() ? Set : NotSet) ); 
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

