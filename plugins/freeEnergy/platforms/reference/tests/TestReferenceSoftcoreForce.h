#ifndef TEST_CUDA_SOFTCORE_H_
#define TEST_CUDA_SOFTCORE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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

/**
 * Utility methods shared across unit tests
 */

#include "../../../tests/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "OpenMM.h"

#include "OpenMMFreeEnergy.h"
#include "openmm/freeEnergyKernels.h"

#include "sfmt/SFMT.h"
#include "openmm/VerletIntegrator.h"

#ifdef OPENMM_SERIALIZE
#include "../../../../../serialization/include/openmm/serialization/SerializationProxy.h"
#include "../../../../../serialization/include/openmm/serialization/SerializationNode.h"
#include "../../../../../serialization/include/openmm/serialization/XmlSerializer.h"
//extern "C" void registerSerializationProxies();
//extern "C" void registerAmoebaSerializationProxies();
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <typeinfo>

extern "C" void registerFreeEnergyCudaKernelFactories();

using namespace OpenMM;
using namespace std;

const double TOL = 1e-5;

static const int NoCutoff          = 0;
static const int CutoffNonPeriodic = 1;
static const int CutoffPeriodic    = 2;

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;
typedef std::vector<StringVector> StringVectorVector;

typedef std::vector<int> IntVector;
typedef IntVector::iterator IntVectorI;
typedef IntVector::const_iterator IntVectorCI;
typedef std::vector<IntVector> IntVectorVector;

typedef std::vector<double> DoubleVector;
typedef DoubleVector::iterator DoubleVectorI;
typedef DoubleVector::const_iterator DoubleVectorCI;
typedef std::vector<DoubleVector> DoubleVectorVector;

// the following are used in parsing parameter file

typedef std::vector<std::string> StringVector;
typedef StringVector::iterator StringVectorI;
typedef StringVector::const_iterator StringVectorCI;

typedef std::vector<StringVector> VectorStringVector;
typedef VectorStringVector::iterator VectorStringVectorI;
typedef VectorStringVector::const_iterator VectorStringVectorCI;

typedef std::vector<std::vector<double> > VectorOfVectors;
typedef VectorOfVectors::iterator VectorOfVectorsI;
typedef VectorOfVectors::const_iterator VectorOfVectorsCI;

typedef std::map< int, int> MapIntInt;
typedef MapIntInt::iterator MapIntIntI;
typedef MapIntInt::const_iterator MapIntIntCI;

typedef std::map< double, int> MapDoubleToInt;
typedef MapDoubleToInt::iterator MapDoubleToIntI;
typedef MapDoubleToInt::const_iterator MapDoubleToIntCI;

typedef std::map< std::string, VectorOfVectors > MapStringVectorOfVectors;
typedef MapStringVectorOfVectors::iterator MapStringVectorOfVectorsI;
typedef MapStringVectorOfVectors::const_iterator MapStringVectorOfVectorsCI;

typedef std::map< std::string, StringVector > MapStringStringVector;
typedef MapStringStringVector::iterator MapStringStringVectorI;
typedef MapStringStringVector::const_iterator MapStringStringVectorCI;

typedef std::map< std::string, std::string > MapStringString;
typedef MapStringString::iterator MapStringStringI;
typedef MapStringString::const_iterator MapStringStringCI;

typedef std::map< std::string, std::string > MapStringToInt;
typedef MapStringToInt::iterator MapStringToIntI;
typedef MapStringToInt::const_iterator MapStringToIntCI;

typedef std::vector< std::map< std::string, std::string > > VectorMapStringString;
typedef VectorMapStringString::iterator VectorMapStringStringI;
typedef VectorMapStringString::const_iterator VectorMapStringStringCI;

typedef std::map< std::string, int > MapStringInt;
typedef MapStringInt::iterator MapStringIntI;
typedef MapStringInt::const_iterator MapStringIntCI;

typedef std::map< std::string,  std::vector<OpenMM::Vec3> > MapStringVec3;
typedef MapStringVec3::iterator MapStringVec3I;
typedef MapStringVec3::const_iterator MapStringVec3CI;

typedef std::map< std::string, double > MapStringToDouble;
typedef MapStringToDouble::iterator MapStringToDoubleI;
typedef MapStringToDouble::const_iterator MapStringToDoubleCI;
typedef std::vector< MapStringToDouble > VectorOfMapStringToDouble;

typedef std::map< std::string, DoubleVector> MapStringToDoubleVector;
typedef MapStringToDoubleVector::iterator MapStringToDoubleVectorI;
typedef MapStringToDoubleVector::const_iterator MapStringToDoubleVectorCI;

typedef std::map< std::string, DoubleVector  > MapStringToDoubleVector;

typedef std::map< std::string, Force*> MapStringForce;
typedef MapStringForce::iterator MapStringForceI;
typedef MapStringForce::const_iterator MapStringForceCI;

typedef std::map< int, IntVector> MapIntIntVector;
typedef MapIntIntVector::const_iterator MapIntIntVectorCI;

typedef std::pair<int, int> IntIntPair;
typedef std::vector<IntIntPair> IntIntPairVector;
typedef IntIntPairVector::iterator IntIntPairVectorI;
typedef IntIntPairVector::const_iterator IntIntPairVectorCI;

typedef std::pair<int, double> IntDoublePair;
typedef std::vector<IntDoublePair> IntDoublePairVector;
typedef IntDoublePairVector::iterator IntDoublePairVectorI;
typedef IntDoublePairVector::const_iterator IntDoublePairVectorCI;

class PositionGenerator {

public:

    enum GenerationMethod {
        /** 
         * Random positions
         */
        Random = 0,
        /** 
         * On grid
         */
        SimpleGrid = 1,
    };  

    PositionGenerator( int numMolecules, int numParticlesPerMolecule, double boxSize );

    ~PositionGenerator( );

    /** 
     * Set logging file reference
     *
     * @param log                       log
     *
     */
     
    void setLog( FILE* log );
    
    /** 
     * Set positions
     *
     * @param method                    method placement
     * @param positions                 output vector of positions
     *
     * @return nonzero if error detected; 0 otherwise
     */
     
    int setPositions( GenerationMethod method, std::vector<Vec3>& positions ) const;
    
    /** 
     * Set positions
     *
     * @param method                    method placement
     * @param sfmt                      input random number generator
     * @param positions                 output vector of positions
     *
     * @return nonzero if error detected; 0 otherwise
     */
     
    int setPositions( GenerationMethod method, OpenMM_SFMT::SFMT& sfmt, std::vector<Vec3>& positions ) const;
    
    /** 
     * Place particles on a grid
     *
     * @param origin                    origin
     * @param boxDimensions             box dimensions
     * @param spacing                   spacing
     * @param sfmt                      input random number generator
     * @param array                     output vector of grid values
     *
     * @return -1 if particles will not fit on grid; 0 if they do
     */
     
    int setParticlesOnGrid( const Vec3& origin, const Vec3& boxDimensions, const Vec3& spacing, 
                            OpenMM_SFMT::SFMT& sfmt, std::vector<Vec3>& array ) const;
    
    /** 
     * Set bond distance
     *
     * @param bondDistance bond distance
     */
     
    void setBondDistance( double bondDistance );
    
    /** 
     * Get bond distance
     *
     * @return bond distance
     */
     
    double getBondDistance( void ) const;

    /** 
     * Get distance
     *
     * @param index1 index of first particle
     * @param index2 index of second particle
     * @param positions particle positions
     *
     * @return distance
     */
     
    double getDistance( int index1, int index2, const std::vector<Vec3>& positions ) const;

    /** 
     * Get distance assumming periodic boundary conditions
     *
     * @param index1 index of first particle
     * @param index2 index of second particle
     * @param positions particle positions
     *
     * @return distance
     */
     
    double getPeriodicDistance( int index1, int index2, const std::vector<Vec3>& positions ) const;

    /** 
     * Get settings
     *
     * @return info string
     */
     
    std::string getSettings( void ) const;

    /** 
     * Get enclosing box
     *
     * @param positions    input vector of positions
     * @param enclosingBox output vector of enclosing box dimensions 
     *
     */
     
    void getEnclosingBox( const std::vector<Vec3>& positions, Vec3 enclosingBox[2] ) const;

    /** 
     * Get sorted distances from particular position
     *
     * @param periodicBoundaryConditions if set, apply PBC
     * @param positionIndex              input position index
     * @param positions                  input vector of positions
     * @param sortVector                 on output sorted IntDoublePairVector
     *
     */
     
    void getSortedDistances( int periodicBoundaryConditions, int positionIndex, const std::vector<Vec3>& positions, IntDoublePairVector& sortVector ) const;

private:

    int _numMolecules;
    int _numParticlesPerMolecule;
    int _numParticles;
    int _seed;

    double _boxSize;
    double _bondDistance;
    Vec3 _origin;
    Vec3 _boxDimensions;
    Vec3 _spacing;

    FILE* _log;
};

PositionGenerator::PositionGenerator( int numMolecules, int numParticlesPerMolecule, double boxSize ) :
               _numMolecules(numMolecules), 
               _seed(0),
               _log(NULL),
               _bondDistance(0.1),
               _numParticlesPerMolecule(numParticlesPerMolecule),
               _numParticles(numMolecules*numParticlesPerMolecule),
               _boxSize(boxSize),
               _boxDimensions(Vec3(boxSize,boxSize,boxSize)),
               _origin(Vec3(0.0,0.0,0.0)) {

    double particlesPerDimension  = pow( static_cast<double>(_numParticles), (1.0/3.0) ); 
    int particlesPerDimensionI    = static_cast<int>(particlesPerDimension+0.999999); 
    double spacingPerDimension    = _boxSize/particlesPerDimension;

    _spacing[0]                   = spacingPerDimension;
    _spacing[1]                   = spacingPerDimension;
    _spacing[2]                   = spacingPerDimension;

}

PositionGenerator::~PositionGenerator( ){};

void PositionGenerator::setBondDistance( double bondDistance ){
    _bondDistance = bondDistance;
}

void PositionGenerator::setLog( FILE* log ){
    _log = log;
}

double PositionGenerator::getBondDistance( void ) const {
    return _bondDistance;
}

double PositionGenerator::getDistance( int index1, int index2, const std::vector<Vec3>& positions ) const {

    Vec3 delta = positions[index2] - positions[index1];
    return sqrt( delta.dot( delta ) );
}

double PositionGenerator::getPeriodicDistance( int index1, int index2, const std::vector<Vec3>& positions ) const {

    Vec3 delta  = positions[index2] - positions[index1];
    if( _boxSize > 0.0 ){
        delta[0]   -= floor(delta[0]/_boxSize+0.5f)*_boxSize;
        delta[1]   -= floor(delta[1]/_boxSize+0.5f)*_boxSize;
        delta[2]   -= floor(delta[2]/_boxSize+0.5f)*_boxSize;
    }
    return sqrt( delta.dot( delta ) );
}

/** 
 * Get positions
 *
 * @param method                    method placement
 * @param positions                 output vector of positions
 *
 * @return nonzero if error detected; 0 otherwise
 */
 
int PositionGenerator::setPositions( GenerationMethod method, std::vector<Vec3>& positions ) const {

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand( _seed, sfmt);
    return setPositions( method, sfmt, positions );
}

/** 
 * Get settings
 *
 * @return info string
 */
 
std::string PositionGenerator::getSettings( void ) const {

    std::stringstream msg;
    msg << "numMolecules            " << _numMolecules            << std::endl;
    msg << "numParticlesPerMolecule " << _numParticlesPerMolecule << std::endl;
    msg << "boxSize                 " << _boxSize                 << std::endl;
    msg << "spacing                 " << _spacing[0]              << std::endl;
    msg << "seed                    " << _seed                    << std::endl;

    return msg.str();
}

/** 
 * Get positions
 *
 * @param method                    method placement
 * @param sfmt                      input random number generator
 * @param positions                 output vector of positions
 *
 * @return nonzero if error detected; 0 otherwise
 */
 
int PositionGenerator::setPositions( GenerationMethod method, OpenMM_SFMT::SFMT& sfmt, std::vector<Vec3>& positions ) const {

    int errorFlag = 0;
    positions.resize( _numParticles );
    if( method == Random ){
        for( unsigned int ii = 0; ii <  static_cast<unsigned int>(_numParticles); ii += _numParticlesPerMolecule ){ 
            positions[ii]    = Vec3(_boxSize*genrand_real2(sfmt), _boxSize*genrand_real2(sfmt), _boxSize*genrand_real2(sfmt));
            for( unsigned int jj = 1; jj < static_cast<unsigned int>(_numParticlesPerMolecule); jj++) { 
                positions[ii+jj]  = positions[ii] + Vec3(_bondDistance*genrand_real2(sfmt), _bondDistance*genrand_real2(sfmt), _bondDistance*genrand_real2(sfmt));
            }
        }
    } else if( method == SimpleGrid ){

        Vec3 origin, boxDimensions, spacing;
        std::stringstream msg;
        if( _numParticlesPerMolecule > 1 && _bondDistance > 0.0 ){
            origin                        = Vec3(_bondDistance, _bondDistance, _bondDistance );
            double particlesPerDimension  = pow( static_cast<double>(_numParticles), (1.0/3.0) ); 
            int particlesPerDimensionI    = static_cast<int>(particlesPerDimension+0.999999); 
            double boxSize                = (_boxSize-2.0*_bondDistance);
            double spacingPerDimension    = boxSize/particlesPerDimension;
            spacing                       = Vec3(spacingPerDimension, spacingPerDimension, spacingPerDimension );
            boxDimensions                 = Vec3(boxSize, boxSize, boxSize );

            msg << "Bond distance " << _bondDistance << std::endl;
            msg << "particlesPerDimension " << particlesPerDimension << std::endl;
            msg << "boxSize " << boxSize << std::endl;
            msg << "spacingPerDimension " << spacingPerDimension << std::endl;

        } else {
           origin                         = _origin;
           spacing                        = _spacing;
           boxDimensions                  = _boxDimensions;

        }
        msg << getSettings() << std::endl;

        if( _log ){
            (void) fprintf( _log, "SimpleGrid %s\n", msg.str().c_str() );
        }

        errorFlag = setParticlesOnGrid( origin, boxDimensions, spacing, sfmt, positions );
    }

    return errorFlag;
}

/** 
 * Place particles on a grid
 *
 * @param origin                    origin
 * @param boxDimensions             box dimensions
 * @param spacing                   spacing
 * @param array                     output vector of grid values
 *
 * @return -1 if particles will not fit on grid; 0 if they do
 */
 
int PositionGenerator::setParticlesOnGrid( const Vec3& origin, const Vec3& boxDimensions, const Vec3& spacing, OpenMM_SFMT::SFMT& sfmt,
                                           std::vector<Vec3>& array ) const {

    Vec3 start(origin);

    if( array.size() != _numParticles ){
        std::stringstream msg;
        msg << "PositionGenerator::setParticlesOnGrid position vector size=" << array.size() << " != numParticles=" << _numParticles;
        msg << getSettings();
        throw OpenMMException( msg.str() );
    }

    // place molecule centers on grid

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(_numParticles); ii += _numParticlesPerMolecule ){
        array[ii]  = Vec3(start);
        bool done  = false;
        for( unsigned int jj = 0; jj < 3 && !done; jj++ ){
            start[jj]  += spacing[jj];
            if( start[jj] > boxDimensions[jj] ){
                start[jj] = origin[jj];
            } else {
                done = true;
            }
        }
        if( !done ){
            std::stringstream msg;
            msg << "PositionGenerator::setParticlesOnGrid error in grid settings";
            throw OpenMMException( msg.str() );
        }
    }

    // add molecule atoms

    Vec3 bondOffset( 0.05, 0.05, 0.05 );
    for( unsigned int ii = 0; ii < static_cast<unsigned int>(_numMolecules); ii++ ){
        int molecularIndex = ii*_numParticlesPerMolecule;
        for( unsigned int jj = 1; jj < static_cast<unsigned int>(_numParticlesPerMolecule); jj++ ){
            array[molecularIndex+jj] = array[molecularIndex] + bondOffset + Vec3(_bondDistance*genrand_real2(sfmt), _bondDistance*genrand_real2(sfmt), _bondDistance*genrand_real2(sfmt));
        }
    }

    return 0;
}

/** 
 * Get enclosing box
 *
 * @param positions    input vector of positions
 * @param enclosingBox output Vec3[2] of minimum enclosing box ranges
 *
 */
 
void PositionGenerator::getEnclosingBox( const std::vector<Vec3>& positions, Vec3 enclosingBox[2] ) const {

    enclosingBox[0][0] = enclosingBox[1][0] = positions[0][0];
    enclosingBox[0][1] = enclosingBox[1][1] = positions[0][1];
    enclosingBox[0][2] = enclosingBox[1][2] = positions[0][2];
 
    for( unsigned int ii = 1; ii < positions.size(); ii++ ){
        if( enclosingBox[0][0] > positions[ii][0] ){
            enclosingBox[0][0] = positions[ii][0];
        }    
        if( enclosingBox[1][0] < positions[ii][0] ){
            enclosingBox[1][0] = positions[ii][0];
        }    
        if( enclosingBox[0][1] > positions[ii][1] ){
            enclosingBox[0][1] = positions[ii][1];
        }    
        if( enclosingBox[1][1] < positions[ii][1] ){
            enclosingBox[1][1] = positions[ii][1];
        }    
        if( enclosingBox[0][2] > positions[ii][2] ){
            enclosingBox[0][2] = positions[ii][2];
        }    
        if( enclosingBox[1][2] < positions[ii][2] ){
            enclosingBox[1][2] = positions[ii][2];
        }    
    }    
 
    return;
}


/** 
 * Predicate for sorting <int,double> pair
 *
 * @param d1 first  IntDoublePair to compare
 * @param d2 second IntDoublePair to compare
 *
 */
 
bool TestIntDoublePair( const IntDoublePair& d1, const IntDoublePair& d2 ){
   return d1.second < d2.second;
}

/** 
 * Get sorted distances from particular position
 *
 * @param periodicBoundaryConditions if set, apply PBC
 * @param positionIndex              input position index
 * @param positions                  input vector of positions
 * @param sortVector                 on output sorted IntDoublePairVector
 *
 */
 
void PositionGenerator::getSortedDistances( int periodicBoundaryConditions, int positionIndex, const std::vector<Vec3>& positions,
                                            IntDoublePairVector& sortVector ) const {

    sortVector.resize( 0 );
    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        if( ii == positionIndex )continue;
        double distance = periodicBoundaryConditions ? getPeriodicDistance( positionIndex, ii, positions) :  getDistance( positionIndex, ii, positions);
        sortVector.push_back( IntDoublePair(ii,sqrt(distance) ) );
    }    

    std::sort( sortVector.begin(), sortVector.end(), TestIntDoublePair );

    return;
}

/**---------------------------------------------------------------------------------------
 *
 * Set string field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setStringFromMap( MapStringString& argumentMap, std::string fieldToCheck, std::string& fieldToSet ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName             = "setStringFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = (*check).second; 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Set int field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setIntFromMap( MapStringString& argumentMap, std::string fieldToCheck, int& fieldToSet ){

// ---------------------------------------------------------------------------------------

   //static const std::string methodName             = "setIntFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = atoi( (*check).second.c_str() ); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Set int field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setIntFromMapStringToDouble( MapStringToDouble& argumentMap, std::string fieldToCheck, int& fieldToSet ){

// ---------------------------------------------------------------------------------------

   MapStringToDoubleCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet   = static_cast<int>(check->second+0.0000001);
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------

 * Set float field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setFloatFromMap( MapStringString& argumentMap, std::string fieldToCheck, float& fieldToSet ){

// ---------------------------------------------------------------------------------------

   static const std::string methodName             = "setFloatFromMap";

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = static_cast<float>(atof( (*check).second.c_str() )); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Set double field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setDoubleFromMap( MapStringString& argumentMap, std::string fieldToCheck, double& fieldToSet ){

// ---------------------------------------------------------------------------------------

   MapStringStringCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = atof( (*check).second.c_str() ); 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Set double field if in map
 * 
 * @param  argumentMap            map to check
 * @param  fieldToCheck           key
 * @param  fieldToSet             field to set
 *
 * @return 1 if argument set, else 0
 *
   --------------------------------------------------------------------------------------- */

static int setDoubleFromMapStringToDouble( MapStringToDouble& argumentMap, std::string fieldToCheck, double& fieldToSet ){

// ---------------------------------------------------------------------------------------

   MapStringToDoubleCI check = argumentMap.find( fieldToCheck );
   if( check != argumentMap.end() ){
      fieldToSet = check->second; 
      return 1;
   }
   return 0;
}

/**---------------------------------------------------------------------------------------
 *
 * Compare forces from two states
 * 
 * @param  state1               state1
 * @param  state2               state2
 * @param  relativeTolerance    relative tolerance
 * @param  log                  if set, output forces
 *
 * @return number of entries with relative difference > tolerance 
 *
   --------------------------------------------------------------------------------------- */

int compareForcesOfTwoStates( State& state1, State& state2, double relativeTolerance,
                              DoubleVector& stats, FILE* log ) {

    int error                             = 0;
    vector<Vec3> force1                   = state1.getForces();
    vector<Vec3> force2                   = state2.getForces();
    double maxRelativeDifference          = -1.0e+30;
    double maxRelativeDifferenceIndex     = -1.0;
    double averageRelativeDifference      = 0.0;
    double count                          = 0.0;
    DoubleVector medians1( force1.size() );
    DoubleVector medians2( force1.size() );
    for( unsigned int ii = 0; ii < force1.size(); ii++ ){

        Vec3 f1                = force1[ii];
        Vec3 f2                = force2[ii];

        double diff            = (f1[0] - f2[0])*(f1[0] - f2[0]) +
                                 (f1[1] - f2[1])*(f1[1] - f2[1]) +
                                 (f1[2] - f2[2])*(f1[2] - f2[2]); 

        double denom1          = sqrt( f1[0]*f1[0] + f1[1]*f1[1] + f1[2]*f1[2] );
        double denom2          = sqrt( f2[0]*f2[0] + f2[1]*f2[1] + f2[2]*f2[2] );
        medians1[ii]            = denom1;
        medians2[ii]            = denom2;
        double relativeDiff;
        if( denom1 > 0.0 || denom2 > 0.0 ){
            relativeDiff = 2.0*sqrt( diff )/(denom1+denom2);
        } else {
            relativeDiff = 0.0;
        }

        if( relativeDiff > maxRelativeDifference ){
            maxRelativeDifference      = relativeDiff;
            maxRelativeDifferenceIndex = static_cast<double>(ii);
        }
        averageRelativeDifference += relativeDiff;
        count                     += 1.0;

        if( relativeDiff > relativeTolerance ){
           error++;
        }
        if( log ){
            (void) fprintf( log, "F %6u %15.7e [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] %15.7e %15.7e %s\n", static_cast<unsigned int>(ii), 
                            relativeDiff, f1[0], f1[1], f1[2], f2[0], f2[1], f2[2], denom1, denom2, (relativeDiff < relativeTolerance ? "":"XXXXXX") );
        }
    }

    if( count > 0.0 ){
        averageRelativeDifference /= count;
    }

    std::sort( medians1.begin(), medians1.end() );
    std::sort( medians2.begin(), medians2.end() );
    double median1 = medians1[medians1.size()/2];
    double median2 = medians2[medians2.size()/2];

    stats.resize( 4 );
    stats[0] = averageRelativeDifference;
    stats[1] = maxRelativeDifference;
    stats[2] = maxRelativeDifferenceIndex;
    stats[3] = median1 < median2 ? median1 : median2;
    
    return error;
}


/** 
 * Get forces in system
 *
 * @param system                   system to serialize
 * @param stringForceVector        output stringForceVector[forceName] = force index
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static void getStringForceMap( System& system, MapStringInt& stringForceVector, FILE* log ){

    // print active forces and relevant parameters

    for( int ii = 0; ii < system.getNumForces(); ii++ ) {

        int hit                 = 0;
        Force& force            = system.getForce(ii);
        if( !hit ){

            try {
               CMAPTorsionForce& castForce = dynamic_cast<CMAPTorsionForce&>(force);
               stringForceVector["CMAPTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               CustomAngleForce& castForce = dynamic_cast<CustomAngleForce&>(force);
               stringForceVector["CustomAngle"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }


        if( !hit ){

            try {
               CustomBondForce& castForce = dynamic_cast<CustomBondForce&>(force);
               stringForceVector["CustomBond"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               CustomExternalForce& castForce = dynamic_cast<CustomExternalForce&>(force);
               stringForceVector["CustomExternal"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               CustomGBForce& castForce = dynamic_cast<CustomGBForce&>(force);
               stringForceVector["CustomGB"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               CustomHbondForce& castForce = dynamic_cast<CustomHbondForce&>(force);
               stringForceVector["CustomHbond"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               CustomNonbondedForce& castForce = dynamic_cast<CustomNonbondedForce&>(force);
               stringForceVector["CustomNonbonded"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }


        if( !hit ){

            try {
               CustomTorsionForce& castForce = dynamic_cast<CustomTorsionForce&>(force);
               stringForceVector["CustomTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }


        if( !hit ){

            try {
               GBSAOBCForce& castForce = dynamic_cast<GBSAOBCForce&>(force);
               stringForceVector["GBSAOBC"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               GBVIForce& castForce = dynamic_cast<GBVIForce&>(force);
               stringForceVector["GBVI"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               HarmonicAngleForce& castForce = dynamic_cast<HarmonicAngleForce&>(force);
               stringForceVector["HarmonicAngle"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }


        if( !hit ){

            try {
               HarmonicBondForce& castForce = dynamic_cast<HarmonicBondForce&>(force);
               stringForceVector["HarmonicBond"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               NonbondedForce& castForce = dynamic_cast<NonbondedForce&>(force);
               stringForceVector["Nonbonded"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               PeriodicTorsionForce& castForce = dynamic_cast<PeriodicTorsionForce&>(force);
               stringForceVector["PeriodicTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               RBTorsionForce& castForce = dynamic_cast<RBTorsionForce&>(force);
               stringForceVector["RBTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               MonteCarloBarostat& castForce = dynamic_cast<MonteCarloBarostat&>(force);
               stringForceVector["MonteCarloBarostat"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AndersenThermostat& castForce = dynamic_cast<AndersenThermostat&>(force);
               stringForceVector["AndersenThermostat"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

#ifdef USE_SOFTCORE
        if( !hit ){

            try {
               GBSAOBCSoftcoreForce& castForce = dynamic_cast<GBSAOBCSoftcoreForce&>(force);
               stringForceVector["GBSAOBCSoftcore"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               GBVISoftcoreForce& castForce = dynamic_cast<GBVISoftcoreForce&>(force);
               stringForceVector["GBVISoftcore"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               NonbondedSoftcoreForce& castForce = dynamic_cast<NonbondedSoftcoreForce&>(force);
               stringForceVector["NonbondedSoftcore"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }
#endif

#ifdef INCLUDE_AMOEBA_FORCES

        if( !hit ){

            try {
               AmoebaHarmonicBondForce& castForce = dynamic_cast<AmoebaHarmonicBondForce&>(force);
               stringForceVector["AmoebaHarmonicBond"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaHarmonicAngleForce& castForce = dynamic_cast<AmoebaHarmonicAngleForce&>(force);
               stringForceVector["AmoebaHarmonicAngle"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaHarmonicInPlaneAngleForce& castForce = dynamic_cast<AmoebaHarmonicInPlaneAngleForce&>(force);
               stringForceVector["AmoebaHarmonicInPlaneAngle"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaMultipoleForce& castForce = dynamic_cast<AmoebaMultipoleForce&>(force);
               stringForceVector["AmoebaMultipole"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaOutOfPlaneBendForce& castForce = dynamic_cast<AmoebaOutOfPlaneBendForce&>(force);
               stringForceVector["AmoebaOutOfPlaneBend"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaPiTorsionForce& castForce = dynamic_cast<AmoebaPiTorsionForce&>(force);
               stringForceVector["AmoebaPiTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaStretchBendForce& castForce = dynamic_cast<AmoebaStretchBendForce&>(force);
               stringForceVector["AmoebaStretchBend"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaTorsionForce& castForce = dynamic_cast<AmoebaTorsionForce&>(force);
               stringForceVector["AmoebaTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaTorsionTorsionForce& castForce = dynamic_cast<AmoebaTorsionTorsionForce&>(force);
               stringForceVector["AmoebaTorsionTorsion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaUreyBradleyForce& castForce = dynamic_cast<AmoebaUreyBradleyForce&>(force);
               stringForceVector["AmoebaUreyBradley"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaVdwForce& castForce = dynamic_cast<AmoebaVdwForce&>(force);
               stringForceVector["AmoebaVdw"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaWcaDispersionForce& castForce = dynamic_cast<AmoebaWcaDispersionForce&>(force);
               stringForceVector["AmoebaWcaDispersion"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaGeneralizedKirkwoodForce& castForce = dynamic_cast<AmoebaGeneralizedKirkwoodForce&>(force);
               stringForceVector["AmoebaGeneralizedKirkwood"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

        if( !hit ){

            try {
               AmoebaTorsionTorsionForce& castForce = dynamic_cast<AmoebaTorsionTorsionForce&>(force);
               stringForceVector["AmoebaTorsionTorsionForce"] = ii;
               hit++;
            } catch( std::bad_cast ){
            }    
        }

#endif

        // COM

        if( !hit ){
    
            try {
               CMMotionRemover& cMMotionRemover = dynamic_cast<CMMotionRemover&>(force);
               hit++;
            } catch( std::bad_cast ){
            }
        }

        if( !hit && log ){
           (void) fprintf( log, "   entry=%2d force not recognized.\n", ii );
        }

    }
}

/**---------------------------------------------------------------------------------------
 *
 * Copy NonbondedSoftcoreForce
 * 
 * @param  nonbondedSoftcoreForce  NonbondedSoftcoreForce to copy
 *
 * @return copy of nonbondedSoftcoreForce
 *
   --------------------------------------------------------------------------------------- */

static NonbondedSoftcoreForce* copyNonbondedSoftcoreForce( const NonbondedSoftcoreForce& nonbondedSoftcoreForce ){

    NonbondedSoftcoreForce* copyNonbondedSoftcoreForce = new NonbondedSoftcoreForce( nonbondedSoftcoreForce );

/*
    copyNonbondedSoftcoreForce->setNonbondedMethod( nonbondedSoftcoreForce.getNonbondedMethod() );
    copyNonbondedSoftcoreForce->setCutoffDistance( nonbondedSoftcoreForce.getCutoffDistance() );
    copyNonbondedSoftcoreForce->setReactionFieldDielectric( nonbondedSoftcoreForce.getReactionFieldDielectric() );

    // particle parameters

    for( unsigned int ii = 0; ii < nonbondedSoftcoreForce.getNumParticles(); ii++ ){

        double charge;
        double sigma;
        double epsilon;
        double softcoreLJLambda;
        nonbondedSoftcoreForce.getParticleParameters(ii, charge, sigma, epsilon, softcoreLJLambda);
        copyNonbondedSoftcoreForce->addParticle( charge, sigma, epsilon, softcoreLJLambda);
    }

    // exceptions

    for( unsigned int ii = 0; ii < nonbondedSoftcoreForce.getNumExceptions(); ii++ ){

        int particle1, particle2;
        double chargeProd;
        double sigma;
        double epsilon;
        double softcoreLJLambda;
        nonbondedSoftcoreForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
        copyNonbondedSoftcoreForce->addException( particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
    }
*/

    return copyNonbondedSoftcoreForce;
}

/**---------------------------------------------------------------------------------------
 *
 * Copy NonbondedForce
 * 
 * @param  nonbondedForce  NonbondedForce to copy
 *
 * @return copy of nonbondedForce
 *
   --------------------------------------------------------------------------------------- */

static NonbondedForce* copyNonbondedForce( const NonbondedForce& nonbondedForce ){

    NonbondedForce* copyNonbondedForce = new NonbondedForce( nonbondedForce );

/*
    copyNonbondedForce->setNonbondedMethod( nonbondedForce.getNonbondedMethod() );
    copyNonbondedForce->setCutoffDistance( nonbondedForce.getCutoffDistance() );
    copyNonbondedForce->setReactionFieldDielectric( nonbondedForce.getReactionFieldDielectric() );

    // particle parameters

    for( unsigned int ii = 0; ii < nonbondedForce.getNumParticles(); ii++ ){

        double charge;
        double sigma;
        double epsilon;
        double softcoreLJLambda;
        nonbondedForce.getParticleParameters(ii, charge, sigma, epsilon, softcoreLJLambda);
        copyNonbondedForce->addParticle( charge, sigma, epsilon, softcoreLJLambda);
    }

    // exceptions

    for( unsigned int ii = 0; ii < nonbondedForce.getNumExceptions(); ii++ ){

        int particle1, particle2;
        double chargeProd;
        double sigma;
        double epsilon;
        double softcoreLJLambda;
        nonbondedForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
        copyNonbondedForce->addException( particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );
    }
*/

    return copyNonbondedForce;
}

/** 
 * Return copy of system (but not forces)
 *
 * @param inputSystem               system to copy
 *
 * @return system copy
 *
 */
 
static void copySystem( const System& inputSystem, System& systemCopy ){

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(inputSystem.getNumParticles()); ii++ ){
        systemCopy.addParticle( inputSystem.getParticleMass( static_cast<int>(ii) ) );
    }

    Vec3 a;
    Vec3 b;
    Vec3 c;
    inputSystem.getDefaultPeriodicBoxVectors( a, b, c );
    systemCopy.setDefaultPeriodicBoxVectors( a, b, c );

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(inputSystem.getNumConstraints()); ii++ ){
        int particle1, particle2;
        double distance;
        inputSystem.getConstraintParameters( ii, particle1, particle2, distance);
        systemCopy.addConstraint( particle1, particle2, distance);
    }

    return;
}

/** 
 * Randomize parameters
 *
 * @param parametersLowerBound      vector of parameter lower bounds
 * @param parametersUpperBound      vector of parameter upper bounds
 * @param sfmt                      SFMT random number generator
 * @param parameters                output vector of randomized parameter values
 *
 */
 
static void randomizeParameters( const std::vector<double>& parametersLowerBound, 
                                 const std::vector<double>& parametersUpperBound,
                                 OpenMM_SFMT::SFMT& sfmt, std::vector<double>& parameters ){

    if( parametersLowerBound.size() != parametersUpperBound.size() ){
        std::stringstream msg;
        msg << " randomizeParameters parametersLowerBound size=" << parametersLowerBound.size() << " != parametersUpperBound size=" << parametersUpperBound.size();
        throw OpenMMException( msg.str() );
    }

    if( parametersLowerBound.size() != parameters.size() ){
        std::stringstream msg;
        msg << " randomizeParameters parametersLowerBound size=" << parametersLowerBound.size() << " != parameter size=" << parameters.size();
        throw OpenMMException( msg.str() );
    }

    for( unsigned int ii = 0; ii < parametersLowerBound.size(); ii++ ){
        parameters[ii] = parametersLowerBound[ii] + (parametersUpperBound[ii] - parametersLowerBound[ii])*(genrand_real2(sfmt));
    }

    return;
}

/** 
 * Randomize Vec3 vector
 *
 * @param average                   mean value
 * @param range                     +/- range
 * @param sfmt                      SFMT random number generator
 * @param array                     output vector of randomized values
 *
 */
 
static void randomizeVec3( double average, double range, 
                           OpenMM_SFMT::SFMT& sfmt, std::vector<Vec3>& array ){

    range *= 2.0;
    for( unsigned int ii = 0; ii < array.size(); ii++ ){
        array[ii] = Vec3( average + range*(genrand_real2(sfmt) - 0.5),
                          average + range*(genrand_real2(sfmt) - 0.5), 
                          average + range*(genrand_real2(sfmt) - 0.5) );
    }
    return;
}

/** 
 * Output contents of MapStringString 
 *
 * @param inputArgumentMap          map to output
 * @param outputStream              output stream
 *
 */
 
static void streamArgumentMap( const MapStringString& inputArgumentMap, std::stringstream& outputStream ){ 

    char buffer[2048];
    for( MapStringStringCI ii = inputArgumentMap.begin(); ii != inputArgumentMap.end(); ii++ ){
        std::string key   = ii->first;
        std::string value = ii->second;
        (void) sprintf( buffer, "      %30s %40s\n", key.c_str(), value.c_str() );
        outputStream << buffer;
    }    

    return;
}

/** 
 * Format argument/value
 *
 * @param buffer                    formatted output
 * @param key                       argument name
 * @param value                     argument value
 * @param format                    format string
 * @param call                      if call > 0, skip key name
 * @param type                      type == 0, then use int value; else double
 *
 */
 
static void formatArgument( char* buffer, const std::string& key, double value, const char* format, int call, int type ){

    // if call > 0, skip key name

    unsigned int index   = 0;
    while( index < key.size() ){
        buffer[index] = call ? ' ' : key[index];
        index++;
    }

    // add blank

    buffer[index++]        = ' ';
    buffer[index]          = static_cast<char>(NULL);

    if( type == 0 ){
        int valueInt       = static_cast<int>(value+0.00001);
        (void) sprintf( buffer + index, format, valueInt );
    } else {
        (void) sprintf( buffer + index, format, value );
    }
    return;
}

/** 
 * Output contents of MapStringString w/ all argument on one line 
 *
 * @param inputArgumentMap          map to output
 * @param exclude                   map of keys to exclude from output
 * @param outputStream              output stream
 *
 */
 
static void streamArgumentMapOneLine( const MapStringToDouble& inputArgumentMap, const MapStringToInt& exclude,
                                      const StringVector& printFirst, int callId, std::stringstream& outputStream ){ 

    char buffer[2048];

    MapStringToInt excludeAll(exclude);

    for( unsigned int ii = 0; ii < printFirst.size(); ii++ ){
        MapStringToDoubleCI iter = inputArgumentMap.find( printFirst[ii] );
        if( iter != inputArgumentMap.end() ){
            std::string key      = iter->first;
            if( exclude.find( key ) == exclude.end() ){
                double      value    = iter->second;

                if( key == "numMolecules" ){
                    formatArgument( buffer, key, value, "%6d ", callId, 0 );
                } else if( key == "nonbondedMethod" ){
                    formatArgument( buffer, key, value, "%1d ", callId, 0 );
                } else if( key == "lambda1" || key == "lambda2" ){
                    formatArgument( buffer, key, value, "%4.2f ", callId, 1 );
                } else if( key == "boxSize" ){
                    formatArgument( buffer, key, value, "%6.2f ", callId, 1 );
                } else {
                    formatArgument( buffer, key, value, "%15.7e ", callId, 1 );
                }
                outputStream << buffer;
                excludeAll[key] = 1;
            }
        }
    }    

    for( MapStringToDoubleCI ii = inputArgumentMap.begin(); ii != inputArgumentMap.end(); ii++ ){
        std::string key      = ii->first;
        if( excludeAll.find( key ) == excludeAll.end() ){
            double      value    = ii->second;
            int valueInt         = static_cast<int>(value+0.00001);
            double valueDouble   = static_cast<double>(valueInt);
            if( key == "numMolecules" ){
                (void) sprintf( buffer, "%s=%6d ", key.c_str(), valueInt );
            } else if( key == "nonbondedMethod" ){
                (void) sprintf( buffer, "%s=%1d ", key.c_str(), valueInt );
            } else if( key == "lambda1" || key == "lambda2" ){
                (void) sprintf( buffer, "%s=%4.2f ", key.c_str(), value );
            } else if( key == "boxSize" ){
                (void) sprintf( buffer, "%s=%6.2f ", key.c_str(), value );
            } else if( valueDouble == value ){
                (void) sprintf( buffer, "%s=%6d ", key.c_str(), valueInt );
            } else {
                (void) sprintf( buffer, "%s=%15.7e ", key.c_str(), value );
            }
            outputStream << buffer;
        }
    }    
    outputStream << std::endl;

    return;
}

/** 
 * Get signature of a MapStringToDouble  object
 *
 * @param inputArgumentMap          map
 * @return signature
 *
 */
 
static double getMapStringToDoubleSignature( const MapStringToDouble& inputArgumentMap ){

    double signature = 0.0;
    double offset    = 0.1;
    for( MapStringToDoubleCI ii = inputArgumentMap.begin(); ii != inputArgumentMap.end(); ii++ ){
        signature           += (offset + ii->second);
        offset              += 0.1;
    }
    return signature;
}

/** 
 * Compare two MapStringToDouble to see if they have the same (key,value) pairs
 *
 * @param inputArgumentMap1 map 1
 * @param inputArgumentMap2 map 2
 *
 * @return true if maps have  same (key,value) pairs; otherwise false
 *
 */
 
static bool compareMapStringToDoubles( const MapStringToDouble& inputArgumentMap1, const MapStringToDouble& inputArgumentMap2 ){

    if( inputArgumentMap1.size() != inputArgumentMap1.size() ){
        return false;
    }
    for( MapStringToDoubleCI ii = inputArgumentMap1.begin(); ii != inputArgumentMap1.end(); ii++ ){
        MapStringToDoubleCI jj = inputArgumentMap2.find( (*ii).first );
        if( jj == inputArgumentMap2.end() || jj->second != ii->second ){
            return false;
        }
    }
    return true;
}

/** 
 * Generate collection of inputArguments maps given
 * list of DoubleVectors for each argument
 *
 * @param inputArguments            map[argumentKey] = vector of double parameter values
 * @param argumentMaps              output vector of generated maps
 *
 */
 
static void generateInputArgumentMapsFromStringVectors( const MapStringToDoubleVector& inputArguments, 
                                                        VectorOfMapStringToDouble& argumentMaps ){

    for( MapStringToDoubleVectorCI ii = inputArguments.begin(); ii != inputArguments.end(); ii++ ){

        std::string  argumentName           = (*ii).first;
        DoubleVector arguments              = (*ii).second;
        unsigned int initialArgumentMapSize = argumentMaps.size();

        // generate signature map for each argument map

        MapDoubleToInt signatures;
        for( unsigned int kk = 0; kk < initialArgumentMapSize; kk++ ){
            double signature      = getMapStringToDoubleSignature( argumentMaps[kk] ); 
            signatures[signature] = 1;
        }

        // for each current argumment map, add a new argument map w/ (key,value)
        // check that no existing map has the same arguments before adding to the 
        // vector of argument maps

        for( unsigned int kk = 0; kk < initialArgumentMapSize; kk++ ){
            for( unsigned int jj = 0; jj < arguments.size(); jj++ ){
               MapStringToDouble inputArgumentMap = MapStringToDouble(argumentMaps[kk]);
               inputArgumentMap[argumentName]     = arguments[jj];
               double signature = getMapStringToDoubleSignature( inputArgumentMap ); 
               if( signatures.find( signature ) == signatures.end() ){
                   argumentMaps.push_back( inputArgumentMap );
               } else {
                   bool match = 0;
                   for( unsigned int mm = 0; mm < initialArgumentMapSize && !match; mm++ ){
                       match = compareMapStringToDoubles( inputArgumentMap, argumentMaps[mm] );
                   }
                   if( !match ){
                       argumentMaps.push_back( inputArgumentMap );
                   }
               }
            }
        }
    }

    return;
}

/** 
 * Predicate for sorting map[string] = double
 *
 * @param d1 first  MapStringToDouble to compare
 * @param d2 second MapStringToDouble to compare
 *
 */
 
bool TestMapSortPredicate( const MapStringToDouble& d1, const MapStringToDouble& d2 ){
    StringVector sortOrder;
    sortOrder.push_back( "numMolecules" );
    sortOrder.push_back( "nonbondedMethod" );
    sortOrder.push_back( "lambda2" );
    sortOrder.push_back( "boxSize" );
    for( unsigned int ii = 0; ii < sortOrder.size(); ii++ ){
        if( d1.find( sortOrder[ii] ) != d1.end() &&
            d2.find( sortOrder[ii] ) != d2.end() ){
           MapStringToDoubleCI d1i = d1.find( sortOrder[ii] );
           MapStringToDoubleCI d2i = d2.find( sortOrder[ii] );
           if( d1i->second != d2i->second ){
               return d1i->second < d2i->second;
           }
        }
    }
    return false;
}


static CustomNonbondedForce* buildCustomNonbondedSoftcoreForce(  const NonbondedSoftcoreForce& nonbondedSoftcoreForce ){

    CustomNonbondedForce* customNonbonded;
    if( nonbondedSoftcoreForce.getNonbondedMethod() == NoCutoff ){

        customNonbonded          = new CustomNonbondedForce("lambda*4*eps*(dem^2-dem)+138.935456*q/r;"
                                                            "q=q1*q2;"
                                                            "dem=1.0/(soft+rsig);"
                                                            "rsig=(r/sigma)^6;"
                                                            "rsig=(r/sigma)^6;"
                                                            "soft=0.5*(1.0-lambda);"
                                                            "sigma=0.5*(sigma1+sigma2);"
                                                            "eps=sqrt(eps1*eps2);"
                                                            "lambda=min(lambda1,lambda2)");

        customNonbonded->setNonbondedMethod( CustomNonbondedForce::NoCutoff );

    } else {

        customNonbonded          = new CustomNonbondedForce("lambda*4*eps*(dem^2-dem)+138.935456*q*(1.0/r+(krf*r*r)-crf);"
                                                            "q=q1*q2;"
                                                            "dem=1.0/(soft+rsig);"
                                                            "rsig=(r/sigma)^6;"
                                                            "rsig=(r/sigma)^6;"
                                                            "soft=0.5*(1.0-lambda);"
                                                            "sigma=0.5*(sigma1+sigma2);"
                                                            "eps=sqrt(eps1*eps2);"
                                                            "lambda=min(lambda1,lambda2)");

        customNonbonded->setCutoffDistance( nonbondedSoftcoreForce.getCutoffDistance() );
        if( nonbondedSoftcoreForce.getNonbondedMethod() == CutoffNonPeriodic ){
            customNonbonded->setNonbondedMethod( CustomNonbondedForce::CutoffNonPeriodic );
        } else {
            customNonbonded->setNonbondedMethod( CustomNonbondedForce::CutoffPeriodic );
        }

        double cutoffDistance           = nonbondedSoftcoreForce.getCutoffDistance();
        double reactionFieldDielectric  = nonbondedSoftcoreForce.getReactionFieldDielectric();

        double eps2                     = (reactionFieldDielectric - 1.0)/(2.0*reactionFieldDielectric+1.0);
        double kValue                   = eps2/(cutoffDistance*cutoffDistance*cutoffDistance);
        customNonbonded->addGlobalParameter("krf", kValue );

        double cValue                   = (1.0/cutoffDistance)*(3.0*reactionFieldDielectric)/(2.0*reactionFieldDielectric + 1.0); 
        customNonbonded->addGlobalParameter("crf", cValue );
    }

    customNonbonded->addPerParticleParameter("q");
    customNonbonded->addPerParticleParameter("sigma");
    customNonbonded->addPerParticleParameter("eps");
    customNonbonded->addPerParticleParameter("lambda");

    vector<double> nonbondedParams(4);
    for( unsigned int ii = 0; ii < static_cast<unsigned int>(nonbondedSoftcoreForce.getNumParticles()); ii++ ){

        double charge;
        double sigma;
        double epsilon;
        double softcoreLJLambda;
        nonbondedSoftcoreForce.getParticleParameters(ii, charge, sigma, epsilon, softcoreLJLambda);

        nonbondedParams[0] = charge;
        nonbondedParams[1] = sigma;
        nonbondedParams[2] = epsilon;
        nonbondedParams[3] = softcoreLJLambda;
        customNonbonded->addParticle( nonbondedParams );
    }

    return customNonbonded;
}

CustomBondForce* buildCustomBondForceForNonbondedExceptions( const NonbondedSoftcoreForce& nonbondedSoftcoreForce ){

    CustomBondForce* customBond;
    if( nonbondedSoftcoreForce.getNonbondedMethod() == NoCutoff ){

        customBond               = new CustomBondForce("lambda*4*eps*(dem^2-dem)+138.935456*q/r;"
                                                       "dem=1.0/(soft+rsig);"
                                                       "rsig=(r/sigma)^6;"
                                                       "soft=0.5*(1.0-lambda)");

    } else {

        customBond               = new CustomBondForce("withinCutoff*(lambda*4*eps*(dem^2-dem)+138.935456*q*(1.0/r+(krf*r*r)-crf));"
                                                       "withinCutoff=step(cutoff-r);"
                                                       "dem=1.0/(soft+rsig);"
                                                       "rsig=(r/sigma)^6;"
                                                       "soft=0.5*(1.0-lambda)");
 

        double cutoffDistance           = nonbondedSoftcoreForce.getCutoffDistance();
        double reactionFieldDielectric  = nonbondedSoftcoreForce.getReactionFieldDielectric();
        double eps2                     = (reactionFieldDielectric - 1.0)/(2.0*reactionFieldDielectric+1.0);
        double kValue                   = eps2/(cutoffDistance*cutoffDistance*cutoffDistance);
        customBond->addGlobalParameter("krf", kValue );

        double cValue                   = (1.0/cutoffDistance)*(3.0*reactionFieldDielectric)/(2.0*reactionFieldDielectric + 1.0); 
        customBond->addGlobalParameter("crf", cValue );
        customBond->addGlobalParameter("cutoff", cutoffDistance );
    }

    customBond->addPerBondParameter("q");
    customBond->addPerBondParameter("sigma");
    customBond->addPerBondParameter("eps");
    customBond->addPerBondParameter("lambda");

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(nonbondedSoftcoreForce.getNumExceptions()); ii++ ){

        int particle1, particle2;
        double chargeProd;
        double sigma;
        double epsilon;
        double softcoreLJLambda;
        nonbondedSoftcoreForce.getExceptionParameters( ii, particle1, particle2, chargeProd, sigma, epsilon, softcoreLJLambda );

        vector<double> bondParams(4);

        bondParams[0] = chargeProd;
        bondParams[1] = sigma;
        bondParams[2] = epsilon;
        bondParams[3] = softcoreLJLambda;
        customBond->addBond( particle1, particle2, bondParams );
    }

    return customBond;
}

/** 
 * Perform comparison of energies/forces for two systems
 *
 * @param system1                  first  system
 * @param system2                  second system
 * @param platform1                first  platform name (Reference, Cuda, OpenCL)
 * @param platform2                second platform name (Reference, Cuda, OpenCL)
 * @param positions                positions
 * @param inputArgumentMap         arguments/flags (relativeTolerance, applyAssert, ...)
 * @param idString                 id string
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
void runSystemComparisonTest( System& system1, System& system2, 
                              const std::string& platform1, const std::string& platform2,
                              const std::vector<Vec3>& positions, MapStringToDouble& inputArgumentMap,
                              const std::string& idString, FILE* log ){

    int applyAssert                      = 0;
    double relativeTolerance             = 1.0e-04;

    setDoubleFromMapStringToDouble( inputArgumentMap, "relativeTolerance",            relativeTolerance );
    setIntFromMapStringToDouble(    inputArgumentMap, "applyAssert",                  applyAssert ) ;

    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);

    if( log ){
        (void) fprintf( log, "System1: particles=%d forces=%d    System2: particles=%d forces=%d\n",
                        system1.getNumParticles(), system1.getNumForces(),
                        system2.getNumParticles(), system2.getNumForces() );
        (void) fprintf( log, "Positions=%u\n", static_cast<unsigned int>(positions.size()) );
        (void) fprintf( log, "Platform1=%s Platform2=%s\n", platform1.c_str(), platform2.c_str() );
        (void) fprintf( log, "relativeTolerance=%8.2e applyAssert=%d\n", relativeTolerance, applyAssert );

        MapStringInt stringForceVector1;
        MapStringInt stringForceVector2;
        getStringForceMap( system1, stringForceVector1, log );
        (void) fprintf( log, "Forces in system 1: [" );
        for( MapStringIntCI ii = stringForceVector1.begin(); ii != stringForceVector1.end(); ii++ ){
            (void) fprintf( log, " %s ", ii->first.c_str() );
        }

        getStringForceMap( system2, stringForceVector2, log );
        (void) fprintf( log, "]\nForces in system 2: [" );
        for( MapStringIntCI ii = stringForceVector2.begin(); ii != stringForceVector2.end(); ii++ ){
            (void) fprintf( log, " %s ", ii->first.c_str() );
        }
        (void) fprintf( log, "]\n" );
    }

    if( system1.getNumParticles() != system2.getNumParticles() ){
        std::stringstream msg;
        msg << "Number of particles for systems to be compared are unequal: " << system1.getNumParticles() << " != " << system2.getNumParticles();
        throw OpenMMException( msg.str() );
    }
 
    if( system1.getNumParticles() != static_cast<int>(positions.size()) ){
        std::stringstream msg;
        msg << "Number of particles for system does not equal size of position array: " << system1.getNumParticles() << " != " << positions.size();
        throw OpenMMException( msg.str() );
    }
 
    Context context1( system1, integrator1, Platform::getPlatformByName( platform1 ));
    context1.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);

    Context context2( system2, integrator2, Platform::getPlatformByName( platform2 ));
    context2.setPositions(positions);

    State state2 = context2.getState(State::Forces | State::Energy);

    double energyDiff = 0.0;
    if( fabs( state1.getPotentialEnergy() ) > 0.0 || fabs( state2.getPotentialEnergy()) > 0.0 ){
        energyDiff = fabs( state1.getPotentialEnergy() - state2.getPotentialEnergy() )/( fabs( state1.getPotentialEnergy() ) + fabs( state2.getPotentialEnergy() ) );
    }

    if( log ){
        DoubleVector stats;
        compareForcesOfTwoStates( state1, state2, relativeTolerance, stats, log );
        (void) fprintf( log, "%s %6d eDff=%15.7e fMx=%15.7e fAvg=%15.7e fMed=%15.7e eCd=%15.7e eRf=%15.7e mxFIdx=%d\n",
                        idString.c_str(), system1.getNumParticles(), energyDiff,
                        stats[1], stats[0], stats[3], state1.getPotentialEnergy(), state2.getPotentialEnergy(), static_cast<int>(stats[2]+0.0001));
        (void) fflush( log );
    }

    if( applyAssert ){
        ASSERT( energyDiff < relativeTolerance );
        for( int ii = 0; ii < system1.getNumParticles(); ii++ ){
    
            Vec3 f1     = state1.getForces()[ii];
            Vec3 f2     = state2.getForces()[ii];
    
            double f1N  = sqrt( (f1[0]*f1[0]) + (f1[1]*f1[1]) + (f1[2]*f1[2]) );
            double f2N  = sqrt( (f2[0]*f2[0]) + (f2[1]*f2[1]) + (f2[2]*f2[2]) );
    
            double diff = (f1[0]-f2[0])*(f1[0]-f2[0]) +
                          (f1[1]-f2[1])*(f1[1]-f2[1]) +
                          (f1[2]-f2[2])*(f1[2]-f2[2]);
            if( f1N > 0.0 || f1N > 0.0 ){
                diff    = 2.0*sqrt( diff )/(f1N + f2N);
            }
            ASSERT( diff < relativeTolerance );
        }
    }

}

/** 
 * Serialize system
 *
 * @param system                   system to serialize
 * @param serializeFileName        file name for xml output
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
void serializeSystem( System& system, const std::string& serializeFileName, FILE* log ){

#ifdef OPENMM_SERIALIZE
    //registerAmoebaSerializationProxies();
    std::stringstream buffer;
    XmlSerializer::serialize<System>(&system, "System", buffer);
    FILE* filePtr = fopen( serializeFileName.c_str(), "w" );
    if( filePtr == NULL ){
        if( log ){
            (void) fprintf( log, "Unable to open xml file %s\n", serializeFileName.c_str() );
            return;
        }
    }
    (void) fprintf( filePtr, "%s", buffer.str().c_str() );
    (void) fclose( filePtr );
    if( log ){
        (void) fprintf( log, "Wrote system to xml file %s\n", serializeFileName.c_str() );
    }
#endif
    return;
}

/** 
 * Output vector of Vec3 to file
 *
 * @param positions                system to serialize
 * @param fileName                 file name for output
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
void serializeVectorOfVec3( const std::vector<Vec3>& positions, std::string fileName, FILE* log ){
#ifdef OPENMM_SERIALIZE
    FILE* filePtr = fopen( fileName.c_str(), "w" );
    if( filePtr == NULL ){
        if( log ){
            (void) fprintf( log, "Unable to open Vec3 file %s\n", fileName.c_str() );
            return;
        }
    }
    (void) fprintf( filePtr, "Positions  %u\n", static_cast<unsigned int>(positions.size()) );
    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        (void) fprintf( filePtr, "%9u %17.10e %17.10e %17.10e\n", ii, positions[ii][0], positions[ii][1], positions[ii][2] );
    }
    (void) fclose( filePtr );
    if( log ){
        (void) fprintf( log, "Wrote to file %s\n", fileName.c_str() );
    }
#endif
    return;
}

/** 
 * Serialize system and positions
 *
 * @param system                   system to serialize
 * @param positions                positions to output
 * @param baseFileName             base file name for xml/txt output
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
void serializeSystemAndPositions( System& system, const std::vector<Vec3>& positions, const std::string& baseFileName, FILE* log ){

    std::stringstream xmlfileName;
    xmlfileName << baseFileName << ".xml";
    serializeSystem( system, xmlfileName.str(), log );

    std::stringstream posfileName;
    posfileName << baseFileName << ".txt";
    serializeVectorOfVec3( positions, posfileName.str(), log );

    return;
}
#endif // TEST_CUDA_SOFTCORE_FORCE_H_
