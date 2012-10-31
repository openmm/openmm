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
 * This tests the OpenCL implementations of GBVIForce, GBSAOBCForce and the softcore versions of
 * these forces.
 */

#define TEST_NONBONDED  0
#define TEST_OBC        1
#define TEST_GBVI       2

#define TEST_CUDA_PLATFORM   0
#define TEST_OPENCL_PLATFORM 1

#include "openmm/GBVIForce.h"
#include "openmm/GBSAOBCForce.h"
#include "openmm/NonbondedForce.h"

#ifdef USE_SOFTCORE
#include "openmm/GBVISoftcoreForce.h"
#include "openmm/GBSAOBCSoftcoreForce.h"
#include "openmm/NonbondedSoftcoreForce.h"
#endif

/**
 * Utility methods shared across unit tests
 */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include "../src/SimTKUtilities/SimTKOpenMMRealType.h"
#include "OpenMM.h"

#if TEST_PLATFORM == TEST_OPENCL_PLATFORM
#include "ReferencePlatform.h"
#include "OpenCLPlatform.h"
#endif

#if TEST_PLATFORM == TEST_CUDA_PLATFORM
#include "ReferencePlatform.h"
#include "CudaPlatform.h"
#endif
    

#ifdef USE_SOFTCORE
#include "OpenMMFreeEnergy.h"
#include "openmm/freeEnergyKernels.h"
#endif

#include "sfmt/SFMT.h"
#include "openmm/VerletIntegrator.h"

#ifdef OPENMM_SERIALIZE
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/XmlSerializer.h"
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <typeinfo>
#include <iomanip>

extern "C" void registerFreeEnergyCudaKernelFactories();

using namespace OpenMM;
using namespace std;

static const int NoCutoff_OpenMMTest          = 0;
static const int CutoffNonPeriodic_OpenMMTest = 1;
static const int CutoffPeriodic_OpenMMTest    = 2;
static const int Ewald_OpenMMTest             = 3;
static const int PME_OpenMMTest               = 4;

static const int ChargeIndex_OpenMMTest       = 0;
static const int SigmaIndex_OpenMMTest        = 1;
static const int EpsIndex_OpenMMTest          = 2;
static const int GammaIndex_OpenMMTest        = 3;
static const int LambdaIndex_OpenMMTest       = 4;

static const int Reference_OpenMMTest         = 0;
static const int Cuda_OpenMMTest              = 1;
static const int OpenCL_OpenMMTest            = 2;

class BondInfo_OpenMMTest {
public:
     BondInfo_OpenMMTest( int particle1, int particle2, double distance );
     int _particle1;
     int _particle2;
     double _distance;
};

BondInfo_OpenMMTest::BondInfo_OpenMMTest( int particle1, int particle2, double distance ){
    _particle1 = particle1;
    _particle2 = particle2;
    _distance  = distance;
}

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

typedef std::vector<std::vector<double> > VectorOfDoubleVectors;
typedef VectorOfDoubleVectors::iterator VectorOfDoubleVectorsI;
typedef VectorOfDoubleVectors::const_iterator VectorOfDoubleVectorsCI;

typedef std::map< int, int> MapIntInt;
typedef MapIntInt::iterator MapIntIntI;
typedef MapIntInt::const_iterator MapIntIntCI;

typedef std::map< double, int> MapDoubleToInt;
typedef MapDoubleToInt::iterator MapDoubleToIntI;
typedef MapDoubleToInt::const_iterator MapDoubleToIntCI;

typedef std::map< std::string, VectorOfDoubleVectors > MapStringVectorOfDoubleVectors;
typedef MapStringVectorOfDoubleVectors::iterator MapStringVectorOfDoubleVectorsI;
typedef MapStringVectorOfDoubleVectors::const_iterator MapStringVectorOfDoubleVectorsCI;

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
     * @param bondDistance              input bond distance
     * @param array                     output vector of grid values
     *
     * @return -1 if particles will not fit on grid; 0 if they do
     */
     
    int setParticlesOnGrid( const Vec3& origin, const Vec3& boxDimensions, const Vec3& spacing, 
                            OpenMM_SFMT::SFMT& sfmt, double bondDistance, std::vector<Vec3>& array ) const;
    
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

    /** 
     * Show min/max distances between all positions
     *
     * @param positions                   input vector of positions
     * @param periodicBoundaryConditions  if set, use PBC in calculating distances
     * @param showIndex                   number of min/maxentries to show
     *
     */
     
    void showMinMaxDistances( const std::vector<Vec3>& positions, 
                              int periodicBoundaryConditions, int showIndex );
    
    /** 
     * Show min/max distances between a single position and all remaining positions 
     *
     * @param positions                   input vector of positions
     * @param periodicBoundaryConditions  if set, use PBC in calculating distances
     * @param showIndex                   number of min/maxentries to show
     * @param positionIndexVector         list of entries to show min/max distances from
     *
     */
     
    void showMinMaxDistances( const std::vector<Vec3>& positions, 
                              int periodicBoundaryConditions, int showIndex,
                              const IntVector& positionIndexVector );
    
    /** 
     * Show distances between positions
     *
     * @param pairs                       particle indcies for which distance is to be reported
     * @param positions                   input vector of positions
     *
     */
     
    void showDistances( const IntIntPairVector& pairs, const std::vector<Vec3>& positions ) const;
    
    /** 
     * Show particles within a specified distance of a given particle
     *
     * @param positions                   input vector of positions
     * @param periodicBoundaryConditions  if set, use PBC in calculating distances
     * @param particleIndex               particle to check
     * @param distanceToCheckFor          distance to check for
     * @param tolerance                   distance tolerance
     *
     */
     
    void showParticlesWithinDistance( const std::vector<Vec3>& positions, 
                                      int periodicBoundaryConditions, unsigned int particleIndex,
                                      double distanceToCheckFor, double tolerance);
    
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
        for( unsigned int ii = 0; ii < static_cast<unsigned int>(_numParticles); ii += _numParticlesPerMolecule ){ 
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
            double particlesPerDimension  = pow( static_cast<double>(_numMolecules), (1.0/3.0) ); 
            int particlesPerDimensionI    = static_cast<int>(particlesPerDimension+0.999999); 
            double boxSize                = _boxSize;
            double spacingPerDimension    = (boxSize-_bondDistance)/(particlesPerDimension+1.0);
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

        errorFlag = setParticlesOnGrid( origin, boxDimensions, spacing, sfmt, _bondDistance, positions );
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
                                           double bondDistance, std::vector<Vec3>& array ) const {

    const double pi  = 3.14159265358979323846;
    const double pi2 = 2.0*pi;

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
            if( (start[jj]+4.0*bondDistance) > boxDimensions[jj] ){
                start[jj] = origin[jj];
            } else {
                done = true;
            }
        }
        if( !done ){
            std::stringstream msg;
            msg << "PositionGenerator::setParticlesOnGrid error in grid settings at molecule=" << ii;
            throw OpenMMException( msg.str() );
        }
    }

    // add molecule atoms

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(_numMolecules); ii++ ){
        int molecularIndex = ii*_numParticlesPerMolecule;
        for( unsigned int jj = 1; jj < static_cast<unsigned int>(_numParticlesPerMolecule); jj++ ){
            double theta             = genrand_real2(sfmt)*pi2;
            double phi               = genrand_real2(sfmt)*pi;
            array[molecularIndex+jj] = array[molecularIndex] + Vec3(_bondDistance*cos(theta)*cos(phi), _bondDistance*cos(theta)*sin(phi), _bondDistance*sin(theta) );
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
 * Show min/max distances between positions
 *
 * @param positions                   input vector of positions
 * @param periodicBoundaryConditions  if set, use PBC in calculating distances
 * @param showIndex                   number of min/maxentries to show
 * @param positionIndexVector         list of entries to show min/max distances from
 *
 */
 
void PositionGenerator::showMinMaxDistances( const std::vector<Vec3>& positions, 
                                             int periodicBoundaryConditions, int showIndex,
                                             const IntVector& positionIndexVector ){

    if( !_log )return;

    Vec3 box[2];
    getEnclosingBox( positions, box );
    (void) fprintf( _log, "Enclosing Box (in A): [%15.7e %15.7e] [%15.7e %15.7e] [%15.7e %15.7e]   [%15.7e %15.7e %15.7e]\n",
                    box[0][0], box[1][0], box[0][1], box[1][1], box[0][2], box[1][2],
                    (box[1][0] - box[0][0]), (box[1][1] - box[0][1]), (box[1][2] - box[0][2]) );

    for( unsigned int ii = 0; ii < positionIndexVector.size(); ii++ ){
        if( positionIndexVector[ii] < static_cast<int>(positions.size()) ){
            int positionIndex = positionIndexVector[ii];
            IntDoublePairVector sortVector;
            getSortedDistances( periodicBoundaryConditions, positionIndex, positions, sortVector );
            (void) fprintf( _log, "Min/max distance from %6d:\n    ", positionIndex );
            for( unsigned int jj = 0; jj < sortVector.size() && jj < static_cast<unsigned int>(showIndex); jj++ ){
                IntDoublePair pair = sortVector[jj];
                (void) fprintf( _log, "[%6d %15.7e] ", pair.first, pair.second);
            }   
            (void) fprintf( _log, "\n    " );
            for( unsigned int jj = (sortVector.size() - showIndex); jj < sortVector.size() && jj >= 0; jj++ ){
                IntDoublePair pair = sortVector[jj];
                (void) fprintf( _log, "[%6d %15.7e] ", pair.first, pair.second);
            }   
            (void) fprintf( _log, "\n" );
        }
    }

    return;
}

/** 
 * Show min/max distances between positions
 *
 * @param positions                   input vector of positions
 * @param periodicBoundaryConditions  if set, use PBC in calculating distances
 * @param showIndex                   number of min/maxentries to show
 *
 */
 
void PositionGenerator::showMinMaxDistances( const std::vector<Vec3>& positions, 
                                             int periodicBoundaryConditions, int showIndex ){

    if( !_log )return;

    Vec3 box[2];
    getEnclosingBox( positions, box );
    (void) fprintf( _log, "Enclosing Box (in A): [%15.7e %15.7e] [%15.7e %15.7e] [%15.7e %15.7e]   [%15.7e %15.7e %15.7e]\n",
                    box[0][0], box[1][0], box[0][1], box[1][1], box[0][2], box[1][2],
                    (box[1][0] - box[0][0]), (box[1][1] - box[0][1]), (box[1][2] - box[0][2]) );

    IntDoublePairVector hitVector;
    double minDistance       = 1.0e+30;
    double minDistanceCutoff = minDistance*1.1;
    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < positions.size(); jj++ ){
            double distance = periodicBoundaryConditions ? getPeriodicDistance( jj, ii, positions) :  
                                                           getDistance( jj, ii, positions);
            if( distance < minDistanceCutoff ){
                if( distance < minDistance ){
                    minDistance        = distance;
                    minDistanceCutoff  = minDistance*1.1;
                }
                hitVector.push_back( IntDoublePair(ii*positions.size()+jj,distance ) );
            }
        }
    }
    std::sort( hitVector.begin(), hitVector.end(), TestIntDoublePair );
            
    (void) fprintf( _log, "Min distances pbc=%d\n", periodicBoundaryConditions );
    for( unsigned int jj = 0; jj < hitVector.size() && jj < static_cast<unsigned int>(showIndex); jj++ ){
        IntDoublePair pair  = hitVector[jj];
        int index           = pair.first;
        int iIndex          = static_cast<int>(index/positions.size());
        int jIndex          = index - iIndex*positions.size();
        (void) fprintf( _log, "   [%6d %6d %15.7e]\n", iIndex, jIndex, pair.second);
    }   
    return;
}

/** 
 * Show particles within a specified distance of a given particle
 *
 * @param positions                   input vector of positions
 * @param periodicBoundaryConditions  if set, use PBC in calculating distances
 * @param particleIndex               particle to check
 * @param distanceToCheckFor          distance to check for
 * @param tolerance                   distance tolerance
 *
 */
 
void PositionGenerator::showParticlesWithinDistance( const std::vector<Vec3>& positions, 
                                                     int periodicBoundaryConditions, unsigned int particleIndex,
                                                     double distanceToCheckFor, double tolerance){

    if( !_log || particleIndex >= positions.size() )return;

    for( unsigned int ii = 0; ii < positions.size(); ii++ ){
        double distance = periodicBoundaryConditions ? getPeriodicDistance( particleIndex, ii, positions) :  
                                                           getDistance( particleIndex, ii, positions);
        double delta    = fabs( distanceToCheckFor - distance );
        if( ii != particleIndex && delta < tolerance ){
            (void) fprintf( _log, "Distance=%15.7e between particles %u %u.\n", distance, particleIndex, ii);
        }
    }

    return;
}

/** 
 * Show distances between positions
 *
 * @param pairs                       particle indcies for which distance is to be reported
 * @param positions                   input vector of positions
 *
 */
 
void PositionGenerator::showDistances( const IntIntPairVector& pairs, const std::vector<Vec3>& positions ) const {

    for( IntIntPairVectorCI ii = pairs.begin(); ii != pairs.end(); ii++ ){
        if( ii->first < static_cast<int>(positions.size()) && ii->second < static_cast<int>(positions.size()) ){
             double d = getDistance( ii->first, ii->second, positions );
             (void) fprintf( _log, "Distance %6d %6d  %15.7e d2=%15.7e\n", ii->first, ii->second,  d, d*d );
        }   
    }   
    return;

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
        sortVector.push_back( IntDoublePair( ii, distance ) );
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
 * Get relative difference between two forces
 * 
 * @param  f1                   force1
 * @param  f2                   force2
 * @param  forceNorm1           output norm of force1
 * @param  forceNorm2           output norm of force2
 * @param  relativeDiff         output relative difference between force norms
 * @param  log                  if set, output forces
 *
 *
   --------------------------------------------------------------------------------------- */

static void getForceRelativeDifference( const Vec3& f1, const Vec3& f2, double& forceNorm1, double& forceNorm2, 
                                        double& relativeDiff, FILE* log ) {

    double diff     = (f1[0] - f2[0])*(f1[0] - f2[0]) +
                      (f1[1] - f2[1])*(f1[1] - f2[1]) +
                      (f1[2] - f2[2])*(f1[2] - f2[2]); 

    forceNorm1      = sqrt( f1[0]*f1[0] + f1[1]*f1[1] + f1[2]*f1[2] );
    forceNorm2      = sqrt( f2[0]*f2[0] + f2[1]*f2[1] + f2[2]*f2[2] );
 
    if( forceNorm1 > 0.0 || forceNorm2 > 0.0 ){
        relativeDiff = 2.0*sqrt( diff )/(forceNorm1+forceNorm2);
    } else {
        relativeDiff = 0.0;
    }

    return;
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
    double averageRelativeDifference      = 0.0;
    double count                          = 0.0;

    DoubleVector medians1( force1.size() );
    DoubleVector medians2( force1.size() );

    IntDoublePairVector relativeDifferences;

    for( unsigned int ii = 0; ii < force1.size(); ii++ ){

        double forceNorm1;
        double forceNorm2;
        double relativeDiff;
        getForceRelativeDifference( force1[ii], force2[ii], forceNorm1, forceNorm2, relativeDiff, log );

        medians1[ii]               = forceNorm1;
        medians2[ii]               = forceNorm2;
 
        relativeDifferences.push_back( IntDoublePair(ii, relativeDiff ) );
        averageRelativeDifference += relativeDiff;
        count                     += 1.0;

        if( relativeDiff > relativeTolerance ){
           error++;
        }

        if( log ){
            (void) fprintf( log, "F %6u %15.7e [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] %15.7e %15.7e %s\n", static_cast<unsigned int>(ii), 
                            relativeDiff, force1[ii][0], force1[ii][1], force1[ii][2], force2[ii][0], force2[ii][1], force2[ii][2],
                            forceNorm1, forceNorm2, (relativeDiff < relativeTolerance ? "":"XXXXXX") );
        }
    }

    // sort relative differences

    std::sort( relativeDifferences.begin(), relativeDifferences.end(), TestIntDoublePair );

    if( log ){
        (void) fprintf( log, "\nEntries w/ largest relative differences.\n" );
        for( unsigned int ii = relativeDifferences.size()-1; ii >= relativeDifferences.size()-10 && ii >= 0; ii-- ){
            double forceNorm1;
            double forceNorm2;
            double relativeDiff;
            int index = relativeDifferences[ii].first;
            getForceRelativeDifference( force1[index], force2[index], forceNorm1, forceNorm2, relativeDiff, log );
            (void) fprintf( log, "Fs %6u %15.7e [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] %15.7e %15.7e %s\n",
                            static_cast<unsigned int>(index), relativeDiff, 
                            force1[index][0], force1[index][1], force1[index][2],
                            force2[index][0], force2[index][1], force2[index][2], 
                            forceNorm1, forceNorm2, (relativeDiff < relativeTolerance ? "":"XXXXXX") );
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
    stats[0]           = averageRelativeDifference;
    IntDoublePair pair = relativeDifferences[relativeDifferences.size()-1];
    stats[1]           = pair.second;
    stats[2]           = static_cast<double>(pair.first);
    stats[3]           = median1 < median2 ? median1 : median2;
    
    return error;
}

/** 
 * Create nonbonded force and set some parameters
 *
 * @param nonbondedMethod          nonbonded method
 * @param cutoffDistance           cutoff distance
 * @param reactionFieldDielectric  reaction field dielectric
 * @param parameterList            list of parameters -- used via addParticle()
 * @param bonds                    list of BondInfo_OpenMMTest containing info for exceptions
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static NonbondedForce* getNonbondedForce( int nonbondedMethod, double cutoffDistance, double reactionFieldDielectric,
                                          VectorOfDoubleVectors& parameterList, std::vector< BondInfo_OpenMMTest >& bonds, FILE* log ){

    
    NonbondedForce* nonbondedForce = new NonbondedForce();
    NonbondedForce::NonbondedMethod method;
    switch( nonbondedMethod ){
        case NoCutoff_OpenMMTest:
            method = NonbondedForce::NoCutoff;
            break;
        case CutoffNonPeriodic_OpenMMTest:
            method = NonbondedForce::CutoffNonPeriodic;
            break;
        case CutoffPeriodic_OpenMMTest:
            method = NonbondedForce::CutoffPeriodic;
            break;
        case Ewald_OpenMMTest:
            method = NonbondedForce::Ewald;
            break;
        case PME_OpenMMTest:
            method = NonbondedForce::PME;
            break;
        default:
            method = NonbondedForce::NoCutoff;
    }
    nonbondedForce->setNonbondedMethod( method );
    nonbondedForce->setCutoffDistance( cutoffDistance );
    nonbondedForce->setReactionFieldDielectric( reactionFieldDielectric );
 
    // load parameters

    for( unsigned int ii = 0; ii < parameterList.size(); ii++ ){
        DoubleVector parameters = parameterList[ii]; 
        nonbondedForce->addParticle( parameters[ChargeIndex_OpenMMTest],  parameters[SigmaIndex_OpenMMTest],  parameters[EpsIndex_OpenMMTest] );
    }

    // add exceptions

    for( unsigned int ii = 0; ii < bonds.size(); ii++ ){
        BondInfo_OpenMMTest bond = bonds[ii]; 
        nonbondedForce->addException( bond._particle1, bond._particle2,  0.0f, 1.0, 0.0f  );
    }

    return nonbondedForce;

}

/** 
 * Create GBVI force and set some parameters
 *
 * @param nonbondedMethod             nonbonded method
 * @param cutoffDistance              cutoff distance
 * @param useQuinticSpline            if set use quintic spline for Born radii scaling; else use no scaling
 * @param quinticLowerLimitFactor     quintic lower limit factor 
 * @param quinticUpperBornRadiusLimit quintic upper Born radius limit
 * @param solventDielectric           solvent dielectric
 * @param soluteDielectric            solute dielectric
 * @param parameterList            list of parameters -- used via addParticle()
 * @param bonds                    list of BondInfo_OpenMMTest containing info for exceptions
 * @param log                         logging file (optional -- may be NULL)
 *
 */
 
static GBVIForce* getGBVIForce( int nonbondedMethod, double cutoffDistance, int useQuinticSpline,
                                double quinticLowerLimitFactor, double quinticUpperBornRadiusLimit,
                                double solventDielectric, double soluteDiecletric,
                                VectorOfDoubleVectors& parameterList, std::vector< BondInfo_OpenMMTest >& bonds, FILE* log ){

    GBVIForce* gbviForce = new GBVIForce();
    GBVIForce::NonbondedMethod method;
    switch( nonbondedMethod ){
        case NoCutoff_OpenMMTest:
            method = GBVIForce::NoCutoff;
            break;
        case CutoffNonPeriodic_OpenMMTest:
            method = GBVIForce::CutoffNonPeriodic;
            break;
        case CutoffPeriodic_OpenMMTest:
            method = GBVIForce::CutoffPeriodic;
            break;
        default:
            method = GBVIForce::NoCutoff;
    }
    gbviForce->setNonbondedMethod( method );
    gbviForce->setCutoffDistance( cutoffDistance );
    gbviForce->setSolventDielectric( solventDielectric );
    gbviForce->setSoluteDielectric( soluteDiecletric );

    if( useQuinticSpline ){
        gbviForce->setBornRadiusScalingMethod( GBVIForce::QuinticSpline );
        gbviForce->setQuinticLowerLimitFactor( quinticLowerLimitFactor );
        gbviForce->setQuinticUpperBornRadiusLimit( quinticUpperBornRadiusLimit );
    } else {
        gbviForce->setBornRadiusScalingMethod( GBVIForce::NoScaling );
    }   

    // load parameters

    for( unsigned int ii = 0; ii < parameterList.size(); ii++ ){
        DoubleVector parameters = parameterList[ii]; 
        gbviForce->addParticle( parameters[ChargeIndex_OpenMMTest],  parameters[SigmaIndex_OpenMMTest],  parameters[GammaIndex_OpenMMTest] );
    }

    // add exceptions

    for( unsigned int ii = 0; ii < bonds.size(); ii++ ){
        BondInfo_OpenMMTest bond = bonds[ii]; 
        gbviForce->addBond( bond._particle1, bond._particle2,  bond._distance);
    }

    return gbviForce;

}

/** 
 * Create OBC force and set some parameters
 *
 * @param nonbondedMethod          nonbonded method
 * @param cutoffDistance           cutoff distance
 * @param solventDielectric        solvent dielectric
 * @param soluteDielectric         solute dielectric
 * @param parameterList            list of parameters -- used via addParticle()
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static GBSAOBCForce* getGBSAOBCForce( int nonbondedMethod, double cutoffDistance, double solventDielectric, double soluteDiecletric,
                                      VectorOfDoubleVectors& parameterList, FILE* log ){

    
    GBSAOBCForce* obcForce = new GBSAOBCForce();
    GBSAOBCForce::NonbondedMethod method;
    switch( nonbondedMethod ){
        case NoCutoff_OpenMMTest:
            method = GBSAOBCForce::NoCutoff;
            break;
        case CutoffNonPeriodic_OpenMMTest:
            method = GBSAOBCForce::CutoffNonPeriodic;
            break;
        case CutoffPeriodic_OpenMMTest:
            method = GBSAOBCForce::CutoffPeriodic;
            break;
        default:
            method = GBSAOBCForce::NoCutoff;
    }
    obcForce->setNonbondedMethod( method );
    obcForce->setCutoffDistance( cutoffDistance );
    obcForce->setSolventDielectric( solventDielectric );
    obcForce->setSoluteDielectric( soluteDiecletric );

    // load parameters

    for( unsigned int ii = 0; ii < parameterList.size(); ii++ ){
        DoubleVector parameters = parameterList[ii]; 
        obcForce->addParticle( parameters[ChargeIndex_OpenMMTest],  parameters[SigmaIndex_OpenMMTest],  parameters[GammaIndex_OpenMMTest] );
    }

    return obcForce;

}

/** 
 * Create nonbonded softcore force and set some parameters
 *
 * @param nonbondedMethod          nonbonded method
 * @param cutoffDistance           cutoff distance
 * @param reactionFieldDielectric  reaction field dielectric
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
#ifdef USE_SOFTCORE
static NonbondedSoftcoreForce* getNonbondedSoftcoreForce( int nonbondedMethod, double cutoffDistance, double reactionFieldDielectric, FILE* log ){

    
    NonbondedSoftcoreForce* nonbondedForce = new NonbondedSoftcoreForce();
    NonbondedSoftcoreForce::NonbondedMethod method;
    switch( nonbondedMethod ){
        case NoCutoff_OpenMMTest:
            method = NonbondedSoftcoreForce::NoCutoff;
            break;
        case CutoffNonPeriodic_OpenMMTest:
            method = NonbondedSoftcoreForce::CutoffNonPeriodic;
            break;
        case CutoffPeriodic_OpenMMTest:
            method = NonbondedSoftcoreForce::CutoffPeriodic;
            break;
        default:
            method = NonbondedSoftcoreForce::NoCutoff;
    }
    nonbondedForce->setNonbondedMethod( method );
    nonbondedForce->setCutoffDistance( cutoffDistance );
    nonbondedForce->setReactionFieldDielectric( reactionFieldDielectric );

    return nonbondedForce;

}

/** 
 * Create GBVI softcore force and set some parameters
 *
 * @param nonbondedMethod             nonbonded method
 * @param cutoffDistance              cutoff distance
 * @param useQuinticSpline            if set use quintic spline for Born radii scaling; else use no scaling
 * @param quinticLowerLimitFactor     quintic lower limit factor 
 * @param quinticUpperBornRadiusLimit quintic upper Born radius limit
 * @param solventDielectric           solvent dielectric
 * @param soluteDielectric            solute dielectric
 * @param log                         logging file (optional -- may be NULL)
 *
 */
 
static GBVISoftcoreForce* getGBVISoftcoreForce( int nonbondedMethod, double cutoffDistance, int useQuinticSpline,
                                                double quinticLowerLimitFactor, double quinticUpperBornRadiusLimit,
                                                double solventDielectric, double soluteDiecletric, FILE* log ){

    GBVISoftcoreForce* gbviForce = new GBVISoftcoreForce();
    GBVISoftcoreForce::NonbondedMethod method;
    switch( gbviMethod ){
        case NoCutoff_OpenMMTest:
            method = GBVISoftcoreForce::NoCutoff;
            break;
        case CutoffNonPeriodic_OpenMMTest:
            method = GBVISoftcoreForce::CutoffNonPeriodic;
            break;
        case CutoffPeriodic_OpenMMTest:
            method = GBVISoftcoreForce::CutoffPeriodic;
            break;
        default:
            method = GBVISoftcoreForce::NoCutoff;
    }
    gbviForce->setNonbondedMethod( method );
    gbviForce->setCutoffDistance( cutoffDistance );
    gbviForce->setSolventDielectric( solventDielectric );
    gbviForce->setSoluteDielectric( soluteDiecletric );

    if( useQuinticSpline ){
        gbviForce->setBornRadiusScalingMethod( GBVISoftcoreForce::QuinticSpline );
        gbviForce->setQuinticLowerLimitFactor( quinticLowerLimitFactor );
        gbviForce->setQuinticUpperBornRadiusLimit( quinticUpperBornRadiusLimit );
    } else {
        gbviForce->setBornRadiusScalingMethod( GBVISoftcoreForce::NoScaling );
    }   

    return gbviForce;

}

/** 
 * Create OBC softcore force and set some parameters
 *
 * @param nonbondedMethod          nonbonded method
 * @param cutoffDistance           cutoff distance
 * @param solventDielectric        solvent dielectric
 * @param soluteDielectric         solute dielectric
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static GBSAOBCSoftcoreForce* getGBSAOBCSoftcoreForce( int nonbondedMethod, double cutoffDistance,
                                                      double solventDielectric, double soluteDiecletric, FILE* log ){

    
    GBSAOBCSoftcoreForce* obcForce = new GBSAOBCSoftcoreForce();
    GBSAOBCSoftcoreForce::NonbondedMethod method;
    switch( gbviMethod ){
        case NoCutoff_OpenMMTest:
            method = GBSAOBCSoftcoreForce::NoCutoff;
            break;
        case CutoffNonPeriodic_OpenMMTest:
            method = GBSAOBCSoftcoreForce::CutoffNonPeriodic;
            break;
        case CutoffPeriodic_OpenMMTest:
            method = GBSAOBCSoftcoreForce::CutoffPeriodic;
            break;
        default:
            method = GBSAOBCSoftcoreForce::NoCutoff;
    }
    obcForce->setNonbondedMethod( method );
    obcForce->setCutoffDistance( cutoffDistance );
    obcForce->setSolventDielectric( solventDielectric );
    obcForce->setSoluteDielectric( soluteDiecletric );

    return obcForce;

}

#endif

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

/** 
 * Get forces in system
 *
 * @param system                   system to serialize
 * @param stringForceVector        output stringForceVector[forceName] = force index
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static Force* copyForce( const Force& force, FILE* log ){

    // print active forces and relevant parameters

    Force* forceCopy = NULL;
    try {
        const CMAPTorsionForce& castForce = dynamic_cast<const CMAPTorsionForce&>(force);
        forceCopy                         = new CMAPTorsionForce( castForce );
    } catch( std::bad_cast ){
    }

    if( forceCopy == NULL ){

        try {
            const CustomAngleForce& castForce = dynamic_cast<const CustomAngleForce&>(force);
            forceCopy                         = new CustomAngleForce( castForce );
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const CustomBondForce& castForce = dynamic_cast<const CustomBondForce&>(force);
           forceCopy                        = new CustomBondForce( castForce );
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const CustomExternalForce& castForce = dynamic_cast<const CustomExternalForce&>(force);
           forceCopy                            = new CustomExternalForce( castForce );
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const CustomGBForce& castForce = dynamic_cast<const CustomGBForce&>(force);
           forceCopy                      = new CustomGBForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }


    if( forceCopy == NULL ){
        try {
           const CustomHbondForce& castForce = dynamic_cast<const CustomHbondForce&>(force);
           forceCopy                         = new CustomHbondForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }


    if( forceCopy == NULL ){
        try {
           const CustomNonbondedForce& castForce = dynamic_cast<const CustomNonbondedForce&>(force);
           forceCopy                             = new CustomNonbondedForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }



    if( forceCopy == NULL ){
        try {
           const CustomTorsionForce& castForce = dynamic_cast<const CustomTorsionForce&>(force);
           forceCopy                           = new CustomTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }



    if( forceCopy == NULL ){
        try {
           const GBSAOBCForce& castForce = dynamic_cast<const GBSAOBCForce&>(force);
           forceCopy                     = new GBSAOBCForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const GBVIForce& castForce = dynamic_cast<const GBVIForce&>(force);
           forceCopy                  = new GBVIForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const HarmonicAngleForce& castForce = dynamic_cast<const HarmonicAngleForce&>(force);
           forceCopy                           = new HarmonicAngleForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){
        try {
           const HarmonicBondForce& castForce = dynamic_cast<const HarmonicBondForce&>(force);
           forceCopy                          = new HarmonicBondForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){
        try {
           const NonbondedForce& castForce = dynamic_cast<const NonbondedForce&>(force);
           forceCopy                       = new NonbondedForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const PeriodicTorsionForce& castForce = dynamic_cast<const PeriodicTorsionForce&>(force);
           forceCopy                             = new PeriodicTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const RBTorsionForce& castForce = dynamic_cast<const RBTorsionForce&>(force);
           forceCopy                       = new RBTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const MonteCarloBarostat& castForce = dynamic_cast<const MonteCarloBarostat&>(force);
           forceCopy                           = new MonteCarloBarostat( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AndersenThermostat& castForce = dynamic_cast<const AndersenThermostat&>(force);
           forceCopy                           = new AndersenThermostat( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

#ifdef USE_SOFTCORE
    if( forceCopy == NULL ){

        try {
           const GBSAOBCSoftcoreForce& castForce = dynamic_cast<const GBSAOBCSoftcoreForce&>(force);
           forceCopy                             = new GBSAOBCSoftcoreForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){
        try {
           const GBVISoftcoreForce& castForce = dynamic_cast<const GBVISoftcoreForce&>(force);
           forceCopy                          = new GBVISoftcoreForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){
        try {
           const NonbondedSoftcoreForce& castForce = dynamic_cast<const NonbondedSoftcoreForce&>(force);
           forceCopy                               = new NonbondedSoftcoreForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }
#endif

#ifdef INCLUDE_AMOEBA_FORCES

    if( forceCopy == NULL ){

        try {
           const AmoebaHarmonicBondForce& castForce = dynamic_cast<const AmoebaHarmonicBondForce&>(force);
           forceCopy                                = new AmoebaHarmonicBondForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaHarmonicAngleForce& castForce = dynamic_cast<const AmoebaHarmonicAngleForce&>(force);
           forceCopy                                 = new AmoebaHarmonicAngleForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaHarmonicInPlaneAngleForce& castForce = dynamic_cast<const AmoebaHarmonicInPlaneAngleForce&>(force);
           forceCopy                                        = new AmoebaHarmonicInPlaneAngleForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaMultipoleForce& castForce = dynamic_cast<const AmoebaMultipoleForce&>(force);
           forceCopy                             = new AmoebaMultipoleForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaOutOfPlaneBendForce& castForce = dynamic_cast<const AmoebaOutOfPlaneBendForce&>(force);
           forceCopy                                  = new AmoebaOutOfPlaneBendForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaPiTorsionForce& castForce = dynamic_cast<const AmoebaPiTorsionForce&>(force);
           forceCopy                             = new AmoebaPiTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaStretchBendForce& castForce = dynamic_cast<const AmoebaStretchBendForce&>(force);
           forceCopy                               = new AmoebaStretchBendForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaTorsionForce& castForce = dynamic_cast<const AmoebaTorsionForce&>(force);
           forceCopy                           = new AmoebaTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaTorsionTorsionForce& castForce = dynamic_cast<const AmoebaTorsionTorsionForce&>(force);
           forceCopy                                  = new AmoebaTorsionTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaUreyBradleyForce& castForce = dynamic_cast<const AmoebaUreyBradleyForce&>(force);
           forceCopy                               = new AmoebaUreyBradleyForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaVdwForce& castForce = dynamic_cast<const AmoebaVdwForce&>(force);
           forceCopy                       = new AmoebaVdwForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaWcaDispersionForce& castForce = dynamic_cast<const AmoebaWcaDispersionForce&>(force);
           forceCopy                                 = new AmoebaWcaDispersionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaGeneralizedKirkwoodForce& castForce = dynamic_cast<const AmoebaGeneralizedKirkwoodForce&>(force);
           forceCopy                                       = new AmoebaGeneralizedKirkwoodForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

    if( forceCopy == NULL ){

        try {
           const AmoebaTorsionTorsionForce& castForce = dynamic_cast<const AmoebaTorsionTorsionForce&>(force);
           forceCopy                                  = new AmoebaTorsionTorsionForce( castForce ); 
        } catch( std::bad_cast ){
        }    
    }

#endif

    if( log && forceCopy == NULL ){
       (void) fprintf( log, " force not recognized.\n" );
    }

    return forceCopy;
}

/** 
 * Return copy of system (but not forces)
 *
 * @param inputSystem               system to copy
 *
 * @return system copy
 *
 */
 
static void copySystem( const System& inputSystem, System& systemCopy, FILE* log ){

    // add particle/mass

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(inputSystem.getNumParticles()); ii++ ){
        systemCopy.addParticle( inputSystem.getParticleMass( static_cast<int>(ii) ) );
    }

    // box

    Vec3 a;
    Vec3 b;
    Vec3 c;
    inputSystem.getDefaultPeriodicBoxVectors( a, b, c );
    systemCopy.setDefaultPeriodicBoxVectors( a, b, c );

    // copy constraints

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(inputSystem.getNumConstraints()); ii++ ){
        int particle1, particle2;
        double distance;
        inputSystem.getConstraintParameters( ii, particle1, particle2, distance);
        systemCopy.addConstraint( particle1, particle2, distance);
    }

    // copy forces

    for( unsigned int ii = 0; ii < static_cast<unsigned int>(inputSystem.getNumForces()); ii++ ){
        systemCopy.addForce( copyForce( inputSystem.getForce(ii), log) );
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
                } else if( key == "cutoffDistance" ){
                    formatArgument( buffer, key, value, "%7.3f ", callId, 1 ); 
                } else if( key == "relativeTolerance" ){
                    formatArgument( buffer, key, value, "%8.1e ", callId, 1 ); 
                } else if( key == "positionPlacementMethod" || key == "applyAssert" || key == "serialize" ){
                    formatArgument( buffer, key, value, "%1d ", callId, 1 ); 
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
            } else if( key == "nonbondedMethod" || key == "positionPlacementMethod" || key == "applyAssert" || key == "serialize" ){
                (void) sprintf( buffer, "%s=%1d ", key.c_str(), valueInt );
            } else if( key == "lambda1" || key == "lambda2" ){
                (void) sprintf( buffer, "%s=%4.2f ", key.c_str(), value );
            } else if( key == "boxSize" || key == "cutoffDistance" ){
                (void) sprintf( buffer, "%s=%6.2f ", key.c_str(), value );
            } else if( key == "relativeTolerance" ){
                (void) sprintf( buffer, "%s=%8.1e ", key.c_str(), value );
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


#ifdef USE_SOFTCORE
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
    for( unsigned int ii = 0; ii < nonbondedSoftcoreForce.getNumParticles(); ii++ ){

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

    for( unsigned int ii = 0; ii < nonbondedSoftcoreForce.getNumExceptions(); ii++ ){

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
#endif

/** 
 * Load plugins
 *
 * @param pluginDirectory          plugin directory; if OPENMM_PLUGIN_DIR use ENV variable
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static int loadPlugins( const std::string& pluginDirectory, std::vector<std::string>& loaded, FILE* log ){

    const char* openmmPluginDirectory;
    int envVariableIsSet = 0;
    if( pluginDirectory.compare( "OPENMM_PLUGIN_DIR") == 0 ){
         openmmPluginDirectory  =  getenv( "OPENMM_PLUGIN_DIR" );
    } else {
         openmmPluginDirectory  = pluginDirectory.c_str();
    }
    try {
        if( openmmPluginDirectory ){
           envVariableIsSet = 1;
           if( log ){
               (void) fprintf( log, "openmmPluginDirectory=%s\n", openmmPluginDirectory );
               (void) fflush( log );
           }    

           loaded = Platform::loadPluginsFromDirectory( openmmPluginDirectory );

           if( log ){
               (void) fprintf( log, "\nLoaded following %u lib(s) from %s:\n", static_cast<unsigned int>(loaded.size()), openmmPluginDirectory ); (void) fflush( log );
               for( unsigned int ii = 0; ii < loaded.size(); ii++ ){
                   (void) fprintf( log, "   %s\n", loaded[ii].c_str() );
               }    
               (void) fprintf( log, "\n" ); (void) fflush( log );
           }
        } else {
           if( log && pluginDirectory == "OPENMM_PLUGIN_DIR" ){
               (void) fprintf( log, "Env variable OPENMM_PLUGIN_DIR is not set.\n" );
               (void) fflush( log );
           }    
        }    

    } catch(const exception& e) { 
        (void) fprintf( log, "Exception: %s\n", e.what() );
        (void) fflush( log );
    }    

    return envVariableIsSet;
}

/** 
 * Set device id
 *
 * @param platform                 platform
 * @param deviceId                 device id
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static void setDeviceId( Platform& platform, int deviceId, FILE* log ){

    std::stringstream deviceIdStr;
    deviceIdStr << deviceId;
    int wasSet = 0;
    if( platform.getName().compare( "Cuda" ) == 0 ){
        platform.setPropertyDefaultValue( "CudaDevice", deviceIdStr.str() );
        wasSet = 1;

    } else if( platform.getName().compare( "OpenCL" ) == 0 ){
        platform.setPropertyDefaultValue( "OpenCLDeviceIndex",  deviceIdStr.str());
        wasSet = 1;
    }

    if( log && wasSet ){
        (void) fprintf( log, "Set deviceId to %d\n", deviceId );
        (void) fflush( log );
    }

    return;
}

/** 
 * Set device id
 *
 * @param platform                 platform
 * @param deviceId                 device id
 * @param log                      logging file (optional -- may be NULL)
 *
 */
 
static void setDeviceIdUsingEnvVariable( Platform& platform, FILE* log ){

    const char* deviceId = getenv( "GPU_DEVICE_ID" );
    if( deviceId == NULL ){
       return;
    }
    int wasSet = 0;
    if( platform.getName().compare( "Cuda" ) == 0 ){
        platform.setPropertyDefaultValue( "CudaDevice", deviceId );
        wasSet = 1;

    } else if( platform.getName().compare( "OpenCL" ) == 0 ){
        platform.setPropertyDefaultValue( "OpenCLDeviceIndex",  deviceId);
        wasSet = 1;
    }

    if( log && wasSet ){
        (void) fprintf( log, "Set deviceId to %s based on env variable GPU_DEVICE_ID setting.\n", deviceId );
        (void) fflush( log );
    }

    return;
}

/** 
 * Get platform name
 *
 * @param platformId               platformId( 0=Reference, 1=Cuda, 2=OpenCL)
 * @param platformName             output platform name
 *
 */
 
static void getPlatformName( int platformId, std::string& platformName ){

    switch( platformId ){
        case Reference_OpenMMTest:
            platformName = "Reference";
            break;
        case Cuda_OpenMMTest:
            platformName = "Cuda";
            break;
        case OpenCL_OpenMMTest:
            platformName = "OpenCL";
            break;
        default:
            platformName = "NA";
            break;
    }
    return;
}

/** 
 * Get lib name
 *
 * @param libPrefix                lib prefix (lib or "")
 * @param libSuffix                lib suffix (.so, .dylib, .dll)
 * @param baseName                 base name
 *
 * @return libname
 *
 */
 
static std::string getLibName( const std::string& libPrefix, const std::string& libSuffix, const std::string& baseName ){

    std::string fullName = libPrefix; 
    fullName.append( baseName );
    fullName.append( libSuffix );
    return fullName;
}

/** 
 * Get nonbonded method name
 *
 * @param nonbondedMethod          nonbonded method flag
 * @return nonbonded method name
 *
 */
 
static std::string getNonbondedMethodName( int nonbondedMethod ){

   switch( nonbondedMethod ){
       case NoCutoff_OpenMMTest:
           return "NoCutoff";
       case CutoffNonPeriodic_OpenMMTest:
           return "CutoffNonPeriodic";
       case CutoffPeriodic_OpenMMTest:
           return "CutoffPeriodic";
       case Ewald_OpenMMTest:
           return "Ewald";
       case PME_OpenMMTest:
           return "PME";
       default:
           return "NA";
    }
}

/** 
 * Check if required libs are available
 *
 * @param requiredLibs   list of required libs
 * @param loadedLibs     list of available libs
 * @param log            optional logging reference
 *
 * @return 1 if all required libs are loaded; else 0
 *
 */
 
static int checkRequiredLibsAreAvailable( const StringVector& requiredLibs, const StringVector& loadedLibs, FILE* log ){

    unsigned int matchCount = 0;
    for( unsigned int kk = 0; kk < requiredLibs.size(); kk++ ){
        unsigned int match = 0;
        for( unsigned int ii = 0; ii < loadedLibs.size() && match == 0; ii++ ){
            if( loadedLibs[ii].compare( requiredLibs[kk] ) == 0 ){
                match = 1;
            }
        }
        if( log && !match ){
            (void) fprintf( log, "Missing lib %s\n", requiredLibs[kk].c_str() );
        }
        matchCount += match;
    }

    int allPresent;
    if( matchCount < requiredLibs.size() ){
        allPresent = 0;
        if( log ){
            (void) fprintf( log, "Aborting tests due to missing libs.\n" );
        }
    } else {
        allPresent = 1;
    }
    return allPresent;
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
                              const std::vector<Vec3>& positions, MapStringToDouble& inputArgumentMap,
                              const std::string& idString, const std::string& precision, FILE* log ){

    int applyAssert                      = 0;
    int platformId1                      = 0;
    int platformId2                      = 0;
    int deviceId1                        = 0;
    int deviceId2                        = 0;
    double relativeTolerance             = 1.0e-04;

    setDoubleFromMapStringToDouble( inputArgumentMap, "relativeTolerance",            relativeTolerance );
    setIntFromMapStringToDouble(    inputArgumentMap, "applyAssert",                  applyAssert ) ;
    setIntFromMapStringToDouble(    inputArgumentMap, "platformId1",                  platformId1 ) ;
    setIntFromMapStringToDouble(    inputArgumentMap, "platformId2",                  platformId2 ) ;
    setIntFromMapStringToDouble(    inputArgumentMap, "deviceId1",                    deviceId1 ) ;
    setIntFromMapStringToDouble(    inputArgumentMap, "deviceId2",                    deviceId2 ) ;

    std::string platformName1;
    std::string platformName2;
    getPlatformName( platformId1, platformName1 );
    getPlatformName( platformId2, platformName2 );

    VerletIntegrator integrator1(0.01);
    VerletIntegrator integrator2(0.01);

    if( log ){
        (void) fprintf( log, "System1: particles=%d forces=%d    System2: particles=%d forces=%d\n",
                        system1.getNumParticles(), system1.getNumForces(),
                        system2.getNumParticles(), system2.getNumForces() );
        (void) fprintf( log, "Positions=%u\n", static_cast<unsigned int>(positions.size()) );
        (void) fprintf( log, "Platform1=%s Platform2=%s\n", platformName1.c_str(), platformName2.c_str() );
        (void) fprintf( log, "deviceId1=%d deviceId2=%d\n", deviceId1, deviceId2 );
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
 
#if TEST_PLATFORM == TEST_OPENCL_PLATFORM
    ReferencePlatform platform1;   
    OpenCLPlatform platform2;
    platform2.setPropertyDefaultValue("OpenCLPrecision", precision);
#elif TEST_PLATFORM == TEST_CUDA_PLATFORM
    ReferencePlatform platform1;   
    CudaPlatform platform2;
    platform2.setPropertyDefaultValue("CudaPrecision", precision);
#else
    Platform& platform1 = Platform::getPlatformByName( platformName1 );
    if( deviceId1 ){
        setDeviceId( platform1, deviceId1, log );
    }
    setDeviceIdUsingEnvVariable( platform1, log );

    Platform& platform2 = Platform::getPlatformByName( platformName2 );
    if( deviceId2 ){
        setDeviceId( platform2, deviceId2, log );
    }
    setDeviceIdUsingEnvVariable( platform2, log );
#endif

    Context context1( system1, integrator1, platform1 );
    context1.setPositions(positions);
    State state1 = context1.getState(State::Forces | State::Energy);

    Context context2( system2, integrator2, platform2 );
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

void runTests( MapStringToDouble& inputArgumentMap, const std::string& precision, FILE* log ){

    double lambda1                       = 1.0;
    double lambda2                       = 1.0;
    int nonbondedMethod                  = 0;
    int numMolecules                     = 1;
    int numParticlesPerMolecule          = 2;
    int useQuinticSpline                 = 1;
    int applyAssert                      = 1;
    int positionPlacementMethod          = 0;
    int serialize                        = 0;
    double boxSize                       = 10.0;
    double relativeTolerance             = 1.0e-04;
    double quinticLowerLimitFactor       = 0.8;
    double quinticUpperBornRadiusLimit   = 2.0;
    std::stringstream baseFileName;

    setDoubleFromMapStringToDouble( inputArgumentMap, "lambda1",                      lambda1 );
    setDoubleFromMapStringToDouble( inputArgumentMap, "lambda2",                      lambda2 );
    setDoubleFromMapStringToDouble( inputArgumentMap, "boxSize",                      boxSize );
    double cutoffDistance                = boxSize*0.4;
    setDoubleFromMapStringToDouble( inputArgumentMap, "cutoffDistance",               cutoffDistance);
    setDoubleFromMapStringToDouble( inputArgumentMap, "relativeTolerance",            relativeTolerance );

    baseFileName << "Nb";
#if IMPLICIT_SOLVENT == TEST_GBVI
    setDoubleFromMapStringToDouble( inputArgumentMap, "quinticLowerLimitFactor",      quinticLowerLimitFactor );
    setDoubleFromMapStringToDouble( inputArgumentMap, "quinticUpperBornRadiusLimit",  quinticUpperBornRadiusLimit );
    baseFileName << "Gbvi";
#endif
#if IMPLICIT_SOLVENT == TEST_OBC
    baseFileName << "Obc";
#endif

    setIntFromMapStringToDouble(    inputArgumentMap, "positionPlacementMethod",      positionPlacementMethod ) ;
    setIntFromMapStringToDouble(    inputArgumentMap, "nonbondedMethod",              nonbondedMethod );
    setIntFromMapStringToDouble(    inputArgumentMap, "numMolecules",                 numMolecules );
    setIntFromMapStringToDouble(    inputArgumentMap, "numParticlesPerMolecule",      numParticlesPerMolecule );
    setIntFromMapStringToDouble(    inputArgumentMap, "serialize",                    serialize );
    setIntFromMapStringToDouble(    inputArgumentMap, "applyAssert",                  applyAssert );
   
    double bondDistance = 0.1;
    double minDistance  = 0.1;
    double cellSize     = 2.0*bondDistance + minDistance;
    double boxLength    = cellSize*pow( static_cast<double>(numMolecules), 0.333333 );
    if( positionPlacementMethod == 1 && boxLength > boxSize ){
        boxSize = boxLength;
        if( log ){
            // (void) fprintf( log, "Updated box size: bL=%6.3f cell=%6.2e bond=%5.2f separation=%5.2f\n", boxLength, cellSize, bondDistance, minDistance );
        }
    }

    if( nonbondedMethod >= 2 && cutoffDistance > boxSize*0.5 ){
        cutoffDistance = boxSize*0.49;
    }

    int numParticles                     = numMolecules*numParticlesPerMolecule;
    int includeGbvi                      = 1;
    double reactionFieldDielectric       = 80.0;

    if( log ){
        double particleDensity = static_cast<double>(numParticles)/(boxSize*boxSize*boxSize);
        double particleCube    = pow( particleDensity, (-1.0/3.0) );
      
        (void) fprintf( log, "\n--------------------------------------------------------------------------------------\n" );
        (void) fprintf( log, "Input arguments\n" );
        (void) fflush( log );
        //(void) fprintf( log, "    includeGbvi                 %d\n", includeGbvi );
        (void) fprintf( log, "    nonbondedMethod             %d\n", nonbondedMethod );
        (void) fprintf( log, "    numParticles                %d\n", numParticles );
        (void) fprintf( log, "    numMolecules                %d\n", numMolecules );
        (void) fprintf( log, "    numParticlesPerMolecule     %d\n", numParticlesPerMolecule );
        (void) fprintf( log, "    positionPlacementMethod     %d\n", positionPlacementMethod);
        (void) fprintf( log, "    boxSize                     %8.3f\n", boxSize );
        (void) fprintf( log, "    cutoffDistance              %15.7e\n", cutoffDistance );
        (void) fprintf( log, "    reactionFieldDielectric     %8.3f\n", reactionFieldDielectric );

#if IMPLICIT_SOLVENT == TEST_GBVI
        (void) fprintf( log, "    useQuinticSpline            %d\n", useQuinticSpline );
        (void) fprintf( log, "    quinticLowerLimitFactor     %8.3f\n", quinticLowerLimitFactor );
        (void) fprintf( log, "    quinticUpperBornRadiusLimit %8.3f\n", quinticUpperBornRadiusLimit );
#endif
#ifdef USE_SOFTCORE
        (void) fprintf( log, "    lambda1                     %8.3f\n", lambda1 );
        (void) fprintf( log, "    lambda2                     %8.3f\n", lambda2 );
#endif
        (void) fprintf( log, "    relativeTolerance           %8.1e\n", relativeTolerance );
        (void) fprintf( log, "    particleDensity             %8.2e\n", particleDensity );
        (void) fprintf( log, "    particleCube                %8.2e\n", particleCube );
    }

    // Create two systems: one with GbviSoftcoreForce NonbondedSoftcoreForce forces, and one using a CustomNonbondedForce, CustomGBVI force to implement the same interaction.

    System standardSystem;
    for (int i = 0; i < numParticles; i++) {
        standardSystem.addParticle(1.0);
    }
    standardSystem.setDefaultPeriodicBoxVectors(Vec3(boxSize, 0, 0), Vec3(0, boxSize, 0), Vec3(0, 0, boxSize));

    double solventDielectric = 78.3; // 1.0 or 1.0e+10
    double soluteDiecletric  = 1.0;

    std::vector<Vec3> positions(numParticles);

    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);

    PositionGenerator positionGenerator( numMolecules, numParticlesPerMolecule, boxSize );
    if( log ){
        positionGenerator.setLog( log );
    }
    if( positionPlacementMethod == 1 ){
        positionGenerator.setPositions( PositionGenerator::SimpleGrid, sfmt, positions );
    } else {
        positionGenerator.setBondDistance( 0.3 );
        positionGenerator.setPositions( PositionGenerator::Random, sfmt, positions );
    }

    // show info on particle positions

    if( log ){
        int periodicBoundaryConditions       = (nonbondedMethod > CutoffNonPeriodic_OpenMMTest) ? 1 : 0;
        int showIndex                        = 5;
        double distanceTolerance             = 1.0e-04;
        IntVector positionIndexVector;
        positionIndexVector.push_back( 0 );
        positionIndexVector.push_back( 5713 );
        positionIndexVector.push_back( 6291 );
        positionIndexVector.push_back( 3191 );
        positionIndexVector.push_back( 3769 );
        positionIndexVector.push_back( static_cast<int>(positions.size())-1 );
        positionGenerator.showMinMaxDistances( positions, periodicBoundaryConditions, showIndex, positionIndexVector);
        positionGenerator.showMinMaxDistances( positions, periodicBoundaryConditions, showIndex );
        positionGenerator.showParticlesWithinDistance( positions, periodicBoundaryConditions, 5713, cutoffDistance, distanceTolerance );
        positionGenerator.showParticlesWithinDistance( positions, periodicBoundaryConditions, 6291, cutoffDistance, distanceTolerance );

        IntIntPairVector pairs;
        pairs.push_back( IntIntPair( 5713, 6291 ) );
        pairs.push_back( IntIntPair( 5713, 3191 ) );
        pairs.push_back( IntIntPair( 5713, 3769 ) );
        pairs.push_back( IntIntPair( 6291, 3191 ) );
        positionGenerator.showDistances( pairs, positions );
    }    

    const int numberOfParameters                        = 5;

    std::vector<double> parameterLowerBound( numberOfParameters, 0.0 );

    double fixedCharge                                  = 0.1;
    parameterLowerBound[ChargeIndex_OpenMMTest]         = fixedCharge;  // charge
    parameterLowerBound[SigmaIndex_OpenMMTest]          = 0.1;          // sigma
    parameterLowerBound[EpsIndex_OpenMMTest]            = 0.5;          // eps
    parameterLowerBound[GammaIndex_OpenMMTest]          = 0.1;          // gamma
    parameterLowerBound[LambdaIndex_OpenMMTest]         = lambda1;      // lambda

    std::vector<double> parameterUpperBound( parameterLowerBound );
    parameterUpperBound[ChargeIndex_OpenMMTest]         = fixedCharge;  // charge
    parameterUpperBound[SigmaIndex_OpenMMTest]          = 0.3;          // sigma
    parameterUpperBound[EpsIndex_OpenMMTest]            = 40.0;         // eps
    parameterUpperBound[GammaIndex_OpenMMTest]          = 40.0;         // gamma

#if IMPLICIT_SOLVENT == TEST_OBC
    parameterLowerBound[GammaIndex_OpenMMTest]          = 0.1;          // overlap factor
    parameterUpperBound[GammaIndex_OpenMMTest]          = 1.5;        
#endif

    std::vector<double> parameters( numberOfParameters );
    VectorOfDoubleVectors parameterList;
    std::vector< BondInfo_OpenMMTest > bonds;

    double charge                                       = fixedCharge;

    for( int ii = 0; ii < numMolecules; ii++) {

        charge                 *= -1.0;
        double lambda           =  ii < (numMolecules/2) ? lambda1 : lambda2;

        randomizeParameters( parameterLowerBound, parameterUpperBound, sfmt, parameters );
        parameters[ChargeIndex_OpenMMTest] = charge;
        parameters[LambdaIndex_OpenMMTest] = lambda;
        parameterList.push_back( parameters );

        int baseParticleIndex                    = ii*numParticlesPerMolecule;
        for( int jj = 1; jj < numParticlesPerMolecule; jj++) {

            // alternate charges

            charge *= -1.0;
            randomizeParameters( parameterLowerBound, parameterUpperBound, sfmt, parameters );
            parameters[ChargeIndex_OpenMMTest] = charge;
            parameters[LambdaIndex_OpenMMTest] = lambda;
            parameterList.push_back( parameters );

            double bondDistance  = positionGenerator.getDistance( baseParticleIndex, baseParticleIndex+jj, positions );
            bonds.push_back( BondInfo_OpenMMTest( baseParticleIndex, baseParticleIndex+jj, bondDistance ) );
        }

        // alternate charge if numParticlesPerMolecule is odd

        if( (numParticlesPerMolecule % 2) ){
            charge *= -1.0;
        }
    }

#ifdef USE_SOFTCORE

    baseFileName << "Softcore";
    NonbondedSoftcoreForce* nonbondedSoftcoreForce   = getNonbondedSoftcoreForce( nonbondedMethod, cutoffDistance, reactionFieldDielectric,
                                                                                  parameterList, bonds, log);

#if IMPLICIT_SOLVENT == TEST_GBVI
    GBVISoftcoreForce* gbviSoftcoreForce             = getGBVISoftcoreForce( nonbondedMethod, nonbondedSoftcoreForce->getCutoffDistance(),
                                                                             useQuinticSpline, quinticLowerLimitFactor, quinticUpperBornRadiusLimit,
                                                                             solventDielectric, soluteDiecletric,
                                                                             parameterList, bonds, log );
#endif

#if IMPLICIT_SOLVENT == TEST_OBC
    GBSAOBCSoftcoreForce* gbviSoftcoreForce          = getGBSAOBCSoftcoreForce( nonbondedMethod, nonbondedSoftcoreForce->getCutoffDistance(),
                                                                                solventDielectric, soluteDiecletric,
                                                                                parameterList, log );
#endif



#else




    NonbondedForce*         nonbondedSoftcoreForce   = getNonbondedForce( nonbondedMethod, cutoffDistance, reactionFieldDielectric,
                                                                          parameterList, bonds, log);
#if IMPLICIT_SOLVENT == TEST_GBVI
    GBVIForce* gbviSoftcoreForce                     = getGBVIForce( nonbondedMethod, nonbondedSoftcoreForce->getCutoffDistance(),
                                                                     useQuinticSpline, quinticLowerLimitFactor, quinticUpperBornRadiusLimit,
                                                                     solventDielectric, soluteDiecletric, parameterList, bonds, log );
#endif

#if IMPLICIT_SOLVENT == TEST_OBC
    GBSAOBCForce* gbviSoftcoreForce                  = getGBSAOBCForce( nonbondedMethod, nonbondedSoftcoreForce->getCutoffDistance(),
                                                                        solventDielectric, soluteDiecletric, parameterList, log );
#endif

#endif

    standardSystem.addForce(nonbondedSoftcoreForce);

#if IMPLICIT_SOLVENT > 0
    if( includeGbvi ){
        standardSystem.addForce(gbviSoftcoreForce);
    }
#endif

    // copy system and forces

    System systemCopy;
    copySystem( standardSystem, systemCopy, log );

    // serialize

    baseFileName  << "_N"     << positions.size();
    baseFileName  << "_Nb"    << nonbondedMethod;
    serializeSystemAndPositions( standardSystem, positions, baseFileName.str(), log);

    // perform comparison

    std::stringstream idString;
    idString << "Nb " << nonbondedMethod << " l2 " << std::fixed << setprecision(2) << lambda2;
    runSystemComparisonTest( standardSystem, systemCopy, positions, inputArgumentMap, idString.str(), precision, log );

}

int main(int argc, char* argv[]) {

    std::string precision;
    if (argc > 1)
        precision = std::string(argv[1]);
    else
        precision = "single";
    
    try {

#ifdef USE_SOFTCORE
        registerFreeEnergyCudaKernelFactories( );
#endif

        VectorOfMapStringToDouble vectorOfMapStringToDouble;
        MapStringToDouble inputArgumentMap;
        MapStringToDoubleVector generativeArgumentMaps;
        //FILE* log = stderr;
        FILE* log = NULL;

        std::vector<std::string> loadedLibs;
        int envVariableIsSet = loadPlugins( "OPENMM_PLUGIN_DIR", loadedLibs, log );

        inputArgumentMap["platformId1"]                     = 0;
        inputArgumentMap["platformId2"]                     = 1;

        inputArgumentMap["deviceId1"]                       = 0;
        inputArgumentMap["deviceId2"]                       = 0;

        inputArgumentMap["lambda2"]                         = 1.0;

        inputArgumentMap["nonbondedMethod"]                 = 0;
        inputArgumentMap["numMolecules"]                    = 10;
        inputArgumentMap["boxSize"]                         = 5.0;
        inputArgumentMap["positionPlacementMethod"]         = 1;
        inputArgumentMap["cutoffDistance"]                  = 0.301*inputArgumentMap["boxSize"];
        //inputArgumentMap["cutoffDistance"]                  = 1.0;
        inputArgumentMap["relativeTolerance"]               = 5.0e-04;
        inputArgumentMap["applyAssert"]                     = 1;
        inputArgumentMap["serialize"]                       = 1;
        inputArgumentMap["numParticlesPerMolecule"]         = 2;

#ifdef USE_SOFTCORE
        DoubleVector lamda2;
        lamda2.push_back( 1.0 );
        lamda2.push_back( 0.5 );
        lamda2.push_back( 0.0 );
        if( lamda2.size() > 0 ){
            generativeArgumentMaps["lambda2"] = lamda2;
            inputArgumentMap["lambda2"]       = lamda2[0];
        }   
#endif

        DoubleVector numberOfMolecules;
        //numberOfMolecules.push_back( 10 );
#if IMPLICIT_SOLVENT != TEST_NONBONDED
        numberOfMolecules.push_back( 100 );
#endif
        numberOfMolecules.push_back( 1000 );
        numberOfMolecules.push_back( 2000 );
        numberOfMolecules.push_back( 4000 );
        //numberOfMolecules.push_back( 8000 );
        if( numberOfMolecules.size() > 0 ){
            generativeArgumentMaps["numMolecules"] = numberOfMolecules;
            inputArgumentMap["numMolecules"]       = numberOfMolecules[0];
        }   

        DoubleVector nonbondedMethod;
        nonbondedMethod.push_back( NoCutoff_OpenMMTest );
        nonbondedMethod.push_back( CutoffNonPeriodic_OpenMMTest );
        nonbondedMethod.push_back( CutoffPeriodic_OpenMMTest );
#if IMPLICIT_SOLVENT == TEST_NONBONDED
        nonbondedMethod.push_back( Ewald_OpenMMTest );
        nonbondedMethod.push_back( PME_OpenMMTest );
#endif
        if( nonbondedMethod.size() > 0 ){
            generativeArgumentMaps["nonbondedMethod"] = nonbondedMethod;
            inputArgumentMap["nonbondedMethod"]       = nonbondedMethod[0];
        }

        DoubleVector platformId2s;
#if TEST_PLATFORM == TEST_OPENCL_PLATFORM
        platformId2s.push_back( OpenCL_OpenMMTest );
#elif TEST_PLATFORM == TEST_CUDA_PLATFORM
        platformId2s.push_back( Cuda_OpenMMTest );
#else
        platformId2s.push_back( Cuda_OpenMMTest );
#endif

        // check that required libs are available for platform to be tested
        // if unavailable, skip tests

        std::string libPrefix = "lib";
        std::string libSuffix = ".so";
#ifdef _MSC_VER
        libPrefix = "";
        libSuffix = ".dll";
#endif
#ifdef __APPLE__
        libSuffix = ".dylib";
#endif

        StringVector requiredLibs;
        for( unsigned int kk = 0; kk < platformId2s.size(); kk++ ){
            if( platformId2s[kk] == OpenCL_OpenMMTest ){
                requiredLibs.push_back( getLibName( libPrefix, libSuffix, "OpenMMOpenCL" ) );
            }
            if( platformId2s[kk] == Cuda_OpenMMTest ){
                requiredLibs.push_back( getLibName( libPrefix, libSuffix, "OpenMMCuda") );
#ifdef USE_SOFTCORE
                requiredLibs.push_back( getLibName( libPrefix, libSuffix, "OpenMMFreeEnergy" ) );
                requiredLibs.push_back( getLibName( libPrefix, libSuffix, "OpenMMFreeEnergyCuda" ) );
#endif
            }
        }

        // if TEST_PLATFORM is not set, then check that required libs are available
     
#if TEST_PLATFORM != TEST_OPENCL_PLATFORM && TEST_PLATFORM != TEST_CUDA_PLATFORM
        envVariableIsSet = checkRequiredLibsAreAvailable( requiredLibs, loadedLibs, log );
        if( envVariableIsSet == 0 && log ){
            (void) fprintf( log, "Aborting tests due to missing libs.\n" );
        }
#else

        // unit test path: force tests to run
        
        envVariableIsSet = 1;
#endif

        if( platformId2s.size() > 0 ){
            generativeArgumentMaps["platformId2"] = platformId2s;
            inputArgumentMap["platformId2"]       = platformId2s[0];
        }

        vectorOfMapStringToDouble.push_back( inputArgumentMap );
        generateInputArgumentMapsFromStringVectors( generativeArgumentMaps, vectorOfMapStringToDouble ); 

        // modify relative tolerance for large systems
        // case: Distance    433    669    1.5000001e+00 d2=  2.2500002e+00 w/ cutoff=1.500

        for( unsigned int kk = 0; kk < vectorOfMapStringToDouble.size(); kk++ ){
            int numMolecules            = 0;
            int numParticlesPerMolecule = 2;
            setIntFromMapStringToDouble( vectorOfMapStringToDouble[kk], "numMolecules",                 numMolecules );
            setIntFromMapStringToDouble( vectorOfMapStringToDouble[kk], "numParticlesPerMolecule",      numParticlesPerMolecule );
            if( numMolecules*numParticlesPerMolecule > 1000 ){
                vectorOfMapStringToDouble[kk]["relativeTolerance"]       = 6.0e-03;
            }
        }

        if( log ){
            MapStringToInt exclude;
            exclude["lambda1"]                 = 1;
            exclude["platformId1"]             = 1;
            exclude["platformId2"]             = 1;
            exclude["deviceId1"]               = 1;
            exclude["deviceId2"]               = 1;
            exclude["numParticlesPerMolecule"] = 1;
            std::stringstream outputStream;
            std::sort( vectorOfMapStringToDouble.begin(), vectorOfMapStringToDouble.end(), TestMapSortPredicate);
            StringVector printOrder;
            printOrder.push_back( "numMolecules" );
            printOrder.push_back( "nonbondedMethod" );
            printOrder.push_back( "lambda2" );
            printOrder.push_back( "boxSize" );
            for( unsigned int kk = 0; kk < vectorOfMapStringToDouble.size(); kk++ ){
                streamArgumentMapOneLine( vectorOfMapStringToDouble[kk], exclude, printOrder, kk, outputStream );
            }
            (void) fprintf( log, "Initial argument maps: %u\n%s", static_cast<unsigned int>(vectorOfMapStringToDouble.size()), outputStream.str().c_str() );
        }

        // run tests

        if( envVariableIsSet ){
            int wasException = 0;
            for( unsigned int kk = 0; kk < vectorOfMapStringToDouble.size() && wasException < 3; kk++ ){
                try {
                    runTests( vectorOfMapStringToDouble[kk], precision, log );
                } catch(const exception& e) {
                    std::stringstream msg;
#if IMPLICIT_SOLVENT == TEST_NONBONDED
                    msg << "Nonbonded";
#elif IMPLICIT_SOLVENT == TEST_OBC
                    msg << "GBSAOBC";
#elif IMPLICIT_SOLVENT == TEST_GBVI
                    msg << "GBVI";
#endif
                    int numMolecules            = 0;
                    int numParticlesPerMolecule = 0;
                    int nonbondedMethod         = 0;
                    setIntFromMapStringToDouble(    vectorOfMapStringToDouble[kk], "numMolecules",                 numMolecules );
                    setIntFromMapStringToDouble(    vectorOfMapStringToDouble[kk], "numParticlesPerMolecule",      numParticlesPerMolecule );
                    setIntFromMapStringToDouble(    vectorOfMapStringToDouble[kk], "nonbondedMethod",              nonbondedMethod);
                    msg << " test: system size=" << numMolecules*numParticlesPerMolecule << " nonbonded method=" << getNonbondedMethodName( nonbondedMethod );
                    msg << " exception: " << e.what() << endl;
                    // msg << "Note cases have been encountered for nonbonded methods with cutoffs where the error was due to particles being within 1.0e-05 of the cutoff." << endl;
                    cout << msg.str();
                    wasException += 1;
                }
            }
            if( wasException ){
                return 1;
            }
        }

    } catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
