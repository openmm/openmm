
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

#include "AmoebaReferenceForce.h"
#include "AmoebaReferenceMultipoleForce.h"
#include <sstream>
#include <algorithm>

#define AMOEBA_DEBUG

AmoebaReferenceMultipoleForce::AmoebaReferenceMultipoleForce( ) : _nonbondedMethod(NoCutoff) {
    initialize();
}

AmoebaReferenceMultipoleForce::AmoebaReferenceMultipoleForce( NonbondedMethod nonbondedMethod ) :
                                                  _nonbondedMethod(nonbondedMethod) {
    initialize();
}

void AmoebaReferenceMultipoleForce::initialize( void ){

    static const RealOpenMM zero                    = 0.0;
    static const RealOpenMM one                     = 1.0;
    static const RealOpenMM pointFour               = 0.4;
    static const RealOpenMM pointEight              = 0.8;

    _electric                                       = 138.9354558456;
    _dielectric                                     = one;

    unsigned int index                              = 0;
    _mScale[index++]                                = zero;
    _mScale[index++]                                = zero;
    _mScale[index++]                                = zero;
    _mScale[index++]                                = pointFour;
    _mScale[index++]                                = pointEight;

    index                                           = 0;
    _dScale[index++]                                = zero;
    _dScale[index++]                                = one;
    _dScale[index++]                                = one;
    _dScale[index++]                                = one;
    _dScale[index++]                                = one;

    index                                           = 0;
    _pScale[index++]                                = zero;
    _pScale[index++]                                = zero;
    _pScale[index++]                                = zero;
    _pScale[index++]                                = one;
    _pScale[index++]                                = one;

    index                                           = 0;
    _uScale[index++]                                = one;
    _uScale[index++]                                = one;
    _uScale[index++]                                = one;
    _uScale[index++]                                = one;
    _uScale[index++]                                = one;

    _mutualInducedDipoleConverged                  = 0;
    _mutualInducedDipoleIterations                 = 0;
    _maximumMutualInducedDipoleIterations          = 100;
    _mutualInducedDipoleEpsilon                    = 1.0e+50;
    _mutualInducedDipoleTargetEpsilon              = 1.0e-04;
    _polarSOR                                      = 0.70;
    _debye                                         = 480.33324;

}

AmoebaReferenceMultipoleForce::NonbondedMethod AmoebaReferenceMultipoleForce::getNonbondedMethod( void ) const {
    return _nonbondedMethod;
}

void AmoebaReferenceMultipoleForce::setNonbondedMethod( AmoebaReferenceMultipoleForce::NonbondedMethod nonbondedMethod ){
    _nonbondedMethod = nonbondedMethod;
}

int AmoebaReferenceMultipoleForce::getMutualInducedDipoleConverged( void ) const {
    return _mutualInducedDipoleConverged;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleConverged( int mutualInducedDipoleConverged ){
    _mutualInducedDipoleConverged = mutualInducedDipoleConverged;
}

int AmoebaReferenceMultipoleForce::getMutualInducedDipoleIterations( void ) const {
    return _mutualInducedDipoleIterations;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleIterations( int mutualInducedDipoleIterations ){
    _mutualInducedDipoleIterations = mutualInducedDipoleIterations;
}

RealOpenMM AmoebaReferenceMultipoleForce::getMutualInducedDipoleEpsilon( void ) const {
    return _mutualInducedDipoleEpsilon;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleEpsilon( RealOpenMM mutualInducedDipoleEpsilon ){
    _mutualInducedDipoleEpsilon = mutualInducedDipoleEpsilon;
}

int AmoebaReferenceMultipoleForce::getMaximumMutualInducedDipoleIterations( void ) const {
    return _maximumMutualInducedDipoleIterations;
}

void AmoebaReferenceMultipoleForce::setMaximumMutualInducedDipoleIterations( int maximumMutualInducedDipoleIterations ){
    _maximumMutualInducedDipoleIterations = maximumMutualInducedDipoleIterations;
}

RealOpenMM AmoebaReferenceMultipoleForce::getMutualInducedDipoleTargetEpsilon( void ) const {
    return _mutualInducedDipoleTargetEpsilon;
}

void AmoebaReferenceMultipoleForce::setMutualInducedDipoleTargetEpsilon( RealOpenMM mutualInducedDipoleTargetEpsilon ){
    _mutualInducedDipoleTargetEpsilon = mutualInducedDipoleTargetEpsilon;
}

void AmoebaReferenceMultipoleForce::getDelta( unsigned int particleI, unsigned int particleJ, RealOpenMM** particlePositions,
                                              RealOpenMM* delta ) const {

    delta[0]  = particlePositions[particleJ][0] - particlePositions[particleI][0];
    delta[1]  = particlePositions[particleJ][1] - particlePositions[particleI][1];
    delta[2]  = particlePositions[particleJ][2] - particlePositions[particleI][2];
}

void AmoebaReferenceMultipoleForce::getDelta( const MultipoleParticleData& particleI,
                                              const MultipoleParticleData& particleJ,
                                              RealOpenMM* delta ) const {

    delta[0]  = particleJ.position[0] - particleI.position[0];
    delta[1]  = particleJ.position[1] - particleI.position[1];
    delta[2]  = particleJ.position[2] - particleI.position[2];
}

void AmoebaReferenceMultipoleForce::loadArrayFromVector( unsigned int particleI, unsigned int offset, const std::vector<RealOpenMM>& vectorToCopy,
                                                         RealOpenMM* loadArray ) const {

    unsigned int vectorOffset = particleI*offset;
    for( unsigned ii = 0; ii < offset; ii++ ){
        loadArray[ii]  = vectorToCopy[vectorOffset+ii];
    }
}

void AmoebaReferenceMultipoleForce::logRealOpenMMVectors( const std::string& header, const VectorOfRealOpenMMVectors& printVector,
                                                          FILE* log, unsigned int itemsPerVector, int maxPrint ) const {

    (void) fprintf( log, "%s", header.c_str() );
    for( unsigned int ii = 0; ii < printVector[0].size()/itemsPerVector; ii++ ){
        (void) fprintf( log, "%5u ", ii );
        for( unsigned int jj = 0; jj < printVector.size(); jj++ ){
            if( itemsPerVector > 1 ){
                (void) fprintf( log, "[" );
            }
            for( unsigned int kk = 0; kk < itemsPerVector; kk++ ){
                (void) fprintf( log, "%15.7e ", printVector[jj][ii*itemsPerVector+kk] );
            }
            if( itemsPerVector > 1 ){
                (void) fprintf( log, "] " );
            }
        }
        (void) fprintf( log, "\n" );
    }
}

void AmoebaReferenceMultipoleForce::logParticleData( const std::string& header, const std::vector<MultipoleParticleData>& particleData,
                                                     unsigned int printFlag, FILE* log, unsigned int maxPrint ) const {

    (void) fprintf( log, "%s", header.c_str() );

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        (void) fprintf( log, "%5u ", particleData[ii].particleIndex );
        if( printFlag & (1 << PARTICLE_POSITION) ){
            (void) fprintf( log, "x[%15.7e %15.7e %15.7e] ", particleData[ii].position[0], particleData[ii].position[1], particleData[ii].position[2] );
        }
        if( printFlag & (1 << PARTICLE_DIPOLE) ){
            (void) fprintf( log, "d[%15.7e %15.7e %15.7e] ", particleData[ii].dipole[0], particleData[ii].dipole[1], particleData[ii].dipole[2] );
        }
        if( printFlag & (1 << PARTICLE_QUADRUPOLE) ){
            (void) fprintf( log, "q[%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e] ",
                                        particleData[ii].quadrupole[QXX], particleData[ii].quadrupole[QXY], particleData[ii].quadrupole[QXZ],
                                        particleData[ii].quadrupole[QXY], particleData[ii].quadrupole[QYY], particleData[ii].quadrupole[QYZ],
                                        particleData[ii].quadrupole[QXZ], particleData[ii].quadrupole[QYZ], particleData[ii].quadrupole[QZZ] );
        }
        if( printFlag & (1 << PARTICLE_FIELD) ){
            (void) fprintf( log, "f[%15.7e %15.7e %15.7e] ", particleData[ii].field[0], particleData[ii].field[1], particleData[ii].field[2] );
        }
        if( printFlag & (1 << PARTICLE_FIELD_POLAR) ){
            (void) fprintf( log, "fp[%15.7e %15.7e %15.7e] ", particleData[ii].fieldPolar[0], particleData[ii].fieldPolar[1], particleData[ii].fieldPolar[2] );
        }
        if( printFlag & (1 << PARTICLE_INDUCED_DIPOLE) ){
            (void) fprintf( log, "i[%15.7e %15.7e %15.7e] ", particleData[ii].inducedDipole[0], particleData[ii].inducedDipole[1], particleData[ii].inducedDipole[2] );
        }
        if( printFlag & (1 << PARTICLE_INDUCED_DIPOLE_POLAR) ){
            (void) fprintf( log, "ip[%15.7e %15.7e %15.7e] ", particleData[ii].inducedDipolePolar[0], particleData[ii].inducedDipolePolar[1], particleData[ii].inducedDipolePolar[2] );
        }
        (void) fprintf( log, "\n" );
        if( (maxPrint > 0) && (ii == maxPrint) && (2*maxPrint < particleData.size()) ){
            ii = particleData.size() - maxPrint;
        }
    }
}

/* 
    Show scaling factors for a particle 
 */ 
void AmoebaReferenceMultipoleForce::showScaleMapForParticle( unsigned int particleI, FILE* log ) const {

    (void) fprintf( log, "Scale map particle %5u maxIndex=%u\n", particleI, _maxScaleIndex[particleI] );

    std::string scaleNames[LAST_SCALE_TYPE_INDEX] = { "D", "P", "M" }; 
    for( unsigned int ii = 0; ii < _scaleMaps[particleI].size(); ii++ ){
        MapIntRealOpenMM scaleMap = _scaleMaps[particleI][ii];
        (void) fprintf( log, "  %s scale ", scaleNames[ii].c_str() );
        for( MapIntRealOpenMMCI jj = scaleMap.begin(); jj != scaleMap.end(); jj++ ){
            //if( jj->first > particleI && jj->second < 1.0 )
            if( jj->second < 1.0 )
            (void) fprintf( log, "%4d=%5.2f ", jj->first, jj->second );
        }
        (void) fprintf( log, "\n" );
    }
    (void) fprintf( log, "\n" );
    (void) fflush( log );
}

void AmoebaReferenceMultipoleForce::setupScaleMaps( const std::vector< std::vector< std::vector<int> > >& multipoleParticleCovalentInfo ){

    /* Setup for scaling maps:

           _scaleMaps[particleIndex][ScaleType] = map, where map[covalentIndex] = scaleFactor 
           _maxScaleIndex[particleIndex]        = max covalent index for particleIndex
 
           multipoleParticleCovalentInfo[ii][jj], jj =0,1,2,3 contains covalent indices (c12, c13, c14, c15)
           multipoleParticleCovalentInfo[ii][jj], jj =4,5,6,7 contains covalent indices (p11, p12, p13, p14)

           only including covalent particles w/ index >= ii
     */

    _scaleMaps.resize( multipoleParticleCovalentInfo.size() );
    _maxScaleIndex.resize( multipoleParticleCovalentInfo.size() );

    for( unsigned int ii = 0; ii < multipoleParticleCovalentInfo.size(); ii++ ){

        _scaleMaps[ii].resize(LAST_SCALE_TYPE_INDEX);
        _maxScaleIndex[ii] = 0;
        const std::vector< std::vector<int> >& covalentInfo = multipoleParticleCovalentInfo[ii];
        const std::vector<int> covalentListP11              = covalentInfo[AmoebaMultipoleForce::PolarizationCovalent11];

        // pScale & mScale

        for( unsigned jj = 0; jj < AmoebaMultipoleForce::PolarizationCovalent11; jj++ ){
            const std::vector<int> covalentList    = covalentInfo[jj];
            for( unsigned int kk = 0; kk < covalentList.size(); kk++ ){
                unsigned int covalentIndex             = static_cast<unsigned int>(covalentList[kk]);
                if( covalentIndex < ii )continue;

                // handle 0.5 factor for p14

                int hit = 0;
                if( jj == AmoebaMultipoleForce::Covalent14 ){
                    for( unsigned int mm = 0; mm < covalentListP11.size() && hit == 0; mm++ ){
                        if( covalentListP11[mm]  == covalentIndex ){
                            hit = 1;
                        }
                    }
                } 
                _scaleMaps[ii][P_SCALE][covalentIndex] = hit ? 0.5*_pScale[jj+1] : _pScale[jj+1];
                _scaleMaps[ii][M_SCALE][covalentIndex] = _mScale[jj+1];
                _maxScaleIndex[ii]                     = _maxScaleIndex[ii] < covalentIndex ? covalentIndex : _maxScaleIndex[ii];
            }
        }

        // dScale & uScale

        for( unsigned jj = AmoebaMultipoleForce::PolarizationCovalent11; jj < covalentInfo.size(); jj++ ){
            const std::vector<int> covalentList = covalentInfo[jj];
            for( unsigned int kk = 0; kk < covalentList.size(); kk++ ){
                unsigned int covalentIndex             = static_cast<unsigned int>(covalentList[kk]);
                if( covalentIndex < ii )continue;
                _scaleMaps[ii][D_SCALE][covalentIndex] = _dScale[jj-4];
                _scaleMaps[ii][U_SCALE][covalentIndex] = _uScale[jj-4];
                _maxScaleIndex[ii]                     = _maxScaleIndex[ii] < covalentIndex ? covalentIndex : _maxScaleIndex[ii];
            }
        }
    }
    //showScaleMapForParticle( 3, stderr );
    //showScaleMapForParticle( 10, stderr );
}

RealOpenMM AmoebaReferenceMultipoleForce::getScaleFactor( unsigned int particleI, unsigned int particleJ, ScaleType scaleType ) const {

    static const RealOpenMM one         = 1.0;

    MapIntRealOpenMM  scaleMap   = _scaleMaps[particleI][scaleType];
    MapIntRealOpenMMCI isPresent = scaleMap.find( particleJ );
    if( isPresent != scaleMap.end() ){
        return isPresent->second;
    } else {
        return one;
    }
}

void AmoebaReferenceMultipoleForce::getDScaleAndPScale( unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale ) const {

    dScale = getScaleFactor( particleI, particleJ, D_SCALE );
    pScale = getScaleFactor( particleI, particleJ, P_SCALE );
}

void AmoebaReferenceMultipoleForce::getScaleFactors( unsigned int particleI, unsigned int particleJ, RealOpenMM* scaleFactors ) const {

    scaleFactors[D_SCALE] = getScaleFactor( particleI, particleJ, D_SCALE );
    scaleFactors[P_SCALE] = getScaleFactor( particleI, particleJ, P_SCALE );
    scaleFactors[M_SCALE] = getScaleFactor( particleI, particleJ, M_SCALE );
    scaleFactors[U_SCALE] = getScaleFactor( particleI, particleJ, U_SCALE );
}

void AmoebaReferenceMultipoleForce::loadParticleData( RealOpenMM** particlePositions, 
                                                      const std::vector<RealOpenMM>& charges,
                                                      const std::vector<RealOpenMM>& dipoles,
                                                      const std::vector<RealOpenMM>& quadrupoles,
                                                      const std::vector<RealOpenMM>& tholes,
                                                      const std::vector<RealOpenMM>& dampingFactors,
                                                      const std::vector<RealOpenMM>& polarity,
                                                      std::vector<MultipoleParticleData>& particleData ) const {
   
    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM zero        = 0.0;

    // ---------------------------------------------------------------------------------------

    for( unsigned int ii = 0; ii < charges.size(); ii++ ){

        particleData[ii].particleIndex        = ii;

        particleData[ii].position[0]          = particlePositions[ii][0];
        particleData[ii].position[1]          = particlePositions[ii][1];
        particleData[ii].position[2]          = particlePositions[ii][2];

        particleData[ii].charge               = charges[ii];

        particleData[ii].dipole[0]            = dipoles[3*ii+0];
        particleData[ii].dipole[1]            = dipoles[3*ii+1];
        particleData[ii].dipole[2]            = dipoles[3*ii+2];

        particleData[ii].quadrupole[QXX]      = quadrupoles[9*ii+0];
        particleData[ii].quadrupole[QXY]      = quadrupoles[9*ii+1];
        particleData[ii].quadrupole[QXZ]      = quadrupoles[9*ii+2];
        particleData[ii].quadrupole[QYY]      = quadrupoles[9*ii+4];
        particleData[ii].quadrupole[QYZ]      = quadrupoles[9*ii+5];
        particleData[ii].quadrupole[QZZ]      = quadrupoles[9*ii+8];

        particleData[ii].thole                = tholes[ii];
        particleData[ii].dampingFactor        = dampingFactors[ii];
        particleData[ii].polarity             = polarity[ii];

        particleData[ii].field[0]             = zero;
        particleData[ii].field[1]             = zero;
        particleData[ii].field[2]             = zero;

        particleData[ii].fieldPolar[0]        = zero;
        particleData[ii].fieldPolar[1]        = zero;
        particleData[ii].fieldPolar[2]        = zero;
    }
}

void AmoebaReferenceMultipoleForce::applyRotationMatrix(       MultipoleParticleData& particleI,
                                                         const MultipoleParticleData& particleZ,
                                                         const MultipoleParticleData& particleX, int axisType ) const {

    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM zero        = 0.0;
    static const RealOpenMM one         = 1.0;
    static const RealOpenMM two         = 2.0;

    static const std::string methodName = "AmoebaReferenceMultipoleForce::applyRotationMatrix";
    static const int Y                  = 1;
    static const int X                  = 0;
    static const int Z                  = 2;
 
    RealOpenMM* vector[3];
    RealOpenMM  rotationMatrix[3][3];
 
    // ---------------------------------------------------------------------------------------

    // handle case where rotation matrix is identity (e.g. single ion)

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom
 
    vector[X]                       = rotationMatrix[0];
    vector[Y]                       = rotationMatrix[1];
    vector[Z]                       = rotationMatrix[2];

    getDelta( particleI, particleZ, vector[Z] );
    getDelta( particleI, particleX, vector[X] );

    AmoebaReferenceForce::normalizeVector3( vector[Z] );
 
    // branch based on axis type
 
    if( axisType == AmoebaMultipoleForce::Bisector ){
 
        // bisector
  
        // dx = dx1 + dx2 (in Tinker code)
       
        AmoebaReferenceForce::normalizeVector3( vector[X] );

        vector[Z][0] += vector[X][0];
        vector[Z][1] += vector[X][1];
        vector[Z][2] += vector[X][2];
       
        AmoebaReferenceForce::normalizeVector3( vector[Z] );
    }
 
    RealOpenMM dot      = vector[Z][0]*vector[X][0] + vector[Z][1]*vector[X][1] + vector[Z][2]*vector[X][2];
    vector[X][0]       -= dot*vector[Z][0];
    vector[X][1]       -= dot*vector[Z][1];
    vector[X][2]       -= dot*vector[Z][2];

    AmoebaReferenceForce::normalizeVector3( vector[X] );
    AmoebaReferenceForce::getCrossProduct( vector[Z], vector[X], vector[Y] );
 
    RealOpenMM labDipole[3];
    for( int ii = 0; ii < 3; ii++ ){
        labDipole[ii] = particleI.dipole[0]*rotationMatrix[0][ii];
        for( int jj = 1; jj < 3; jj++ ){
            labDipole[ii] += particleI.dipole[jj]*rotationMatrix[jj][ii];
        }
    }
    particleI.dipole[0] = labDipole[0];
    particleI.dipole[1] = labDipole[1];
    particleI.dipole[2] = labDipole[2];
 
    RealOpenMM mPole[3][3];
    RealOpenMM rPole[3][3] = { { zero, zero, zero },
                               { zero, zero, zero },
                               { zero, zero, zero } };

    mPole[0][0] = particleI.quadrupole[QXX];
    mPole[0][1] = particleI.quadrupole[QXY];
    mPole[0][2] = particleI.quadrupole[QXZ];

    mPole[1][0] = particleI.quadrupole[QXY];
    mPole[1][1] = particleI.quadrupole[QYY];
    mPole[1][2] = particleI.quadrupole[QYZ];

    mPole[2][0] = particleI.quadrupole[QXZ];
    mPole[2][1] = particleI.quadrupole[QYZ];
    mPole[2][2] = particleI.quadrupole[QZZ];
 
    for( int ii = 0; ii < 3; ii++ ){
       for( int jj = ii; jj < 3; jj++ ){
          for( int kk = 0; kk < 3; kk++ ){
             for( int mm = 0; mm < 3; mm++ ){
                 rPole[ii][jj] += rotationMatrix[kk][ii]*rotationMatrix[mm][jj]*mPole[kk][mm];
             }
          }
       }
    }
 
    particleI.quadrupole[QXX] = rPole[0][0];
    particleI.quadrupole[QXY] = rPole[0][1];
    particleI.quadrupole[QXZ] = rPole[0][2];

    particleI.quadrupole[QYY] = rPole[1][1];
    particleI.quadrupole[QYZ] = rPole[1][2];
    particleI.quadrupole[QZZ] = rPole[2][2];

    return;

}

void AmoebaReferenceMultipoleForce::getAndScaleInverseRs( RealOpenMM dampI, RealOpenMM dampJ,
                                                          RealOpenMM tholeI, RealOpenMM tholeJ,
                                                          RealOpenMM r, std::vector<RealOpenMM>& rrI ) const {

    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM zero        = 0.0;
    static const RealOpenMM one         = 1.0;
    static const RealOpenMM two         = 2.0;
    static const RealOpenMM three       = 3.0;
    static const RealOpenMM five        = 5.0;
    static const RealOpenMM fifty       = 50.0;

    static const std::string methodName = "AmoebaReferenceMultipoleForce::scaleInverseRs";

    // ---------------------------------------------------------------------------------------
 
    RealOpenMM rI             =  one/r;
    RealOpenMM r2I            =  rI*rI;
 
    rrI[0]                    = rI*r2I;
    RealOpenMM constantFactor = 3.0;
    for( unsigned int ii  = 1; ii < rrI.size(); ii++ ){ 
       rrI[ii]         = constantFactor*rrI[ii-1]*r2I;
       constantFactor += 2.0;
    }
 
    RealOpenMM damp      = dampI*dampJ;
    if( damp != zero ){
        RealOpenMM pgamma    = tholeI < tholeJ ? tholeI : tholeJ;
        RealOpenMM ratio     = (r/damp);
                   ratio     = ratio*ratio*ratio;
                   damp      = -pgamma*ratio;

        if( damp > -fifty ){ 
            RealOpenMM dampExp   = EXP( damp );

            rrI[0]              *= one - dampExp;
            rrI[1]              *= one - ( one - damp )*dampExp;
            if( rrI.size() > 2 ){
                rrI[2]          *= one - ( one - damp + (0.6*damp*damp))*dampExp;
            }
       }
   }
}

void AmoebaReferenceMultipoleForce::calculateFixedEFieldPairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                                                 RealOpenMM dScale, RealOpenMM pScale ) const {

    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM zero        = 0.0;
    static const RealOpenMM one         = 1.0;
    static const RealOpenMM two         = 2.0;

    //static const std::string methodName = "AmoebaReferenceMultipoleForce::calculateFixedEFieldPairIxn";

    // ---------------------------------------------------------------------------------------
 
    // get deltaR, R2, and R between 2 atoms
 
    RealOpenMM deltaR[3];
    getDelta( particleI, particleJ, deltaR );

    RealOpenMM r      = AmoebaReferenceForce::getNorm3( deltaR );
    std::vector<RealOpenMM> rrI(3);
 
    // get scaling factors, if needed
  
    getAndScaleInverseRs( particleI.dampingFactor, particleJ.dampingFactor, particleI.thole, particleJ.thole, r, rrI );

    RealOpenMM rr3    = rrI[0];
    RealOpenMM rr5    = rrI[1];
    RealOpenMM rr7    = rrI[2];
    RealOpenMM rr5_2  = rr5*two;

    // field at particle I due multipoles at particle J

    RealOpenMM qDotDelta[3];
    qDotDelta[0]                  = deltaR[0]*particleJ.quadrupole[QXX] + deltaR[1]*particleJ.quadrupole[QXY] + deltaR[2]*particleJ.quadrupole[QXZ];
    qDotDelta[1]                  = deltaR[0]*particleJ.quadrupole[QXY] + deltaR[1]*particleJ.quadrupole[QYY] + deltaR[2]*particleJ.quadrupole[QYZ];
    qDotDelta[2]                  = deltaR[0]*particleJ.quadrupole[QXZ] + deltaR[1]*particleJ.quadrupole[QYZ] + deltaR[2]*particleJ.quadrupole[QZZ];

    RealOpenMM dipoleDelta        = AmoebaReferenceForce::getDotProduct3( particleJ.dipole, deltaR ); 
    RealOpenMM qdpoleDelta        = AmoebaReferenceForce::getDotProduct3( qDotDelta, deltaR ); 
    RealOpenMM factor             = rr3*particleJ.charge - rr5*dipoleDelta + rr7*qdpoleDelta;
 
    RealOpenMM field[3];
    field[0]                      = deltaR[0]*factor + rr3*particleJ.dipole[0] - rr5_2*qDotDelta[0];
    field[1]                      = deltaR[1]*factor + rr3*particleJ.dipole[1] - rr5_2*qDotDelta[1];
    field[2]                      = deltaR[2]*factor + rr3*particleJ.dipole[2] - rr5_2*qDotDelta[2];

    for( unsigned int ii = 0; ii < 3; ii++ ){
        particleI.field[ii]       -= dScale*field[ii];
        particleI.fieldPolar[ii]  -= pScale*field[ii];
    }
 
    // field at particle J due multipoles at particle I

    qDotDelta[0]                  = deltaR[0]*particleI.quadrupole[QXX] + deltaR[1]*particleI.quadrupole[QXY] + deltaR[2]*particleI.quadrupole[QXZ];
    qDotDelta[1]                  = deltaR[0]*particleI.quadrupole[QXY] + deltaR[1]*particleI.quadrupole[QYY] + deltaR[2]*particleI.quadrupole[QYZ];
    qDotDelta[2]                  = deltaR[0]*particleI.quadrupole[QXZ] + deltaR[1]*particleI.quadrupole[QYZ] + deltaR[2]*particleI.quadrupole[QZZ];

    dipoleDelta                   = AmoebaReferenceForce::getDotProduct3( particleI.dipole, deltaR ); 
    qdpoleDelta                   = AmoebaReferenceForce::getDotProduct3( qDotDelta, deltaR ); 
    factor                        = rr3*particleI.charge + rr5*dipoleDelta + rr7*qdpoleDelta;
 
    field[0]                      = deltaR[0]*factor - rr3*particleI.dipole[0] - rr5_2*qDotDelta[0];
    field[1]                      = deltaR[1]*factor - rr3*particleI.dipole[1] - rr5_2*qDotDelta[1];
    field[2]                      = deltaR[2]*factor - rr3*particleI.dipole[2] - rr5_2*qDotDelta[2];

    for( unsigned int ii = 0; ii < 3; ii++ ){
        particleJ.field[ii]       += dScale*field[ii];
        particleJ.fieldPolar[ii]  += pScale*field[ii];
    }
 
    return;
}

void AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxn( MultipoleParticleData& particleI, 
                                                                   MultipoleParticleData& particleJ,
                                                                   std::vector<RealOpenMM>& field,
                                                                   std::vector<RealOpenMM>& fieldPolar ) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM zero        = 0.0;
    static const RealOpenMM one         = 1.0;
    static const RealOpenMM two         = 2.0;

    static const std::string methodName = "AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxn";

    // ---------------------------------------------------------------------------------------

   // get deltaR, R2, and R between 2 atoms

    RealOpenMM deltaR[3];
    getDelta( particleI, particleJ, deltaR );

    RealOpenMM r         =  AmoebaReferenceForce::getNorm3( deltaR );
    std::vector<RealOpenMM> rrI(2);
  
    getAndScaleInverseRs( particleI.dampingFactor, particleJ.dampingFactor,
                          particleI.thole, particleJ.thole, r, rrI );
 
    RealOpenMM rr3                  = -rrI[0];
    RealOpenMM rr5                  =  rrI[1];
 
    unsigned int indexOffset        = 3*particleI.particleIndex;
    RealOpenMM dDotDelta            = rr5*AmoebaReferenceForce::getDotProduct3( particleJ.inducedDipole, deltaR );
    field[indexOffset+0]           += rr3*particleJ.inducedDipole[0] + dDotDelta*deltaR[0];
    field[indexOffset+1]           += rr3*particleJ.inducedDipole[1] + dDotDelta*deltaR[1];
    field[indexOffset+2]           += rr3*particleJ.inducedDipole[2] + dDotDelta*deltaR[2];

    dDotDelta                       = rr5*AmoebaReferenceForce::getDotProduct3( particleJ.inducedDipolePolar, deltaR );
    fieldPolar[indexOffset+0]      += rr3*particleJ.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fieldPolar[indexOffset+1]      += rr3*particleJ.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fieldPolar[indexOffset+2]      += rr3*particleJ.inducedDipolePolar[2] + dDotDelta*deltaR[2];
 
    dDotDelta                       = rr5*AmoebaReferenceForce::getDotProduct3( particleI.inducedDipole, deltaR );
    indexOffset                     = 3*particleJ.particleIndex;
    field[indexOffset+0]           += rr3*particleI.inducedDipole[0] + dDotDelta*deltaR[0];
    field[indexOffset+1]           += rr3*particleI.inducedDipole[1] + dDotDelta*deltaR[1];
    field[indexOffset+2]           += rr3*particleI.inducedDipole[2] + dDotDelta*deltaR[2];
 
    dDotDelta                       = rr5*AmoebaReferenceForce::getDotProduct3( particleI.inducedDipolePolar, deltaR );
    fieldPolar[indexOffset+0]      += rr3*particleI.inducedDipolePolar[0] + dDotDelta*deltaR[0];
    fieldPolar[indexOffset+1]      += rr3*particleI.inducedDipolePolar[1] + dDotDelta*deltaR[1];
    fieldPolar[indexOffset+2]      += rr3*particleI.inducedDipolePolar[2] + dDotDelta*deltaR[2];
 
    return;

}

void AmoebaReferenceMultipoleForce::calculateInducedDipoleField( std::vector<MultipoleParticleData>& particleData,
                                                                 std::vector<RealOpenMM>& field,
                                                                 std::vector<RealOpenMM>& fieldPolar ) const {


    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceMultipoleForce::calculateInducedDipoleField";

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;
    static const RealOpenMM  two           = 2.0;

	 FILE* log                              = stderr;

    // ---------------------------------------------------------------------------------------

    std::fill( field.begin(),      field.end(),      zero );
    std::fill( fieldPolar.begin(), fieldPolar.end(), zero );

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){
            calculateInducedDipolePairIxn( particleData[ii], particleData[jj], field, fieldPolar);
        }
    }

#ifdef AMOEBA_DEBUG
    if( 0 && log ){
        std::string header    = "Induced fields\n";
        VectorOfRealOpenMMVectors printVector(2);
        printVector[0]        = field;
        printVector[1]        = fieldPolar;
        int maxPrint          = 10;
        logRealOpenMMVectors( header, printVector, log, 3, maxPrint );
        exit(0);
    }
#endif

    return;
}

void AmoebaReferenceMultipoleForce::updateInducedDipole(  MultipoleParticleData& particleI,
                                                         const std::vector<RealOpenMM>& field,
                                                         const std::vector<RealOpenMM>& fieldPolar,
                                                         RealOpenMM& epsilonD, RealOpenMM& epsilonP ) const {

    unsigned int particleOffset = particleI.particleIndex*3;
    for( unsigned int ii = 0; ii < 3; ii++ ){

        RealOpenMM oldValue               = particleI.inducedDipole[ii];
        RealOpenMM newValue               = particleI.field[ii] + particleI.polarity*field[particleOffset+ii];

        RealOpenMM delta                  = newValue - oldValue;
                   newValue               = oldValue + _polarSOR*delta;
        particleI.inducedDipole[ii]       = newValue;

        // polarSOR*polarSOR factor included in final epsilon calculation

        epsilonD                         += delta*delta;

        oldValue                          = particleI.inducedDipolePolar[ii];
        newValue                          = particleI.fieldPolar[ii] + particleI.polarity*fieldPolar[particleOffset+ii];

        delta                             = newValue - oldValue;
        newValue                          = oldValue + _polarSOR*delta;
        particleI.inducedDipolePolar[ii]  = newValue;

        // polarSOR*polarSOR factor included in final epsilon calculation

        epsilonP                         += delta*delta;
    }
}

void AmoebaReferenceMultipoleForce::calculateNoCutoffInducedDipoles( std::vector<MultipoleParticleData>& particleData ){

    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceMultipoleForce::calculateNoCutoffInducedDipoles";

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;
    static const RealOpenMM  two           = 2.0;

	 FILE* log                              = stderr;

    // ---------------------------------------------------------------------------------------

    // initialize inducedDipole to fixedE_Field

    unsigned int numParticles = particleData.size();
    for( unsigned int ii = 0; ii < numParticles; ii++ ){

        particleData[ii].field[0]                 *= particleData[ii].polarity; 
        particleData[ii].field[1]                 *= particleData[ii].polarity; 
        particleData[ii].field[2]                 *= particleData[ii].polarity; 

        particleData[ii].fieldPolar[0]            *= particleData[ii].polarity; 
        particleData[ii].fieldPolar[1]            *= particleData[ii].polarity; 
        particleData[ii].fieldPolar[2]            *= particleData[ii].polarity; 

        particleData[ii].inducedDipole[0]          = particleData[ii].field[0]; 
        particleData[ii].inducedDipole[1]          = particleData[ii].field[1]; 
        particleData[ii].inducedDipole[2]          = particleData[ii].field[2]; 

        particleData[ii].inducedDipolePolar[0]     = particleData[ii].fieldPolar[0]; 
        particleData[ii].inducedDipolePolar[1]     = particleData[ii].fieldPolar[1]; 
        particleData[ii].inducedDipolePolar[2]     = particleData[ii].fieldPolar[2]; 
    }

    std::vector<RealOpenMM> field( numParticles*3 );
    std::vector<RealOpenMM> fieldPolar( numParticles*3 );

    bool done                 = false;
    setMutualInducedDipoleConverged( false );
    unsigned int iteration    = 0;
    RealOpenMM currentEpsilon =  1.0e+50;
    while( !done ){

        calculateInducedDipoleField( particleData, field, fieldPolar );   

        RealOpenMM epsilonDirect    = zero;
        RealOpenMM epsilonPolar     = zero;

        for( unsigned int ii = 0; ii < numParticles; ii++ ){
            updateInducedDipole( particleData[ii], field, fieldPolar, epsilonDirect, epsilonPolar );
        }    

        RealOpenMM epsilon = epsilonDirect > epsilonPolar ? epsilonDirect : epsilonPolar;
                   epsilon = _polarSOR*_debye*SQRT( epsilon/( static_cast<RealOpenMM>(numParticles) ) );

fprintf( stderr, "MIb %3u eps=%15.7e %15.7e %15.7e\n", iteration, epsilonDirect, epsilonPolar, epsilon );

        if( epsilon < getMutualInducedDipoleTargetEpsilon() ){
            setMutualInducedDipoleConverged( true );
            done = true;
        }
        if( currentEpsilon < epsilon || static_cast<int>(iteration) >= getMaximumMutualInducedDipoleIterations() ){
            done = true;
        }
        currentEpsilon = epsilon;
        iteration++;
    }
    setMutualInducedDipoleEpsilon( currentEpsilon );
    setMutualInducedDipoleIterations( iteration );

#ifdef AMOEBA_DEBUG
    if( log ){
        std::stringstream header;
        header << "MutualInducedDipoles:";
        header << " converged=" << getMutualInducedDipoleConverged();
        header << " iter="      << getMutualInducedDipoleIterations();
        header << " eps="       << getMutualInducedDipoleEpsilon();
        header << "\n";
        unsigned int printFlag    =  (1 << PARTICLE_INDUCED_DIPOLE ) | (1 << PARTICLE_INDUCED_DIPOLE_POLAR);
        int maxPrint              = 10;
        logParticleData( header.str(), particleData, printFlag, log, maxPrint );
    }
#endif

    return;
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateNoCutoffElectrostaticPairIxn( const MultipoleParticleData& particleI,
                                                                                 const MultipoleParticleData& particleK,
                                                                                 RealOpenMM* scalingFactors, RealOpenMM** forces,
                                                                                 std::vector<Vec3>& torque ) const { 

/*
off2
*/
    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceMultipoleForce::calculateNoCutoffElectrostatic";

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;
    static const RealOpenMM  two           = 2.0;

    FILE* log                              = stderr;

    // ---------------------------------------------------------------------------------------

    RealOpenMM e,ei;
    RealOpenMM damp,expdamp;
    RealOpenMM pdi,pti,pgamma;
    RealOpenMM scale3,scale3i;
    RealOpenMM scale5,scale5i;
    RealOpenMM scale7,scale7i;
    RealOpenMM temp3,temp5,temp7;
    RealOpenMM psc3,psc5,psc7;
    RealOpenMM dsc3,dsc5,dsc7;
    RealOpenMM xr,yr,zr;
    //RealOpenMM xix,yix,zix;
    //RealOpenMM xiy,yiy,ziy;
    //RealOpenMM xiz,yiz,ziz;
    //RealOpenMM xkx,ykx,zkx;
    //RealOpenMM xky,yky,zky;
    //RealOpenMM xkz,ykz,zkz;
    RealOpenMM r,r2,rr1,rr3;
    RealOpenMM rr5,rr7,rr9,rr11;
    //RealOpenMM vxx,vyy,vzz;
    //RealOpenMM vyx,vzx,vzy;
    RealOpenMM qi[10];
    RealOpenMM qk[10];
/*
    RealOpenMM frcxi[4],frcxk[4];
    RealOpenMM frcyi[4],frcyk[4];
    RealOpenMM frczi[4],frczk[4];
*/
    RealOpenMM fridmp[4],findmp[4];
    RealOpenMM ftm2[4],ftm2i[4];
    RealOpenMM ttm2[4],ttm3[4];
    RealOpenMM ttm2i[4],ttm3i[4];
    RealOpenMM dixdk[4];
    RealOpenMM dixuk[4],dkxui[4];
    RealOpenMM dixukp[4],dkxuip[4];
    RealOpenMM uixqkr[4],ukxqir[4];
    RealOpenMM uixqkrp[4],ukxqirp[4];
    RealOpenMM qiuk[4],qkui[4];
    RealOpenMM qiukp[4],qkuip[4];
    RealOpenMM rxqiuk[4],rxqkui[4];
    RealOpenMM rxqiukp[4],rxqkuip[4];
    RealOpenMM qidk[4],qkdi[4];
    RealOpenMM qir[4],qkr[4];
    RealOpenMM qiqkr[4],qkqir[4];
    RealOpenMM qixqk[4],rxqir[4];
    RealOpenMM dixr[4],dkxr[4];
    RealOpenMM dixqkr[4],dkxqir[4];
    RealOpenMM rxqkr[4],qkrxqir[4];
    RealOpenMM rxqikr[4],rxqkir[4];
    RealOpenMM rxqidk[4],rxqkdi[4];
    RealOpenMM ddsc3[4],ddsc5[4];
    RealOpenMM ddsc7[4];
    RealOpenMM gl[10],gli[8],glip[8];
    RealOpenMM sc[11],sci[9],scip[9];
    RealOpenMM gf[8],gfi[7],gti[7];
 
    RealOpenMM delta[3];
    getDelta( particleI, particleK, delta );
    xr = delta[0];
    yr = delta[1];
    zr = delta[2];
//  call image (xr,yr,zr);
    r2 = AmoebaReferenceForce::getNormSquared3( delta );
//    if( r2 > off2 )return zero;

    // set conversion factor, cutoff and switching coefficients

    RealOpenMM f = _electric /_dielectric;

    // call switch ('MPOLE');

    // set scale factors for permanent multipole and induced terms

    pdi   = particleI.dampingFactor;
    pti   = particleI.thole;
    qi[0] = particleI.quadrupole[QXX];
    qi[1] = particleI.quadrupole[QXY];
    qi[2] = particleI.quadrupole[QXZ];
    qi[3] = particleI.quadrupole[QXY];
    qi[4] = particleI.quadrupole[QYY];
    qi[5] = particleI.quadrupole[QYZ];
    qi[6] = particleI.quadrupole[QXZ];
    qi[7] = particleI.quadrupole[QYZ];
    qi[8] = particleI.quadrupole[QZZ];

    // set interaction scaling coefficients for connected atoms;

/*
            kz = zaxis(k);
            kx = xaxis(k);

            usek = (use(kk) .or. use(kz) .or. use(kx));
            proceed = .true.;
            if( use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0);
            if( .not. use_intra)  proceed = .true.;
            if( proceed)  proceed = (usei .or. usek);
            if( .not. proceed)  goto 10;
*/
    qk[0] = particleK.quadrupole[QXX];
    qk[1] = particleK.quadrupole[QXY];
    qk[2] = particleK.quadrupole[QXZ];
    qk[3] = particleK.quadrupole[QXY];
    qk[4] = particleK.quadrupole[QYY];
    qk[5] = particleK.quadrupole[QYZ];
    qk[6] = particleK.quadrupole[QXZ];
    qk[7] = particleK.quadrupole[QYZ];
    qk[8] = particleK.quadrupole[QZZ];
 
    // apply Thole polarization damping to scale factors

    r = sqrt(r2);
    rr1 = 1.0 / r;
    rr3 = rr1 / r2;
    rr5 = 3.0 * rr3 / r2;
    rr7 = 5.0 * rr5 / r2;
    rr9 = 7.0 * rr7 / r2;
    rr11 = 9.0 * rr9 / r2;
    scale3 = 1.0;
    scale5 = 1.0;
    scale7 = 1.0;
    for( unsigned int jj = 0; jj < 3; jj++ ){
        ddsc3[jj] = 0.0;
        ddsc5[jj] = 0.0;
        ddsc7[jj] = 0.0;
    }
    damp = particleI.dampingFactor*particleK.dampingFactor;
    if( damp != zero ){
        pgamma = particleI.thole < particleK.thole ? particleI.thole : particleK.thole;
        RealOpenMM ratio = (r/damp);
        damp   = -pgamma * ratio*ratio*ratio;
        if( damp > -50.0){
            expdamp = exp(damp);
            scale3 = 1.0 - expdamp;
            scale5 = 1.0 - (1.0-damp)*expdamp;
            scale7 = 1.0 - (1.0-damp+0.6*damp*damp)*expdamp;
            temp3 = -3.0 * damp * expdamp / r2;
            temp5 = -damp;
            temp7 = -0.2 - 0.6*damp;
            ddsc3[0] = temp3 * xr;
            ddsc3[1] = temp3 * yr;
            ddsc3[2] = temp3 * zr;
            ddsc5[0] = temp5 * ddsc3[0];
            ddsc5[1] = temp5 * ddsc3[1];
            ddsc5[2] = temp5 * ddsc3[2];
            ddsc7[0] = temp7 * ddsc5[0];
            ddsc7[1] = temp7 * ddsc5[1];
            ddsc7[2] = temp7 * ddsc5[2];
       }
   }
   scale3i = scale3 * scalingFactors[U_SCALE];
   scale5i = scale5 * scalingFactors[U_SCALE];
   scale7i = scale7 * scalingFactors[U_SCALE];

   dsc3 = scale3 * scalingFactors[D_SCALE];
   dsc5 = scale5 * scalingFactors[D_SCALE];
   dsc7 = scale7 * scalingFactors[D_SCALE];

   psc3 = scale3 * scalingFactors[P_SCALE];
   psc5 = scale5 * scalingFactors[P_SCALE];
   psc7 = scale7 * scalingFactors[P_SCALE];

    // construct necessary auxiliary vectors

    dixdk[0] = particleI.dipole[1]*particleK.dipole[2] - particleI.dipole[2]*particleK.dipole[1];
    dixdk[1] = particleI.dipole[2]*particleK.dipole[0] - particleI.dipole[0]*particleK.dipole[2];
    dixdk[2] = particleI.dipole[0]*particleK.dipole[1] - particleI.dipole[1]*particleK.dipole[0];

    dixuk[0] = particleI.dipole[1]*particleK.inducedDipole[2] - particleI.dipole[2]*particleK.inducedDipole[1];
    dixuk[1] = particleI.dipole[2]*particleK.inducedDipole[0] - particleI.dipole[0]*particleK.inducedDipole[2];
    dixuk[2] = particleI.dipole[0]*particleK.inducedDipole[1] - particleI.dipole[1]*particleK.inducedDipole[0];

    dkxui[0] = particleK.dipole[1]*particleI.inducedDipole[2] - particleK.dipole[2]*particleI.inducedDipole[1];
    dkxui[1] = particleK.dipole[2]*particleI.inducedDipole[0] - particleK.dipole[0]*particleI.inducedDipole[2];
    dkxui[2] = particleK.dipole[0]*particleI.inducedDipole[1] - particleK.dipole[1]*particleI.inducedDipole[0];

    dixukp[0] = particleI.dipole[1]*particleK.inducedDipolePolar[2] - particleI.dipole[2]*particleK.inducedDipolePolar[1];
    dixukp[1] = particleI.dipole[2]*particleK.inducedDipolePolar[0] - particleI.dipole[0]*particleK.inducedDipolePolar[2];
    dixukp[2] = particleI.dipole[0]*particleK.inducedDipolePolar[1] - particleI.dipole[1]*particleK.inducedDipolePolar[0];

    dkxuip[0] = particleK.dipole[1]*particleI.inducedDipolePolar[2] - particleK.dipole[2]*particleI.inducedDipolePolar[1];
    dkxuip[1] = particleK.dipole[2]*particleI.inducedDipolePolar[0] - particleK.dipole[0]*particleI.inducedDipolePolar[2];
    dkxuip[2] = particleK.dipole[0]*particleI.inducedDipolePolar[1] - particleK.dipole[1]*particleI.inducedDipolePolar[0];

    dixr[0] = particleI.dipole[1]*zr -particleI.dipole[2]*yr;
    dixr[1] = particleI.dipole[2]*xr -particleI.dipole[0]*zr;
    dixr[2] = particleI.dipole[0]*yr -particleI.dipole[1]*xr;

    dkxr[0] = particleK.dipole[1]*zr - particleK.dipole[2]*yr;
    dkxr[1] = particleK.dipole[2]*xr - particleK.dipole[0]*zr;
    dkxr[2] = particleK.dipole[0]*yr - particleK.dipole[1]*xr;

    qir[0] = qi[0]*xr + qi[3]*yr + qi[6]*zr;
    qir[1] = qi[1]*xr + qi[4]*yr + qi[7]*zr;
    qir[2] = qi[2]*xr + qi[5]*yr + qi[8]*zr;
    qkr[0] = qk[0]*xr + qk[3]*yr + qk[6]*zr;
    qkr[1] = qk[1]*xr + qk[4]*yr + qk[7]*zr;
    qkr[2] = qk[2]*xr + qk[5]*yr + qk[8]*zr;

    qiqkr[0] = qi[0]*qkr[0] + qi[3]*qkr[1] + qi[6]*qkr[2];
    qiqkr[1] = qi[1]*qkr[0] + qi[4]*qkr[1] + qi[7]*qkr[2];
    qiqkr[2] = qi[2]*qkr[0] + qi[5]*qkr[1] + qi[8]*qkr[2];

    qkqir[0] = qk[0]*qir[0] + qk[3]*qir[1] + qk[6]*qir[2];
    qkqir[1] = qk[1]*qir[0] + qk[4]*qir[1] + qk[7]*qir[2];
    qkqir[2] = qk[2]*qir[0] + qk[5]*qir[1] + qk[8]*qir[2];

    qixqk[0] = qi[1]*qk[2] + qi[4]*qk[5] + qi[7]*qk[8] - qi[2]*qk[1] - qi[5]*qk[4] - qi[8]*qk[7];
    qixqk[1] = qi[2]*qk[0] + qi[5]*qk[3] + qi[8]*qk[6] - qi[0]*qk[2] - qi[3]*qk[5] - qi[6]*qk[8];
    qixqk[2] = qi[0]*qk[1] + qi[3]*qk[4] + qi[6]*qk[7] - qi[1]*qk[0] - qi[4]*qk[3] - qi[7]*qk[6];

    rxqir[0] = yr*qir[2] - zr*qir[1];
    rxqir[1] = zr*qir[0] - xr*qir[2];
    rxqir[2] = xr*qir[1] - yr*qir[0];

    rxqkr[0] = yr*qkr[2] - zr*qkr[1];
    rxqkr[1] = zr*qkr[0] - xr*qkr[2];
    rxqkr[2] = xr*qkr[1] - yr*qkr[0];

    rxqikr[0] = yr*qiqkr[2] - zr*qiqkr[1];
    rxqikr[1] = zr*qiqkr[0] - xr*qiqkr[2];
    rxqikr[2] = xr*qiqkr[1] - yr*qiqkr[0];

    rxqkir[0] = yr*qkqir[2] - zr*qkqir[1];
    rxqkir[1] = zr*qkqir[0] - xr*qkqir[2];
    rxqkir[2] = xr*qkqir[1] - yr*qkqir[0];

    qkrxqir[0] = qkr[1]*qir[2] - qkr[2]*qir[1];
    qkrxqir[1] = qkr[2]*qir[0] - qkr[0]*qir[2];
    qkrxqir[2] = qkr[0]*qir[1] - qkr[1]*qir[0];

    qidk[0] = qi[0]*particleK.dipole[0] + qi[3]*particleK.dipole[1] + qi[6]*particleK.dipole[2];
    qidk[1] = qi[1]*particleK.dipole[0] + qi[4]*particleK.dipole[1] + qi[7]*particleK.dipole[2];
    qidk[2] = qi[2]*particleK.dipole[0] + qi[5]*particleK.dipole[1] + qi[8]*particleK.dipole[2];

    qkdi[0] = qk[0]*particleI.dipole[0] + qk[3]*particleI.dipole[1] + qk[6]*particleI.dipole[2];
    qkdi[1] = qk[1]*particleI.dipole[0] + qk[4]*particleI.dipole[1] + qk[7]*particleI.dipole[2];
    qkdi[2] = qk[2]*particleI.dipole[0] + qk[5]*particleI.dipole[1] + qk[8]*particleI.dipole[2];

    qiuk[0] = qi[0]*particleK.inducedDipole[0] + qi[3]*particleK.inducedDipole[1] + qi[6]*particleK.inducedDipole[2];
    qiuk[1] = qi[1]*particleK.inducedDipole[0] + qi[4]*particleK.inducedDipole[1] + qi[7]*particleK.inducedDipole[2];
    qiuk[2] = qi[2]*particleK.inducedDipole[0] + qi[5]*particleK.inducedDipole[1] + qi[8]*particleK.inducedDipole[2];

    qkui[0] = qk[0]*particleI.inducedDipole[0] + qk[3]*particleI.inducedDipole[1] + qk[6]*particleI.inducedDipole[2];
    qkui[1] = qk[1]*particleI.inducedDipole[0] + qk[4]*particleI.inducedDipole[1] + qk[7]*particleI.inducedDipole[2];
    qkui[2] = qk[2]*particleI.inducedDipole[0] + qk[5]*particleI.inducedDipole[1] + qk[8]*particleI.inducedDipole[2];

    qiukp[0] = qi[0]*particleK.inducedDipolePolar[0] + qi[3]*particleK.inducedDipolePolar[1] + qi[6]*particleK.inducedDipolePolar[2];
    qiukp[1] = qi[1]*particleK.inducedDipolePolar[0] + qi[4]*particleK.inducedDipolePolar[1] + qi[7]*particleK.inducedDipolePolar[2];
    qiukp[2] = qi[2]*particleK.inducedDipolePolar[0] + qi[5]*particleK.inducedDipolePolar[1] + qi[8]*particleK.inducedDipolePolar[2];

    qkuip[0] = qk[0]*particleI.inducedDipolePolar[0] + qk[3]*particleI.inducedDipolePolar[1] + qk[6]*particleI.inducedDipolePolar[2];
    qkuip[1] = qk[1]*particleI.inducedDipolePolar[0] + qk[4]*particleI.inducedDipolePolar[1] + qk[7]*particleI.inducedDipolePolar[2];
    qkuip[2] = qk[2]*particleI.inducedDipolePolar[0] + qk[5]*particleI.inducedDipolePolar[1] + qk[8]*particleI.inducedDipolePolar[2];

    dixqkr[0] = particleI.dipole[1]*qkr[2] -particleI.dipole[2]*qkr[1];
    dixqkr[1] = particleI.dipole[2]*qkr[0] -particleI.dipole[0]*qkr[2];
    dixqkr[2] = particleI.dipole[0]*qkr[1] -particleI.dipole[1]*qkr[0];

    dkxqir[0] = particleK.dipole[1]*qir[2] - particleK.dipole[2]*qir[1];
    dkxqir[1] = particleK.dipole[2]*qir[0] - particleK.dipole[0]*qir[2];
    dkxqir[2] = particleK.dipole[0]*qir[1] - particleK.dipole[1]*qir[0];

    uixqkr[0] = particleI.inducedDipole[1]*qkr[2] - particleI.inducedDipole[2]*qkr[1];
    uixqkr[1] = particleI.inducedDipole[2]*qkr[0] - particleI.inducedDipole[0]*qkr[2];
    uixqkr[2] = particleI.inducedDipole[0]*qkr[1] - particleI.inducedDipole[1]*qkr[0];

    ukxqir[0] = particleK.inducedDipole[1]*qir[2] - particleK.inducedDipole[2]*qir[1];
    ukxqir[1] = particleK.inducedDipole[2]*qir[0] - particleK.inducedDipole[0]*qir[2];
    ukxqir[2] = particleK.inducedDipole[0]*qir[1] - particleK.inducedDipole[1]*qir[0];

    uixqkrp[0] = particleI.inducedDipolePolar[1]*qkr[2] - particleI.inducedDipolePolar[2]*qkr[1];
    uixqkrp[1] = particleI.inducedDipolePolar[2]*qkr[0] - particleI.inducedDipolePolar[0]*qkr[2];
    uixqkrp[2] = particleI.inducedDipolePolar[0]*qkr[1] - particleI.inducedDipolePolar[1]*qkr[0];

    ukxqirp[0] = particleK.inducedDipolePolar[1]*qir[2] - particleK.inducedDipolePolar[2]*qir[1];
    ukxqirp[1] = particleK.inducedDipolePolar[2]*qir[0] - particleK.inducedDipolePolar[0]*qir[2];
    ukxqirp[2] = particleK.inducedDipolePolar[0]*qir[1] - particleK.inducedDipolePolar[1]*qir[0];

    rxqidk[0] = yr*qidk[2] - zr*qidk[1];
    rxqidk[1] = zr*qidk[0] - xr*qidk[2];
    rxqidk[2] = xr*qidk[1] - yr*qidk[0];

    rxqkdi[0] = yr*qkdi[2] - zr*qkdi[1];
    rxqkdi[1] = zr*qkdi[0] - xr*qkdi[2];
    rxqkdi[2] = xr*qkdi[1] - yr*qkdi[0];

    rxqiuk[0] = yr*qiuk[2] - zr*qiuk[1];
    rxqiuk[1] = zr*qiuk[0] - xr*qiuk[2];
    rxqiuk[2] = xr*qiuk[1] - yr*qiuk[0];

    rxqkui[0] = yr*qkui[2] - zr*qkui[1];
    rxqkui[1] = zr*qkui[0] - xr*qkui[2];
    rxqkui[2] = xr*qkui[1] - yr*qkui[0];

    rxqiukp[0] = yr*qiukp[2] - zr*qiukp[1];
    rxqiukp[1] = zr*qiukp[0] - xr*qiukp[2];
    rxqiukp[2] = xr*qiukp[1] - yr*qiukp[0];

    rxqkuip[0] = yr*qkuip[2] - zr*qkuip[1];
    rxqkuip[1] = zr*qkuip[0] - xr*qkuip[2];
    rxqkuip[2] = xr*qkuip[1] - yr*qkuip[0];

    // calculate scalar products for permanent components

    sc[1] = particleI.dipole[0]*particleK.dipole[0] +particleI.dipole[1]*particleK.dipole[1] +particleI.dipole[2]*particleK.dipole[2];
    sc[2] = particleI.dipole[0]*xr +particleI.dipole[1]*yr +particleI.dipole[2]*zr;
    sc[3] = particleK.dipole[0]*xr + particleK.dipole[1]*yr + particleK.dipole[2]*zr;
    sc[4] = qir[0]*xr + qir[1]*yr + qir[2]*zr;
    sc[5] = qkr[0]*xr + qkr[1]*yr + qkr[2]*zr;
    sc[6] = qir[0]*particleK.dipole[0] + qir[1]*particleK.dipole[1] + qir[2]*particleK.dipole[2];
    sc[7] = qkr[0]*particleI.dipole[0] + qkr[1]*particleI.dipole[1] + qkr[2]*particleI.dipole[2];
    sc[8] = qir[0]*qkr[0] + qir[1]*qkr[1] + qir[2]*qkr[2];
    sc[9] = qi[0]*qk[0] + qi[1]*qk[1] + qi[2]*qk[2] +
             qi[3]*qk[3] + qi[4]*qk[4] + qi[5]*qk[5] +
             qi[6]*qk[6] + qi[7]*qk[7] + qi[8]*qk[8];
 
    // calculate scalar products for induced components

    sci[0] = particleI.inducedDipole[0]*particleK.dipole[0] + particleI.inducedDipole[1]*particleK.dipole[1] + particleI.inducedDipole[2]*particleK.dipole[2] +particleI.dipole[0]*particleK.inducedDipole[0] +particleI.dipole[1]*particleK.inducedDipole[1] +particleI.dipole[2]*particleK.inducedDipole[2];
    sci[1] = particleI.inducedDipole[0]*particleK.inducedDipole[0] + particleI.inducedDipole[1]*particleK.inducedDipole[1] + particleI.inducedDipole[2]*particleK.inducedDipole[2];
    sci[2] = particleI.inducedDipole[0]*xr + particleI.inducedDipole[1]*yr + particleI.inducedDipole[2]*zr;
    sci[3] = particleK.inducedDipole[0]*xr + particleK.inducedDipole[1]*yr + particleK.inducedDipole[2]*zr;
    sci[6] = qir[0]*particleK.inducedDipole[0] + qir[1]*particleK.inducedDipole[1] + qir[2]*particleK.inducedDipole[2];
    sci[7] = qkr[0]*particleI.inducedDipole[0] + qkr[1]*particleI.inducedDipole[1] + qkr[2]*particleI.inducedDipole[2];
    scip[0] = particleI.inducedDipolePolar[0]*particleK.dipole[0] + particleI.inducedDipolePolar[1]*particleK.dipole[1] + particleI.inducedDipolePolar[2]*particleK.dipole[2] +particleI.dipole[0]*particleK.inducedDipolePolar[0] +particleI.dipole[1]*particleK.inducedDipolePolar[1] +particleI.dipole[2]*particleK.inducedDipolePolar[2];
    scip[1] = particleI.inducedDipole[0]*particleK.inducedDipolePolar[0]+particleI.inducedDipole[1]*particleK.inducedDipolePolar[1] + particleI.inducedDipole[2]*particleK.inducedDipolePolar[2]+particleI.inducedDipolePolar[0]*particleK.inducedDipole[0] + particleI.inducedDipolePolar[1]*particleK.inducedDipole[1]+particleI.inducedDipolePolar[2]*particleK.inducedDipole[2];
    scip[2] = particleI.inducedDipolePolar[0]*xr + particleI.inducedDipolePolar[1]*yr + particleI.inducedDipolePolar[2]*zr;
    scip[3] = particleK.inducedDipolePolar[0]*xr + particleK.inducedDipolePolar[1]*yr + particleK.inducedDipolePolar[2]*zr;
    scip[6] = qir[0]*particleK.inducedDipolePolar[0] + qir[1]*particleK.inducedDipolePolar[1] + qir[2]*particleK.inducedDipolePolar[2];
    scip[7] = qkr[0]*particleI.inducedDipolePolar[0] + qkr[1]*particleI.inducedDipolePolar[1] + qkr[2]*particleI.inducedDipolePolar[2];

    // calculate the gl functions for permanent components

    gl[0] = particleI.charge*particleK.charge;
    gl[1] = particleK.charge*sc[2] - particleI.charge*sc[3];
    gl[2] = particleI.charge*sc[5] + particleK.charge*sc[4] - sc[2]*sc[3];
    gl[3] = sc[2]*sc[5] - sc[3]*sc[4];
    gl[4] = sc[4]*sc[5];
    gl[5] = -4.0 * sc[8];
    gl[6] = sc[1];
    gl[7] = 2.0 * (sc[6]-sc[7]);
    gl[8] = 2.0 * sc[9];

    // calculate the gl functions for induced components

    gli[0] = particleK.charge*sci[2] - particleI.charge*sci[3];
    gli[1] = -sc[2]*sci[3] - sci[2]*sc[3];
    gli[2] = sci[2]*sc[5] - sci[3]*sc[4];
    gli[5] = sci[0];
    gli[6] = 2.0 * (sci[6]-sci[7]);
    glip[0] = particleK.charge*scip[2] - particleI.charge*scip[3];
    glip[1] = -sc[2]*scip[3] - scip[2]*sc[3];
    glip[2] = scip[2]*sc[5] - scip[3]*sc[4];
    glip[5] = scip[0];
    glip[6] = 2.0 * (scip[6]-scip[7]);

    // compute the energy contributions for this interaction

    e = rr1*gl[0] + rr3*(gl[1]+gl[6]) + rr5*(gl[2]+gl[7]+gl[8]) + rr7*(gl[3]+gl[5]) + rr9*gl[4];
    ei = 0.5*(rr3*(gli[0]+gli[5])*psc3 + rr5*(gli[1]+gli[6])*psc5 + rr7*gli[2]*psc7);
    e = f * scalingFactors[M_SCALE] * e;
    ei = f * ei;
    RealOpenMM energy = e + ei;

    // intermediate variables for the permanent components

    gf[0] = rr3*gl[0] + rr5*(gl[1]+gl[6]) + rr7*(gl[2]+gl[7]+gl[8]) + rr9*(gl[3]+gl[5]) + rr11*gl[4];
    gf[1] = -particleK.charge*rr3 + sc[3]*rr5 - sc[5]*rr7;
    gf[2] = particleI.charge*rr3 + sc[2]*rr5 + sc[4]*rr7;
    gf[3] = 2.0 * rr5;
    gf[4] = 2.0 * (-particleK.charge*rr5+sc[3]*rr7-sc[5]*rr9);
    gf[5] = 2.0 * (-particleI.charge*rr5-sc[2]*rr7-sc[4]*rr9);
    gf[6] = 4.0 * rr7;

    // intermediate variables for the induced components

    gfi[0] = 0.5 * rr5 * ((gli[0]+gli[5])*psc3 + (glip[0]+glip[5])*dsc3 + scip[1]*scale3i)
           + 0.5 * rr7 * ((gli[6]+gli[1])*psc5 + (glip[6]+glip[1])*dsc5 - (sci[2]*scip[3]+scip[2]*sci[3])*scale5i)
           + 0.5 * rr9 * (gli[2]*psc7+glip[2]*dsc7);

    gfi[1] = -rr3*particleK.charge + rr5*sc[3] - rr7*sc[5];
    gfi[2] =  rr3*particleI.charge + rr5*sc[2] + rr7*sc[4];
    gfi[3] = 2.0 * rr5;
    gfi[4] = rr7 * (sci[3]*psc7+scip[3]*dsc7);
    gfi[5] = -rr7 * (sci[2]*psc7+scip[2]*dsc7);

    // get the permanent force components

    ftm2[0] = gf[0]*xr + gf[1]*particleI.dipole[0] + gf[2]*particleK.dipole[0]
                 + gf[3]*(qkdi[0]-qidk[0]) + gf[4]*qir[0]
                 + gf[5]*qkr[0] + gf[6]*(qiqkr[0]+qkqir[0]);
    ftm2[1] = gf[0]*yr + gf[1]*particleI.dipole[1] + gf[2]*particleK.dipole[1]
                 + gf[3]*(qkdi[1]-qidk[1]) + gf[4]*qir[1]
                 + gf[5]*qkr[1] + gf[6]*(qiqkr[1]+qkqir[1]);
    ftm2[2] = gf[0]*zr + gf[1]*particleI.dipole[2] + gf[2]*particleK.dipole[2]
                 + gf[3]*(qkdi[2]-qidk[2]) + gf[4]*qir[2]
                 + gf[5]*qkr[2] + gf[6]*(qiqkr[2]+qkqir[2]);


    // get the induced force components

    ftm2i[0] = gfi[0]*xr + 0.5*
      (- rr3*particleK.charge*(particleI.inducedDipole[0]*psc3+particleI.inducedDipolePolar[0]*dsc3)
       + rr5*sc[3]*(particleI.inducedDipole[0]*psc5+particleI.inducedDipolePolar[0]*dsc5)
       - rr7*sc[5]*(particleI.inducedDipole[0]*psc7+particleI.inducedDipolePolar[0]*dsc7))
       + (rr3*particleI.charge*(particleK.inducedDipole[0]*psc3+particleK.inducedDipolePolar[0]*dsc3)
       + rr5*sc[2]*(particleK.inducedDipole[0]*psc5+particleK.inducedDipolePolar[0]*dsc5)
       + rr7*sc[4]*(particleK.inducedDipole[0]*psc7+particleK.inducedDipolePolar[0]*dsc7))*0.5
       + rr5*scale5i*(sci[3]*particleI.inducedDipolePolar[0]+scip[3]*particleI.inducedDipole[0]
       + sci[2]*particleK.inducedDipolePolar[0]+scip[2]*particleK.inducedDipole[0])*0.5
       + 0.5*(sci[3]*psc5+scip[3]*dsc5)*rr5*particleI.dipole[0]
       + 0.5*(sci[2]*psc5+scip[2]*dsc5)*rr5*particleK.dipole[0]
       + 0.5*gfi[3]*((qkui[0]-qiuk[0])*psc5
       + (qkuip[0]-qiukp[0])*dsc5)
       + gfi[4]*qir[0] + gfi[5]*qkr[0];

    ftm2i[1] = gfi[0]*yr + 0.5*
      (- rr3*particleK.charge*(particleI.inducedDipole[1]*psc3+particleI.inducedDipolePolar[1]*dsc3)
       + rr5*sc[3]*(particleI.inducedDipole[1]*psc5+particleI.inducedDipolePolar[1]*dsc5)
       - rr7*sc[5]*(particleI.inducedDipole[1]*psc7+particleI.inducedDipolePolar[1]*dsc7))
       + (rr3*particleI.charge*(particleK.inducedDipole[1]*psc3+particleK.inducedDipolePolar[1]*dsc3)
       + rr5*sc[2]*(particleK.inducedDipole[1]*psc5+particleK.inducedDipolePolar[1]*dsc5)
       + rr7*sc[4]*(particleK.inducedDipole[1]*psc7+particleK.inducedDipolePolar[1]*dsc7))*0.5
       + rr5*scale5i*(sci[3]*particleI.inducedDipolePolar[1]+scip[3]*particleI.inducedDipole[1]
       + sci[2]*particleK.inducedDipolePolar[1]+scip[2]*particleK.inducedDipole[1])*0.5
       + 0.5*(sci[3]*psc5+scip[3]*dsc5)*rr5*particleI.dipole[1]
       + 0.5*(sci[2]*psc5+scip[2]*dsc5)*rr5*particleK.dipole[1]
       + 0.5*gfi[3]*((qkui[1]-qiuk[1])*psc5
       + (qkuip[1]-qiukp[1])*dsc5)
       + gfi[4]*qir[1] + gfi[5]*qkr[1];

    ftm2i[2] = gfi[0]*zr  + 0.5*
      (- rr3*particleK.charge*(particleI.inducedDipole[2]*psc3+particleI.inducedDipolePolar[2]*dsc3)
       + rr5*sc[3]*(particleI.inducedDipole[2]*psc5+particleI.inducedDipolePolar[2]*dsc5)
       - rr7*sc[5]*(particleI.inducedDipole[2]*psc7+particleI.inducedDipolePolar[2]*dsc7))
       + (rr3*particleI.charge*(particleK.inducedDipole[2]*psc3+particleK.inducedDipolePolar[2]*dsc3)
       + rr5*sc[2]*(particleK.inducedDipole[2]*psc5+particleK.inducedDipolePolar[2]*dsc5)
       + rr7*sc[4]*(particleK.inducedDipole[2]*psc7+particleK.inducedDipolePolar[2]*dsc7))*0.5
       + rr5*scale5i*(sci[3]*particleI.inducedDipolePolar[2]+scip[3]*particleI.inducedDipole[2]
       + sci[2]*particleK.inducedDipolePolar[2]+scip[2]*particleK.inducedDipole[2])*0.5
       + 0.5*(sci[3]*psc5+scip[3]*dsc5)*rr5*particleI.dipole[2]
       + 0.5*(sci[2]*psc5+scip[2]*dsc5)*rr5*particleK.dipole[2]
       + 0.5*gfi[3]*((qkui[2]-qiuk[2])*psc5
       + (qkuip[2]-qiukp[2])*dsc5)
       + gfi[4]*qir[2] + gfi[5]*qkr[2];

    // account for partially excluded induced interactions

    temp3 = 0.5 * rr3 * ((gli[0]+gli[5])*scalingFactors[P_SCALE] +(glip[0]+glip[5])*scalingFactors[D_SCALE]);
    temp5 = 0.5 * rr5 * ((gli[1]+gli[6])*scalingFactors[P_SCALE] +(glip[1]+glip[6])*scalingFactors[D_SCALE]);
    temp7 = 0.5 * rr7 * (gli[2]*scalingFactors[P_SCALE] +glip[2]*scalingFactors[D_SCALE]);

    fridmp[0] = temp3*ddsc3[0] + temp5*ddsc5[0] + temp7*ddsc7[0];
    fridmp[1] = temp3*ddsc3[1] + temp5*ddsc5[1] + temp7*ddsc7[1];
    fridmp[2] = temp3*ddsc3[2] + temp5*ddsc5[2] + temp7*ddsc7[2];

    // find some scaling terms for induced-induced force

    temp3 = 0.5 * rr3 * scalingFactors[U_SCALE] * scip[1];
    temp5 = -0.5 * rr5 * scalingFactors[U_SCALE]
               * (sci[2]*scip[3]+scip[2]*sci[3]);
    findmp[0] = temp3*ddsc3[0] + temp5*ddsc5[0];
    findmp[1] = temp3*ddsc3[1] + temp5*ddsc5[1];
    findmp[2] = temp3*ddsc3[2] + temp5*ddsc5[2];

    // modify induced force for partially excluded interactions

    ftm2i[0] = ftm2i[0] - fridmp[0] - findmp[0];
    ftm2i[1] = ftm2i[1] - fridmp[1] - findmp[1];
    ftm2i[2] = ftm2i[2] - fridmp[2] - findmp[2];

    // correction to convert mutual to direct polarization force

/*
    if( poltyp .eq. 'DIRECT'){
       gfd = 0.5 * (rr5*scip[1]*scale3i
             - rr7*(scip[2]*sci[3]+sci[2]*scip[3])*scale5i);
       temp5 = 0.5 * rr5 * scale5i;
       fdir[0] = gfd*xr + temp5
                    * (sci[3]*particleI.inducedDipolePolar[0]+scip[3]*particleI.inducedDipole[0]
                      +sci[2]*particleK.inducedDipolePolar[0]+scip[2]*particleK.inducedDipole[0]);
       fdir[1] = gfd*yr + temp5
                    * (sci[3]*particleI.inducedDipolePolar[1]+scip[3]*particleI.inducedDipole[1]
                      +sci[2]*particleK.inducedDipolePolar[1]+scip[2]*particleK.inducedDipole[1]);
       fdir[2] = gfd*zr + temp5
                    * (sci[3]*particleI.inducedDipolePolar[2]+scip[3]*particleI.inducedDipole[2]
                      +sci[2]*particleK.inducedDipolePolar[2]+scip[2]*particleK.inducedDipole[2]);
       ftm2i[0] = ftm2i[0] - fdir[0] + findmp[0];
       ftm2i[1] = ftm2i[1] - fdir[1] + findmp[1];
       ftm2i[2] = ftm2i[2] - fdir[2] + findmp[2];
    }
*/

    // intermediate terms for induced torque on multipoles

    gti[1] = 0.5 * rr5 * (sci[3]*psc5+scip[3]*dsc5);
    gti[2] = 0.5 * rr5 * (sci[2]*psc5+scip[2]*dsc5);
    gti[3] = gfi[3];
    gti[4] = gfi[4];
    gti[5] = gfi[5];

    // get the permanent torque components

    ttm2[0] = -rr3*dixdk[0] + gf[1]*dixr[0] - gf[4]*rxqir[0]
      + gf[3]*(dixqkr[0]+dkxqir[0]+rxqidk[0]-2.0*qixqk[0])
      - gf[6]*(rxqikr[0]+qkrxqir[0]);

    ttm2[1] = -rr3*dixdk[1] + gf[1]*dixr[1] - gf[4]*rxqir[1]
      + gf[3]*(dixqkr[1]+dkxqir[1]+rxqidk[1]-2.0*qixqk[1])
      - gf[6]*(rxqikr[1]+qkrxqir[1]);

    ttm2[2] = -rr3*dixdk[2] + gf[1]*dixr[2] - gf[4]*rxqir[2]
      + gf[3]*(dixqkr[2]+dkxqir[2]+rxqidk[2]-2.0*qixqk[2])
      - gf[6]*(rxqikr[2]+qkrxqir[2]);

    ttm3[0] = rr3*dixdk[0] + gf[2]*dkxr[0] - gf[5]*rxqkr[0]
      - gf[3]*(dixqkr[0]+dkxqir[0]+rxqkdi[0]-2.0*qixqk[0])
      - gf[6]*(rxqkir[0]-qkrxqir[0]);

    ttm3[1] = rr3*dixdk[1] + gf[2]*dkxr[1] - gf[5]*rxqkr[1]
      - gf[3]*(dixqkr[1]+dkxqir[1]+rxqkdi[1]-2.0*qixqk[1])
      - gf[6]*(rxqkir[1]-qkrxqir[1]);

    ttm3[2] = rr3*dixdk[2] + gf[2]*dkxr[2] - gf[5]*rxqkr[2]
      - gf[3]*(dixqkr[2]+dkxqir[2]+rxqkdi[2]-2.0*qixqk[2])
      - gf[6]*(rxqkir[2]-qkrxqir[2]);

    // get the induced torque components
   
    ttm2i[0] = -rr3*(dixuk[0]*psc3+dixukp[0]*dsc3)*0.5
      + gti[1]*dixr[0] + gti[3]*((ukxqir[0]+rxqiuk[0])*psc5
      +(ukxqirp[0]+rxqiukp[0])*dsc5)*0.5 - gti[4]*rxqir[0];

    ttm2i[1] = -rr3*(dixuk[1]*psc3+dixukp[1]*dsc3)*0.5
      + gti[1]*dixr[1] + gti[3]*((ukxqir[1]+rxqiuk[1])*psc5
      +(ukxqirp[1]+rxqiukp[1])*dsc5)*0.5 - gti[4]*rxqir[1];

    ttm2i[2] = -rr3*(dixuk[2]*psc3+dixukp[2]*dsc3)*0.5
      + gti[1]*dixr[2] + gti[3]*((ukxqir[2]+rxqiuk[2])*psc5
      +(ukxqirp[2]+rxqiukp[2])*dsc5)*0.5 - gti[4]*rxqir[2];

    ttm3i[0] = -rr3*(dkxui[0]*psc3+dkxuip[0]*dsc3)*0.5
      + gti[2]*dkxr[0] - gti[3]*((uixqkr[0]+rxqkui[0])*psc5
      +(uixqkrp[0]+rxqkuip[0])*dsc5)*0.5 - gti[5]*rxqkr[0];

    ttm3i[1] = -rr3*(dkxui[1]*psc3+dkxuip[1]*dsc3)*0.5
      + gti[2]*dkxr[1] - gti[3]*((uixqkr[1]+rxqkui[1])*psc5
      +(uixqkrp[1]+rxqkuip[1])*dsc5)*0.5 - gti[5]*rxqkr[1];

    ttm3i[2] = -rr3*(dkxui[2]*psc3+dkxuip[2]*dsc3)*0.5
      + gti[2]*dkxr[2] - gti[3]*((uixqkr[2]+rxqkui[2])*psc5
      +(uixqkrp[2]+rxqkuip[2])*dsc5)*0.5 - gti[5]*rxqkr[2];

    // handle the case where scaling is used

    for( unsigned int jj = 0; jj < 3; jj++ ){
       ftm2[jj]   = f * ftm2[jj] * scalingFactors[M_SCALE];
       ftm2i[jj]  = f * ftm2i[jj];
       ttm2[jj]   = f * ttm2[jj] * scalingFactors[M_SCALE];
       ttm2i[jj]  = f * ttm2i[jj];
       ttm3[jj]   = f * ttm3[jj] * scalingFactors[M_SCALE];
       ttm3i[jj]  = f * ttm3i[jj];
    }

    // increment forces

    for( unsigned int jj = 0; jj < 3; jj++ ){

        forces[particleI.particleIndex][jj] -= ftm2[jj] + ftm2i[jj];
        torque[particleI.particleIndex][jj] += ttm2[jj] + ttm2i[jj];

        forces[particleK.particleIndex][jj] += ftm2[jj] + ftm2i[jj];
        torque[particleK.particleIndex][jj] += ttm3[jj] + ttm3i[jj];
    }

    return energy;
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateNoCutoffElectrostatic( std::vector<MultipoleParticleData>& particleData, RealOpenMM** forces ) const {


    // ---------------------------------------------------------------------------------------

    static const std::string methodName    = "AmoebaReferenceMultipoleForce::calculateNoCutoffElectrostatic: ";

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;

    static const int debug                 = 1;
    FILE* log                              = stderr;

    // ---------------------------------------------------------------------------------------

    // initialize forces/energy and scaleing factors

    std::vector<Vec3> torques( particleData.size() );
    for( unsigned int ii = 0; ii <  particleData.size(); ii++ ){
        torques[ii]   = Vec3( zero, zero, zero );
    }

    if( debug ){
        if( log ){
            (void) fprintf( log, "\n%s Zeroing forces\n", methodName.c_str() );
        }
        for( unsigned int ii = 0; ii <  particleData.size(); ii++ ){
            forces[ii][0] = zero;
            forces[ii][1] = zero;
            forces[ii][2] = zero;
        }
    }
    RealOpenMM energy = zero;

    RealOpenMM scaleFactors[LAST_SCALE_TYPE_INDEX];
    for( unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++ ){
        scaleFactors[kk] = one;
    }   

    // main loop over particle pairs

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){

            if( jj <= _maxScaleIndex[ii] ){
                getScaleFactors( ii, jj, scaleFactors);
            }

            energy += calculateNoCutoffElectrostaticPairIxn( particleData[ii], particleData[jj], scaleFactors, forces, torques );

            if( jj <= _maxScaleIndex[ii] ){
                for( unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++ ){
                    scaleFactors[kk] = one;
                }
            }
        }
    }

    if( log ){
        RealOpenMM conversion = 1.0/4.184;
        (void) fprintf( log, "Force & torques energy=%15.7e\n", 10.0*energy*conversion );
        conversion *= -0.1;
        unsigned int maxPrint = 10;
        for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
            (void) fprintf( log, "%6u [%15.7e %15.7e %15.7e] [%15.7e %15.7e %15.7e]\n", ii,
                            conversion*forces[ii][0],   conversion*forces[ii][1],  conversion*forces[ii][2], 
                            conversion*torques[ii][0], conversion*torques[ii][1], conversion*torques[ii][2] );
            if( ii == maxPrint && 2*maxPrint < particleData.size() ){
                ii = particleData.size() - maxPrint;
            }
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateNoCutoffForceAndEnergy( unsigned int numParticles,
                                                                           RealOpenMM** particlePositions, 
                                                                           const std::vector<RealOpenMM>& charges,
                                                                           const std::vector<RealOpenMM>& dipoles,
                                                                           const std::vector<RealOpenMM>& quadrupoles,
                                                                           const std::vector<RealOpenMM>& tholes,
                                                                           const std::vector<RealOpenMM>& dampingFactors,
                                                                           const std::vector<RealOpenMM>& polarity,
                                                                           const std::vector<int>& axisTypes,
                                                                           const std::vector<int>& multipoleAtomId1s,
                                                                           const std::vector<int>& multipoleAtomId2s,
                                                                           RealOpenMM** forces ){


    // ---------------------------------------------------------------------------------------

    //static const std::string methodName = "AmoebaReferenceMultipoleForce::calculateNoCutoffForceAndEnergy";

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;
    static const RealOpenMM  two           = 2.0;

	 FILE* log                              = stderr;

    // ---------------------------------------------------------------------------------------

    // get lab frame dipole and quadrupoles

    std::vector<MultipoleParticleData> particleData( numParticles );
    loadParticleData( particlePositions, charges, dipoles, quadrupoles,
                      tholes, dampingFactors, polarity, particleData );

    for( unsigned int ii = 0; ii < numParticles; ii++ ){
        if( multipoleAtomId1s[ii] >= 0 && multipoleAtomId2s[ii] >= 0 ){
            applyRotationMatrix( particleData[ii], particleData[multipoleAtomId1s[ii]], particleData[multipoleAtomId2s[ii]], axisTypes[ii] );
        }
    }

#ifdef AMOEBA_DEBUG
    if( log ){
        std::string header = "labFrameQuadrupole and labFrameDipole\n";
        unsigned int printFlag    =  (1 << PARTICLE_DIPOLE) | (1 << PARTICLE_QUADRUPOLE);
        int maxPrint              = 10;
        logParticleData( header, particleData, printFlag, log, maxPrint );
    }
#endif

   // particleData[].field & particleData[].fieldPolar zeroed in loadParticleData 

    for( unsigned int ii = 0; ii < numParticles; ii++ ){
        for( unsigned int jj = ii+1; jj < numParticles; jj++ ){

            // field[0]: field at site ii due site jj's charge/dipole/quadrupole
            // field[1]: field at site jj due site ii's charge/dipole/quadrupole

            // if site jj is less than max covalent scaling index then get/apply scaling constants
            // otherwise add unmodified field and fieldPolar to particle fields 

            RealOpenMM dScale, pScale;
            if( jj <= _maxScaleIndex[ii] ){
                getDScaleAndPScale( ii, jj, dScale, pScale );
            } else {
                dScale = pScale = one;
            }
            calculateFixedEFieldPairIxn( particleData[ii], particleData[jj], dScale, pScale );
        }
    }

#ifdef AMOEBA_DEBUG
    if( log ){
        std::string header        = "Fixed fields\n";
        unsigned int printFlag    =  (1<<PARTICLE_FIELD) | (1<<PARTICLE_FIELD_POLAR);
        int maxPrint              = 10;
        logParticleData( header, particleData, printFlag, log, maxPrint );
    }
#endif

    // get induced dipoles

    calculateNoCutoffInducedDipoles( particleData );

    // check if induced dipoles converged

    if( !getMutualInducedDipoleConverged() ){
        std::stringstream message;
        message << "Induced dipoles did not converge: ";
        message << " iterations="      << getMutualInducedDipoleIterations();
        message << " eps="             << getMutualInducedDipoleEpsilon();
        throw OpenMMException(message.str());
    }

    RealOpenMM energy = calculateNoCutoffElectrostatic( particleData, forces );
    return energy;

}

RealOpenMM AmoebaReferenceMultipoleForce::calculateForceAndEnergy( int numParticles,
                                                         RealOpenMM** particlePositions, 
                                                         const std::vector<RealOpenMM>& charges,
                                                         const std::vector<RealOpenMM>& dipoles,
                                                         const std::vector<RealOpenMM>& quadrupoles,
                                                         const std::vector<RealOpenMM>& tholes,
                                                         const std::vector<RealOpenMM>& dampingFactors,
                                                         const std::vector<RealOpenMM>& polarity,
                                                         const std::vector<int>& axisTypes,
                                                         const std::vector<int>& multipoleAtomId1s,
                                                         const std::vector<int>& multipoleAtomId2s,
                                                         const std::vector< std::vector< std::vector<int> > >& multipoleParticleCovalentInfo,
                                                         RealOpenMM** forces ){

    
    setupScaleMaps( multipoleParticleCovalentInfo );
    if( getNonbondedMethod() == NoCutoff || 1 ){

        return calculateNoCutoffForceAndEnergy( static_cast<unsigned int>(numParticles), particlePositions, charges, dipoles, quadrupoles, tholes, dampingFactors,
                                                polarity, axisTypes, multipoleAtomId1s, multipoleAtomId2s, forces );


    }    
}    

