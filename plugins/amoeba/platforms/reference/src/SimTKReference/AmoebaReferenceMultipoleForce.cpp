
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

#include "AmoebaReferenceMultipoleForce.h"
#include <algorithm>

using std::vector;
using OpenMM::RealVec;

#undef AMOEBA_DEBUG

AmoebaReferenceMultipoleForce::AmoebaReferenceMultipoleForce( ) :
                                                   _nonbondedMethod(NoCutoff),
                                                   _numParticles(0), 
                                                   _electric(138.9354558456),
                                                   _dielectric(1.0),
                                                   _mutualInducedDipoleConverged(0),
                                                   _mutualInducedDipoleIterations(0),
                                                   _maximumMutualInducedDipoleIterations(100),
                                                   _mutualInducedDipoleEpsilon(1.0e+50),
                                                   _mutualInducedDipoleTargetEpsilon(1.0e-04),
                                                   _polarSOR(0.55),
                                                   _debye(48.033324) {
    initialize();
}

AmoebaReferenceMultipoleForce::AmoebaReferenceMultipoleForce( NonbondedMethod nonbondedMethod ) :
                                                   _nonbondedMethod(NoCutoff),
                                                   _numParticles(0), 
                                                   _electric(138.9354558456),
                                                   _dielectric(1.0),
                                                   _mutualInducedDipoleConverged(0),
                                                   _mutualInducedDipoleIterations(0),
                                                   _maximumMutualInducedDipoleIterations(100),
                                                   _mutualInducedDipoleEpsilon(1.0e+50),
                                                   _mutualInducedDipoleTargetEpsilon(1.0e-04),
                                                   _polarSOR(0.55),
                                                   _debye(48.033324) {
    initialize();
}

void AmoebaReferenceMultipoleForce::initialize( void ){

    unsigned int index                              = 0;
    _mScale[index++]                                = 0.0;
    _mScale[index++]                                = 0.0;
    _mScale[index++]                                = 0.0;
    _mScale[index++]                                = 0.4;
    _mScale[index++]                                = 0.8;

    index                                           = 0;
    _dScale[index++]                                = 0.0;
    _dScale[index++]                                = 1.0;
    _dScale[index++]                                = 1.0;
    _dScale[index++]                                = 1.0;
    _dScale[index++]                                = 1.0;

    index                                           = 0;
    _pScale[index++]                                = 0.0;
    _pScale[index++]                                = 0.0;
    _pScale[index++]                                = 0.0;
    _pScale[index++]                                = 1.0;
    _pScale[index++]                                = 1.0;

    index                                           = 0;
    _uScale[index++]                                = 1.0;
    _uScale[index++]                                = 1.0;
    _uScale[index++]                                = 1.0;
    _uScale[index++]                                = 1.0;
    _uScale[index++]                                = 1.0;

    return;
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

void AmoebaReferenceMultipoleForce::setupScaleMaps( const std::vector< std::vector< std::vector<int> > >& multipoleParticleCovalentInfo ){

    /* Setup for scaling maps:
     *
     *     _scaleMaps[particleIndex][ScaleType] = map, where map[covalentIndex] = scaleFactor 
     *     _maxScaleIndex[particleIndex]        = max covalent index for particleIndex
     *
     *     multipoleParticleCovalentInfo[ii][jj], jj =0,1,2,3 contains covalent indices (c12, c13, c14, c15)
     *     multipoleParticleCovalentInfo[ii][jj], jj =4,5,6,7 contains covalent indices (p11, p12, p13, p14)
     *
     *     only including covalent particles w/ index >= ii
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

    //showScaleMapForParticle( 2, stderr );
    //showScaleMapForParticle( 10, stderr );

    return;
}

void AmoebaReferenceMultipoleForce::showScaleMapForParticle( unsigned int particleI, FILE* log ) const {

#ifdef AMOEBA_DEBUG
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
#endif
}

RealOpenMM AmoebaReferenceMultipoleForce::getMultipoleScaleFactor( unsigned int particleI, unsigned int particleJ, ScaleType scaleType ) const {

    MapIntRealOpenMM  scaleMap   = _scaleMaps[particleI][scaleType];
    MapIntRealOpenMMCI isPresent = scaleMap.find( particleJ );
    if( isPresent != scaleMap.end() ){
        return isPresent->second;
    } else {
        return 1.0;
    }
}

void AmoebaReferenceMultipoleForce::getDScaleAndPScale( unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale ) const {
    dScale = getMultipoleScaleFactor( particleI, particleJ, D_SCALE );
    pScale = getMultipoleScaleFactor( particleI, particleJ, P_SCALE );
}

void AmoebaReferenceMultipoleForce::getMultipoleScaleFactors( unsigned int particleI, unsigned int particleJ, std::vector<RealOpenMM>& scaleFactors ) const {
    scaleFactors[D_SCALE] = getMultipoleScaleFactor( particleI, particleJ, D_SCALE );
    scaleFactors[P_SCALE] = getMultipoleScaleFactor( particleI, particleJ, P_SCALE );
    scaleFactors[M_SCALE] = getMultipoleScaleFactor( particleI, particleJ, M_SCALE );
    scaleFactors[U_SCALE] = getMultipoleScaleFactor( particleI, particleJ, U_SCALE );
}

RealOpenMM AmoebaReferenceMultipoleForce::normalizeRealVec( RealVec& vectorToNormalize ) const {
    RealOpenMM norm = SQRT( vectorToNormalize.dot( vectorToNormalize ) );
    if( norm > 0.0 ){
        vectorToNormalize *= (1.0/norm);
    }
    return norm;
}

void AmoebaReferenceMultipoleForce::initializeRealOpenMMVector( vector<RealOpenMM>& vectorToInitialize ) const {
    RealOpenMM zero = 0.0;
    vectorToInitialize.resize( _numParticles );
    std::fill( vectorToInitialize.begin(), vectorToInitialize.end(), zero );
}

void AmoebaReferenceMultipoleForce::initializeRealVecVector( vector<RealVec>& vectorToInitialize ) const {
    vectorToInitialize.resize( _numParticles );
    RealVec zeroVec( 0.0, 0.0, 0.0 );
    std::fill( vectorToInitialize.begin(), vectorToInitialize.end(), zeroVec );
}

void AmoebaReferenceMultipoleForce::copyRealVecVector( const std::vector<OpenMM::RealVec>& inputVector, std::vector<OpenMM::RealVec>& outputVector ) const { 
    outputVector.resize( inputVector.size() );
    for( unsigned int ii = 0; ii < inputVector.size(); ii++ ){
        outputVector[ii] = inputVector[ii];
    }   
    return;
}

void AmoebaReferenceMultipoleForce::loadParticleData( const std::vector<RealVec>& particlePositions,
                                                      const std::vector<RealOpenMM>& charges,
                                                      const std::vector<RealOpenMM>& dipoles,
                                                      const std::vector<RealOpenMM>& quadrupoles,
                                                      const std::vector<RealOpenMM>& tholes,
                                                      const std::vector<RealOpenMM>& dampingFactors,
                                                      const std::vector<RealOpenMM>& polarity,
                                                      std::vector<MultipoleParticleData>& particleData ) const {
   
    for( unsigned int ii = 0; ii < _numParticles; ii++ ){

        particleData[ii].particleIndex        = ii;

        particleData[ii].position             = particlePositions[ii];
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

    }
}

void AmoebaReferenceMultipoleForce::zeroFixed_E_Fields( void ){
    initializeRealVecVector( _fixed_E_Field );
    initializeRealVecVector( _fixed_E_FieldPolar );
}

void AmoebaReferenceMultipoleForce::checkChiralCenterAtParticle( MultipoleParticleData& particleI, int axisType,
                                                                 MultipoleParticleData& particleZ, MultipoleParticleData& particleX, 
                                                                 MultipoleParticleData& particleY ) const {

    if( axisType == AmoebaMultipoleForce::ZThenX ){
        return;
    }
        
    RealVec deltaAD   = particleI.position - particleY.position;
    RealVec deltaBD   = particleZ.position - particleY.position;
    RealVec deltaCD   = particleX.position - particleY.position;

    RealVec deltaC    = deltaBD.cross( deltaCD );
    RealOpenMM volume = deltaC.dot( deltaAD );

    if( volume < 0.0 ){
        particleI.dipole[1]         *= -1.0; // pole(3,i)
        particleI.quadrupole[QXY]   *= -1.0; // pole(6,i)  && pole(8,i)
        particleI.quadrupole[QYZ]   *= -1.0; // pole(10,i) && pole(12,i)
    }
    return;

}

void AmoebaReferenceMultipoleForce::checkChiral( std::vector<MultipoleParticleData>& particleData,
                                                 const std::vector<int>& multipoleAtomXs,
                                                 const std::vector<int>& multipoleAtomYs,
                                                 const std::vector<int>& multipoleAtomZs,
                                                 const std::vector<int>& axisTypes ) const {

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        if( multipoleAtomYs[ii] > -1 ){
            checkChiralCenterAtParticle( particleData[ii], axisTypes[ii], particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]], particleData[multipoleAtomYs[ii]] );
        }
    }
    return;
}

void AmoebaReferenceMultipoleForce::applyRotationMatrixToParticle(       MultipoleParticleData& particleI,
                                                                   const MultipoleParticleData& particleZ,
                                                                   const MultipoleParticleData& particleX,
                                                                         MultipoleParticleData* particleY,
                                                                         int axisType ) const {

    // handle case where rotation matrix is identity (e.g. single ion)

    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom

    RealVec vectorY;
    RealVec vectorZ = particleZ.position - particleI.position;
    RealVec vectorX = particleX.position - particleI.position;

    normalizeRealVec( vectorZ );
 
    // branch based on axis type
 
    if( axisType == AmoebaMultipoleForce::Bisector ){
 
        // bisector
  
        // dx = dx1 + dx2 (in TINKER code)
       
        normalizeRealVec( vectorX );
        vectorZ      += vectorX;
        normalizeRealVec( vectorZ );
       
    } else if( axisType == AmoebaMultipoleForce::ZBisect ){
 
        // z-bisect
  
        // dx = dx1 + dx2 (in TINKER code)
       
        normalizeRealVec( vectorX );

        vectorY  = particleY->position - particleI.position;
        normalizeRealVec( vectorY );

        vectorX += vectorY;
        normalizeRealVec( vectorX );
       
    } else if( axisType == AmoebaMultipoleForce::ThreeFold ){
 
        // 3-fold
  
        // dx = dx1 + dx2 + dx3 (in TINKER code)
       
        normalizeRealVec( vectorX );

        vectorY   = particleY->position - particleI.position;
        normalizeRealVec( vectorY );

        vectorZ  += vectorX +  vectorY;
        normalizeRealVec( vectorZ );
       
    } else if( axisType == AmoebaMultipoleForce::ZOnly ){
 
        // z-only
  
        vectorX = RealVec( 0.1, 0.1, 0.1 );

    }
 
    RealOpenMM dot      = vectorZ.dot( vectorX );
    vectorX            -= vectorZ*dot;

    normalizeRealVec( vectorX );
    vectorY = vectorZ.cross( vectorX );

    RealVec rotationMatrix[3];
    rotationMatrix[0] = vectorX;
    rotationMatrix[1] = vectorY;
    rotationMatrix[2] = vectorZ; 

    RealVec labDipole;
    for( int ii = 0; ii < 3; ii++ ){
        labDipole[ii] = particleI.dipole[0]*rotationMatrix[0][ii];
        for( int jj = 1; jj < 3; jj++ ){
            labDipole[ii] += particleI.dipole[jj]*rotationMatrix[jj][ii];
        }
    }
    particleI.dipole = labDipole;
 
    RealOpenMM mPole[3][3];
    RealOpenMM rPole[3][3] = { { 0.0, 0.0, 0.0 },
                               { 0.0, 0.0, 0.0 },
                               { 0.0, 0.0, 0.0 } };

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

void AmoebaReferenceMultipoleForce::applyRotationMatrix( std::vector<MultipoleParticleData>& particleData,
                                                         const std::vector<int>& multipoleAtomXs,
                                                         const std::vector<int>& multipoleAtomYs,
                                                         const std::vector<int>& multipoleAtomZs,
                                                         const std::vector<int>& axisTypes ) const {
    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        if( multipoleAtomZs[ii] >= 0 && multipoleAtomXs[ii] >= 0 ){
            applyRotationMatrixToParticle( particleData[ii], particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]],
                                           multipoleAtomYs[ii] > -1 ? &particleData[multipoleAtomYs[ii]] : NULL, axisTypes[ii] );
        }
    }

    return;
}

void AmoebaReferenceMultipoleForce::getAndScaleInverseRs( RealOpenMM dampI, RealOpenMM dampJ,
                                                          RealOpenMM tholeI, RealOpenMM tholeJ,
                                                          RealOpenMM r, std::vector<RealOpenMM>& rrI ) const {

    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM one         = 1.0;
    static const RealOpenMM two         = 2.0;
    static const RealOpenMM three       = 3.0;
    static const RealOpenMM five        = 5.0;
    static const RealOpenMM fifty       = 50.0;

    // ---------------------------------------------------------------------------------------
 
    RealOpenMM rI             =  one/r;
    RealOpenMM r2I            =  rI*rI;
 
    rrI[0]                    = rI*r2I;
    RealOpenMM constantFactor = three;
    for( unsigned int ii  = 1; ii < rrI.size(); ii++ ){ 
       rrI[ii]         = constantFactor*rrI[ii-1]*r2I;
       constantFactor += two;
    }
 
    RealOpenMM damp      = dampI*dampJ;
    if( damp != 0.0 ){
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
    return;
}

void AmoebaReferenceMultipoleForce::calculateFixedEFieldPairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                                                 RealOpenMM dScale, RealOpenMM pScale ) {

    // ---------------------------------------------------------------------------------------
 
    static const RealOpenMM zero        = 0.0;
    static const RealOpenMM one         = 1.0;
    static const RealOpenMM two         = 2.0;

    // ---------------------------------------------------------------------------------------
 
    if( particleI.particleIndex == particleJ.particleIndex )return;

    RealVec deltaR    = particleJ.position - particleI.position;
    RealOpenMM r      = SQRT( deltaR.dot( deltaR ) );
    std::vector<RealOpenMM> rrI(3);
 
    // get scaling factors, if needed
  
    getAndScaleInverseRs( particleI.dampingFactor, particleJ.dampingFactor, particleI.thole, particleJ.thole, r, rrI );

    RealOpenMM rr3    = rrI[0];
    RealOpenMM rr5    = rrI[1];
    RealOpenMM rr7    = rrI[2];
    RealOpenMM rr5_2  = rr5*two;

    // field at particle I due multipoles at particle J

    RealVec qDotDelta;
    qDotDelta[0]                          = deltaR[0]*particleJ.quadrupole[QXX] + deltaR[1]*particleJ.quadrupole[QXY] + deltaR[2]*particleJ.quadrupole[QXZ];
    qDotDelta[1]                          = deltaR[0]*particleJ.quadrupole[QXY] + deltaR[1]*particleJ.quadrupole[QYY] + deltaR[2]*particleJ.quadrupole[QYZ];
    qDotDelta[2]                          = deltaR[0]*particleJ.quadrupole[QXZ] + deltaR[1]*particleJ.quadrupole[QYZ] + deltaR[2]*particleJ.quadrupole[QZZ];

    RealOpenMM dipoleDelta                = particleJ.dipole.dot( deltaR ); 
    RealOpenMM qdpoleDelta                = qDotDelta.dot( deltaR ); 
    RealOpenMM factor                     = rr3*particleJ.charge - rr5*dipoleDelta + rr7*qdpoleDelta;

    RealVec field                         = deltaR*factor + particleJ.dipole*rr3 - qDotDelta*rr5_2;

    unsigned int particleIndex            = particleI.particleIndex;
    _fixed_E_Field[particleIndex]        -= field*dScale;
    _fixed_E_FieldPolar[particleIndex]   -= field*pScale;
 
    // field at particle J due multipoles at particle I

    qDotDelta[0]                          = deltaR[0]*particleI.quadrupole[QXX] + deltaR[1]*particleI.quadrupole[QXY] + deltaR[2]*particleI.quadrupole[QXZ];
    qDotDelta[1]                          = deltaR[0]*particleI.quadrupole[QXY] + deltaR[1]*particleI.quadrupole[QYY] + deltaR[2]*particleI.quadrupole[QYZ];
    qDotDelta[2]                          = deltaR[0]*particleI.quadrupole[QXZ] + deltaR[1]*particleI.quadrupole[QYZ] + deltaR[2]*particleI.quadrupole[QZZ];

    dipoleDelta                           = particleI.dipole.dot( deltaR ); 
    qdpoleDelta                           = qDotDelta.dot( deltaR ); 
    factor                                = rr3*particleI.charge + rr5*dipoleDelta + rr7*qdpoleDelta;
 
    field                                 = deltaR*factor - particleI.dipole*rr3 - qDotDelta*rr5_2;
    particleIndex                         = particleJ.particleIndex;
    _fixed_E_Field[particleIndex]        += field*dScale;
    _fixed_E_FieldPolar[particleIndex]   += field*pScale;
 
    return;
}

void AmoebaReferenceMultipoleForce::calculateFixedEField( vector<MultipoleParticleData>& particleData ){

    // calculate fixed electric fields

    zeroFixed_E_Fields();
    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        for( unsigned int jj = ii; jj < _numParticles; jj++ ){

            // field[0]: field at site ii due site jj's charge/dipole/quadrupole
            // field[1]: field at site jj due site ii's charge/dipole/quadrupole

            // if site jj is less than max covalent scaling index then get/apply scaling constants
            // otherwise add unmodified field and fieldPolar to particle fields 

            RealOpenMM dScale, pScale;
            if( jj <= _maxScaleIndex[ii] ){
                getDScaleAndPScale( ii, jj, dScale, pScale );
            } else {
                dScale = pScale = 1.0;
            }
            calculateFixedEFieldPairIxn( particleData[ii], particleData[jj], dScale, pScale );
        }
    }
    return;
}

void AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxn( unsigned int particleI, 
                                                                   unsigned int particleJ,
                                                                   RealOpenMM rr3,
                                                                   RealOpenMM rr5,
                                                                   const RealVec& deltaR,
                                                                   const std::vector<RealVec>& inducedDipole,
                                                                   std::vector<RealVec>& field ) const {

    RealOpenMM dDotDelta            = rr5*(inducedDipole[particleJ].dot( deltaR ) );
    field[particleI]               += inducedDipole[particleJ]*rr3 + deltaR*dDotDelta;
    dDotDelta                       = rr5*(inducedDipole[particleI].dot( deltaR ) );
    field[particleJ]               += inducedDipole[particleI]*rr3 + deltaR*dDotDelta;

}

void AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxns( const MultipoleParticleData& particleI, 
                                                                    const MultipoleParticleData& particleJ,
                                                                    std::vector<UpdateInducedDipoleField>& updateInducedDipoleFields ){

    if( particleI.particleIndex == particleJ.particleIndex )return;

    RealVec deltaR       = particleJ.position - particleI.position;
    RealOpenMM r         =  SQRT( deltaR.dot( deltaR ) );
    std::vector<RealOpenMM> rrI(2);
  
    getAndScaleInverseRs( particleI.dampingFactor, particleJ.dampingFactor,
                          particleI.thole, particleJ.thole, r, rrI );
 
    RealOpenMM rr3       = -rrI[0];
    RealOpenMM rr5       =  rrI[1];
 
    for( unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++ ){
        calculateInducedDipolePairIxn( particleI.particleIndex, particleJ.particleIndex, rr3, rr5, deltaR,
                                       *(updateInducedDipoleFields[ii].inducedDipoles), updateInducedDipoleFields[ii].inducedDipoleField );
    }
    return;

}

RealOpenMM AmoebaReferenceMultipoleForce::updateInducedDipoleFields( const std::vector<MultipoleParticleData>& particleData,
                                                                     std::vector<UpdateInducedDipoleField>& updateInducedDipoleFields){

    // (1) zero fields
    // (2) calculate induced dipole pair ixns 
    // (3) update induced dipoles based on pair ixns and calculate/return convergence factor, maxEpsilon

    RealVec zeroVec( 0.0, 0.0, 0.0 );
    for( unsigned int ii = 0; ii < updateInducedDipoleFields.size(); ii++ ){
        std::fill( updateInducedDipoleFields[ii].inducedDipoleField.begin(), updateInducedDipoleFields[ii].inducedDipoleField.end(), zeroVec );
    }

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii; jj < particleData.size(); jj++ ){
            calculateInducedDipolePairIxns( particleData[ii], particleData[jj], updateInducedDipoleFields );
        }
    }

    RealOpenMM maxEpsilon = 0.0;
    for( unsigned int kk = 0; kk < updateInducedDipoleFields.size(); kk++ ){
        RealOpenMM epsilon = updateInducedDipole( particleData,
                                                  *(updateInducedDipoleFields[kk].fixed_E_Field),
                                                    updateInducedDipoleFields[kk].inducedDipoleField,
                                                  *(updateInducedDipoleFields[kk].inducedDipoles) );
   
        maxEpsilon = epsilon > maxEpsilon ? epsilon : maxEpsilon;
    }

    return maxEpsilon;
}

RealOpenMM AmoebaReferenceMultipoleForce::updateInducedDipole( const std::vector<MultipoleParticleData>& particleData,
                                                               const std::vector<RealVec>& fixed_E_Field,
                                                               const std::vector<RealVec>& inducedDipoleField,
                                                               std::vector<RealVec>& inducedDipole ){

    RealOpenMM epsilon                    = 0.0;
    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        RealVec    oldValue               = inducedDipole[ii];
        RealVec    newValue               = fixed_E_Field[ii] + inducedDipoleField[ii]*particleData[ii].polarity;
        RealVec delta                     = newValue - oldValue;
        inducedDipole[ii]                 = oldValue + delta*_polarSOR;
        epsilon                          += delta.dot( delta );
    }
    return epsilon;
}

void AmoebaReferenceMultipoleForce::convergeInduceDipoles( const std::vector<MultipoleParticleData>& particleData,
                                                           std::vector<UpdateInducedDipoleField>& updateInducedDipoleField){

    bool done                 = false;
    setMutualInducedDipoleConverged( false );
    unsigned int iteration    = 0;
    RealOpenMM currentEpsilon = 1.0e+50;

    // loop until (1) induced dipoles are converged or
    //            (2) iterations == max iterations or
    //            (3) convergence factor (spsilon) increases 

    while( !done ){

        RealOpenMM epsilon = updateInducedDipoleFields( particleData, updateInducedDipoleField);   
                   epsilon = _polarSOR*_debye*SQRT( epsilon/( static_cast<RealOpenMM>(_numParticles) ) );

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

    return;
}

void AmoebaReferenceMultipoleForce::calculateInducedDipoles( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                             const std::vector<MultipoleParticleData>& particleData ){
    // initialize inducedDipoles
    // if polarization type is direct, then return after initializing; otherwise 
    // converge induced dipoles.

    _inducedDipole.resize( _numParticles );
    _inducedDipolePolar.resize( _numParticles );

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        _fixed_E_Field[ii]      *= particleData[ii].polarity;
        _fixed_E_FieldPolar[ii] *= particleData[ii].polarity; 
        _inducedDipole[ii]       = _fixed_E_Field[ii];
        _inducedDipolePolar[ii]  = _fixed_E_FieldPolar[ii];
    }

    if( polarizationType == AmoebaReferenceMultipoleForce::Direct ){
        setMutualInducedDipoleConverged( true );
        return;
    }

    // UpdateInducedDipoleField contains induced dipole, fixed E-fields and field due to other induced dipoles at each site

    std::vector<UpdateInducedDipoleField> updateInducedDipoleField;
    updateInducedDipoleField.push_back( UpdateInducedDipoleField( &_fixed_E_Field,       &_inducedDipole ) );
    updateInducedDipoleField.push_back( UpdateInducedDipoleField( &_fixed_E_FieldPolar,  &_inducedDipolePolar ) );

    convergeInduceDipoles( particleData, updateInducedDipoleField );

    return;
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateElectrostaticPairIxn( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                         const MultipoleParticleData& particleI,
                                                                         const MultipoleParticleData& particleK,
                                                                         const std::vector<RealOpenMM>& scalingFactors,
                                                                         std::vector<RealVec>& forces,
                                                                         std::vector<RealVec>& torque ) const {

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;
    static const RealOpenMM  two           = 2.0;

    RealOpenMM temp3,temp5,temp7;
    RealOpenMM gl[9],gli[7],glip[7];
    RealOpenMM sc[10],sci[8],scip[8];
    RealOpenMM gf[7],gfi[6],gti[6];
 
    unsigned int iIndex = particleI.particleIndex;
    unsigned int kIndex = particleK.particleIndex;

    RealVec delta       = particleK.position - particleI.position;
    RealOpenMM r2       = delta.dot( delta );

    // set conversion factor, cutoff and switching coefficients

    RealOpenMM f        = _electric/_dielectric;

    // set scale factors for permanent multipole and induced terms

    RealOpenMM pdi      = particleI.dampingFactor;
    RealOpenMM pti      = particleI.thole;

    // apply Thole polarization damping to scale factors

    RealOpenMM r        = SQRT(r2);
    RealOpenMM rr1      = 1.0/r;
    RealOpenMM rr3      = rr1/r2;
    RealOpenMM rr5      = 3.0*rr3/r2;
    RealOpenMM rr7      = 5.0*rr5/r2;
    RealOpenMM rr9      = 7.0*rr7/r2;
    RealOpenMM rr11     = 9.0*rr9/r2;

    RealOpenMM scale3   = 1.0;
    RealOpenMM scale5   = 1.0;
    RealOpenMM scale7   = 1.0;

    RealVec ddsc3( 0.0, 0.0, 0.0 );
    RealVec ddsc5( 0.0, 0.0, 0.0 );
    RealVec ddsc7( 0.0, 0.0, 0.0 );

    RealOpenMM damp = particleI.dampingFactor*particleK.dampingFactor;
    if( damp != zero ){
        RealOpenMM pgamma = particleI.thole < particleK.thole ? particleI.thole : particleK.thole;
        RealOpenMM ratio  = (r/damp);
        damp              = -pgamma * ratio*ratio*ratio;
        if( damp > -50.0){
            RealOpenMM expdamp = EXP(damp);
            scale3 = 1.0 - expdamp;
            scale5 = 1.0 - (1.0-damp)*expdamp;
            scale7 = 1.0 - (1.0-damp+0.6*damp*damp)*expdamp;
            temp3 = -3.0 * damp * expdamp / r2;
            temp5 = -damp;
            temp7 = -0.2 - 0.6*damp;

            ddsc3 = delta*temp3;
            ddsc5 = ddsc3*temp5;
            ddsc7 = ddsc5*temp7;

        }
    }
    RealOpenMM scale3i = scale3*scalingFactors[U_SCALE];
    RealOpenMM scale5i = scale5*scalingFactors[U_SCALE];
    RealOpenMM scale7i = scale7*scalingFactors[U_SCALE];
 
    RealOpenMM dsc3    = scale3*scalingFactors[D_SCALE];
    RealOpenMM dsc5    = scale5*scalingFactors[D_SCALE];
    RealOpenMM dsc7    = scale7*scalingFactors[D_SCALE];
 
    RealOpenMM psc3    = scale3*scalingFactors[P_SCALE];
    RealOpenMM psc5    = scale5*scalingFactors[P_SCALE];
    RealOpenMM psc7    = scale7*scalingFactors[P_SCALE];
 
    // construct necessary auxiliary vectors

    RealVec dixdk      = particleI.dipole.cross( particleK.dipole );
    RealVec dixuk      = particleI.dipole.cross( _inducedDipole[kIndex] );
    RealVec dkxui      = particleK.dipole.cross( _inducedDipole[iIndex] );
    RealVec dixukp     = particleI.dipole.cross( _inducedDipolePolar[kIndex] );
    RealVec dkxuip     = particleK.dipole.cross( _inducedDipolePolar[iIndex] );
    RealVec dixr       = particleI.dipole.cross( delta );
    RealVec dkxr       = particleK.dipole.cross( delta );

    RealVec qir;
    qir[0] = particleI.quadrupole[QXX]*delta[0] + particleI.quadrupole[QXY]*delta[1] + particleI.quadrupole[QXZ]*delta[2];
    qir[1] = particleI.quadrupole[QXY]*delta[0] + particleI.quadrupole[QYY]*delta[1] + particleI.quadrupole[QYZ]*delta[2];
    qir[2] = particleI.quadrupole[QXZ]*delta[0] + particleI.quadrupole[QYZ]*delta[1] + particleI.quadrupole[QZZ]*delta[2];

    RealVec qkr;
    qkr[0] = particleK.quadrupole[QXX]*delta[0] + particleK.quadrupole[QXY]*delta[1] + particleK.quadrupole[QXZ]*delta[2];
    qkr[1] = particleK.quadrupole[QXY]*delta[0] + particleK.quadrupole[QYY]*delta[1] + particleK.quadrupole[QYZ]*delta[2];
    qkr[2] = particleK.quadrupole[QXZ]*delta[0] + particleK.quadrupole[QYZ]*delta[1] + particleK.quadrupole[QZZ]*delta[2];

    RealVec qiqkr;
    qiqkr[0] = particleI.quadrupole[QXX]*qkr[0] + particleI.quadrupole[QXY]*qkr[1] + particleI.quadrupole[QXZ]*qkr[2];
    qiqkr[1] = particleI.quadrupole[QXY]*qkr[0] + particleI.quadrupole[QYY]*qkr[1] + particleI.quadrupole[QYZ]*qkr[2];
    qiqkr[2] = particleI.quadrupole[QXZ]*qkr[0] + particleI.quadrupole[QYZ]*qkr[1] + particleI.quadrupole[QZZ]*qkr[2];

    RealVec qkqir;
    qkqir[0] = particleK.quadrupole[QXX]*qir[0] + particleK.quadrupole[QXY]*qir[1] + particleK.quadrupole[QXZ]*qir[2];
    qkqir[1] = particleK.quadrupole[QXY]*qir[0] + particleK.quadrupole[QYY]*qir[1] + particleK.quadrupole[QYZ]*qir[2];
    qkqir[2] = particleK.quadrupole[QXZ]*qir[0] + particleK.quadrupole[QYZ]*qir[1] + particleK.quadrupole[QZZ]*qir[2];

    RealVec qixqk;
    qixqk[0] = particleI.quadrupole[QXY]*particleK.quadrupole[QXZ] +
               particleI.quadrupole[QYY]*particleK.quadrupole[QYZ] +
               particleI.quadrupole[QYZ]*particleK.quadrupole[QZZ] -
               particleI.quadrupole[QXZ]*particleK.quadrupole[QXY] -
               particleI.quadrupole[QYZ]*particleK.quadrupole[QYY] -
               particleI.quadrupole[QZZ]*particleK.quadrupole[QYZ];

    qixqk[1] = particleI.quadrupole[QXZ]*particleK.quadrupole[QXX] +
               particleI.quadrupole[QYZ]*particleK.quadrupole[QXY] +
               particleI.quadrupole[QZZ]*particleK.quadrupole[QXZ] -
               particleI.quadrupole[QXX]*particleK.quadrupole[QXZ] -
               particleI.quadrupole[QXY]*particleK.quadrupole[QYZ] -
               particleI.quadrupole[QXZ]*particleK.quadrupole[QZZ];

    qixqk[2] = particleI.quadrupole[QXX]*particleK.quadrupole[QXY] +
               particleI.quadrupole[QXY]*particleK.quadrupole[QYY] +
               particleI.quadrupole[QXZ]*particleK.quadrupole[QYZ] -
               particleI.quadrupole[QXY]*particleK.quadrupole[QXX] -
               particleI.quadrupole[QYY]*particleK.quadrupole[QXY] -
               particleI.quadrupole[QYZ]*particleK.quadrupole[QXZ];

    RealVec rxqir   = delta.cross( qir );
    RealVec rxqkr   = delta.cross( qkr );
    RealVec rxqikr  = delta.cross( qiqkr );
    RealVec rxqkir  = delta.cross( qkqir );
    RealVec qkrxqir = qkr.cross( qir );

    RealVec qidk,qkdi;
    qidk[0] = particleI.quadrupole[QXX]*particleK.dipole[0] + particleI.quadrupole[QXY]*particleK.dipole[1] + particleI.quadrupole[QXZ]*particleK.dipole[2];
    qidk[1] = particleI.quadrupole[QXY]*particleK.dipole[0] + particleI.quadrupole[QYY]*particleK.dipole[1] + particleI.quadrupole[QYZ]*particleK.dipole[2];
    qidk[2] = particleI.quadrupole[QXZ]*particleK.dipole[0] + particleI.quadrupole[QYZ]*particleK.dipole[1] + particleI.quadrupole[QZZ]*particleK.dipole[2];

    qkdi[0] = particleK.quadrupole[QXX]*particleI.dipole[0] + particleK.quadrupole[QXY]*particleI.dipole[1] + particleK.quadrupole[QXZ]*particleI.dipole[2];
    qkdi[1] = particleK.quadrupole[QXY]*particleI.dipole[0] + particleK.quadrupole[QYY]*particleI.dipole[1] + particleK.quadrupole[QYZ]*particleI.dipole[2];
    qkdi[2] = particleK.quadrupole[QXZ]*particleI.dipole[0] + particleK.quadrupole[QYZ]*particleI.dipole[1] + particleK.quadrupole[QZZ]*particleI.dipole[2];

    RealVec qiuk,qkui;
    qiuk[0] = particleI.quadrupole[QXX]*_inducedDipole[kIndex][0] + particleI.quadrupole[QXY]*_inducedDipole[kIndex][1] + particleI.quadrupole[QXZ]*_inducedDipole[kIndex][2];
    qiuk[1] = particleI.quadrupole[QXY]*_inducedDipole[kIndex][0] + particleI.quadrupole[QYY]*_inducedDipole[kIndex][1] + particleI.quadrupole[QYZ]*_inducedDipole[kIndex][2];
    qiuk[2] = particleI.quadrupole[QXZ]*_inducedDipole[kIndex][0] + particleI.quadrupole[QYZ]*_inducedDipole[kIndex][1] + particleI.quadrupole[QZZ]*_inducedDipole[kIndex][2];

    qkui[0] = particleK.quadrupole[QXX]*_inducedDipole[iIndex][0] + particleK.quadrupole[QXY]*_inducedDipole[iIndex][1] + particleK.quadrupole[QXZ]*_inducedDipole[iIndex][2];
    qkui[1] = particleK.quadrupole[QXY]*_inducedDipole[iIndex][0] + particleK.quadrupole[QYY]*_inducedDipole[iIndex][1] + particleK.quadrupole[QYZ]*_inducedDipole[iIndex][2];
    qkui[2] = particleK.quadrupole[QXZ]*_inducedDipole[iIndex][0] + particleK.quadrupole[QYZ]*_inducedDipole[iIndex][1] + particleK.quadrupole[QZZ]*_inducedDipole[iIndex][2];

    RealVec qiukp,qkuip;
    qiukp[0] = particleI.quadrupole[QXX]*_inducedDipolePolar[kIndex][0] + particleI.quadrupole[QXY]*_inducedDipolePolar[kIndex][1] + particleI.quadrupole[QXZ]*_inducedDipolePolar[kIndex][2];
    qiukp[1] = particleI.quadrupole[QXY]*_inducedDipolePolar[kIndex][0] + particleI.quadrupole[QYY]*_inducedDipolePolar[kIndex][1] + particleI.quadrupole[QYZ]*_inducedDipolePolar[kIndex][2];
    qiukp[2] = particleI.quadrupole[QXZ]*_inducedDipolePolar[kIndex][0] + particleI.quadrupole[QYZ]*_inducedDipolePolar[kIndex][1] + particleI.quadrupole[QZZ]*_inducedDipolePolar[kIndex][2];

    qkuip[0] = particleK.quadrupole[QXX]*_inducedDipolePolar[iIndex][0] + particleK.quadrupole[QXY]*_inducedDipolePolar[iIndex][1] + particleK.quadrupole[QXZ]*_inducedDipolePolar[iIndex][2];
    qkuip[1] = particleK.quadrupole[QXY]*_inducedDipolePolar[iIndex][0] + particleK.quadrupole[QYY]*_inducedDipolePolar[iIndex][1] + particleK.quadrupole[QYZ]*_inducedDipolePolar[iIndex][2];
    qkuip[2] = particleK.quadrupole[QXZ]*_inducedDipolePolar[iIndex][0] + particleK.quadrupole[QYZ]*_inducedDipolePolar[iIndex][1] + particleK.quadrupole[QZZ]*_inducedDipolePolar[iIndex][2];

    RealVec dixqkr   = particleI.dipole.cross( qkr );
    RealVec dkxqir   = particleK.dipole.cross( qir );
    RealVec uixqkr   = _inducedDipole[iIndex].cross( qkr );
    RealVec ukxqir   = _inducedDipole[kIndex].cross( qir );
    RealVec uixqkrp  = _inducedDipolePolar[iIndex].cross( qkr );
    RealVec ukxqirp  = _inducedDipolePolar[kIndex].cross( qir );
    RealVec rxqidk   = delta.cross( qidk );
    RealVec rxqkdi   = delta.cross( qkdi );
    RealVec rxqiuk   = delta.cross( qiuk );
    RealVec rxqkui   = delta.cross( qkui );
    RealVec rxqiukp  = delta.cross( qiukp );
    RealVec rxqkuip  = delta.cross( qkuip );

    // calculate scalar products for permanent components

    sc[1] = particleI.dipole.dot( particleK.dipole );
    sc[2] = particleI.dipole.dot( delta );
    sc[3] = particleK.dipole.dot( delta );
    sc[4] = qir.dot(delta );
    sc[5] = qkr.dot( delta );
    sc[6] = qir.dot(particleK.dipole );
    sc[7] = qkr.dot( particleI.dipole );
    sc[8] = qir.dot( qkr );
    sc[9] = particleI.quadrupole[QXX]*particleK.quadrupole[QXX] + particleI.quadrupole[QXY]*particleK.quadrupole[QXY] + particleI.quadrupole[QXZ]*particleK.quadrupole[QXZ] +
            particleI.quadrupole[QXY]*particleK.quadrupole[QXY] + particleI.quadrupole[QYY]*particleK.quadrupole[QYY] + particleI.quadrupole[QYZ]*particleK.quadrupole[QYZ] +
            particleI.quadrupole[QXZ]*particleK.quadrupole[QXZ] + particleI.quadrupole[QYZ]*particleK.quadrupole[QYZ] + particleI.quadrupole[QZZ]*particleK.quadrupole[QZZ];
 
    // calculate scalar products for induced components

    sci[0] = _inducedDipole[iIndex][0]*particleK.dipole[0] + _inducedDipole[iIndex][1]*particleK.dipole[1] + _inducedDipole[iIndex][2]*particleK.dipole[2] +particleI.dipole[0]*_inducedDipole[kIndex][0] +particleI.dipole[1]*_inducedDipole[kIndex][1] +particleI.dipole[2]*_inducedDipole[kIndex][2];
    sci[1] = _inducedDipole[iIndex][0]*_inducedDipole[kIndex][0] + _inducedDipole[iIndex][1]*_inducedDipole[kIndex][1] + _inducedDipole[iIndex][2]*_inducedDipole[kIndex][2];
    sci[2] = _inducedDipole[iIndex][0]*delta[0] + _inducedDipole[iIndex][1]*delta[1] + _inducedDipole[iIndex][2]*delta[2];
    sci[3] = _inducedDipole[kIndex][0]*delta[0] + _inducedDipole[kIndex][1]*delta[1] + _inducedDipole[kIndex][2]*delta[2];
    sci[6] = qir[0]*_inducedDipole[kIndex][0] + qir[1]*_inducedDipole[kIndex][1] + qir[2]*_inducedDipole[kIndex][2];
    sci[7] = qkr[0]*_inducedDipole[iIndex][0] + qkr[1]*_inducedDipole[iIndex][1] + qkr[2]*_inducedDipole[iIndex][2];
    scip[0] = _inducedDipolePolar[iIndex][0]*particleK.dipole[0] + _inducedDipolePolar[iIndex][1]*particleK.dipole[1] + _inducedDipolePolar[iIndex][2]*particleK.dipole[2] +particleI.dipole[0]*_inducedDipolePolar[kIndex][0] +particleI.dipole[1]*_inducedDipolePolar[kIndex][1] +particleI.dipole[2]*_inducedDipolePolar[kIndex][2];
    scip[1] = _inducedDipole[iIndex][0]*_inducedDipolePolar[kIndex][0]+_inducedDipole[iIndex][1]*_inducedDipolePolar[kIndex][1] + _inducedDipole[iIndex][2]*_inducedDipolePolar[kIndex][2]+_inducedDipolePolar[iIndex][0]*_inducedDipole[kIndex][0] + _inducedDipolePolar[iIndex][1]*_inducedDipole[kIndex][1]+_inducedDipolePolar[iIndex][2]*_inducedDipole[kIndex][2];
    scip[2] = _inducedDipolePolar[iIndex][0]*delta[0] + _inducedDipolePolar[iIndex][1]*delta[1] + _inducedDipolePolar[iIndex][2]*delta[2];
    scip[3] = _inducedDipolePolar[kIndex][0]*delta[0] + _inducedDipolePolar[kIndex][1]*delta[1] + _inducedDipolePolar[kIndex][2]*delta[2];
    scip[6] = qir[0]*_inducedDipolePolar[kIndex][0] + qir[1]*_inducedDipolePolar[kIndex][1] + qir[2]*_inducedDipolePolar[kIndex][2];
    scip[7] = qkr[0]*_inducedDipolePolar[iIndex][0] + qkr[1]*_inducedDipolePolar[iIndex][1] + qkr[2]*_inducedDipolePolar[iIndex][2];

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

    RealOpenMM energy = rr1*gl[0] + rr3*(gl[1]+gl[6]) + rr5*(gl[2]+gl[7]+gl[8]) + rr7*(gl[3]+gl[5]) + rr9*gl[4];
    energy           *= scalingFactors[M_SCALE];
    energy           += 0.5*(rr3*(gli[0]+gli[5])*psc3 + rr5*(gli[1]+gli[6])*psc5 + rr7*gli[2]*psc7);
    energy           *= f;

    // intermediate variables for the permanent components

    gf[0] = rr3*gl[0] + rr5*(gl[1]+gl[6]) + rr7*(gl[2]+gl[7]+gl[8]) + rr9*(gl[3]+gl[5]) + rr11*gl[4];
    gf[1] = -particleK.charge*rr3 + sc[3]*rr5 - sc[5]*rr7;
    gf[2] =  particleI.charge*rr3 + sc[2]*rr5 + sc[4]*rr7;
    gf[3] = 2.0*rr5;
    gf[4] = 2.0*(-particleK.charge*rr5+sc[3]*rr7-sc[5]*rr9);
    gf[5] = 2.0*(-particleI.charge*rr5-sc[2]*rr7-sc[4]*rr9);
    gf[6] = 4.0*rr7;

    // intermediate variables for the induced components

    gfi[0] = 0.5 * rr5 * ((gli[0]+gli[5])*psc3 + (glip[0]+glip[5])*dsc3 + scip[1]*scale3i)
           + 0.5 * rr7 * ((gli[6]+gli[1])*psc5 + (glip[6]+glip[1])*dsc5 - (sci[2]*scip[3]+scip[2]*sci[3])*scale5i)
           + 0.5 * rr9 * (gli[2]*psc7+glip[2]*dsc7);

    gfi[1] = -rr3*particleK.charge + rr5*sc[3] - rr7*sc[5];
    gfi[2] =  rr3*particleI.charge + rr5*sc[2] + rr7*sc[4];
    gfi[3] = 2.0*rr5;
    gfi[4] = rr7*(sci[3]*psc7+scip[3]*dsc7);
    gfi[5] = -rr7*(sci[2]*psc7+scip[2]*dsc7);

    // get the permanent force components

    RealVec ftm2 = delta*gf[0] + particleI.dipole*gf[1] + particleK.dipole*gf[2] + (qkdi -qidk)*gf[3] + qir*gf[4] +
                   qkr*gf[5] + (qiqkr+qkqir)*gf[6];

    // get the induced force components

    RealVec ftm2i  = delta*gfi[0] + qir*gfi[4] + qkr*gfi[5];
            ftm2i += ( 

         (_inducedDipole[iIndex]*psc3 + _inducedDipolePolar[iIndex]*dsc3)*(-rr3*particleK.charge) +
         (_inducedDipole[iIndex]*psc5 + _inducedDipolePolar[iIndex]*dsc5)*( rr5*sc[3]) +
         (_inducedDipole[iIndex]*psc7 + _inducedDipolePolar[iIndex]*dsc7)*(-rr7*sc[5]) +
         (_inducedDipole[kIndex]*psc3 + _inducedDipolePolar[kIndex]*dsc3)*(rr3*particleI.charge) +
         (_inducedDipole[kIndex]*psc5 + _inducedDipolePolar[kIndex]*dsc5)*(rr5*sc[2]) +
         (_inducedDipole[kIndex]*psc7 + _inducedDipolePolar[kIndex]*dsc7)*(rr7*sc[4]) +
         (_inducedDipolePolar[iIndex]*sci[3] + _inducedDipole[iIndex]*scip[3] +
          _inducedDipolePolar[kIndex]*sci[2] + _inducedDipole[kIndex]*scip[2])*(rr5*scale5i) +
         particleI.dipole*((sci[3]*psc5  + scip[3]*dsc5)*rr5) +
         particleK.dipole*((sci[2]*psc5  + scip[2]*dsc5)*rr5) +
         ((qkui - qiuk)*psc5 + (qkuip - qiukp)*dsc5)*(gfi[3]))*0.5;

    // account for partially excluded induced interactions

    temp3 = rr3 * ((gli[0]+gli[5])*scalingFactors[P_SCALE] +(glip[0]+glip[5])*scalingFactors[D_SCALE]);
    temp5 = rr5 * ((gli[1]+gli[6])*scalingFactors[P_SCALE] +(glip[1]+glip[6])*scalingFactors[D_SCALE]);
    temp7 = rr7 * (gli[2]*scalingFactors[P_SCALE] +glip[2]*scalingFactors[D_SCALE]);

    RealVec fridmp,findmp;
    fridmp = (ddsc3*temp3 + ddsc5*temp5 + ddsc7*temp7);

    // find some scaling terms for induced-induced force

    temp3 =  rr3*scalingFactors[U_SCALE]*scip[1];
    temp5 = -rr5*scalingFactors[U_SCALE]*(sci[2]*scip[3]+scip[2]*sci[3]);

    findmp = (ddsc3*temp3 + ddsc5*temp5);

    // modify induced force for partially excluded interactions

    ftm2i -= ( fridmp + findmp )*0.5;

    // correction to convert mutual to direct polarization force

    if( polarizationType == AmoebaReferenceMultipoleForce::Direct ){
       RealOpenMM gfd   = (rr5*scip[1]*scale3i - rr7*(scip[2]*sci[3]+sci[2]*scip[3])*scale5i);
       temp5            = rr5*scale5i;

       RealVec fdir;
       fdir = delta*gfd + (_inducedDipolePolar[iIndex]*sci[3] +
                           _inducedDipole[iIndex]*scip[3] +
                           _inducedDipolePolar[kIndex]*sci[2] +
                           _inducedDipole[kIndex]*scip[2])*temp5;

       ftm2i += ( findmp - fdir )*0.5;
    }

    // intermediate terms for induced torque on multipoles

    gti[1] = 0.5*rr5*(sci[3]*psc5+scip[3]*dsc5);
    gti[2] = 0.5*rr5*(sci[2]*psc5+scip[2]*dsc5);
    gti[3] = gfi[3];
    gti[4] = gfi[4];
    gti[5] = gfi[5];

    // get the permanent torque components

    RealVec ttm2  =  dixdk*(-rr3) + dixr*gf[1] - rxqir*gf[4] +
                     (dixqkr + dkxqir + rxqidk - qixqk*two)*gf[3] -
                     (rxqikr + qkrxqir)*gf[6];

    RealVec ttm3  =  dixdk*rr3 + dkxr*gf[2] - rxqkr*gf[5] -
                     (dixqkr + dkxqir + rxqkdi - qixqk*2.0)*gf[3] -
                     (rxqkir - qkrxqir)*gf[6];

    // get the induced torque components
   
    RealVec ttm2i = (dixuk*psc3 + dixukp*dsc3)*(0.5*(-rr3)) +
                     dixr*gti[1] + 
                     ((ukxqir+rxqiuk)*psc5 + (ukxqirp+rxqiukp)*dsc5)*(0.5*gti[3]) - 
                     rxqir*gti[4];

    RealVec ttm3i = (dkxui*psc3 + dkxuip*dsc3)*(0.5*(-rr3)) +
                     dkxr*gti[2] - 
                    ((uixqkr + rxqkui)*psc5 + (uixqkrp + rxqkuip)*dsc5)*(0.5*gti[3]) - 
                      rxqkr*gti[5];

    // increment forces and torques
    // remove factor of f from torques and add back in?

    RealVec force   = ftm2*scalingFactors[M_SCALE] + ftm2i;
            force  *= f;

    forces[iIndex] -= force;
    forces[kIndex] += force;

    torque[iIndex] += ( ttm2*scalingFactors[M_SCALE] + ttm2i )*f;
    torque[kIndex] += ( ttm3*scalingFactors[M_SCALE] + ttm3i )*f;

    return energy;
}

void AmoebaReferenceMultipoleForce::mapTorqueToForceForParticle( const MultipoleParticleData& particleI,
                                                                 const MultipoleParticleData& particleU,
                                                                 const MultipoleParticleData& particleV,
                                                                       MultipoleParticleData* particleW,
                                                                       int axisType, const Vec3& torque,
                                                                       std::vector<RealVec>& forces ) const {
 
    // ---------------------------------------------------------------------------------------
 
    static const int U                  = 0;
    static const int V                  = 1;
    static const int W                  = 2;
    static const int R                  = 3;
    static const int S                  = 4;
    static const int UV                 = 5;
    static const int UW                 = 6;
    static const int VW                 = 7;
    static const int UR                 = 8;
    static const int US                 = 9;
    static const int VS                 = 10;
    static const int WS                 = 11;
    static const int LastVectorIndex    = 12;
    
    static const int X                  = 0;
    static const int Y                  = 1;
    static const int Z                  = 2;
    static const int I                  = 3;
    
    RealOpenMM norms[LastVectorIndex];
    RealOpenMM angles[LastVectorIndex][2];

    // ---------------------------------------------------------------------------------------
 
    // get coordinates of this atom and the z & x axis atoms
    // compute the vector between the atoms and 1/sqrt(d2), d2 is distance between
    // this atom and the axis atom

    if( axisType == AmoebaMultipoleForce::NoAxisType ){
        return;
    }

    RealVec vectorU = particleU.position - particleI.position;
    norms[U]        = normalizeRealVec( vectorU );

    RealVec vectorV = particleV.position - particleI.position;
    norms[V]        = normalizeRealVec( vectorV );

    RealVec vectorW;
    if( particleW && (axisType == AmoebaMultipoleForce::ZBisect || axisType == AmoebaMultipoleForce::ThreeFold) ){
         vectorW = particleW->position - particleI.position;
    } else {
         vectorW = vectorU.cross( vectorV );
    }
    norms[W]  = normalizeRealVec( vectorW );
 
    RealVec vectorUV, vectorUW, vectorVW;
    vectorUV = vectorV.cross( vectorU );
    vectorUW = vectorW.cross( vectorU );
    vectorVW = vectorW.cross( vectorV );
    
    norms[UV]                     = normalizeRealVec( vectorUV );
    norms[UW]                     = normalizeRealVec( vectorUW );
    norms[VW]                     = normalizeRealVec( vectorVW );

    // angles[][0] is cosine of angle
    // angles[][1] is sine   of angle

    angles[UV][0]                 = vectorU.dot( vectorV );
    angles[UV][1]                 = SQRT( 1.0 - angles[UV][0]*angles[UV][0]);
    
    angles[UW][0]                 = vectorU.dot( vectorW );
    angles[UW][1]                 = SQRT( 1.0 - angles[UW][0]*angles[UW][0]);

    angles[VW][0]                 = vectorV.dot( vectorW );
    angles[VW][1]                 = SQRT( 1.0 - angles[VW][0]*angles[VW][0]);

    RealVec dphi;
    dphi[U]                       = vectorU.dot( torque );
    dphi[V]                       = vectorV.dot( torque );
    dphi[W]                       = vectorW.dot( torque );
    dphi                         *= -1.0;

    // branch based on axis type
 
    if( axisType == AmoebaMultipoleForce::ZThenX || axisType == AmoebaMultipoleForce::Bisector ){
 
        RealOpenMM factor1;
        RealOpenMM factor2;
        RealOpenMM factor3;
        RealOpenMM factor4;
        RealOpenMM half = 0.5;
    
        factor1                 =  dphi[V]/(norms[U]*angles[UV][1]);
        factor2                 =  dphi[W]/(norms[U]);
        factor3                 = -dphi[U]/(norms[V]*angles[UV][1]);
    
        if( axisType == AmoebaMultipoleForce::Bisector ){ 
            factor2    *= half;
            factor4     = half*dphi[W]/(norms[V]);
        } else {
            factor4     = 0.0;
        }
 
        for( int ii = 0; ii < 3; ii++ ){
            double forceU                                        =  vectorUV[ii]*factor1 + factor2*vectorUW[ii];
            forces[particleU.particleIndex][ii]                 -=  forceU;

            double forceV                                        =  vectorUV[ii]*factor3 + factor4*vectorVW[ii];
            forces[particleV.particleIndex][ii]                 -=  forceV;

            forces[particleI.particleIndex][ii]                 +=  (forceU + forceV);
        }

    } else if( axisType == AmoebaMultipoleForce::ZBisect ){

        RealVec vectorR           = vectorV + vectorW; 
        RealVec vectorS           = vectorU.cross( vectorR );

        norms[R]                  = normalizeRealVec( vectorR );
        norms[S]                  = normalizeRealVec( vectorS );

        RealVec vectorUR          =  vectorR.cross( vectorU );
        RealVec vectorUS          =  vectorS.cross( vectorU );
        RealVec vectorVS          =  vectorS.cross( vectorV );
        RealVec vectorWS          =  vectorS.cross( vectorW );

        norms[UR]                 = normalizeRealVec( vectorUR );
        norms[US]                 = normalizeRealVec( vectorUS );
        norms[VS]                 = normalizeRealVec( vectorVS );
        norms[WS]                 = normalizeRealVec( vectorWS );

        angles[UR][0]             = vectorU.dot( vectorR );
        angles[UR][1]             = SQRT( 1.0 - angles[UR][0]*angles[UR][0]);

        angles[US][0]             = vectorU.dot( vectorS );
        angles[US][1]             = SQRT( 1.0 - angles[US][0]*angles[US][0]);

        angles[VS][0]             = vectorV.dot( vectorS );
        angles[VS][1]             = SQRT( 1.0 - angles[VS][0]*angles[VS][0]);

        angles[WS][0]             = vectorW.dot( vectorS );
        angles[WS][1]             = SQRT( 1.0 - angles[WS][0]*angles[WS][0]);
 
        RealVec t1                = vectorV - vectorS*angles[VS][0];
        RealVec t2                = vectorW - vectorS*angles[WS][0];

        RealOpenMM notUsed        = normalizeRealVec( t1 );
              notUsed             = normalizeRealVec( t2 );

        RealOpenMM ut1cos         = vectorU.dot( t1 );
        RealOpenMM ut1sin         = SQRT( 1.0 - ut1cos*ut1cos);

        RealOpenMM ut2cos         = vectorU.dot( t2 );
        RealOpenMM ut2sin         = SQRT( 1.0 - ut2cos*ut2cos);

        RealOpenMM dphiR          = vectorR.dot( torque )*(-1.0);
        RealOpenMM dphiS          = vectorS.dot( torque )*(-1.0);

        RealOpenMM factor1        = dphiR/(norms[U]*angles[UR][1]);
        RealOpenMM factor2        = dphiS/(norms[U]);
        RealOpenMM factor3        = dphi[U]/(norms[V]*(ut1sin+ut2sin));
        RealOpenMM factor4        = dphi[U]/(norms[W]*(ut1sin+ut2sin));

        RealVec forceU            =  vectorUR*factor1 + vectorUS*factor2;
        forces[particleU.particleIndex]        -= forceU;

        RealVec forceV            = (vectorS*angles[VS][1] - t1*angles[VS][0])*factor3;
        forces[particleV.particleIndex]        -= forceV;

        RealVec forceW            = (vectorS*angles[WS][1] - t2*angles[WS][0])*factor4;
        forces[particleW->particleIndex]       -= forceW;

        forces[particleI.particleIndex]        += (forceU + forceV + forceW);

    } else if( axisType == AmoebaMultipoleForce::ThreeFold ){

        // 3-fold

        for( int ii = 0; ii < 3; ii++ ){

            RealOpenMM du =  vectorUW[ii]*dphi[W]/(norms[U]*angles[UW][1]) +
                             vectorUV[ii]*dphi[V]/(norms[U]*angles[UV][1]) -
                             vectorUW[ii]*dphi[U]/(norms[U]*angles[UW][1]) -
                             vectorUV[ii]*dphi[U]/(norms[U]*angles[UV][1]);

            RealOpenMM dv =  vectorVW[ii]*dphi[W]/(norms[V]*angles[VW][1]) -
                             vectorUV[ii]*dphi[U]/(norms[V]*angles[UV][1]) -
                             vectorVW[ii]*dphi[V]/(norms[V]*angles[VW][1]) +
                             vectorUV[ii]*dphi[V]/(norms[V]*angles[UV][1]);

            RealOpenMM dw = -vectorUW[ii]*dphi[U]/(norms[W]*angles[UW][1]) -
                             vectorVW[ii]*dphi[V]/(norms[W]*angles[VW][1]) +
                             vectorUW[ii]*dphi[W]/(norms[W]*angles[UW][1]) +
                             vectorVW[ii]*dphi[W]/(norms[W]*angles[VW][1]);

            du /= 3.0;
            dv /= 3.0;
            dw /= 3.0;

            forces[particleU.particleIndex][ii] -= du;
            forces[particleV.particleIndex][ii] -= dv;
            if( particleW )
                forces[particleW->particleIndex][ii] -= dw;
            forces[particleI.particleIndex][ii] += (du + dv + dw);
        }

    } else if( axisType == AmoebaMultipoleForce::ZOnly ){

        // z-only

        for( int ii = 0; ii < 3; ii++ ){
            RealOpenMM du                               = vectorUV[ii]*dphi[V]/(norms[U]*angles[UV][1]);
            forces[particleU.particleIndex][ii]        -= du;
            forces[particleI.particleIndex][ii]        += du;
        }
    }
 
    return;
 
}

void AmoebaReferenceMultipoleForce::mapTorqueToForce( std::vector<MultipoleParticleData>& particleData,
                                                      const std::vector<int>& multipoleAtomXs,
                                                      const std::vector<int>& multipoleAtomYs,
                                                      const std::vector<int>& multipoleAtomZs,
                                                      const std::vector<int>& axisTypes,
                                                      std::vector<RealVec>& torques,
                                                      std::vector<RealVec>& forces ) const {

    // map torques to forces

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        if( axisTypes[ii] != AmoebaMultipoleForce::NoAxisType ){
             mapTorqueToForceForParticle( particleData[ii],
                                          particleData[multipoleAtomZs[ii]], particleData[multipoleAtomXs[ii]],
                                          multipoleAtomYs[ii] > -1 ? &particleData[multipoleAtomYs[ii]] : NULL,
                                          axisTypes[ii], torques[ii], forces ); 
        }
    }
    return;
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateElectrostatic( std::vector<MultipoleParticleData>& particleData,
                                                                  const std::vector<int>& axisTypes,
                                                                  const std::vector<int>& multipoleAtomZs,
                                                                  const std::vector<int>& multipoleAtomXs,
                                                                  const std::vector<int>& multipoleAtomYs,
                                                                  AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                  std::vector<RealVec>& torques,
                                                                  std::vector<RealVec>& forces ) const {

    // ---------------------------------------------------------------------------------------

    static const RealOpenMM  zero          = 0.0;
    static const RealOpenMM  one           = 1.0;

    // ---------------------------------------------------------------------------------------

    RealOpenMM energy = zero;
    std::vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for( unsigned int kk = 0; kk < scaleFactors.size(); kk++ ){
        scaleFactors[kk] = one;
    }   

    // main loop over particle pairs

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){

            if( jj <= _maxScaleIndex[ii] ){
                getMultipoleScaleFactors( ii, jj, scaleFactors);
            }

            energy += calculateElectrostaticPairIxn( polarizationType, particleData[ii], particleData[jj], scaleFactors, forces, torques );

            if( jj <= _maxScaleIndex[ii] ){
                for( unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++ ){
                    scaleFactors[kk] = one;
                }
            }
        }
    }

    return energy;
}

RealOpenMM AmoebaReferenceMultipoleForce::calculateForceAndEnergy( const std::vector<RealVec>& particlePositions,
                                                                   const std::vector<RealOpenMM>& charges,
                                                                   const std::vector<RealOpenMM>& dipoles,
                                                                   const std::vector<RealOpenMM>& quadrupoles,
                                                                   const std::vector<RealOpenMM>& tholes,
                                                                   const std::vector<RealOpenMM>& dampingFactors,
                                                                   const std::vector<RealOpenMM>& polarity,
                                                                   const std::vector<int>& axisTypes,
                                                                   const std::vector<int>& multipoleAtomZs,
                                                                   const std::vector<int>& multipoleAtomXs,
                                                                   const std::vector<int>& multipoleAtomYs,
                                                                   const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo,
                                                                   AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                   std::vector<RealVec>& forces ){


    _numParticles = particlePositions.size();

    // load particle parameters into vecotr of MultipoleParticleData
    // and check for inverted chiral centers

    std::vector<MultipoleParticleData> particleData( _numParticles );
    loadParticleData( particlePositions, charges, dipoles, quadrupoles,
                      tholes, dampingFactors, polarity, particleData );

    checkChiral( particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes );

    // apply rotation matrix to get lab frame dipole and quadrupoles

    applyRotationMatrix( particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes );

    // calculate fixed electric fields

    setupScaleMaps( multipoleAtomCovalentInfo );
    calculateFixedEField( particleData );

    // get induced dipoles

    calculateInducedDipoles( polarizationType, particleData );

    // check if induced dipoles converged

    if( !getMutualInducedDipoleConverged() ){
        std::stringstream message;
        message << "Induced dipoles did not converge: ";
        message << " iterations="      << getMutualInducedDipoleIterations();
        message << " eps="             << getMutualInducedDipoleEpsilon();
        throw OpenMMException(message.str());
    }

    // calculate electrostatic ixns

    std::vector<RealVec> torques;
    initializeRealVecVector( torques );
    RealOpenMM energy = calculateElectrostatic( particleData, axisTypes, multipoleAtomZs, multipoleAtomXs, multipoleAtomYs, polarizationType, torques, forces );

    // map torques to forces

    mapTorqueToForce( particleData, multipoleAtomXs, multipoleAtomYs, multipoleAtomZs, axisTypes, torques, forces );

    return energy;
}

AmoebaReferenceMultipoleForce::UpdateInducedDipoleField::UpdateInducedDipoleField( std::vector<OpenMM::RealVec>* inputFixed_E_Field, std::vector<OpenMM::RealVec>* inputInducedDipoles ) { 
    fixed_E_Field  = inputFixed_E_Field;
    inducedDipoles = inputInducedDipoles;
    inducedDipoleField.resize( fixed_E_Field->size() );
}   

AmoebaReferenceGeneralizedKirkwoodMultipoleForce::AmoebaReferenceGeneralizedKirkwoodMultipoleForce( AmoebaReferenceGeneralizedKirkwoodForce* amoebaReferenceGeneralizedKirkwoodForce ) :
               AmoebaReferenceMultipoleForce(NoCutoff),
               _amoebaReferenceGeneralizedKirkwoodForce(amoebaReferenceGeneralizedKirkwoodForce),
               _gkc(2.455) {

    const RealOpenMM zero         = 0.0;
    const RealOpenMM one          = 1.0;
    const RealOpenMM two          = 2.0;
    const RealOpenMM three        = 3.0;

    RealOpenMM solventDielectric  = _amoebaReferenceGeneralizedKirkwoodForce->getSolventDielectric();

    _fc                           =       (one - solventDielectric)/(solventDielectric);
    _fd                           =   two*(one - solventDielectric)/(one+two*solventDielectric);
    _fq                           = three*(one - solventDielectric)/(two+three*solventDielectric);

    _amoebaReferenceGeneralizedKirkwoodForce->getGrycukBornRadii( _bornRadii );
    _amoebaReferenceGeneralizedKirkwoodForce->getAtomicRadii( _atomicRadii );
    _amoebaReferenceGeneralizedKirkwoodForce->getScaleFactors( _scaledRadii );

    _includeCavityTerm            = _amoebaReferenceGeneralizedKirkwoodForce->getIncludeCavityTerm();
    _probeRadius                  = _amoebaReferenceGeneralizedKirkwoodForce->getProbeRadius();
    _surfaceAreaFactor            = _amoebaReferenceGeneralizedKirkwoodForce->getSurfaceAreaFactor();
    _dielectricOffset             = _amoebaReferenceGeneralizedKirkwoodForce->getDielectricOffset();

    for( unsigned int ii = 0; ii < _scaledRadii.size(); ii++ ){
        _scaledRadii[ii] *= _atomicRadii[ii];
    }
} 

AmoebaReferenceGeneralizedKirkwoodMultipoleForce::~AmoebaReferenceGeneralizedKirkwoodMultipoleForce( ){
     delete _amoebaReferenceGeneralizedKirkwoodForce;
};
 
int AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getIncludeCavityTerm( void ) const {
     return _includeCavityTerm;
};
 
RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getProbeRadius( void ) const {
     return _probeRadius;
};
 
RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getSurfaceAreaFactor( void ) const {
     return _surfaceAreaFactor;
};
 
RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::getDielectricOffset( void ) const {
     return _dielectricOffset;
};
 
void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::zeroFixed_E_Fields( void ){
    this->AmoebaReferenceMultipoleForce::zeroFixed_E_Fields( );
    initializeRealVecVector( _gkField );
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateFixedEFieldPairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                                                                    RealOpenMM dScale, RealOpenMM pScale ) {
    static const RealOpenMM zero         = 0.0;
    static const RealOpenMM one          = 1.0;
    static const RealOpenMM two          = 2.0;
    static const RealOpenMM three        = 3.0;
    static const RealOpenMM fifteen      = 15.0;
    static const RealOpenMM oneHalf      = 0.5;

    // ---------------------------------------------------------------------------------------
 
    this->AmoebaReferenceMultipoleForce::calculateFixedEFieldPairIxn( particleI, particleJ, dScale, pScale );

    // get deltaR, R2, and R between 2 atoms
 
    RealVec deltaR               = particleJ.position - particleI.position;
    RealOpenMM r                 = SQRT( deltaR.dot( deltaR ) );
    RealOpenMM rb2               = _bornRadii[particleI.particleIndex]*_bornRadii[particleJ.particleIndex];

    RealOpenMM ci                = particleI.charge;

    RealOpenMM uxi               = particleI.dipole[0];
    RealOpenMM uyi               = particleI.dipole[1];
    RealOpenMM uzi               = particleI.dipole[2];
    
    RealOpenMM qxxi              = particleI.quadrupole[QXX];
    RealOpenMM qxyi              = particleI.quadrupole[QXY];
    RealOpenMM qxzi              = particleI.quadrupole[QXZ];
    RealOpenMM qyyi              = particleI.quadrupole[QYY];
    RealOpenMM qyzi              = particleI.quadrupole[QYZ];
    RealOpenMM qzzi              = particleI.quadrupole[QZZ];
    
    RealOpenMM xr                = deltaR[0];
    RealOpenMM yr                = deltaR[1];
    RealOpenMM zr                = deltaR[2];
    
    RealOpenMM ck                = particleJ.charge;
    
    RealOpenMM xr2               = xr*xr;
    RealOpenMM yr2               = yr*yr;
    RealOpenMM zr2               = zr*zr;
    RealOpenMM r2                = xr2 + yr2 + zr2;
    
    RealOpenMM uxk               = particleJ.dipole[0];
    RealOpenMM uyk               = particleJ.dipole[1];
    RealOpenMM uzk               = particleJ.dipole[2];
    
    RealOpenMM qxxk              = particleJ.quadrupole[QXX];
    RealOpenMM qxyk              = particleJ.quadrupole[QXY];
    RealOpenMM qxzk              = particleJ.quadrupole[QXZ];
    RealOpenMM qyyk              = particleJ.quadrupole[QYY];
    RealOpenMM qyzk              = particleJ.quadrupole[QYZ];
    RealOpenMM qzzk              = particleJ.quadrupole[QZZ];
    
    RealOpenMM expterm           = EXP(-r2/(_gkc*rb2));
    RealOpenMM expc              = expterm/_gkc;
    RealOpenMM dexpc             = -two/(_gkc*rb2);
    RealOpenMM gf2               = one/(r2+rb2*expterm);
    RealOpenMM gf                = SQRT(gf2);
    RealOpenMM gf3               = gf2*gf;
    RealOpenMM gf5               = gf3*gf2;
    RealOpenMM gf7               = gf5*gf2;

    // reaction potential auxiliary terms

    RealOpenMM a[4][4];
    a[0][0]                      = gf;
    a[1][0]                      = -gf3;
    a[2][0]                      = three*gf5;
    a[3][0]                      = -fifteen*gf7;

    // reaction potential gradient auxiliary terms

    RealOpenMM expc1             = one - expc;
    a[0][1]                      = expc1*a[1][0];
    a[1][1]                      = expc1*a[2][0];
    a[2][1]                      = expc1*a[3][0];

    // dipole second reaction potential gradient auxiliary term

    RealOpenMM expcdexpc         = -expc*dexpc;
    a[1][2]                      = expc1*a[2][1] + expcdexpc*a[2][0];

    // multiply the auxillary terms by dielectric functions;

    a[0][1]                      = _fc*a[0][1];
    a[1][0]                      = _fd*a[1][0];
    a[1][1]                      = _fd*a[1][1];
    a[1][2]                      = _fd*a[1][2];
    a[2][0]                      = _fq*a[2][0];
    a[2][1]                      = _fq*a[2][1];

    // unweighted dipole reaction potential tensor

    RealOpenMM gux[11],guy[11],guz[11];
    gux[1]                       = xr*a[1][0];
    guy[1]                       = yr*a[1][0];
    guz[1]                       = zr*a[1][0];

    // unweighted reaction potential gradient tensor

    RealOpenMM gc[5];
    gc[2]                        = xr*a[0][1];
    gc[3]                        = yr*a[0][1];
    gc[4]                        = zr*a[0][1];

    gux[2]                       = a[1][0] + xr2*a[1][1];
    gux[3]                       = xr*yr*a[1][1];
    gux[4]                       = xr*zr*a[1][1];
    guy[2]                       = gux[3];
    guy[3]                       = a[1][0] + yr2*a[1][1];
    guy[4]                       = yr*zr*a[1][1];
    guz[2]                       = gux[4];
    guz[3]                       = guy[4];
    guz[4]                       = a[1][0] + zr2*a[1][1];

    RealOpenMM gqxx[5],gqxy[5];
    RealOpenMM gqxz[5],gqyy[5];
    RealOpenMM gqyz[5],gqzz[5];

    gqxx[2]                      = xr*(two*a[2][0]+xr2*a[2][1]);
    gqxx[3]                      = yr*xr2*a[2][1];
    gqxx[4]                      = zr*xr2*a[2][1];

    gqyy[2]                      = xr*yr2*a[2][1];
    gqyy[3]                      = yr*(two*a[2][0]+yr2*a[2][1]);
    gqyy[4]                      = zr*yr2*a[2][1];

    gqzz[2]                      = xr*zr2*a[2][1];
    gqzz[3]                      = yr*zr2*a[2][1];
    gqzz[4]                      = zr*(two*a[2][0]+zr2*a[2][1]);

    gqxy[2]                      = yr*(a[2][0]+xr2*a[2][1]);
    gqxy[3]                      = xr*(a[2][0]+yr2*a[2][1]);
    gqxy[4]                      = zr*xr*yr*a[2][1];

    gqxz[2]                      = zr*(a[2][0]+xr2*a[2][1]);
    gqxz[3]                      = gqxy[4];
    gqxz[4]                      = xr*(a[2][0]+zr2*a[2][1]);

    gqyz[2]                      = gqxy[4];
    gqyz[3]                      = zr*(a[2][0]+yr2*a[2][1]);
    gqyz[4]                      = yr*(a[2][0]+zr2*a[2][1]);

    // unweighted dipole second reaction potential gradient tensor

    gux[5]                       = xr*(three*a[1][1]+xr2*a[1][2]);
    gux[6]                       = yr*(a[1][1]+xr2*a[1][2]);
    gux[7]                       = zr*(a[1][1]+xr2*a[1][2]);

    gux[8]                       = xr*(a[1][1]+yr2*a[1][2]);
    gux[9]                       = zr*xr*yr*a[1][2];
    gux[10]                      = xr*(a[1][1]+zr2*a[1][2]);

    guy[5]                       = yr*(a[1][1]+xr2*a[1][2]);
    guy[6]                       = xr*(a[1][1]+yr2*a[1][2]);
    guy[7]                       = gux[9];

    guy[8]                       = yr*(three*a[1][1]+yr2*a[1][2]);
    guy[9]                       = zr*(a[1][1]+yr2*a[1][2]);
    guy[10]                      = yr*(a[1][1]+zr2*a[1][2]);
    guz[5]                       = zr*(a[1][1]+xr2*a[1][2]);
    guz[6]                       = gux[9];
    guz[7]                       = xr*(a[1][1]+zr2*a[1][2]);
    guz[8]                       = zr*(a[1][1]+yr2*a[1][2]);
    guz[9]                       = yr*(a[1][1]+zr2*a[1][2]);
    guz[10]                      = zr*(three*a[1][1]+zr2*a[1][2]);

    // generalized Kirkwood permanent reaction field

    RealVec fid;
    RealVec fjd;
    fid[0] = uxk*gux[2] + uyk*gux[3] + uzk*gux[4]
                                   + oneHalf*(ck*gux[1] + qxxk*gux[5]
                                   + qyyk*gux[8] + qzzk*gux[10]
                                   + two*(qxyk*gux[6]+qxzk*gux[7]
                                   + qyzk*gux[9]))
                                   + oneHalf*(ck*gc[2] + qxxk*gqxx[2]
                                   + qyyk*gqyy[2] + qzzk*gqzz[2]
                                   + two*(qxyk*gqxy[2]+qxzk*gqxz[2]
                                   + qyzk*gqyz[2]));

    fid[1] = uxk*guy[2] + uyk*guy[3] + uzk*guy[4]
                                   + oneHalf*(ck*guy[1] + qxxk*guy[5]
                                   + qyyk*guy[8] + qzzk*guy[10]
                                   + two*(qxyk*guy[6]+qxzk*guy[7]
                                   + qyzk*guy[9]))
                                   + oneHalf*(ck*gc[3] + qxxk*gqxx[3]
                                   + qyyk*gqyy[3] + qzzk*gqzz[3]
                                   + two*(qxyk*gqxy[3]+qxzk*gqxz[3]
                                   + qyzk*gqyz[3]));

    fid[2] = uxk*guz[2] + uyk*guz[3] + uzk*guz[4]
                                   + oneHalf*(ck*guz[1] + qxxk*guz[5]
                                   + qyyk*guz[8] + qzzk*guz[10]
                                   + two*(qxyk*guz[6]+qxzk*guz[7]
                                   + qyzk*guz[9]))
                                   + oneHalf*(ck*gc[4] + qxxk*gqxx[4]
                                   + qyyk*gqyy[4] + qzzk*gqzz[4]
                                   + two*(qxyk*gqxy[4]+qxzk*gqxz[4]
                                   + qyzk*gqyz[4]));

    fjd[0] = uxi*gux[2] + uyi*gux[3] + uzi*gux[4]
                                   - oneHalf*(ci*gux[1] + qxxi*gux[5]
                                   + qyyi*gux[8] + qzzi*gux[10]
                                   + two*(qxyi*gux[6]+qxzi*gux[7]
                                   + qyzi*gux[9]))
                                   - oneHalf*(ci*gc[2] + qxxi*gqxx[2]
                                   + qyyi*gqyy[2] + qzzi*gqzz[2]
                                   + two*(qxyi*gqxy[2]+qxzi*gqxz[2]
                                   + qyzi*gqyz[2]));

    fjd[1] = uxi*guy[2] + uyi*guy[3] + uzi*guy[4]
                                   - oneHalf*(ci*guy[1] + qxxi*guy[5]
                                   + qyyi*guy[8] + qzzi*guy[10]
                                   + two*(qxyi*guy[6]+qxzi*guy[7]
                                   + qyzi*guy[9]))
                                   - oneHalf*(ci*gc[3]      + qxxi*gqxx[3]
                                   + qyyi*gqyy[3] + qzzi*gqzz[3]
                                   + two*(qxyi*gqxy[3]+qxzi*gqxz[3]
                                   + qyzi*gqyz[3]));

    fjd[2] = uxi*guz[2] + uyi*guz[3] + uzi*guz[4]
                                   - oneHalf*(ci*guz[1] + qxxi*guz[5]
                                   + qyyi*guz[8] + qzzi*guz[10]
                                   + two*(qxyi*guz[6]+qxzi*guz[7]
                                   + qyzi*guz[9]))
                                   - oneHalf*(ci*gc[4] + qxxi*gqxx[4]
                                   + qyyi*gqyy[4] + qzzi*gqzz[4]
                                   + two*(qxyi*gqxy[4]+qxzi*gqxz[4]
                                   + qyzi*gqyz[4]));

    _gkField[particleI.particleIndex] += fid;
    if( particleI.particleIndex != particleJ.particleIndex ){
        _gkField[particleJ.particleIndex] += fjd;
    }
 
    return;
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateInducedDipolePairGkIxn( const MultipoleParticleData& particleI,
                                                                                        const MultipoleParticleData& particleJ,
                                                                                        const std::vector<OpenMM::RealVec>& inputFields,
                                                                                        std::vector<OpenMM::RealVec>& outputFields ) const {

    static const RealOpenMM one    = 1.0;
    static const RealOpenMM three  = 3.0;

    RealOpenMM a[3][3];
    
    RealVec deltaR               = particleJ.position - particleI.position;
    RealOpenMM r                 = SQRT( deltaR.dot( deltaR ) );

    RealOpenMM xr                = deltaR[0];
    RealOpenMM yr                = deltaR[1];
    RealOpenMM zr                = deltaR[2];

    RealOpenMM xr2               = xr*xr;
    RealOpenMM yr2               = yr*yr;
    RealOpenMM zr2               = zr*zr;

    unsigned int iIndex          = particleI.particleIndex;
    unsigned int jIndex          = particleJ.particleIndex;
    RealOpenMM rb2               = _bornRadii[iIndex]*_bornRadii[jIndex];

    RealOpenMM r2                = xr2 + yr2 + zr2;
    RealOpenMM expterm           = EXP(-r2/(_gkc*rb2));
    RealOpenMM expc              = expterm /_gkc; 

    RealOpenMM gf2               = one/(r2+rb2*expterm);
    RealOpenMM gf                = SQRT(gf2);
    RealOpenMM gf3               = gf2*gf;
    RealOpenMM gf5               = gf3*gf2;

    RealVec duis                 = inputFields[iIndex];
    RealVec duks                 = inputFields[jIndex];

    // reaction potential auxiliary terms
 
    a[1][0]                      = -gf3;
    a[2][0]                      = three*gf5;

    // reaction potential gradient auxiliary terms

    RealOpenMM expc1             = one - expc;
    a[1][1]                      = expc1*a[2][0];
 
    // unweighted dipole reaction potential gradient tensor

    RealVec gux, guy, guz;
    gux[0]                       = (a[1][0] + xr2*a[1][1]);
    gux[1]                       = xr*yr*a[1][1];
    gux[2]                       = xr*zr*a[1][1];
 
    guy[0]                       = gux[1];
    guy[1]                       = (a[1][0] + yr2*a[1][1]);
    guy[2]                       = yr*zr*a[1][1];
 
    guz[0]                       = gux[2];
    guz[1]                       = guy[2];
    guz[2]                       = (a[1][0] + zr2*a[1][1]);
 
    outputFields[iIndex][0]     += _fd*duks.dot( gux );
    outputFields[iIndex][1]     += _fd*duks.dot( guy );
    outputFields[iIndex][2]     += _fd*duks.dot( guz );

    // skip i == j, i.e., do not include contribution twice

    if( iIndex !=jIndex ){
        outputFields[jIndex][0] += _fd*duis.dot( gux );
        outputFields[jIndex][1] += _fd*duis.dot( guy );
        outputFields[jIndex][2] += _fd*duis.dot( guz );
    }

    return;
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateInducedDipolePairIxns( const MultipoleParticleData& particleI,
                                                                                       const MultipoleParticleData& particleJ,
                                                                                       std::vector<UpdateInducedDipoleField>& updateInducedDipoleFields ){

    AmoebaReferenceMultipoleForce::calculateInducedDipolePairIxns( particleI, particleJ, updateInducedDipoleFields );

    // include GK contribution

    for( unsigned int ii = 2; ii < updateInducedDipoleFields.size(); ii++ ){
        calculateInducedDipolePairGkIxn( particleI, particleJ, *(updateInducedDipoleFields[ii].inducedDipoles), updateInducedDipoleFields[ii].inducedDipoleField );
    }

    return;
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateInducedDipoles( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                                const std::vector<MultipoleParticleData>& particleData ){

    // initialize inducedDipoles

    _inducedDipole.resize( _numParticles );
    _inducedDipolePolar.resize( _numParticles );
    _inducedDipoleS.resize( _numParticles );
    _inducedDipolePolarS.resize( _numParticles );
    vector<RealVec> gkFieldPolar( _numParticles );

    for( unsigned int ii = 0; ii < _numParticles; ii++ ){

        _fixed_E_Field[ii]            *= particleData[ii].polarity;
        _fixed_E_FieldPolar[ii]       *= particleData[ii].polarity; 
        _gkField[ii]                  *= particleData[ii].polarity; 

        _inducedDipole[ii]             =  _fixed_E_Field[ii];
        _inducedDipolePolar[ii]        =  _fixed_E_FieldPolar[ii];

        _inducedDipoleS[ii]            = (_fixed_E_Field[ii]       + _gkField[ii]);
        _inducedDipolePolarS[ii]       = (_fixed_E_FieldPolar[ii]  + _gkField[ii]);

        _gkField[ii]                   = _inducedDipoleS[ii];
         gkFieldPolar[ii]              = _inducedDipolePolarS[ii];
    }

    if( polarizationType == AmoebaReferenceMultipoleForce::Direct ){
        setMutualInducedDipoleConverged( true );
        return;
    }

    std::vector<UpdateInducedDipoleField> updateInducedDipoleField;
    updateInducedDipoleField.push_back( UpdateInducedDipoleField( &_fixed_E_Field,       &_inducedDipole ) );
    updateInducedDipoleField.push_back( UpdateInducedDipoleField( &_fixed_E_FieldPolar,  &_inducedDipolePolar ) );
    updateInducedDipoleField.push_back( UpdateInducedDipoleField( &_gkField,             &_inducedDipoleS ) );
    updateInducedDipoleField.push_back( UpdateInducedDipoleField( &gkFieldPolar,         &_inducedDipolePolarS ) );

    convergeInduceDipoles( particleData, updateInducedDipoleField );

    return;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateKirkwoodPairIxn( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                                       const MultipoleParticleData& particleI,
                                                                                       const MultipoleParticleData& particleJ,
                                                                                       std::vector<RealVec>& forces,
                                                                                       std::vector<RealVec>& torques,
                                                                                       std::vector<RealOpenMM>& dBorn ) const {

    static const RealOpenMM one      = 1.0;
    static const RealOpenMM two      = 2.0;
    static const RealOpenMM three    = 3.0;
    static const RealOpenMM five     = 5.0;
    static const RealOpenMM six      = 6.0;
    static const RealOpenMM nine     = 9.0;
    static const RealOpenMM twelve   = 12.0;
    static const RealOpenMM oneHalf  = 0.5;

    RealOpenMM e,ei;
    RealOpenMM xr,yr,zr;
    RealOpenMM xr2,yr2,zr2;
    RealOpenMM sxi,syi,szi;
    RealOpenMM sxk,syk,szk;
    RealOpenMM r2,rb2;
    RealOpenMM dedx,dedy,dedz;
    RealOpenMM drbi;
    RealOpenMM drbk;
    RealOpenMM dpdx,dpdy,dpdz;
    RealOpenMM dpbi;
    RealOpenMM dpbk;
    RealOpenMM fc,fd,fq;
    RealOpenMM expterm;
    RealOpenMM gf,gf2,gf3,gf5;
    RealOpenMM gf7,gf9,gf11;
    RealOpenMM expc,dexpc;
    RealOpenMM expc1,expcdexpc;
    RealOpenMM expcr,dexpcr;
    RealOpenMM dgfdr;
    RealOpenMM esym,ewi,ewk;
    RealOpenMM desymdx,dewidx,dewkdx;
    RealOpenMM desymdy,dewidy,dewkdy;
    RealOpenMM desymdz,dewidz,dewkdz;
    RealOpenMM dsumdr,desymdr;
    RealOpenMM dewidr,dewkdr;
    RealOpenMM dsymdr;
    RealOpenMM esymi,ewii,ewki;
    RealOpenMM dpwidx,dpwkdx;
    RealOpenMM dpsymdy,dpwidy,dpwkdy;
    RealOpenMM dpsymdz,dpwidz,dpwkdz;
    RealOpenMM dwipdr,dwkpdr;
    RealOpenMM duvdr;

    unsigned int iIndex   = particleI.particleIndex;
    unsigned int jIndex   = particleJ.particleIndex;

    // set the bulk dielectric constant to the water value

    fc           = _electric*_fc;
    fd           = _electric*_fd;
    fq           = _electric*_fq;

    sxi          = _inducedDipoleS[iIndex][0] + _inducedDipolePolarS[iIndex][0];
    syi          = _inducedDipoleS[iIndex][1] + _inducedDipolePolarS[iIndex][1];
    szi          = _inducedDipoleS[iIndex][2] + _inducedDipolePolarS[iIndex][2];

    // decide whether to compute the current interaction
 
    RealVec deltaR = particleJ.position - particleI.position;
    RealOpenMM r = SQRT( deltaR.dot( deltaR ) ); 

    xr           = deltaR[0];
    yr           = deltaR[1];
    zr           = deltaR[2];

    xr2          = xr*xr;
    yr2          = yr*yr;
    zr2          = zr*zr;
    r2           = xr2 + yr2 + zr2;

    sxk          = _inducedDipoleS[jIndex][0] + _inducedDipolePolarS[jIndex][0];
    syk          = _inducedDipoleS[jIndex][1] + _inducedDipolePolarS[jIndex][1];
    szk          = _inducedDipoleS[jIndex][2] + _inducedDipolePolarS[jIndex][2];

    rb2          = _bornRadii[iIndex]*_bornRadii[jIndex];

    expterm      = EXP(-r2/(_gkc*rb2));
    expc         = expterm/_gkc;
    expcr        = r2*expterm/(_gkc*_gkc*rb2*rb2);
    dexpc        = -two / (_gkc*rb2);
    dexpcr       = two / (_gkc*rb2*rb2);
    dgfdr        = oneHalf * expterm * (one+r2/(rb2*_gkc));
    gf2          = one / (r2+rb2*expterm);
    gf           = SQRT(gf2);
    gf3          = gf2*gf;
    gf5          = gf3*gf2;
    gf7          = gf5*gf2;
    gf9          = gf7*gf2;
    gf11         = gf9*gf2;

    // reaction potential auxiliary terms;

    RealOpenMM a00      =        gf;
    RealOpenMM a10      =       -gf3;
    RealOpenMM a20      =  three*gf5;
    RealOpenMM a30      =  -15.0*gf7;
    RealOpenMM a40      =  105.0*gf9;
    RealOpenMM a50      = -945.0*gf11;

    // Born radii derivatives of reaction potential auxiliary terms;

    RealOpenMM b00      = dgfdr*a10;
    RealOpenMM b10      = dgfdr*a20;
    RealOpenMM b20      = dgfdr*a30;
    RealOpenMM b30      = dgfdr*a40;
    RealOpenMM b40      = dgfdr*a50;

    // reaction potential gradient auxiliary terms;

    expc1               = one - expc;
    RealOpenMM a01      = expc1*a10;
    RealOpenMM a11      = expc1*a20;
    RealOpenMM a21      = expc1*a30;
    RealOpenMM a31      = expc1*a40;
    RealOpenMM a41      = expc1*a50;

    // Born radii derivs of reaction potential gradient auxiliary terms;

    RealOpenMM b01      = b10 - expcr*a10 - expc*b10;
    RealOpenMM b11      = b20 - expcr*a20 - expc*b20;
    RealOpenMM b21      = b30 - expcr*a30 - expc*b30;
    RealOpenMM b31      = b40 - expcr*a40 - expc*b40;

    // 2nd reaction potential gradient auxiliary terms;

    expcdexpc           = -expc*dexpc;
    RealOpenMM a02      = expc1*a11 + expcdexpc*a10;
    RealOpenMM a12      = expc1*a21 + expcdexpc*a20;
    RealOpenMM a22      = expc1*a31 + expcdexpc*a30;
    RealOpenMM a32      = expc1*a41 + expcdexpc*a40;

    // Born radii derivatives of the 2nd reaction potential
    // gradient auxiliary terms

     RealOpenMM b02     = b11 - (expcr*(a11 + dexpc*a10)
                           + expc*(b11 + dexpcr*a10+dexpc*b10));

     RealOpenMM b12     = b21 - (expcr*(a21 + dexpc*a20)
                           +   expc*(b21 + dexpcr*a20+dexpc*b20));

     RealOpenMM b22     = b31 - (expcr*(a31 + dexpc*a30)
                           +   expc*(b31 + dexpcr*a30+dexpc*b30));

    // 3rd reaction potential gradient auxiliary terms

    expcdexpc           = two*expcdexpc;
    RealOpenMM a03      = expc1*a12 + expcdexpc*a11;
    RealOpenMM a13      = expc1*a22 + expcdexpc*a21;
    RealOpenMM a23      = expc1*a32 + expcdexpc*a31;

    expcdexpc      = -expc*dexpc*dexpc;
    a03            = a03 + expcdexpc*a10;
    a13            = a13 + expcdexpc*a20;
    a23            = a23 + expcdexpc*a30;

    // multiply the auxillary terms by their dieletric functions;

    a00            *= fc;
    a01            *= fc;
    a02            *= fc;
    a03            *= fc;

    b00            *= fc;
    b01            *= fc;
    b02            *= fc;

    a10            *= fd;
    a11            *= fd;
    a12            *= fd;
    a13            *= fd;

    b10            *= fd;
    b11            *= fd;
    b12            *= fd;

    a20            *= fq;
    a21            *= fq;
    a22            *= fq;
    a23            *= fq;

    b20            *= fq;
    b21            *= fq;
    b22            *= fq;

    // unweighted reaction potential tensor

    RealOpenMM gc1        = a00;

    RealOpenMM gux1       = xr*a10;
    RealOpenMM guy1       = yr*a10;
    RealOpenMM guz1       = zr*a10;

    RealOpenMM gqxx1      = xr2*a20;
    RealOpenMM gqyy1      = yr2*a20;
    RealOpenMM gqzz1      = zr2*a20;

    RealOpenMM gqxy1      = xr*yr*a20;
    RealOpenMM gqxz1      = xr*zr*a20;
    RealOpenMM gqyz1      = yr*zr*a20;

    // Born radii derivs of unweighted reaction potential tensor

    RealOpenMM gc21       = b00;

    RealOpenMM gux21      = xr*b10;
    RealOpenMM guy21      = yr*b10;
    RealOpenMM guz21      = zr*b10;

    RealOpenMM gqxx21     = xr2*b20;
    RealOpenMM gqyy21     = yr2*b20;
    RealOpenMM gqzz21     = zr2*b20;

    RealOpenMM gqxy21     = xr*yr*b20;
    RealOpenMM gqxz21     = xr*zr*b20;
    RealOpenMM gqyz21     = yr*zr*b20;

    // unweighted reaction potential gradient tensor;

    RealOpenMM gc2        = xr*a01;
    RealOpenMM gc3        = yr*a01;
    RealOpenMM gc4        = zr*a01;

    RealOpenMM gux2       = a10 + xr2*a11;
    RealOpenMM gux3       = xr*yr*a11;
    RealOpenMM gux4       = xr*zr*a11;

    RealOpenMM guy2       = gux3;
    RealOpenMM guy3       = a10 + yr2*a11;
    RealOpenMM guy4       = yr*zr*a11;
    RealOpenMM guz2       = gux4;
    RealOpenMM guz3       = guy4;
    RealOpenMM guz4       = a10 + zr2*a11;
    RealOpenMM gqxx2      = xr*(two*a20+xr2*a21);
    RealOpenMM gqxx3      = yr*xr2*a21;
    RealOpenMM gqxx4      = zr*xr2*a21;
    RealOpenMM gqyy2      = xr*yr2*a21;
    RealOpenMM gqyy3      = yr*(two*a20+yr2*a21);
    RealOpenMM gqyy4      = zr*yr2*a21;
    RealOpenMM gqzz2      = xr*zr2*a21;
    RealOpenMM gqzz3      = yr*zr2*a21;
    RealOpenMM gqzz4      = zr*(two*a20+zr2*a21);
    RealOpenMM gqxy2      = yr*(a20+xr2*a21);
    RealOpenMM gqxy3      = xr*(a20+yr2*a21);
    RealOpenMM gqxy4      = zr*xr*yr*a21;
    RealOpenMM gqxz2      = zr*(a20+xr2*a21);
    RealOpenMM gqxz3      = gqxy4;
    RealOpenMM gqxz4      = xr*(a20+zr2*a21);
    RealOpenMM gqyz2      = gqxy4;
    RealOpenMM gqyz3      = zr*(a20+yr2*a21);
    RealOpenMM gqyz4      = yr*(a20+zr2*a21);

    // Born derivs of the unweighted reaction potential gradient tensor

    RealOpenMM gc22       = xr*b01;
    RealOpenMM gc23       = yr*b01;
    RealOpenMM gc24       = zr*b01;
    RealOpenMM gux22      = b10 + xr2*b11;
    RealOpenMM gux23      = xr * yr * b11;
    RealOpenMM gux24      = xr*zr*b11;
    RealOpenMM guy22      = gux23;
    RealOpenMM guy23      = b10 + yr2*b11;
    RealOpenMM guy24      = yr*zr*b11;
    RealOpenMM guz22      = gux24;
    RealOpenMM guz23      = guy24;
    RealOpenMM guz24      = b10 + zr2*b11;
    RealOpenMM gqxx22     = xr*(two*b20+xr2*b21);
    RealOpenMM gqxx23     = yr*xr2*b21;
    RealOpenMM gqxx24     = zr*xr2*b21;
    RealOpenMM gqyy22     = xr*yr2*b21;
    RealOpenMM gqyy23     = yr*(two*b20+yr2*b21);
    RealOpenMM gqyy24     = zr*yr2*b21;
    RealOpenMM gqzz22     = xr*zr2*b21;
    RealOpenMM gqzz23     = yr*zr2*b21;
    RealOpenMM gqzz24     = zr*(two*b20 + zr2*b21);
    RealOpenMM gqxy22     = yr*(b20+xr2*b21);
    RealOpenMM gqxy23     = xr*(b20+yr2*b21);
    RealOpenMM gqxy24     = zr*xr*yr*b21;
    RealOpenMM gqxz22     = zr*(b20+xr2*b21);
    RealOpenMM gqxz23     = gqxy24;
    RealOpenMM gqxz24     = xr*(b20+zr2*b21);
    RealOpenMM gqyz22     = gqxy24;
    RealOpenMM gqyz23     = zr*(b20+yr2*b21);
    RealOpenMM gqyz24     = yr*(b20+zr2*b21);

    // unweighted 2nd reaction potential gradient tensor;

    RealOpenMM gc5        = a01 + xr2*a02;
    RealOpenMM gc6        = xr*yr*a02;
    RealOpenMM gc7        = xr*zr*a02;
    RealOpenMM gc8        = a01 + yr2*a02;
    RealOpenMM gc9        = yr*zr*a02;
    RealOpenMM gc10       = a01 + zr2*a02;
    RealOpenMM gux5       = xr*(three*a11+xr2*a12);
    RealOpenMM gux6       = yr*(a11+xr2*a12);
    RealOpenMM gux7       = zr*(a11+xr2*a12);
    RealOpenMM gux8       = xr*(a11+yr2*a12);
    RealOpenMM gux9       = zr*xr*yr*a12;
    RealOpenMM gux10      = xr*(a11+zr2*a12);
    RealOpenMM guy5       = yr*(a11+xr2*a12);
    RealOpenMM guy6       = xr*(a11+yr2*a12);
    RealOpenMM guy7       = gux9;
    RealOpenMM guy8       = yr*(three*a11+yr2*a12);
    RealOpenMM guy9       = zr*(a11+yr2*a12);
    RealOpenMM guy10      = yr*(a11+zr2*a12);
    RealOpenMM guz5       = zr*(a11+xr2*a12);
    RealOpenMM guz6       = gux9;
    RealOpenMM guz7       = xr*(a11+zr2*a12);
    RealOpenMM guz8       = zr*(a11+yr2*a12);
    RealOpenMM guz9       = yr*(a11+zr2*a12);
    RealOpenMM guz10      = zr*(three*a11+zr2*a12);
    RealOpenMM gqxx5      = two*a20 + xr2*(five*a21+xr2*a22);
    RealOpenMM gqxx6      = yr*xr*(two*a21+xr2*a22);
    RealOpenMM gqxx7      = zr*xr*(two*a21+xr2*a22);
    RealOpenMM gqxx8      = xr2*(a21+yr2*a22);
    RealOpenMM gqxx9      = zr*yr*xr2*a22;
    RealOpenMM gqxx10     = xr2*(a21+zr2*a22);
    RealOpenMM gqyy5      = yr2*(a21+xr2*a22);
    RealOpenMM gqyy6      = xr*yr*(two*a21+yr2*a22);
    RealOpenMM gqyy7      = xr*zr*yr2*a22;
    RealOpenMM gqyy8      = two*a20 + yr2*(five*a21+yr2*a22);
    RealOpenMM gqyy9      = yr*zr*(two*a21+yr2*a22);
    RealOpenMM gqyy10     = yr2*(a21+zr2*a22);
    RealOpenMM gqzz5      = zr2*(a21+xr2*a22);
    RealOpenMM gqzz6      = xr*yr*zr2*a22;
    RealOpenMM gqzz7      = xr*zr*(two*a21+zr2*a22);
    RealOpenMM gqzz8      = zr2*(a21+yr2*a22);
    RealOpenMM gqzz9      = yr*zr*(two*a21+zr2*a22);
    RealOpenMM gqzz10     = two*a20 + zr2*(five*a21+zr2*a22);
    RealOpenMM gqxy5      = xr*yr*(three*a21+xr2*a22);
    RealOpenMM gqxy6      = a20 + (xr2+yr2)*a21 + xr2*yr2*a22;
    RealOpenMM gqxy7      = zr*yr*(a21+xr2*a22);
    RealOpenMM gqxy8      = xr*yr*(three*a21+yr2*a22);
    RealOpenMM gqxy9      = zr*xr*(a21+yr2*a22);
    RealOpenMM gqxy10     = xr*yr*(a21+zr2*a22);
    RealOpenMM gqxz5      = xr*zr*(three*a21+xr2*a22);
    RealOpenMM gqxz6      = yr*zr*(a21+xr2*a22);
    RealOpenMM gqxz7      = a20 + (xr2+zr2)*a21 + xr2*zr2*a22;
    RealOpenMM gqxz8      = xr*zr*(a21+yr2*a22);
    RealOpenMM gqxz9      = xr*yr*(a21+zr2*a22);
    RealOpenMM gqxz10     = xr*zr*(three*a21+zr2*a22);
    RealOpenMM gqyz5      = zr*yr*(a21+xr2*a22);
    RealOpenMM gqyz6      = xr*zr*(a21+yr2*a22);
    RealOpenMM gqyz7      = xr*yr*(a21+zr2*a22);
    RealOpenMM gqyz8      = yr*zr*(three*a21+yr2*a22);
    RealOpenMM gqyz9      = a20 + (yr2+zr2)*a21 + yr2*zr2*a22;
    RealOpenMM gqyz10     = yr*zr*(three*a21+zr2*a22);

    // Born radii derivatives of the unweighted 2nd reaction;
    // potential gradient tensor;

    RealOpenMM gc25       = b01 + xr2*b02;
    RealOpenMM gc26       = xr*yr*b02;
    RealOpenMM gc27       = xr*zr*b02;
    RealOpenMM gc28       = b01 + yr2*b02;
    RealOpenMM gc29       = yr*zr*b02;
    RealOpenMM gc30       = b01 + zr2*b02;
    RealOpenMM gux25      = xr*(three*b11+xr2*b12);
    RealOpenMM gux26      = yr*(b11+xr2*b12);
    RealOpenMM gux27      = zr*(b11+xr2*b12);
    RealOpenMM gux28      = xr*(b11+yr2*b12);
    RealOpenMM gux29      = zr*xr*yr*b12;
    RealOpenMM gux30      = xr*(b11+zr2*b12);
    RealOpenMM guy25      = yr*(b11+xr2*b12);
    RealOpenMM guy26      = xr*(b11+yr2*b12);
    RealOpenMM guy27      = gux29;
    RealOpenMM guy28      = yr*(three*b11+yr2*b12);
    RealOpenMM guy29      = zr*(b11+yr2*b12);
    RealOpenMM guy30      = yr*(b11+zr2*b12);
    RealOpenMM guz25      = zr*(b11+xr2*b12);
    RealOpenMM guz26      = gux29;
    RealOpenMM guz27      = xr*(b11+zr2*b12);
    RealOpenMM guz28      = zr*(b11+yr2*b12);
    RealOpenMM guz29      = yr*(b11+zr2*b12);
    RealOpenMM guz30      = zr*(three*b11+zr2*b12);
    RealOpenMM gqxx25     = two*b20 + xr2*(five*b21+xr2*b22);
    RealOpenMM gqxx26     = yr*xr*(two*b21+xr2*b22);
    RealOpenMM gqxx27     = zr*xr*(two*b21+xr2*b22);
    RealOpenMM gqxx28     = xr2*(b21+yr2*b22);
    RealOpenMM gqxx29     = zr*yr*xr2*b22;
    RealOpenMM gqxx30     = xr2*(b21+zr2*b22);
    RealOpenMM gqyy25     = yr2*(b21+xr2*b22);
    RealOpenMM gqyy26     = xr*yr*(two*b21+yr2*b22);
    RealOpenMM gqyy27     = xr*zr*yr2*b22;
    RealOpenMM gqyy28     = two*b20 + yr2*(five*b21+yr2*b22);
    RealOpenMM gqyy29     = yr*zr*(two*b21+yr2*b22);
    RealOpenMM gqyy30     = yr2*(b21+zr2*b22);
    RealOpenMM gqzz25     = zr2*(b21+xr2*b22);
    RealOpenMM gqzz26     = xr*yr*zr2*b22;
    RealOpenMM gqzz27     = xr*zr*(two*b21+zr2*b22);
    RealOpenMM gqzz28     = zr2*(b21+yr2*b22);
    RealOpenMM gqzz29     = yr*zr*(two*b21+zr2*b22);
    RealOpenMM gqzz30     = two*b20 + zr2*(five*b21+zr2*b22);
    RealOpenMM gqxy25     = xr*yr*(three*b21 + xr2*b22);
    RealOpenMM gqxy26     = b20 + (xr2+yr2)*b21 + xr2*yr2*b22;
    RealOpenMM gqxy27     = zr*yr*(b21+xr2*b22);
    RealOpenMM gqxy28     = xr*yr*(three*b21+yr2*b22);
    RealOpenMM gqxy29     = zr*xr*(b21+yr2*b22);
    RealOpenMM gqxy30     = xr*yr*(b21+zr2*b22);
    RealOpenMM gqxz25     = xr*zr*(three*b21+xr2*b22);
    RealOpenMM gqxz26     = yr*zr*(b21+xr2*b22);
    RealOpenMM gqxz27     = b20 + (xr2+zr2)*b21 + xr2*zr2*b22;
    RealOpenMM gqxz28     = xr*zr*(b21+yr2*b22);
    RealOpenMM gqxz29     = xr*yr*(b21+zr2*b22);
    RealOpenMM gqxz30     = xr*zr*(three*b21+zr2*b22);
    RealOpenMM gqyz25     = zr*yr*(b21+xr2*b22);
    RealOpenMM gqyz26     = xr*zr*(b21+yr2*b22);
    RealOpenMM gqyz27     = xr*yr*(b21+zr2*b22);
    RealOpenMM gqyz28     = yr*zr*(three*b21+yr2*b22);
    RealOpenMM gqyz29     = b20 + (yr2+zr2)*b21 + yr2*zr2*b22;
    RealOpenMM gqyz30     = yr*zr*(three*b21+zr2*b22);

    // unweighted 3rd reaction potential gradient tensor;

    RealOpenMM gc11       = xr*(three*a02+xr2*a03);
    RealOpenMM gc12       = yr*(a02+xr2*a03);
    RealOpenMM gc13       = zr*(a02+xr2*a03);
    RealOpenMM gc14       = xr*(a02+yr2*a03);
    RealOpenMM gc15       = xr*yr*zr*a03;
    RealOpenMM gc16       = xr*(a02+zr2*a03);
    RealOpenMM gc17       = yr*(three*a02+yr2*a03);
    RealOpenMM gc18       = zr*(a02+yr2*a03);
    RealOpenMM gc19       = yr*(a02+zr2*a03);
    RealOpenMM gc20       = zr*(three*a02+zr2*a03);
    RealOpenMM gux11      = three*a11 + xr2*(six*a12+xr2*a13);
    RealOpenMM gux12      = xr*yr*(three*a12+xr2*a13);
    RealOpenMM gux13      = xr*zr*(three*a12+xr2*a13);
    RealOpenMM gux14      = a11 + (xr2+yr2)*a12 + xr2*yr2*a13;
    RealOpenMM gux15      = yr*zr*(a12+xr2*a13);
    RealOpenMM gux16      = a11 + (xr2+zr2)*a12 + xr2*zr2*a13;
    RealOpenMM gux17      = xr*yr*(three*a12+yr2*a13);
    RealOpenMM gux18      = xr*zr*(a12+yr2*a13);
    RealOpenMM gux19      = xr*yr*(a12+zr2*a13);
    RealOpenMM gux20      = xr*zr*(three*a12+zr2*a13);
    RealOpenMM guy11      = gux12;
    RealOpenMM guy12      = gux14;
    RealOpenMM guy13      = gux15;
    RealOpenMM guy14      = gux17;
    RealOpenMM guy15      = gux18;
    RealOpenMM guy16      = gux19;
    RealOpenMM guy17      = three*a11 + yr2*(six*a12+yr2*a13);
    RealOpenMM guy18      = yr*zr*(three*a12+yr2*a13);
    RealOpenMM guy19      = a11 + (yr2+zr2)*a12 + yr2*zr2*a13;
    RealOpenMM guy20      = yr*zr*(three*a12+zr2*a13);
    RealOpenMM guz11      = gux13;
    RealOpenMM guz12      = gux15;
    RealOpenMM guz13      = gux16;
    RealOpenMM guz14      = gux18;
    RealOpenMM guz15      = gux19;
    RealOpenMM guz16      = gux20;
    RealOpenMM guz17      = guy18;
    RealOpenMM guz18      = guy19;
    RealOpenMM guz19      = guy20;
    RealOpenMM guz20      = three*a11 + zr2*(six*a12+zr2*a13);

    RealOpenMM gqxx11     = xr*(twelve*a21+xr2*(nine*a22 + xr2*a23));
    RealOpenMM gqxx12     = yr*(two*a21+xr2*(five*a22  + xr2*a23));
    RealOpenMM gqxx13     = zr*(two*a21+xr2*(five*a22  + xr2*a23));
    RealOpenMM gqxx14     = xr*(two*a21+yr2*two*a22   +xr2*(a22+yr2*a23));
    RealOpenMM gqxx15     = xr*yr*zr*(two*a22+xr2*a23);
    RealOpenMM gqxx16     = xr*(two*a21+zr2*two*a22 +xr2*(a22+zr2*a23));
    RealOpenMM gqxx17     = yr*xr2*(three*a22+yr2*a23);
    RealOpenMM gqxx18     = zr*xr2*(a22+yr2*a23);
    RealOpenMM gqxx19     = yr*xr2*(a22+zr2*a23);
    RealOpenMM gqxx20     = zr*xr2*(three*a22+zr2*a23);

    RealOpenMM gqxy11     = yr*(three*a21+xr2*(six*a22 +xr2*a23));
    RealOpenMM gqxy12     = xr*(three*(a21+yr2*a22) +xr2*(a22+yr2*a23));
    RealOpenMM gqxy13     = xr*yr*zr*(three*a22+xr2*a23);
    RealOpenMM gqxy14     = yr*(three*(a21+xr2*a22) +yr2*(a22+xr2*a23));
    RealOpenMM gqxy15     = zr*(a21+(yr2+xr2)*a22 +yr2*xr2*a23);
    RealOpenMM gqxy16     = yr*(a21+(xr2+zr2)*a22 +xr2*zr2*a23);
    RealOpenMM gqxy17     = xr*(three*(a21+yr2*a22) +yr2*(three*a22+yr2*a23));
    RealOpenMM gqxy18     = xr*yr*zr*(three*a22+yr2*a23);
    RealOpenMM gqxy19     = xr*(a21+(yr2+zr2)*a22 +yr2*zr2*a23);
    RealOpenMM gqxy20     = xr*yr*zr*(three*a22+zr2*a23);
    RealOpenMM gqxz11     = zr*(three*a21+xr2*(six*a22 +xr2*a23));
    RealOpenMM gqxz12     = xr*yr*zr*(three*a22+xr2*a23);
    RealOpenMM gqxz13     = xr*(three*(a21+zr2*a22) +xr2*(a22+zr2*a23));
    RealOpenMM gqxz14     = zr*(a21+(xr2+yr2)*a22 +xr2*yr2*a23);
    RealOpenMM gqxz15     = yr*(a21+(xr2+zr2)*a22 +zr2*xr2*a23);
    RealOpenMM gqxz16     = zr*(three*(a21+xr2*a22) +zr2*(a22+xr2*a23));
    RealOpenMM gqxz17     = xr*yr*zr*(three*a22+yr2*a23);
    RealOpenMM gqxz18     = xr*(a21+(zr2+yr2)*a22 +zr2*yr2*a23);
    RealOpenMM gqxz19     = xr*yr*zr*(three*a22+zr2*a23);
    RealOpenMM gqxz20     = xr*(three*a21+zr2*(six*a22 +zr2*a23));
    RealOpenMM gqyy11     = xr*yr2*(three*a22+xr2*a23);
    RealOpenMM gqyy12     = yr*(two*a21+xr2*two*a22 +yr2*(a22+xr2*a23));
    RealOpenMM gqyy13     = zr*yr2*(a22+xr2*a23);
    RealOpenMM gqyy14     = xr*(two*a21+yr2*(five*a22 +yr2*a23));
    RealOpenMM gqyy15     = xr*yr*zr*(two*a22+yr2*a23);
    RealOpenMM gqyy16     = xr*yr2*(a22+zr2*a23);
    RealOpenMM gqyy17     = yr*(twelve*a21+yr2*(nine*a22 +yr2*a23));
    RealOpenMM gqyy18     = zr*(two*a21+yr2*(five*a22 +yr2*a23));
    RealOpenMM gqyy19     = yr*(two*a21+zr2*two*a22 +yr2*(a22+zr2*a23));
    RealOpenMM gqyy20     = zr*yr2*(three*a22+zr2*a23);
    RealOpenMM gqyz11     = xr*yr*zr*(three*a22+xr2*a23);
    RealOpenMM gqyz12     = zr*(a21+(xr2+yr2)*a22 +xr2*yr2*a23);
    RealOpenMM gqyz13     = yr*(a21+(xr2+zr2)*a22 +xr2*zr2*a23);
    RealOpenMM gqyz14     = xr*yr*zr*(three*a22+yr2*a23);
    RealOpenMM gqyz15     = xr*(a21+(yr2+zr2)*a22 +yr2*zr2*a23);
    RealOpenMM gqyz16     = xr*yr*zr*(three*a22+zr2*a23);
    RealOpenMM gqyz17     = zr*(three*a21+yr2*(six*a22 +yr2*a23));
    RealOpenMM gqyz18     = yr*(three*(a21+zr2*a22) +yr2*(a22+zr2*a23));
    RealOpenMM gqyz19     = zr*(three*(a21+yr2*a22) +zr2*(a22+yr2*a23));
    RealOpenMM gqyz20     = yr*(three*a21+zr2*(six*a22 +zr2*a23));
    RealOpenMM gqzz11     = xr*zr2*(three*a22+xr2*a23);
    RealOpenMM gqzz12     = yr*(zr2*a22+xr2*(zr2*a23));
    RealOpenMM gqzz13     = zr*(two*a21+xr2*two*a22 +zr2*(a22+xr2*a23));
    RealOpenMM gqzz14     = xr*zr2*(a22+yr2*a23);
    RealOpenMM gqzz15     = xr*yr*zr*(two*a22+zr2*a23);
    RealOpenMM gqzz16     = xr*(two*a21+zr2*(five*a22 +zr2*a23));
    RealOpenMM gqzz17     = yr*zr2*(three*a22+yr2*a23);
    RealOpenMM gqzz18     = zr*(two*a21+yr2*two*a22 +zr2*(a22+yr2*a23));
    RealOpenMM gqzz19     = yr*(two*a21+zr2*(five*a22 +zr2*a23));
    RealOpenMM gqzz20     = zr*(twelve*a21+zr2*(nine*a22 +zr2*a23));

    // electrostatic solvation energy of the permanent multipoles
    // in their own GK reaction potential

    esym = particleI.charge*particleJ.charge*gc1 - (particleI.dipole[0]*(particleJ.dipole[0]*gux2+particleJ.dipole[1]*guy2+particleJ.dipole[2]*guz2)
                                   +  particleI.dipole[1]*(particleJ.dipole[0]*gux3+particleJ.dipole[1]*guy3+particleJ.dipole[2]*guz3)
                                   +  particleI.dipole[2]*(particleJ.dipole[0]*gux4+particleJ.dipole[1]*guy4+particleJ.dipole[2]*guz4));

    ewi =  particleI.charge*(particleJ.dipole[0]*gc2+particleJ.dipole[1]*gc3+particleJ.dipole[2]*gc4)
          -particleJ.charge*(particleI.dipole[0]*gux1+particleI.dipole[1]*guy1+particleI.dipole[2]*guz1)
           +particleI.charge*(particleJ.quadrupole[QXX]*gc5+particleJ.quadrupole[QYY]*gc8+particleJ.quadrupole[QZZ]*gc10
              +two*(particleJ.quadrupole[QXY]*gc6+particleJ.quadrupole[QXZ]*gc7+particleJ.quadrupole[QYZ]*gc9))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx1+particleI.quadrupole[QYY]*gqyy1+particleI.quadrupole[QZZ]*gqzz1
              +two*(particleI.quadrupole[QXY]*gqxy1+particleI.quadrupole[QXZ]*gqxz1+particleI.quadrupole[QYZ]*gqyz1))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux5+particleJ.quadrupole[QYY]*gux8+particleJ.quadrupole[QZZ]*gux10
              +two*(particleJ.quadrupole[QXY]*gux6+particleJ.quadrupole[QXZ]*gux7+particleJ.quadrupole[QYZ]*gux9))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy5+particleJ.quadrupole[QYY]*guy8+particleJ.quadrupole[QZZ]*guy10
              +two*(particleJ.quadrupole[QXY]*guy6+particleJ.quadrupole[QXZ]*guy7+particleJ.quadrupole[QYZ]*guy9))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz5+particleJ.quadrupole[QYY]*guz8+particleJ.quadrupole[QZZ]*guz10
              +two*(particleJ.quadrupole[QXY]*guz6+particleJ.quadrupole[QXZ]*guz7+particleJ.quadrupole[QYZ]*guz9))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
              +two*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
              +two*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
              +two*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqxx8+particleJ.quadrupole[QZZ]*gqxx10
              +two*(particleJ.quadrupole[QXY]*gqxx6+particleJ.quadrupole[QXZ]*gqxx7+particleJ.quadrupole[QYZ]*gqxx9))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy5+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqyy10
              +two*(particleJ.quadrupole[QXY]*gqyy6+particleJ.quadrupole[QXZ]*gqyy7+particleJ.quadrupole[QYZ]*gqyy9))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz5+particleJ.quadrupole[QYY]*gqzz8+particleJ.quadrupole[QZZ]*gqzz10
              +two*(particleJ.quadrupole[QXY]*gqzz6+particleJ.quadrupole[QXZ]*gqzz7+particleJ.quadrupole[QYZ]*gqzz9))
              + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy5+particleJ.quadrupole[QYY]*gqxy8+particleJ.quadrupole[QZZ]*gqxy10
              +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxy7+particleJ.quadrupole[QYZ]*gqxy9))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz5+particleJ.quadrupole[QYY]*gqxz8+particleJ.quadrupole[QZZ]*gqxz10
              +two*(particleJ.quadrupole[QXY]*gqxz6+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqxz9))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz5+particleJ.quadrupole[QYY]*gqyz8+particleJ.quadrupole[QZZ]*gqyz10
              +two*(particleJ.quadrupole[QXY]*gqyz6+particleJ.quadrupole[QXZ]*gqyz7+particleJ.quadrupole[QYZ]*gqyz9)));

    ewk = particleI.charge*(particleJ.dipole[0]*gux1+particleJ.dipole[1]*guy1+particleJ.dipole[2]*guz1)
                      -particleJ.charge*(particleI.dipole[0]*gc2+particleI.dipole[1]*gc3+particleI.dipole[2]*gc4)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx1+particleJ.quadrupole[QYY]*gqyy1+particleJ.quadrupole[QZZ]*gqzz1
              +two*(particleJ.quadrupole[QXY]*gqxy1+particleJ.quadrupole[QXZ]*gqxz1+particleJ.quadrupole[QYZ]*gqyz1))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc5+particleI.quadrupole[QYY]*gc8+particleI.quadrupole[QZZ]*gc10
              +two*(particleI.quadrupole[QXY]*gc6+particleI.quadrupole[QXZ]*gc7+particleI.quadrupole[QYZ]*gc9))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
              +two*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
              +two*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
              +two*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux5+particleI.quadrupole[QYY]*gux8+particleI.quadrupole[QZZ]*gux10
              +two*(particleI.quadrupole[QXY]*gux6+particleI.quadrupole[QXZ]*gux7+particleI.quadrupole[QYZ]*gux9))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy5+particleI.quadrupole[QYY]*guy8+particleI.quadrupole[QZZ]*guy10
              +two*(particleI.quadrupole[QXY]*guy6+particleI.quadrupole[QXZ]*guy7+particleI.quadrupole[QYZ]*guy9))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz5+particleI.quadrupole[QYY]*guz8+particleI.quadrupole[QZZ]*guz10
              +two*(particleI.quadrupole[QXY]*guz6+particleI.quadrupole[QXZ]*guz7+particleI.quadrupole[QYZ]*guz9))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
              +two*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
              +two*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
              +two*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
              +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
              +two*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
              +two*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9)));

    desymdx = particleI.charge*particleJ.charge*gc2 - (particleI.dipole[0]*(particleJ.dipole[0]*gux5+particleJ.dipole[1]*guy5+particleJ.dipole[2]*guz5)
                              +  particleI.dipole[1]*(particleJ.dipole[0]*gux6+particleJ.dipole[1]*guy6+particleJ.dipole[2]*guz6)
                              +  particleI.dipole[2]*(particleJ.dipole[0]*gux7+particleJ.dipole[1]*guy7+particleJ.dipole[2]*guz7));

    dewidx = particleI.charge*(particleJ.dipole[0]*gc5+particleJ.dipole[1]*gc6+particleJ.dipole[2]*gc7)
                      -particleJ.charge*(particleI.dipole[0]*gux2+particleI.dipole[1]*guy2+particleI.dipole[2]*guz2)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc11+particleJ.quadrupole[QYY]*gc14+particleJ.quadrupole[QZZ]*gc16
              +two*(particleJ.quadrupole[QXY]*gc12+particleJ.quadrupole[QXZ]*gc13+particleJ.quadrupole[QYZ]*gc15))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
              +two*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux11+particleJ.quadrupole[QYY]*gux14+particleJ.quadrupole[QZZ]*gux16
              +two*(particleJ.quadrupole[QXY]*gux12+particleJ.quadrupole[QXZ]*gux13+particleJ.quadrupole[QYZ]*gux15))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy11+particleJ.quadrupole[QYY]*guy14+particleJ.quadrupole[QZZ]*guy16
              +two*(particleJ.quadrupole[QXY]*guy12+particleJ.quadrupole[QXZ]*guy13+particleJ.quadrupole[QYZ]*guy15))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz11+particleJ.quadrupole[QYY]*guz14+particleJ.quadrupole[QZZ]*guz16
              +two*(particleJ.quadrupole[QXY]*guz12+particleJ.quadrupole[QXZ]*guz13+particleJ.quadrupole[QYZ]*guz15))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqyy5+particleI.quadrupole[QZZ]*gqzz5
              +two*(particleI.quadrupole[QXY]*gqxy5+particleI.quadrupole[QXZ]*gqxz5+particleI.quadrupole[QYZ]*gqyz5))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
              +two*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
              +two*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx11+particleJ.quadrupole[QYY]*gqxx14+particleJ.quadrupole[QZZ]*gqxx16
              +two*(particleJ.quadrupole[QXY]*gqxx12+particleJ.quadrupole[QXZ]*gqxx13+particleJ.quadrupole[QYZ]*gqxx15))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy11+particleJ.quadrupole[QYY]*gqyy14+particleJ.quadrupole[QZZ]*gqyy16
              +two*(particleJ.quadrupole[QXY]*gqyy12+particleJ.quadrupole[QXZ]*gqyy13+particleJ.quadrupole[QYZ]*gqyy15))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz11+particleJ.quadrupole[QYY]*gqzz14+particleJ.quadrupole[QZZ]*gqzz16
              +two*(particleJ.quadrupole[QXY]*gqzz12+particleJ.quadrupole[QXZ]*gqzz13+particleJ.quadrupole[QYZ]*gqzz15))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy11+particleJ.quadrupole[QYY]*gqxy14+particleJ.quadrupole[QZZ]*gqxy16
              +two*(particleJ.quadrupole[QXY]*gqxy12+particleJ.quadrupole[QXZ]*gqxy13+particleJ.quadrupole[QYZ]*gqxy15))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz11+particleJ.quadrupole[QYY]*gqxz14+particleJ.quadrupole[QZZ]*gqxz16
              +two*(particleJ.quadrupole[QXY]*gqxz12+particleJ.quadrupole[QXZ]*gqxz13+particleJ.quadrupole[QYZ]*gqxz15))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz11+particleJ.quadrupole[QYY]*gqyz14+particleJ.quadrupole[QZZ]*gqyz16
              +two*(particleJ.quadrupole[QXY]*gqyz12+particleJ.quadrupole[QXZ]*gqyz13+particleJ.quadrupole[QYZ]*gqyz15)));

    dewkdx = particleI.charge*(particleJ.dipole[0]*gux2+particleJ.dipole[1]*guy2+particleJ.dipole[2]*guz2)
                      -particleJ.charge*(particleI.dipole[0]*gc5+particleI.dipole[1]*gc6+particleI.dipole[2]*gc7)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
              +two*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc11+particleI.quadrupole[QYY]*gc14+particleI.quadrupole[QZZ]*gc16
              +two*(particleI.quadrupole[QXY]*gc12+particleI.quadrupole[QXZ]*gc13+particleI.quadrupole[QYZ]*gc15))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
              +two*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
              +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
              +two*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux11+particleI.quadrupole[QYY]*gux14+particleI.quadrupole[QZZ]*gux16
              +two*(particleI.quadrupole[QXY]*gux12+particleI.quadrupole[QXZ]*gux13+particleI.quadrupole[QYZ]*gux15))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy11+particleI.quadrupole[QYY]*guy14+particleI.quadrupole[QZZ]*guy16
              +two*(particleI.quadrupole[QXY]*guy12+particleI.quadrupole[QXZ]*guy13+particleI.quadrupole[QYZ]*guy15))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz11+particleI.quadrupole[QYY]*guz14+particleI.quadrupole[QZZ]*guz16
              +two*(particleI.quadrupole[QXY]*guz12+particleI.quadrupole[QXZ]*guz13+particleI.quadrupole[QYZ]*guz15))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx11+particleJ.quadrupole[QYY]*gqyy11+particleJ.quadrupole[QZZ]*gqzz11
              +two*(particleJ.quadrupole[QXY]*gqxy11+particleJ.quadrupole[QXZ]*gqxz11+particleJ.quadrupole[QYZ]*gqyz11))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx14+particleJ.quadrupole[QYY]*gqyy14+particleJ.quadrupole[QZZ]*gqzz14
              +two*(particleJ.quadrupole[QXY]*gqxy14+particleJ.quadrupole[QXZ]*gqxz14+particleJ.quadrupole[QYZ]*gqyz14))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx16+particleJ.quadrupole[QYY]*gqyy16+particleJ.quadrupole[QZZ]*gqzz16
              +two*(particleJ.quadrupole[QXY]*gqxy16+particleJ.quadrupole[QXZ]*gqxz16+particleJ.quadrupole[QYZ]*gqyz16))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx12+particleJ.quadrupole[QYY]*gqyy12+particleJ.quadrupole[QZZ]*gqzz12
              +two*(particleJ.quadrupole[QXY]*gqxy12+particleJ.quadrupole[QXZ]*gqxz12+particleJ.quadrupole[QYZ]*gqyz12))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx13+particleJ.quadrupole[QYY]*gqyy13+particleJ.quadrupole[QZZ]*gqzz13
              +two*(particleJ.quadrupole[QXY]*gqxy13+particleJ.quadrupole[QXZ]*gqxz13+particleJ.quadrupole[QYZ]*gqyz13))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx15+particleJ.quadrupole[QYY]*gqyy15+particleJ.quadrupole[QZZ]*gqzz15
              +two*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqyz15)));

    dedx = desymdx + oneHalf*(dewidx + dewkdx);

    desymdy = particleI.charge*particleJ.charge*gc3
                           - (particleI.dipole[0]*(particleJ.dipole[0]*gux6+particleJ.dipole[1]*guy6+particleJ.dipole[2]*guz6)
                             +particleI.dipole[1]*(particleJ.dipole[0]*gux8+particleJ.dipole[1]*guy8+particleJ.dipole[2]*guz8)
                             +particleI.dipole[2]*(particleJ.dipole[0]*gux9+particleJ.dipole[1]*guy9+particleJ.dipole[2]*guz9));

    dewidy = particleI.charge*(particleJ.dipole[0]*gc6+particleJ.dipole[1]*gc8+particleJ.dipole[2]*gc9)
                      -particleJ.charge*(particleI.dipole[0]*gux3+particleI.dipole[1]*guy3+particleI.dipole[2]*guz3)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc12+particleJ.quadrupole[QYY]*gc17+particleJ.quadrupole[QZZ]*gc19
              +two*(particleJ.quadrupole[QXY]*gc14+particleJ.quadrupole[QXZ]*gc15+particleJ.quadrupole[QYZ]*gc18))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
              +two*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux12+particleJ.quadrupole[QYY]*gux17+particleJ.quadrupole[QZZ]*gux19
              +two*(particleJ.quadrupole[QXY]*gux14+particleJ.quadrupole[QXZ]*gux15+particleJ.quadrupole[QYZ]*gux18))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy12+particleJ.quadrupole[QYY]*guy17+particleJ.quadrupole[QZZ]*guy19
              +two*(particleJ.quadrupole[QXY]*guy14+particleJ.quadrupole[QXZ]*guy15+particleJ.quadrupole[QYZ]*guy18))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz12+particleJ.quadrupole[QYY]*guz17+particleJ.quadrupole[QZZ]*guz19
              +two*(particleJ.quadrupole[QXY]*guz14+particleJ.quadrupole[QXZ]*guz15+particleJ.quadrupole[QYZ]*guz18))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
              +two*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx8+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqzz8
              +two*(particleI.quadrupole[QXY]*gqxy8+particleI.quadrupole[QXZ]*gqxz8+particleI.quadrupole[QYZ]*gqyz8))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
              +two*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx12+particleJ.quadrupole[QYY]*gqxx17+particleJ.quadrupole[QZZ]*gqxx19
              +two*(particleJ.quadrupole[QXY]*gqxx14+particleJ.quadrupole[QXZ]*gqxx15+particleJ.quadrupole[QYZ]*gqxx18))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy12+particleJ.quadrupole[QYY]*gqyy17+particleJ.quadrupole[QZZ]*gqyy19
              +two*(particleJ.quadrupole[QXY]*gqyy14+particleJ.quadrupole[QXZ]*gqyy15+particleJ.quadrupole[QYZ]*gqyy18))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz12+particleJ.quadrupole[QYY]*gqzz17+particleJ.quadrupole[QZZ]*gqzz19
              +two*(particleJ.quadrupole[QXY]*gqzz14+particleJ.quadrupole[QXZ]*gqzz15+particleJ.quadrupole[QYZ]*gqzz18))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy12+particleJ.quadrupole[QYY]*gqxy17+particleJ.quadrupole[QZZ]*gqxy19
              +two*(particleJ.quadrupole[QXY]*gqxy14+particleJ.quadrupole[QXZ]*gqxy15+particleJ.quadrupole[QYZ]*gqxy18))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz12+particleJ.quadrupole[QYY]*gqxz17+particleJ.quadrupole[QZZ]*gqxz19
              +two*(particleJ.quadrupole[QXY]*gqxz14+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqxz18))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz12+particleJ.quadrupole[QYY]*gqyz17+particleJ.quadrupole[QZZ]*gqyz19
              +two*(particleJ.quadrupole[QXY]*gqyz14+particleJ.quadrupole[QXZ]*gqyz15+particleJ.quadrupole[QYZ]*gqyz18)));

    dewkdy = particleI.charge*(particleJ.dipole[0]*gux3+particleJ.dipole[1]*guy3+particleJ.dipole[2]*guz3)
                      -particleJ.charge*(particleI.dipole[0]*gc6+particleI.dipole[1]*gc8+particleI.dipole[2]*gc9)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
              +two*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc12+particleI.quadrupole[QYY]*gc17+particleI.quadrupole[QZZ]*gc19
              +two*(particleI.quadrupole[QXY]*gc14+particleI.quadrupole[QXZ]*gc15+particleI.quadrupole[QYZ]*gc18))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
              +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
              +two*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
              +two*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux12+particleI.quadrupole[QYY]*gux17+particleI.quadrupole[QZZ]*gux19
              +two*(particleI.quadrupole[QXY]*gux14+particleI.quadrupole[QXZ]*gux15+particleI.quadrupole[QYZ]*gux18))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy12+particleI.quadrupole[QYY]*guy17+particleI.quadrupole[QZZ]*guy19
              +two*(particleI.quadrupole[QXY]*guy14+particleI.quadrupole[QXZ]*guy15+particleI.quadrupole[QYZ]*guy18))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz12+particleI.quadrupole[QYY]*guz17+particleI.quadrupole[QZZ]*guz19
              +two*(particleI.quadrupole[QXY]*guz14+particleI.quadrupole[QXZ]*guz15+particleI.quadrupole[QYZ]*guz18))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx12+particleJ.quadrupole[QYY]*gqyy12+particleJ.quadrupole[QZZ]*gqzz12
              +two*(particleJ.quadrupole[QXY]*gqxy12+particleJ.quadrupole[QXZ]*gqxz12+particleJ.quadrupole[QYZ]*gqyz12))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx17+particleJ.quadrupole[QYY]*gqyy17+particleJ.quadrupole[QZZ]*gqzz17
              +two*(particleJ.quadrupole[QXY]*gqxy17+particleJ.quadrupole[QXZ]*gqxz17+particleJ.quadrupole[QYZ]*gqyz17))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx19+particleJ.quadrupole[QYY]*gqyy19+particleJ.quadrupole[QZZ]*gqzz19
              +two*(particleJ.quadrupole[QXY]*gqxy19+particleJ.quadrupole[QXZ]*gqxz19+particleJ.quadrupole[QYZ]*gqyz19))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx14+particleJ.quadrupole[QYY]*gqyy14+particleJ.quadrupole[QZZ]*gqzz14
              +two*(particleJ.quadrupole[QXY]*gqxy14+particleJ.quadrupole[QXZ]*gqxz14+particleJ.quadrupole[QYZ]*gqyz14))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx15+particleJ.quadrupole[QYY]*gqyy15+particleJ.quadrupole[QZZ]*gqzz15
              +two*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqyz15))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx18+particleJ.quadrupole[QYY]*gqyy18+particleJ.quadrupole[QZZ]*gqzz18
              +two*(particleJ.quadrupole[QXY]*gqxy18+particleJ.quadrupole[QXZ]*gqxz18+particleJ.quadrupole[QYZ]*gqyz18)));

    dedy = desymdy + oneHalf*(dewidy + dewkdy);

    desymdz = particleI.charge*particleJ.charge*gc4
                           - (particleI.dipole[0]*(particleJ.dipole[0]*gux7+particleJ.dipole[1]*guy7+particleJ.dipole[2]*guz7)
                             +particleI.dipole[1]*(particleJ.dipole[0]*gux9+particleJ.dipole[1]*guy9+particleJ.dipole[2]*guz9)
                             +particleI.dipole[2]*(particleJ.dipole[0]*gux10+particleJ.dipole[1]*guy10+particleJ.dipole[2]*guz10));

    dewidz = particleI.charge*(particleJ.dipole[0]*gc7+particleJ.dipole[1]*gc9+particleJ.dipole[2]*gc10)
                      -particleJ.charge*(particleI.dipole[0]*gux4+particleI.dipole[1]*guy4+particleI.dipole[2]*guz4)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc13+particleJ.quadrupole[QYY]*gc18+particleJ.quadrupole[QZZ]*gc20
              +two*(particleJ.quadrupole[QXY]*gc15+particleJ.quadrupole[QXZ]*gc16+particleJ.quadrupole[QYZ]*gc19))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
              +two*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux13+particleJ.quadrupole[QYY]*gux18+particleJ.quadrupole[QZZ]*gux20
              +two*(particleJ.quadrupole[QXY]*gux15+particleJ.quadrupole[QXZ]*gux16+particleJ.quadrupole[QYZ]*gux19))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy13+particleJ.quadrupole[QYY]*guy18+particleJ.quadrupole[QZZ]*guy20
              +two*(particleJ.quadrupole[QXY]*guy15+particleJ.quadrupole[QXZ]*guy16+particleJ.quadrupole[QYZ]*guy19))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz13+particleJ.quadrupole[QYY]*guz18+particleJ.quadrupole[QZZ]*guz20
              +two*(particleJ.quadrupole[QXY]*guz15+particleJ.quadrupole[QXZ]*guz16+particleJ.quadrupole[QYZ]*guz19))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
              +two*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
              +two*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx10+particleI.quadrupole[QYY]*gqyy10+particleI.quadrupole[QZZ]*gqzz10
              +two*(particleI.quadrupole[QXY]*gqxy10+particleI.quadrupole[QXZ]*gqxz10+particleI.quadrupole[QYZ]*gqyz10))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx13+particleJ.quadrupole[QYY]*gqxx18+particleJ.quadrupole[QZZ]*gqxx20
              +two*(particleJ.quadrupole[QXY]*gqxx15+particleJ.quadrupole[QXZ]*gqxx16+particleJ.quadrupole[QYZ]*gqxx19))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy13+particleJ.quadrupole[QYY]*gqyy18+particleJ.quadrupole[QZZ]*gqyy20
              +two*(particleJ.quadrupole[QXY]*gqyy15+particleJ.quadrupole[QXZ]*gqyy16+particleJ.quadrupole[QYZ]*gqyy19))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz13+particleJ.quadrupole[QYY]*gqzz18+particleJ.quadrupole[QZZ]*gqzz20
              +two*(particleJ.quadrupole[QXY]*gqzz15+particleJ.quadrupole[QXZ]*gqzz16+particleJ.quadrupole[QYZ]*gqzz19))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy13+particleJ.quadrupole[QYY]*gqxy18+particleJ.quadrupole[QZZ]*gqxy20
              +two*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxy16+particleJ.quadrupole[QYZ]*gqxy19))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz13+particleJ.quadrupole[QYY]*gqxz18+particleJ.quadrupole[QZZ]*gqxz20
              +two*(particleJ.quadrupole[QXY]*gqxz15+particleJ.quadrupole[QXZ]*gqxz16+particleJ.quadrupole[QYZ]*gqxz19))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz13+particleJ.quadrupole[QYY]*gqyz18+particleJ.quadrupole[QZZ]*gqyz20
              +two*(particleJ.quadrupole[QXY]*gqyz15+particleJ.quadrupole[QXZ]*gqyz16+particleJ.quadrupole[QYZ]*gqyz19)));

    dewkdz = particleI.charge*(particleJ.dipole[0]*gux4+particleJ.dipole[1]*guy4+particleJ.dipole[2]*guz4)
                      -particleJ.charge*(particleI.dipole[0]*gc7+particleI.dipole[1]*gc9+particleI.dipole[2]*gc10)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
              +two*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc13+particleI.quadrupole[QYY]*gc18+particleI.quadrupole[QZZ]*gc20
              +two*(particleI.quadrupole[QXY]*gc15+particleI.quadrupole[QXZ]*gc16+particleI.quadrupole[QYZ]*gc19))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
              +two*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
              +two*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
              +two*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux13+particleI.quadrupole[QYY]*gux18+particleI.quadrupole[QZZ]*gux20
              +two*(particleI.quadrupole[QXY]*gux15+particleI.quadrupole[QXZ]*gux16+particleI.quadrupole[QYZ]*gux19))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy13+particleI.quadrupole[QYY]*guy18+particleI.quadrupole[QZZ]*guy20
              +two*(particleI.quadrupole[QXY]*guy15+particleI.quadrupole[QXZ]*guy16+particleI.quadrupole[QYZ]*guy19))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz13+particleI.quadrupole[QYY]*guz18+particleI.quadrupole[QZZ]*guz20
              +two*(particleI.quadrupole[QXY]*guz15+particleI.quadrupole[QXZ]*guz16+particleI.quadrupole[QYZ]*guz19))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx13+particleJ.quadrupole[QYY]*gqyy13+particleJ.quadrupole[QZZ]*gqzz13
              +two*(particleJ.quadrupole[QXY]*gqxy13+particleJ.quadrupole[QXZ]*gqxz13+particleJ.quadrupole[QYZ]*gqyz13))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx18+particleJ.quadrupole[QYY]*gqyy18+particleJ.quadrupole[QZZ]*gqzz18
              +two*(particleJ.quadrupole[QXY]*gqxy18+particleJ.quadrupole[QXZ]*gqxz18+particleJ.quadrupole[QYZ]*gqyz18))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx20+particleJ.quadrupole[QYY]*gqyy20+particleJ.quadrupole[QZZ]*gqzz20
              +two*(particleJ.quadrupole[QXY]*gqxy20+particleJ.quadrupole[QXZ]*gqxz20+particleJ.quadrupole[QYZ]*gqyz20))
       + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx15+particleJ.quadrupole[QYY]*gqyy15+particleJ.quadrupole[QZZ]*gqzz15
              +two*(particleJ.quadrupole[QXY]*gqxy15+particleJ.quadrupole[QXZ]*gqxz15+particleJ.quadrupole[QYZ]*gqyz15))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx16+particleJ.quadrupole[QYY]*gqyy16+particleJ.quadrupole[QZZ]*gqzz16
              +two*(particleJ.quadrupole[QXY]*gqxy16+particleJ.quadrupole[QXZ]*gqxz16+particleJ.quadrupole[QYZ]*gqyz16))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx19+particleJ.quadrupole[QYY]*gqyy19+particleJ.quadrupole[QZZ]*gqzz19
              +two*(particleJ.quadrupole[QXY]*gqxy19+particleJ.quadrupole[QXZ]*gqxz19+particleJ.quadrupole[QYZ]*gqyz19)));

    dedz = desymdz + oneHalf*(dewidz + dewkdz);

    desymdr = particleI.charge*particleJ.charge*gc21
                           - (particleI.dipole[0]*(particleJ.dipole[0]*gux22+particleJ.dipole[1]*guy22+particleJ.dipole[2]*guz22)
                             +particleI.dipole[1]*(particleJ.dipole[0]*gux23+particleJ.dipole[1]*guy23+particleJ.dipole[2]*guz23)
                             +particleI.dipole[2]*(particleJ.dipole[0]*gux24+particleJ.dipole[1]*guy24+particleJ.dipole[2]*guz24));

    dewidr = particleI.charge*(particleJ.dipole[0]*gc22+particleJ.dipole[1]*gc23+particleJ.dipole[2]*gc24)
                      -particleJ.charge*(particleI.dipole[0]*gux21+particleI.dipole[1]*guy21+particleI.dipole[2]*guz21)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gc25+particleJ.quadrupole[QYY]*gc28+particleJ.quadrupole[QZZ]*gc30
              +two*(particleJ.quadrupole[QXY]*gc26+particleJ.quadrupole[QXZ]*gc27+particleJ.quadrupole[QYZ]*gc29))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gqxx21+particleI.quadrupole[QYY]*gqyy21+particleI.quadrupole[QZZ]*gqzz21
              +two*(particleI.quadrupole[QXY]*gqxy21+particleI.quadrupole[QXZ]*gqxz21+particleI.quadrupole[QYZ]*gqyz21))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gux25+particleJ.quadrupole[QYY]*gux28+particleJ.quadrupole[QZZ]*gux30
              +two*(particleJ.quadrupole[QXY]*gux26+particleJ.quadrupole[QXZ]*gux27+particleJ.quadrupole[QYZ]*gux29))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*guy25+particleJ.quadrupole[QYY]*guy28+particleJ.quadrupole[QZZ]*guy30
              +two*(particleJ.quadrupole[QXY]*guy26+particleJ.quadrupole[QXZ]*guy27+particleJ.quadrupole[QYZ]*guy29))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*guz25+particleJ.quadrupole[QYY]*guz28+particleJ.quadrupole[QZZ]*guz30
              +two*(particleJ.quadrupole[QXY]*guz26+particleJ.quadrupole[QXZ]*guz27+particleJ.quadrupole[QYZ]*guz29))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gqxx22+particleI.quadrupole[QYY]*gqyy22+particleI.quadrupole[QZZ]*gqzz22
              +two*(particleI.quadrupole[QXY]*gqxy22+particleI.quadrupole[QXZ]*gqxz22+particleI.quadrupole[QYZ]*gqyz22))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*gqxx23+particleI.quadrupole[QYY]*gqyy23+particleI.quadrupole[QZZ]*gqzz23
              +two*(particleI.quadrupole[QXY]*gqxy23+particleI.quadrupole[QXZ]*gqxz23+particleI.quadrupole[QYZ]*gqyz23))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*gqxx24+particleI.quadrupole[QYY]*gqyy24+particleI.quadrupole[QZZ]*gqzz24
              +two*(particleI.quadrupole[QXY]*gqxy24+particleI.quadrupole[QXZ]*gqxz24+particleI.quadrupole[QYZ]*gqyz24))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx25+particleJ.quadrupole[QYY]*gqxx28+particleJ.quadrupole[QZZ]*gqxx30
              +two*(particleJ.quadrupole[QXY]*gqxx26+particleJ.quadrupole[QXZ]*gqxx27+particleJ.quadrupole[QYZ]*gqxx29))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqyy25+particleJ.quadrupole[QYY]*gqyy28+particleJ.quadrupole[QZZ]*gqyy30
              +two*(particleJ.quadrupole[QXY]*gqyy26+particleJ.quadrupole[QXZ]*gqyy27+particleJ.quadrupole[QYZ]*gqyy29))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqzz25+particleJ.quadrupole[QYY]*gqzz28+particleJ.quadrupole[QZZ]*gqzz30
              +two*(particleJ.quadrupole[QXY]*gqzz26+particleJ.quadrupole[QXZ]*gqzz27+particleJ.quadrupole[QYZ]*gqzz29))
              + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxy25+particleJ.quadrupole[QYY]*gqxy28+particleJ.quadrupole[QZZ]*gqxy30
              +two*(particleJ.quadrupole[QXY]*gqxy26+particleJ.quadrupole[QXZ]*gqxy27+particleJ.quadrupole[QYZ]*gqxy29))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxz25+particleJ.quadrupole[QYY]*gqxz28+particleJ.quadrupole[QZZ]*gqxz30
              +two*(particleJ.quadrupole[QXY]*gqxz26+particleJ.quadrupole[QXZ]*gqxz27+particleJ.quadrupole[QYZ]*gqxz29))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqyz25+particleJ.quadrupole[QYY]*gqyz28+particleJ.quadrupole[QZZ]*gqyz30
              +two*(particleJ.quadrupole[QXY]*gqyz26+particleJ.quadrupole[QXZ]*gqyz27+particleJ.quadrupole[QYZ]*gqyz29)));

    dewkdr = particleI.charge*(particleJ.dipole[0]*gux21+particleJ.dipole[1]*guy21+particleJ.dipole[2]*guz21)
                      -particleJ.charge*(particleI.dipole[0]*gc22+particleI.dipole[1]*gc23+particleI.dipole[2]*gc24)
                 +particleI.charge*(particleJ.quadrupole[QXX]*gqxx21+particleJ.quadrupole[QYY]*gqyy21+particleJ.quadrupole[QZZ]*gqzz21
              +two*(particleJ.quadrupole[QXY]*gqxy21+particleJ.quadrupole[QXZ]*gqxz21+particleJ.quadrupole[QYZ]*gqyz21))
                 +particleJ.charge*(particleI.quadrupole[QXX]*gc25+particleI.quadrupole[QYY]*gc28+particleI.quadrupole[QZZ]*gc30
              +two*(particleI.quadrupole[QXY]*gc26+particleI.quadrupole[QXZ]*gc27+particleI.quadrupole[QYZ]*gc29))
               - particleI.dipole[0]*(particleJ.quadrupole[QXX]*gqxx22+particleJ.quadrupole[QYY]*gqyy22+particleJ.quadrupole[QZZ]*gqzz22
              +two*(particleJ.quadrupole[QXY]*gqxy22+particleJ.quadrupole[QXZ]*gqxz22+particleJ.quadrupole[QYZ]*gqyz22))
               - particleI.dipole[1]*(particleJ.quadrupole[QXX]*gqxx23+particleJ.quadrupole[QYY]*gqyy23+particleJ.quadrupole[QZZ]*gqzz23
              +two*(particleJ.quadrupole[QXY]*gqxy23+particleJ.quadrupole[QXZ]*gqxz23+particleJ.quadrupole[QYZ]*gqyz23))
               - particleI.dipole[2]*(particleJ.quadrupole[QXX]*gqxx24+particleJ.quadrupole[QYY]*gqyy24+particleJ.quadrupole[QZZ]*gqzz24
              +two*(particleJ.quadrupole[QXY]*gqxy24+particleJ.quadrupole[QXZ]*gqxz24+particleJ.quadrupole[QYZ]*gqyz24))
               + particleJ.dipole[0]*(particleI.quadrupole[QXX]*gux25+particleI.quadrupole[QYY]*gux28+particleI.quadrupole[QZZ]*gux30
              +two*(particleI.quadrupole[QXY]*gux26+particleI.quadrupole[QXZ]*gux27+particleI.quadrupole[QYZ]*gux29))
               + particleJ.dipole[1]*(particleI.quadrupole[QXX]*guy25+particleI.quadrupole[QYY]*guy28+particleI.quadrupole[QZZ]*guy30
              +two*(particleI.quadrupole[QXY]*guy26+particleI.quadrupole[QXZ]*guy27+particleI.quadrupole[QYZ]*guy29))
               + particleJ.dipole[2]*(particleI.quadrupole[QXX]*guz25+particleI.quadrupole[QYY]*guz28+particleI.quadrupole[QZZ]*guz30
              +two*(particleI.quadrupole[QXY]*guz26+particleI.quadrupole[QXZ]*guz27+particleI.quadrupole[QYZ]*guz29))
              + particleI.quadrupole[QXX]*(particleJ.quadrupole[QXX]*gqxx25+particleJ.quadrupole[QYY]*gqyy25+particleJ.quadrupole[QZZ]*gqzz25
              +two*(particleJ.quadrupole[QXY]*gqxy25+particleJ.quadrupole[QXZ]*gqxz25+particleJ.quadrupole[QYZ]*gqyz25))
              + particleI.quadrupole[QYY]*(particleJ.quadrupole[QXX]*gqxx28+particleJ.quadrupole[QYY]*gqyy28+particleJ.quadrupole[QZZ]*gqzz28
              +two*(particleJ.quadrupole[QXY]*gqxy28+particleJ.quadrupole[QXZ]*gqxz28+particleJ.quadrupole[QYZ]*gqyz28))
              + particleI.quadrupole[QZZ]*(particleJ.quadrupole[QXX]*gqxx30+particleJ.quadrupole[QYY]*gqyy30+particleJ.quadrupole[QZZ]*gqzz30
              +two*(particleJ.quadrupole[QXY]*gqxy30+particleJ.quadrupole[QXZ]*gqxz30+particleJ.quadrupole[QYZ]*gqyz30))
              + two*(particleI.quadrupole[QXY]*(particleJ.quadrupole[QXX]*gqxx26+particleJ.quadrupole[QYY]*gqyy26+particleJ.quadrupole[QZZ]*gqzz26
              +two*(particleJ.quadrupole[QXY]*gqxy26+particleJ.quadrupole[QXZ]*gqxz26+particleJ.quadrupole[QYZ]*gqyz26))
              + particleI.quadrupole[QXZ]*(particleJ.quadrupole[QXX]*gqxx27+particleJ.quadrupole[QYY]*gqyy27+particleJ.quadrupole[QZZ]*gqzz27
              +two*(particleJ.quadrupole[QXY]*gqxy27+particleJ.quadrupole[QXZ]*gqxz27+particleJ.quadrupole[QYZ]*gqyz27))
              + particleI.quadrupole[QYZ]*(particleJ.quadrupole[QXX]*gqxx29+particleJ.quadrupole[QYY]*gqyy29+particleJ.quadrupole[QZZ]*gqzz29
              +two*(particleJ.quadrupole[QXY]*gqxy29+particleJ.quadrupole[QXZ]*gqxz29+particleJ.quadrupole[QYZ]*gqyz29)));

    dsumdr = desymdr + oneHalf*(dewidr + dewkdr);
RealOpenMM dsumdrXXX = dsumdr;
    drbi = _bornRadii[jIndex]*dsumdr;
    drbk = _bornRadii[iIndex]*dsumdr;

    // torque on permanent dipoles due to permanent reaction field

    RealOpenMM trq1   = 0.0;
    RealOpenMM trq2   = 0.0;
    RealOpenMM trq3   = 0.0;

    RealOpenMM trq_k1 = 0.0;
    RealOpenMM trq_k2 = 0.0;
    RealOpenMM trq_k3 = 0.0;

    if ( xr != 0.0 || yr != 0.0 || zr != 0.0 )
    {

        RealOpenMM fid1 = particleJ.dipole[0]*gux2 + particleJ.dipole[1]*gux3 + particleJ.dipole[2]*gux4
                + oneHalf*(particleJ.charge*gux1+particleJ.quadrupole[QXX]*gux5+particleJ.quadrupole[QYY]*gux8+particleJ.quadrupole[QZZ]*gux10
                      +two*(particleJ.quadrupole[QXY]*gux6+particleJ.quadrupole[QXZ]*gux7+particleJ.quadrupole[QYZ]*gux9)
                      +particleJ.charge*gc2+particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
                      +two*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2));

        RealOpenMM fid2 = particleJ.dipole[0]*guy2 + particleJ.dipole[1]*guy3 + particleJ.dipole[2]*guy4
                + oneHalf*(particleJ.charge*guy1+particleJ.quadrupole[QXX]*guy5+particleJ.quadrupole[QYY]*guy8+particleJ.quadrupole[QZZ]*guy10
                      +two*(particleJ.quadrupole[QXY]*guy6+particleJ.quadrupole[QXZ]*guy7+particleJ.quadrupole[QYZ]*guy9)
                      +particleJ.charge*gc3+particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
                      +two*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3));

        RealOpenMM fid3 = particleJ.dipole[0]*guz2 + particleJ.dipole[1]*guz3 + particleJ.dipole[2]*guz4
                + oneHalf*(particleJ.charge*guz1+particleJ.quadrupole[QXX]*guz5+particleJ.quadrupole[QYY]*guz8+particleJ.quadrupole[QZZ]*guz10
                      +two*(particleJ.quadrupole[QXY]*guz6+particleJ.quadrupole[QXZ]*guz7+particleJ.quadrupole[QYZ]*guz9)
                      +particleJ.charge*gc4+particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
                      +two*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4));

        RealOpenMM fkd1 = particleI.dipole[0]*gux2 + particleI.dipole[1]*gux3 + particleI.dipole[2]*gux4
                - oneHalf*(particleI.charge*gux1+particleI.quadrupole[QXX]*gux5+particleI.quadrupole[QYY]*gux8+particleI.quadrupole[QZZ]*gux10
                      +two*(particleI.quadrupole[QXY]*gux6+particleI.quadrupole[QXZ]*gux7+particleI.quadrupole[QYZ]*gux9)
                      +particleI.charge*gc2+particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
                      +two*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2));

        RealOpenMM fkd2 = particleI.dipole[0]*guy2 + particleI.dipole[1]*guy3 + particleI.dipole[2]*guy4
                - oneHalf*(particleI.charge*guy1+particleI.quadrupole[QXX]*guy5+particleI.quadrupole[QYY]*guy8+particleI.quadrupole[QZZ]*guy10
                      +two*(particleI.quadrupole[QXY]*guy6+particleI.quadrupole[QXZ]*guy7+particleI.quadrupole[QYZ]*guy9)
                      +particleI.charge*gc3+particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
                      +two*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3));

        RealOpenMM fkd3 = particleI.dipole[0]*guz2 + particleI.dipole[1]*guz3 + particleI.dipole[2]*guz4
                - oneHalf*(particleI.charge*guz1+particleI.quadrupole[QXX]*guz5+particleI.quadrupole[QYY]*guz8+particleI.quadrupole[QZZ]*guz10
                      +two*(particleI.quadrupole[QXY]*guz6+particleI.quadrupole[QXZ]*guz7+particleI.quadrupole[QYZ]*guz9)
                      +particleI.charge*gc4+particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
                      +two*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4));

        trq1    = particleI.dipole[1]*fid3 - particleI.dipole[2]*fid2;
        trq2    = particleI.dipole[2]*fid1 - particleI.dipole[0]*fid3;
        trq3    = particleI.dipole[0]*fid2 - particleI.dipole[1]*fid1;

        trq_k1  = particleJ.dipole[1]*fkd3 - particleJ.dipole[2]*fkd2;
        trq_k2  = particleJ.dipole[2]*fkd1 - particleJ.dipole[0]*fkd3;
        trq_k3  = particleJ.dipole[0]*fkd2 - particleJ.dipole[1]*fkd1;

        // torque on quadrupoles due to permanent reaction field gradient

        RealOpenMM fidg11 =
                - oneHalf*(particleJ.charge*gqxx1+particleJ.dipole[0]*gqxx2+particleJ.dipole[1]*gqxx3+particleJ.dipole[2]*gqxx4
                      +particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqxx8+particleJ.quadrupole[QZZ]*gqxx10
                      +two*(particleJ.quadrupole[QXY]*gqxx6+particleJ.quadrupole[QXZ]*gqxx7+particleJ.quadrupole[QYZ]*gqxx9)
                      +particleJ.charge*gc5+particleJ.dipole[0]*gux5+particleJ.dipole[1]*guy5+particleJ.dipole[2]*guz5
                      +particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
                      +two*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5));

        RealOpenMM fidg12 =
                - oneHalf*(particleJ.charge*gqxy1+particleJ.dipole[0]*gqxy2+particleJ.dipole[1]*gqxy3+particleJ.dipole[2]*gqxy4
                      +particleJ.quadrupole[QXX]*gqxy5+particleJ.quadrupole[QYY]*gqxy8+particleJ.quadrupole[QZZ]*gqxy10
                      +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxy7+particleJ.quadrupole[QYZ]*gqxy9)
                      +particleJ.charge*gc6+particleJ.dipole[0]*gux6+particleJ.dipole[1]*guy6+particleJ.dipole[2]*guz6
                      +particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
                      +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6));

        RealOpenMM fidg13 =
                - oneHalf*(particleJ.charge*gqxz1+particleJ.dipole[0]*gqxz2+particleJ.dipole[1]*gqxz3+particleJ.dipole[2]*gqxz4
                      +particleJ.quadrupole[QXX]*gqxz5+particleJ.quadrupole[QYY]*gqxz8+particleJ.quadrupole[QZZ]*gqxz10
                      +two*(particleJ.quadrupole[QXY]*gqxz6+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqxz9)
                      +particleJ.charge*gc7+particleJ.dipole[0]*gux7+particleJ.dipole[1]*guy7+particleJ.dipole[2]*guz7
                      +particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
                      +two*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7));

        RealOpenMM fidg22 =
                - oneHalf*(particleJ.charge*gqyy1+particleJ.dipole[0]*gqyy2+particleJ.dipole[1]*gqyy3+particleJ.dipole[2]*gqyy4
                      +particleJ.quadrupole[QXX]*gqyy5+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqyy10
                      +two*(particleJ.quadrupole[QXY]*gqyy6+particleJ.quadrupole[QXZ]*gqyy7+particleJ.quadrupole[QYZ]*gqyy9)
                      +particleJ.charge*gc8+particleJ.dipole[0]*gux8+particleJ.dipole[1]*guy8+particleJ.dipole[2]*guz8
                      +particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
                      +two*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8));

        RealOpenMM fidg23 =
                - oneHalf*(particleJ.charge*gqyz1+particleJ.dipole[0]*gqyz2+particleJ.dipole[1]*gqyz3+particleJ.dipole[2]*gqyz4
                      +particleJ.quadrupole[QXX]*gqyz5+particleJ.quadrupole[QYY]*gqyz8+particleJ.quadrupole[QZZ]*gqyz10
                      +two*(particleJ.quadrupole[QXY]*gqyz6+particleJ.quadrupole[QXZ]*gqyz7+particleJ.quadrupole[QYZ]*gqyz9)
                      +particleJ.charge*gc9+particleJ.dipole[0]*gux9+particleJ.dipole[1]*guy9+particleJ.dipole[2]*guz9
                      +particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
                      +two*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9));

        RealOpenMM fidg33 =
                - oneHalf*(particleJ.charge*gqzz1+particleJ.dipole[0]*gqzz2+particleJ.dipole[1]*gqzz3+particleJ.dipole[2]*gqzz4
                      +particleJ.quadrupole[QXX]*gqzz5+particleJ.quadrupole[QYY]*gqzz8+particleJ.quadrupole[QZZ]*gqzz10
                      +two*(particleJ.quadrupole[QXY]*gqzz6+particleJ.quadrupole[QXZ]*gqzz7+particleJ.quadrupole[QYZ]*gqzz9)
                      +particleJ.charge*gc10+particleJ.dipole[0]*gux10+particleJ.dipole[1]*guy10+particleJ.dipole[2]*guz10
                      +particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
                   +two*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10));

        RealOpenMM fidg21 = fidg12;
        RealOpenMM fidg31 = fidg13;
        RealOpenMM fidg32 = fidg23;

        RealOpenMM fkdg11 =
                - oneHalf*(particleI.charge*gqxx1-particleI.dipole[0]*gqxx2-particleI.dipole[1]*gqxx3-particleI.dipole[2] *gqxx4
                      +particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqxx8+particleI.quadrupole[QZZ]*gqxx10
                      +two*(particleI.quadrupole[QXY]*gqxx6+particleI.quadrupole[QXZ]*gqxx7+particleI.quadrupole[QYZ]*gqxx9)
                      +particleI.charge*gc5-particleI.dipole[0]*gux5-particleI.dipole[1]*guy5-particleI.dipole[2]*guz5
                      +particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqyy5+particleI.quadrupole[QZZ]*gqzz5
                      +two*(particleI.quadrupole[QXY]*gqxy5+particleI.quadrupole[QXZ]*gqxz5+particleI.quadrupole[QYZ]*gqyz5));

        RealOpenMM fkdg12 =
                - oneHalf*(particleI.charge*gqxy1-particleI.dipole[0]*gqxy2-particleI.dipole[1]*gqxy3-particleI.dipole[2]*gqxy4
                      +particleI.quadrupole[QXX]*gqxy5+particleI.quadrupole[QYY]*gqxy8+particleI.quadrupole[QZZ]*gqxy10
                      +two*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxy7+particleI.quadrupole[QYZ]*gqxy9)
                      +particleI.charge*gc6-particleI.dipole[0]*gux6-particleI.dipole[1]*guy6-particleI.dipole[2]*guz6
                      +particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
                      +two*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6));

        RealOpenMM fkdg13 =
                - oneHalf*(particleI.charge*gqxz1-particleI.dipole[0]*gqxz2-particleI.dipole[1]*gqxz3-particleI.dipole[2]*gqxz4
                      +particleI.quadrupole[QXX]*gqxz5+particleI.quadrupole[QYY]*gqxz8+particleI.quadrupole[QZZ]*gqxz10
                      +two*(particleI.quadrupole[QXY]*gqxz6+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqxz9)
                      +particleI.charge*gc7-particleI.dipole[0]*gux7-particleI.dipole[1]*guy7-particleI.dipole[2]*guz7
                      +particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
                      +two*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7));

        RealOpenMM fkdg22 =
                - oneHalf*(particleI.charge*gqyy1-particleI.dipole[0]*gqyy2-particleI.dipole[1]*gqyy3-particleI.dipole[2]*gqyy4
                      +particleI.quadrupole[QXX]*gqyy5+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqyy10
                      +two*(particleI.quadrupole[QXY]*gqyy6+particleI.quadrupole[QXZ]*gqyy7+particleI.quadrupole[QYZ]*gqyy9)
                      +particleI.charge*gc8-particleI.dipole[0]*gux8-particleI.dipole[1]*guy8-particleI.dipole[2]*guz8
                      +particleI.quadrupole[QXX]*gqxx8+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqzz8
                      +two*(particleI.quadrupole[QXY]*gqxy8+particleI.quadrupole[QXZ]*gqxz8+particleI.quadrupole[QYZ]*gqyz8));

        RealOpenMM fkdg23 =
                - oneHalf*(particleI.charge*gqyz1-particleI.dipole[0]*gqyz2-particleI.dipole[1]*gqyz3-particleI.dipole[2]*gqyz4
                      +particleI.quadrupole[QXX]*gqyz5+particleI.quadrupole[QYY]*gqyz8+particleI.quadrupole[QZZ]*gqyz10
                      +two*(particleI.quadrupole[QXY]*gqyz6+particleI.quadrupole[QXZ]*gqyz7+particleI.quadrupole[QYZ]*gqyz9)
                      +particleI.charge*gc9-particleI.dipole[0]*gux9-particleI.dipole[1]*guy9-particleI.dipole[2]*guz9
                      +particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
                      +two*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9));
        RealOpenMM fkdg33 =
                - oneHalf*(particleI.charge*gqzz1-particleI.dipole[0]*gqzz2-particleI.dipole[1]*gqzz3-particleI.dipole[2]*gqzz4
                      +particleI.quadrupole[QXX]*gqzz5+particleI.quadrupole[QYY]*gqzz8+particleI.quadrupole[QZZ]*gqzz10
                      +two*(particleI.quadrupole[QXY]*gqzz6+particleI.quadrupole[QXZ]*gqzz7+particleI.quadrupole[QYZ]*gqzz9)
                      +particleI.charge*gc10-particleI.dipole[0]*gux10-particleI.dipole[1]*guy10-particleI.dipole[2]*guz10
                      +particleI.quadrupole[QXX]*gqxx10+particleI.quadrupole[QYY]*gqyy10+particleI.quadrupole[QZZ]*gqzz10
                    +two*(particleI.quadrupole[QXY]*gqxy10+particleI.quadrupole[QXZ]*gqxz10+particleI.quadrupole[QYZ]*gqyz10));

        RealOpenMM fkdg21 = fkdg12;
        RealOpenMM fkdg31 = fkdg13;
        RealOpenMM fkdg32 = fkdg23;

        trq1   += two* (particleI.quadrupole[QXY]*fidg13+particleI.quadrupole[QYY]*fidg23+particleI.quadrupole[QYZ]*fidg33
                           -particleI.quadrupole[QXZ]*fidg12-particleI.quadrupole[QYZ]*fidg22-particleI.quadrupole[QZZ]*fidg32);

        trq2   += two*(particleI.quadrupole[QXZ]*fidg11+particleI.quadrupole[QYZ]*fidg21+particleI.quadrupole[QZZ]*fidg31
                         -particleI.quadrupole[QXX]*fidg13-particleI.quadrupole[QXY]*fidg23-particleI.quadrupole[QXZ]*fidg33);

        trq3   += two*(particleI.quadrupole[QXX]*fidg12+particleI.quadrupole[QXY]*fidg22+particleI.quadrupole[QXZ]*fidg32
                         -particleI.quadrupole[QXY]*fidg11-particleI.quadrupole[QYY]*fidg21-particleI.quadrupole[QYZ]*fidg31);

        trq_k1 += two*
                          (particleJ.quadrupole[QXY]*fkdg13+particleJ.quadrupole[QYY]*fkdg23+particleJ.quadrupole[QYZ]*fkdg33
                          -particleJ.quadrupole[QXZ]*fkdg12-particleJ.quadrupole[QYZ]*fkdg22-particleJ.quadrupole[QZZ]*fkdg32);

        trq_k2 += two*
                          (particleJ.quadrupole[QXZ]*fkdg11+particleJ.quadrupole[QYZ]*fkdg21+particleJ.quadrupole[QZZ]*fkdg31
                          -particleJ.quadrupole[QXX]*fkdg13-particleJ.quadrupole[QXY]*fkdg23-particleJ.quadrupole[QXZ]*fkdg33);

        trq_k3 += two*
                          (particleJ.quadrupole[QXX]*fkdg12+particleJ.quadrupole[QXY]*fkdg22+particleJ.quadrupole[QXZ]*fkdg32
                          -particleJ.quadrupole[QXY]*fkdg11-particleJ.quadrupole[QYY]*fkdg21-particleJ.quadrupole[QYZ]*fkdg31);
    }

    // electrostatic solvation energy of the permanent multipoles in
    // the GK reaction potential of the induced dipoles

    esymi =              -particleI.dipole[0]*(_inducedDipoleS[jIndex][0]*gux2+_inducedDipoleS[jIndex][1]*guy2+_inducedDipoleS[jIndex][2]*guz2)
                        - particleI.dipole[1]*(_inducedDipoleS[jIndex][0]*gux3+_inducedDipoleS[jIndex][1]*guy3+_inducedDipoleS[jIndex][2]*guz3)
                        - particleI.dipole[2]*(_inducedDipoleS[jIndex][0]*gux4+_inducedDipoleS[jIndex][1]*guy4+_inducedDipoleS[jIndex][2]*guz4)
                        - particleJ.dipole[0]*(_inducedDipoleS[iIndex][0]*gux2+_inducedDipoleS[iIndex][1]*guy2+_inducedDipoleS[iIndex][2]*guz2)
                        - particleJ.dipole[1]*(_inducedDipoleS[iIndex][0]*gux3+_inducedDipoleS[iIndex][1]*guy3+_inducedDipoleS[iIndex][2]*guz3)
                        - particleJ.dipole[2]*(_inducedDipoleS[iIndex][0]*gux4+_inducedDipoleS[iIndex][1]*guy4+_inducedDipoleS[iIndex][2]*guz4);

    ewii = particleI.charge*(_inducedDipoleS[jIndex][0]*gc2+_inducedDipoleS[jIndex][1]*gc3+_inducedDipoleS[jIndex][2]*gc4)
                      - particleJ.charge*(_inducedDipoleS[iIndex][0]*gux1+_inducedDipoleS[iIndex][1]*guy1+_inducedDipoleS[iIndex][2]*guz1)
                      - _inducedDipoleS[iIndex][0]*(particleJ.quadrupole[QXX]*gux5+particleJ.quadrupole[QYY]*gux8+particleJ.quadrupole[QZZ]*gux10
                     +two*(particleJ.quadrupole[QXY]*gux6+particleJ.quadrupole[QXZ]*gux7+particleJ.quadrupole[QYZ]*gux9))
                      - _inducedDipoleS[iIndex][1]*(particleJ.quadrupole[QXX]*guy5+particleJ.quadrupole[QYY]*guy8+particleJ.quadrupole[QZZ]*guy10
                     +two*(particleJ.quadrupole[QXY]*guy6+particleJ.quadrupole[QXZ]*guy7+particleJ.quadrupole[QYZ]*guy9))
                      - _inducedDipoleS[iIndex][2]*(particleJ.quadrupole[QXX]*guz5+particleJ.quadrupole[QYY]*guz8+particleJ.quadrupole[QZZ]*guz10
                     +two*(particleJ.quadrupole[QXY]*guz6+particleJ.quadrupole[QXZ]*guz7+particleJ.quadrupole[QYZ]*guz9))
                      + _inducedDipoleS[jIndex][0]*(particleI.quadrupole[QXX]*gqxx2+particleI.quadrupole[QYY]*gqyy2+particleI.quadrupole[QZZ]*gqzz2
                     +two*(particleI.quadrupole[QXY]*gqxy2+particleI.quadrupole[QXZ]*gqxz2+particleI.quadrupole[QYZ]*gqyz2))
                      + _inducedDipoleS[jIndex][1]*(particleI.quadrupole[QXX]*gqxx3+particleI.quadrupole[QYY]*gqyy3+particleI.quadrupole[QZZ]*gqzz3
                     +two*(particleI.quadrupole[QXY]*gqxy3+particleI.quadrupole[QXZ]*gqxz3+particleI.quadrupole[QYZ]*gqyz3))
                      + _inducedDipoleS[jIndex][2]*(particleI.quadrupole[QXX]*gqxx4+particleI.quadrupole[QYY]*gqyy4+particleI.quadrupole[QZZ]*gqzz4
                     +two*(particleI.quadrupole[QXY]*gqxy4+particleI.quadrupole[QXZ]*gqxz4+particleI.quadrupole[QYZ]*gqyz4));

    ewki = particleI.charge*(_inducedDipoleS[jIndex][0]*gux1+_inducedDipoleS[jIndex][1]*guy1+_inducedDipoleS[jIndex][2]*guz1)
                      - particleJ.charge*(_inducedDipoleS[iIndex][0]*gc2+_inducedDipoleS[iIndex][1]*gc3+_inducedDipoleS[iIndex][2]*gc4)
                      - _inducedDipoleS[iIndex][0]*(particleJ.quadrupole[QXX]*gqxx2+particleJ.quadrupole[QYY]*gqyy2+particleJ.quadrupole[QZZ]*gqzz2
                     +two*(particleJ.quadrupole[QXY]*gqxy2+particleJ.quadrupole[QXZ]*gqxz2+particleJ.quadrupole[QYZ]*gqyz2))
                      - _inducedDipoleS[iIndex][1]*(particleJ.quadrupole[QXX]*gqxx3+particleJ.quadrupole[QYY]*gqyy3+particleJ.quadrupole[QZZ]*gqzz3
                     +two*(particleJ.quadrupole[QXY]*gqxy3+particleJ.quadrupole[QXZ]*gqxz3+particleJ.quadrupole[QYZ]*gqyz3))
                      - _inducedDipoleS[iIndex][2]*(particleJ.quadrupole[QXX]*gqxx4+particleJ.quadrupole[QYY]*gqyy4+particleJ.quadrupole[QZZ]*gqzz4
                     +two*(particleJ.quadrupole[QXY]*gqxy4+particleJ.quadrupole[QXZ]*gqxz4+particleJ.quadrupole[QYZ]*gqyz4))
                      + _inducedDipoleS[jIndex][0]*(particleI.quadrupole[QXX]*gux5+particleI.quadrupole[QYY]*gux8+particleI.quadrupole[QZZ]*gux10
                     +two*(particleI.quadrupole[QXY]*gux6+particleI.quadrupole[QXZ]*gux7+particleI.quadrupole[QYZ]*gux9))
                      + _inducedDipoleS[jIndex][1]*(particleI.quadrupole[QXX]*guy5+particleI.quadrupole[QYY]*guy8+particleI.quadrupole[QZZ]*guy10
                     +two*(particleI.quadrupole[QXY]*guy6+particleI.quadrupole[QXZ]*guy7+particleI.quadrupole[QYZ]*guy9))
                      + _inducedDipoleS[jIndex][2]*(particleI.quadrupole[QXX]*guz5+particleI.quadrupole[QYY]*guz8+particleI.quadrupole[QZZ]*guz10
                     +two*(particleI.quadrupole[QXY]*guz6+particleI.quadrupole[QXZ]*guz7+particleI.quadrupole[QYZ]*guz9));

    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles

    RealOpenMM dpsymdx =       - particleI.dipole[0]*(sxk*gux5+syk*guy5+szk*guz5)
                          - particleI.dipole[1]*(sxk*gux6+syk*guy6+szk*guz6)
                          - particleI.dipole[2]*(sxk*gux7+syk*guy7+szk*guz7)

                          - particleJ.dipole[0]*(sxi*gux5+syi*guy5+szi*guz5)
                          - particleJ.dipole[1]*(sxi*gux6+syi*guy6+szi*guz6)
                          - particleJ.dipole[2]*(sxi*gux7+syi*guy7+szi*guz7);

    dpwidx = particleI.charge*(sxk*gc5+syk*gc6+szk*gc7)
                        - particleJ.charge*(sxi*gux2+syi*guy2+szi*guz2)
                      - sxi*(particleJ.quadrupole[QXX]*gux11+particleJ.quadrupole[QYY]*gux14+particleJ.quadrupole[QZZ]*gux16
                     +two*(particleJ.quadrupole[QXY]*gux12+particleJ.quadrupole[QXZ]*gux13+particleJ.quadrupole[QYZ]*gux15))
                      - syi*(particleJ.quadrupole[QXX]*guy11+particleJ.quadrupole[QYY]*guy14+particleJ.quadrupole[QZZ]*guy16
                     +two*(particleJ.quadrupole[QXY]*guy12+particleJ.quadrupole[QXZ]*guy13+particleJ.quadrupole[QYZ]*guy15))
                      - szi*(particleJ.quadrupole[QXX]*guz11+particleJ.quadrupole[QYY]*guz14+particleJ.quadrupole[QZZ]*guz16
                     +two*(particleJ.quadrupole[QXY]*guz12+particleJ.quadrupole[QXZ]*guz13+particleJ.quadrupole[QYZ]*guz15))
                      + sxk*(particleI.quadrupole[QXX]*gqxx5+particleI.quadrupole[QYY]*gqyy5+particleI.quadrupole[QZZ]*gqzz5
                     +two*(particleI.quadrupole[QXY]*gqxy5+particleI.quadrupole[QXZ]*gqxz5+particleI.quadrupole[QYZ]*gqyz5))
                      + syk*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
                     +two*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
                      + szk*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
                     +two*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7));

    dpwkdx = particleI.charge*(sxk*gux2+syk*guy2+szk*guz2)
                        - particleJ.charge*(sxi*gc5+syi*gc6+szi*gc7)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx5+particleJ.quadrupole[QYY]*gqyy5+particleJ.quadrupole[QZZ]*gqzz5
                     +two*(particleJ.quadrupole[QXY]*gqxy5+particleJ.quadrupole[QXZ]*gqxz5+particleJ.quadrupole[QYZ]*gqyz5))
                      - syi*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
                     +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
                      - szi*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
                     +two*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
                      + sxk*(particleI.quadrupole[QXX]*gux11+particleI.quadrupole[QYY]*gux14+particleI.quadrupole[QZZ]*gux16
                     +two*(particleI.quadrupole[QXY]*gux12+particleI.quadrupole[QXZ]*gux13+particleI.quadrupole[QYZ]*gux15))
                      + syk*(particleI.quadrupole[QXX]*guy11+particleI.quadrupole[QYY]*guy14+particleI.quadrupole[QZZ]*guy16
                     +two*(particleI.quadrupole[QXY]*guy12+particleI.quadrupole[QXZ]*guy13+particleI.quadrupole[QYZ]*guy15))
                      + szk*(particleI.quadrupole[QXX]*guz11+particleI.quadrupole[QYY]*guz14+particleI.quadrupole[QZZ]*guz16
                     +two*(particleI.quadrupole[QXY]*guz12+particleI.quadrupole[QXZ]*guz13+particleI.quadrupole[QYZ]*guz15));

    dpdx = oneHalf*(dpsymdx + oneHalf*(dpwidx + dpwkdx));

    dpsymdy = -particleI.dipole[0]*(sxk*gux6+syk*guy6+szk*guz6)
                          - particleI.dipole[1]*(sxk*gux8+syk*guy8+szk*guz8)
                          - particleI.dipole[2]*(sxk*gux9+syk*guy9+szk*guz9)
                          - particleJ.dipole[0]*(sxi*gux6+syi*guy6+szi*guz6)
                          - particleJ.dipole[1]*(sxi*gux8+syi*guy8+szi*guz8)
                          - particleJ.dipole[2]*(sxi*gux9+syi*guy9+szi*guz9);

    dpwidy = particleI.charge*(sxk*gc6+syk*gc8+szk*gc9)
                        - particleJ.charge*(sxi*gux3+syi*guy3+szi*guz3)
                         - sxi*(particleJ.quadrupole[QXX]*gux12+particleJ.quadrupole[QYY]*gux17+particleJ.quadrupole[QZZ]*gux19
                        +two*(particleJ.quadrupole[QXY]*gux14+particleJ.quadrupole[QXZ]*gux15+particleJ.quadrupole[QYZ]*gux18))
                         - syi*(particleJ.quadrupole[QXX]*guy12+particleJ.quadrupole[QYY]*guy17+particleJ.quadrupole[QZZ]*guy19
                        +two*(particleJ.quadrupole[QXY]*guy14+particleJ.quadrupole[QXZ]*guy15+particleJ.quadrupole[QYZ]*guy18))
                         - szi*(particleJ.quadrupole[QXX]*guz12+particleJ.quadrupole[QYY]*guz17+particleJ.quadrupole[QZZ]*guz19
                        +two*(particleJ.quadrupole[QXY]*guz14+particleJ.quadrupole[QXZ]*guz15+particleJ.quadrupole[QYZ]*guz18))
                         + sxk*(particleI.quadrupole[QXX]*gqxx6+particleI.quadrupole[QYY]*gqyy6+particleI.quadrupole[QZZ]*gqzz6
                        +two*(particleI.quadrupole[QXY]*gqxy6+particleI.quadrupole[QXZ]*gqxz6+particleI.quadrupole[QYZ]*gqyz6))
                         + syk*(particleI.quadrupole[QXX]*gqxx8+particleI.quadrupole[QYY]*gqyy8+particleI.quadrupole[QZZ]*gqzz8
                        +two*(particleI.quadrupole[QXY]*gqxy8+particleI.quadrupole[QXZ]*gqxz8+particleI.quadrupole[QYZ]*gqyz8))
                         + szk*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
                        +two*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9));

    dpwkdy = particleI.charge*(sxk*gux3+syk*guy3+szk*guz3)
                        - particleJ.charge*(sxi*gc6+syi*gc8+szi*gc9)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx6+particleJ.quadrupole[QYY]*gqyy6+particleJ.quadrupole[QZZ]*gqzz6
                     +two*(particleJ.quadrupole[QXY]*gqxy6+particleJ.quadrupole[QXZ]*gqxz6+particleJ.quadrupole[QYZ]*gqyz6))
                      - syi*(particleJ.quadrupole[QXX]*gqxx8+particleJ.quadrupole[QYY]*gqyy8+particleJ.quadrupole[QZZ]*gqzz8
                     +two*(particleJ.quadrupole[QXY]*gqxy8+particleJ.quadrupole[QXZ]*gqxz8+particleJ.quadrupole[QYZ]*gqyz8))
                      - szi*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
                     +two*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
                      + sxk*(particleI.quadrupole[QXX]*gux12+particleI.quadrupole[QYY]*gux17+particleI.quadrupole[QZZ]*gux19
                     +two*(particleI.quadrupole[QXY]*gux14+particleI.quadrupole[QXZ]*gux15+particleI.quadrupole[QYZ]*gux18))
                      + syk*(particleI.quadrupole[QXX]*guy12+particleI.quadrupole[QYY]*guy17+particleI.quadrupole[QZZ]*guy19
                     +two*(particleI.quadrupole[QXY]*guy14+particleI.quadrupole[QXZ]*guy15+particleI.quadrupole[QYZ]*guy18))
                      + szk*(particleI.quadrupole[QXX]*guz12+particleI.quadrupole[QYY]*guz17+particleI.quadrupole[QZZ]*guz19
                     +two*(particleI.quadrupole[QXY]*guz14+particleI.quadrupole[QXZ]*guz15+particleI.quadrupole[QYZ]*guz18));

    dpdy    = oneHalf*(dpsymdy + oneHalf*(dpwidy + dpwkdy));

    dpsymdz = -particleI.dipole[0]*(sxk*gux7+syk*guy7+szk*guz7)
                          - particleI.dipole[1]*(sxk*gux9+syk*guy9+szk*guz9)
                          - particleI.dipole[2]*(sxk*gux10+syk*guy10+szk*guz10)
                          - particleJ.dipole[0]*(sxi*gux7+syi*guy7+szi*guz7)
                          - particleJ.dipole[1]*(sxi*gux9+syi*guy9+szi*guz9)
                          - particleJ.dipole[2]*(sxi*gux10+syi*guy10+szi*guz10);

    dpwidz = particleI.charge*(sxk*gc7+syk*gc9+szk*gc10)
                        - particleJ.charge*(sxi*gux4+syi*guy4+szi*guz4)
                      - sxi*(particleJ.quadrupole[QXX]*gux13+particleJ.quadrupole[QYY]*gux18+particleJ.quadrupole[QZZ]*gux20
                     +two*(particleJ.quadrupole[QXY]*gux15+particleJ.quadrupole[QXZ]*gux16+particleJ.quadrupole[QYZ]*gux19))
                      - syi*(particleJ.quadrupole[QXX]*guy13+particleJ.quadrupole[QYY]*guy18+particleJ.quadrupole[QZZ]*guy20
                     +two*(particleJ.quadrupole[QXY]*guy15+particleJ.quadrupole[QXZ]*guy16+particleJ.quadrupole[QYZ]*guy19))
                      - szi*(particleJ.quadrupole[QXX]*guz13+particleJ.quadrupole[QYY]*guz18+particleJ.quadrupole[QZZ]*guz20
                     +two*(particleJ.quadrupole[QXY]*guz15+particleJ.quadrupole[QXZ]*guz16+particleJ.quadrupole[QYZ]*guz19))
                      + sxk*(particleI.quadrupole[QXX]*gqxx7+particleI.quadrupole[QYY]*gqyy7+particleI.quadrupole[QZZ]*gqzz7
                     +two*(particleI.quadrupole[QXY]*gqxy7+particleI.quadrupole[QXZ]*gqxz7+particleI.quadrupole[QYZ]*gqyz7))
                      + syk*(particleI.quadrupole[QXX]*gqxx9+particleI.quadrupole[QYY]*gqyy9+particleI.quadrupole[QZZ]*gqzz9
                     +two*(particleI.quadrupole[QXY]*gqxy9+particleI.quadrupole[QXZ]*gqxz9+particleI.quadrupole[QYZ]*gqyz9))
                      + szk*(particleI.quadrupole[QXX]*gqxx10+particleI.quadrupole[QYY]*gqyy10+particleI.quadrupole[QZZ]*gqzz10
                     +two*(particleI.quadrupole[QXY]*gqxy10+particleI.quadrupole[QXZ]*gqxz10+particleI.quadrupole[QYZ]*gqyz10));

    dpwkdz = particleI.charge*(sxk*gux4+syk*guy4+szk*guz4)
                        - particleJ.charge*(sxi*gc7+syi*gc9+szi*gc10)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx7+particleJ.quadrupole[QYY]*gqyy7+particleJ.quadrupole[QZZ]*gqzz7
                     +two*(particleJ.quadrupole[QXY]*gqxy7+particleJ.quadrupole[QXZ]*gqxz7+particleJ.quadrupole[QYZ]*gqyz7))
                      - syi*(particleJ.quadrupole[QXX]*gqxx9+particleJ.quadrupole[QYY]*gqyy9+particleJ.quadrupole[QZZ]*gqzz9
                     +two*(particleJ.quadrupole[QXY]*gqxy9+particleJ.quadrupole[QXZ]*gqxz9+particleJ.quadrupole[QYZ]*gqyz9))
                      - szi*(particleJ.quadrupole[QXX]*gqxx10+particleJ.quadrupole[QYY]*gqyy10+particleJ.quadrupole[QZZ]*gqzz10
                     +two*(particleJ.quadrupole[QXY]*gqxy10+particleJ.quadrupole[QXZ]*gqxz10+particleJ.quadrupole[QYZ]*gqyz10))
                      + sxk*(particleI.quadrupole[QXX]*gux13+particleI.quadrupole[QYY]*gux18+particleI.quadrupole[QZZ]*gux20
                     +two*(particleI.quadrupole[QXY]*gux15+particleI.quadrupole[QXZ]*gux16+particleI.quadrupole[QYZ]*gux19))
                      + syk*(particleI.quadrupole[QXX]*guy13+particleI.quadrupole[QYY]*guy18+particleI.quadrupole[QZZ]*guy20
                     +two*(particleI.quadrupole[QXY]*guy15+particleI.quadrupole[QXZ]*guy16+particleI.quadrupole[QYZ]*guy19))
                      + szk*(particleI.quadrupole[QXX]*guz13+particleI.quadrupole[QYY]*guz18+particleI.quadrupole[QZZ]*guz20
                     +two*(particleI.quadrupole[QXY]*guz15+particleI.quadrupole[QXZ]*guz16+particleI.quadrupole[QYZ]*guz19));

    dpdz = oneHalf*(dpsymdz + oneHalf*(dpwidz + dpwkdz));

    // effective radii chain rule terms for the
    // electrostatic solvation free energy gradient of the permanent
    // multipoles in the reaction potential of the induced dipoles

    dsymdr = -particleI.dipole[0]*(sxk*gux22+syk*guy22+szk*guz22)
                          - particleI.dipole[1]*(sxk*gux23+syk*guy23+szk*guz23)
                          - particleI.dipole[2]*(sxk*gux24+syk*guy24+szk*guz24)
                          - particleJ.dipole[0]*(sxi*gux22+syi*guy22+szi*guz22)
                          - particleJ.dipole[1]*(sxi*gux23+syi*guy23+szi*guz23)
                          - particleJ.dipole[2]*(sxi*gux24+syi*guy24+szi*guz24);

    dwipdr = particleI.charge*(sxk*gc22+syk*gc23+szk*gc24)
                         - particleJ.charge*(sxi*gux21+syi*guy21+szi*guz21)
                      - sxi*(particleJ.quadrupole[QXX]*gux25+particleJ.quadrupole[QYY]*gux28+particleJ.quadrupole[QZZ]*gux30
                     +two*(particleJ.quadrupole[QXY]*gux26+particleJ.quadrupole[QXZ]*gux27+particleJ.quadrupole[QYZ]*gux29))
                      - syi*(particleJ.quadrupole[QXX]*guy25+particleJ.quadrupole[QYY]*guy28+particleJ.quadrupole[QZZ]*guy30
                     +two*(particleJ.quadrupole[QXY]*guy26+particleJ.quadrupole[QXZ]*guy27+particleJ.quadrupole[QYZ]*guy29))
                      - szi*(particleJ.quadrupole[QXX]*guz25+particleJ.quadrupole[QYY]*guz28+particleJ.quadrupole[QZZ]*guz30
                     +two*(particleJ.quadrupole[QXY]*guz26+particleJ.quadrupole[QXZ]*guz27+particleJ.quadrupole[QYZ]*guz29))
                      + sxk*(particleI.quadrupole[QXX]*gqxx22+particleI.quadrupole[QYY]*gqyy22+particleI.quadrupole[QZZ]*gqzz22
                     +two*(particleI.quadrupole[QXY]*gqxy22+particleI.quadrupole[QXZ]*gqxz22+particleI.quadrupole[QYZ]*gqyz22))
                      + syk*(particleI.quadrupole[QXX]*gqxx23+particleI.quadrupole[QYY]*gqyy23+particleI.quadrupole[QZZ]*gqzz23
                     +two*(particleI.quadrupole[QXY]*gqxy23+particleI.quadrupole[QXZ]*gqxz23+particleI.quadrupole[QYZ]*gqyz23))
                      + szk*(particleI.quadrupole[QXX]*gqxx24+particleI.quadrupole[QYY]*gqyy24+particleI.quadrupole[QZZ]*gqzz24
                     +two*(particleI.quadrupole[QXY]*gqxy24+particleI.quadrupole[QXZ]*gqxz24+particleI.quadrupole[QYZ]*gqyz24));

    dwkpdr = particleI.charge*(sxk*gux21+syk*guy21+szk*guz21)
                         - particleJ.charge*(sxi*gc22+syi*gc23+szi*gc24)
                      - sxi*(particleJ.quadrupole[QXX]*gqxx22+particleJ.quadrupole[QYY]*gqyy22+particleJ.quadrupole[QZZ]*gqzz22
                     +two*(particleJ.quadrupole[QXY]*gqxy22+particleJ.quadrupole[QXZ]*gqxz22+particleJ.quadrupole[QYZ]*gqyz22))
                      - syi*(particleJ.quadrupole[QXX]*gqxx23+particleJ.quadrupole[QYY]*gqyy23+particleJ.quadrupole[QZZ]*gqzz23
                     +two*(particleJ.quadrupole[QXY]*gqxy23+particleJ.quadrupole[QXZ]*gqxz23+particleJ.quadrupole[QYZ]*gqyz23))
                      - szi*(particleJ.quadrupole[QXX]*gqxx24+particleJ.quadrupole[QYY]*gqyy24+particleJ.quadrupole[QZZ]*gqzz24
                     +two*(particleJ.quadrupole[QXY]*gqxy24+particleJ.quadrupole[QXZ]*gqxz24+particleJ.quadrupole[QYZ]*gqyz24))
                      + sxk*(particleI.quadrupole[QXX]*gux25+particleI.quadrupole[QYY]*gux28+particleI.quadrupole[QZZ]*gux30
                     +two*(particleI.quadrupole[QXY]*gux26+particleI.quadrupole[QXZ]*gux27+particleI.quadrupole[QYZ]*gux29))
                      + syk*(particleI.quadrupole[QXX]*guy25+particleI.quadrupole[QYY]*guy28+particleI.quadrupole[QZZ]*guy30
                     +two*(particleI.quadrupole[QXY]*guy26+particleI.quadrupole[QXZ]*guy27+particleI.quadrupole[QYZ]*guy29))
                      + szk*(particleI.quadrupole[QXX]*guz25+particleI.quadrupole[QYY]*guz28+particleI.quadrupole[QZZ]*guz30
                     +two*(particleI.quadrupole[QXY]*guz26+particleI.quadrupole[QXZ]*guz27+particleI.quadrupole[QYZ]*guz29));

    dsumdr = dsymdr + oneHalf*(dwipdr + dwkpdr);
    dpbi = oneHalf*_bornRadii[jIndex]*dsumdr;
    dpbk = oneHalf*_bornRadii[iIndex]*dsumdr;

    // mutual polarization electrostatic solvation free energy gradient

    if( polarizationType == AmoebaReferenceMultipoleForce::Mutual ){

        dpdx = dpdx - oneHalf *
                           (_inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux5+_inducedDipolePolarS[jIndex][1]*gux6+_inducedDipolePolarS[jIndex][2]*gux7)
                           +_inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy5+_inducedDipolePolarS[jIndex][1]*guy6+_inducedDipolePolarS[jIndex][2]*guy7)
                           +_inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz5+_inducedDipolePolarS[jIndex][1]*guz6+_inducedDipolePolarS[jIndex][2]*guz7)
                           +_inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux5+_inducedDipolePolarS[iIndex][1]*gux6+_inducedDipolePolarS[iIndex][2]*gux7)
                           +_inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy5+_inducedDipolePolarS[iIndex][1]*guy6+_inducedDipolePolarS[iIndex][2]*guy7)
                           +_inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz5+_inducedDipolePolarS[iIndex][1]*guz6+_inducedDipolePolarS[iIndex][2]*guz7));

        dpdy = dpdy - oneHalf *
                           (_inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux6+_inducedDipolePolarS[jIndex][1]*gux8+_inducedDipolePolarS[jIndex][2]*gux9)
                           +_inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy6+_inducedDipolePolarS[jIndex][1]*guy8+_inducedDipolePolarS[jIndex][2]*guy9)
                           +_inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz6+_inducedDipolePolarS[jIndex][1]*guz8+_inducedDipolePolarS[jIndex][2]*guz9)
                           +_inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux6+_inducedDipolePolarS[iIndex][1]*gux8+_inducedDipolePolarS[iIndex][2]*gux9)
                           +_inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy6+_inducedDipolePolarS[iIndex][1]*guy8+_inducedDipolePolarS[iIndex][2]*guy9)
                           +_inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz6+_inducedDipolePolarS[iIndex][1]*guz8+_inducedDipolePolarS[iIndex][2]*guz9));

        dpdz = dpdz - oneHalf *
                           (_inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux7+_inducedDipolePolarS[jIndex][1]*gux9+_inducedDipolePolarS[jIndex][2]*gux10)
                           +_inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy7+_inducedDipolePolarS[jIndex][1]*guy9+_inducedDipolePolarS[jIndex][2]*guy10)
                           +_inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz7+_inducedDipolePolarS[jIndex][1]*guz9+_inducedDipolePolarS[jIndex][2]*guz10)
                           +_inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux7+_inducedDipolePolarS[iIndex][1]*gux9+_inducedDipolePolarS[iIndex][2]*gux10)
                           +_inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy7+_inducedDipolePolarS[iIndex][1]*guy9+_inducedDipolePolarS[iIndex][2]*guy10)
                           +_inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz7+_inducedDipolePolarS[iIndex][1]*guz9+_inducedDipolePolarS[iIndex][2]*guz10));

        duvdr = _inducedDipoleS[iIndex][0]*(_inducedDipolePolarS[jIndex][0]*gux22+_inducedDipolePolarS[jIndex][1]*gux23+_inducedDipolePolarS[jIndex][2]*gux24)
                            + _inducedDipoleS[iIndex][1]*(_inducedDipolePolarS[jIndex][0]*guy22+_inducedDipolePolarS[jIndex][1]*guy23+_inducedDipolePolarS[jIndex][2]*guy24)
                            + _inducedDipoleS[iIndex][2]*(_inducedDipolePolarS[jIndex][0]*guz22+_inducedDipolePolarS[jIndex][1]*guz23+_inducedDipolePolarS[jIndex][2]*guz24)
                            + _inducedDipoleS[jIndex][0]*(_inducedDipolePolarS[iIndex][0]*gux22+_inducedDipolePolarS[iIndex][1]*gux23+_inducedDipolePolarS[iIndex][2]*gux24)
                            + _inducedDipoleS[jIndex][1]*(_inducedDipolePolarS[iIndex][0]*guy22+_inducedDipolePolarS[iIndex][1]*guy23+_inducedDipolePolarS[iIndex][2]*guy24)
                            + _inducedDipoleS[jIndex][2]*(_inducedDipolePolarS[iIndex][0]*guz22+_inducedDipolePolarS[iIndex][1]*guz23+_inducedDipolePolarS[iIndex][2]*guz24);

        dpbi = dpbi - oneHalf*_bornRadii[jIndex]*duvdr;
        dpbk = dpbk - oneHalf*_bornRadii[iIndex]*duvdr;
    }

    // torque due to induced reaction field on permanent dipoles

    RealOpenMM fid1 = oneHalf*(sxk*gux2+syk*guy2+szk*guz2);
    RealOpenMM fid2 = oneHalf*(sxk*gux3+syk*guy3+szk*guz3);
    RealOpenMM fid3 = oneHalf*(sxk*gux4+syk*guy4+szk*guz4);
    RealOpenMM fkd1 = oneHalf*(sxi*gux2+syi*guy2+szi*guz2);
    RealOpenMM fkd2 = oneHalf*(sxi*gux3+syi*guy3+szi*guz3);
    RealOpenMM fkd3 = oneHalf*(sxi*gux4+syi*guy4+szi*guz4);

    RealOpenMM trqi1   = particleI.dipole[1]*fid3 - particleI.dipole[2]*fid2;
    RealOpenMM trqi2   = particleI.dipole[2]*fid1 - particleI.dipole[0]*fid3;
    RealOpenMM trqi3   = particleI.dipole[0]*fid2 - particleI.dipole[1]*fid1;

    RealOpenMM trqi_k1 = particleJ.dipole[1]*fkd3 - particleJ.dipole[2]*fkd2;
    RealOpenMM trqi_k2 = particleJ.dipole[2]*fkd1 - particleJ.dipole[0]*fkd3;
    RealOpenMM trqi_k3 = particleJ.dipole[0]*fkd2 - particleJ.dipole[1]*fkd1;

    // torque due to induced reaction field gradient on quadrupoles;

    RealOpenMM fidg11 = -0.25 *
                              ( (sxk*gqxx2+syk*gqxx3+szk*gqxx4)
                              + (sxk*gux5+syk*guy5+szk*guz5));

    RealOpenMM fidg12 = -0.25 *
                              ( (sxk*gqxy2+syk*gqxy3+szk*gqxy4)
                              + (sxk*gux6+syk*guy6+szk*guz6));

    RealOpenMM fidg13 = -0.25 *
                              ( (sxk*gqxz2+syk*gqxz3+szk*gqxz4)
                              + (sxk*gux7+syk*guy7+szk*guz7));

    RealOpenMM fidg22 = -0.25 *
                              ( (sxk*gqyy2+syk*gqyy3+szk*gqyy4)
                              + (sxk*gux8+syk*guy8+szk*guz8));

    RealOpenMM fidg23 = -0.25 *
                              ( (sxk*gqyz2+syk*gqyz3+szk*gqyz4)
                              + (sxk*gux9+syk*guy9+szk*guz9));

    RealOpenMM fidg33 = -0.25 *
                              ( (sxk*gqzz2+syk*gqzz3+szk*gqzz4)
                              + (sxk*gux10+syk*guy10+szk*guz10));

    RealOpenMM fidg21 = fidg12;
    RealOpenMM fidg31 = fidg13;
    RealOpenMM fidg32 = fidg23;

    RealOpenMM fkdg11 = 0.25 *
                              ( (sxi*gqxx2+syi*gqxx3+szi*gqxx4)
                              + (sxi*gux5+syi*guy5+szi*guz5));

    RealOpenMM fkdg12 = 0.25 *
                              ( (sxi*gqxy2+syi*gqxy3+szi*gqxy4)
                              + (sxi*gux6+syi*guy6+szi*guz6));
    RealOpenMM fkdg13 = 0.25 *
                              ( (sxi*gqxz2+syi*gqxz3+szi*gqxz4)
                              + (sxi*gux7+syi*guy7+szi*guz7));
    RealOpenMM fkdg22 = 0.25 *
                              ( (sxi*gqyy2+syi*gqyy3+szi*gqyy4)
                              + (sxi*gux8+syi*guy8+szi*guz8));
    RealOpenMM fkdg23 = 0.25 *
                              ( (sxi*gqyz2+syi*gqyz3+szi*gqyz4)
                              + (sxi*gux9+syi*guy9+szi*guz9));
    RealOpenMM fkdg33 = 0.25 *
                              ( (sxi*gqzz2+syi*gqzz3+szi*gqzz4)
                              + (sxi*gux10+syi*guy10+szi*guz10));
    RealOpenMM fkdg21 = fkdg12;
    RealOpenMM fkdg31 = fkdg13;
    RealOpenMM fkdg32 = fkdg23;

    trqi1 += two*(particleI.quadrupole[QXY]*fidg13+particleI.quadrupole[QYY]*fidg23+particleI.quadrupole[QYZ]*fidg33
                        -particleI.quadrupole[QXZ]*fidg12-particleI.quadrupole[QYZ]*fidg22-particleI.quadrupole[QZZ]*fidg32);
    trqi2 += two*(particleI.quadrupole[QXZ]*fidg11+particleI.quadrupole[QYZ]*fidg21+particleI.quadrupole[QZZ]*fidg31
                        -particleI.quadrupole[QXX]*fidg13-particleI.quadrupole[QXY]*fidg23-particleI.quadrupole[QXZ]*fidg33);

    trqi3 += two*(particleI.quadrupole[QXX]*fidg12+particleI.quadrupole[QXY]*fidg22+particleI.quadrupole[QXZ]*fidg32
                        -particleI.quadrupole[QXY]*fidg11-particleI.quadrupole[QYY]*fidg21-particleI.quadrupole[QYZ]*fidg31);

    trqi_k1 += two*
                        (particleJ.quadrupole[QXY]*fkdg13+particleJ.quadrupole[QYY]*fkdg23+particleJ.quadrupole[QYZ]*fkdg33
                        -particleJ.quadrupole[QXZ]*fkdg12-particleJ.quadrupole[QYZ]*fkdg22-particleJ.quadrupole[QZZ]*fkdg32);

    trqi_k2 += two*
                        (particleJ.quadrupole[QXZ]*fkdg11+particleJ.quadrupole[QYZ]*fkdg21+particleJ.quadrupole[QZZ]*fkdg31
                        -particleJ.quadrupole[QXX]*fkdg13-particleJ.quadrupole[QXY]*fkdg23-particleJ.quadrupole[QXZ]*fkdg33);

    trqi_k3 += two*
                        (particleJ.quadrupole[QXX]*fkdg12+particleJ.quadrupole[QXY]*fkdg22+particleJ.quadrupole[QXZ]*fkdg32
                        -particleJ.quadrupole[QXY]*fkdg11-particleJ.quadrupole[QYY]*fkdg21-particleJ.quadrupole[QYZ]*fkdg31);

    // total permanent and induced energies for this interaction

    e                        = esym + oneHalf*(ewi+ewk);
    ei                       = oneHalf*(esymi + oneHalf*(ewii+ewki));
    RealOpenMM energy        = e + ei;
    energy                   = iIndex == jIndex ? 0.5*energy : energy;

    forces[iIndex][0]       += (dedx + dpdx);
    forces[iIndex][1]       += (dedy + dpdy);
    forces[iIndex][2]       += (dedz + dpdz);

    forces[jIndex][0]       -= (dedx + dpdx);
    forces[jIndex][1]       -= (dedy + dpdy);
    forces[jIndex][2]       -= (dedz + dpdz);

    torques[iIndex][0]      += (trq1 + trqi1);
    torques[iIndex][1]      += (trq2 + trqi2);
    torques[iIndex][2]      += (trq3 + trqi3);

    dBorn[iIndex]           += drbi + dpbi;
    if( iIndex != jIndex ){

        torques[jIndex][0]  += (trq_k1 + trqi_k1);
        torques[jIndex][1]  += (trq_k2 + trqi_k2);
        torques[jIndex][2]  += (trq_k3 + trqi_k3);

        dBorn[jIndex]       += drbk + dpbk;
    }

    return (energy);
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateElectrostatic( std::vector<MultipoleParticleData>& particleData,
                                                                                     const std::vector<int>& axisTypes,
                                                                                     const std::vector<int>& multipoleAtomZs,
                                                                                     const std::vector<int>& multipoleAtomXs,
                                                                                     const std::vector<int>& multipoleAtomYs,
                                                                                     AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                                     std::vector<RealVec>& torques,
                                                                                     std::vector<RealVec>& forces ) const {

    RealOpenMM energy = AmoebaReferenceMultipoleForce::calculateElectrostatic( particleData, axisTypes,
                                                                               multipoleAtomZs, multipoleAtomXs, multipoleAtomYs,
                                                                               polarizationType, torques, forces );

    vector<RealOpenMM> dBorn;
    initializeRealOpenMMVector( dBorn );

    // Kirkwood loop over particle pairs

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii; jj < particleData.size(); jj++ ){
            energy += calculateKirkwoodPairIxn( polarizationType, particleData[ii], particleData[jj], forces, torques, dBorn );
        }
    }

    // cavity term

    if( getIncludeCavityTerm() ){
        energy += calculateCavityTermEnergyAndForces( dBorn );
    }

    // apply Born chain rule; skip diagonal terms since these make no contribution to forces

    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){
            calculateGrycukChainRulePairIxn( particleData[ii], particleData[jj], dBorn, forces );
            calculateGrycukChainRulePairIxn( particleData[jj], particleData[ii], dBorn, forces );
        }
    }

    // correct vacuum to SCRF derivatives (ediff1 in TINKER)

    vector<RealOpenMM> scaleFactors(LAST_SCALE_TYPE_INDEX);
    for( unsigned int kk = 0; kk < scaleFactors.size(); kk++ ){
        scaleFactors[kk] = 1.0;
    }   

    RealOpenMM eDiffEnergy = 0.0;
    for( unsigned int ii = 0; ii < particleData.size(); ii++ ){
        for( unsigned int jj = ii+1; jj < particleData.size(); jj++ ){

            if( jj <= _maxScaleIndex[ii] ){
                getMultipoleScaleFactors( ii, jj, scaleFactors);
            }

            eDiffEnergy += calculateKirkwoodEDiffPairIxn( polarizationType, particleData[ii], particleData[jj],
                                                          scaleFactors[P_SCALE], scaleFactors[D_SCALE], forces, torques );

            if( jj <= _maxScaleIndex[ii] ){
                for( unsigned int kk = 0; kk < LAST_SCALE_TYPE_INDEX; kk++ ){
                    scaleFactors[kk] = 1.0;
                }
            }
        }
    }
    energy += (_electric/_dielectric)*eDiffEnergy;

    return energy;
}

void AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateGrycukChainRulePairIxn( const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                                                        const std::vector<RealOpenMM>& dBorn, std::vector<RealVec>& forces ) const {

    static const RealOpenMM zero        = 0.0; 
    static const RealOpenMM one         = 1.0; 
    static const RealOpenMM two         = 2.0; 
    static const RealOpenMM three       = 3.0; 
    static const RealOpenMM four        = 4.0; 
    static const RealOpenMM five        = 5.0; 
    static const RealOpenMM six         = 6.0; 
    static const RealOpenMM nine        = 9.0; 
    static const RealOpenMM seventeen   = 17.0; 
    static const RealOpenMM oneThird    = 1.0/3.0; 
    static const RealOpenMM oneFourth   = 0.25; 
    static const RealOpenMM pi          = 3.1415926535897;

    unsigned int iIndex                 = particleI.particleIndex;
    unsigned int jIndex                 = particleJ.particleIndex;

    RealOpenMM pi43                     = four*oneThird*pi;

    RealOpenMM lik, uik;
    RealOpenMM lik4, uik4;
    RealOpenMM factor                   = -POW(pi,oneThird)*powf(six,(two*oneThird))/nine;
    RealOpenMM term                     = pi43/(_bornRadii[iIndex]*_bornRadii[iIndex]*_bornRadii[iIndex]);
          term                          = factor/POW( term, (four*oneThird) );

    RealVec deltaR                      = particleJ.position - particleI.position;

    RealOpenMM sk                       = _scaledRadii[jIndex];
    RealOpenMM sk2                      = sk*sk;
    RealOpenMM r2                       = deltaR.dot( deltaR );
    RealOpenMM r                        = SQRT(r2);
    RealOpenMM de                       = zero;

    if( (_atomicRadii[iIndex] + r) < sk ){
        RealOpenMM uik4;
        uik                = sk - r;
        uik4               = uik*uik;
        uik4               = uik4*uik4;
        de                 = -four*pi/uik4;
    }

    if( (_atomicRadii[iIndex] + r) < sk){
        lik          = sk - r;
        lik4         = lik*lik;
        lik4         = lik4*lik4;
        de          += oneFourth*pi*(sk2-four*sk*r+seventeen*r2)/ (r2*lik4);
    } else if( r < (_atomicRadii[iIndex] +sk) ){
        lik          = _atomicRadii[iIndex];
        lik4         = lik*lik;
        lik4         = lik4*lik4;
        de          += oneFourth*pi*(two*_atomicRadii[iIndex]*_atomicRadii[iIndex]-sk2-r2)/ (r2*lik4);
    } else {
        lik          = r - sk;
        lik4         = lik*lik;
        lik4         = lik4*lik4;
        de          += oneFourth*pi*(sk2-four*sk*r+r2)/ (r2*lik4);
    }
    uik                = r + sk;
    uik4               = uik*uik;
    uik4               = uik4*uik4;
    
    de                -= oneFourth*pi*(sk2+four*sk*r+r2)/ (r2*uik4);
    RealOpenMM dbr     = term*de/r;
          de           = dbr*dBorn[iIndex];

    forces[iIndex]    -= deltaR*de;
    forces[jIndex]    += deltaR*de;

    return;
}

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateCavityTermEnergyAndForces( std::vector<RealOpenMM>& dBorn ) const {

    RealOpenMM energy       = 0.0;
    for( unsigned int ii = 0; ii < _numParticles; ii++ ){
        RealOpenMM r        = _atomicRadii[ii] + _dielectricOffset + _probeRadius;
        RealOpenMM ratio    = ( (_atomicRadii[ii] + _dielectricOffset)/_bornRadii[ii]);
        RealOpenMM saTerm   = _surfaceAreaFactor*r*r*POW( ratio, 6.0 );
        dBorn[ii]          += saTerm/_bornRadii[ii];
        energy             += saTerm;
    }
    return (energy/-6.0);
}

/*****************************************************************************

       ediff1 correct vacuum to SCRF derivatives

       calculates the energy and derivatives of polarizing
       the vacuum induced dipoles to their SCRF polarized values

******************************************************************************/

RealOpenMM AmoebaReferenceGeneralizedKirkwoodMultipoleForce::calculateKirkwoodEDiffPairIxn( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                                            const MultipoleParticleData& particleI,
                                                                                            const MultipoleParticleData& particleJ,
                                                                                            RealOpenMM pscale, RealOpenMM dscale,
                                                                                            std::vector<RealVec>& forces,
                                                                                            std::vector<RealVec>& torques ) const {

    static const RealOpenMM zero        = 0.0; 
    static const RealOpenMM one         = 1.0; 
    static const RealOpenMM two         = 2.0; 
    static const RealOpenMM three       = 3.0; 
    static const RealOpenMM four        = 4.0; 
    static const RealOpenMM five        = 5.0; 
    static const RealOpenMM six         = 6.0; 
    static const RealOpenMM seven       = 7.0; 
    static const RealOpenMM nine        = 9.0; 
    static const RealOpenMM seventeen   = 17.0; 
    static const RealOpenMM oneHalf     = 0.5;
    static const RealOpenMM oneThird    = 1.0/3.0; 
    static const RealOpenMM oneFourth   = 0.25; 
    static const RealOpenMM pi          = 3.1415926535897;

    static const RealOpenMM uscale      = 1.0;

    RealOpenMM scale3,scale5;
    RealOpenMM scale7;
    RealOpenMM psc3,psc5,psc7;
    RealOpenMM dsc3,dsc5,dsc7;
    RealOpenMM scale3i,scale5i;
    RealOpenMM r,rr1,rr3;
    RealOpenMM rr5,rr7,rr9;
    RealOpenMM pgamma;

    unsigned int iIndex   = particleI.particleIndex;
    unsigned int jIndex   = particleJ.particleIndex;

    // deltaR

    RealVec deltaR        = particleJ.position - particleI.position;
    RealOpenMM r2         = deltaR.dot( deltaR );

    RealOpenMM xr         = deltaR[0];
    RealOpenMM yr         = deltaR[1];
    RealOpenMM zr         = deltaR[2];

	 r                     = SQRT(r2);
    rr1                   = one / r;
    rr3                   = rr1 / r2;
    rr5                   = three*rr3 / r2;
    rr7                   = five*rr5 / r2;
    rr9                   = seven*rr7 / r2;

    scale3                = one;
    scale5                = one;
    scale7                = one;

    RealOpenMM ddsc3_1    = zero;
    RealOpenMM ddsc3_2    = zero;
    RealOpenMM ddsc3_3    = zero;

    RealOpenMM ddsc5_1    = zero;
    RealOpenMM ddsc5_2    = zero;
    RealOpenMM ddsc5_3    = zero;

    RealOpenMM ddsc7_1    = zero;
    RealOpenMM ddsc7_2    = zero;
    RealOpenMM ddsc7_3    = zero;

    // apply Thole polarization damping to scale factors
 
    RealOpenMM damp = particleI.dampingFactor*particleJ.dampingFactor;
    if( damp != zero ){
        pgamma          = particleJ.thole > particleI.thole ? particleI.thole : particleJ.thole;
        RealOpenMM ratio     = (r/damp);
        damp            = -pgamma*ratio*ratio*ratio;
        if( damp > -50.0 ){
            RealOpenMM dampE  = expf( damp );
            RealOpenMM damp2  = damp*damp;
            scale3       = one - dampE;
            scale5       = one - (one - damp)*dampE;
            scale7       = one - (one - damp + 0.6*damp2)*dampE;
            //scale9       = one - (one - damp + (18.0f*damp2 - (9.0f*damp*damp2))/3five)*dampE;

            ddsc3_1     = -three*damp*expf(damp)*xr/r2;
            ddsc3_2     = -three*damp*expf(damp)*yr/r2;
            ddsc3_3     = -three*damp*expf(damp)*zr/r2;

            ddsc5_1     = -damp*ddsc3_1;
            ddsc5_2     = -damp*ddsc3_2;
            ddsc5_3     = -damp*ddsc3_3;

            ddsc7_1     = (-0.2-0.6*damp)*ddsc5_1;
            ddsc7_2     = (-0.2-0.6*damp)*ddsc5_2;
            ddsc7_3     = (-0.2-0.6*damp)*ddsc5_3;
        }
    }

    scale3i             = scale3*uscale;
    scale5i             = scale5*uscale;

    dsc3                = scale3*dscale;
    dsc5                = scale5*dscale;
    dsc7                = scale7*dscale;

    psc3                = scale3*pscale;
    psc5                = scale5*pscale;
    psc7                = scale7*pscale;
 
    // construct auxiliary vectors for permanent terms
 
    RealOpenMM dixr1             = particleI.dipole[1]*zr - particleI.dipole[2]*yr;
    RealOpenMM dixr2             = particleI.dipole[2]*xr - particleI.dipole[0]*zr;
    RealOpenMM dixr3             = particleI.dipole[0]*yr - particleI.dipole[1]*xr;

    RealOpenMM dkxr1             = particleJ.dipole[1]*zr - particleJ.dipole[2]*yr;
    RealOpenMM dkxr2             = particleJ.dipole[2]*xr - particleJ.dipole[0]*zr;
    RealOpenMM dkxr3             = particleJ.dipole[0]*yr - particleJ.dipole[1]*xr;

    RealOpenMM qir1              = particleI.quadrupole[QXX]*xr + particleI.quadrupole[QXY]*yr + particleI.quadrupole[QXZ]*zr;
    RealOpenMM qir2              = particleI.quadrupole[QXY]*xr + particleI.quadrupole[QYY]*yr + particleI.quadrupole[QYZ]*zr;
    RealOpenMM qir3              = particleI.quadrupole[QXZ]*xr + particleI.quadrupole[QYZ]*yr + particleI.quadrupole[QZZ]*zr;

    RealOpenMM qkr1              = particleJ.quadrupole[QXX]*xr + particleJ.quadrupole[QXY]*yr + particleJ.quadrupole[QXZ]*zr;
    RealOpenMM qkr2              = particleJ.quadrupole[QXY]*xr + particleJ.quadrupole[QYY]*yr + particleJ.quadrupole[QYZ]*zr;
    RealOpenMM qkr3              = particleJ.quadrupole[QXZ]*xr + particleJ.quadrupole[QYZ]*yr + particleJ.quadrupole[QZZ]*zr;

    RealOpenMM rxqir1            = yr*qir3 - zr*qir2;
    RealOpenMM rxqir2            = zr*qir1 - xr*qir3;
    RealOpenMM rxqir3            = xr*qir2 - yr*qir1;

    RealOpenMM rxqkr1            = yr*qkr3 - zr*qkr2;
    RealOpenMM rxqkr2            = zr*qkr1 - xr*qkr3;
    RealOpenMM rxqkr3            = xr*qkr2 - yr*qkr1;

    // get intermediate variables for permanent energy terms
 
    RealOpenMM sc3               = particleI.dipole[0]*xr  + particleI.dipole[1]*yr  + particleI.dipole[2]*zr;
    RealOpenMM sc4               = particleJ.dipole[0]*xr  + particleJ.dipole[1]*yr  + particleJ.dipole[2]*zr;
    RealOpenMM sc5               = qir1*xr + qir2*yr + qir3*zr;
    RealOpenMM sc6               = qkr1*xr + qkr2*yr + qkr3*zr;
 
    // construct auxiliary vectors for induced terms
 
    RealOpenMM dixuk1            = particleI.dipole[1]*_inducedDipoleS[jIndex][2] - particleI.dipole[2]*_inducedDipoleS[jIndex][1];
    RealOpenMM dixuk2            = particleI.dipole[2]*_inducedDipoleS[jIndex][0] - particleI.dipole[0]*_inducedDipoleS[jIndex][2];
    RealOpenMM dixuk3            = particleI.dipole[0]*_inducedDipoleS[jIndex][1] - particleI.dipole[1]*_inducedDipoleS[jIndex][0];

    RealOpenMM dkxui1            = particleJ.dipole[1]*_inducedDipoleS[iIndex][2] - particleJ.dipole[2]*_inducedDipoleS[iIndex][1];
    RealOpenMM dkxui2            = particleJ.dipole[2]*_inducedDipoleS[iIndex][0] - particleJ.dipole[0]*_inducedDipoleS[iIndex][2];
    RealOpenMM dkxui3            = particleJ.dipole[0]*_inducedDipoleS[iIndex][1] - particleJ.dipole[1]*_inducedDipoleS[iIndex][0];

    RealOpenMM dixukp1           = particleI.dipole[1]*_inducedDipolePolarS[jIndex][2] - particleI.dipole[2]*_inducedDipolePolarS[jIndex][1];
    RealOpenMM dixukp2           = particleI.dipole[2]*_inducedDipolePolarS[jIndex][0] - particleI.dipole[0]*_inducedDipolePolarS[jIndex][2];
    RealOpenMM dixukp3           = particleI.dipole[0]*_inducedDipolePolarS[jIndex][1] - particleI.dipole[1]*_inducedDipolePolarS[jIndex][0];

    RealOpenMM dkxuip1           = particleJ.dipole[1]*_inducedDipolePolarS[iIndex][2] - particleJ.dipole[2]*_inducedDipolePolarS[iIndex][1];
    RealOpenMM dkxuip2           = particleJ.dipole[2]*_inducedDipolePolarS[iIndex][0] - particleJ.dipole[0]*_inducedDipolePolarS[iIndex][2];
    RealOpenMM dkxuip3           = particleJ.dipole[0]*_inducedDipolePolarS[iIndex][1] - particleJ.dipole[1]*_inducedDipolePolarS[iIndex][0];

    RealOpenMM qiuk1             = particleI.quadrupole[QXX]*_inducedDipoleS[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipoleS[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipoleS[jIndex][2];
    RealOpenMM qiuk2             = particleI.quadrupole[QXY]*_inducedDipoleS[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipoleS[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipoleS[jIndex][2];
    RealOpenMM qiuk3             = particleI.quadrupole[QXZ]*_inducedDipoleS[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipoleS[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipoleS[jIndex][2];

    RealOpenMM qkui1             = particleJ.quadrupole[QXX]*_inducedDipoleS[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipoleS[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipoleS[iIndex][2];
    RealOpenMM qkui2             = particleJ.quadrupole[QXY]*_inducedDipoleS[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipoleS[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipoleS[iIndex][2];
    RealOpenMM qkui3             = particleJ.quadrupole[QXZ]*_inducedDipoleS[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipoleS[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipoleS[iIndex][2];

    RealOpenMM qiukp1            = particleI.quadrupole[QXX]*_inducedDipolePolarS[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipolePolarS[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipolePolarS[jIndex][2];
    RealOpenMM qiukp2            = particleI.quadrupole[QXY]*_inducedDipolePolarS[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipolePolarS[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipolePolarS[jIndex][2];
    RealOpenMM qiukp3            = particleI.quadrupole[QXZ]*_inducedDipolePolarS[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipolePolarS[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipolePolarS[jIndex][2];

    RealOpenMM qkuip1            = particleJ.quadrupole[QXX]*_inducedDipolePolarS[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipolePolarS[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipolePolarS[iIndex][2];
    RealOpenMM qkuip2            = particleJ.quadrupole[QXY]*_inducedDipolePolarS[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipolePolarS[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipolePolarS[iIndex][2];
    RealOpenMM qkuip3            = particleJ.quadrupole[QXZ]*_inducedDipolePolarS[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipolePolarS[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipolePolarS[iIndex][2];

    RealOpenMM uixqkr1           = _inducedDipoleS[iIndex][1]*qkr3 - _inducedDipoleS[iIndex][2]*qkr2;
    RealOpenMM uixqkr2           = _inducedDipoleS[iIndex][2]*qkr1 - _inducedDipoleS[iIndex][0]*qkr3;
    RealOpenMM uixqkr3           = _inducedDipoleS[iIndex][0]*qkr2 - _inducedDipoleS[iIndex][1]*qkr1;

    RealOpenMM ukxqir1           = _inducedDipoleS[jIndex][1]*qir3 - _inducedDipoleS[jIndex][2]*qir2;
    RealOpenMM ukxqir2           = _inducedDipoleS[jIndex][2]*qir1 - _inducedDipoleS[jIndex][0]*qir3;
    RealOpenMM ukxqir3           = _inducedDipoleS[jIndex][0]*qir2 - _inducedDipoleS[jIndex][1]*qir1;

    RealOpenMM uixqkrp1          = _inducedDipolePolarS[iIndex][1]*qkr3 - _inducedDipolePolarS[iIndex][2]*qkr2;
    RealOpenMM uixqkrp2          = _inducedDipolePolarS[iIndex][2]*qkr1 - _inducedDipolePolarS[iIndex][0]*qkr3;
    RealOpenMM uixqkrp3          = _inducedDipolePolarS[iIndex][0]*qkr2 - _inducedDipolePolarS[iIndex][1]*qkr1;

    RealOpenMM ukxqirp1          = _inducedDipolePolarS[jIndex][1]*qir3 - _inducedDipolePolarS[jIndex][2]*qir2;
    RealOpenMM ukxqirp2          = _inducedDipolePolarS[jIndex][2]*qir1 - _inducedDipolePolarS[jIndex][0]*qir3;
    RealOpenMM ukxqirp3          = _inducedDipolePolarS[jIndex][0]*qir2 - _inducedDipolePolarS[jIndex][1]*qir1;

    RealOpenMM rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    RealOpenMM rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    RealOpenMM rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    RealOpenMM rxqkui1           = yr*qkui3 - zr*qkui2;
    RealOpenMM rxqkui2           = zr*qkui1 - xr*qkui3;
    RealOpenMM rxqkui3           = xr*qkui2 - yr*qkui1;

    RealOpenMM rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    RealOpenMM rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    RealOpenMM rxqiukp3          = xr*qiukp2 - yr*qiukp1;

    RealOpenMM rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    RealOpenMM rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    RealOpenMM rxqkuip3          = xr*qkuip2 - yr*qkuip1;
 
    // get intermediate variables for induction energy terms

    RealOpenMM sci1              = _inducedDipoleS[iIndex][0]*particleJ.dipole[0] + _inducedDipoleS[iIndex][1]*particleJ.dipole[1] +
                              _inducedDipoleS[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipoleS[jIndex][0] +
                              particleI.dipole[1]*_inducedDipoleS[jIndex][1] + particleI.dipole[2]*_inducedDipoleS[jIndex][2];

    RealOpenMM sci3              = _inducedDipoleS[iIndex][0]*xr + _inducedDipoleS[iIndex][1]*yr + _inducedDipoleS[iIndex][2]*zr;
    RealOpenMM sci4              = _inducedDipoleS[jIndex][0]*xr + _inducedDipoleS[jIndex][1]*yr + _inducedDipoleS[jIndex][2]*zr;

    RealOpenMM sci7              = qir1*_inducedDipoleS[jIndex][0] + qir2*_inducedDipoleS[jIndex][1] + qir3*_inducedDipoleS[jIndex][2];
    RealOpenMM sci8              = qkr1*_inducedDipoleS[iIndex][0] + qkr2*_inducedDipoleS[iIndex][1] + qkr3*_inducedDipoleS[iIndex][2];

    RealOpenMM scip1             = _inducedDipolePolarS[iIndex][0]*particleJ.dipole[0] + _inducedDipolePolarS[iIndex][1]*particleJ.dipole[1] +
                              _inducedDipolePolarS[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipolePolarS[jIndex][0] +
                              particleI.dipole[1]*_inducedDipolePolarS[jIndex][1] + particleI.dipole[2]*_inducedDipolePolarS[jIndex][2];

    RealOpenMM scip2             = _inducedDipoleS[iIndex][0]*_inducedDipolePolarS[jIndex][0] + _inducedDipoleS[iIndex][1]*_inducedDipolePolarS[jIndex][1] +
                              _inducedDipoleS[iIndex][2]*_inducedDipolePolarS[jIndex][2] + _inducedDipolePolarS[iIndex][0]*_inducedDipoleS[jIndex][0] +
                              _inducedDipolePolarS[iIndex][1]*_inducedDipoleS[jIndex][1] + _inducedDipolePolarS[iIndex][2]*_inducedDipoleS[jIndex][2];

    RealOpenMM scip3             = _inducedDipolePolarS[iIndex][0]*xr + _inducedDipolePolarS[iIndex][1]*yr + _inducedDipolePolarS[iIndex][2]*zr;
    RealOpenMM scip4             = _inducedDipolePolarS[jIndex][0]*xr + _inducedDipolePolarS[jIndex][1]*yr + _inducedDipolePolarS[jIndex][2]*zr;

    RealOpenMM scip7             = qir1*_inducedDipolePolarS[jIndex][0] + qir2*_inducedDipolePolarS[jIndex][1] + qir3*_inducedDipolePolarS[jIndex][2];
    RealOpenMM scip8             = qkr1*_inducedDipolePolarS[iIndex][0] + qkr2*_inducedDipolePolarS[iIndex][1] + qkr3*_inducedDipolePolarS[iIndex][2];

    // calculate the gl functions for potential energy

    RealOpenMM gli1              = particleJ.charge*sci3 - particleI.charge*sci4;
    RealOpenMM gli2              = -sc3*sci4 - sci3*sc4;
    RealOpenMM gli3              = sci3*sc6 - sci4*sc5;
    RealOpenMM gli6              = sci1;
    RealOpenMM gli7              = two*(sci7-sci8);
    RealOpenMM glip1             = particleJ.charge*scip3 - particleI.charge*scip4;
    RealOpenMM glip2             = -sc3*scip4 - scip3*sc4;
    RealOpenMM glip3             = scip3*sc6 - scip4*sc5;
    RealOpenMM glip6             = scip1;
    RealOpenMM glip7             = two*(scip7-scip8);

    // get the permanent multipole and induced energies;

    RealOpenMM energy            = oneHalf*(rr3*(gli1+gli6)*psc3 + rr5*(gli2+gli7)*psc5 + rr7*gli3*psc7);

    // intermediate variables for the induced-permanent terms

    RealOpenMM gfi1              = oneHalf*rr5*((gli1+gli6)*psc3 +
                                        (glip1+glip6)*dsc3+scip2*scale3i) +
                              oneHalf*rr7*((gli7+gli2)*psc5 +
                                        (glip7+glip2)*dsc5 -
                              (sci3*scip4+scip3*sci4)*scale5i) +
                              oneHalf*rr9*(gli3*psc7+glip3*dsc7);

    RealOpenMM gfi4              = two*rr5;
    RealOpenMM gfi5              = rr7*(sci4*psc7+scip4*dsc7);
    RealOpenMM gfi6              = -rr7*(sci3*psc7+scip3*dsc7);

    // get the induced force;
 
    RealOpenMM ftm2i1            = gfi1*xr + oneHalf*
                          (- rr3*particleJ.charge*(_inducedDipoleS[iIndex][0]*psc3+_inducedDipolePolarS[iIndex][0]*dsc3)
                           + rr5*sc4*(_inducedDipoleS[iIndex][0]*psc5+_inducedDipolePolarS[iIndex][0]*dsc5)
                           - rr7*sc6*(_inducedDipoleS[iIndex][0]*psc7+_inducedDipolePolarS[iIndex][0]*dsc7))
                           +(rr3*particleI.charge*(_inducedDipoleS[jIndex][0]*psc3+_inducedDipolePolarS[jIndex][0]*dsc3)
                           + rr5*sc3*(_inducedDipoleS[jIndex][0]*psc5+_inducedDipolePolarS[jIndex][0]*dsc5)
                           + rr7*sc5*(_inducedDipoleS[jIndex][0]*psc7+_inducedDipolePolarS[jIndex][0]*dsc7))*oneHalf
                           + rr5*scale5i*(sci4*_inducedDipolePolarS[iIndex][0]+scip4*_inducedDipoleS[iIndex][0]
                           + sci3*_inducedDipolePolarS[jIndex][0]+scip3*_inducedDipoleS[jIndex][0])*oneHalf
                           + oneHalf*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[0]
                           + oneHalf*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[0]
                           + oneHalf*gfi4*((qkui1-qiuk1)*psc5
                           + (qkuip1-qiukp1)*dsc5)
                           + gfi5*qir1 + gfi6*qkr1;
 
    RealOpenMM ftm2i2            = gfi1*yr + oneHalf*
                          (- rr3*particleJ.charge*(_inducedDipoleS[iIndex][1]*psc3+_inducedDipolePolarS[iIndex][1]*dsc3)
                           + rr5*sc4*(_inducedDipoleS[iIndex][1]*psc5+_inducedDipolePolarS[iIndex][1]*dsc5)
                           - rr7*sc6*(_inducedDipoleS[iIndex][1]*psc7+_inducedDipolePolarS[iIndex][1]*dsc7))
                           +(rr3*particleI.charge*(_inducedDipoleS[jIndex][1]*psc3+_inducedDipolePolarS[jIndex][1]*dsc3)
                           + rr5*sc3*(_inducedDipoleS[jIndex][1]*psc5+_inducedDipolePolarS[jIndex][1]*dsc5)
                           + rr7*sc5*(_inducedDipoleS[jIndex][1]*psc7+_inducedDipolePolarS[jIndex][1]*dsc7))*oneHalf
                           + rr5*scale5i*(sci4*_inducedDipolePolarS[iIndex][1]+scip4*_inducedDipoleS[iIndex][1]
                           + sci3*_inducedDipolePolarS[jIndex][1]+scip3*_inducedDipoleS[jIndex][1])*oneHalf
                           + oneHalf*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[1]
                           + oneHalf*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[1]
                           + oneHalf*gfi4*((qkui2-qiuk2)*psc5
                           + (qkuip2-qiukp2)*dsc5)
                           + gfi5*qir2 + gfi6*qkr2;

    RealOpenMM ftm2i3            = gfi1*zr  + oneHalf*
                          (- rr3*particleJ.charge*(_inducedDipoleS[iIndex][2]*psc3+_inducedDipolePolarS[iIndex][2]*dsc3)
                           + rr5*sc4*(_inducedDipoleS[iIndex][2]*psc5+_inducedDipolePolarS[iIndex][2]*dsc5)
                           - rr7*sc6*(_inducedDipoleS[iIndex][2]*psc7+_inducedDipolePolarS[iIndex][2]*dsc7))
                           +(rr3*particleI.charge*(_inducedDipoleS[jIndex][2]*psc3+_inducedDipolePolarS[jIndex][2]*dsc3)
                           + rr5*sc3*(_inducedDipoleS[jIndex][2]*psc5+_inducedDipolePolarS[jIndex][2]*dsc5)
                           + rr7*sc5*(_inducedDipoleS[jIndex][2]*psc7+_inducedDipolePolarS[jIndex][2]*dsc7))*oneHalf
                           + rr5*scale5i*(sci4*_inducedDipolePolarS[iIndex][2]+scip4*_inducedDipoleS[iIndex][2]
                           + sci3*_inducedDipolePolarS[jIndex][2]+scip3*_inducedDipoleS[jIndex][2])*oneHalf
                           + oneHalf*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[2]
                           + oneHalf*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[2]
                           + oneHalf*gfi4*((qkui3-qiuk3)*psc5
                           + (qkuip3-qiukp3)*dsc5)
                           + gfi5*qir3 + gfi6*qkr3;
 
    // intermediate values needed for partially excluded interactions

    RealOpenMM fridmp1           = oneHalf*(rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_1
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_1
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_1);

    RealOpenMM fridmp2           = oneHalf*(rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_2
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_2
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_2);

    RealOpenMM fridmp3           = oneHalf*(rr3*((gli1+gli6)*pscale
                          +(glip1+glip6)*dscale)*ddsc3_3
                          + rr5*((gli2+gli7)*pscale
                          +(glip2+glip7)*dscale)*ddsc5_3
                          + rr7*(gli3*pscale+glip3*dscale)*ddsc7_3);

    // get the induced-induced derivative terms
 
    RealOpenMM findmp1           = oneHalf*uscale*(scip2*rr3*ddsc3_1
                          - rr5*ddsc5_1*(sci3*scip4+scip3*sci4));

    RealOpenMM findmp2           = oneHalf*uscale*(scip2*rr3*ddsc3_2
                          - rr5*ddsc5_2*(sci3*scip4+scip3*sci4));

    RealOpenMM findmp3           = oneHalf*uscale*(scip2*rr3*ddsc3_3
                          - rr5*ddsc5_3*(sci3*scip4+scip3*sci4));

    // handle of scaling for partially excluded interactions

    ftm2i1           -= fridmp1 + findmp1;
    ftm2i2           -= fridmp2 + findmp2;
    ftm2i3           -= fridmp3 + findmp3;

    // correction to convert mutual to direct polarization force

    if( polarizationType == AmoebaReferenceMultipoleForce::Direct ){
        RealOpenMM gfd      = oneHalf*(rr5*scip2*scale3i - rr7*(scip3*sci4+sci3*scip4)*scale5i);
        RealOpenMM fdir1    = gfd*xr + oneHalf*rr5*scale5i* (sci4*_inducedDipolePolarS[iIndex][0]+scip4*_inducedDipoleS[iIndex][0] + sci3*_inducedDipolePolarS[jIndex][0]+scip3*_inducedDipoleS[jIndex][0]);
        RealOpenMM fdir2    = gfd*yr + oneHalf*rr5*scale5i* (sci4*_inducedDipolePolarS[iIndex][1]+scip4*_inducedDipoleS[iIndex][1] + sci3*_inducedDipolePolarS[jIndex][1]+scip3*_inducedDipoleS[jIndex][1]);
        RealOpenMM fdir3    = gfd*zr + oneHalf*rr5*scale5i* (sci4*_inducedDipolePolarS[iIndex][2]+scip4*_inducedDipoleS[iIndex][2] + sci3*_inducedDipolePolarS[jIndex][2]+scip3*_inducedDipoleS[jIndex][2]);
        ftm2i1        -= fdir1 - findmp1;
        ftm2i2        -= fdir2 - findmp2;
        ftm2i3        -= fdir3 - findmp3;
     
    }

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and k
 
    RealOpenMM gti2              = oneHalf*(sci4*psc5+scip4*dsc5)*rr5;
    RealOpenMM gti3              = oneHalf*(sci3*psc5+scip3*dsc5)*rr5;
    RealOpenMM gti4              = gfi4;
    RealOpenMM gti5              = gfi5;
    RealOpenMM gti6              = gfi6;

    // calculate the induced torque components

    RealOpenMM ttm2i1            = -rr3*(dixuk1*psc3+dixukp1*dsc3)*oneHalf
                          + gti2*dixr1 + gti4*((ukxqir1+rxqiuk1)*psc5
                          +(ukxqirp1+rxqiukp1)*dsc5)*oneHalf - gti5*rxqir1;

    RealOpenMM ttm2i2            = -rr3*(dixuk2*psc3+dixukp2*dsc3)*oneHalf
                          + gti2*dixr2 + gti4*((ukxqir2+rxqiuk2)*psc5
                          +(ukxqirp2+rxqiukp2)*dsc5)*oneHalf - gti5*rxqir2;

    RealOpenMM ttm2i3            = -rr3*(dixuk3*psc3+dixukp3*dsc3)*oneHalf
                          + gti2*dixr3 + gti4*((ukxqir3+rxqiuk3)*psc5
                          +(ukxqirp3+rxqiukp3)*dsc5)*oneHalf - gti5*rxqir3;

    RealOpenMM ttm3i1            = -rr3*(dkxui1*psc3+dkxuip1*dsc3)*oneHalf
                          + gti3*dkxr1 - gti4*((uixqkr1+rxqkui1)*psc5
                          +(uixqkrp1+rxqkuip1)*dsc5)*oneHalf - gti6*rxqkr1;

    RealOpenMM ttm3i2            = -rr3*(dkxui2*psc3+dkxuip2*dsc3)*oneHalf
                          + gti3*dkxr2 - gti4*((uixqkr2+rxqkui2)*psc5
                          +(uixqkrp2+rxqkuip2)*dsc5)*oneHalf - gti6*rxqkr2;

    RealOpenMM ttm3i3            = -rr3*(dkxui3*psc3+dkxuip3*dsc3)*oneHalf
                          + gti3*dkxr3 - gti4*((uixqkr3+rxqkui3)*psc5
                          +(uixqkrp3+rxqkuip3)*dsc5)*oneHalf - gti6*rxqkr3;
 
    // update force and torque on site k

    RealVec force, torqueI, torqueJ;
    force[0]      = -ftm2i1;
    force[1]      = -ftm2i2;
    force[2]      = -ftm2i3;

    torqueI[0]    = ttm2i1;
    torqueI[1]    = ttm2i2;
    torqueI[2]    = ttm2i3;

    torqueJ[0]    = ttm3i1;
    torqueJ[1]    = ttm3i2;
    torqueJ[2]    = ttm3i3;

    // construct auxiliary vectors for induced terms

    dixuk1            = particleI.dipole[1]*_inducedDipole[jIndex][2] - particleI.dipole[2]*_inducedDipole[jIndex][1];
    dixuk2            = particleI.dipole[2]*_inducedDipole[jIndex][0] - particleI.dipole[0]*_inducedDipole[jIndex][2];
    dixuk3            = particleI.dipole[0]*_inducedDipole[jIndex][1] - particleI.dipole[1]*_inducedDipole[jIndex][0];

    dkxui1            = particleJ.dipole[1]*_inducedDipole[iIndex][2] - particleJ.dipole[2]*_inducedDipole[iIndex][1];
    dkxui2            = particleJ.dipole[2]*_inducedDipole[iIndex][0] - particleJ.dipole[0]*_inducedDipole[iIndex][2];
    dkxui3            = particleJ.dipole[0]*_inducedDipole[iIndex][1] - particleJ.dipole[1]*_inducedDipole[iIndex][0];

    dixukp1           = particleI.dipole[1]*_inducedDipolePolar[jIndex][2] - particleI.dipole[2]*_inducedDipolePolar[jIndex][1];
    dixukp2           = particleI.dipole[2]*_inducedDipolePolar[jIndex][0] - particleI.dipole[0]*_inducedDipolePolar[jIndex][2];
    dixukp3           = particleI.dipole[0]*_inducedDipolePolar[jIndex][1] - particleI.dipole[1]*_inducedDipolePolar[jIndex][0];

    dkxuip1           = particleJ.dipole[1]*_inducedDipolePolar[iIndex][2] - particleJ.dipole[2]*_inducedDipolePolar[iIndex][1];
    dkxuip2           = particleJ.dipole[2]*_inducedDipolePolar[iIndex][0] - particleJ.dipole[0]*_inducedDipolePolar[iIndex][2];
    dkxuip3           = particleJ.dipole[0]*_inducedDipolePolar[iIndex][1] - particleJ.dipole[1]*_inducedDipolePolar[iIndex][0];

    qiuk1             = particleI.quadrupole[QXX]*_inducedDipole[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipole[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipole[jIndex][2];
    qiuk2             = particleI.quadrupole[QXY]*_inducedDipole[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipole[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipole[jIndex][2];
    qiuk3             = particleI.quadrupole[QXZ]*_inducedDipole[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipole[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipole[jIndex][2];

    qkui1             = particleJ.quadrupole[QXX]*_inducedDipole[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipole[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipole[iIndex][2];
    qkui2             = particleJ.quadrupole[QXY]*_inducedDipole[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipole[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipole[iIndex][2];
    qkui3             = particleJ.quadrupole[QXZ]*_inducedDipole[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipole[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipole[iIndex][2];

    qiukp1            = particleI.quadrupole[QXX]*_inducedDipolePolar[jIndex][0] + particleI.quadrupole[QXY]*_inducedDipolePolar[jIndex][1] + particleI.quadrupole[QXZ]*_inducedDipolePolar[jIndex][2];
    qiukp2            = particleI.quadrupole[QXY]*_inducedDipolePolar[jIndex][0] + particleI.quadrupole[QYY]*_inducedDipolePolar[jIndex][1] + particleI.quadrupole[QYZ]*_inducedDipolePolar[jIndex][2];
    qiukp3            = particleI.quadrupole[QXZ]*_inducedDipolePolar[jIndex][0] + particleI.quadrupole[QYZ]*_inducedDipolePolar[jIndex][1] + particleI.quadrupole[QZZ]*_inducedDipolePolar[jIndex][2];

    qkuip1            = particleJ.quadrupole[QXX]*_inducedDipolePolar[iIndex][0] + particleJ.quadrupole[QXY]*_inducedDipolePolar[iIndex][1] + particleJ.quadrupole[QXZ]*_inducedDipolePolar[iIndex][2];
    qkuip2            = particleJ.quadrupole[QXY]*_inducedDipolePolar[iIndex][0] + particleJ.quadrupole[QYY]*_inducedDipolePolar[iIndex][1] + particleJ.quadrupole[QYZ]*_inducedDipolePolar[iIndex][2];
    qkuip3            = particleJ.quadrupole[QXZ]*_inducedDipolePolar[iIndex][0] + particleJ.quadrupole[QYZ]*_inducedDipolePolar[iIndex][1] + particleJ.quadrupole[QZZ]*_inducedDipolePolar[iIndex][2];

    uixqkr1           = _inducedDipole[iIndex][1]*qkr3 - _inducedDipole[iIndex][2]*qkr2;
    uixqkr2           = _inducedDipole[iIndex][2]*qkr1 - _inducedDipole[iIndex][0]*qkr3;
    uixqkr3           = _inducedDipole[iIndex][0]*qkr2 - _inducedDipole[iIndex][1]*qkr1;

    ukxqir1           = _inducedDipole[jIndex][1]*qir3 - _inducedDipole[jIndex][2]*qir2;
    ukxqir2           = _inducedDipole[jIndex][2]*qir1 - _inducedDipole[jIndex][0]*qir3;
    ukxqir3           = _inducedDipole[jIndex][0]*qir2 - _inducedDipole[jIndex][1]*qir1;

    uixqkrp1          = _inducedDipolePolar[iIndex][1]*qkr3 - _inducedDipolePolar[iIndex][2]*qkr2;
    uixqkrp2          = _inducedDipolePolar[iIndex][2]*qkr1 - _inducedDipolePolar[iIndex][0]*qkr3;
    uixqkrp3          = _inducedDipolePolar[iIndex][0]*qkr2 - _inducedDipolePolar[iIndex][1]*qkr1;

    ukxqirp1          = _inducedDipolePolar[jIndex][1]*qir3 - _inducedDipolePolar[jIndex][2]*qir2;
    ukxqirp2          = _inducedDipolePolar[jIndex][2]*qir1 - _inducedDipolePolar[jIndex][0]*qir3;
    ukxqirp3          = _inducedDipolePolar[jIndex][0]*qir2 - _inducedDipolePolar[jIndex][1]*qir1;

    rxqiuk1           = yr*qiuk3 - zr*qiuk2;
    rxqiuk2           = zr*qiuk1 - xr*qiuk3;
    rxqiuk3           = xr*qiuk2 - yr*qiuk1;

    rxqkui1           = yr*qkui3 - zr*qkui2;
    rxqkui2           = zr*qkui1 - xr*qkui3;
    rxqkui3           = xr*qkui2 - yr*qkui1;

    rxqiukp1          = yr*qiukp3 - zr*qiukp2;
    rxqiukp2          = zr*qiukp1 - xr*qiukp3;
    rxqiukp3          = xr*qiukp2 - yr*qiukp1;

    rxqkuip1          = yr*qkuip3 - zr*qkuip2;
    rxqkuip2          = zr*qkuip1 - xr*qkuip3;
    rxqkuip3          = xr*qkuip2 - yr*qkuip1;

    // get intermediate variables for induction energy terms

    sci1              = _inducedDipole[iIndex][0]*particleJ.dipole[0] + _inducedDipole[iIndex][1]*particleJ.dipole[1]
                          + _inducedDipole[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipole[jIndex][0]
                          + particleI.dipole[1]*_inducedDipole[jIndex][1] + particleI.dipole[2]*_inducedDipole[jIndex][2];

    sci3              = _inducedDipole[iIndex][0]*xr + _inducedDipole[iIndex][1]*yr + _inducedDipole[iIndex][2]*zr;
    sci4              = _inducedDipole[jIndex][0]*xr + _inducedDipole[jIndex][1]*yr + _inducedDipole[jIndex][2]*zr;

    sci7              = qir1*_inducedDipole[jIndex][0] + qir2*_inducedDipole[jIndex][1] + qir3*_inducedDipole[jIndex][2];
    sci8              = qkr1*_inducedDipole[iIndex][0] + qkr2*_inducedDipole[iIndex][1] + qkr3*_inducedDipole[iIndex][2];

    scip1             = _inducedDipolePolar[iIndex][0]*particleJ.dipole[0] + _inducedDipolePolar[iIndex][1]*particleJ.dipole[1] + _inducedDipolePolar[iIndex][2]*particleJ.dipole[2] + particleI.dipole[0]*_inducedDipolePolar[jIndex][0] + particleI.dipole[1]*_inducedDipolePolar[jIndex][1] + particleI.dipole[2]*_inducedDipolePolar[jIndex][2];
    scip2             = _inducedDipole[iIndex][0]*_inducedDipolePolar[jIndex][0]+_inducedDipole[iIndex][1]*_inducedDipolePolar[jIndex][1] + _inducedDipole[iIndex][2]*_inducedDipolePolar[jIndex][2]+_inducedDipolePolar[iIndex][0]*_inducedDipole[jIndex][0] + _inducedDipolePolar[iIndex][1]*_inducedDipole[jIndex][1]+_inducedDipolePolar[iIndex][2]*_inducedDipole[jIndex][2];

    scip3             = _inducedDipolePolar[iIndex][0]*xr + _inducedDipolePolar[iIndex][1]*yr + _inducedDipolePolar[iIndex][2]*zr;
    scip4             = _inducedDipolePolar[jIndex][0]*xr + _inducedDipolePolar[jIndex][1]*yr + _inducedDipolePolar[jIndex][2]*zr;

    scip7             = qir1*_inducedDipolePolar[jIndex][0] + qir2*_inducedDipolePolar[jIndex][1] + qir3*_inducedDipolePolar[jIndex][2];
    scip8             = qkr1*_inducedDipolePolar[iIndex][0] + qkr2*_inducedDipolePolar[iIndex][1] + qkr3*_inducedDipolePolar[iIndex][2];

    // calculate the gl functions for potential energy

    gli1              = particleJ.charge*sci3 - particleI.charge*sci4;
    gli2              = -sc3*sci4 - sci3*sc4;
    gli3              = sci3*sc6 - sci4*sc5;
    gli6              = sci1;
    gli7              = two*(sci7-sci8);

    glip1             = particleJ.charge*scip3 - particleI.charge*scip4;
    glip2             = -sc3*scip4 - scip3*sc4;
    glip3             = scip3*sc6 - scip4*sc5;
    glip6             = scip1;
    glip7             = two*(scip7-scip8);

    // get the permanent multipole and induced energies

    energy           += -oneHalf*(rr3*(gli1+gli6)*psc3 + rr5*(gli2+gli7)*psc5 + rr7*gli3*psc7);

    // intermediate variables for the induced-permanent terms;

    gfi1             = oneHalf*rr5*((gli1+gli6)*psc3 + (glip1+glip6)*dsc3+scip2*scale3i)
                         + oneHalf*rr7*((gli7+gli2)*psc5 +(glip7+glip2)*dsc5
                         -(sci3*scip4+scip3*sci4)*scale5i)
                         + oneHalf*rr9*(gli3*psc7+glip3*dsc7);

    gfi4             = two*rr5;
    gfi5             = rr7*(sci4*psc7+scip4*dsc7);
    gfi6             = -rr7*(sci3*psc7+scip3*dsc7);

    // get the induced force

    ftm2i1           = gfi1*xr + oneHalf*
                        (- rr3*particleJ.charge*(_inducedDipole[iIndex][0]*psc3+_inducedDipolePolar[iIndex][0]*dsc3)
                         + rr5*sc4*(_inducedDipole[iIndex][0]*psc5+_inducedDipolePolar[iIndex][0]*dsc5)
                         - rr7*sc6*(_inducedDipole[iIndex][0]*psc7+_inducedDipolePolar[iIndex][0]*dsc7))
                         +(rr3*particleI.charge*(_inducedDipole[jIndex][0]*psc3+_inducedDipolePolar[jIndex][0]*dsc3)
                         + rr5*sc3*(_inducedDipole[jIndex][0]*psc5+_inducedDipolePolar[jIndex][0]*dsc5)
                         + rr7*sc5*(_inducedDipole[jIndex][0]*psc7+_inducedDipolePolar[jIndex][0]*dsc7))*oneHalf
                         + rr5*scale5i*(sci4*_inducedDipolePolar[iIndex][0]+scip4*_inducedDipole[iIndex][0]
                         + sci3*_inducedDipolePolar[jIndex][0]+scip3*_inducedDipole[jIndex][0])*oneHalf
                         + oneHalf*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[0]
                         + oneHalf*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[0]
                         + oneHalf*gfi4*((qkui1-qiuk1)*psc5
                         + (qkuip1-qiukp1)*dsc5)
                         + gfi5*qir1 + gfi6*qkr1;

    ftm2i2           = gfi1*yr + oneHalf*
                        (- rr3*particleJ.charge*(_inducedDipole[iIndex][1]*psc3+_inducedDipolePolar[iIndex][1]*dsc3)
                         + rr5*sc4*(_inducedDipole[iIndex][1]*psc5+_inducedDipolePolar[iIndex][1]*dsc5)
                         - rr7*sc6*(_inducedDipole[iIndex][1]*psc7+_inducedDipolePolar[iIndex][1]*dsc7))
                         +(rr3*particleI.charge*(_inducedDipole[jIndex][1]*psc3+_inducedDipolePolar[jIndex][1]*dsc3)
                         + rr5*sc3*(_inducedDipole[jIndex][1]*psc5+_inducedDipolePolar[jIndex][1]*dsc5)
                         + rr7*sc5*(_inducedDipole[jIndex][1]*psc7+_inducedDipolePolar[jIndex][1]*dsc7))*oneHalf
                         + rr5*scale5i*(sci4*_inducedDipolePolar[iIndex][1]+scip4*_inducedDipole[iIndex][1]
                         + sci3*_inducedDipolePolar[jIndex][1]+scip3*_inducedDipole[jIndex][1])*oneHalf
                         + oneHalf*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[1]
                         + oneHalf*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[1]
                         + oneHalf*gfi4*((qkui2-qiuk2)*psc5
                         + (qkuip2-qiukp2)*dsc5)
                         + gfi5*qir2 + gfi6*qkr2;

    ftm2i3           = gfi1*zr  + oneHalf*
                         (- rr3*particleJ.charge*(_inducedDipole[iIndex][2]*psc3+_inducedDipolePolar[iIndex][2]*dsc3)
                         + rr5*sc4*(_inducedDipole[iIndex][2]*psc5+_inducedDipolePolar[iIndex][2]*dsc5)
                         - rr7*sc6*(_inducedDipole[iIndex][2]*psc7+_inducedDipolePolar[iIndex][2]*dsc7))
                         +(rr3*particleI.charge*(_inducedDipole[jIndex][2]*psc3+_inducedDipolePolar[jIndex][2]*dsc3)
                         + rr5*sc3*(_inducedDipole[jIndex][2]*psc5+_inducedDipolePolar[jIndex][2]*dsc5)
                         + rr7*sc5*(_inducedDipole[jIndex][2]*psc7+_inducedDipolePolar[jIndex][2]*dsc7))*oneHalf
                         + rr5*scale5i*(sci4*_inducedDipolePolar[iIndex][2]+scip4*_inducedDipole[iIndex][2]
                         + sci3*_inducedDipolePolar[jIndex][2]+scip3*_inducedDipole[jIndex][2])*oneHalf
                         + oneHalf*(sci4*psc5+scip4*dsc5)*rr5*particleI.dipole[2]
                         + oneHalf*(sci3*psc5+scip3*dsc5)*rr5*particleJ.dipole[2]
                         + oneHalf*gfi4*((qkui3-qiuk3)*psc5
                         + (qkuip3-qiukp3)*dsc5)
                         + gfi5*qir3 + gfi6*qkr3;

    // intermediate values needed for partially excluded interactions

    fridmp1          = oneHalf*(rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_1
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_1
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_1);

    fridmp2          = oneHalf*(rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_2
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_2
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_2);

    fridmp3          = oneHalf*(rr3*((gli1+gli6)*pscale
                         +(glip1+glip6)*dscale)*ddsc3_3
                         + rr5*((gli2+gli7)*pscale
                         +(glip2+glip7)*dscale)*ddsc5_3
                         + rr7*(gli3*pscale+glip3*dscale)*ddsc7_3);

    // get the induced-induced derivative terms;

    findmp1          = oneHalf*uscale*(scip2*rr3*ddsc3_1
                         - rr5*ddsc5_1*(sci3*scip4+scip3*sci4));

    findmp2          = oneHalf*uscale*(scip2*rr3*ddsc3_2
                         - rr5*ddsc5_2*(sci3*scip4+scip3*sci4));

    findmp3          = oneHalf*uscale*(scip2*rr3*ddsc3_3
                         - rr5*ddsc5_3*(sci3*scip4+scip3*sci4));

    // handle of scaling for partially excluded interactions;

    ftm2i1           = ftm2i1 - fridmp1 - findmp1;
    ftm2i2           = ftm2i2 - fridmp2 - findmp2;
    ftm2i3           = ftm2i3 - fridmp3 - findmp3;

    // correction to convert mutual to direct polarization force;

    if( polarizationType == AmoebaReferenceMultipoleForce::Direct ){

        RealOpenMM gfd    = oneHalf*(rr5*scip2*scale3i- rr7*(scip3*sci4+sci3*scip4)*scale5i);
        RealOpenMM fdir1  = gfd*xr + oneHalf*rr5*scale5i* (sci4*_inducedDipolePolar[iIndex][0]+scip4*_inducedDipole[iIndex][0] + sci3*_inducedDipolePolar[jIndex][0]+scip3*_inducedDipole[jIndex][0]);
        RealOpenMM fdir2  = gfd*yr + oneHalf*rr5*scale5i* (sci4*_inducedDipolePolar[iIndex][1]+scip4*_inducedDipole[iIndex][1] + sci3*_inducedDipolePolar[jIndex][1]+scip3*_inducedDipole[jIndex][1]);
        RealOpenMM fdir3  = gfd*zr + oneHalf*rr5*scale5i* (sci4*_inducedDipolePolar[iIndex][2]+scip4*_inducedDipole[iIndex][2] + sci3*_inducedDipolePolar[jIndex][2]+scip3*_inducedDipole[jIndex][2]);
        ftm2i1      -= fdir1 - findmp1;
        ftm2i2      -= fdir2 - findmp2;
        ftm2i3      -= fdir3 - findmp3;
    }

    // now perform the torque calculation
    // intermediate terms for torque between multipoles i and k

    gti2             = oneHalf*(sci4*psc5+scip4*dsc5)*rr5;
    gti3             = oneHalf*(sci3*psc5+scip3*dsc5)*rr5;
    gti4             = gfi4;
    gti5             = gfi5;
    gti6             = gfi6;

    // calculate the induced torque components

    torqueI[0]      -= -rr3*(dixuk1*psc3+dixukp1*dsc3)*oneHalf
                         + gti2*dixr1 + gti4*((ukxqir1+rxqiuk1)*psc5
                         +(ukxqirp1+rxqiukp1)*dsc5)*oneHalf - gti5*rxqir1;

    torqueI[1]      -= -rr3*(dixuk2*psc3+dixukp2*dsc3)*oneHalf
                         + gti2*dixr2 + gti4*((ukxqir2+rxqiuk2)*psc5
                         +(ukxqirp2+rxqiukp2)*dsc5)*oneHalf - gti5*rxqir2;

    torqueI[2]      -= -rr3*(dixuk3*psc3+dixukp3*dsc3)*oneHalf
                         + gti2*dixr3 + gti4*((ukxqir3+rxqiuk3)*psc5
                         +(ukxqirp3+rxqiukp3)*dsc5)*oneHalf - gti5*rxqir3;

    torqueJ[0]      -= -rr3*(dkxui1*psc3+dkxuip1*dsc3)*oneHalf
                         + gti3*dkxr1 - gti4*((uixqkr1+rxqkui1)*psc5
                         +(uixqkrp1+rxqkuip1)*dsc5)*oneHalf - gti6*rxqkr1;

    torqueJ[1]      -= -rr3*(dkxui2*psc3+dkxuip2*dsc3)*oneHalf
                         + gti3*dkxr2 - gti4*((uixqkr2+rxqkui2)*psc5
                         +(uixqkrp2+rxqkuip2)*dsc5)*oneHalf - gti6*rxqkr2;

    torqueJ[2]      -= -rr3*(dkxui3*psc3+dkxuip3*dsc3)*oneHalf
                         + gti3*dkxr3 - gti4*((uixqkr3+rxqkui3)*psc5
                         +(uixqkrp3+rxqkuip3)*dsc5)*oneHalf - gti6*rxqkr3;

    // update force and torque on site k;

    force[0]       += ftm2i1;
    force[1]       += ftm2i2;
    force[2]       += ftm2i3;

    force          *= (_electric/_dielectric);
    forces[iIndex] += force;
    forces[jIndex] -= force;

    torqueI         *= (_electric/_dielectric);
    torqueJ         *= (_electric/_dielectric);
    torques[iIndex] += torqueI;
    torques[jIndex] += torqueJ;

    return energy;
}
