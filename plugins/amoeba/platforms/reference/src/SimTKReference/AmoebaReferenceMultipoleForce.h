
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

#ifndef __AmoebaReferenceMultipoleForce_H__
#define __AmoebaReferenceMultipoleForce_H__

#include "SimTKUtilities/SimTKOpenMMRealType.h"
#include "openmm/Vec3.h"
#include "AmoebaMultipoleForce.h"
#include <string>
#include <vector>
#include <map>

typedef std::map< unsigned int, RealOpenMM> MapIntRealOpenMM;
typedef MapIntRealOpenMM::iterator MapIntRealOpenMMI;
typedef MapIntRealOpenMM::const_iterator MapIntRealOpenMMCI;

typedef std::vector<std::vector<RealOpenMM> > VectorOfRealOpenMMVectors;
typedef VectorOfRealOpenMMVectors::iterator VectorOfRealOpenMMVectorsI;
typedef VectorOfRealOpenMMVectors::const_iterator VectorOfRealOpenMMVectorsCI;

using namespace OpenMM;

// ---------------------------------------------------------------------------------------

class AmoebaReferenceMultipoleForce {

public:

    /** 
     * This is an enumeration of the different methods that may be used for handling long range Multipole forces.
     */
    enum NonbondedMethod {

        /**
         * No cutoff is applied to the interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */

        NoCutoff = 0,

        /**
         * Interactions beyond the cutoff distance are ignored.  
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.  
         */
        CutoffPeriodic = 2
    };
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceMultipoleForce( );
 
    /**---------------------------------------------------------------------------------------
       
       Constructor
       
       --------------------------------------------------------------------------------------- */
 
    AmoebaReferenceMultipoleForce( NonbondedMethod nonbondedMethod );
 
    /**---------------------------------------------------------------------------------------
       
       Destructor
       
       --------------------------------------------------------------------------------------- */
 
    ~AmoebaReferenceMultipoleForce( ){};
 
    /**---------------------------------------------------------------------------------------
    
       Get nonbonded method
    
       @return nonbonded method
    
       --------------------------------------------------------------------------------------- */
    
    NonbondedMethod getNonbondedMethod( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Set nonbonded method
    
       @param nonbonded method
    
       --------------------------------------------------------------------------------------- */
    
    void setNonbondedMethod( NonbondedMethod nonbondedMethod );

    /**---------------------------------------------------------------------------------------
    
       Get flag indicating if mutual induced dipoles are converged
    
       @return nonzero if converged
    
       --------------------------------------------------------------------------------------- */
    
    int getMutualInducedDipoleConverged( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Get number of iterations used in computing mutual induced dipoles
    
       @return number of iterations
    
       --------------------------------------------------------------------------------------- */
    
    int getMutualInducedDipoleIterations( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Get the final epsilon for mutual induced dipoles
    
       @return epsilon
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getMutualInducedDipoleEpsilon( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Set the target epsilon for converging mutual induced dipoles
    
       @param  target epsilon for converging mutual induced dipoles
    
       --------------------------------------------------------------------------------------- */
    
    void setMutualInducedDipoleTargetEpsilon( RealOpenMM targetEpsilon );

    /**---------------------------------------------------------------------------------------
    
       Get the target epsilon for converging mutual induced dipoles
    
       @return target epsilon for converging mutual induced dipoles
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getMutualInducedDipoleTargetEpsilon( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Set the maximum number of iterations to be executed in converging mutual induced dipoles
    
       @param  maximum number of iterations to be executed in converging mutual induced dipoles
    
       --------------------------------------------------------------------------------------- */
    
    void setMaximumMutualInducedDipoleIterations( int maximumMutualInducedDipoleIterations );

    /**---------------------------------------------------------------------------------------
    
       Get the maximum number of iterations to be executed in converging mutual induced dipoles
    
       @return maximum number of iterations to be executed in converging mutual induced dipoles
    
       --------------------------------------------------------------------------------------- */
    
    int getMaximumMutualInducedDipoleIterations( void ) const;

    /**---------------------------------------------------------------------------------------
    
       Calculate Amoeba Hal vdw ixns
    
       @param numParticles            number of particles
       @param particlePositions       Cartesian coordinates of particles
       @param indexIVs                position index for associated reducing particle
       @param indexClasses            class index for combining sigmas/epsilons (not currently used)
       @param sigmas                  particle sigmas 
       @param epsilons                particle epsilons
       @param reductions              particle reduction factors
       @param vdwExclusions           particle exclusions
       @param forces                  add forces to this vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculateForceAndEnergy( int numParticles, RealOpenMM** particlePositions,
                                        const std::vector<RealOpenMM>& charges,
                                        const std::vector<RealOpenMM>& dipoles,
                                        const std::vector<RealOpenMM>& quadrupoles,
                                        const std::vector<RealOpenMM>& tholes,
                                        const std::vector<RealOpenMM>& dampingFactors,
                                        const std::vector<RealOpenMM>& polarity,
                                        const std::vector<int>& axisTypes,
                                        const std::vector<int>& multipoleAtomId1s,
                                        const std::vector<int>& multipoleAtomId2s,
                                        const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo,
                                        RealOpenMM** forces );


         
private:

    struct MultipoleParticleData {
        unsigned int particleIndex;    
        RealOpenMM position[3];
        RealOpenMM charge;
        RealOpenMM dipole[3];
        RealOpenMM quadrupole[6];
        RealOpenMM thole;
        RealOpenMM dampingFactor;
        RealOpenMM polarity;
        RealOpenMM field[3];
        RealOpenMM fieldPolar[3];
        RealOpenMM inducedDipole[3];
        RealOpenMM inducedDipolePolar[3];
    };
    enum MultipoleParticleDataEnum { PARTICLE_POSITION, PARTICLE_CHARGE, PARTICLE_DIPOLE, PARTICLE_QUADRUPOLE,
                                     PARTICLE_THOLE, PARTICLE_DAMPING_FACTOR, PARTICLE_POLARITY, PARTICLE_FIELD, 
                                     PARTICLE_FIELD_POLAR, PARTICLE_INDUCED_DIPOLE, PARTICLE_INDUCED_DIPOLE_POLAR };

    NonbondedMethod _nonbondedMethod;

    RealOpenMM _electric;
    RealOpenMM _dielectric;

    enum ScaleType { D_SCALE, P_SCALE, M_SCALE, U_SCALE, LAST_SCALE_TYPE_INDEX };
    std::vector<  std::vector< MapIntRealOpenMM > > _scaleMaps;
    std::vector<unsigned int> _maxScaleIndex;
    RealOpenMM _dScale[5];
    RealOpenMM _pScale[5];
    RealOpenMM _mScale[5];
    RealOpenMM _uScale[5];

    int _mutualInducedDipoleConverged;
    int _mutualInducedDipoleIterations;
    int _maximumMutualInducedDipoleIterations;
    RealOpenMM  _mutualInducedDipoleEpsilon;
    RealOpenMM  _mutualInducedDipoleTargetEpsilon;
    RealOpenMM  _polarSOR;
    RealOpenMM  _debye;

    enum QuadrupoleIndices { QXX, QXY, QXZ, QYY, QYZ, QZZ };

    /**---------------------------------------------------------------------------------------
    
       Helper constructor method to centralize initialization of objects
    
       --------------------------------------------------------------------------------------- */
    
    void initialize( void );

    void loadParticleData( RealOpenMM** particlePositions, 
                                                      const std::vector<RealOpenMM>& charges,
                                                      const std::vector<RealOpenMM>& dipoles,
                                                      const std::vector<RealOpenMM>& quadrupoles,
                                                      const std::vector<RealOpenMM>& tholes,
                                                      const std::vector<RealOpenMM>& dampingFactors,
                                                      const std::vector<RealOpenMM>& polarity,
                                                      std::vector<MultipoleParticleData>& particleData ) const;

    /**---------------------------------------------------------------------------------------
    
       Set flag indicating if mutual induced dipoles are converged
    
       @param nonzero if converged
    
       --------------------------------------------------------------------------------------- */
    
    void setMutualInducedDipoleConverged( int iterations );

    /**---------------------------------------------------------------------------------------
    
       Set number of iterations used in computing mutual induced dipoles
    
       @param  number of iterations
    
       --------------------------------------------------------------------------------------- */
    
    void setMutualInducedDipoleIterations( int iterations );

    /**---------------------------------------------------------------------------------------
    
       Set the final epsilon for mutual induced dipoles
    
       @param epsilon
    
       --------------------------------------------------------------------------------------- */
    
    void setMutualInducedDipoleEpsilon( RealOpenMM epsilon );

    /**---------------------------------------------------------------------------------------
    
       Get delta between positions of two particles
    
       @param  particleI            index of particleI
       @param  particleJ            index of particleJ
       @param  delta                output: delta[0] = pos[particleJ][0] - pos[particleI][0], ...
    
       --------------------------------------------------------------------------------------- */
    
    void getDelta( unsigned int particleI, unsigned int particleJ, RealOpenMM** particlePositions, RealOpenMM* delta ) const;

    /**---------------------------------------------------------------------------------------
    
       Get delta between positions of two particles
    
       @param  particleI            index of particleI
       @param  particleJ            index of particleJ
       @param  delta                output: delta[0] = pos[particleJ][0] - pos[particleI][0], ...
    
       --------------------------------------------------------------------------------------- */
    
    void getDelta( const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, RealOpenMM* delta ) const;

    /**---------------------------------------------------------------------------------------
    
       Load values from vector (dipole, quadrupoles vectors for example) into an array
    
       @param  particleI            index of particleI whose values are to be loaded
       @param  offset               offset to use (3 for dipole, 9 for quadrupole)
       @param  vectorToCopy         vector of values to be loaded
       @param  loadArray            on return contains loaded values

       --------------------------------------------------------------------------------------- */
    
    void loadArrayFromVector( unsigned int particleI, unsigned int offset, const std::vector<RealOpenMM>& vectorToCopy, RealOpenMM* loadArray ) const;

    /**---------------------------------------------------------------------------------------
    
       Add fields valus to vector
    
       @param  particleIOffset     index of particleI times offset (3,9 for dipole. quadrupole) whose values are to be updated
       @param  sign                if > 0, then add the field values; otherwise subtract the values
       @param  field               field values
       @param  vectorToAddTo       vector whose of values are to be incremented by the values in field

       --------------------------------------------------------------------------------------- */
    
    void addField( unsigned int particleIOffset, int sign, RealOpenMM field[3], std::vector<RealOpenMM>& vectorToAddTo ) const;

    void logRealOpenMMVectors( const std::string& header, const VectorOfRealOpenMMVectors& printVector,
                               FILE* log = stderr, unsigned int itemsPerVector = 1, int maxPrint = -1 ) const;

    void logParticleData( const std::string& header, const std::vector<MultipoleParticleData>& particleData,
                          unsigned int printFlag, FILE* log, unsigned int maxPrint ) const;

    void showScaleMapForParticle( unsigned int particleI, FILE* log ) const;

    /**---------------------------------------------------------------------------------------
    
       Setup scale factors given covalent info
    
       @param  multipoleAtomCovalentInfo vector of vectors containing the covalent info

       --------------------------------------------------------------------------------------- */
    
    void setupScaleMaps( const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo );

    /**---------------------------------------------------------------------------------------
    
       Get scale factor for particleI & particleJ
    
       @param  particleI           index of particleI whose scale factor is to be retreived
       @param  particleJ           index of particleJ whose scale factor is to be retreived
       @param  scaleType           scale type (D_SCALE, P_SCALE, M_SCAL)

       @return scaleFactor 

       --------------------------------------------------------------------------------------- */
    
    RealOpenMM getScaleFactor(  unsigned int particleI, unsigned int particleJ, ScaleType scaleType      ) const;
    void       getScaleFactors( unsigned int particleI, unsigned int particleJ, RealOpenMM* scaleFactors ) const;

    /**---------------------------------------------------------------------------------------
    
       Scale field for particleI & particleJ ixn
    
       @param  particleI           index of particleI 
       @param  particleJ           index of particleJ
       @param  field               field to scale

       --------------------------------------------------------------------------------------- */
    
    void getDScaleAndPScale( unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale ) const;
    
    void getAndScaleInverseRs( RealOpenMM dampI, RealOpenMM dampJ, RealOpenMM tholeI, RealOpenMM tholeJ,
                               RealOpenMM r, std::vector<RealOpenMM>& rrI ) const;

    /**---------------------------------------------------------------------------------------
    
       Apply roatation matrix to molecular dipole/quadrupoles to get corresponding lab frame values
       for particle I
    
       @param  particleI            particleI data
       @param  particleJ            particleI data
       @param  axisType             axis type
    
       --------------------------------------------------------------------------------------- */
    
    void applyRotationMatrix(       MultipoleParticleData& particleI,
                              const MultipoleParticleData& particleZ,
                              const MultipoleParticleData& particleX, int axisType ) const;

    void calculateFixedEFieldPairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                      RealOpenMM dScale, RealOpenMM pScale ) const;

    void calculateInducedDipolePairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                        std::vector<RealOpenMM>& field, std::vector<RealOpenMM>& fieldPolar ) const;

    void updateInducedDipole( MultipoleParticleData& particleI,
                              const std::vector<RealOpenMM>& field,
                              const std::vector<RealOpenMM>& fieldPolar,
                              RealOpenMM& epsilonD, RealOpenMM& epsilonP ) const;

    void calculateNoCutoffInducedDipoles( std::vector<MultipoleParticleData>& particleData );

    void calculateInducedDipoleField( std::vector<MultipoleParticleData>& particleData,
                                      std::vector<RealOpenMM>& field,
                                      std::vector<RealOpenMM>& fieldPolar ) const;

    RealOpenMM calculateNoCutoffElectrostaticPairIxn( const MultipoleParticleData& particleI,
                                                      const MultipoleParticleData& particleK,
                                                      RealOpenMM* scalingFactors, RealOpenMM** forces, std::vector<Vec3>& torque ) const;

    RealOpenMM calculateNoCutoffElectrostatic( std::vector<MultipoleParticleData>& particleData, RealOpenMM** forces ) const;

    /**---------------------------------------------------------------------------------------
    
       Calculate Multipole ixns w/ no cutoff
    
       @param numParticles            number of particles
       @param particlePositions       Cartesian coordinates of particles
       @param indexIVs                position index for associated reducing particle
       @param indexClasses            class index for combining sigmas/epsilons (not currently used)
       @param sigmas                  particle sigmas 
       @param epsilons                particle epsilons
       @param reductions              particle reduction factors
       @param vdwExclusions           particle exclusions
       @param forces                  add forces to this vector
    
       @return energy
    
       --------------------------------------------------------------------------------------- */
    
    RealOpenMM calculateNoCutoffForceAndEnergy( unsigned int numParticles, RealOpenMM** particlePositions,
                                                const std::vector<RealOpenMM>& charges,
                                                const std::vector<RealOpenMM>& dipoles,
                                                const std::vector<RealOpenMM>& quadrupoles,
                                                const std::vector<RealOpenMM>& tholes,
                                                const std::vector<RealOpenMM>& dampingFactors,
                                                const std::vector<RealOpenMM>& polarity,
                                                const std::vector<int>& axisTypes,
                                                const std::vector<int>& multipoleAtomId1s,
                                                const std::vector<int>& multipoleAtomId2s,
                                                RealOpenMM** forces );

};

// ---------------------------------------------------------------------------------------

#endif // _AmoebaReferenceMultipoleForce___
