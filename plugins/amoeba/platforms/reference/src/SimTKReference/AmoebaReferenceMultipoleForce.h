
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

#include "SimTKUtilities/RealVec.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "AmoebaReferenceGeneralizedKirkwoodForce.h"
#include <map>

typedef std::map< unsigned int, RealOpenMM> MapIntRealOpenMM;
typedef MapIntRealOpenMM::iterator MapIntRealOpenMMI;
typedef MapIntRealOpenMM::const_iterator MapIntRealOpenMMCI;

using namespace OpenMM;

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

    enum PolarizationType {

        /** 
         * Mutual polarization
         */
        Mutual = 0,

        /** 
         * Direct polarization
         */
        Direct = 1 
    };  

    /**
     * Constructor
     * 
     */
    AmoebaReferenceMultipoleForce( );
 
    /**
     * Constructor
     * 
     * @param nonbondedMethod nonbonded method
     */
    AmoebaReferenceMultipoleForce( NonbondedMethod nonbondedMethod );
 
    /**
     * Destructor
     * 
     */
    virtual ~AmoebaReferenceMultipoleForce( ){};
 
    /**
     * Get nonbonded method
     * 
     * @return nonbonded method
     */
    NonbondedMethod getNonbondedMethod( void ) const;

    /**
     * Set nonbonded method
     * 
     * @param nonbondedMethod nonbonded method
     */
    void setNonbondedMethod( NonbondedMethod nonbondedMethod );

    /**
     * Get flag indicating if mutual induced dipoles are converged
     *
     * @return nonzero if converged
     *
     */
    int getMutualInducedDipoleConverged( void ) const;

    /**
     * Get number of iterations used in computing mutual induced dipoles
     *
     * @return number of iterations
     *
     */
    int getMutualInducedDipoleIterations( void ) const;

    /**
     * Get the final epsilon for mutual induced dipoles
     *
     *  @return epsilon
     *
     */
    RealOpenMM getMutualInducedDipoleEpsilon( void ) const;

    /**
     * Set the target epsilon for converging mutual induced dipoles
     *
     * @param targetEpsilon target epsilon for converging mutual induced dipoles
     *
     */
    void setMutualInducedDipoleTargetEpsilon( RealOpenMM targetEpsilon );

    /**
     * Get the target epsilon for converging mutual induced dipoles
     *
     * @return target epsilon for converging mutual induced dipoles
     *
     */
    RealOpenMM getMutualInducedDipoleTargetEpsilon( void ) const;

    /**
     * Set the maximum number of iterations to be executed in converging mutual induced dipoles
     *
     * @param maximumMutualInducedDipoleIterations maximum number of iterations to be executed in converging mutual induced dipoles
     *
     */
    void setMaximumMutualInducedDipoleIterations( int maximumMutualInducedDipoleIterations );

    /**
     * Get the maximum number of iterations to be executed in converging mutual induced dipoles
     *
     * @return maximum number of iterations to be executed in converging mutual induced dipoles
     * 
     */
    int getMaximumMutualInducedDipoleIterations( void ) const;

    /**
     * Calculate force and energy       
     *
     * @param numParticles            number of particles
     * @param particlePositions       Cartesian coordinates of particles
     * @param forces                  add forces to this vector
     *
     * @return energy
     */
    RealOpenMM calculateForceAndEnergy( const std::vector<OpenMM::RealVec>& particlePositions,
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
                                        std::vector<OpenMM::RealVec>& forces );

protected:

    enum MultipoleParticleDataEnum { PARTICLE_POSITION, PARTICLE_CHARGE, PARTICLE_DIPOLE, PARTICLE_QUADRUPOLE,
                                     PARTICLE_THOLE, PARTICLE_DAMPING_FACTOR, PARTICLE_POLARITY, PARTICLE_FIELD, 
                                     PARTICLE_FIELD_POLAR, GK_FIELD, PARTICLE_INDUCED_DIPOLE, PARTICLE_INDUCED_DIPOLE_POLAR };

    enum QuadrupoleIndices { QXX, QXY, QXZ, QYY, QYZ, QZZ };

    /* 
     * Particle parameters and coordinates
     */
    class MultipoleParticleData {
        public:
            unsigned int particleIndex;    
            RealVec position;
            RealOpenMM charge;
            RealVec dipole;
            RealOpenMM quadrupole[6];
            RealOpenMM thole;
            RealOpenMM dampingFactor;
            RealOpenMM polarity;
    };

    /* 
     * Helper class used in calculating induced dipoles
     */
    class UpdateInducedDipoleField {
        public:
            UpdateInducedDipoleField( std::vector<OpenMM::RealVec>* inputFixed_E_Field, std::vector<OpenMM::RealVec>* inputInducedDipoles );
            std::vector<OpenMM::RealVec>* fixed_E_Field;
            std::vector<OpenMM::RealVec>* inducedDipoles;
            std::vector<OpenMM::RealVec> inducedDipoleField;
    };

    unsigned int _numParticles;

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

    std::vector<RealVec> _fixed_E_Field;
    std::vector<RealVec> _fixed_E_FieldPolar;
    std::vector<RealVec> _inducedDipole;
    std::vector<RealVec> _inducedDipolePolar;

    int _mutualInducedDipoleConverged;
    int _mutualInducedDipoleIterations;
    int _maximumMutualInducedDipoleIterations;
    RealOpenMM  _mutualInducedDipoleEpsilon;
    RealOpenMM  _mutualInducedDipoleTargetEpsilon;
    RealOpenMM  _polarSOR;
    RealOpenMM  _debye;

    /**
     * Helper constructor method to centralize initialization of objects
     *
     */
    void initialize( void );

    /**
     * Load particle data
     *
     * @param particlePositions   particle coordinates
     * @param charges             charges
     * @param dipoles             dipoles
     * @param quadrupoles         quadrupoles
     * @param tholes              Thole parameters
     * @param dampingFactors      dampming factors
     * @param polarity            polarity
     * @param particleData        output data struct
     *
     */
    void loadParticleData( const std::vector<OpenMM::RealVec>& particlePositions, 
                           const std::vector<RealOpenMM>& charges,
                           const std::vector<RealOpenMM>& dipoles,
                           const std::vector<RealOpenMM>& quadrupoles,
                           const std::vector<RealOpenMM>& tholes,
                           const std::vector<RealOpenMM>& dampingFactors,
                           const std::vector<RealOpenMM>& polarity,
                           std::vector<MultipoleParticleData>& particleData ) const;

    /**
     * Calculate fixed electric fields
     *
     * @param particleData vector particle data
     * 
     */
    virtual void calculateFixedEField( vector<MultipoleParticleData>& particleData );

    /**
     * Set flag indicating if mutual induced dipoles are converged
     * 
     * @param converged nonzero if converged
     *
     */
    void setMutualInducedDipoleConverged( int converged );

    /**
     * Set number of iterations used in computing mutual induced dipoles
     * 
     * @param  number of iterations
     * 
     */
    void setMutualInducedDipoleIterations( int iterations );

    /**
     * Set the final epsilon for mutual induced dipoles
     * 
     * @param epsilon
     *
     */
    void setMutualInducedDipoleEpsilon( RealOpenMM epsilon );

    /**
     * Setup scale factors given covalent info
     *
     * @param  multipoleAtomCovalentInfo vector of vectors containing the covalent info
     *
     */
    void setupScaleMaps( const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo );

    /**
     * Get multipole scale factor for particleI & particleJ
     * 
     * @param  particleI           index of particleI whose scale factor is to be retrieved
     * @param  particleJ           index of particleJ whose scale factor is to be retrieved
     * @param  scaleType           scale type (D_SCALE, P_SCALE, M_SCAL)
     *
     * @return scaleFactor 
     */
    RealOpenMM getMultipoleScaleFactor( unsigned int particleI, unsigned int particleJ, ScaleType scaleType ) const;

    /**
     * Get scale factor for particleI & particleJ
     * 
     * @param  particleI           index of particleI whose scale factor is to be retrieved
     * @param  particleJ           index of particleJ whose scale factor is to be retrieved
     * @param  scaleType           scale type (D_SCALE, P_SCALE, M_SCAL)
     *
     * @return array of scaleFactors 
     */
    void getMultipoleScaleFactors( unsigned int particleI, unsigned int particleJ, std::vector<RealOpenMM>& scaleFactors ) const;

    /**
     * Get p- and d-scale factors for particleI & particleJ ixn
     *
     * @param  particleI           index of particleI 
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */
    void getDScaleAndPScale( unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale ) const;
    
    /**
     * Calculate damped powers of 1/r
     *
     * @param  particleI           index of particleI 
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */
    void getAndScaleInverseRs( RealOpenMM dampI, RealOpenMM dampJ, RealOpenMM tholeI, RealOpenMM tholeJ,
                               RealOpenMM r, std::vector<RealOpenMM>& rrI ) const;

    /**
     * Check if multipoles at chiral site should be inverted
     *
     * @param  particleI            particleI data 
     * @param  axisType             axis type
     * @param  particleZ            z-axis particle to particleI
     * @param  particleX            x-axis particle to particleI
     * @param  particleY            y-axis particle to particleI
     *
     */
    void checkChiralCenterAtParticle( MultipoleParticleData& particleI, int axisType, MultipoleParticleData& particleZ,
                                      MultipoleParticleData& particleX, MultipoleParticleData& particleY ) const;

    /**
     * Invert multipole moments (dipole[Y], quadrupole[XY] and quadrupole[YZ]) if chiral center inverted 
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param multipoleAtomXs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomZs         vector of y-particle indices used to map molecular frame to lab frame
     * @param axisType                axis type
     */
    void checkChiral( std::vector<MultipoleParticleData>& particleData, 
                              const std::vector<int>& multipoleAtomXs,
                              const std::vector<int>& multipoleAtomYs,
                              const std::vector<int>& multipoleAtomZs,
                              const std::vector<int>& axisTypes ) const;
    /**
     * Apply rotation matrix to molecular dipole/quadrupoles to get corresponding lab frame values
     * for particle I
     * 
     * @param  particleI            particleI data
     * @param  particleJ            particleI data
     * @param  axisType             axis type
     */
    void applyRotationMatrixToParticle(       MultipoleParticleData& particleI,
                                        const MultipoleParticleData& particleZ,
                                        const MultipoleParticleData& particleX,
                                              MultipoleParticleData* particleY, int axisType ) const;

    /**
     * Apply rotation matrix to molecular dipole/quadrupoles to get corresponding lab frame values
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     *                                dipole and quadrupole entries are modified
     * @param multipoleAtomXs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomZs         vector of y-particle indices used to map molecular frame to lab frame
     * @param axisType                axis type
     */
    void applyRotationMatrix( std::vector<MultipoleParticleData>& particleData, 
                              const std::vector<int>& multipoleAtomXs,
                              const std::vector<int>& multipoleAtomYs,
                              const std::vector<int>& multipoleAtomZs,
                              const std::vector<int>& axisTypes ) const;
    /**
     * Zero fixed E-fields
     */
    virtual void zeroFixed_E_Fields( void );

    /**
     * Calculate electric field at particle I due fixed multipoles at particle J and vice versa
     * (field at particle J due fixed multipoles at particle I)
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    virtual void calculateFixedEFieldPairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                              RealOpenMM dScale, RealOpenMM pScale );

    /**
     * Calculate field at particle I due induced dipole at particle J and vice versa
     * (field at particle J due induced dipole at particle I)
     * 
     * @param particleI               index of particle I
     * @param particleJ               index of particle J
     * @param rr3                     damped 1/r^3 factor
     * @param rr5                     damped 1/r^5 factor
     * @param delta                   delta of particle positions: particleJ.x - particleI.x, ...
     * @param inducedDipole           vector of induced dipoles
     * @param field                   vector of induced dipole fields
     */
    void calculateInducedDipolePairIxn( unsigned int particleI, unsigned int particleJ,
                                        RealOpenMM rr3, RealOpenMM rr5, const RealVec& delta,
                                        const std::vector<RealVec>& inducedDipole,
                                        std::vector<RealVec>& field ) const;

    /**
     * Calculate fields due induced dipoles at each site
     * 
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleField containing input induced dipoles and output fields
     */
    virtual void calculateInducedDipolePairIxns( const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                 std::vector<UpdateInducedDipoleField>& updateInducedDipoleFields );

    /**
     * Converge induced dipoles
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleField containing input induced dipoles and output fields
     */
    void convergeInduceDipoles( const std::vector<MultipoleParticleData>& particleData,
                                std::vector<UpdateInducedDipoleField>& calculateInducedDipoleField );

    /**
     * Update fields due to induced dipoles for each particle
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleField containing input induced dipoles and output fields
     */
    RealOpenMM updateInducedDipoleFields( const std::vector<MultipoleParticleData>& particleData,
                                          std::vector<UpdateInducedDipoleField>& calculateInducedDipoleField);

    /**
     * Update induced dipole for a particle given updated induced dipole field at the site
     * 
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param fixed_E_Field             fields due fixed multipoles at each site
     * @param inducedDipoleField        fields due induced dipoles at each site
     * @param inducedDipoles            output vector of updated induced dipoles
     */
    RealOpenMM updateInducedDipole( const std::vector<MultipoleParticleData>& particleI,
                                    const std::vector<RealVec>& fixed_E_Field,
                                    const std::vector<RealVec>& inducedDipoleField,
                                    std::vector<RealVec>& inducedDipoles);

    /**
     * Calculate induced dipoles
     * 
     * @param polarizationType  if 'Direct' polariztion, only initial induced dipoles calculated
     *                          if 'Mutual' polariztion, induced dipoles converged to specified tolerance
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    virtual void calculateInducedDipoles( AmoebaReferenceMultipoleForce::PolarizationType polarizationType, const std::vector<MultipoleParticleData>& particleData );

    /**
     * Calculate electrostatic interaction between particles I and K
     * 
     * @param polarizationType  if 'Direct' polariztion, only initial induced dipoles calculated
     *                          if 'Mutual' polariztion, induced dipoles converged to specified tolerance
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param scalingFactors    scaling factors for interaction
     * @param forces            vector of particle forces to be updated
     * @param torques           vector of particle torques to be updated
     */
    RealOpenMM calculateElectrostaticPairIxn( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                              const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                              const std::vector<RealOpenMM>& scalingFactors, std::vector<OpenMM::RealVec>& forces, std::vector<RealVec>& torque ) const;

    /**
     * Map torque to forces
     * 
     * @param particleI               particle whose torque is to be mapped
     * @param particleU               particle1 of lab frame for particleI 
     * @param particleV               particle2 of lab frame for particleI 
     * @param particleW               particle3 of lab frame for particleI 
     * @param axisType                axis type (Bisector/Z-then-X, ...)
     * @param torque                  torque on particle I
     */
    void mapTorqueToForceForParticle( const MultipoleParticleData& particleI,
                                      const MultipoleParticleData& particleU,
                                      const MultipoleParticleData& particleV,
                                            MultipoleParticleData* particleW,
                                      int axisType, const Vec3& torque, std::vector<OpenMM::RealVec>& forces ) const;

    /**
     * Map torques to forces
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param multipoleAtomZs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomXs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of y-particle indices used to map molecular frame to lab frame
     * @param axisType                vector of axis types (Bisector/Z-then-X, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    void mapTorqueToForce( std::vector<MultipoleParticleData>& particleData, 
                           const std::vector<int>& multipoleAtomXs,
                           const std::vector<int>& multipoleAtomYs,
                           const std::vector<int>& multipoleAtomZs,
                           const std::vector<int>& axisTypes,
                           std::vector<OpenMM::RealVec>& torques,
                           std::vector<OpenMM::RealVec>& forces ) const;

    /**
     * Calculate electrostatic forces
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param axisType                vector of axis types (Bisector/Z-then-X, ...) for particles
     * @param multipoleAtomZs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomXs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of y-particle indices used to map molecular frame to lab frame
     * @param polarizationType        polarization type (direct or mutual)
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    virtual RealOpenMM calculateElectrostatic( std::vector<MultipoleParticleData>& particleData, 
                                               const std::vector<int>& axisTypes,
                                               const std::vector<int>& multipoleAtomZs,
                                               const std::vector<int>& multipoleAtomXs,
                                               const std::vector<int>& multipoleAtomYs,
                                               AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                               std::vector<OpenMM::RealVec>& torques,
                                               std::vector<OpenMM::RealVec>& forces ) const;

    /**
     * Normalize a RealVec
     *
     * @param vectorToNormalize vector to normalize
     *
     * @return norm of vector on input
     * 
     */
    RealOpenMM normalizeRealVec( RealVec& vectorToNormalize ) const;

    /**
     * Initialize vector of RealOpenMM (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeRealOpenMMVector( vector<RealOpenMM>& vectorToInitialize ) const;

    /**
     * Initialize vector of RealVec (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeRealVecVector( vector<RealVec>& vectorToInitialize ) const;

    /**
     * Copy vector of RealVec
     *
     * @param inputVector  vector to copy
     * @param outputVector output vector
     * 
     */
    void copyRealVecVector( const std::vector<OpenMM::RealVec>& inputVector, std::vector<OpenMM::RealVec>& outputVector ) const;

 void showScaleMapForParticle( unsigned int particleI, FILE* log ) const;
};

class AmoebaReferenceGeneralizedKirkwoodMultipoleForce : public AmoebaReferenceMultipoleForce {

public:

    /**
     * Constructor
     * 
     */
    AmoebaReferenceGeneralizedKirkwoodMultipoleForce( AmoebaReferenceGeneralizedKirkwoodForce* amoebaReferenceGeneralizedKirkwoodForce );
 
    /**
     * Destructor
     * 
     */
    ~AmoebaReferenceGeneralizedKirkwoodMultipoleForce( );
 
    /**
     * Get flag signalling whether cavity term is to be included
     *
     * @return flag
     *
     */
    int getIncludeCavityTerm( void ) const;

    /**
     * Get probe radius
     *
     * @return probe radius
     *
     */
    RealOpenMM getProbeRadius( void ) const;

    /**
     * Get surface area factor
     *
     * @return surface area factor
     *
     */
    RealOpenMM getSurfaceAreaFactor( void ) const;

    /**
     * Get dielectric offset
     *
     * @return dielectric offset
     *
     */
    RealOpenMM getDielectricOffset( void ) const;

private:

    AmoebaReferenceGeneralizedKirkwoodForce* _amoebaReferenceGeneralizedKirkwoodForce;

    RealOpenMM _gkc;
    RealOpenMM _fc;
    RealOpenMM _fd;
    RealOpenMM _fq;

    std::vector<RealOpenMM> _atomicRadii;
    std::vector<RealOpenMM> _scaledRadii;
    std::vector<RealOpenMM> _bornRadii;
    std::vector<RealOpenMM> _bornForce;

    std::vector<RealVec> _gkField;
    std::vector<RealVec> _inducedDipoleS;
    std::vector<RealVec> _inducedDipolePolarS;

    int _includeCavityTerm;
    RealOpenMM _probeRadius;
    RealOpenMM _surfaceAreaFactor;
    RealOpenMM _dielectricOffset;

    /**
     * Zero fixed E-fields
     *
     */
    void zeroFixed_E_Fields( void );

    /**
     * Calculate electric field at particle I due fixed multipoles at particle J and vice versa
     * (field at particle J due fixed multipoles at particle I)
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    void calculateFixedEFieldPairIxn( MultipoleParticleData& particleI, MultipoleParticleData& particleJ,
                                      RealOpenMM dScale, RealOpenMM pScale );

    /**
     * Calculate induced dipoles
     * 
     * @param polarizationType  if 'Direct' polariztion, only initial induced dipoles calculated
     *                          if 'Mutual' polariztion, induced dipoles converged to specified tolerance
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void calculateInducedDipoles( AmoebaReferenceMultipoleForce::PolarizationType polarizationType, const std::vector<MultipoleParticleData>& particleData );

    /**
     * Calculate fields due induced dipoles at each site
     * 
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleField containing input induced dipoles and output fields
     */
    void calculateInducedDipolePairIxns( const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                         std::vector<UpdateInducedDipoleField>& updateInducedDipoleFields );

    /**
     * Calculate electrostatic forces
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param axisType                vector of axis types (Bisector/Z-then-X, ...) for particles
     * @param multipoleAtomZs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomXs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of y-particle indices used to map molecular frame to lab frame
     * @param polarizationType        polarization type (direct or mutual)
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    RealOpenMM calculateElectrostatic( std::vector<MultipoleParticleData>& particleData, 
                                       const std::vector<int>& axisTypes,
                                       const std::vector<int>& multipoleAtomZs,
                                       const std::vector<int>& multipoleAtomXs,
                                       const std::vector<int>& multipoleAtomYs,
                                       AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                       std::vector<OpenMM::RealVec>& torques,
                                       std::vector<OpenMM::RealVec>& forces ) const;


    /**
     * Calculate GK field at particle I due induced dipole at particle J and vice versa
     * (field at particle J due induced dipole at particle I)
     * 
     * @param particleI               index of particle I
     * @param particleJ               index of particle J
     * @param field                   vector of induced dipole fields
     * @param fieldPolar              vector of induced dipole polar fields
     */
    void calculateInducedDipolePairGkIxn( const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                          const std::vector<RealVec>& field, std::vector<RealVec>& fieldPolar ) const;

    /**
     * Calculate Kirkwood interaction
     * 
     * @param polarizationType        polarization type ( direct or mutual )
     * @param particleI               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param forces                  add Kirkwood force to forces
     * @param torques                 add Kirkwood torque to torques
     * @param dBorn                   chain-rule factor
     *
     * @return energy
     */
    RealOpenMM calculateKirkwoodPairIxn( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                         const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                         std::vector<RealVec>& forces, 
                                         std::vector<RealVec>& torques,
                                         std::vector<RealOpenMM>& dBorn ) const;

    /**
     * Calculate Grycuk 'chain-rule' force
     * 
     * @param particleI               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dBorn                   chain-rule Born force factor
     * @param forces                  add Kirkwood force to forces
     *
     */
    void calculateGrycukChainRulePairIxn( const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                          const std::vector<RealOpenMM>& dBorn, std::vector<RealVec>& forces ) const;

    /**
     * Calculate TINKER's ACE approximation to non-polar cavity term 
     * 
     * @param forces                  add Kirkwood force to forces
     *
     * @return energy
     *
     */
    RealOpenMM calculateCavityTermEnergyAndForces( std::vector<RealOpenMM>& dBorn ) const;

    /**
     * Correct vacuum to SCRF derivatives
     * 
     * @param polarizationType        polarization type ( direct or mutual )
     * @param particleI               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param pscale                  p-scale factor
     * @param dscale                  d-scale factor
     * @param forces                  force accumulator
     * @param torques                 torque accumulator
     *
     * @return energy
     */
    RealOpenMM calculateKirkwoodEDiffPairIxn( AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                              const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                              RealOpenMM pscale, RealOpenMM dscale, 
                                              std::vector<RealVec>& forces, std::vector<RealVec>& torques ) const;

};

#endif // _AmoebaReferenceMultipoleForce___
