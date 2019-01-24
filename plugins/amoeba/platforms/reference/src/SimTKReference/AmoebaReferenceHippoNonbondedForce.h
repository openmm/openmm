/* Portions copyright (c) 2006-2018 Stanford University and Simbios.
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

#ifndef __AmoebaReferenceHippoNonbondedForce_H__
#define __AmoebaReferenceHippoNonbondedForce_H__

#include "openmm/HippoNonbondedForce.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <map>
#include <utility>
#include <vector>
#include "fftpack.h"
#include <complex>

namespace OpenMM {

/**
 * 3-dimensional int vector
 */
class HippoIntVec {
public:
    /**
     * Create a HippoIntVec whose elements are all 0.
     */
    HippoIntVec() {
        data[0] = data[1] = data[2] = 0;
    }
    /**
     * Create a HippoIntVec with specified x, y, z, w components.
     */
    HippoIntVec(int x, int y, int z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    int operator[](int index) const {
        assert(index >= 0 && index < 3);
        return data[index];
    }
    int& operator[](int index) {
        assert(index >= 0 && index < 3);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    HippoIntVec operator+() const {
        return HippoIntVec(*this);
    }

    // plus
    HippoIntVec operator+(const HippoIntVec& rhs) const {
        const HippoIntVec& lhs = *this;
        return HippoIntVec(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }

    HippoIntVec& operator+=(const HippoIntVec& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    HippoIntVec& operator-=(const HippoIntVec& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

private:
    int data[3];
};

/**
 * 4-dimensional double vector
 */
class HippoDouble4 {
public:
    /**
     * Create a HippoDouble4 whose elements are all 0.
     */
    HippoDouble4() {
        data[0] = data[1] = data[2] = data[3] = 0.0;
    }
    /**
     * Create a HippoDouble4 with specified x, y, z, w components.
     */
    HippoDouble4(double x, double y, double z, double w) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }
    double operator[](int index) const {
        assert(index >= 0 && index < 4);
        return data[index];
    }
    double& operator[](int index) {
        assert(index >= 0 && index < 4);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    HippoDouble4 operator+() const {
        return HippoDouble4(*this);
    }

    // plus
    HippoDouble4 operator+(const HippoDouble4& rhs) const {
        const HippoDouble4& lhs = *this;
        return HippoDouble4(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],lhs[3] + rhs[3]);
    }

    HippoDouble4& operator+=(const HippoDouble4& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        data[3] += rhs[3];
        return *this;
    }

    HippoDouble4& operator-=(const HippoDouble4& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        data[3] -= rhs[3];
        return *this;
    }

private:
    double data[4];
};

using namespace OpenMM;

class AmoebaReferenceHippoNonbondedForce {

   /**
    * AmoebaReferenceHippoNonbondedForce is base class for MultipoleForce calculations
    * AmoebaReferencePmeHippoNonbondedForce is derived class for PME calculations
    *
    * Below is a outline of the sequence of methods called to evaluate the force and energy 
    * for each scenario: Generalized Kirkwood (GK) and PME.
    *
    * If 'virtual' appears before the method name, the method is overridden in one or more of the derived classes.
    *
    * calculateForceAndEnergy()                            calculate forces  and energy
    *
    *    setup()                                           rotate molecular multipole moments to lab frame
    *                                                      setup scaling maps and calculate induced dipoles (see calculateInducedDipoles below)
    *
    *    virtual calculateElectrostatic()                  calculate forces and torques
    *
    *                                                      GK case includes the following calls:
    *
    *                                                          AmoebaReferenceHippoNonbondedForce::calculateElectrostatic()
    *                                                               loop over particle pairs: calculateElectrostaticPairIxn()
    *
    *                                                          TINKER's egk1a: calculateKirkwoodPairIxn()
    *
    *                                                          SASA force and energy: calculateCavityTermEnergyAndForces()
    *
    *                                                          TINKER's born1 (Born chain rule term): loop over particle pairs: calculateGrycukChainRulePairIxn()
    *
    *                                                          TINKER's ediff1: loop over particle pairs: calculateKirkwoodEDiffPairIxn()
    *                                                       
    *                                                      PME case includes the following calls:
    *
    *                                                          reciprocal [computeReciprocalSpaceInducedDipoleForceAndEnergy(),
    *                                                                      computeReciprocalSpaceFixedMultipoleForceAndEnergy]
    *
    *                                                          direct space calculations [calculatePmeDirectElectrostaticPairIxn()]
    *
    *                                                          self-energy [calculatePmeSelfEnergy()]
    *
    *                                                          torques [calculatePmeSelfTorque()]
    *
    *    mapTorqueToForce()                                map torques to forces
    * 
    * setup()
    *    loadParticleData()                                load particle data (polarity, multipole moments, Thole factors, ...)
    *    checkChiral()                                     if needed, invert multipole moments at chiral centers
    *    applyRotationMatrix()                             rotate molecular multipole moments to lab frame
    *    setupScaleMaps()                                  setup scaling maps
    *    calculateInducedDipoles()                         calculate induced dipoles
    * 
    * 
    * virtual calculateInducedDipoles()                    calculate induced dipoles:
    *                                                          field at each site due to fixed multipoles first calculated
    *                                                          if polarization type == Direct,
    *                                                          initial induced dipoles are calculated, but are not converged.
    *                                                          if polarization type == Mutual, then loop until
    *                                                          induce dipoles converge.
    *                                                       For GK, include gkField in setup
    *                                                       For PME, base class method is used
    *
    * 
    *     virtual zeroFixedMultipoleFields()                zero fixed multipole vectors; for GK includes zeroing of gkField vector
    * 
    *     virtual calculateFixedMultipoleField()            calculate fixed multipole field -- particle pair loop 
    *                                                       gkField also calculated for GK
    *                                                       for PME, reciprocal, direct space (particle pair loop) and self terms calculated
    *                                                       
    * 
    *         virtual calculateFixedMultipoleFieldPairIxn() pair ixn for fixed multipole
    *                                                       gkField also calculated for GK
    *                                                       for PME, direct space ixn calculated here
    * 
    *     virtual initializeInducedDipoles()                initialize induced dipoles; for PME, calculateReciprocalSpaceInducedDipoleField()
    *                                                       called in case polarization type == Direct
    *
    *     convergeInduceDipoles()                           loop until induced dipoles converge
    *
    *         updateInducedDipoleFields()                   update fields at each site due other induced dipoles
    *
    *           virtual calculateInducedDipoleFields()      calculate induced dipole field at each site by looping over particle pairs 
    *                                                       for PME includes reciprocal space calculation calculateReciprocalSpaceInducedDipoleField(), 
    *                                                       direct space calculateDirectInducedDipolePairIxns() and self terms
    *
    *              virtual calculateInducedDipolePairIxns() field at particle i due particle j's induced dipole and vice versa; for GK includes GK field
    */

public:

    /**
     * Constructor
     * 
     */
    AmoebaReferenceHippoNonbondedForce(const HippoNonbondedForce& force);

    /**
     * Destructor
     * 
     */
    virtual ~AmoebaReferenceHippoNonbondedForce() {};
 
    /**
     * Get nonbonded method.
     * 
     * @return nonbonded method
     */
    HippoNonbondedForce::NonbondedMethod getNonbondedMethod() const;

    /**
     * Set the coefficients for the µ_0, µ_1, µ_2, µ_n terms in the extrapolation
     * theory algorithm for induced dipoles
     *
     * @param optCoefficients a vector whose mth entry specifies the coefficient for µ_m
     *
     */
    void setExtrapolationCoefficients(const std::vector<double> &coefficients);

    /**
     * Calculate force and energy.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param forces                    add forces to this vector
     *
     * @return energy
     */
    double calculateForceAndEnergy(const std::vector<OpenMM::Vec3>& particlePositions,
                                   std::vector<OpenMM::Vec3>& forces);

    /**
     * Calculate particle induced dipoles.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param outputMultipoleMoments    output multipole moments
     */
    void calculateInducedDipoles(const std::vector<OpenMM::Vec3>& particlePositions,
                                 std::vector<Vec3>& outputInducedDipoles);

    /**
     * Calculate particle permanent dipoles rotated in the lab frame.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param outputMultipoleMoments    output multipole moments
     */

    void calculateLabFramePermanentDipoles(const std::vector<Vec3>& particlePositions,
                                           std::vector<Vec3>& outputRotatedPermanentDipoles);

    /**
     * Calculate particle total dipoles.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param outputMultipoleMoments    output multipole moments
     */


    void calculateTotalDipoles(const std::vector<Vec3>& particlePositions,
                               std::vector<Vec3>& outputRotatedPermanentDipoles);

    /**
     * Calculate electrostatic potential at a set of grid points.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param input grid                input grid points to compute potential
     * @param outputPotential           output electrostatic potential
     */
    void calculateElectrostaticPotential(const std::vector<OpenMM::Vec3>& particlePositions,
                                         const std::vector<Vec3>& inputGrid,
                                         std::vector<double>& outputPotential);

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
            int particleIndex, axisType, multipoleAtomX, multipoleAtomY, multipoleAtomZ;    
            Vec3 position;
            Vec3 dipole;
            double quadrupole[6];
            Vec3 sphericalDipole;
            double sphericalQuadrupole[5];
            double coreCharge, valenceCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
    };
    
    /**
     * Particle parameters transformed into fractional coordinates
     */
    class TransformedMultipole {
    public:
        double charge;
        Vec3 dipole;
        double quadrupole[6];
    };
    
    /**
     * Information defining an exception.
     */
    class Exception {
    public:
        int particle1, particle2;
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale;
    };

    /* 
     * Helper class used in calculating induced dipoles
     */
    struct UpdateInducedDipoleFieldStruct {
            UpdateInducedDipoleFieldStruct(std::vector<OpenMM::Vec3>& inputFixed_E_Field, std::vector<OpenMM::Vec3>& inputInducedDipoles, std::vector<std::vector<Vec3> >& extrapolatedDipoles, std::vector<std::vector<double> >& extrapolatedDipoleFieldGradient);
            std::vector<OpenMM::Vec3>* fixedMultipoleField;
            std::vector<OpenMM::Vec3>* inducedDipoles;
            std::vector<std::vector<Vec3> >* extrapolatedDipoles;
            std::vector<std::vector<double> >* extrapolatedDipoleFieldGradient;
            std::vector<OpenMM::Vec3> inducedDipoleField;
            std::vector<std::vector<double> > inducedDipoleFieldGradient;
    };

    unsigned int _numParticles;

    HippoNonbondedForce::NonbondedMethod _nonbondedMethod;

    double _electric;

    enum ScaleType { D_SCALE, P_SCALE, M_SCALE, U_SCALE, LAST_SCALE_TYPE_INDEX };
    std::map<std::pair<int, int>, Exception> exceptions;

    std::vector<MultipoleParticleData> particleData;
    std::vector<TransformedMultipole> _transformed;
    std::vector<Vec3> _fixedMultipoleField;
    std::vector<Vec3> _fixedMultipoleFieldPolar;
    std::vector<Vec3> _inducedDipole;
    std::vector<Vec3> _inducedDipolePolar;
    std::vector<std::vector<Vec3> > _ptDipoleP;
    std::vector<std::vector<Vec3> > _ptDipoleD;
    std::vector<std::vector<double> > _ptDipoleFieldGradientP;
    std::vector<std::vector<double> > _ptDipoleFieldGradientD;

    int _maxPTOrder;
    std::vector<double>  _extrapolationCoefficients;
    std::vector<double>  _extPartCoefficients;

    /**
     * Load particle data.
     *
     * @param particlePositions   particle coordinates
     *
     */
    void loadParticleData(const std::vector<OpenMM::Vec3>& particlePositions);

    /**
     * Calculate fixed multipole fields.
     *
     * @param particleData vector of particle data
     * 
     */
    virtual void calculateFixedMultipoleField();


    /**
     * Get scale factor for particleI & particleJ
     * 
     * @param  particleI           index of particleI whose scale factor is to be retrieved
     * @param  particleJ           index of particleJ whose scale factor is to be retrieved
     * @param  scaleType           scale type (D_SCALE, P_SCALE, M_SCAL)
     *
     * @return array of scaleFactors 
     */
    void getMultipoleScaleFactors(unsigned int particleI, unsigned int particleJ, std::vector<double>& scaleFactors) const;

    /**
     * Get p- and d-scale factors for particleI & particleJ ixn
     *
     * @param  particleI           index of particleI 
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */
    void getDScaleAndPScale(unsigned int particleI, unsigned int particleJ, double& dScale, double& pScale) const;
    
    /**
     * Calculate damped powers of 1/r.
     *
     * @param  particleI           index of particleI 
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */
    void getAndScaleInverseRs(double dampI, double dampJ, double tholeI, double tholeJ,
                              double r, std::vector<double>& rrI) const;
    /**
     * Compute the damping factors for the field due to the fixed multipole moments of a particle.
     * 
     * @param particle     parameters for the particle
     * @param r            the distance from the particle at which to calculate the damping factors
     * @param fdamp3       outputs the damping factor for the r^-3 term
     * @param fdamp5       outputs the damping factor for the r^-5 term
     * @param fdamp7       outputs the damping factor for the r^-7 term
     */
    void computeDirectFieldDampingFactors(const MultipoleParticleData& particle, double r, double& fdamp3, double& fdamp5, double& fdamp7);
    /**
     * Compute the damping factors for the field at one particle due to the induced dipole of another particle.
     * 
     * @param particleI    parameters for the first particle
     * @param particleJ    parameters for the second particle
     * @param r            the distance between the two particles
     * @param fdamp3       outputs the damping factor for the r^-3 term
     * @param fdamp5       outputs the damping factor for the r^-5 term
     */
    void computeMutualFieldDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r, double& fdamp3, double& fdamp5);
    /**
     * Compute the damping factors for the overlap of two particles.
     * 
     * @param particleI    parameters for the first particle
     * @param particleJ    parameters for the second particle
     * @param r            the distance between the two particles
     * @param fdampI1      outputs the damping factor for the r^-1 term for the valence charge of particle I
     * @param fdampI3      outputs the damping factor for the r^-3 term for the valence charge of particle I
     * @param fdampI5      outputs the damping factor for the r^-5 term for the valence charge of particle I
     * @param fdampI7      outputs the damping factor for the r^-7 term for the valence charge of particle I
     * @param fdampI9      outputs the damping factor for the r^-9 term for the valence charge of particle I
     * @param fdampJ1      outputs the damping factor for the r^-1 term for the valence charge of particle J
     * @param fdampJ3      outputs the damping factor for the r^-3 term for the valence charge of particle J
     * @param fdampJ5      outputs the damping factor for the r^-5 term for the valence charge of particle J
     * @param fdampJ7      outputs the damping factor for the r^-7 term for the valence charge of particle J
     * @param fdampJ9      outputs the damping factor for the r^-9 term for the valence charge of particle J
     * @param fdampIJ1     outputs the damping factor for the r^-1 term for the overlap between valence charges
     * @param fdampIJ3     outputs the damping factor for the r^-3 term for the overlap between valence charges
     * @param fdampIJ5     outputs the damping factor for the r^-5 term for the overlap between valence charges
     * @param fdampIJ7     outputs the damping factor for the r^-7 term for the overlap between valence charges
     * @param fdampIJ9     outputs the damping factor for the r^-9 term for the overlap between valence charges
     * @param fdampIJ11    outputs the damping factor for the r^-9 term for the overlap between valence charges
     */
    void computeOverlapDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r,
            double& fdampI1, double& fdampI3, double& fdampI5, double& fdampI7, double& fdampI9,
            double& fdampJ1, double& fdampJ3, double& fdampJ5, double& fdampJ7, double& fdampJ9,
            double& fdampIJ1, double& fdampIJ3, double& fdampIJ5, double& fdampIJ7, double& fdampIJ9, double& fdampIJ11);
    /**
     * Check if multipoles at chiral site should be inverted.
     *
     * @param  particleI            particleI data 
     * @param  axisType             axis type
     * @param  particleZ            z-axis particle to particleI
     * @param  particleX            x-axis particle to particleI
     * @param  particleY            y-axis particle to particleI
     *
     */
    void checkChiralCenterAtParticle(MultipoleParticleData& particleI, int axisType, MultipoleParticleData& particleZ,
                                     MultipoleParticleData& particleX, MultipoleParticleData& particleY);

    /**
     * Invert multipole moments (dipole[Y], quadrupole[XY] and quadrupole[YZ]) if chiral center inverted.
     */
    void checkChiral();
    /**
     * Apply rotation matrix to molecular dipole/quadrupoles to get corresponding lab frame values
     * for particle I.
     * 
     * @param  particleI            particleI data
     * @param  particleJ            particleI data
     * @param  axisType             axis type
     */
    void applyRotationMatrixToParticle(      MultipoleParticleData& particleI,
                                       const MultipoleParticleData* particleZ,
                                       const MultipoleParticleData* particleX,
                                             MultipoleParticleData* particleY, int axisType) const;

    /**
     * Forms the rotation matrix for the quasi-internal coordinate system,
     * which is the rotation matrix that describes the orientation of the
     * internuclear vector for a given pair (I,J) in lab frame.
     *
     * @param particleI             particleI position
     * @param particleJ             particleJ position
     * @param deltaR                the internuclear vector, corrected for periodic boundary conditions
     * @param r                     the bond length between atoms I and J
     * @param rotationmatrix        the output rotation matrix for a 3-vector
     */
    void formQIRotationMatrix(const Vec3& iPosition,
                              const Vec3& jPosition,
                              const Vec3 &deltaR,
                              double r,
                              double (&rotationMatrix)[3][3]) const;


    /**
     * Constructs a rotation matrix for spherical harmonic quadrupoles, using the dipole rotation matrix.
     *
     * @param D1                    The input spherical harmonic dipole rotation matrix
     * @param D2                    The output spherical harmonic quadrupole rotation matrix
     */
     void buildSphericalQuadrupoleRotationMatrix(const double (&D1)[3][3], double (&D2)[5][5]) const;

     /**
      * Constructs a rotation matrix for spherical harmonic quadrupoles, using the dipole rotation matrix.
      * Only the m={0,1c,1s} terms are constructed; these are the only terms needed to evaluate the field.
      *
      * @param D1                    The input spherical harmonic dipole rotation matrix
      * @param D2                    The output spherical harmonic quadrupole rotation matrix
      */
      void buildPartialSphericalQuadrupoleRotationMatrix(const double (&D1)[3][3], double (&D2)[3][5]) const;

    /**
     * Apply rotation matrix to molecular dipole/quadrupoles to get corresponding lab frame values.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     *                                dipole and quadrupole entries are modified
     * @param multipoleAtomXs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomZs         vector of y-particle indices used to map molecular frame to lab frame
     * @param axisType                axis type
     */
    void applyRotationMatrix();
    /**
     * Zero fixed multipole fields.
     */
    virtual void zeroFixedMultipoleFields();

    /**
     * Calculate electric field at particle I due fixed multipoles at particle J and vice versa
     * (field at particle J due fixed multipoles at particle I).
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    virtual void calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                     double dScale, double pScale);

    /**
     * Initialize induced dipoles
     *
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    virtual void initializeInducedDipoles(std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields); 

    /**
     * Calculate field at particle I due induced dipole at particle J and vice versa
     * (field at particle J due induced dipole at particle I).
     * 
     * @param particleI               index of particle I
     * @param particleJ               index of particle J
     * @param rr3                     damped 1/r^3 factor
     * @param rr5                     damped 1/r^5 factor
     * @param delta                   delta of particle positions: particleJ.x - particleI.x, ...
     * @param inducedDipole           vector of induced dipoles
     * @param field                   vector of induced dipole fields
     */
    void calculateInducedDipolePairIxn(unsigned int particleI, unsigned int particleJ,
                                       double rr3, double rr5, const Vec3& delta,
                                       const std::vector<Vec3>& inducedDipole,
                                       std::vector<Vec3>& field) const;

    /**
     * Calculate fields due induced dipoles at each site.
     *
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    virtual void calculateInducedDipolePairIxns(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields);

    /**
     * Calculate induced dipole fields.
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    virtual void calculateInducedDipoleFields(const std::vector<MultipoleParticleData>& particleData,
                                              std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields);
    /**
     * Calculated induced dipoles using extrapolated perturbation theory.
     *
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void convergeInduceDipolesByExtrapolation(const std::vector<MultipoleParticleData>& particleData,
                                              std::vector<UpdateInducedDipoleFieldStruct>& calculateInducedDipoleField);

    /**
     * Calculate induced dipoles.
     */
    virtual void calculateInducedDipoles();

    /**
     * Setup: 
     *        if needed invert multipole moments at chiral centers
     *        rotate molecular multipole moments to lab frame 
     *        setup scaling maps and 
     *        calculate induced dipoles (see calculateInducedDipoles below)
     *
     * @param particlePositions         Cartesian coordinates of particles
     */
    void setup(const std::vector<OpenMM::Vec3>& particlePositions);

    /**
     * Calculate electrostatic interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param scalingFactors    scaling factors for interaction
     * @param forces            vector of particle forces to be updated
     * @param torque            vector of particle torques to be updated
     */
    double calculateElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                         const std::vector<double>& scalingFactors, std::vector<OpenMM::Vec3>& forces, std::vector<Vec3>& torque) const;

    /**
     * Map particle torque to force.
     * 
     * @param particleI               particle whose torque is to be mapped
     * @param particleU               particle1 of lab frame for particleI 
     * @param particleV               particle2 of lab frame for particleI 
     * @param particleW               particle3 of lab frame for particleI 
     * @param axisType                axis type (Bisector/Z-then-X, ...)
     * @param torque                  torque on particle I
     * @param forces                  vector of particle forces to be updated
     */
    void mapTorqueToForceForParticle(const MultipoleParticleData& particleI,
                                     const MultipoleParticleData& particleU,
                                     const MultipoleParticleData& particleV,
                                           MultipoleParticleData* particleW,
                                     int axisType, const Vec3& torque, std::vector<OpenMM::Vec3>& forces);

    /**
     * Map torques to forces.
     * 
     * @param torques                 output torques
     * @param forces                  output forces 
     */
    void mapTorqueToForce(std::vector<OpenMM::Vec3>& torques,
                          std::vector<OpenMM::Vec3>& forces);

    /**
     * Calculate electrostatic forces
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    virtual double calculateElectrostatic(std::vector<OpenMM::Vec3>& torques,
                                          std::vector<OpenMM::Vec3>& forces);

    /**
     * Normalize a Vec3
     *
     * @param vectorToNormalize vector to normalize
     *
     * @return norm of vector on input
     * 
     */
    double normalizeVec3(Vec3& vectorToNormalize) const;

    /**
     * Initialize vector of double (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeRealOpenMMVector(std::vector<double>& vectorToInitialize) const;

    /**
     * Initialize vector of Vec3 (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeVec3Vector(std::vector<Vec3>& vectorToInitialize) const;

    /**
     * Copy vector of Vec3
     *
     * @param inputVector  vector to copy
     * @param outputVector output vector
     * 
     */
    void copyVec3Vector(const std::vector<OpenMM::Vec3>& inputVector, std::vector<OpenMM::Vec3>& outputVector) const;

    /**
     * Calculate potential at grid point due to a particle
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param gridPoint               grid point
     *
     * @return potential at grid point
     * 
     */
    double calculateElectrostaticPotentialForParticleGridPoint(const MultipoleParticleData& particleI, const Vec3& gridPoint) const;

    /**
     * Apply periodic boundary conditions to difference in positions
     * 
     * @param deltaR  difference in particle positions; modified on output after applying PBC
     * 
     */
    virtual void getPeriodicDelta(Vec3& deltaR) const {};
};


class AmoebaReferencePmeHippoNonbondedForce : public AmoebaReferenceHippoNonbondedForce {

public:

    /**
     * Constructor
     * 
     */
    AmoebaReferencePmeHippoNonbondedForce(const HippoNonbondedForce& force, const System& system);
 
    /**
     * Destructor
     * 
     */
    ~AmoebaReferencePmeHippoNonbondedForce();
 
    /**
     * Get cutoff distance.
     *
     * @return cutoff distance
     *
     */
    double getCutoffDistance() const;

    /**
     * Set cutoff distance.
     *
     * @return cutoff distance
     *
     */
    void setCutoffDistance(double cutoffDistance);

    /**
     * Get alpha used in Ewald summation.
     *
     * @return alpha
     *
     */
    double getAlphaEwald() const;

    /**
     * Get alpha used in dispersion Ewald summation.
     *
     * @return alpha
     *
     */
    double getDispersionAlphaEwald() const;

    /**
     * Get PME grid dimensions.
     *
     * @param pmeGridDimensions contains PME grid dimensions upon return

     *
     */
    void getPmeGridDimensions(std::vector<int>& pmeGridDimensions) const;

    /**
     * Get PME grid dimensions.
     *
     * @param pmeGridDimensions contains PME grid dimensions upon return

     *
     */
    void getDispersionPmeGridDimensions(std::vector<int>& pmeGridDimensions) const;

    /**
     * Set PME grid dimensions.
     *
     * @param pmeGridDimensions input PME grid dimensions 
     *
     */
    void setPmeGridDimensions(std::vector<int>& pmeGridDimensions);

    /**
     * Set PME grid dimensions.
     *
     * @param pmeGridDimensions input PME grid dimensions 
     *
     */
    void setDispersionPmeGridDimensions(std::vector<int>& pmeGridDimensions);

    /**
     * Set periodic box size.
     *
     * @param vectors    the vectors defining the periodic box
     */
     void setPeriodicBoxSize(OpenMM::Vec3* vectors);

private:

    static const int AMOEBA_PME_ORDER;
    static const double SQRT_PI;

    double _alphaEwald, _dalphaEwald;
    double _cutoffDistance;
    double _cutoffDistanceSquared;

    Vec3 _recipBoxVectors[3];
    Vec3 _periodicBoxVectors[3];

    int _totalGridSize;
    HippoIntVec _pmeGridDimensions, _dpmeGridDimensions;

    fftpack_t   _fftplan;

    unsigned int _pmeGridSize;
    t_complex* _pmeGrid;
 
    std::vector<double> _pmeBsplineModuli[3];
    std::vector<HippoDouble4> _thetai[3];
    std::vector<HippoIntVec> _iGrid;
    std::vector<double> _phi;
    std::vector<double> _phid;
    std::vector<double> _phip;
    std::vector<double> _phidp;
    std::vector<HippoDouble4> _pmeBsplineTheta;
    std::vector<HippoDouble4> _pmeBsplineDtheta;

    /**
     * Resize PME arrays.
     * 
     */
    void resizePmeArrays();

    /**
     * Zero Pme grid.
     */
    void initializePmeGrid();

    /**
     * Modify input vector of differences in particle positions for periodic boundary conditions.
     * 
     * @param delta                   input vector of difference in particle positions; on output adjusted for
     *                                periodic boundary conditions
     */
    void getPeriodicDelta(Vec3& deltaR) const;

    /**
     * Calculate damped inverse distances.
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     * @param dampedDInverseDistances damped inverse distances (drr3,drr5,drr7 in udirect2a() in TINKER)
     * @param dampedPInverseDistances damped inverse distances (prr3,prr5,prr7 in udirect2a() in TINKER)
     */
    void getDampedInverseDistances(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                   double dscale, double pscale, double r,
                                   Vec3& dampedDInverseDistances, Vec3& dampedPInverseDistances) const;
    
    /**
     * Initialize B-spline moduli.
     * 
     */
    void initializeBSplineModuli();

    /**
     * Calculate direct-space field at site I due fixed multipoles at site J and vice versa.
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    void calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                             double dscale, double pscale);
    
    /**
     * Calculate fixed multipole fields.
     *
     * @param particleData vector particle data
     * 
     */
    void calculateFixedMultipoleField();

    /**
     * This is called from computeAmoebaBsplines().  It calculates the spline coefficients for a single atom along a single axis.
     * 
     * @param thetai output spline coefficients
     * @param w offset from grid point
     */
    void computeBSplinePoint(std::vector<HippoDouble4>& thetai, double w);
    
    /**
     * Compute bspline coefficients.
     *
     * @param particleData   vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void computeAmoebaBsplines(const std::vector<MultipoleParticleData>& particleData);

    /**
     * Transform multipoles from cartesian coordinates to fractional coordinates.
     */
    void transformMultipolesToFractionalCoordinates(const std::vector<MultipoleParticleData>& particleData);

    /**
     * Transform potential from fractional coordinates to cartesian coordinates.
     */
    void transformPotentialToCartesianCoordinates(const std::vector<double>& fphi, std::vector<double>& cphi) const;

    /**
     * Spread fixed multipoles onto PME grid.
     * 
     * @param particleData vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void spreadFixedMultipolesOntoGrid(const std::vector<MultipoleParticleData>& particleData);

    /**
     * Perform reciprocal convolution.
     * 
     */
    void performAmoebaReciprocalConvolution();

    /**
     * Compute reciprocal potential due fixed multipoles at each particle site.
     * 
     */
    void computeFixedPotentialFromGrid(void);

    /**
     * Compute reciprocal potential due fixed multipoles at each particle site.
     * 
     */
    void computeInducedPotentialFromGrid();

    /**
     * Calculate reciprocal space energy and force due to fixed multipoles.
     * 
     * @param particleData    vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces          upon return updated vector of forces
     * @param torques         upon return updated vector of torques
     *
     * @return energy
     */
    double computeReciprocalSpaceFixedMultipoleForceAndEnergy(const std::vector<MultipoleParticleData>& particleData,
                                                              std::vector<Vec3>& forces, std::vector<Vec3>& torques) const;

    /**
     * Set reciprocal space fixed multipole fields.
     * 
     */
    void recordFixedMultipoleField();

    /**
     * Compute the potential due to the reciprocal space PME calculation for induced dipoles.
     *
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void calculateReciprocalSpaceInducedDipoleField(std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields);

    /**
     * Calculate field at particleI due to induced dipole at particle J and vice versa.
     *
     * @param iIndex        particle I index
     * @param jIndex        particle J index
     * @param preFactor1    first factor used in calculating field
     * @param preFactor2    second factor used in calculating field
     * @param delta         delta in particle positions after adjusting for periodic boundary conditions
     * @param inducedDipole vector of induced dipoles
     * @param field         vector of field at each particle due induced dipole of other particles
     */
    void calculateDirectInducedDipolePairIxn(unsigned int iIndex, unsigned int jIndex,
                                             double preFactor1, double preFactor2, const Vec3& delta,
                                             const std::vector<Vec3>& inducedDipole,
                                             std::vector<Vec3>& field) const;

    /**
     * Calculate direct space field at particleI due to induced dipole at particle J and vice versa for
     * inducedDipole and inducedDipolePolar.
     * 
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void calculateDirectInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                              const MultipoleParticleData& particleJ,
                                              std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields);

    /**
     * Initialize induced dipoles
     *
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void initializeInducedDipoles(std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields); 

    /**
     * Spread induced dipoles onto grid.
     *
     * @param inputInducedDipole      induced dipole value
     * @param inputInducedDipolePolar induced dipole polar value
     */
    void spreadInducedDipolesOnGrid(const std::vector<Vec3>& inputInducedDipole,
                                    const std::vector<Vec3>& inputInducedDipolePolar);

    /**
     * Calculate induced dipole fields.
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void calculateInducedDipoleFields(const std::vector<MultipoleParticleData>& particleData,
                                      std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields);

    /**
     * Set reciprocal space induced dipole fields. 
     *
     * @param field       reciprocal space output induced dipole field value at each site
     * @param fieldPolar  reciprocal space output induced dipole polar field value at each site
     * 
     */
    void recordInducedDipoleField(std::vector<Vec3>& field, std::vector<Vec3>& fieldPolar);

    /**
     * Compute Pme self energy.
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     */
    double calculatePmeSelfEnergy(const std::vector<MultipoleParticleData>& particleData) const;

    /**
     * Compute the self torques.
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 vector of torques
     */
    void calculatePmeSelfTorque(const std::vector<MultipoleParticleData>& particleData, std::vector<Vec3>& torques) const;

    /**
     * Calculate direct space electrostatic interaction between particles I and J.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param scalingFactors    scaling factors for interaction
     * @param forces            vector of particle forces to be updated
     * @param torques           vector of particle torques to be updated
     */
    double calculatePmeDirectElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                  const std::vector<double>& scalingFactors,
                                                  std::vector<Vec3>& forces, std::vector<Vec3>& torques) const;

    /**
     * Calculate reciprocal space energy/force/torque for dipole interaction.
     * 
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces            vector of particle forces to be updated
     * @param torques           vector of particle torques to be updated
     */
     double computeReciprocalSpaceInducedDipoleForceAndEnergy(const std::vector<MultipoleParticleData>& particleData,
                                                              std::vector<Vec3>& forces, std::vector<Vec3>& torques) const;

    /**
     * Calculate electrostatic forces.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    double calculateElectrostatic(const std::vector<MultipoleParticleData>& particleData, 
                                  std::vector<OpenMM::Vec3>& torques,
                                  std::vector<OpenMM::Vec3>& forces);

};

} // namespace OpenMM

#endif // _AmoebaReferenceHippoNonbondedForce___
