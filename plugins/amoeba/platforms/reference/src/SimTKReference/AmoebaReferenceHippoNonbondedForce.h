/* Portions copyright (c) 2006-2019 Stanford University and Simbios.
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
#include <array>
#include <map>
#include <utility>
#include <vector>
#include "fftpack.h"
#include <complex>

namespace OpenMM {

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

protected:

    enum QuadrupoleIndices { QXX, QXY, QXZ, QYY, QYZ, QZZ };

    /* 
     * Particle parameters and coordinates
     */
    class MultipoleParticleData {
        public:
            int index, axisType, multipoleAtomX, multipoleAtomY, multipoleAtomZ;    
            Vec3 position;
            Vec3 dipole, localDipole, qiDipole, qiInducedDipole;
            double quadrupole[6], localQuadrupole[6], qiQuadrupole[6];
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
        double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
    };

    unsigned int _numParticles;
    HippoNonbondedForce::NonbondedMethod _nonbondedMethod;
    double _electric, _cutoffDistance, _cutoffDistanceSquared, _switchingDistance;
    std::map<std::pair<int, int>, Exception> exceptions;
    std::vector<MultipoleParticleData> particleData;
    std::vector<TransformedMultipole> _transformed;
    std::vector<Vec3> _fixedMultipoleField;
    std::vector<Vec3> _inducedDipole;
    std::vector<Vec3> _inducedDipoleField;
    std::vector<std::vector<Vec3> > _ptDipoleD;

    int _maxPTOrder;
    std::vector<double> _extrapolationCoefficients;
    std::vector<double> _extPartCoefficients;

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
     * Compute the damping factors for the field due to the fixed multipole moments of a particle.
     * 
     * @param particle     parameters for the particle
     * @param r            the distance from the particle at which to calculate the damping factors
     * @param fdamp3       outputs the damping factor for the r^-3 term
     * @param fdamp5       outputs the damping factor for the r^-5 term
     * @param fdamp7       outputs the damping factor for the r^-7 term
     */
    void computeDirectFieldDampingFactors(const MultipoleParticleData& particle, double r, double& fdamp3, double& fdamp5, double& fdamp7) const;
    /**
     * Compute the damping factors for the field at one particle due to the induced dipole of another particle.
     * 
     * @param particleI    parameters for the first particle
     * @param particleJ    parameters for the second particle
     * @param r            the distance between the two particles
     * @param fdamp3       outputs the damping factor for the r^-3 term
     * @param fdamp5       outputs the damping factor for the r^-5 term
     */
    void computeMutualFieldDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r, double& fdamp3, double& fdamp5) const;
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
            double& fdampIJ1, double& fdampIJ3, double& fdampIJ5, double& fdampIJ7, double& fdampIJ9, double& fdampIJ11) const;
    /**
     * Compute the damping factors used for dispersion.
     * 
     * @param particleI    parameters for the first particle
     * @param particleJ    parameters for the second particle
     * @param r            the distance between the two particles
     * @param fdamp        outputs the damping factor
     * @param ddamp        outputs the derivative of the damping factor
     */
    void computeDispersionDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r, double& fdamp, double& ddamp) const;
    /**
     * Compute the damping factors used for Pauli repulsion.
     * 
     * @param particleI    parameters for the first particle
     * @param particleJ    parameters for the second particle
     * @param r            the distance between the two particles
     * @param fdamp1       outputs the damping factor for the r^-1 term
     * @param fdamp3       outputs the damping factor for the r^-3 term
     * @param fdamp5       outputs the damping factor for the r^-5 term
     * @param fdamp7       outputs the damping factor for the r^-7 term
     * @param fdamp9       outputs the damping factor for the r^-9 term
     * @param fdamp11       outputs the damping factor for the r^-11 term
     */
    void computeRepulsionDampingFactors(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ, double r, double& fdamp1, double& fdamp3,
                                        double& fdamp5, double& fdamp7, double& fdamp9, double& fdamp11) const;
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
     * @param deltaR                the internuclear vector, corrected for periodic boundary conditions
     * @param r                     the distance between atoms I and J
     * @param rotationmatrix        the output rotation matrix for a 3-vector
     */
    void formQIRotationMatrix(const Vec3 &deltaR, double r, double (&rotationMatrix)[3][3]) const;

    /**
     * Rotate a vector from the lab frame to the quasi-internal coordinate system.
     *
     * @param v      the vector to rotate
     * @param mat    the rotation matrix computed by formQIRotationMatrix()
     * @return the rotated vector
     */
    Vec3 rotateVectorToQI(const Vec3 v, const double (&mat)[3][3]) const;

    /**
     * Rotate a vector from the quasi-internal coordinate system back to the lab frame.
     *
     * @param v      the vector to rotate
     * @param mat    the rotation matrix computed by formQIRotationMatrix()
     * @return the rotated vector
     */
    Vec3 rotateVectorFromQI(const Vec3 v, const double (&mat)[3][3]) const;

    /**
     * Rotate a quadrupole from the lab frame to the quasi-internal coordinate system.
     *
     * @param q      the quadrupole to rotate
     * @param r      the rotated quadrupole to stored into this
     * @param mat    the rotation matrix computed by formQIRotationMatrix()
     */
    void rotateQuadrupoleToQI(const double (&q)[6], double (&r)[6], const double (&mat)[3][3]) const;

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
     * Calculate electric field at particle I due fixed multipoles at particle J and vice versa
     * (field at particle J due fixed multipoles at particle I).
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     */
    virtual void calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ);

    /**
     * Initialize induced dipoles
     */
    virtual void initializeInducedDipoles(); 

    /**
     * Calculate fields due induced dipoles at each site.
     *
     * @param particleI     positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ     positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     */
    void calculateInducedDipolePairIxns(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ);

    /**
     * Calculate induced dipole fields.
     * 
     * @param particleData    vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    virtual void calculateInducedDipoleFields(const std::vector<MultipoleParticleData>& particleData, int optOrder);
    /**
     * Calculated induced dipoles using extrapolated perturbation theory.
     *
     * @param particleData   vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void convergeInduceDipolesByExtrapolation(const std::vector<MultipoleParticleData>& particleData);

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
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueI           the torque to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueK           the torque to apply (in the quasi-internal frame) to particle K should be added to this
     */
    virtual double calculateElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                                 double r, Vec3& force, Vec3& torqueI, Vec3& torqueK) const;

    /**
     * Calculate electrostatic interactions involving induced dipoles on particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param deltaR            the displacement between the two particles (in the lab frame)
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueI           the torque to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueK           the torque to apply (in the quasi-internal frame) to particle K should be added to this
     * @param labForce          the force to apply (in the lab frame) to particle I should be added to this
     */
    virtual void calculateInducedDipolePairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                               Vec3 deltaR, double r, Vec3& force, Vec3& torqueI, Vec3& torqueK, Vec3& labForce) const;

    /**
     * Calculate dispersion interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     */
    virtual double calculateDispersionPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                              double r, Vec3& force) const;

    /**
     * Calculate the Pauli repulsion interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueI           the torque to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueK           the torque to apply (in the quasi-internal frame) to particle K should be added to this
     */
    double calculateRepulsionPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                     double r, Vec3& force, Vec3& torqueI, Vec3& torqueK) const;

    /**
     * Calculate the charge transfer interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     */
    double calculateChargeTransferPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                          double r, Vec3& force) const;

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
     * Calculate the forces and energy
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    virtual double calculateInteractions(std::vector<OpenMM::Vec3>& torques,
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
     * Initialize vector of Vec3 (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeVec3Vector(std::vector<Vec3>& vectorToInitialize) const;

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

    Vec3 _recipBoxVectors[3];
    Vec3 _periodicBoxVectors[3];

    int _totalGridSize;
    int _pmeGridDimensions[3];
    int _dpmeGridDimensions[3];

    fftpack_t   _fftplan;

    std::vector<t_complex> _pmeGrid;
 
    std::vector<double> _pmeBsplineModuli[3];
    std::vector<HippoDouble4> _thetai[3];
    std::vector<std::array<int, 3> > _iGrid;
    std::vector<double> _phi;
    std::vector<double> _phidp;
    std::vector<std::vector<double> > optPhi;
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
     * Initialize B-spline moduli.
     * 
     */
    void initializeBSplineModuli();

    /**
     * Calculate direct-space field at site I due fixed multipoles at site J and vice versa.
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     */
    void calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ);
    
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
     */
    void calculateReciprocalSpaceInducedDipoleField();

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
    void calculateDirectInducedDipolePairIxn(int iIndex, int jIndex,
                                             double preFactor1, double preFactor2, const Vec3& delta,
                                             const std::vector<Vec3>& inducedDipole,
                                             std::vector<Vec3>& field) const;

    /**
     * Calculate direct space field at particleI due to induced dipole at particle J and vice versa for
     * inducedDipole and inducedDipolePolar.
     * 
     * @param particleI    positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ    positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     */
    void calculateDirectInducedDipolePairIxns(const MultipoleParticleData& particleI,
                                              const MultipoleParticleData& particleJ);

    /**
     * Initialize induced dipoles
     */
    void initializeInducedDipoles(); 

    /**
     * Spread induced dipoles onto grid.
     *
     * @param inputInducedDipole      induced dipole value
     */
    void spreadInducedDipolesOnGrid(const std::vector<Vec3>& inputInducedDipole);

    /**
     * Calculate induced dipole fields.
     * 
     * @param particleData   vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void calculateInducedDipoleFields(const std::vector<MultipoleParticleData>& particleData, int optOrder);

    /**
     * Set reciprocal space induced dipole fields. 
     *
     * @param field       reciprocal space output induced dipole field value at each site
     * 
     */
    void recordInducedDipoleField(std::vector<Vec3>& field);

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
     * Calculate direct space electrostatic interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueI           the torque to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueK           the torque to apply (in the quasi-internal frame) to particle K should be added to this
     */
    double calculateElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                         double r, Vec3& force, Vec3& torqueI, Vec3& torqueK) const;

    /**
     * Calculate electrostatic interactions involving induced dipoles on particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param deltaR            the displacement between the two particles (in the lab frame)
     * @param r                 the distance between the two particles
     * @param force             the force to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueI           the torque to apply (in the quasi-internal frame) to particle I should be added to this
     * @param torqueK           the torque to apply (in the quasi-internal frame) to particle K should be added to this
     * @param labForce          the force to apply (in the lab frame) to particle I should be added to this
     */
    void calculateInducedDipolePairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                       Vec3 deltaR, double r, Vec3& force, Vec3& torqueI, Vec3& torqueK, Vec3& labForce) const;

    /**
     * Calculate dispersion interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param deltaR            the displacement between the two particles
     * @param r                 the distance between the two particles
     * @param force             the force to apply to particle I should be added to this
     */
    double calculateDispersionPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                      double r, Vec3& force) const;

    /**
     * Calculate reciprocal space force/torque for dipole interaction.
     * 
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces            vector of particle forces to be updated
     * @param torques           vector of particle torques to be updated
     */
     void computeReciprocalSpaceInducedDipoleForce(const std::vector<MultipoleParticleData>& particleData,
                                                   std::vector<Vec3>& forces, std::vector<Vec3>& torques) const;

    /**
     * Calculate reciprocal space energy and force due to dispersion.
     * 
     * @param particleData    vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces          upon return updated vector of forces
     *
     * @return energy
     */
    double computeReciprocalSpaceDispersionForceAndEnergy(const std::vector<MultipoleParticleData>& particleData, std::vector<Vec3>& forces) const;

    /**
     * Calculate the forces and energy.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    double calculateInteractions(std::vector<OpenMM::Vec3>& torques,
                                 std::vector<OpenMM::Vec3>& forces);

};

} // namespace OpenMM

#endif // _AmoebaReferenceHippoNonbondedForce___
