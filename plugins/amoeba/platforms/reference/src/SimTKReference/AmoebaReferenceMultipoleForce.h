/* Portions copyright (c) 2006-2015 Stanford University and Simbios.
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

#include "RealVec.h"
#include "openmm/AmoebaMultipoleForce.h"
#include "AmoebaReferenceGeneralizedKirkwoodForce.h"
#include <map>
#include "fftpack.h"
#include <complex>

namespace OpenMM {

typedef std::map< unsigned int, RealOpenMM> MapIntRealOpenMM;
typedef MapIntRealOpenMM::iterator MapIntRealOpenMMI;
typedef MapIntRealOpenMM::const_iterator MapIntRealOpenMMCI;


// A few useful constants for the spherical harmonic multipole code.
const RealOpenMM oneThird = 1.0/3.0;
const RealOpenMM twoThirds = 2.0/3.0;
const RealOpenMM fourThirds = 4.0/3.0;
const RealOpenMM fourSqrtOneThird = 4.0/sqrt(3.0);
const RealOpenMM sqrtFourThirds = 2.0/sqrt(3.0);
const RealOpenMM sqrtOneThird = 1.0/sqrt(3.0);
const RealOpenMM sqrtThree = sqrt(3.0);
const RealOpenMM oneNinth = 1.0/9.0;
const RealOpenMM fourOverFortyFive = 4.0/45.0;
const RealOpenMM fourOverFifteen = 4.0/15.0;


/**
 * 2-dimensional int vector
 */
class int2 {
public:
    /**
     * Create a int2 whose elements are all 0.
     */
    int2() {
        data[0] = data[1] = 0;
    }
    /**
     * Create a int2 with specified x, y components.
     */
    int2(int x, int y) {
        data[0] = x;
        data[1] = y;
    }
    int operator[](int index) const {
        assert(index >= 0 && index < 2);
        return data[index];
    }
    int& operator[](int index) {
        assert(index >= 0 && index < 2);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    int2 operator+() const {
        return int2(*this);
    }

    // plus
    int2 operator+(const int2& rhs) const {
        const int2& lhs = *this;
        return int2(lhs[0] + rhs[0], lhs[1] + rhs[1]);
    }

    int2& operator+=(const int2& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        return *this;
    }

    int2& operator-=(const int2& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        return *this;
    }

private:
    int data[2];
};

/**
 * 3-dimensional int vector
 */
class IntVec {
public:
    /**
     * Create a IntVec whose elements are all 0.
     */
    IntVec() {
        data[0] = data[1] = data[2] = 0;
    }
    /**
     * Create a IntVec with specified x, y, z, w components.
     */
    IntVec(int x, int y, int z) {
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
    IntVec operator+() const {
        return IntVec(*this);
    }

    // plus
    IntVec operator+(const IntVec& rhs) const {
        const IntVec& lhs = *this;
        return IntVec(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
    }

    IntVec& operator+=(const IntVec& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        return *this;
    }

    IntVec& operator-=(const IntVec& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        return *this;
    }

private:
    int data[3];
};

/**
 * 4-dimensional RealOpenMM vector
 */
class RealOpenMM4 {
public:
    /**
     * Create a RealOpenMM4 whose elements are all 0.
     */
    RealOpenMM4() {
        data[0] = data[1] = data[2] = data[3] = 0.0;
    }
    /**
     * Create a RealOpenMM4 with specified x, y, z, w components.
     */
    RealOpenMM4(RealOpenMM x, RealOpenMM y, RealOpenMM z, RealOpenMM w) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
    }
    RealOpenMM operator[](int index) const {
        assert(index >= 0 && index < 4);
        return data[index];
    }
    RealOpenMM& operator[](int index) {
        assert(index >= 0 && index < 4);
        return data[index];
    }

    // Arithmetic operators

    // unary plus
    RealOpenMM4 operator+() const {
        return RealOpenMM4(*this);
    }

    // plus
    RealOpenMM4 operator+(const RealOpenMM4& rhs) const {
        const RealOpenMM4& lhs = *this;
        return RealOpenMM4(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2],lhs[3] + rhs[3]);
    }

    RealOpenMM4& operator+=(const RealOpenMM4& rhs) {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];
        data[3] += rhs[3];
        return *this;
    }

    RealOpenMM4& operator-=(const RealOpenMM4& rhs) {
        data[0] -= rhs[0];
        data[1] -= rhs[1];
        data[2] -= rhs[2];
        data[3] -= rhs[3];
        return *this;
    }

private:
    RealOpenMM data[4];
};

using namespace OpenMM;

class AmoebaReferenceMultipoleForce {

   /**
    * AmoebaReferenceMultipoleForce is base class for MultipoleForce calculations
    * AmoebaReferenceGeneralizedKirkwoodMultipoleForce  is derived class for Generalized Kirkwood calculations
    * AmoebaReferencePmeMultipoleForce is derived class for PME calculations
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
    *                                                          AmoebaReferenceMultipoleForce::calculateElectrostatic()
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
     * This is an enumeration of the different methods that may be used for handling long range Multipole forces.
     */
    enum NonbondedMethod {

        /**
         * No cutoff is applied to the interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,

       /** 
         * Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        PME = 1
    };

    enum PolarizationType {

        /** 
         * Mutual polarization
         */
        Mutual = 0,

        /**
         * Direct polarization
         */
        Direct = 1,

        /**
         * Extrapolated perturbation theory
         */
        Extrapolated = 2
    };

    /**
     * Constructor
     * 
     */
    AmoebaReferenceMultipoleForce();
 
    /**
     * Constructor
     * 
     * @param nonbondedMethod nonbonded method
     */
    AmoebaReferenceMultipoleForce(NonbondedMethod nonbondedMethod);
 
    /**
     * Destructor
     * 
     */
    virtual ~AmoebaReferenceMultipoleForce() {};
 
    /**
     * Get nonbonded method.
     * 
     * @return nonbonded method
     */
    NonbondedMethod getNonbondedMethod() const;

    /**
     * Set nonbonded method.
     * 
     * @param nonbondedMethod nonbonded method
     */
    void setNonbondedMethod(NonbondedMethod nonbondedMethod);

    /**
     * Get polarization type.
     * 
     * @return polarization type
     */
    PolarizationType getPolarizationType() const;

    /**
     * Set polarization type.
     * 
     * @param  polarizationType polarization type
     */
    void setPolarizationType(PolarizationType polarizationType);

    /**
     * Get flag indicating if mutual induced dipoles are converged.
     *
     * @return nonzero if converged
     *
     */
    int getMutualInducedDipoleConverged() const;

    /**
     * Get the number of iterations used in computing mutual induced dipoles.
     *
     * @return number of iterations
     *
     */
    int getMutualInducedDipoleIterations() const;

    /**
     * Get the final epsilon for mutual induced dipoles.
     *
     *  @return epsilon
     *
     */
    RealOpenMM getMutualInducedDipoleEpsilon() const;

    /**
     * Set the coefficients for the µ_0, µ_1, µ_2, µ_n terms in the extrapolation
     * theory algorithm for induced dipoles
     *
     * @param optCoefficients a vector whose mth entry specifies the coefficient for µ_m
     *
     */
    void setExtrapolationCoefficients(const std::vector<RealOpenMM> &coefficients);

    /**
     * Set the target epsilon for converging mutual induced dipoles.
     *
     * @param targetEpsilon target epsilon for converging mutual induced dipoles
     *
     */
    void setMutualInducedDipoleTargetEpsilon(RealOpenMM targetEpsilon);

    /**
     * Get the target epsilon for converging mutual induced dipoles.
     *
     * @return target epsilon for converging mutual induced dipoles
     *
     */
    RealOpenMM getMutualInducedDipoleTargetEpsilon() const;

    /**
     * Set the maximum number of iterations to be executed in converging mutual induced dipoles.
     *
     * @param maximumMutualInducedDipoleIterations maximum number of iterations to be executed in converging mutual induced dipoles
     *
     */
    void setMaximumMutualInducedDipoleIterations(int maximumMutualInducedDipoleIterations);

    /**
     * Get the maximum number of iterations to be executed in converging mutual induced dipoles.
     *
     * @return maximum number of iterations to be executed in converging mutual induced dipoles
     * 
     */
    int getMaximumMutualInducedDipoleIterations() const;

    /**
     * Calculate force and energy.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            damping factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param forces                    add forces to this vector
     *
     * @return energy
     */
    RealOpenMM calculateForceAndEnergy(const std::vector<OpenMM::RealVec>& particlePositions,
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
                                       std::vector<OpenMM::RealVec>& forces);

    /**
     * Calculate particle induced dipoles.
     *
     * @param masses                    particle masses
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param outputMultipoleMoments    output multipole moments
     */
    void calculateInducedDipoles(const std::vector<OpenMM::RealVec>& particlePositions,
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
                                 std::vector<RealVec>& outputInducedDipoles);

    /**
     * Calculate particle permanent dipoles rotated in the lab frame.
     *
     * @param masses                    particle masses
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param outputMultipoleMoments    output multipole moments
     */

    void calculateLabFramePermanentDipoles(const std::vector<RealVec>& particlePositions,
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
                                           const std::vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                           std::vector<RealVec>& outputRotatedPermanentDipoles);

    /**
     * Calculate particle total dipoles.
     *
     * @param masses                    particle masses
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param outputMultipoleMoments    output multipole moments
     */


    void calculateTotalDipoles(const std::vector<RealVec>& particlePositions,
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
                                           const std::vector< vector< vector<int> > >& multipoleAtomCovalentInfo,
                                           std::vector<RealVec>& outputRotatedPermanentDipoles);



    /**
     * Calculate system multipole moments.
     *
     * @param masses                    particle masses
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param outputMultipoleMoments    output multipole moments
     */
    void calculateAmoebaSystemMultipoleMoments(const std::vector<RealOpenMM>& masses,
                                               const std::vector<OpenMM::RealVec>& particlePositions,
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
                                               std::vector<RealOpenMM>& outputMultipoleMoments);

    /**
     * Calculate electrostatic potential at a set of grid points.
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param input grid                input grid points to compute potential
     * @param outputPotential           output electrostatic potential
     */
    void calculateElectrostaticPotential(const std::vector<OpenMM::RealVec>& particlePositions,
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
                                         const std::vector<RealVec>& inputGrid,
                                         std::vector<RealOpenMM>& outputPotential);

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
            RealVec sphericalDipole;
            RealOpenMM sphericalQuadrupole[5];
            RealOpenMM thole;
            RealOpenMM dampingFactor;
            RealOpenMM polarity;
    };
    
    /**
     * Particle parameters transformed into fractional coordinates
     */
    class TransformedMultipole {
    public:
        RealOpenMM charge;
        RealVec dipole;
        RealOpenMM quadrupole[6];
    };

    /* 
     * Helper class used in calculating induced dipoles
     */
    struct UpdateInducedDipoleFieldStruct {
            UpdateInducedDipoleFieldStruct(std::vector<OpenMM::RealVec>& inputFixed_E_Field, std::vector<OpenMM::RealVec>& inputInducedDipoles, std::vector<std::vector<RealVec> >& extrapolatedDipoles, std::vector<std::vector<RealOpenMM> >& extrapolatedDipoleFieldGradient);
            std::vector<OpenMM::RealVec>* fixedMultipoleField;
            std::vector<OpenMM::RealVec>* inducedDipoles;
            std::vector<std::vector<RealVec> >* extrapolatedDipoles;
            std::vector<std::vector<RealOpenMM> >* extrapolatedDipoleFieldGradient;
            std::vector<OpenMM::RealVec> inducedDipoleField;
            std::vector<std::vector<RealOpenMM> > inducedDipoleFieldGradient;
    };

    unsigned int _numParticles;

    NonbondedMethod _nonbondedMethod;
    PolarizationType _polarizationType;

    RealOpenMM _electric;
    RealOpenMM _dielectric;

    enum ScaleType { D_SCALE, P_SCALE, M_SCALE, U_SCALE, LAST_SCALE_TYPE_INDEX };
    std::vector<  std::vector< MapIntRealOpenMM > > _scaleMaps;
    std::vector<unsigned int> _maxScaleIndex;
    RealOpenMM _dScale[5];
    RealOpenMM _pScale[5];
    RealOpenMM _mScale[5];
    RealOpenMM _uScale[5];

    std::vector<TransformedMultipole> _transformed;
    std::vector<RealVec> _fixedMultipoleField;
    std::vector<RealVec> _fixedMultipoleFieldPolar;
    std::vector<RealVec> _inducedDipole;
    std::vector<RealVec> _inducedDipolePolar;
    std::vector<std::vector<RealVec> > _ptDipoleP;
    std::vector<std::vector<RealVec> > _ptDipoleD;
    std::vector<std::vector<RealOpenMM> > _ptDipoleFieldGradientP;
    std::vector<std::vector<RealOpenMM> > _ptDipoleFieldGradientD;

    int _mutualInducedDipoleConverged;
    int _mutualInducedDipoleIterations;
    int _maximumMutualInducedDipoleIterations;
    int _maxPTOrder;
    std::vector<RealOpenMM>  _extrapolationCoefficients;
    std::vector<RealOpenMM>  _extPartCoefficients;
    RealOpenMM  _mutualInducedDipoleEpsilon;
    RealOpenMM  _mutualInducedDipoleTargetEpsilon;
    RealOpenMM  _polarSOR;
    RealOpenMM  _debye;

    /**
     * Helper constructor method to centralize initialization of objects.
     *
     */
    void initialize();

    /**
     * Load particle data.
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
    void loadParticleData(const std::vector<OpenMM::RealVec>& particlePositions, 
                          const std::vector<RealOpenMM>& charges,
                          const std::vector<RealOpenMM>& dipoles,
                          const std::vector<RealOpenMM>& quadrupoles,
                          const std::vector<RealOpenMM>& tholes,
                          const std::vector<RealOpenMM>& dampingFactors,
                          const std::vector<RealOpenMM>& polarity,
                          std::vector<MultipoleParticleData>& particleData) const;

    /**
     * Calculate fixed multipole fields.
     *
     * @param particleData vector of particle data
     * 
     */
    virtual void calculateFixedMultipoleField(const vector<MultipoleParticleData>& particleData);

    /**
     * Set flag indicating if mutual induced dipoles are converged.
     * 
     * @param converged nonzero if converged
     *
     */
    void setMutualInducedDipoleConverged(int converged);

    /**
     * Set number of iterations used in computing mutual induced dipoles.
     * 
     * @param  number of iterations
     * 
     */
    void setMutualInducedDipoleIterations(int iterations);

    /**
     * Set the final epsilon for mutual induced dipoles.
     * 
     * @param epsilon
     *
     */
    void setMutualInducedDipoleEpsilon(RealOpenMM epsilon);

    /**
     * Setup scale factors given covalent info.
     *
     * @param  multipoleAtomCovalentInfo vector of vectors containing the covalent info
     *
     */
    void setupScaleMaps(const std::vector< std::vector< std::vector<int> > >& multipoleAtomCovalentInfo);

    /**
     * Get multipole scale factor for particleI & particleJ
     * 
     * @param  particleI           index of particleI whose scale factor is to be retrieved
     * @param  particleJ           index of particleJ whose scale factor is to be retrieved
     * @param  scaleType           scale type (D_SCALE, P_SCALE, M_SCAL)
     *
     * @return scaleFactor 
     */
    RealOpenMM getMultipoleScaleFactor(unsigned int particleI, unsigned int particleJ, ScaleType scaleType) const;

    /**
     * Get scale factor for particleI & particleJ
     * 
     * @param  particleI           index of particleI whose scale factor is to be retrieved
     * @param  particleJ           index of particleJ whose scale factor is to be retrieved
     * @param  scaleType           scale type (D_SCALE, P_SCALE, M_SCAL)
     *
     * @return array of scaleFactors 
     */
    void getMultipoleScaleFactors(unsigned int particleI, unsigned int particleJ, std::vector<RealOpenMM>& scaleFactors) const;

    /**
     * Get p- and d-scale factors for particleI & particleJ ixn
     *
     * @param  particleI           index of particleI 
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */
    void getDScaleAndPScale(unsigned int particleI, unsigned int particleJ, RealOpenMM& dScale, RealOpenMM& pScale) const;
    
    /**
     * Calculate damped powers of 1/r.
     *
     * @param  particleI           index of particleI 
     * @param  particleJ           index of particleJ
     * @param  dScale              output d-scale factor
     * @param  pScale              output p-scale factor
     */
    void getAndScaleInverseRs(RealOpenMM dampI, RealOpenMM dampJ, RealOpenMM tholeI, RealOpenMM tholeJ,
                               RealOpenMM r, std::vector<RealOpenMM>& rrI) const;

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
                                     MultipoleParticleData& particleX, MultipoleParticleData& particleY) const;

    /**
     * Invert multipole moments (dipole[Y], quadrupole[XY] and quadrupole[YZ]) if chiral center inverted.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param multipoleAtomXs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomZs         vector of y-particle indices used to map molecular frame to lab frame
     * @param axisType                axis type
     */
    void checkChiral(std::vector<MultipoleParticleData>& particleData, 
                     const std::vector<int>& multipoleAtomXs,
                     const std::vector<int>& multipoleAtomYs,
                     const std::vector<int>& multipoleAtomZs,
                     const std::vector<int>& axisTypes) const;
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
    void formQIRotationMatrix(const RealVec& iPosition,
                              const RealVec& jPosition,
                              const RealVec &deltaR,
                              RealOpenMM r,
                              RealOpenMM (&rotationMatrix)[3][3]) const;


    /**
     * Constructs a rotation matrix for spherical harmonic quadrupoles, using the dipole rotation matrix.
     *
     * @param D1                    The input spherical harmonic dipole rotation matrix
     * @param D2                    The output spherical harmonic quadrupole rotation matrix
     */
     void buildSphericalQuadrupoleRotationMatrix(const RealOpenMM (&D1)[3][3], RealOpenMM (&D2)[5][5]) const;

     /**
      * Constructs a rotation matrix for spherical harmonic quadrupoles, using the dipole rotation matrix.
      * Only the m={0,1c,1s} terms are constructed; these are the only terms needed to evaluate the field.
      *
      * @param D1                    The input spherical harmonic dipole rotation matrix
      * @param D2                    The output spherical harmonic quadrupole rotation matrix
      */
      void buildPartialSphericalQuadrupoleRotationMatrix(const RealOpenMM (&D1)[3][3], RealOpenMM (&D2)[3][5]) const;

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
    void applyRotationMatrix(std::vector<MultipoleParticleData>& particleData, 
                             const std::vector<int>& multipoleAtomXs,
                             const std::vector<int>& multipoleAtomYs,
                             const std::vector<int>& multipoleAtomZs,
                             const std::vector<int>& axisTypes) const;
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
                                                     RealOpenMM dScale, RealOpenMM pScale);

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
                                       RealOpenMM rr3, RealOpenMM rr5, const RealVec& delta,
                                       const std::vector<RealVec>& inducedDipole,
                                       std::vector<RealVec>& field) const;

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
     * Converge induced dipoles.
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void convergeInduceDipolesBySOR(const std::vector<MultipoleParticleData>& particleData,
                                    std::vector<UpdateInducedDipoleFieldStruct>& calculateInducedDipoleField);
    /**
     * Converge induced dipoles.
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void convergeInduceDipolesByDIIS(const std::vector<MultipoleParticleData>& particleData,
                                     std::vector<UpdateInducedDipoleFieldStruct>& calculateInducedDipoleField);
    
    /**
     * Use DIIS to compute the weighting coefficients for the new induced dipoles.
     * 
     * @param prevErrors    the vector of errors from previous iterations
     * @param coefficients  the coefficients will be stored into this
     */
    void computeDIISCoefficients(const std::vector<std::vector<RealVec> >& prevErrors, std::vector<RealOpenMM>& coefficients) const;

    /**
     * Update fields due to induced dipoles for each particle.
     * 
     * @param particleData              vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    RealOpenMM updateInducedDipoleFields(const std::vector<MultipoleParticleData>& particleData,
                                         std::vector<UpdateInducedDipoleFieldStruct>& calculateInducedDipoleField);

    /**
     * Update induced dipole for a particle given updated induced dipole field at the site.
     * 
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param fixedMultipoleField       fields due fixed multipoles at each site
     * @param inducedDipoleField        fields due induced dipoles at each site
     * @param inducedDipoles            output vector of updated induced dipoles
     */
    RealOpenMM updateInducedDipole(const std::vector<MultipoleParticleData>& particleI,
                                   const std::vector<RealVec>& fixedMultipoleField,
                                   const std::vector<RealVec>& inducedDipoleField,
                                   std::vector<RealVec>& inducedDipoles);

    /**
     * Calculate induced dipoles.
     * 
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    virtual void calculateInducedDipoles(const std::vector<MultipoleParticleData>& particleData);

    /**
     * Setup: 
     *        if needed invert multipole moments at chiral centers
     *        rotate molecular multipole moments to lab frame 
     *        setup scaling maps and 
     *        calculate induced dipoles (see calculateInducedDipoles below)
     *
     * @param particlePositions         Cartesian coordinates of particles
     * @param charges                   scalar charges for each particle
     * @param dipoles                   molecular frame dipoles for each particle
     * @param quadrupoles               molecular frame quadrupoles for each particle
     * @param tholes                    Thole factors for each particle
     * @param dampingFactors            dampling factors for each particle
     * @param polarity                  polarity for each particle
     * @param axisTypes                 axis type (Z-then-X, ...) for each particle
     * @param multipoleAtomZs           indicies of particle specifying the molecular frame z-axis for each particle
     * @param multipoleAtomXs           indicies of particle specifying the molecular frame x-axis for each particle
     * @param multipoleAtomYs           indicies of particle specifying the molecular frame y-axis for each particle
     * @param multipoleAtomCovalentInfo covalent info needed to set scaling factors
     * @param particleData              output vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     *
     */
    void setup(const std::vector<OpenMM::RealVec>& particlePositions,
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
               std::vector<MultipoleParticleData>& particleData);

    /**
     * Calculate electrostatic interaction between particles I and K.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleK         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle K
     * @param scalingFactors    scaling factors for interaction
     * @param forces            vector of particle forces to be updated
     * @param torque            vector of particle torques to be updated
     */
    RealOpenMM calculateElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleK,
                                             const std::vector<RealOpenMM>& scalingFactors, std::vector<OpenMM::RealVec>& forces, std::vector<RealVec>& torque) const;

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
                                     int axisType, const Vec3& torque, std::vector<OpenMM::RealVec>& forces) const;

    /**
     * Map torques to forces.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param multipoleAtomZs         vector of z-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomXs         vector of x-particle indices used to map molecular frame to lab frame
     * @param multipoleAtomYs         vector of y-particle indices used to map molecular frame to lab frame
     * @param axisType                vector of axis types (Bisector/Z-then-X, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     */
    void mapTorqueToForce(std::vector<MultipoleParticleData>& particleData, 
                          const std::vector<int>& multipoleAtomXs,
                          const std::vector<int>& multipoleAtomYs,
                          const std::vector<int>& multipoleAtomZs,
                          const std::vector<int>& axisTypes,
                          std::vector<OpenMM::RealVec>& torques,
                          std::vector<OpenMM::RealVec>& forces) const;

    /**
     * Calculate electrostatic forces
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    virtual RealOpenMM calculateElectrostatic(const std::vector<MultipoleParticleData>& particleData, 
                                              std::vector<OpenMM::RealVec>& torques,
                                              std::vector<OpenMM::RealVec>& forces);

    /**
     * Normalize a RealVec
     *
     * @param vectorToNormalize vector to normalize
     *
     * @return norm of vector on input
     * 
     */
    RealOpenMM normalizeRealVec(RealVec& vectorToNormalize) const;

    /**
     * Initialize vector of RealOpenMM (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeRealOpenMMVector(vector<RealOpenMM>& vectorToInitialize) const;

    /**
     * Initialize vector of RealVec (size=numParticles)
     *
     * @param vectorToInitialize vector to initialize
     * 
     */
    void initializeRealVecVector(vector<RealVec>& vectorToInitialize) const;

    /**
     * Copy vector of RealVec
     *
     * @param inputVector  vector to copy
     * @param outputVector output vector
     * 
     */
    void copyRealVecVector(const std::vector<OpenMM::RealVec>& inputVector, std::vector<OpenMM::RealVec>& outputVector) const;

    /**
     * Calculate potential at grid point due to a particle
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param gridPoint               grid point
     *
     * @return potential at grid point
     * 
     */
    RealOpenMM calculateElectrostaticPotentialForParticleGridPoint(const MultipoleParticleData& particleI, const RealVec& gridPoint) const;

    /**
     * Apply periodic boundary conditions to difference in positions
     * 
     * @param deltaR  difference in particle positions; modified on output after applying PBC
     * 
     */
    virtual void getPeriodicDelta(RealVec& deltaR) const {};
};

class AmoebaReferenceGeneralizedKirkwoodMultipoleForce : public AmoebaReferenceMultipoleForce {

public:

    /**
     * Constructor
     * 
     */
    AmoebaReferenceGeneralizedKirkwoodMultipoleForce(AmoebaReferenceGeneralizedKirkwoodForce* amoebaReferenceGeneralizedKirkwoodForce);
 
    /**
     * Destructor
     * 
     */
    ~AmoebaReferenceGeneralizedKirkwoodMultipoleForce();
 
    /**
     * Get flag signalling whether cavity term is to be included.
     *
     * @return flag
     *
     */
    int getIncludeCavityTerm() const;

    /**
     * Get probe radius.
     *
     * @return probe radius
     *
     */
    RealOpenMM getProbeRadius() const;

    /**
     * Get surface area factor.
     *
     * @return surface area factor
     *
     */
    RealOpenMM getSurfaceAreaFactor() const;

    /**
     * Get dielectric offset.
     *
     * @return dielectric offset
     *
     */
    RealOpenMM getDielectricOffset() const;

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
    std::vector<std::vector<RealVec> > _ptDipolePS;
    std::vector<std::vector<RealVec> > _ptDipoleDS;
    std::vector<std::vector<RealOpenMM> > _ptDipoleFieldGradientPS;
    std::vector<std::vector<RealOpenMM> > _ptDipoleFieldGradientDS;

    int _includeCavityTerm;
    RealOpenMM _probeRadius;
    RealOpenMM _surfaceAreaFactor;
    RealOpenMM _dielectricOffset;

    /**
     * Zero fixed multipole fields.
     *
     */
    void zeroFixedMultipoleFields();

    /**
     * Calculate electric field at particle I due fixed multipoles at particle J and vice versa
     * (field at particle J due fixed multipoles at particle I).
     * 
     * @param particleI               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dScale                  d-scale value for i-j interaction
     * @param pScale                  p-scale value for i-j interaction
     */
    void calculateFixedMultipoleFieldPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                             RealOpenMM dScale, RealOpenMM pScale);

    /**
     * Calculate induced dipoles.
     * 
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void calculateInducedDipoles(const std::vector<MultipoleParticleData>& particleData);

    /**
     * Calculate fields due induced dipoles at each site.
     * 
     * @param particleI                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ                 positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param updateInducedDipoleFields vector of UpdateInducedDipoleFieldStruct containing input induced dipoles and output fields
     */
    void calculateInducedDipolePairIxns(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                        std::vector<UpdateInducedDipoleFieldStruct>& updateInducedDipoleFields);

    /**
     * Calculate electrostatic forces and torques.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    RealOpenMM calculateElectrostatic(const std::vector<MultipoleParticleData>& particleData, 
                                      std::vector<OpenMM::RealVec>& torques,
                                      std::vector<OpenMM::RealVec>& forces);

    /**
     * Calculate GK field at particle I due induced dipole at particle J and vice versa
     * (field at particle J due induced dipole at particle I).
     * 
     * @param particleI               index of particle I
     * @param particleJ               index of particle J
     * @param field                   vector of induced dipole fields
     * @param fieldPolar              vector of induced dipole polar fields
     */
    void calculateInducedDipolePairGkIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                         const std::vector<RealVec>& field, std::vector<RealVec>& fieldPolar) const;

    /**
     * Calculate Kirkwood interaction.
     * 
     * @param particleI               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param forces                  add Kirkwood force to forces
     * @param torques                 add Kirkwood torque to torques
     * @param dBorn                   chain-rule factor
     *
     * @return energy
     */
    RealOpenMM calculateKirkwoodPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                        std::vector<RealVec>& forces, 
                                        std::vector<RealVec>& torques,
                                        std::vector<RealOpenMM>& dBorn) const;

    /**
     * Calculate Grycuk 'chain-rule' force.
     * 
     * @param particleI               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param dBorn                   chain-rule Born force factor
     * @param forces                  add Kirkwood force to forces
     *
     */
    void calculateGrycukChainRulePairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                         const std::vector<RealOpenMM>& dBorn, std::vector<RealVec>& forces) const;

    /**
     * Calculate TINKER's ACE approximation to non-polar cavity term. 
     * 
     * @param dBorn add ACE force to chain rule force
     *
     * @return ACE energy
     *
     */
    RealOpenMM calculateCavityTermEnergyAndForces(std::vector<RealOpenMM>& dBorn) const;

    /**
     * Correct vacuum to SCRF derivatives (TINKER's ediff1()).
     * 
     * @param particleI               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ               particle parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param pscale                  p-scale factor
     * @param dscale                  d-scale factor
     * @param forces                  force accumulator
     * @param torques                 torque accumulator
     *
     * @return energy
     */
    RealOpenMM calculateKirkwoodEDiffPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                             RealOpenMM pscale, RealOpenMM dscale, 
                                             std::vector<RealVec>& forces, std::vector<RealVec>& torques) const;

};

class AmoebaReferencePmeMultipoleForce : public AmoebaReferenceMultipoleForce {

public:

    /**
     * Constructor
     * 
     */
    AmoebaReferencePmeMultipoleForce();
 
    /**
     * Destructor
     * 
     */
    ~AmoebaReferencePmeMultipoleForce();
 
    /**
     * Get cutoff distance.
     *
     * @return cutoff distance
     *
     */
    RealOpenMM getCutoffDistance() const;

    /**
     * Set cutoff distance.
     *
     * @return cutoff distance
     *
     */
    void setCutoffDistance(RealOpenMM cutoffDistance);

    /**
     * Get alpha used in Ewald summation.
     *
     * @return alpha
     *
     */
    RealOpenMM getAlphaEwald() const;

    /**
     * Set alpha used in Ewald summation.
     *
     * @return alpha
     *
     */
    void setAlphaEwald(RealOpenMM alphaEwald);

    /**
     * Get PME grid dimensions.
     *
     * @param pmeGridDimensions contains PME grid dimensions upon return

     *
     */
    void getPmeGridDimensions(std::vector<int>& pmeGridDimensions) const;

    /**
     * Set PME grid dimensions.
     *
     * @param pmeGridDimensions input PME grid dimensions 
     *
     */
    void setPmeGridDimensions(std::vector<int>& pmeGridDimensions);

    /**
     * Set periodic box size.
     *
     * @param vectors    the vectors defining the periodic box
     */
     void setPeriodicBoxSize(OpenMM::RealVec* vectors);

private:

    static const int AMOEBA_PME_ORDER;
    static const RealOpenMM SQRT_PI;

    RealOpenMM _alphaEwald;
    RealOpenMM _cutoffDistance;
    RealOpenMM _cutoffDistanceSquared;

    RealVec _recipBoxVectors[3];
    RealVec _periodicBoxVectors[3];

    int _totalGridSize;
    IntVec _pmeGridDimensions;

    fftpack_t   _fftplan;

    unsigned int _pmeGridSize;
    t_complex* _pmeGrid;
 
    std::vector<RealOpenMM> _pmeBsplineModuli[3];
    std::vector<RealOpenMM4> _thetai[3];
    std::vector<IntVec> _iGrid;
    std::vector<RealOpenMM> _phi;
    std::vector<RealOpenMM> _phid;
    std::vector<RealOpenMM> _phip;
    std::vector<RealOpenMM> _phidp;
    std::vector<RealOpenMM4> _pmeBsplineTheta;
    std::vector<RealOpenMM4> _pmeBsplineDtheta;

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
    void getPeriodicDelta(RealVec& deltaR) const;

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
                                   RealOpenMM dscale, RealOpenMM pscale, RealOpenMM r,
                                   RealVec& dampedDInverseDistances, RealVec& dampedPInverseDistances) const;
    
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
                                             RealOpenMM dscale, RealOpenMM pscale);
    
    /**
     * Calculate fixed multipole fields.
     *
     * @param particleData vector particle data
     * 
     */
    void calculateFixedMultipoleField(const vector<MultipoleParticleData>& particleData);

    /**
     * This is called from computeAmoebaBsplines().  It calculates the spline coefficients for a single atom along a single axis.
     * 
     * @param thetai output spline coefficients
     * @param w offset from grid point
     */
    void computeBSplinePoint(std::vector<RealOpenMM4>& thetai, RealOpenMM w);
    
    /**
     * Compute bspline coefficients.
     *
     * @param particleData   vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void computeAmoebaBsplines(const std::vector<MultipoleParticleData>& particleData);

    /**
     * Transform multipoles from cartesian coordinates to fractional coordinates.
     */
    void transformMultipolesToFractionalCoordinates(const vector<MultipoleParticleData>& particleData);

    /**
     * Transform potential from fractional coordinates to cartesian coordinates.
     */
    void transformPotentialToCartesianCoordinates(const std::vector<RealOpenMM>& fphi, std::vector<RealOpenMM>& cphi) const;

    /**
     * Spread fixed multipoles onto PME grid.
     * 
     * @param particleData vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     */
    void spreadFixedMultipolesOntoGrid(const vector<MultipoleParticleData>& particleData);

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
    RealOpenMM computeReciprocalSpaceFixedMultipoleForceAndEnergy(const std::vector<MultipoleParticleData>& particleData,
                                                                  std::vector<RealVec>& forces, std::vector<RealVec>& torques) const;

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
                                             RealOpenMM preFactor1, RealOpenMM preFactor2, const RealVec& delta,
                                             const std::vector<RealVec>& inducedDipole,
                                             std::vector<RealVec>& field) const;

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
    void spreadInducedDipolesOnGrid(const std::vector<RealVec>& inputInducedDipole,
                                    const std::vector<RealVec>& inputInducedDipolePolar);

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
    void recordInducedDipoleField(vector<RealVec>& field, vector<RealVec>& fieldPolar);

    /**
     * Compute Pme self energy.
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     */
    RealOpenMM calculatePmeSelfEnergy(const std::vector<MultipoleParticleData>& particleData) const;

    /**
     * Compute the self torques.
     *
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 vector of torques
     */
    void calculatePmeSelfTorque(const std::vector<MultipoleParticleData>& particleData, std::vector<RealVec>& torques) const;

    /**
     * Calculate direct space electrostatic interaction between particles I and J.
     * 
     * @param particleI         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle I
     * @param particleJ         positions and parameters (charge, labFrame dipoles, quadrupoles, ...) for particle J
     * @param scalingFactors    scaling factors for interaction
     * @param forces            vector of particle forces to be updated
     * @param torques           vector of particle torques to be updated
     */
    RealOpenMM calculatePmeDirectElectrostaticPairIxn(const MultipoleParticleData& particleI, const MultipoleParticleData& particleJ,
                                                      const std::vector<RealOpenMM>& scalingFactors,
                                                      std::vector<RealVec>& forces, std::vector<RealVec>& torques) const;

    /**
     * Calculate reciprocal space energy/force/torque for dipole interaction.
     * 
     * @param polarizationType  if 'Direct' polariztion, only initial induced dipoles calculated
     *                          if 'Mutual' polariztion, induced dipoles converged to specified tolerance
     * @param particleData      vector of particle positions and parameters (charge, labFrame dipoles, quadrupoles, ...)
     * @param forces            vector of particle forces to be updated
     * @param torques           vector of particle torques to be updated
     */
     RealOpenMM computeReciprocalSpaceInducedDipoleForceAndEnergy(AmoebaReferenceMultipoleForce::PolarizationType polarizationType,
                                                                  const std::vector<MultipoleParticleData>& particleData,
                                                                  std::vector<RealVec>& forces, std::vector<RealVec>& torques) const;

    /**
     * Calculate electrostatic forces.
     * 
     * @param particleData            vector of parameters (charge, labFrame dipoles, quadrupoles, ...) for particles
     * @param torques                 output torques
     * @param forces                  output forces 
     *
     * @return energy
     */
    RealOpenMM calculateElectrostatic(const std::vector<MultipoleParticleData>& particleData, 
                                      std::vector<OpenMM::RealVec>& torques,
                                      std::vector<OpenMM::RealVec>& forces);

};

} // namespace OpenMM

#endif // _AmoebaReferenceMultipoleForce___
