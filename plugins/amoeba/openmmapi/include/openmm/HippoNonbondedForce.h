#ifndef OPENMM_HIPPO_NONBONDED_FORCE_H_
#define OPENMM_HIPPO_NONBONDED_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2019 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Mark Friedrichs                                    *
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

#include "openmm/Force.h"
#include "internal/windowsExportAmoeba.h"
#include "openmm/Vec3.h"
#include <map>
#include <utility>
#include <vector>

namespace OpenMM {

/**
 * This class implements all nonbonded interactions in the HIPPO force field: electrostatics,
 * induction, charge transfer, dispersion, and repulsion.  Although some of these are
 * conceptually distinct, they share parameters in common and are most efficiently computed
 * together.  For example, the same multipole definitions are used for both electrostatics
 * and Pauli repulsion.  Therefore, all of them are computed by a single Force object.
 *
 * To use it, create a HippoNonbondedForce object, then call addParticle() once for each particle.  After
 * an entry has been added, you can modify its force field parameters by calling setParticleParameters().
 * This will have no effect on Contexts that already exist unless you call updateParametersInContext().
 * 
 * You also can specify "exceptions", particular pairs of particles whose interactions should be
 * reduced or completely omitted.  Call addException() to define exceptions.
 */

class OPENMM_EXPORT_AMOEBA HippoNonbondedForce : public Force {
public:

    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        PME = 1
    };

    enum ParticleAxisTypes { ZThenX = 0, Bisector = 1, ZBisect = 2, ThreeFold = 3, ZOnly = 4, NoAxisType = 5 };

    /**
     * Create a HippoNonbondedForce.
     */
    HippoNonbondedForce();
    /**
     * Get the number of particles in the potential function.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of exceptions.
     */
    int getNumExceptions() const {
        return exceptions.size();
    }
    /**
     * Get the method used for handling long-range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long-range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);
    /**
     * Get the distance at which the switching function begins to reduce the repulsion and charge transfer interactions.  This must be
     * less than the cutoff distance.
     */
    double getSwitchingDistance() const;
    /**
     * Set the distance at which the switching function begins to reduce the repulsion and charge transfer interactions.  This must be
     * less than the cutoff distance.
     */
    void setSwitchingDistance(double distance);
    /**
     * Get the coefficients for the mu_0, mu_1, mu_2, ..., mu_n terms in the extrapolation
     * algorithm for induced dipoles.
     */
    const std::vector<double>& getExtrapolationCoefficients() const;
    /**
     * Set the coefficients for the mu_0, mu_1, mu_2, ..., mu_n terms in the extrapolation
     * algorithm for induced dipoles.
     *
     * @param coefficients      a vector whose mth entry specifies the coefficient for mu_m.  The length of this
     *                          vector determines how many iterations are performed.
     */
    void setExtrapolationCoefficients(const std::vector<double> &coefficients);
    /**
     * Get the parameters to use for PME calculations.  If alpha is 0 (the default), these parameters are
     * ignored and instead their values are chosen based on the Ewald error tolerance.
     *
     * @param[out] alpha   the separation parameter
     * @param[out] nx      the number of grid points along the X axis
     * @param[out] ny      the number of grid points along the Y axis
     * @param[out] nz      the number of grid points along the Z axis
     */
    void getPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Get the parameters to use for dispersion PME calculations.  If alpha is 0 (the default),
     * these parameters are ignored and instead their values are chosen based on the Ewald error tolerance.
     *
     * @param[out] alpha   the separation parameter
     * @param[out] nx      the number of dispersion grid points along the X axis
     * @param[out] ny      the number of dispersion grid points along the Y axis
     * @param[out] nz      the number of dispersion grid points along the Z axis
     */
    void getDPMEParameters(double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Set the parameters to use for PME calculations.  If alpha is 0 (the default), these parameters are
     * ignored and instead their values are chosen based on the Ewald error tolerance.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void setPMEParameters(double alpha, int nx, int ny, int nz);
    /**
     * Set the parameters to use for dispersion PME calculations.  If alpha is 0 (the default),
     * these parameters are ignored and instead their values are chosen based on the Ewald error tolerance.
     *
     * @param alpha   the separation parameter
     * @param nx      the number of grid points along the X axis
     * @param ny      the number of grid points along the Y axis
     * @param nz      the number of grid points along the Z axis
     */
    void setDPMEParameters(double alpha, int nx, int ny, int nz);
    /**
     * Get the parameters being used for PME in a particular Context.  Because some platforms have restrictions
     * on the allowed grid sizes, the values that are actually used may be slightly different from those
     * specified with setPmeGridDimensions(), or the standard values calculated based on the Ewald error tolerance.
     * See the manual for details.
     *
     * @param context      the Context for which to get the parameters
     * @param[out] alpha   the separation parameter
     * @param[out] nx      the number of grid points along the X axis
     * @param[out] ny      the number of grid points along the Y axis
     * @param[out] nz      the number of grid points along the Z axis
     */
    void getPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Get the parameters being used for dispersion PME in a particular Context.  Because some
     * platforms have restrictions on the allowed grid sizes, the values that are actually used may be slightly different
     * from those specified with setPMEParameters(), or the standard values calculated based on the Ewald error tolerance.
     * See the manual for details.
     *
     * @param context      the Context for which to get the parameters
     * @param[out] alpha   the separation parameter
     * @param[out] nx      the number of grid points along the X axis
     * @param[out] ny      the number of grid points along the Y axis
     * @param[out] nz      the number of grid points along the Z axis
     */
    void getDPMEParametersInContext(const Context& context, double& alpha, int& nx, int& ny, int& nz) const;
    /**
     * Add the nonbonded force parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param charge            the particle's charge
     * @param dipole            the particle's molecular dipole (vector of size 3)
     * @param quadrupole        the particle's molecular quadrupole (vector of size 9)
     * @param coreCharge        the charge of the atomic core
     * @param alpha             controls the width of the particle's electron density
     * @param epsilon           sets the magnitude of charge transfer
     * @param damping           sets the length scale for charge transfer
     * @param c6                the coefficient of the dispersion interaction
     * @param pauliK            the coefficient of the Pauli repulsion interaction
     * @param pauliQ            the charge used in computing the Pauli repulsion interaction
     * @param pauliAlpha        the width of the particle's electron density for computing the Pauli repulsion interaction
     * @param polarizability    atomic polarizability
     * @param axisType          the particle's axis type
     * @param multipoleAtomZ    index of first atom used in defining the local coordinate system for multipoles
     * @param multipoleAtomX    index of second atom used in defining the local coordinate system for multipoles
     * @param multipoleAtomY    index of third atom used in defining the local coordinate system for multipoles
     * 
     * @return the index of the particle that was added
     */
    int addParticle(double charge, const std::vector<double>& dipole, const std::vector<double>& quadrupole, double coreCharge,
                    double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                    double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY);
    /**
     * Get the nonbonded force parameters for a particle.
     *
     * @param index             the index of the particle for which to get parameters
     * @param charge            the particle's charge
     * @param dipole            the particle's molecular dipole (vector of size 3)
     * @param quadrupole        the particle's molecular quadrupole (vector of size 9)
     * @param coreCharge        the charge of the atomic core
     * @param alpha             controls the width of the particle's electron density
     * @param epsilon           sets the magnitude of charge transfer
     * @param damping           sets the length scale for charge transfer
     * @param c6                the coefficient of the dispersion interaction
     * @param pauliK            the coefficient of the Pauli repulsion interaction
     * @param pauliQ            the charge used in computing the Pauli repulsion interaction
     * @param pauliAlpha        the width of the particle's electron density for computing the Pauli repulsion interaction
     * @param polarizability    atomic polarizability
     * @param axisType          the particle's axis type
     * @param multipoleAtomZ    index of first atom used in defining the local coordinate system for multipoles
     * @param multipoleAtomX    index of second atom used in defining the local coordinate system for multipoles
     * @param multipoleAtomY    index of third atom used in defining the local coordinate system for multipoles
     */
    void getParticleParameters(int index, double& charge, std::vector<double>& dipole, std::vector<double>& quadrupole, double& coreCharge,
                               double& alpha, double& epsilon, double& damping, double& c6, double& pauliK, double& pauliQ, double& pauliAlpha,
                               double& polarizability, int& axisType, int& multipoleAtomZ, int& multipoleAtomX, int& multipoleAtomY) const;
    /**
     * Set the nonbonded force parameters for a particle.
     *
     * @param index             the index of the particle for which to set parameters
     * @param charge            the particle's charge
     * @param dipole            the particle's molecular dipole (vector of size 3)
     * @param quadrupole        the particle's molecular quadrupole (vector of size 9)
     * @param coreCharge        the charge of the atomic core
     * @param alpha             controls the width of the particle's electron density
     * @param epsilon           sets the magnitude of charge transfer
     * @param damping           sets the length scale for charge transfer
     * @param c6                the coefficient of the dispersion interaction
     * @param pauliK            the coefficient of the Pauli repulsion interaction
     * @param pauliQ            the charge used in computing the Pauli repulsion interaction
     * @param pauliAlpha        the width of the particle's electron density for computing the Pauli repulsion interaction
     * @param polarizability    atomic polarizability
     * @param axisType          the particle's axis type
     * @param multipoleAtomZ    index of first atom used in defining the local coordinate system for multipoles
     * @param multipoleAtomX    index of second atom used in defining the local coordinate system for multipoles
     * @param multipoleAtomY    index of third atom used in defining the local coordinate system for multipoles
     */
    void setParticleParameters(int index, double charge, const std::vector<double>& dipole, const std::vector<double>& quadrupole, double coreCharge,
                               double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                               double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY);
    /**
     * Add an interaction to the list of exceptions that should be calculated differently from other interactions.
     * If all scale factors are set to 0, this will cause the interaction to be completely omitted from
     * force and energy calculations.
     *
     * @param particle1                  the index of the first particle involved in the interaction
     * @param particle2                  the index of the second particle involved in the interaction
     * @param multipoleMultipoleScale    the factor by which to scale the Coulomb interaction between fixed multipoles
     * @param dipoleMultipoleScale       the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
     * @param dipoleDipoleScale          the factor by which to scale the Coulomb interaction between induced dipoles
     * @param dispersionScale            the factor by which to scale the dispersion interaction
     * @param repulsionScale             the factor by which to scale the Pauli repulsion
     * @param chargeTransferScale        the factor by which to scale the charge transfer interaction
     * @param replace                    determines the behavior if there is already an exception for the same two particles.  If true, the existing one is replaced.
     *                                   If false, an exception is thrown.
     * @return the index of the exception that was added
     */
    int addException(int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale,
                     double dispersionScale, double repulsionScale, double chargeTransferScale, bool replace = false);
    /**
     * Get the scale factors for an interaction that should be calculated differently from others.
     *
     * @param index                      the index of the interaction for which to get parameters
     * @param particle1                  the index of the first particle involved in the interaction
     * @param particle2                  the index of the second particle involved in the interaction
     * @param multipoleMultipoleScale    the factor by which to scale the Coulomb interaction between fixed multipoles
     * @param dipoleMultipoleScale       the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
     * @param dipoleDipoleScale          the factor by which to scale the Coulomb interaction between induced dipoles
     * @param dispersionScale            the factor by which to scale the dispersion interaction
     * @param repulsionScale             the factor by which to scale the Pauli repulsion
     * @param chargeTransferScale        the factor by which to scale the charge transfer interaction
     */
    void getExceptionParameters(int index, int& particle1, int& particle2, double& multipoleMultipoleScale, double& dipoleMultipoleScale, double& dipoleDipoleScale,
                                double& dispersionScale, double& repulsionScale, double& chargeTransferScale) const;
    /**
     * Set the scale factors for an interaction that should be calculated differently from others.
     * If all scale factors are set to 0, this will cause the interaction to be completely omitted from
     * force and energy calculations.
     *
     * @param index                      the index of the interaction for which to set parameters
     * @param particle1                  the index of the first particle involved in the interaction
     * @param particle2                  the index of the second particle involved in the interaction
     * @param multipoleMultipoleScale    the factor by which to scale the Coulomb interaction between fixed multipoles
     * @param dipoleMultipoleScale       the factor by which to scale the Coulomb interaction between an induced dipole and a fixed multipole
     * @param dipoleDipoleScale          the factor by which to scale the Coulomb interaction between induced dipoles
     * @param dispersionScale            the factor by which to scale the dispersion interaction
     * @param repulsionScale             the factor by which to scale the Pauli repulsion
     * @param chargeTransferScale        the factor by which to scale the charge transfer interaction
     */
    void setExceptionParameters(int index, int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale,
                                double dispersionScale, double repulsionScale, double chargeTransferScale);
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     *
     * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
     */
    double getEwaldErrorTolerance() const;
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the grid dimensions and separation (alpha)
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     *
     * This can be overridden by explicitly setting an alpha parameter and grid dimensions to use.
     */
    void setEwaldErrorTolerance(double tol);
    /**
     * Get the fixed dipole moments of all particles in the global reference frame.
     *
     * @param context         the Context for which to get the fixed dipoles
     * @param[out] dipoles    the fixed dipole moment of particle i is stored into the i'th element
     */
    void getLabFramePermanentDipoles(Context& context, std::vector<Vec3>& dipoles);
    /**
     * Get the induced dipole moments of all particles.
     *
     * @param context         the Context for which to get the induced dipoles
     * @param[out] dipoles    the induced dipole moment of particle i is stored into the i'th element
     */
    void getInducedDipoles(Context& context, std::vector<Vec3>& dipoles);
    /**
     * Update the particle and exception parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext() to
     * copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the parameters of particles and exceptions.
     * All other aspects of the Force (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be
     * changed by reinitializing the Context.  Furthermore, only the scale factors for an exception  can be changed; the
     * pair of particles involved in the exception cannot change.  Finally, this method cannot be used to add new
     * particles or exceptions, only to change the parameters of existing ones.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if nonbondedMethod uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == HippoNonbondedForce::PME;
    }
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    class ExceptionInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance, switchingDistance;
    double ewaldErrorTol;
    double alpha, dalpha;
    int nx, ny, nz, dnx, dny, dnz;
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
    std::map<std::pair<int, int>, int> exceptionMap;
    std::vector<double> extrapolationCoefficients;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class HippoNonbondedForce::ParticleInfo {
public:
    int axisType, multipoleAtomZ, multipoleAtomX, multipoleAtomY;
    double charge, coreCharge, alpha, epsilon, damping, c6, pauliK, pauliQ, pauliAlpha, polarizability;
    std::vector<double> dipole, quadrupole;
    ParticleInfo() {
        axisType = multipoleAtomZ = multipoleAtomX = multipoleAtomY = -1;
        charge = coreCharge = alpha = epsilon = damping = c6 = pauliK = pauliQ = pauliAlpha = polarizability = 0.0;
        dipole.resize(3);
        quadrupole.resize(9);
    }
    ParticleInfo(double charge, const std::vector<double>& dipole, const std::vector<double>& quadrupole, double coreCharge,
                 double alpha, double epsilon, double damping, double c6, double pauliK, double pauliQ, double pauliAlpha,
                 double polarizability, int axisType, int multipoleAtomZ, int multipoleAtomX, int multipoleAtomY) :
        charge(charge), dipole(dipole), quadrupole(quadrupole), coreCharge(coreCharge), alpha(alpha), epsilon(epsilon),
        damping(damping), c6(c6), pauliK(pauliK), pauliQ(pauliQ), pauliAlpha(pauliAlpha), polarizability(polarizability),
        axisType(axisType), multipoleAtomZ(multipoleAtomZ), multipoleAtomX(multipoleAtomX), multipoleAtomY(multipoleAtomY) {
    }
};

/**
 * This is an internal class used to record information about an exception.
 * @private
 */
class HippoNonbondedForce::ExceptionInfo {
public:
    int particle1, particle2;
    double multipoleMultipoleScale, dipoleMultipoleScale, dipoleDipoleScale, dispersionScale, repulsionScale, chargeTransferScale;
    ExceptionInfo() {
        particle1 = particle2 = -1;
        multipoleMultipoleScale = dipoleMultipoleScale = dipoleDipoleScale = dispersionScale = repulsionScale = chargeTransferScale = 0.0;
    }
    ExceptionInfo(int particle1, int particle2, double multipoleMultipoleScale, double dipoleMultipoleScale, double dipoleDipoleScale, double dispersionScale, double repulsionScale, double chargeTransferScale) :
        particle1(particle1), particle2(particle2), multipoleMultipoleScale(multipoleMultipoleScale), dipoleMultipoleScale(dipoleMultipoleScale),
        dipoleDipoleScale(dipoleDipoleScale), dispersionScale(dispersionScale), repulsionScale(repulsionScale), chargeTransferScale(chargeTransferScale) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_HIPPO_NONBONDED_FORCE_H_*/
