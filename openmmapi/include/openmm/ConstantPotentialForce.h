#ifndef OPENMM_CONSTANTPOTENTIALFORCE_H_
#define OPENMM_CONSTANTPOTENTIALFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Evan Pretti                                        *
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

#include "Context.h"
#include "Force.h"
#include <map>
#include <set>
#include <utility>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements the periodic finite field constant potential method.
 * This is described in Dufils et al., Phys. Rev. Lett. 123, 195501 (2019).
 * This class implements a standard Coulomb force using the particle mesh Ewald
 * method, but it can also solve for the magnitudes of Gaussian charge
 * distributions on conducting electrodes held at fixed potentials and apply
 * external electric fields.  An implementation of a semiclassical Thomas-Fermi
 * model described in Scalfi et al., J. Chem. Phys. 153, 174704 (2020) is also
 * included.  Unlike NonbondedForce, Lennard-Jones interactions are not computed
 * by this force, and should be specified in another force if they are desired.
 *
 * To use this class, create a ConstantPotentialForce object, then call
 * addParticle() once for each particle in the System to define its parameters.
 * The number of particles for which you define nonbonded parameters must be
 * exactly equal to the number of particles in the System, or else an exception
 * will be thrown when you try to create a Context.  After a particle has been
 * added, you can modify its force field parameters by calling
 * setParticleParameters().  This will have no effect on Contexts that already
 * exist unless you call updateParametersInContext().
 *
 * ConstantPotentialForce also lets you specify "exceptions", particular pairs
 * of particles whose interactions should be computed based on different
 * parameters than those defined for the individual particles.  This can be used
 * to completely exclude certain interactions from the force calculation, or to
 * alter how they interact with each other.
 *
 * Many molecular force fields omit Coulomb and Lennard-Jones interactions
 * between particles separated by one or two bonds, while using modified
 * parameters for those separated by three bonds (known as "1-4 interactions").
 * This class provides a convenience method for this case called
 * createExceptionsFromBonds().  You pass to it a list of bonds and the scale
 * factors to use for 1-4 interactions.  It identifies all pairs of particles
 * which are separated by 1, 2, or 3 bonds, then automatically creates
 * exceptions for them.
 *
 * To treat a group of particles as an electrode (such that the charges
 * specified for them using addParticle() or setParticleParameters() will be
 * ignored, and charges will be solved for over the course of the simulation
 * instead), call addElectrode() with a set of particle indices.  After an
 * electrode has been added, you can modify its parameters by calling
 * setElectrodeParameters().  A constraint on the total charge of the system can
 * be enabled with setUseChargeConstraint() and setChargeConstraintTarget(), and
 * an external field can be applied by using setExternalField() to specify a
 * non-zero electric field strength.  Once a Context has been created, calling
 * getCharges() with the Context will return a vector of the actual current
 * charges on each particle, including the solved fluctuating charges on the
 * electrode particles.
 */

class OPENMM_EXPORT ConstantPotentialForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for
     * solving for electrode charges.
     */
    enum ConstantPotentialMethod {
        /**
         * A conjugate gradient method is used to iteratively solve for the
         * electrode charges at each simulation step.
         */
        CG = 0,
        /**
         * A capacitance matrix is precomputed at the start of the simulation
         * and is used to directly calculate the electrode charges at each
         * simulation step.  This method can only be used if all electrode
         * particles have fixed positions and the periodic box does not change.
         */
        Matrix = 1,
    };
    /**
     * Create a ConstantPotentialForce.
     */
    ConstantPotentialForce();
    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of special interactions that should be calculated differently from other interactions.
     */
    int getNumExceptions() const {
        return exceptions.size();
    }
    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the reciprocal space cutoff and separation
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     *
     * For PME calculations, if setPMEParameters() is used to set alpha to something other than 0,
     * this value is ignored.
     */
    double getEwaldErrorTolerance() const;
    /**
     * Set the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the reciprocal space cutoff and separation
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     *
     * For PME calculations, if setPMEParameters() is used to set alpha to something other than 0,
     * this value is ignored.
     */
    void setEwaldErrorTolerance(double tol);
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
     * Get the parameters being used for PME in a particular Context.  Because some platforms have restrictions
     * on the allowed grid sizes, the values that are actually used may be slightly different from those
     * specified with setPMEParameters(), or the standard values calculated based on the Ewald error tolerance.
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
     * Add the nonbonded force parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param charge    the charge of the particle, measured in units of the proton charge
     * @return the index of the particle that was added
     */
    int addParticle(double charge);
    /**
     * Get the nonbonded force parameters for a particle.
     *
     * @param index          the index of the particle for which to get parameters
     * @param[out] charge    the charge of the particle, measured in units of the proton charge
     */
    void getParticleParameters(int index, double& charge) const;
    /**
     * Set the nonbonded force parameters for a particle.
     *
     * @param index     the index of the particle for which to set parameters
     * @param charge    the charge of the particle, measured in units of the proton charge
     */
    void setParticleParameters(int index, double charge);
    /**
     * Add an interaction to the list of exceptions that should be calculated
     * differently from other interactions.  If chargeProd is equal to 0, this
     * will cause the interaction to be completely omitted from force and energy
     * calculations.
     * 
     * Cutoffs are never applied to exceptions.  That is because they are
     * primarily used for 1-4 interactions, which are really a type of bonded
     * interaction and are parametrized together with the other bonded
     * interactions.
     *
     * In many cases, you can use createExceptionsFromBonds() rather than adding
     * each exception explicitly.
     *
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
     * @param replace    determines the behavior if there is already an exception for the same two particles.  If true, the existing one is replaced.  If false,
     *                   an exception is thrown.
     * @return the index of the exception that was added
     */
    int addException(int particle1, int particle2, double chargeProd, bool replace = false);
    /**
     * Get the force field parameters for an interaction that should be calculated differently from others.
     *
     * @param index           the index of the interaction for which to get parameters
     * @param[out] particle1  the index of the first particle involved in the interaction
     * @param[out] particle2  the index of the second particle involved in the interaction
     * @param[out] chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
     */
    void getExceptionParameters(int index, int& particle1, int& particle2, double& chargeProd) const;
    /**
     * Set the force field parameters for an interaction that should be
     * calculated differently from others.  If chargeProd is equal to 0, this
     * will cause the interaction to be completely omitted from force and energy
     * calculations.
     * 
     * Cutoffs are never applied to exceptions.  That is because they are
     * primarily used for 1-4 interactions, which are really a type of bonded
     * interaction and are parametrized together with the other bonded
     * interactions.
     *
     * @param index      the index of the interaction for which to get parameters
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
     */
    void setExceptionParameters(int index, int particle1, int particle2, double chargeProd);
    /**
     * Identify exceptions based on the molecular topology.  Particles which are separated by one or two bonds are set
     * to not interact at all, while pairs of particles separated by three bonds (known as "1-4 interactions") have
     * their Coulomb interactions reduced by a fixed factor.
     *
     * @param bonds           the set of bonds based on which to construct exceptions.  Each element specifies the indices of
     *                        two particles that are bonded to each other.
     * @param coulomb14Scale  pairs of particles separated by three bonds will have the strength of their Coulomb interaction
     *                        multiplied by this factor
     */
    void createExceptionsFromBonds(const std::vector<std::pair<int, int> >& bonds, double coulomb14Scale);
    /**
     * Update particle, exception, and electrode parameters in a Context to
     * match those stored in this Force object.  This method provides an
     * efficient method to update certain parameters in an existing Context
     * without needing to reinitialize it.  Simply call setParticleParameters()
     * and setExceptionParameters() to modify this object's parameters, then
     * call updateParametersInContext() to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is
     * the parameters of particles, exceptions, and electrodes, as well as the
     * target total charge for the charge constraint.  All other aspects of the
     * Force (the constant potential method, the cutoff distance, etc.) are
     * unaffected and can only be changed by reinitializing the Context.
     * Furthermore, only the chargeProd value of an exception can be changed;
     * the pair of particles involved in the exception cannot change.
     * Similarly, for electrodes, the set of particles involved in the electrode
     * cannot be updated with this method.  Finally, this method cannot be used
     * to add new particles, exceptions, or electrodes, only to change the
     * parameters of existing ones.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return true;
    }
    /**
     * Get whether periodic boundary conditions should be applied to exceptions.
     * Usually this is not appropriate, because exceptions are normally used to
     * represent bonded interactions (1-2, 1-3, and 1-4 pairs), but there are
     * situations when it does make sense.  For example, you may want to
     * simulate an infinite chain where one end of a molecule is bonded to the
     * opposite end of the next periodic copy.  Note that cutoffs are never
     * applied to exceptions, again because they are normally used to represent
     * bonded interactions.
     */
    bool getExceptionsUsePeriodicBoundaryConditions() const;
    /**
     * Set whether periodic boundary conditions should be applied to exceptions.
     * Usually this is not appropriate, because exceptions are normally used to
     * represent bonded interactions (1-2, 1-3, and 1-4 pairs), but there are
     * situations when it does make sense.  For example, you may want to
     * simulate an infinite chain where one end of a molecule is bonded to the
     * opposite end of the next periodic copy.  Note that cutoffs are never
     * applied to exceptions, again because they are normally used to represent
     * bonded interactions.
     */
    void setExceptionsUsePeriodicBoundaryConditions(bool periodic);
    /**
     * Get the method used for calculating electrode charges.
     */
    ConstantPotentialMethod getConstantPotentialMethod() const;
    /**
     * Set the method used for calculating electrode charges.
     */
    void setConstantPotentialMethod(ConstantPotentialMethod method);
    /**
     * Get whether or not to use a preconditioner when solving for electrode
     * charges with the conjugate gradient method.  This option has no effect
     * when the matrix solver is in use.
     */
    bool getUsePreconditioner() const;
    /**
     * Set whether or not to use a preconditioner when solving for electrode
     * charges with the conjugate gradient method.  This option has no effect
     * when the matrix solver is in use.
     */
    void setUsePreconditioner(bool use);
    /**
     * Get the tolerance, in units of kJ/mol per proton charge, used for the
     * conjugate gradient method of calculating electrode charges.  The method
     * will iterate until the RMS error in the electrode potentials (gradient of
     * the energy with respect to the electrode particle charges) no longer
     * exceeds this tolerance.
     */
    double getCGErrorTolerance() const;
    /**
     * Set the tolerance, in units of kJ/mol per proton charge, used for the
     * conjugate gradient method of calculating electrode charges.  The method
     * will iterate until the RMS error in the electrode potentials (gradient of
     * the energy with respect to the electrode particle charges) no longer
     * exceeds this tolerance.
     */
    void setCGErrorTolerance(double tol);
    /**
     * Get the number of electrodes that have been added.
     */
    int getNumElectrodes() const {
        return electrodes.size();
    }
    /**
     * Add a new electrode from a set of particles.  The specified particles
     * will have their charges solved for such that they are held at the
     * specified electric potential.  Particles in a system may belong to at
     * most a single electrode.
     *
     * @param electrodeParticles  the indices of the particles to be included in this electrode
     * @param potential           the electric potential of the particles in this electrode, measured in kJ/mol per proton charge
     * @param gaussianWidth       the width of the Gaussian charge distribution assigned to particles in this electrode, measured in nm
     * @param thomasFermiScale    the square of the Thomas-Fermi length divided by the Thomas-Fermi Voronoi volume, measured in 1/nm
     * @return                    the index of the electrode added
     */
    int addElectrode(const std::set<int>& electrodeParticles, double potential, double gaussianWidth, double thomasFermiScale);
    /**
     * Get the parameters for an electrode.
     *
     * @param index                    the index of the electrode for which to get parameters
     * @param[out] electrodeParticles  the indices of the particles to be included in this electrode
     * @param[out] potential           the electric potential of the particles in this electrode, measured in kJ/mol per proton charge
     * @param[out] gaussianWidth       the width of the Gaussian charge distribution assigned to particles in this electrode, measured in nm
     * @param[out] thomasFermiScale    the square of the Thomas-Fermi length divided by the Thomas-Fermi Voronoi volume, measured in 1/nm
     */
    void getElectrodeParameters(int index, std::set<int>& electrodeParticles, double& potential, double& gaussianWidth, double& thomasFermiScale) const;
    /**
     * Set the parameters for an electrode.
     *
     * @param index               the index of the electrode for which to set parameters
     * @param electrodeParticles  the indices of the particles to be included in this electrode
     * @param potential           the electric potential of the particles in this electrode, measured in kJ/mol per proton charge
     * @param gaussianWidth       the width of the Gaussian charge distribution assigned to particles in this electrode, measured in nm
     * @param thomasFermiScale    the square of the Thomas-Fermi length divided by the Thomas-Fermi Voronoi volume, measured in 1/nm
     */
    void setElectrodeParameters(int index, const std::set<int>& electrodeParticles, double potential, double gaussianWidth, double thomasFermiScale);
    /**
     * Get whether or not to apply a constraint to hold the total charge of the
     * system constant.
     */
    bool getUseChargeConstraint() const;
    /**
     * Set whether or not to apply a constraint to hold the total charge of the
     * system constant.
     */
    void setUseChargeConstraint(bool use);
    /**
     * Get the desired charge, in units of the proton charge, at which to hold
     * the system if the total charge constraint is active.  This includes the
     * (fluctuating) charges on all electrode particles as well as the (fixed)
     * charges on all non-electrode particles.
     */
    double getChargeConstraintTarget() const;
    /**
     * Get the desired charge, in units of the proton charge, at which to hold
     * the system if the total charge constraint is active.  This includes the
     * (fluctuating) charges on all electrode particles as well as the (fixed)
     * charges on all non-electrode particles.
     */
    void setChargeConstraintTarget(double charge);
    /**
     * Get the external electric field strength, measured in kJ/mol/nm per
     * proton charge.
     */
    void getExternalField(Vec3& field) const;
    /**
     * Set the external electric field strength, measured in kJ/mol/nm per
     * proton charge.
     */
    void setExternalField(const Vec3& field);
    /**
     * Get the charges on all particles: for non-electrode particles, these will
     * simply be the fixed charges set, while for electrode particles, they will
     * be the current charges solved for by the constant potential method.
     */
    void getCharges(Context& context, std::vector<double>& charges) const;
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    class ExceptionInfo;
    class ElectrodeInfo;
    ConstantPotentialMethod constantPotentialMethod;
    double cutoffDistance, ewaldErrorTol, alpha, cgErrorTol, chargeTarget;
    Vec3 externalField;
    bool exceptionsUsePeriodic, useChargeConstraint, usePreconditioner;
    int nx, ny, nz;
    void addExclusionsToSet(const std::vector<std::set<int> >& bonded12, std::set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const;
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
    std::vector<ElectrodeInfo> electrodes;
    std::map<std::pair<int, int>, int> exceptionMap;
    mutable int numContexts, firstChangedParticle, lastChangedParticle, firstChangedException, lastChangedException, firstChangedElectrode, lastChangedElectrode;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class ConstantPotentialForce::ParticleInfo {
public:
    double charge;
    ParticleInfo() {
        charge = 0.0;
    }
    ParticleInfo(double charge) :
        charge(charge) {
    }
};

/**
 * This is an internal class used to record information about an exception.
 * @private
 */
class ConstantPotentialForce::ExceptionInfo {
public:
    int particle1, particle2;
    double chargeProd;
    ExceptionInfo() {
        particle1 = particle2 = -1;
        chargeProd = 0.0;
    }
    ExceptionInfo(int particle1, int particle2, double chargeProd) :
        particle1(particle1), particle2(particle2), chargeProd(chargeProd) {
    }
};

/**
 * This is an internal class used to record information about an electrode.
 * @private
 */
class ConstantPotentialForce::ElectrodeInfo {
public:
    std::set<int> particles;
    double potential;
    double gaussianWidth;
    double thomasFermiScale;
    ElectrodeInfo() {
        potential = 0.0;
        gaussianWidth = 0.0;
        thomasFermiScale = 0.0;
    }
    ElectrodeInfo(const std::set<int>& particles, double potential, double gaussianWidth, double thomasFermiScale) :
        particles(particles), potential(potential), gaussianWidth(gaussianWidth), thomasFermiScale(thomasFermiScale) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_CONSTANTPOTENTIALFORCE_H_*/
