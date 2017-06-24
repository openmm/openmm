#ifndef OPENMM_GAYBERNEFORCE_H_
#define OPENMM_GAYBERNEFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2016 Stanford University and the Authors.           *
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

#include "Context.h"
#include "Force.h"
#include <map>
#include <utility>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements the Gay-Berne anisotropic potential.  This is similar to a Lennard-Jones potential,
 * but it represents the particles as ellipsoids rather than point particles.  In addition to the standard
 * sigma and epsilon parameters, each particle has three widths sx, sy, and sz that give the diameter of the
 * ellipsoid along each axis.  It also has three scale factors ex, ey, and ez that scale the strength
 * of the interaction along each axis.  You can think of this force as a Lennard-Jones interaction computed
 * based on the distance between the nearest points on two ellipsoids.  The scale factors act as multipliers
 * for epsilon along each axis, so the strength of the interaction along the ellipsoid's x axis is multiplied by
 * ex, and likewise for the other axes.  If two particles each have all their widths set to sigma and all their
 * scale factors set to 1, the interaction simplifies to a standard Lennard-Jones force between point particles.
 *
 * The orientation of a particle's ellipsoid is determined based on the positions of two other particles.
 * The vector to the first particle sets the direction of the x axis.  The vector to the second particle
 * (after subtracting out any x component) sets the direction of the y axis.  If the ellipsoid is axially
 * symmetric (sy=sz and ey=ez), you can omit the second particle and define only an x axis direction.
 * If the ellipsoid is a sphere (all three widths and all three scale factors are equal), both particles
 * can be omitted.
 *
 * To determine the values of sigma and epsilon for an interaction, this class uses Lorentz-Berthelot
 * combining rules: it takes the arithmetic mean of the sigmas and the geometric mean of the epsilons for
 * the two interacting particles.  You also can specify "exceptions", particular pairs of particles for
 * which different values should be used.
 *
 * To use this class, create a GayBerneForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define parameters must be exactly
 * equal to the number of particles in the System, or else an exception will be thrown when you try to
 * create a Context.  After a particle has been added, you can modify its force field parameters by calling
 * setParticleParameters().  This will have no effect on Contexts that already exist unless you call
 * updateParametersInContext().
 *
 * When using a cutoff, by default interactions are sharply truncated at the cutoff distance.  Optionally
 * you can instead use a switching function to make the interaction smoothly go to zero over a finite
 * distance range.  To enable this, call setUseSwitchingFunction().  You must also call setSwitchingDistance()
 * to specify the distance at which the interaction should begin to decrease.  The switching distance must be
 * less than the cutoff distance.
 */

class OPENMM_EXPORT GayBerneForce : public Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
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
    /**
     * Create a GayBerneForce.
     */
    GayBerneForce();
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
     * Get the method used for handling long range interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long range interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Get the cutoff distance (in nm) being used for interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);
    /**
     * Get whether a switching function is applied to the interaction.  If the nonbonded method is set
     * to NoCutoff, this option is ignored.
     */
    bool getUseSwitchingFunction() const;
    /**
     * Set whether a switching function is applied to the interaction.  If the nonbonded method is set
     * to NoCutoff, this option is ignored.
     */
    void setUseSwitchingFunction(bool use);
    /**
     * Get the distance at which the switching function begins to reduce the interaction.  This must be
     * less than the cutoff distance.
     */
    double getSwitchingDistance() const;
    /**
     * Set the distance at which the switching function begins to reduce the interaction.  This must be
     * less than the cutoff distance.
     */
    void setSwitchingDistance(double distance);
    /**
     * Add the parameters for a particle.  This should be called once for each particle in the System.
     * When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param sigma     the sigma parameter (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon   the epsilon parameter (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     * @param xparticle the index of the particle whose position defines the ellipsoid's x axis, or -1 if the ellipsoid is a sphere
     * @param yparticle the index of the particle whose position defines the ellipsoid's y axis, or -1 if the ellipsoid is axially symmetric
     * @param sx        the diameter of the ellipsoid along its x axis
     * @param sy        the diameter of the ellipsoid along its y axis
     * @param sz        the diameter of the ellipsoid along its z axis
     * @param ex        the factor by which epsilon is scaled along the ellipsoid's x axis
     * @param ey        the factor by which epsilon is scaled along the ellipsoid's y axis
     * @param ez        the factor by which epsilon is scaled along the ellipsoid's z axis
     * @return the index of the particle that was added
     */
    int addParticle(double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez);
    /**
     * Get the parameters for a particle.
     *
     * @param index          the index of the particle for which to get parameters
     * @param[out] sigma     the sigma parameter (corresponding to the van der Waals radius of the particle), measured in nm
     * @param[out] epsilon   the epsilon parameter (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     * @param[out] xparticle the index of the particle whose position defines the ellipsoid's x axis, or -1 if the ellipsoid is a sphere
     * @param[out] yparticle the index of the particle whose position defines the ellipsoid's y axis, or -1 if the ellipsoid is axially symmetric
     * @param[out] sx        the diameter of the ellipsoid along its x axis
     * @param[out] sy        the diameter of the ellipsoid along its y axis
     * @param[out] sz        the diameter of the ellipsoid along its z axis
     * @param[out] ex        the factor by which epsilon is scaled along the ellipsoid's x axis
     * @param[out] ey        the factor by which epsilon is scaled along the ellipsoid's y axis
     * @param[out] ez        the factor by which epsilon is scaled along the ellipsoid's z axis
     */
    void getParticleParameters(int index, double& sigma, double& epsilon, int& xparticle, int& yparticle, double& sx, double& sy, double& sz, double& ex, double& ey, double& ez) const;
    /**
     * Set the parameters for a particle.
     *
     * @param index     the index of the particle for which to set parameters
     * @param sigma     the sigma parameter (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon   the epsilon parameter (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     * @param xparticle the index of the particle whose position defines the ellipsoid's x axis, or -1 if the ellipsoid is a sphere
     * @param yparticle the index of the particle whose position defines the ellipsoid's y axis, or -1 if the ellipsoid is axially symmetric
     * @param sx        the diameter of the ellipsoid along its x axis
     * @param sy        the diameter of the ellipsoid along its y axis
     * @param sz        the diameter of the ellipsoid along its z axis
     * @param ex        the factor by which epsilon is scaled along the ellipsoid's x axis
     * @param ey        the factor by which epsilon is scaled along the ellipsoid's y axis
     * @param ez        the factor by which epsilon is scaled along the ellipsoid's z axis
     */
    void setParticleParameters(int index, double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez);
    /**
     * Add an interaction to the list of exceptions that should be calculated differently from other interactions.  If
     * epsilon is equal to 0, this will cause the interaction to be completely omitted from force and energy calculations.
     *
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param sigma      the sigma parameter (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon    the epsilon parameter (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     * @param replace    determines the behavior if there is already an exception for the same two particles.  If true, the existing one is replaced.  If false,
     *                   an exception is thrown.
     * @return the index of the exception that was added
     */
    int addException(int particle1, int particle2, double sigma, double epsilon, bool replace = false);
    /**
     * Get the force field parameters for an interaction that should be calculated differently from others.
     *
     * @param index           the index of the interaction for which to get parameters
     * @param[out] particle1  the index of the first particle involved in the interaction
     * @param[out] particle2  the index of the second particle involved in the interaction
     * @param[out] sigma      the sigma parameter (corresponding to the van der Waals radius of the particle), measured in nm
     * @param[out] epsilon    the epsilon parameter (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     */
    void getExceptionParameters(int index, int& particle1, int& particle2, double& sigma, double& epsilon) const;
    /**
     * Set the force field parameters for an interaction that should be calculated differently from others.  If
     * epsilon is equal to 0, this will cause the interaction to be completely omitted from force and energy calculations.
     *
     * @param index      the index of the interaction for which to get parameters
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param sigma      the sigma parameter (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon    the epsilon parameter (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     */
    void setExceptionParameters(int index, int particle1, int particle2, double sigma, double epsilon);
    /**
     * Update the particle and exception parameters in a Context to match those stored in this Force object.  This method
     * provides an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() and setExceptionParameters() to modify this object's parameters, then call
     * updateParametersInContext() to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the parameters of particles and exceptions.
     * All other aspects of the Force (the nonbonded method, the cutoff distance, etc.) are unaffected and can only be
     * changed by reinitializing the Context.  Furthermore, only the sigma and epsilon values of an exception can be
     * changed; the pair of particles involved in the exception cannot change.  Likewise, the xparticle and yparticle
     * defining the orientation of an ellipse cannot be changed.  Finally, this method cannot be used to add new
     * particles or exceptions, only to change the parameters of existing ones.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == GayBerneForce::CutoffPeriodic;
    }
protected:
    ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    class ExceptionInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance, switchingDistance;
    bool useSwitchingFunction;
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
    std::map<std::pair<int, int>, int> exceptionMap;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class GayBerneForce::ParticleInfo {
public:
    int xparticle, yparticle;
    double sigma, epsilon, sx, sy, sz, ex, ey, ez;
    ParticleInfo() {
        xparticle = yparticle = -1;
        sigma = epsilon = sx = sy = sz = ex = ey = ez = 0.0;
    }
    ParticleInfo(double sigma, double epsilon, int xparticle, int yparticle, double sx, double sy, double sz, double ex, double ey, double ez) :
        sigma(sigma), epsilon(epsilon), xparticle(xparticle), yparticle(yparticle), sx(sx), sy(sy), sz(sz), ex(ex), ey(ey), ez(ez) {
    }
};

/**
 * This is an internal class used to record information about an exception.
 * @private
 */
class GayBerneForce::ExceptionInfo {
public:
    int particle1, particle2;
    double sigma, epsilon;
    ExceptionInfo() {
        particle1 = particle2 = -1;
        sigma = epsilon = 0.0;
    }
    ExceptionInfo(int particle1, int particle2, double sigma, double epsilon) :
        particle1(particle1), particle2(particle2), sigma(sigma), epsilon(epsilon) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_GAYBERNEFORCE_H_*/
