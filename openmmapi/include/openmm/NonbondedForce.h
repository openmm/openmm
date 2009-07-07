#ifndef OPENMM_NONBONDEDFORCE_H_
#define OPENMM_NONBONDEDFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2009 Stanford University and the Authors.      *
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

#include "Force.h"
#include "Vec3.h"
#include <map>
#include <set>
#include <vector>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements nonbonded interactions between particles, including a Coulomb force to represent
 * electrostatics and a Lennard-Jones force to represent van der Waals interactions.  It optionally supports
 * periodic boundary conditions and cutoffs for long range interactions.  Lennard-Jones interactions are
 * calculated with the Lorentz-Bertelot combining rule: it uses the arithmetic mean of the sigmas and the
 * geometric mean of the epsilons for the two interacting particles.
 *
 * To use this class, create a NonbondedForce object, then call addParticle() once for each particle in the
 * System to define its parameters.  The number of particles for which you define nonbonded parameters must
 * be exactly equal to the number of particles in the System, or else an exception will be thrown when you
 * try to create an OpenMMContext.  After a particle has been added, you can modify its force field parameters
 * by calling setParticleParameters().
 *
 * NonbondedForce also lets you specify "exceptions", particular pairs of particles whose interactions should be
 * computed based on different parameters than those defined for the individual particles.  This can be used to
 * completely exclude certain interactions from the force calculation, or to alter how they interact with each other.
 *
 * Many molecular force fields omit Coulomb and Lennard-Jones interactions between particles separated by one
 * or two bonds, while using modified parameters for those separated by three bonds (known as "1-4 interactions").
 * This class provides a convenience method for this case called createExceptionsFromBonds().  You pass to it
 * a list of bonds and the scale factors to use for 1-4 interactions.  It identifies all pairs of particles which
 * are separated by 1, 2, or 3 bonds, then automatically creates exceptions for them.
 */

class OPENMM_EXPORT NonbondedForce : public Force {
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
         * Interactions beyond the cutoff distance are ignored.  Coulomb interactions closer than the cutoff distance
         * are modified using the reaction field method.
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.  Coulomb interactions closer than the
         * cutoff distance are modified using the reaction field method.
         */
        CutoffPeriodic = 2,
        /**
         * Periodic boundary conditions are used, and Ewald summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        Ewald = 3,
        /**
         * Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the interaction of each particle
         * with all periodic copies of every other particle.
         */
        PME = 4
    };
    /**
     * Create a NonbondedForce.
     */
    NonbondedForce();
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
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     */
    void setCutoffDistance(double distance);
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the reciprocal space cutoff and separation
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     */
    double getEwaldErrorTolerance() const;
    /**
     * Get the error tolerance for Ewald summation.  This corresponds to the fractional error in the forces
     * which is acceptable.  This value is used to select the reciprocal space cutoff and separation
     * parameter so that the average error level will be less than the tolerance.  There is not a
     * rigorous guarantee that all forces on all atoms will be less than the tolerance, however.
     */
    void setEwaldErrorTolerance(double tol);
    /**
     * Get the vectors which define the axes of the periodic box (measured in nm).  If the NonbondedMethod
     * in use does not use periodic boundary conditions, these values will have no effect.
     *
     * Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
     * x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
     *
     * @param a      on exit, this contains the vector defining the first edge of the periodic box
     * @param b      on exit, this contains the vector defining the second edge of the periodic box
     * @param c      on exit, this contains the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Set the vectors which define the axes of the periodic box (measured in nm).  If the NonbondedMethod
     * in use does not use periodic boundary conditions, these values will have no effect.
     *
     * Currently, only rectangular boxes are supported.  This means that a, b, and c must be aligned with the
     * x, y, and z axes respectively.  Future releases may support arbitrary triclinic boxes.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(Vec3 a, Vec3 b, Vec3 c);
    /**
     * Add the nonbonded force parameters for a particle.  This should be called once for each particle
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     * For calculating the Lennard-Jones interaction between two particles, the arithmetic mean of the sigmas
     * and the geometric mean of the epsilons for the two interacting particles is used (the Lorentz-Bertelot
     * combining rule).
     *
     * @param charge    the charge of the particle, measured in units of the proton charge
     * @param sigma     the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon   the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     * @return the index of the particle that was added
     */
    int addParticle(double charge, double sigma, double epsilon);
    /**
     * Get the nonbonded force parameters for a particle.
     *
     * @param index     the index of the particle for which to get parameters
     * @param charge    the charge of the particle, measured in units of the proton charge
     * @param sigma     the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon   the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     */
    void getParticleParameters(int index, double& charge, double& sigma, double& epsilon) const;
    /**
     * Set the nonbonded force parameters for a particle.  When calculating the Lennard-Jones interaction between two particles,
     * it uses the arithmetic mean of the sigmas and the geometric mean of the epsilons for the two interacting particles
     * (the Lorentz-Bertelot combining rule).
     *
     * @param index     the index of the particle for which to set parameters
     * @param charge    the charge of the particle, measured in units of the proton charge
     * @param sigma     the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon   the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     */
    void setParticleParameters(int index, double charge, double sigma, double epsilon);
    /**
     * Add an interaction to the list of exceptions that should be calculated differently from other interactions.
     * If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from
     * force and energy calculations.
     *
     * In many cases, you can use createExceptionsFromBonds() rather than adding each exception explicitly.
     *
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
     * @param sigma      the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon    the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     * @return the index of the exception that was added
     */
    int addException(int particle1, int particle2, double chargeProd, double sigma, double epsilon);
    /**
     * Get the force field parameters for an interaction that should be calculated differently from others.
     * 
     * @param index      the index of the interaction for which to get parameters
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
     * @param sigma      the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon    the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     */
    void getExceptionParameters(int index, int& particle1, int& particle2, double& chargeProd, double& sigma, double& epsilon) const;
    /**
     * Set the force field parameters for an interaction that should be calculated differently from others.
     * If chargeProd and epsilon are both equal to 0, this will cause the interaction to be completely omitted from
     * force and energy calculations.
     * 
     * @param index      the index of the interaction for which to get parameters
     * @param particle1  the index of the first particle involved in the interaction
     * @param particle2  the index of the second particle involved in the interaction
     * @param chargeProd the scaled product of the atomic charges (i.e. the strength of the Coulomb interaction), measured in units of the proton charge squared
     * @param sigma      the sigma parameter of the Lennard-Jones potential (corresponding to the van der Waals radius of the particle), measured in nm
     * @param epsilon    the epsilon parameter of the Lennard-Jones potential (corresponding to the well depth of the van der Waals interaction), measured in kJ/mol
     */
    void setExceptionParameters(int index, int particle1, int particle2, double chargeProd, double sigma, double epsilon);
    /**
     * Identify exceptions based on the molecular topology.  Particles which are separated by one or two bonds are set
     * to not interact at all, while pairs of particles separated by three bonds (known as "1-4 interactions") have
     * their Coulomb and Lennard-Jones interactions reduced by a fixed factor.
     *
     * @param bonds           the set of bonds based on which to construct exceptions.  Each element specifies the indices of
     *                        two particles that are bonded to each other.
     * @param coulomb14Scale  pairs of particles separated by three bonds will have the strength of their Coulomb interaction
     *                        multiplied by this factor
     * @param lj14Scale       pairs of particles separated by three bonds will have the strength of their Lennard-Jones interaction
     *                        multiplied by this factor
     */
    void createExceptionsFromBonds(const std::vector<std::pair<int, int> >& bonds, double coulomb14Scale, double lj14Scale);
protected:
    ForceImpl* createImpl();
private:
    class ParticleInfo;
    class ExceptionInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance, ewaldErrorTol;
    Vec3 periodicBoxVectors[3];
    void addExclusionsToSet(const std::vector<std::set<int> >& bonded12, std::set<int>& exclusions, int baseParticle, int fromParticle, int currentLevel) const;


// Retarded visual studio compiler complains about being unable to 
// export private stl class members.
// This stanza explains that it should temporarily shut up.
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable:4251)
#endif
    std::vector<ParticleInfo> particles;
    std::vector<ExceptionInfo> exceptions;
#if defined(_MSC_VER)
#pragma warning(pop)
#endif

};

class NonbondedForce::ParticleInfo {
public:
    double charge, sigma, epsilon;
    ParticleInfo() {
        charge = sigma = epsilon = 0.0;
    }
    ParticleInfo(double charge, double sigma, double epsilon) :
        charge(charge), sigma(sigma), epsilon(epsilon) {
    }
};

class NonbondedForce::ExceptionInfo {
public:
    int particle1, particle2;
    double chargeProd, sigma, epsilon;
    ExceptionInfo() {
        particle1 = particle2 = -1;
        chargeProd = sigma = epsilon = 0.0;
    }
    ExceptionInfo(int particle1, int particle2, double chargeProd, double sigma, double epsilon) :
        particle1(particle1), particle2(particle2), chargeProd(chargeProd), sigma(sigma), epsilon(epsilon) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_NONBONDEDFORCE_H_*/
