#ifndef OPENMM_ATMFORCE_H_
#define OPENMM_ATMFORCE_H_

/* -------------------------------------------------------------------------- *
 *                    OpenMM's Alchemical Transfer Force                      *
 * -------------------------------------------------------------------------- *
 * This is a Force of the OpenMM molecular simulation toolkit                 *
 * that implements the Alchemical Transfer Potential                          *
 * for absolute and relative binding free energy estimation                   *
 * (https://doi.org/10.1021/acs.jcim.1c01129). The code is derived from the   *
 * ATMMetaForce plugin                                                        *
 * https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin               *
 * with support from the National Science Foundation CAREER 1750511           *
 *                                                                            *
 * Portions copyright (c) 2021-2023 by the Authors                            *
 * Authors: Emilio Gallicchio                                                 *
 * Contributors: Peter Eastman                                                *
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
#include "internal/AssertionUtilities.h"
#include <openmm/Vec3.h>
#include <vector>
#include <string>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * 
 * The ATMForce class implements the Alchemical Transfer Method (ATM) for OpenMM.
 * ATM is used to compute
 * the binding free energies of molecular complexes and of other equilibrium processes.
 * ATM and its implementation are described in the open access article:
 *
 *    Solmaz Azimi, Sheenam Khuttan, Joe Z. Wu, Rajat K. Pal, and Emilio  Gallicchio. 
 *    Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method. 
 *    J. Chem. Inf. Model.  62, 309 (2022)
 *    https://doi.org/10.1021/acs.jcim.1c01129
 *
 * Refer to the publication above for a detailed description of the ATM method and the parameters used in this API 
 * and please cite it to support our work if you use this software in your research.
 *
 * The ATMForce implements an arbitrary potential energy function that depends on the potential
 * energies (u0 and u1) of the system before and after a set of atoms are displaced by a specified amount.
 * For example, to displace a molecule from the solvent bulk to a receptor binding site to simulate 
 * a binding process.
 * The potential energy function typically depends on one or more parameters that are dialed to implement
 * alchemical transformations.
 *
 * To use this class, create a ATMForce object, passing an algebraic expression to the
 * constructor that defines the potential energy. This expression can be any combination
 * of the variables u0 and u1. Then call addGlobalParameter() to define the parameters on which the potential energy expression depends.
 * The values of global parameters may be modified during a simulation by calling Context::setParameter().
 * Next, call addForce() to add Force objects that define the terms of the potential energy function
 * that change upon displacement. Finally, call addParticle() to specify the displacement applied to
 * each particle. Displacements can be changed by calling setParticleParameters(). As any per-particle parameters, 
 * changes in displacements take effect only after calling updateParametersInContext().
 *
 * As an example, the following code creates a ATMForce based on the change in energy of
 * two particles when the second particle is displaced by 1 nm in the x direction.
 * The energy change is dialed using an alchemical parameter Lambda, which in this case is set to 1/2:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *    ATMForce *atmforce = new ATMForce("u0 + Lambda*(u1 - u0)");
 *    atm->addGlobalParameter("Lambda", 0.5);
 *    atm->addParticle( Vec3(0, 0, 0));
 *    atm->addParticle( Vec3(1, 0, 0));
 *    CustomBondForce* force = new CustomBondForce("0.5*r^2");
 *    atm->addForce(force);
 * \endverbatim
 *
 * Expressions may involve the operators + (add), - (subtract), * (multiply), / (divide), and ^ (power), and the following
 * functions: sqrt, exp, log, sin, cos, sec, csc, tan, cot, asin, acos, atan, atan2, sinh, cosh, tanh, erf, erfc, min, max, abs, floor, ceil, step, delta,
 * select.  All trigonometric functions
 * are defined in radians, and log is the natural logarithm.  step(x) = 0 if x is less than 0, 1 otherwise.  delta(x) = 1 if x is 0, 0 otherwise.
 * select(x,y,z) = z if x = 0, y otherwise.
 *
 * If instead of the energy expression the ATMForce constructor specifies the values of a series of parameters,
 * the default energy expression is used:
   \verbatim
   select(step(Direction), u0, u1) + ((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + Lambda2*usc + W0;
     usc = select(step(u-Ubcore), (Umax-Ubcore)*fsc+Ubcore, u), u);
     fsc = (z^Acore-1)/(z^Acore+1);
     z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;
     y = (u-Ubcore)/(Umax-Ubcore);
     u = Direction*(u1-u0)"
   \endverbatim
 * which is the same as the soft-core softplus alchemical potential energy function in the Azimi et al. paper above.
 *
 * The ATMForce is then added to the System as any other Force
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: cpp
 *
 *  system.addForce(atmforce);
 * \endverbatim
 *
 * after which it will be used for energy/force evaluations for molecular dynamics and energy optimization.
 *
 */

class OPENMM_EXPORT ATMForce : public OpenMM::Force {
public:
    /**
     * Create an ATMForce object. 
     *
     * @param energy   an algebraic expression giving the energy of the system as a function
     *                 of u0 and u1, the energies before and after displacement
     */
    explicit ATMForce(const std::string& energy);
    /**
     * Create an ATMForce object with the default energy expression.
     *
     * @param lambda1    the Lambda1 parameter of the softplus alchemical potential (dimensionless)
     * @param lambda2    the Lambda2 parameter of the softplus alchemical potential (dimensionless)
     * @param alpha      the Alpha   parameter of the softplus alchemical potential (kJ/mol)^-1
     * @param uh         the Uh      parameter of the softplus alchemical potential (kJ/mol)
     * @param w0         the W0      parameter of the softplus alchemical potential (kJ/mol)
     * @param umax       the Umax    parameter of the softcore perturbation energy  (kJ/mol)
     * @param ubcore     the Ubcore  parameter of the softcore perturbation energy  (kJ/mol)
     * @param acore      the Acore   parameter of the softcore perturbation energy  (dimensionless)
     * @param direction  the Direction parameter (dimensionless)
     *
     * The parameters provided in this constructor are added to OpenMM's Context as global parameters.
     * Their values can be changed by calling setParameter() on the Context using the parameter
     * names defined by the Lambda1(), Lambda2(), etc. methods below. 
     * For example: Context.setParameter(ATMForce::Lambda1(), 1.0)
     *
     * @return An ATMForce object
     */
    ATMForce(double lambda1, double lambda2, double alpha, double uh, double w0, double umax, double ubcore, double acore, double direction);
    ~ATMForce();

    /**
     * Get the number of particles managed by ATMForce
     *
     * This should be the same number of particles as the System
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of Forces included in the ATMForce.
     */
    int getNumForces() const {
        return forces.size();
    }
    /**
     * Get the number of global parameters that the interaction depends on.
     */
    int getNumGlobalParameters() const {
        return globalParameters.size();
    }
    /**
     * Get the algebraic expression that gives the energy of the system
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the energy of the system
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Add a Force that will be computed by the ATM orce.
     *
     * @param force  the Force to the be added, which should have been created on the heap with the
     *               "new" operator.  The ATMForce takes over ownership of it, and deletes the Force when the
     *               ATMForce itself is deleted.
     * @return       The index within ATMForce of the force that was added
     */
    int addForce(Force* force);
    /**
     * return the force from index
     */
    Force& getForce(int index) const;
    /**
     * Add a particle to the force.
     *
     * Normally, all of the particles in the OpenMM's System should be added to the ATMForce
     * in the same order as they appear in the System.
     *
     * @param displacement1    the displacement of the particle for the target state in nm
     * @param displacement0    the displacement of the particle for the initial state in nm
     * @return                 the index of the particle that was added
     */
    int addParticle(const Vec3& displacement1, const Vec3& displacement0=Vec3());
    /**
     * Get the parameters for a particle
     * 
     * @param index           the index in the force for the particle for which to get parameters
     * @param displacement1   the the displacement of the particle for the target state in nm
     * @param displacement0   the the displacement of the particle for the initial state in nm
     */
    void getParticleParameters(int index, Vec3& displacement1, Vec3& displacement0) const;
    /**
     * Set the parameters for a particle
     * 
     * @param index           the index in the force of the particle for which to set parameters
     * @param displacement1   the displacement of the particle for the target state in nm
     * @param displacement0   the displacement of the particle for the initial state in nm
     */
    void setParticleParameters(int index, const Vec3& displacement1, const Vec3& displacement0=Vec3());
    /**
     * Add a new global parameter that the interaction may depend on.  The default value provided to
     * this method is the initial value of the parameter in newly created Contexts.  You can change
     * the value at any time by calling setParameter() on the Context.
     *
     * @param name             the name of the parameter
     * @param defaultValue     the default value of the parameter
     * @return the index of the parameter that was added
     */
    int addGlobalParameter(const std::string& name, double defaultValue);
    /**
     * Get the name of a global parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getGlobalParameterName(int index) const;
    /**
     * Set the name of a global parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setGlobalParameterName(int index, const std::string& name);
    /**
     * Get the default value of a global parameter.
     *
     * @param index     the index of the parameter for which to get the default value
     * @return the parameter default value
     */
    double getGlobalParameterDefaultValue(int index) const;
    /**
     * Set the default value of a global parameter.
     *
     * @param index          the index of the parameter for which to set the default value
     * @param defaultValue   the default value of the parameter
     */
    void setGlobalParameterDefaultValue(int index, double defaultValue);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method 
     * should be called after updating parameters with setParticleParameters() to copy them over to the Context.
     * The only information this method updates is the values of per-particle parameters.  The number of particles
     * cannot be changed.
     */
    void updateParametersInContext(Context& context);
    /**
     * Always returns False. Included for compatibility.
     */
    bool usesPeriodicBoundaryConditions() const {
        return false; //the non-bonded force with PBC is in the system so it would be queried correctly
    }
    /**
     * Returns the current perturbation energy calculated by the ATMForce
     *
     * The perturbation energy is U2(x) - U1(x) (for direction = 1) or U1(x) - U2(x) (for direction = -1),
     * as further modified by the soft-core function.
     */
    void getPerturbationEnergy(const Context& context, double& u1, double& u0, double& energy) const;
    /**
     * Returns the name of the global parameter corresponding to lambda1
     */
    static const std::string& Lambda1() {
        static const std::string key = "Lambda1";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to lambda2
     */
    static const std::string& Lambda2() {
        static const std::string key = "Lambda2";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to lambda2
     */
    static const std::string& Alpha() {
        static const std::string key = "Alpha";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to uh
     */
    static const std::string& Uh() {
        static const std::string key = "Uh";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to w0
     */
    static const std::string& W0() {
        static const std::string key = "W0";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to umax
     */
    static const std::string& Umax() {
        static const std::string key = "Umax";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to ubcore
     */
    static const std::string& Ubcore() {
        static const std::string key = "Ubcore";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to acore
     */
    static const std::string& Acore() {
        static const std::string key = "Acore";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to direction
     */
    static const std::string& Direction() {
        static const std::string key = "Direction";
        return key;
    }

protected:
  ForceImpl* createImpl() const;
private:
    class ParticleInfo;
    class GlobalParameterInfo;
    std::string energyExpression;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<Force *> forces;
    std::vector<ParticleInfo> particles;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class ATMForce::ParticleInfo {
 public:
  int index;
  Vec3 displacement1;
  Vec3 displacement0;
  ParticleInfo() {
    index = -1;
    displacement1 = Vec3();
    displacement0 = Vec3();
  }
  ParticleInfo( int index ) : index(index) {
    displacement1 = Vec3();
    displacement0 = Vec3();
  }
  ParticleInfo(int index, Vec3 displacement1, Vec3 displacement0) :
    index(index), displacement1(displacement1), displacement0(displacement0) {
  }
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class ATMForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

} // namespace OpenMM

#endif /*OPENMM_ATMFORCE_H_*/
