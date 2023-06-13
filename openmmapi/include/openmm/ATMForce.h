#ifndef OPENMM_ATMFORCE_H_
#define OPENMM_ATMFORCE_H_

/* -------------------------------------------------------------------------- *
 *                    OpenMM's Alchemical Tranfer Force                       *
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

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include <vector>
#include <string>
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * 
 * The ATMForce class implements the Alchemical Transfer Method (ATM) for OpenMM
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
 */

class OPENMM_EXPORT ATMForce : public OpenMM::Force {
public:
    /**
     * Create an ATMForce object. 
     *
     * @param lambda1    the lambda1 parameter of the softplus alchemical potential (dimensionless)
     * @param lambda2    the lambda2 parameter of the softplus alchemical potential (dimensionless)
     * @param alpha      the alpha   parameter of the softplus alchemical potential (kJ/mol)^-1
     * @param u0         the u0      parameter of the softplus alchemical potential (kJ/mol)
     * @param w0         the w0      parameter of the softplus alchemical potential (kJ/mol)
     * @param direction  the direction parameter (dimensionless)
     *
     * @param umax       the umax    parameter of the softcore perturbation energy  (kJ/mol)
     * @param ubcore     the uc      parameter of the softcore perturbation energy  (kJ/mol)
     * @param acore      the a       parameter of the softcore perturbation energy  (dimensionless)
     *
     * @param VariableForceGroups  the list of force groups managed by the ATMForce
     *
     * The parameters provided in this constructor are added to OpenMM's Context as global parameters.
     * Their values can be changed by calling setParameter() on the Context using the parameter
     * names defined by the Lambda1(), Lambda2(), etc. methods below. 
     * For example: Context.setParameter(ATMForce::Lambda1(), 1.0)
     *
     * @return An ATMForce object
     */
  ATMForce(double lambda1, double lambda2, double alpha, double u0, double w0, double umax, double ubcore, double acore, double direction,
	   const std::vector<int>& VariableForceGroups): 
     defaultLambda1(lambda1), defaultLambda2(lambda2), defaultAlpha(alpha), defaultU0(u0), defaultW0(w0),
       defaultUmax(umax), defaultUbcore(ubcore), defaultAcore(acore), defaultDirection(direction), VariableForceGroups(VariableForceGroups) {
     }

     /**
     * Get the number of particles managed by ATMForce
     *
     * This should be the same number of particles as the OpenMM's System
     */
    int getNumParticles() const {
        return particles.size();
    }
    
    /**
     * Add a particle to the force.
     *
     * Normally, all of the particles in the OpenMM's System should be added to the ATMForce
     * in the same order as they appear in the System.
     *
     * @param particle    the index of the particle
     * @param dx, dy, dz  the displacement vector in nm
     * @return the index of the particle that was added
     */
    int addParticle(int particle, double dx, double dy, double dz);

    /**
     * Get the parameters for a particle
     * 
     * @param index      the index in the force for the particle for which to get parameters
     * @param particle   the index of the particle
     * @param dx, dy, dz the coordinates of the displacement vector in nm
     */
    void getParticleParameters(int index, int& particle, double& dx, double &dy, double &dz) const;

    /**
     * Set the parameters for a particle
     * 
     * @param index      the index in the force of the particle for which to set parameters
     * @param particle   the particle associated with this index
     * @param dx, dy, dz the coordinates of the displacement vector in nm
     */
    void setParticleParameters(int index, int particle, double dx, double dy, double dz);
    
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method 
     * should be called after updating parameters with setParticleParameters() to copy them over to the Context.
     * The only information this method updates is the values of per-particle parameters.  The number of particles
     * cannot be changed.
     */
    void updateParametersInContext(OpenMM::Context& context);


    /* TO DO */
    
    /**
     * return the number of Forces included in the ATM Force

    int getNumForces() const {
      return forces.size();
    }
     */
    
    /**
     * add a Force to the ATM Force.
     *
     * @param Force  The Force to the be added. ATMForce internally uses a copy of this Force.
     *               Subsequent modifications of the force will have no effect on the ATMForce.
     *               It is the responsibility of the caller to delete this Force if it is not
     *               needed.
     * @return       The index within ATMForce of the force that was added

    int addForce(Force* force);

     */

    /**
     * return the force from index

    Force& getForce(int index){
      ASSERT_VALID_INDEX(index, forces); 
      return *forces[index];
    }

     */

    
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
    double getPerturbationEnergy(const OpenMM::Context& context) const;


    /**
     * Returns the name of the global parameter corresponding to lambda1
     */
    static const std::string& Lambda1() {
        static const std::string key = "ATMLambda1";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to lambda2
     */
    static const std::string& Lambda2() {
        static const std::string key = "ATMLambda2";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to lambda2
     */
    static const std::string& Alpha() {
        static const std::string key = "ATMAlpha";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to u0
     */
    static const std::string& U0() {
        static const std::string key = "ATMU0";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to w0
     */
    static const std::string& W0() {
        static const std::string key = "ATMW0";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to umax
     */
    static const std::string& Umax() {
        static const std::string key = "ATMUmax";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to ubcore
     */
    static const std::string& Ubcore() {
        static const std::string key = "ATMUbcore";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to acore
     */
    static const std::string& Acore() {
        static const std::string key = "ATMAcore";
        return key;
    }

    /**
     * Returns the name of the global parameter corresponding to direction
     */
    static const std::string& Direction() {
        static const std::string key = "ATMDirection";
        return key;
    }

    /**
     *  Get the value of the lambda1 parameter passed in the constructor
     */
    double getDefaultLambda1() const {
        return defaultLambda1;
    }

    /**
     *  Get the value of the lambda2 parameter passed in the constructor
     */
    double getDefaultLambda2() const {
        return defaultLambda2;
    }

    /**
     *  Get the value of the alpha parameter passed in the constructor
     */
    double getDefaultAlpha() const {
        return defaultAlpha;
    }

    /**
     *  Get the value of the u0 parameter passed in the constructor
     */
    double getDefaultU0() const {
        return defaultU0;
    }

    /**
     *  Get the value of the w0 parameter passed in the constructor
     */
    double getDefaultW0() const {
        return defaultW0;
    }

    /**
     *  Get the value of the umax parameter passed in the constructor
     */
    double getDefaultUmax() const {
        return defaultUmax;
    }

    /**
     *  Get the value of the ubcore parameter passed in the constructor
     */
    double getDefaultUbcore() const {
        return defaultUbcore;
    }

    /**
     *  Get the value of the acore parameter passed in the constructor
     */
    double getDefaultAcore() const {
        return defaultAcore;
    }

    /**
     *  Get the value of the direction parameter passed in the constructor
     */
    double getDefaultDirection() const {
        return defaultDirection;
    }

    /**
     *  Get the list of the group ids of the Forces managed by the ATMForce 
     */
    const std::vector<int>& getVariableForceGroups() const {
      return VariableForceGroups;
    }

protected:
  OpenMM::ForceImpl* createImpl() const;
private:
    //  std::vector<Force *> forces;

    class ParticleInfo;
    std::vector<ParticleInfo> particles;

    //softplus parameters
    double defaultLambda1, defaultLambda2, defaultAlpha, defaultU0, defaultW0;
    //soft core parameters
    double defaultUmax, defaultUbcore, defaultAcore;
    //alchemical direction parameter
    double defaultDirection;

    //the forces that are recalculated after the coordinate transformation
    std::vector<int> VariableForceGroups;

};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class ATMForce::ParticleInfo {
 public:
  int particle;
  double dx, dy, dz;
  ParticleInfo() {
    particle = -1;
    dx = dy = dz = 0.0;
  }
  ParticleInfo(int particle) : particle(particle) {
    dx = dy = dz = 0.0;
  }
  ParticleInfo(int particle, double dx, double dy, double dz) : particle(particle), dx(dx), dy(dy), dz(dz) {
  }
};
 
} // namespace OpenMM

#endif /*OPENMM_ATMFORCE_H_*/
