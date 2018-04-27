
/* Portions copyright (c) 2011-2016 Stanford University and Simbios.
 * Contributors: Peter Eastman
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

#ifndef __ReferenceCustomDynamics_H__
#define __ReferenceCustomDynamics_H__

#include "ReferenceDynamics.h"
#include "openmm/CustomIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/CustomIntegratorUtilities.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/VectorExpression.h"
#include "lepton/CompiledExpression.h"

#include <map>
#include <string>
#include <vector>

namespace OpenMM {

class ReferenceCustomDynamics : public ReferenceDynamics {
private:

    class DerivFunction;
    const OpenMM::CustomIntegrator& integrator;
    std::vector<double> inverseMasses;
    std::vector<OpenMM::Vec3> sumBuffer, oldPos;
    std::vector<OpenMM::CustomIntegrator::ComputationType> stepType;
    std::vector<std::string> stepVariable;
    std::vector<std::vector<Lepton::CompiledExpression> > stepExpressions;
    std::vector<std::vector<VectorExpression> > stepVectorExpressions;
    std::vector<CustomIntegratorUtilities::Comparison> comparisons;
    std::vector<bool> invalidatesForces, needsForces, needsEnergy, computeBothForceAndEnergy;
    std::vector<int> forceGroupFlags, blockEnd;
    std::map<std::string, double> energyParamDerivs;
    Lepton::CompiledExpression kineticEnergyExpression;
    bool kineticEnergyNeedsForce;
    CompiledExpressionSet expressionSet;
    double x, v, m, f, energy, gaussian, uniform;
    int xIndex, vIndex;
    std::vector<int> perDofVariableIndex, stepVariableIndex;
    std::vector<double> perDofVariable;

    void initialize(OpenMM::ContextImpl& context, std::vector<double>& masses, std::map<std::string, double>& globals);
    
    Lepton::ExpressionTreeNode replaceDerivFunctions(const Lepton::ExpressionTreeNode& node, OpenMM::ContextImpl& context);
    
    void computePerDof(int numberOfAtoms, std::vector<OpenMM::Vec3>& results, const std::vector<OpenMM::Vec3>& atomCoordinates,
                  const std::vector<OpenMM::Vec3>& velocities, const std::vector<OpenMM::Vec3>& forces, const std::vector<double>& masses,
                  const std::vector<std::vector<OpenMM::Vec3> >& perDof, const Lepton::CompiledExpression& expression);
    
    void computePerParticle(int numberOfAtoms, std::vector<OpenMM::Vec3>& results, const std::vector<OpenMM::Vec3>& atomCoordinates,
                  const std::vector<OpenMM::Vec3>& velocities, const std::vector<OpenMM::Vec3>& forces, const std::vector<double>& masses,
                  const std::vector<std::vector<OpenMM::Vec3> >& perDof, const std::map<std::string, double>& globals, const VectorExpression& expression);
    
    void recordChangedParameters(OpenMM::ContextImpl& context, std::map<std::string, double>& globals);

    bool evaluateCondition(int step);
      
public:

      /**---------------------------------------------------------------------------------------
      
         Constructor

         @param numberOfAtoms  number of atoms
         @param integrator     the integrator definition to use
      
         --------------------------------------------------------------------------------------- */

       ReferenceCustomDynamics(int numberOfAtoms, const OpenMM::CustomIntegrator& integrator);

      /**---------------------------------------------------------------------------------------
      
         Destructor
      
         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomDynamics();

      /**---------------------------------------------------------------------------------------
      
         Update
      
         @param context             the context this integrator is updating
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param globals             a map containing values of global variables
         @param perDof              the values of per-DOF variables
         @param forcesAreValid      whether the current forces are valid or need to be recomputed
         @param tolerance           the constraint tolerance
      
         --------------------------------------------------------------------------------------- */
     
      void update(OpenMM::ContextImpl& context, int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                  std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses,
                  std::map<std::string, double>& globals, std::vector<std::vector<OpenMM::Vec3> >& perDof, bool& forcesAreValid, double tolerance);
      
      /**---------------------------------------------------------------------------------------
      
         Compute the kinetic energy of the system.
      
         @param context             the context this integrator is updating
         @param numberOfAtoms       number of atoms
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param globals             a map containing values of global variables
         @param perDof              the values of per-DOF variables
         @param forcesAreValid      whether the current forces are valid or need to be recomputed

         --------------------------------------------------------------------------------------- */

       double computeKineticEnergy(OpenMM::ContextImpl& context, int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                                   std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses,
                                   std::map<std::string, double>& globals, std::vector<std::vector<OpenMM::Vec3> >& perDof, bool& forcesAreValid);
};

} // namespace OpenMM

#endif // __ReferenceCustomDynamics_H__
