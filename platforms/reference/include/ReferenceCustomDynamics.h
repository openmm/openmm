
/* Portions copyright (c) 2011 Stanford University and Simbios.
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
#include "lepton/ExpressionProgram.h"

#include <map>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------------------

class ReferenceCustomDynamics : public ReferenceDynamics {
private:

    const OpenMM::CustomIntegrator& integrator;
    std::vector<RealOpenMM> inverseMasses;
    std::vector<OpenMM::RealVec> sumBuffer, oldPos;
    std::vector<OpenMM::CustomIntegrator::ComputationType> stepType;
    std::vector<std::string> stepVariable, forceName, energyName;
    std::vector<Lepton::ExpressionProgram> stepExpression;
    std::vector<bool> invalidatesForces, needsForces, needsEnergy;
    std::vector<int> forceGroup;
    RealOpenMM energy;
    Lepton::ExpressionProgram kineticEnergyExpression;
    bool kineticEnergyNeedsForce;
    
    void computePerDof(int numberOfAtoms, std::vector<OpenMM::RealVec>& results, const std::vector<OpenMM::RealVec>& atomCoordinates,
                  const std::vector<OpenMM::RealVec>& velocities, const std::vector<OpenMM::RealVec>& forces, const std::vector<RealOpenMM>& masses,
                  const std::map<std::string, RealOpenMM>& globals, const std::vector<std::vector<OpenMM::RealVec> >& perDof,
                  const Lepton::ExpressionProgram& expression, const std::string& forceName);
    
    void recordChangedParameters(OpenMM::ContextImpl& context, std::map<std::string, RealOpenMM>& globals);
      
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
     
      void update(OpenMM::ContextImpl& context, int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                  std::vector<OpenMM::RealVec>& velocities, std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& masses,
                  std::map<std::string, RealOpenMM>& globals, std::vector<std::vector<OpenMM::RealVec> >& perDof, bool& forcesAreValid, RealOpenMM tolerance);
      
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

       double computeKineticEnergy(OpenMM::ContextImpl& context, int numberOfAtoms, std::vector<OpenMM::RealVec>& atomCoordinates,
                                   std::vector<OpenMM::RealVec>& velocities, std::vector<OpenMM::RealVec>& forces, std::vector<RealOpenMM>& masses,
                                   std::map<std::string, RealOpenMM>& globals, std::vector<std::vector<OpenMM::RealVec> >& perDof, bool& forcesAreValid);
};

// ---------------------------------------------------------------------------------------

#endif // __ReferenceCustomDynamics_H__
