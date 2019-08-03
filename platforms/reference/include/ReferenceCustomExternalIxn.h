
/* Portions copyright (c) 2009-2018 Stanford University and Simbios.
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

#ifndef __ReferenceCustomExternalIxn_H__
#define __ReferenceCustomExternalIxn_H__

#include "ReferenceCustomExternalIxn.h"
#include "openmm/Vec3.h"
#include "lepton/CompiledExpression.h"

namespace OpenMM {

class ReferenceCustomExternalIxn {

   private:
      Lepton::CompiledExpression energyExpression;
      Lepton::CompiledExpression forceExpressionX;
      Lepton::CompiledExpression forceExpressionY;
      Lepton::CompiledExpression forceExpressionZ;
      std::vector<double*> energyParams;
      std::vector<double*> forceXParams;
      std::vector<double*> forceYParams;
      std::vector<double*> forceZParams;
      double *energyX, *energyY, *energyZ;
      double *forceXX, *forceXY, *forceXZ;
      double *forceYX, *forceYY, *forceYZ;
      double *forceZX, *forceZY, *forceZZ;
      int numParameters;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomExternalIxn(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpressionX,
                              const Lepton::CompiledExpression& forceExpressionY, const Lepton::CompiledExpression& forceExpressionZ,
                              const std::vector<std::string>& parameterNames);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomExternalIxn();

       /**---------------------------------------------------------------------------------------
      
         Set the values of all global parameters.
      
         --------------------------------------------------------------------------------------- */
      
       void setGlobalParameters(std::map<std::string, double> parameters);

      /**---------------------------------------------------------------------------------------

         Calculate Custom External Force

         @param atomIndex        the index of the atom to apply the force to
         @param atomCoordinates  atom coordinates
         @param parameters       parameter values
         @param forces           force array (forces added)
         @param energy           energy is added to this

         --------------------------------------------------------------------------------------- */

      void calculateForce(int atomIndex, std::vector<OpenMM::Vec3>& atomCoordinates,
                          std::vector<double>& parameters, std::vector<OpenMM::Vec3>& forces, double* energy) const;


};

} // namespace OpenMM

#endif // _ReferenceCustomBondIxn___
