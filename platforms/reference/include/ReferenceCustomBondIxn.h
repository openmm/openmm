
/* Portions copyright (c) 2009-2013 Stanford University and Simbios.
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

#ifndef __ReferenceCustomBondIxn_H__
#define __ReferenceCustomBondIxn_H__

#include "ReferenceBondIxn.h"
#include "lepton/CompiledExpression.h"

// ---------------------------------------------------------------------------------------

class ReferenceCustomBondIxn : public ReferenceBondIxn {

   private:
      Lepton::CompiledExpression energyExpression;
      Lepton::CompiledExpression forceExpression;
      std::vector<double*> energyParams;
      std::vector<double*> forceParams;
      double* energyR;
      double* forceR;
      int numParameters;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomBondIxn(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpression,
                              const std::vector<std::string>& parameterNames, std::map<std::string, double> globalParameters);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomBondIxn( );

      /**---------------------------------------------------------------------------------------

         Calculate Custom Bond Ixn

         @param atomIndices      two bond indices
         @param atomCoordinates  atom coordinates
         @param parameters       parameter values
         @param forces           force array (forces added)
         @param totalEnergy      if not null, the energy will be added to this

         --------------------------------------------------------------------------------------- */

      void calculateBondIxn( int* atomIndices, std::vector<OpenMM::RealVec>& atomCoordinates,
                            RealOpenMM* parameters, std::vector<OpenMM::RealVec>& forces,
                            RealOpenMM* totalEnergy ) const;


};

// ---------------------------------------------------------------------------------------

#endif // _ReferenceCustomBondIxn___
