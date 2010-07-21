
/* Portions copyright (c) 2009 Stanford University and Simbios.
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
#include "lepton/ExpressionProgram.h"

// ---------------------------------------------------------------------------------------

class ReferenceCustomExternalIxn {

   private:
      Lepton::ExpressionProgram energyExpression;
      Lepton::ExpressionProgram forceExpressionX;
      Lepton::ExpressionProgram forceExpressionY;
      Lepton::ExpressionProgram forceExpressionZ;
      std::vector<std::string> paramNames;
      std::map<std::string, double> globalParameters;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomExternalIxn(const Lepton::ExpressionProgram& energyExpression, const Lepton::ExpressionProgram& forceExpressionX,
                              const Lepton::ExpressionProgram& forceExpressionY, const Lepton::ExpressionProgram& forceExpressionZ,
                              const std::vector<std::string>& parameterNames, std::map<std::string, double> globalParameters);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomExternalIxn( );

      /**---------------------------------------------------------------------------------------

         Calculate Custom External Force

         @param atomIndex        the index of the atom to apply the force to
         @param atomCoordinates  atom coordinates
         @param parameters       parameter values
         @param forces           force array (forces added)
         @param energy           energy is added to this

         --------------------------------------------------------------------------------------- */

      void calculateForce( int atomIndex, RealOpenMM** atomCoordinates,
                            RealOpenMM* parameters, RealOpenMM** forces, RealOpenMM* energy ) const;


};

// ---------------------------------------------------------------------------------------

#endif // _ReferenceCustomBondIxn___
