
/* Portions copyright (c) 2009-2010 Stanford University and Simbios.
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

#ifndef __ReferenceCustomHbondIxn_H__
#define __ReferenceCustomHbondIxn_H__

#include "ReferenceBondIxn.h"
#include "lepton/ExpressionProgram.h"
#include <map>
#include <vector>

// ---------------------------------------------------------------------------------------

class ReferenceCustomHbondIxn : public ReferenceBondIxn {

   private:

      bool cutoff;
      bool periodic;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;
      Lepton::ExpressionProgram energyExpression;
      Lepton::ExpressionProgram rForceExpression;
      Lepton::ExpressionProgram thetaForceExpression;
      Lepton::ExpressionProgram psiForceExpression;
      Lepton::ExpressionProgram chiForceExpression;
      std::vector<std::string> donorParamNames, acceptorParamNames;
      std::vector<std::pair<int, int> > donorAtoms, acceptorAtoms;

      /**---------------------------------------------------------------------------------------

         Calculate custom interaction between a donor and an acceptor

         @param donor            the index of the donor
         @param acceptor         the index of the acceptor
         @param atomCoordinates  atom coordinates
         @param variables        the values of variables that may appear in expressions
         @param forces           force array (forces added)
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn(int donor, int acceptor, RealOpenMM** atomCoordinates,
                           std::map<std::string, double>& variables, RealOpenMM** forces,
                           RealOpenMM* totalEnergy) const;

      static RealOpenMM computeAngle(RealOpenMM* vec1, RealOpenMM* vec2, RealOpenMM sign);


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomHbondIxn(const std::vector<std::pair<int, int> >& donorAtoms, const std::vector<std::pair<int, int> >& acceptorAtoms,
                               const Lepton::ExpressionProgram& energyExpression, const Lepton::ExpressionProgram& rForceExpression,
                               const Lepton::ExpressionProgram& thetaForceExpression, const Lepton::ExpressionProgram& psiForceExpression,
                               const Lepton::ExpressionProgram& chiForceExpression, const std::vector<std::string>& donorParameterNames,
                               const std::vector<std::string>& acceptorParameterNames);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomHbondIxn();

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(RealOpenMM distance);

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param boxSize             the X, Y, and Z widths of the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(RealOpenMM* boxSize);

      /**---------------------------------------------------------------------------------------

         Calculate custom hbond interaction

         @param atomCoordinates    atom coordinates
         @param donorParameters    donor parameters values       donorParameters[donorIndex][parameterIndex]
         @param acceptorParameters acceptor parameters values    acceptorParameters[acceptorIndex][parameterIndex]
         @param exclusions         exclusion indices             exclusions[donorIndex][acceptorToExcludeIndex]
                                   exclusions[donorIndex][0] = number of exclusions
                                   exclusions[donorIndex][no.-1] = indices of acceptors to excluded from
                                   interacting w/ donor donorIndex
         @param globalParameters   the values of global parameters
         @param forces             force array (forces added)
         @param totalEnergy        total energy

         --------------------------------------------------------------------------------------- */

      void calculatePairIxn(RealOpenMM** atomCoordinates, RealOpenMM** donorParameters, RealOpenMM** acceptorParameters,
                            int** exclusions, const std::map<std::string, double>& globalParameters,
                            RealOpenMM** forces, RealOpenMM* totalEnergy) const;

// ---------------------------------------------------------------------------------------

};

#endif // __ReferenceCustomHbondIxn_H__
