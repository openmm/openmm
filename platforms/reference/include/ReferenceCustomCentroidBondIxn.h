
/* Portions copyright (c) 2009-2016 Stanford University and Simbios.
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

#ifndef __ReferenceCustomCentroidBondIxn_H__
#define __ReferenceCustomCentroidBondIxn_H__

#include "ReferenceBondIxn.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <vector>

namespace OpenMM {

class ReferenceCustomCentroidBondIxn : public ReferenceBondIxn {

   private:

      class PositionTermInfo;
      class DistanceTermInfo;
      class AngleTermInfo;
      class DihedralTermInfo;
      std::vector<std::vector<int> > groupAtoms;
      std::vector<std::vector<double> > normalizedWeights;
      std::vector<std::vector<int> > bondGroups;
      Lepton::ExpressionProgram energyExpression;
      std::vector<std::string> bondParamNames;
      std::vector<PositionTermInfo> positionTerms;
      std::vector<DistanceTermInfo> distanceTerms;
      std::vector<AngleTermInfo> angleTerms;
      std::vector<DihedralTermInfo> dihedralTerms;
      bool usePeriodic;
      RealVec boxVectors[3];


      /**---------------------------------------------------------------------------------------

         Calculate custom interaction for one bond

         @param bond             the index of the bond
         @param groupCenters     group center coordinates
         @param variables        the values of variables that may appear in expressions
         @param forces           force array (forces added)
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn(int bond, std::vector<OpenMM::RealVec>& groupCenters,
                           std::map<std::string, double>& variables, std::vector<OpenMM::RealVec>& forces,
                           RealOpenMM* totalEnergy) const;

      void computeDelta(int group1, int group2, RealOpenMM* delta, std::vector<OpenMM::RealVec>& groupCenters) const;

      static RealOpenMM computeAngle(RealOpenMM* vec1, RealOpenMM* vec2);


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomCentroidBondIxn(int numGroupsPerBond, const std::vector<std::vector<int> >& groupAtoms,
                               const std::vector<std::vector<double> >& normalizedWeights, const std::vector<std::vector<int> >& bondGroups, const Lepton::ParsedExpression& energyExpression,
                               const std::vector<std::string>& bondParameterNames, const std::map<std::string, std::vector<int> >& distances,
                               const std::map<std::string, std::vector<int> >& angles, const std::map<std::string, std::vector<int> >& dihedrals);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomCentroidBondIxn();

       /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.
      
         @param vectors    the vectors defining the periodic box
      
         --------------------------------------------------------------------------------------- */
      
       void setPeriodic(OpenMM::RealVec* vectors);

      /**---------------------------------------------------------------------------------------

         Get the list of groups in each bond.

         --------------------------------------------------------------------------------------- */

       const std::vector<std::vector<int> >& getBondGroups() const {
           return bondGroups;
       }

      /**---------------------------------------------------------------------------------------

         Calculate custom compound bond interaction

         @param atomCoordinates    atom coordinates
         @param bondParameters     bond parameters values       bondParameters[bondIndex][parameterIndex]
         @param globalParameters   the values of global parameters
         @param forces             force array (forces added)
         @param totalEnergy        total energy

         --------------------------------------------------------------------------------------- */

      void calculatePairIxn(std::vector<OpenMM::RealVec>& atomCoordinates, RealOpenMM** bondParameters,
                            const std::map<std::string, double>& globalParameters,
                            std::vector<OpenMM::RealVec>& forces, RealOpenMM* totalEnergy) const;

// ---------------------------------------------------------------------------------------

};

class ReferenceCustomCentroidBondIxn::PositionTermInfo {
public:
    std::string name;
    int group, component;
    Lepton::ExpressionProgram forceExpression;
    PositionTermInfo(const std::string& name, int group, int component, const Lepton::ExpressionProgram& forceExpression) :
            name(name), group(group), component(component), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCentroidBondIxn::DistanceTermInfo {
public:
    std::string name;
    int g1, g2;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta[ReferenceForce::LastDeltaRIndex];
    DistanceTermInfo(const std::string& name, const std::vector<int>& groups, const Lepton::ExpressionProgram& forceExpression) :
            name(name), g1(groups[0]), g2(groups[1]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCentroidBondIxn::AngleTermInfo {
public:
    std::string name;
    int g1, g2, g3;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta1[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM delta2[ReferenceForce::LastDeltaRIndex];
    AngleTermInfo(const std::string& name, const std::vector<int>& groups, const Lepton::ExpressionProgram& forceExpression) :
            name(name), g1(groups[0]), g2(groups[1]), g3(groups[2]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCentroidBondIxn::DihedralTermInfo {
public:
    std::string name;
    int g1, g2, g3, g4;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta1[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM delta2[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM delta3[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM cross1[3];
    mutable RealOpenMM cross2[3];
    DihedralTermInfo(const std::string& name, const std::vector<int>& groups, const Lepton::ExpressionProgram& forceExpression) :
            name(name), g1(groups[0]), g2(groups[1]), g3(groups[2]), g4(groups[3]), forceExpression(forceExpression) {
    }
};

} // namespace OpenMM

#endif // __ReferenceCustomCentroidBondIxn_H__
