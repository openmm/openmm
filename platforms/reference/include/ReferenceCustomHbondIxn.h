
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
#include "lepton/ParsedExpression.h"
#include <map>
#include <set>
#include <vector>

namespace OpenMM {

class ReferenceCustomHbondIxn : public ReferenceBondIxn {

   private:

      class DistanceTermInfo;
      class AngleTermInfo;
      class DihedralTermInfo;
      bool cutoff;
      bool periodic;
      OpenMM::Vec3 periodicBoxVectors[3];
      double cutoffDistance;
      std::vector<std::vector<int> > donorAtoms, acceptorAtoms;
      Lepton::ExpressionProgram energyExpression;
      std::vector<std::string> donorParamNames, acceptorParamNames;
      std::vector<DistanceTermInfo> distanceTerms;
      std::vector<AngleTermInfo> angleTerms;
      std::vector<DihedralTermInfo> dihedralTerms;

      /**---------------------------------------------------------------------------------------

         Calculate custom interaction between a donor and an acceptor

         @param donor            the index of the donor
         @param acceptor         the index of the acceptor
         @param atomCoordinates  atom coordinates
         @param variables        the values of variables that may appear in expressions
         @param forces           force array (forces added)
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn(int donor, int acceptor, std::vector<OpenMM::Vec3>& atomCoordinates,
                           std::map<std::string, double>& variables, std::vector<OpenMM::Vec3>& forces,
                           double* totalEnergy) const;

      void computeDelta(int atom1, int atom2, double* delta, std::vector<OpenMM::Vec3>& atomCoordinates) const;

      static double computeAngle(double* vec1, double* vec2);


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomHbondIxn(const std::vector<std::vector<int> >& donorAtoms, const std::vector<std::vector<int> >& acceptorAtoms,
                               const Lepton::ParsedExpression& energyExpression, const std::vector<std::string>& donorParameterNames,
                               const std::vector<std::string>& acceptorParameterNames, const std::map<std::string, std::vector<int> >& distances,
                               const std::map<std::string, std::vector<int> >& angles, const std::map<std::string, std::vector<int> >& dihedrals);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomHbondIxn();

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(double distance);

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param vectors    the vectors defining the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(OpenMM::Vec3* vectors);

      /**---------------------------------------------------------------------------------------

         Get the list of atoms for each donor group.

         --------------------------------------------------------------------------------------- */

      const std::vector<std::vector<int> >& getDonorAtoms() const {
          return donorAtoms;
      }

      /**---------------------------------------------------------------------------------------

         Get the list of atoms for each acceptor group.

         --------------------------------------------------------------------------------------- */

      const std::vector<std::vector<int> >& getAcceptorAtoms() const {
          return acceptorAtoms;
      }

      /**---------------------------------------------------------------------------------------

         Calculate custom hbond interaction

         @param atomCoordinates    atom coordinates
         @param donorParameters    donor parameters values       donorParameters[donorIndex][parameterIndex]
         @param acceptorParameters acceptor parameters values    acceptorParameters[acceptorIndex][parameterIndex]
         @param exclusions         exclusion indices
                                   exclusions[donorIndex] contains the list of excluded acceptors for that donor
         @param globalParameters   the values of global parameters
         @param forces             force array (forces added)
         @param totalEnergy        total energy

         --------------------------------------------------------------------------------------- */

      void calculatePairIxn(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& donorParameters, std::vector<std::vector<double> >& acceptorParameters,
                            std::vector<std::set<int> >& exclusions, const std::map<std::string, double>& globalParameters,
                            std::vector<OpenMM::Vec3>& forces, double* totalEnergy) const;

// ---------------------------------------------------------------------------------------

};

class ReferenceCustomHbondIxn::DistanceTermInfo {
public:
    std::string name;
    int p1, p2;
    Lepton::ExpressionProgram forceExpression;
    mutable double delta[ReferenceForce::LastDeltaRIndex];
    DistanceTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomHbondIxn::AngleTermInfo {
public:
    std::string name;
    int p1, p2, p3;
    Lepton::ExpressionProgram forceExpression;
    mutable double delta1[ReferenceForce::LastDeltaRIndex];
    mutable double delta2[ReferenceForce::LastDeltaRIndex];
    AngleTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomHbondIxn::DihedralTermInfo {
public:
    std::string name;
    int p1, p2, p3, p4;
    Lepton::ExpressionProgram forceExpression;
    mutable double delta1[ReferenceForce::LastDeltaRIndex];
    mutable double delta2[ReferenceForce::LastDeltaRIndex];
    mutable double delta3[ReferenceForce::LastDeltaRIndex];
    mutable double cross1[3];
    mutable double cross2[3];
    DihedralTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), p4(atoms[3]), forceExpression(forceExpression) {
    }
};

} // namespace OpenMM

#endif // __ReferenceCustomHbondIxn_H__
