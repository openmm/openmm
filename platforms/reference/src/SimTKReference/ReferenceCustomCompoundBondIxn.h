
/* Portions copyright (c) 2009-2012 Stanford University and Simbios.
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

#ifndef __ReferenceCustomCompoundBondIxn_H__
#define __ReferenceCustomCompoundBondIxn_H__

#include "ReferenceBondIxn.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <vector>

// ---------------------------------------------------------------------------------------

class ReferenceCustomCompoundBondIxn : public ReferenceBondIxn {

   private:

      class ParticleTermInfo;
      class DistanceTermInfo;
      class AngleTermInfo;
      class DihedralTermInfo;
      std::vector<std::vector<int> > bondAtoms;
      Lepton::ExpressionProgram energyExpression;
      std::vector<std::string> bondParamNames;
      std::vector<ParticleTermInfo> particleTerms;
      std::vector<DistanceTermInfo> distanceTerms;
      std::vector<AngleTermInfo> angleTerms;
      std::vector<DihedralTermInfo> dihedralTerms;


      /**---------------------------------------------------------------------------------------

         Calculate custom interaction for one bond

         @param bond             the index of the bond
         @param atomCoordinates  atom coordinates
         @param variables        the values of variables that may appear in expressions
         @param forces           force array (forces added)
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn(int bond, std::vector<OpenMM::RealVec>& atomCoordinates,
                           std::map<std::string, double>& variables, std::vector<OpenMM::RealVec>& forces,
                           RealOpenMM* totalEnergy) const;

      void computeDelta(int atom1, int atom2, RealOpenMM* delta, std::vector<OpenMM::RealVec>& atomCoordinates) const;

      static RealOpenMM computeAngle(RealOpenMM* vec1, RealOpenMM* vec2);


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomCompoundBondIxn(int numParticlesPerBond, const std::vector<std::vector<int> >& bondAtoms, const Lepton::ParsedExpression& energyExpression,
                               const std::vector<std::string>& bondParameterNames, const std::map<std::string, std::vector<int> >& distances,
                               const std::map<std::string, std::vector<int> >& angles, const std::map<std::string, std::vector<int> >& dihedrals);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomCompoundBondIxn();

      /**---------------------------------------------------------------------------------------

         Get the list atoms in each bond.

         --------------------------------------------------------------------------------------- */

       const std::vector<std::vector<int> >& getBondAtoms() const {
           return bondAtoms;
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

class ReferenceCustomCompoundBondIxn::ParticleTermInfo {
public:
    std::string name;
    int atom, component;
    Lepton::ExpressionProgram forceExpression;
    ParticleTermInfo(const std::string& name, int atom, int component, const Lepton::ExpressionProgram& forceExpression) :
            name(name), atom(atom), component(component), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCompoundBondIxn::DistanceTermInfo {
public:
    std::string name;
    int p1, p2;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta[ReferenceForce::LastDeltaRIndex];
    DistanceTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCompoundBondIxn::AngleTermInfo {
public:
    std::string name;
    int p1, p2, p3;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta1[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM delta2[ReferenceForce::LastDeltaRIndex];
    AngleTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCompoundBondIxn::DihedralTermInfo {
public:
    std::string name;
    int p1, p2, p3, p4;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta1[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM delta2[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM delta3[ReferenceForce::LastDeltaRIndex];
    mutable RealOpenMM cross1[3];
    mutable RealOpenMM cross2[3];
    DihedralTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), p4(atoms[3]), forceExpression(forceExpression) {
    }
};

#endif // __ReferenceCustomCompoundBondIxn_H__
