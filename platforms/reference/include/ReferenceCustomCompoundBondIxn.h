
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

#ifndef __ReferenceCustomCompoundBondIxn_H__
#define __ReferenceCustomCompoundBondIxn_H__

#include "ReferenceBondIxn.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <vector>

namespace OpenMM {

class ReferenceCustomCompoundBondIxn : public ReferenceBondIxn {

   private:

      class ParticleTermInfo;
      class DistanceTermInfo;
      class AngleTermInfo;
      class DihedralTermInfo;
      std::vector<std::vector<int> > bondAtoms;
      CompiledExpressionSet expressionSet;
      Lepton::CompiledExpression energyExpression;
      std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
      std::vector<int> bondParamIndex;
      std::vector<ParticleTermInfo> particleTerms;
      std::vector<DistanceTermInfo> distanceTerms;
      std::vector<AngleTermInfo> angleTerms;
      std::vector<DihedralTermInfo> dihedralTerms;
      int numParameters;
      bool usePeriodic;
      Vec3 boxVectors[3];


      /**---------------------------------------------------------------------------------------

         Calculate custom interaction for one bond

         @param bond             the index of the bond
         @param atomCoordinates  atom coordinates
         @param forces           force array (forces added)
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn(int bond, std::vector<OpenMM::Vec3>& atomCoordinates,
                           std::vector<OpenMM::Vec3>& forces, double* totalEnergy, double* energyParamDerivs);

      void computeDelta(int atom1, int atom2, double* delta, std::vector<OpenMM::Vec3>& atomCoordinates) const;

      static double computeAngle(double* vec1, double* vec2);


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomCompoundBondIxn(int numParticlesPerBond, const std::vector<std::vector<int> >& bondAtoms, const Lepton::ParsedExpression& energyExpression,
                               const std::vector<std::string>& bondParameterNames, const std::map<std::string, std::vector<int> >& distances,
                               const std::map<std::string, std::vector<int> >& angles, const std::map<std::string, std::vector<int> >& dihedrals,
                               const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomCompoundBondIxn();

       /**---------------------------------------------------------------------------------------
      
         Set the force to use periodic boundary conditions.
      
         @param vectors    the vectors defining the periodic box
      
         --------------------------------------------------------------------------------------- */
      
       void setPeriodic(OpenMM::Vec3* vectors);

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

      void calculatePairIxn(std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& bondParameters,
                            const std::map<std::string, double>& globalParameters,
                            std::vector<OpenMM::Vec3>& forces, double* totalEnergy, double* energyParamDerivs);

// ---------------------------------------------------------------------------------------

};

class ReferenceCustomCompoundBondIxn::ParticleTermInfo {
public:
    std::string name;
    int atom, component, index;
    Lepton::CompiledExpression forceExpression;
    ParticleTermInfo(const std::string& name, int atom, int component, const Lepton::CompiledExpression& forceExpression) :
            name(name), atom(atom), component(component), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCompoundBondIxn::DistanceTermInfo {
public:
    std::string name;
    int p1, p2, index;
    Lepton::CompiledExpression forceExpression;
    mutable double delta[ReferenceForce::LastDeltaRIndex];
    DistanceTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::CompiledExpression& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCompoundBondIxn::AngleTermInfo {
public:
    std::string name;
    int p1, p2, p3, index;
    Lepton::CompiledExpression forceExpression;
    mutable double delta1[ReferenceForce::LastDeltaRIndex];
    mutable double delta2[ReferenceForce::LastDeltaRIndex];
    AngleTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::CompiledExpression& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomCompoundBondIxn::DihedralTermInfo {
public:
    std::string name;
    int p1, p2, p3, p4, index;
    Lepton::CompiledExpression forceExpression;
    mutable double delta1[ReferenceForce::LastDeltaRIndex];
    mutable double delta2[ReferenceForce::LastDeltaRIndex];
    mutable double delta3[ReferenceForce::LastDeltaRIndex];
    mutable double cross1[3];
    mutable double cross2[3];
    DihedralTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::CompiledExpression& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), p3(atoms[2]), p4(atoms[3]), forceExpression(forceExpression) {
    }
};

} // namespace OpenMM

#endif // __ReferenceCustomCompoundBondIxn_H__
