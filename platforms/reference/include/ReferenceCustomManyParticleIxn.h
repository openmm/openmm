
/* Portions copyright (c) 2009-2014 Stanford University and Simbios.
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

#ifndef __ReferenceCustomManyParticleIxn_H__
#define __ReferenceCustomManyParticleIxn_H__

#include "ReferenceBondIxn.h"
#include "openmm/CustomManyParticleForce.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/ParsedExpression.h"
#include <map>
#include <set>
#include <vector>

namespace OpenMM {

class ReferenceCustomManyParticleIxn {

   private:

      class ParticleTermInfo;
      class DistanceTermInfo;
      class AngleTermInfo;
      class DihedralTermInfo;
      int numParticlesPerSet, numPerParticleParameters, numTypes;
      bool useCutoff, usePeriodic, centralParticleMode;
      RealOpenMM cutoffDistance;
      OpenMM::RealVec periodicBoxVectors[3];
      Lepton::ExpressionProgram energyExpression;
      std::vector<std::vector<std::string> > particleParamNames;
      std::vector<std::set<int> > exclusions;
      std::vector<int> particleTypes;
      std::vector<int> orderIndex;
      std::vector<std::vector<int> > particleOrder;
      std::vector<ParticleTermInfo> particleTerms;
      std::vector<DistanceTermInfo> distanceTerms;
      std::vector<AngleTermInfo> angleTerms;
      std::vector<DihedralTermInfo> dihedralTerms;

      void loopOverInteractions(std::vector<int>& particles, int loopIndex, std::vector<OpenMM::RealVec>& atomCoordinates,
                                RealOpenMM** particleParameters, std::map<std::string, double>& variables, std::vector<OpenMM::RealVec>& forces,
                                RealOpenMM* totalEnergy) const;

      /**---------------------------------------------------------------------------------------

         Calculate custom interaction for one set of particles

         @param particles          the indices of the particles
         @param atomCoordinates    atom coordinates
         @param particleParameters particle parameter values (particleParameters[particleIndex][parameterIndex])
         @param variables          the values of variables that may appear in expressions
         @param forces             force array (forces added)
         @param totalEnergy        total energy

         --------------------------------------------------------------------------------------- */

      void calculateOneIxn(const std::vector<int>& particles, std::vector<OpenMM::RealVec>& atomCoordinates,
                           RealOpenMM** particleParameters, std::map<std::string, double>& variables, std::vector<OpenMM::RealVec>& forces,
                           RealOpenMM* totalEnergy) const;

      void computeDelta(int atom1, int atom2, RealOpenMM* delta, std::vector<OpenMM::RealVec>& atomCoordinates) const;

      static RealOpenMM computeAngle(RealOpenMM* vec1, RealOpenMM* vec2);


   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomManyParticleIxn(const OpenMM::CustomManyParticleForce& force);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomManyParticleIxn();

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(RealOpenMM distance);

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param vectors    the vectors defining the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(OpenMM::RealVec* vectors);

      /**---------------------------------------------------------------------------------------

         Calculate the interaction

         @param atomCoordinates    atom coordinates
         @param particleParameters particle parameter values (particleParameters[particleIndex][parameterIndex])
         @param globalParameters   the values of global parameters
         @param forces             force array (forces added)
         @param totalEnergy        total energy

         --------------------------------------------------------------------------------------- */

      void calculateIxn(std::vector<OpenMM::RealVec>& atomCoordinates, RealOpenMM** particleParameters,
                        const std::map<std::string, double>& globalParameters,
                        std::vector<OpenMM::RealVec>& forces, RealOpenMM* totalEnergy) const;

// ---------------------------------------------------------------------------------------

};

class ReferenceCustomManyParticleIxn::ParticleTermInfo {
public:
    std::string name;
    int atom, component;
    Lepton::ExpressionProgram forceExpression;
    ParticleTermInfo(const std::string& name, int atom, int component, const Lepton::ExpressionProgram& forceExpression) :
            name(name), atom(atom), component(component), forceExpression(forceExpression) {
    }
};

class ReferenceCustomManyParticleIxn::DistanceTermInfo {
public:
    std::string name;
    int p1, p2;
    Lepton::ExpressionProgram forceExpression;
    mutable RealOpenMM delta[ReferenceForce::LastDeltaRIndex];
    DistanceTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::ExpressionProgram& forceExpression) :
            name(name), p1(atoms[0]), p2(atoms[1]), forceExpression(forceExpression) {
    }
};

class ReferenceCustomManyParticleIxn::AngleTermInfo {
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

class ReferenceCustomManyParticleIxn::DihedralTermInfo {
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

} // namespace OpenMM

#endif // __ReferenceCustomManyParticleIxn_H__
