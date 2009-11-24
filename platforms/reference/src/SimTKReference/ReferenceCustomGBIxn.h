
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

#ifndef __ReferenceCustomGBIxn_H__
#define __ReferenceCustomGBIxn_H__

#include "ReferenceNeighborList.h"
#include "lepton/ExpressionProgram.h"
#include "openmm/CustomGBForce.h"
#include <map>
#include <set>
#include <vector>

// ---------------------------------------------------------------------------------------

class ReferenceCustomGBIxn {

   private:

      bool cutoff;
      bool periodic;
      const OpenMM::NeighborList* neighborList;
      RealOpenMM periodicBoxSize[3];
      RealOpenMM cutoffDistance;
      std::vector<Lepton::ExpressionProgram> valueExpressions;
      std::vector<std::vector<Lepton::ExpressionProgram> > valueDerivExpressions;
      std::vector<std::string> valueNames;
      std::vector<OpenMM::CustomGBForce::ComputationType> valueTypes;
      std::vector<Lepton::ExpressionProgram> energyExpressions;
      std::vector<std::vector<Lepton::ExpressionProgram> > energyDerivExpressions;
      std::vector<std::string> paramNames;
      std::vector<OpenMM::CustomGBForce::ComputationType> energyTypes;
      std::vector<std::string> particleParamNames;
      std::vector<std::string> particleValueNames;
      struct ComputedValue {
          RealOpenMM value;
          RealOpenMM gradient[3];
          ComputedValue() {
              value = (RealOpenMM) 0;
              gradient[0] = (RealOpenMM) 0;
              gradient[1] = (RealOpenMM) 0;
              gradient[2] = (RealOpenMM) 0;
          }
      };

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn between two atoms

         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][parameterIndex]
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateSingleParticleValue(int index, int numAtoms, std::vector<std::vector<ComputedValue> >& values,
                                        const std::map<std::string, double>& globalParameters, RealOpenMM** atomParameters) const;

      void calculateParticlePairValue(int index, int numAtoms, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
                                      std::vector<std::vector<ComputedValue> >& values,
                                      const std::map<std::string, double>& globalParameters,
                                      const std::vector<std::set<int> >& exclusions, bool useExclusions) const;

      void calculateOnePairValue(int index, int atom1, int atom2, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
                                 const std::map<std::string, double>& globalParameters,
                                 std::vector<std::vector<ComputedValue> >& values) const;

      void calculateSingleParticleEnergyTerm(int index, int numAtoms, const std::vector<std::vector<ComputedValue> >& values,
                                        const std::map<std::string, double>& globalParameters, RealOpenMM** atomParameters,
                                        RealOpenMM** forces, RealOpenMM* totalEnergy) const;

      void calculateParticlePairEnergyTerm(int index, int numAtoms, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
                                      const std::vector<std::vector<ComputedValue> >& values,
                                      const std::map<std::string, double>& globalParameters,
                                      const std::vector<std::set<int> >& exclusions, bool useExclusions,
                                      RealOpenMM** forces, RealOpenMM* totalEnergy) const;

      void calculateOnePairEnergyTerm(int index, int atom1, int atom2, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters,
                                 const std::map<std::string, double>& globalParameters,
                                 const std::vector<std::vector<ComputedValue> >& values,
                                 RealOpenMM** forces, RealOpenMM* totalEnergy) const;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomGBIxn(const std::vector<Lepton::ExpressionProgram>& valueExpressions,
                            const std::vector<std::vector<Lepton::ExpressionProgram> >& valueDerivExpressions,
                            const std::vector<std::string>& valueNames,
                            const std::vector<OpenMM::CustomGBForce::ComputationType>& valueTypes,
                            const std::vector<Lepton::ExpressionProgram>& energyExpressions,
                            const std::vector<std::vector<Lepton::ExpressionProgram> > energyDerivExpressions,
                            const std::vector<OpenMM::CustomGBForce::ComputationType>& energyTypes,
                            const std::vector<std::string>& parameterNames);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomGBIxn( );

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         @return ReferenceForce::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int setUseCutoff( RealOpenMM distance, const OpenMM::NeighborList& neighbors );

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param boxSize             the X, Y, and Z widths of the periodic box

         @return ReferenceForce::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int setPeriodic( RealOpenMM* boxSize );

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices                      exclusions[atomIndex][atomToExcludeIndex]
                                 exclusions[atomIndex][0] = number of exclusions
                                 exclusions[atomIndex][1-no.] = atom indices of atoms to excluded from
                                 interacting w/ atom atomIndex
         @param fixedParameters  non atom parameters (not currently used)
         @param globalParameters the values of global parameters
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy

         @return ReferenceForce::DefaultReturn

         --------------------------------------------------------------------------------------- */

      int calculateIxn(int numberOfAtoms, RealOpenMM** atomCoordinates, RealOpenMM** atomParameters, const std::vector<std::set<int> >& exclusions,
                       std::map<std::string, double>& globalParameters, RealOpenMM** forces,
                       RealOpenMM* energyByAtom, RealOpenMM* totalEnergy) const;

// ---------------------------------------------------------------------------------------

};

#endif // __ReferenceCustomGBIxn_H__
