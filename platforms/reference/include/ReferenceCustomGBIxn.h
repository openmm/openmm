
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

#ifndef __ReferenceCustomGBIxn_H__
#define __ReferenceCustomGBIxn_H__

#include "ReferenceNeighborList.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/CustomGBForce.h"
#include <map>
#include <set>
#include <vector>

namespace OpenMM {

class ReferenceCustomGBIxn {

   private:

      bool cutoff;
      bool periodic;
      const OpenMM::NeighborList* neighborList;
      OpenMM::Vec3 periodicBoxVectors[3];
      double cutoffDistance;
      CompiledExpressionSet expressionSet;
      std::vector<Lepton::CompiledExpression> valueExpressions;
      std::vector<std::vector<Lepton::CompiledExpression> > valueDerivExpressions;
      std::vector<std::vector<Lepton::CompiledExpression> > valueGradientExpressions;
      std::vector<std::vector<Lepton::CompiledExpression> > valueParamDerivExpressions;
      std::vector<OpenMM::CustomGBForce::ComputationType> valueTypes;
      std::vector<Lepton::CompiledExpression> energyExpressions;
      std::vector<std::vector<Lepton::CompiledExpression> > energyDerivExpressions;
      std::vector<std::vector<Lepton::CompiledExpression> > energyGradientExpressions;
      std::vector<std::vector<Lepton::CompiledExpression> > energyParamDerivExpressions;
      std::vector<OpenMM::CustomGBForce::ComputationType> energyTypes;
      std::vector<int> paramIndex;
      std::vector<int> valueIndex;
      std::vector<int> particleParamIndex;
      std::vector<int> particleValueIndex;
      int rIndex, xIndex, yIndex, zIndex;
      std::vector<std::vector<double> > values, dEdV;
      std::vector<std::vector<std::vector<double> > > dValuedParam;

      /**---------------------------------------------------------------------------------------

         Calculate a computed value of type SingleParticle

         @param index            the index of the value to compute
         @param numAtoms         number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]

         --------------------------------------------------------------------------------------- */

      void calculateSingleParticleValue(int index, int numAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters);

      /**---------------------------------------------------------------------------------------

         Calculate a computed value that is based on particle pairs

         @param index            the index of the value to compute
         @param numAtoms         number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param exclusions       exclusions[i] is the set of excluded indices for atom i
         @param useExclusions    specifies whether to use exclusions

         --------------------------------------------------------------------------------------- */

      void calculateParticlePairValue(int index, int numAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters,
                                      const std::vector<std::set<int> >& exclusions, bool useExclusions);

      /**---------------------------------------------------------------------------------------

         Evaluate a single atom pair as part of calculating a computed value

         @param index            the index of the value to compute
         @param atom1            the index of the first atom in the pair
         @param atom2            the index of the second atom in the pair
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]

         --------------------------------------------------------------------------------------- */

      void calculateOnePairValue(int index, int atom1, int atom2, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters);

      /**---------------------------------------------------------------------------------------

         Calculate an energy term of type SingleParticle

         @param index            the index of the value to compute
         @param numAtoms         number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param forces           forces on atoms are added to this
         @param totalEnergy      the energy contribution is added to this

         --------------------------------------------------------------------------------------- */

      void calculateSingleParticleEnergyTerm(int index, int numAtoms, std::vector<OpenMM::Vec3>& atomCoordinates,
                        std::vector<std::vector<double> >& atomParameters, std::vector<OpenMM::Vec3>& forces, double* totalEnergy, double* energyParamDerivs);

      /**---------------------------------------------------------------------------------------

         Calculate an energy term that is based on particle pairs

         @param index            the index of the term to compute
         @param numAtoms         number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param exclusions       exclusions[i] is the set of excluded indices for atom i
         @param useExclusions    specifies whether to use exclusions
         @param forces           forces on atoms are added to this
         @param totalEnergy      the energy contribution is added to this

         --------------------------------------------------------------------------------------- */

      void calculateParticlePairEnergyTerm(int index, int numAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters,
                                      const std::vector<std::set<int> >& exclusions, bool useExclusions,
                                      std::vector<OpenMM::Vec3>& forces, double* totalEnergy, double* energyParamDerivs);

      /**---------------------------------------------------------------------------------------

         Evaluate a single atom pair as part of calculating an energy term

         @param index            the index of the term to compute
         @param atom1            the index of the first atom in the pair
         @param atom2            the index of the second atom in the pair
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param forces           forces on atoms are added to this
         @param totalEnergy      the energy contribution is added to this

         --------------------------------------------------------------------------------------- */

      void calculateOnePairEnergyTerm(int index, int atom1, int atom2, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters,
                                 std::vector<OpenMM::Vec3>& forces, double* totalEnergy, double* energyParamDerivs);

      /**---------------------------------------------------------------------------------------

         Apply the chain rule to compute forces on atoms

         @param numAtoms         number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param exclusions       exclusions[i] is the set of excluded indices for atom i
         @param forces           forces on atoms are added to this

         --------------------------------------------------------------------------------------- */

      void calculateChainRuleForces(int numAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters,
                                      const std::vector<std::set<int> >& exclusions, std::vector<OpenMM::Vec3>& forces, double* energyParamDerivs);

      /**---------------------------------------------------------------------------------------

         Evaluate a single atom pair as part of applying the chain rule

         @param atom1            the index of the first atom in the pair
         @param atom2            the index of the second atom in the pair
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param forces           forces on atoms are added to this
         @param isExcluded       specifies whether this is an excluded pair

         --------------------------------------------------------------------------------------- */

      void calculateOnePairChainRule(int atom1, int atom2, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters,
                                 std::vector<OpenMM::Vec3>& forces, bool isExcluded);

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       ReferenceCustomGBIxn(const std::vector<Lepton::CompiledExpression>& valueExpressions,
                            const std::vector<std::vector<Lepton::CompiledExpression> > valueDerivExpressions,
                            const std::vector<std::vector<Lepton::CompiledExpression> > valueGradientExpressions,
                            const std::vector<std::vector<Lepton::CompiledExpression> > valueParamDerivExpressions,
                            const std::vector<std::string>& valueNames,
                            const std::vector<OpenMM::CustomGBForce::ComputationType>& valueTypes,
                            const std::vector<Lepton::CompiledExpression>& energyExpressions,
                            const std::vector<std::vector<Lepton::CompiledExpression> > energyDerivExpressions,
                            const std::vector<std::vector<Lepton::CompiledExpression> > energyGradientExpressions,
                            const std::vector<std::vector<Lepton::CompiledExpression> > energyParamDerivExpressions,
                            const std::vector<OpenMM::CustomGBForce::ComputationType>& energyTypes,
                            const std::vector<std::string>& parameterNames);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceCustomGBIxn();

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(double distance, const OpenMM::NeighborList& neighbors);

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param vectors    the vectors defining the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(OpenMM::Vec3* vectors);

      /**---------------------------------------------------------------------------------------

         Calculate custom GB ixn

         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param atomParameters   atomParameters[atomIndex][paramterIndex]
         @param exclusions       exclusions[i] is the set of excluded indices for atom i
         @param globalParameters the values of global parameters
         @param forces           force array (forces added)
         @param totalEnergy      total energy

         --------------------------------------------------------------------------------------- */

      void calculateIxn(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters, const std::vector<std::set<int> >& exclusions,
                       std::map<std::string, double>& globalParameters, std::vector<OpenMM::Vec3>& forces, double* totalEnergy, double* energyParamDerivs);

// ---------------------------------------------------------------------------------------

};

} // namespace OpenMM

#endif // __ReferenceCustomGBIxn_H__
