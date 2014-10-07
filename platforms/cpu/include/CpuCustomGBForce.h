
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

#ifndef OPENMM_CPU_CUSTOM_GB_FORCE_H__
#define OPENMM_CPU_CUSTOM_GB_FORCE_H__

#include "CompiledExpressionSet.h"
#include "CpuNeighborList.h"
#include "lepton/CompiledExpression.h"
#include "openmm/CustomGBForce.h"
#include "openmm/internal/vectorize.h"
#include <map>
#include <set>
#include <vector>

namespace OpenMM {

class CpuCustomGBForce {
private:
    bool cutoff;
    bool periodic;
    const CpuNeighborList* neighborList;
    float periodicBoxSize[3];
    float cutoffDistance;
    CompiledExpressionSet expressionSet;
    std::vector<Lepton::CompiledExpression> valueExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueDerivExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueGradientExpressions;
    std::vector<std::string> valueNames;
    std::vector<int> valueIndex;
    std::vector<CustomGBForce::ComputationType> valueTypes;
    std::vector<Lepton::CompiledExpression> energyExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyDerivExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyGradientExpressions;
    std::vector<std::string> paramNames;
    std::vector<int> paramIndex;
    std::vector<CustomGBForce::ComputationType> energyTypes;
    std::vector<int> particleParamIndex;
    std::vector<int> particleValueIndex;
    int xindex, yindex, zindex, rindex;
    // Workspace vectors
    std::vector<std::vector<float> > values, dEdV;
    std::vector<float> dVdR1, dVdR2, dVdX, dVdY, dVdZ;

    /**
     * Calculate a computed value of type SingleParticle
     * 
     * @param index            the index of the value to compute
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param values           the vector to store computed values into
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     */

    void calculateSingleParticleValue(int index, int numAtoms, float* posq, std::vector<std::vector<float> >& values,
                                      RealOpenMM** atomParameters);

    /**
     * Calculate a computed value that is based on particle pairs
     * 
     * @param index            the index of the value to compute
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param values           the vector to store computed values into
     * @param exclusions       exclusions[i] is the set of excluded indices for atom i
     * @param useExclusions    specifies whether to use exclusions
     */

    void calculateParticlePairValue(int index, int numAtoms, float* posq, RealOpenMM** atomParameters,
                                    std::vector<std::vector<float> >& values,
                                    const std::vector<std::set<int> >& exclusions, bool useExclusions, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Evaluate a single atom pair as part of calculating a computed value
     * 
     * @param index            the index of the value to compute
     * @param atom1            the index of the first atom in the pair
     * @param atom2            the index of the second atom in the pair
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param values           the vector to store computed values into
     */

    void calculateOnePairValue(int index, int atom1, int atom2, float* posq, RealOpenMM** atomParameters,
                               std::vector<std::vector<float> >& values, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Calculate an energy term of type SingleParticle
     * 
     * @param index            the index of the value to compute
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param values           the vector containing computed values
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param forces           forces on atoms are added to this
     * @param totalEnergy      the energy contribution is added to this
     * @param dEdV             the derivative of energy with respect to computed values is stored in this
     */

    void calculateSingleParticleEnergyTerm(int index, int numAtoms, float* posq, const std::vector<std::vector<float> >& values,
                                      RealOpenMM** atomParameters, float* forces,
                                      double* totalEnergy, std::vector<std::vector<float> >& dEdV);

    /**
     * Calculate an energy term that is based on particle pairs
     * 
     * @param index            the index of the term to compute
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param values           the vector containing computed values
     * @param exclusions       exclusions[i] is the set of excluded indices for atom i
     * @param useExclusions    specifies whether to use exclusions
     * @param forces           forces on atoms are added to this
     * @param totalEnergy      the energy contribution is added to this
     * @param dEdV             the derivative of energy with respect to computed values is stored in this
     */

    void calculateParticlePairEnergyTerm(int index, int numAtoms, float* posq, RealOpenMM** atomParameters,
                                    const std::vector<std::vector<float> >& values,
                                    const std::vector<std::set<int> >& exclusions, bool useExclusions,
                                    float* forces, double* totalEnergy, std::vector<std::vector<float> >& dEdV,
                                    const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Evaluate a single atom pair as part of calculating an energy term
     * 
     * @param index            the index of the term to compute
     * @param atom1            the index of the first atom in the pair
     * @param atom2            the index of the second atom in the pair
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param values           the vector containing computed values
     * @param forces           forces on atoms are added to this
     * @param totalEnergy      the energy contribution is added to this
     * @param dEdV             the derivative of energy with respect to computed values is stored in this
     */

    void calculateOnePairEnergyTerm(int index, int atom1, int atom2, float* posq, RealOpenMM** atomParameters,
                               const std::vector<std::vector<float> >& values,
                               float* forces, double* totalEnergy, std::vector<std::vector<float> >& dEdV,
                               const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Apply the chain rule to compute forces on atoms
     * 
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param values           the vector containing computed values
     * @param exclusions       exclusions[i] is the set of excluded indices for atom i
     * @param forces           forces on atoms are added to this
     * @param dEdV             the derivative of energy with respect to computed values is stored in this
     */

    void calculateChainRuleForces(int numAtoms, float* posq, RealOpenMM** atomParameters,
                                    const std::vector<std::vector<float> >& values,
                                    const std::vector<std::set<int> >& exclusions,
                                    float* forces, std::vector<std::vector<float> >& dEdV,
                                    const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Evaluate a single atom pair as part of applying the chain rule
     * 
     * @param atom1            the index of the first atom in the pair
     * @param atom2            the index of the second atom in the pair
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param values           the vector containing computed values
     * @param forces           forces on atoms are added to this
     * @param dEdV             the derivative of energy with respect to computed values is stored in this
     * @param isExcluded       specifies whether this is an excluded pair
     */

    void calculateOnePairChainRule(int atom1, int atom2, float* posq, RealOpenMM** atomParameters,
                               const std::vector<std::vector<float> >& values,
                               float* forces, std::vector<std::vector<float> >& dEdV,
                               bool isExcluded, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Compute the displacement and squared distance between two points, optionally using
     * periodic boundary conditions.
     */
    void getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const;

public:

    /**
     * Construct a new CpuCustomGBForce.
     */

     CpuCustomGBForce(int numAtoms, const std::vector<Lepton::CompiledExpression>& valueExpressions,
                          const std::vector<std::vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                          const std::vector<std::vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                          const std::vector<std::string>& valueNames,
                          const std::vector<CustomGBForce::ComputationType>& valueTypes,
                          const std::vector<Lepton::CompiledExpression>& energyExpressions,
                          const std::vector<std::vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                          const std::vector<std::vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                          const std::vector<CustomGBForce::ComputationType>& energyTypes,
                          const std::vector<std::string>& parameterNames);

     ~CpuCustomGBForce();

    /**
     * Set the force to use a cutoff.
     * 
     * @param distance            the cutoff distance
     * @param neighbors           the neighbor list to use
     */

    void setUseCutoff(float distance, const CpuNeighborList& neighbors);

    /**
     * Set the force to use periodic boundary conditions.  This requires that a cutoff has
     * already been set, and the smallest side of the periodic box is at least twice the cutoff
     * distance.
     * 
     * @param boxSize             the X, Y, and Z widths of the periodic box
     */

    void setPeriodic(RealVec& boxSize);

    /**
     * Calculate custom GB ixn
     * 
     * @param numberOfAtoms    number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param exclusions       exclusions[i] is the set of excluded indices for atom i
     * @param globalParameters the values of global parameters
     * @param forces           force array (forces added)
     * @param totalEnergy      total energy
     */

    void calculateIxn(int numberOfAtoms, float* posq, RealOpenMM** atomParameters, const std::vector<std::set<int> >& exclusions,
                     std::map<std::string, double>& globalParameters, float* forces, double* totalEnergy);
};

} // namespace OpenMM

#endif // OPENMM_CPU_CUSTOM_GB_FORCE_H__
