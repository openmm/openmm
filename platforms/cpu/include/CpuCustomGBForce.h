
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
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
#include <map>
#include <set>
#include <vector>

namespace OpenMM {

class CpuCustomGBForce {
private:
    class ComputeForceTask;
    class ThreadData;

    bool cutoff;
    bool periodic;
    const CpuNeighborList* neighborList;
    float periodicBoxSize[3];
    float cutoffDistance, cutoffDistance2;
    const std::vector<std::set<int> > exclusions;
    std::vector<std::string> valueNames;
    std::vector<CustomGBForce::ComputationType> valueTypes;
    std::vector<std::string> paramNames;
    std::vector<CustomGBForce::ComputationType> energyTypes;
    ThreadPool& threads;
    std::vector<ThreadData*> threadData;
    std::vector<double> threadEnergy;
    // Workspace vectors
    std::vector<std::vector<float> > values, dEdV;
    // The following variables are used to make information accessible to the individual threads.
    int numberOfAtoms;
    float* posq;
    RealOpenMM** atomParameters;
    const std::map<std::string, double>* globalParameters;
    std::vector<AlignedArray<float> >* threadForce;
    bool includeForce, includeEnergy;
    void* atomicCounter;
    
    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeForce(ThreadPool& threads, int threadIndex);

    /**
     * Calculate a computed value that is based on particle pairs
     * 
     * @param index            the index of the value to compute
     * @param data             workspace for the current thread
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param useExclusions    specifies whether to use exclusions
     */

    void calculateParticlePairValue(int index, ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters,
                                    bool useExclusions, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Evaluate a single atom pair as part of calculating a computed value
     * 
     * @param index            the index of the value to compute
     * @param atom1            the index of the first atom in the pair
     * @param atom2            the index of the second atom in the pair
     * @param data             workspace for the current thread
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     */

    void calculateOnePairValue(int index, int atom1, int atom2, ThreadData& data, float* posq, RealOpenMM** atomParameters,
                               std::vector<float>& valueArray, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Calculate an energy term of type SingleParticle
     * 
     * @param index            the index of the value to compute
     * @param data             workspace for the current thread
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param forces           forces on atoms are added to this
     * @param totalEnergy      the energy contribution is added to this
     */

    void calculateSingleParticleEnergyTerm(int index, ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters, float* forces, double& totalEnergy);

    /**
     * Calculate an energy term that is based on particle pairs
     * 
     * @param index            the index of the term to compute
     * @param data             workspace for the current thread
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param useExclusions    specifies whether to use exclusions
     * @param forces           forces on atoms are added to this
     * @param totalEnergy      the energy contribution is added to this
     */

    void calculateParticlePairEnergyTerm(int index, ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters,
                                    bool useExclusions, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Evaluate a single atom pair as part of calculating an energy term
     * 
     * @param index            the index of the term to compute
     * @param atom1            the index of the first atom in the pair
     * @param atom2            the index of the second atom in the pair
     * @param data             workspace for the current thread
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param forces           forces on atoms are added to this
     * @param totalEnergy      the energy contribution is added to this
     */

    void calculateOnePairEnergyTerm(int index, int atom1, int atom2, ThreadData& data, float* posq, RealOpenMM** atomParameters,
                               float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Apply the chain rule to compute forces on atoms
     * 
     * @param data             workspace for the current thread
     * @param numAtoms         number of atoms
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param forces           forces on atoms are added to this
     */

    void calculateChainRuleForces(ThreadData& data, int numAtoms, float* posq, RealOpenMM** atomParameters,
                                    float* forces, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Evaluate a single atom pair as part of applying the chain rule
     * 
     * @param atom1            the index of the first atom in the pair
     * @param atom2            the index of the second atom in the pair
     * @param data             workspace for the current thread
     * @param posq             atom coordinates
     * @param atomParameters   atomParameters[atomIndex][paramterIndex]
     * @param forces           forces on atoms are added to this
     * @param isExcluded       specifies whether this is an excluded pair
     */

    void calculateOnePairChainRule(int atom1, int atom2, ThreadData& data, float* posq, RealOpenMM** atomParameters,
                               float* forces, bool isExcluded, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Compute the displacement and squared distance between two points, optionally using
     * periodic boundary conditions.
     */
    void getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, bool periodic, const fvec4& boxSize, const fvec4& invBoxSize) const;

public:

    /**
     * Construct a new CpuCustomGBForce.
     */

     CpuCustomGBForce(int numAtoms, const std::vector<std::set<int> >& exclusions,
                        const std::vector<Lepton::CompiledExpression>& valueExpressions,
                        const std::vector<std::vector<Lepton::CompiledExpression> >& valueDerivExpressions,
                        const std::vector<std::vector<Lepton::CompiledExpression> >& valueGradientExpressions,
                        const std::vector<std::string>& valueNames,
                        const std::vector<CustomGBForce::ComputationType>& valueTypes,
                        const std::vector<Lepton::CompiledExpression>& energyExpressions,
                        const std::vector<std::vector<Lepton::CompiledExpression> >& energyDerivExpressions,
                        const std::vector<std::vector<Lepton::CompiledExpression> >& energyGradientExpressions,
                        const std::vector<CustomGBForce::ComputationType>& energyTypes,
                        const std::vector<std::string>& parameterNames, ThreadPool& threads);

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
     * @param globalParameters the values of global parameters
     * @param forces           force array (forces added)
     * @param totalEnergy      total energy
     */

    void calculateIxn(int numberOfAtoms, float* posq, RealOpenMM** atomParameters,
                     std::map<std::string, double>& globalParameters, std::vector<AlignedArray<float> >& threadForce, bool includeForce, bool includeEnergy, double& totalEnergy);
};

class CpuCustomGBForce::ThreadData {
public:
    ThreadData(int numAtoms, int numThreads, int threadIndex,
               const std::vector<Lepton::CompiledExpression>& valueExpressions,
               const std::vector<std::vector<Lepton::CompiledExpression> >& valueDerivExpressions,
               const std::vector<std::vector<Lepton::CompiledExpression> >& valueGradientExpressions,
               const std::vector<std::string>& valueNames,
               const std::vector<Lepton::CompiledExpression>& energyExpressions,
               const std::vector<std::vector<Lepton::CompiledExpression> >& energyDerivExpressions,
               const std::vector<std::vector<Lepton::CompiledExpression> >& energyGradientExpressions,
               const std::vector<std::string>& parameterNames);
    CompiledExpressionSet expressionSet;
    std::vector<Lepton::CompiledExpression> valueExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueDerivExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > valueGradientExpressions;
    std::vector<int> valueIndex;
    std::vector<Lepton::CompiledExpression> energyExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyDerivExpressions;
    std::vector<std::vector<Lepton::CompiledExpression> > energyGradientExpressions;
    std::vector<int> paramIndex;
    std::vector<int> particleParamIndex;
    std::vector<int> particleValueIndex;
    int xindex, yindex, zindex, rindex;
    int firstAtom, lastAtom;
    // Workspace vectors
    std::vector<float> value0, dVdR1, dVdR2, dVdX, dVdY, dVdZ;
    std::vector<std::vector<float> > dEdV;
};

} // namespace OpenMM

#endif // OPENMM_CPU_CUSTOM_GB_FORCE_H__
