
/* Portions copyright (c) 2009-2022 Stanford University and Simbios.
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

#ifndef OPENMM_CPU_CUSTOM_NONBONDED_FORCE_H__
#define OPENMM_CPU_CUSTOM_NONBONDED_FORCE_H__

#include "AlignedArray.h"
#include "CpuNeighborList.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
#include "lepton/CompiledVectorExpression.h"
#include "lepton/ParsedExpression.h"
#include <atomic>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {

class CpuCustomNonbondedForce {
public:
      /**---------------------------------------------------------------------------------------

         Constructor

         --------------------------------------------------------------------------------------- */

       CpuCustomNonbondedForce(ThreadPool& threads, const CpuNeighborList& neighbors);

       void initialize(const Lepton::ParsedExpression& energyExpression, const Lepton::ParsedExpression& forceExpression,
                       const std::vector<std::string>& parameterNames, const std::vector<std::set<int> >& exclusions,
                       const std::vector<Lepton::ParsedExpression> energyParamDerivExpressions,
                       const std::vector<std::string>& computedValueNames, const std::vector<Lepton::ParsedExpression> computedValueExpressions);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       virtual ~CpuCustomNonbondedForce();

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(double distance);

      /**---------------------------------------------------------------------------------------

         Restrict the force to a list of interaction groups.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         --------------------------------------------------------------------------------------- */

      void setInteractionGroups(const std::vector<std::pair<std::set<int>, std::set<int> > >& groups);

      /**---------------------------------------------------------------------------------------
      
         Set the force to use a switching function.
      
         @param distance            the switching distance
      
         --------------------------------------------------------------------------------------- */
      
      void setUseSwitchingFunction(double distance);

      /**---------------------------------------------------------------------------------------

         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.

         @param periodicBoxVectors    the vectors defining the periodic box

         --------------------------------------------------------------------------------------- */

      void setPeriodic(Vec3* periodicBoxVectors);

      /**---------------------------------------------------------------------------------------

         Calculate custom pair ixn

         @param numberOfAtoms    number of atoms
         @param posq             atom coordinates in float format
         @param atomCoordinates  atom coordinates
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param globalParameters the values of global parameters
         @param forces           force array (forces added)
         @param totalEnergy      total energy
         @param threads          the thread pool to use

         --------------------------------------------------------------------------------------- */

    void calculatePairIxn(int numberOfAtoms, float* posq, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<std::vector<double> >& atomParameters,
                          const std::map<std::string, double>& globalParameters, std::vector<AlignedArray<float> >& threadForce,
                          bool includeForce, bool includeEnergy, double& totalEnergy, double* energyParamDerivs);
protected:
    class ThreadData;

    bool cutoff;
    bool useSwitch;
    bool periodic;
    bool triclinic;
    bool useInteractionGroups;
    const CpuNeighborList* neighborList;
    float recipBoxSize[3];
    Vec3 periodicBoxVectors[3];
    AlignedArray<fvec4> periodicBoxVec4;
    double cutoffDistance, switchingDistance;
    ThreadPool& threads;
    std::vector<std::set<int> > exclusions;
    std::vector<ThreadData*> threadData;
    std::vector<std::string> paramNames, computedValueNames;
    std::vector<std::pair<int, int> > groupInteractions;
    std::vector<double> threadEnergy;
    std::vector<std::vector<double> > atomComputedValues;
    // The following variables are used to make information accessible to the individual threads.
    int numberOfAtoms;
    float* posq;
    Vec3 const* atomCoordinates;
    std::vector<double>* atomParameters;        
    const std::map<std::string, double>* globalParameters;
    std::vector<AlignedArray<float> >* threadForce;
    bool includeForce, includeEnergy;
    std::atomic<int> atomicCounter;

    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeForce(ThreadPool& threads, int threadIndex);

    /**
     * Calculate the interaction between two atoms.
     * 
     * @param atom1            the index of the first atom
     * @param atom2            the index of the second atom
     * @param data             workspace for the current thread
     * @param forces           force array (forces added)
     * @param totalEnergy      total energy
     * @param boxSize          the size of the periodic box
     * @param invBoxSize       the inverse size of the periodic box
     */
    void calculateOneIxn(int atom1, int atom2, ThreadData& data, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Calculate all the interactions for one block of atoms.
     * 
     * @param data            workspace for the current thread
     * @param blockIndex      the index of the atom block
     * @param forces          force array (forces added)
     * @param totalEnergy     total energy
     * @param boxSize         the size of the periodic box
     * @param invBoxSize       the inverse size of the periodic box
     */
    virtual void calculateBlockIxn(ThreadData& data, int blockIndex, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize) = 0;

    /**
     * Compute the displacement and squared distance between two points, optionally using
     * periodic boundary conditions.
     */
    void getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const;
};

class CpuCustomNonbondedForce::ThreadData {
public:
    ThreadData(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledVectorExpression& energyVecExpression,
            const Lepton::CompiledExpression& forceExpression, const Lepton::CompiledVectorExpression& forceVecExpression, const std::vector<std::string>& parameterNames,
            const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions, const std::vector<std::string>& computedValueNames,
            const std::vector<Lepton::CompiledExpression> computedValueExpressions, std::vector<std::vector<double> >& atomComputedValues);
    Lepton::CompiledExpression energyExpression, forceExpression;
    Lepton::CompiledVectorExpression energyVecExpression, forceVecExpression;
    std::vector<Lepton::CompiledExpression> computedValueExpressions, energyParamDerivExpressions;
    CompiledExpressionSet expressionSet;
    std::vector<double> particleParam, computedValues;
    std::vector<float> rvec, vecParticle1Params, vecParticle2Params, vecParticle1Values, vecParticle2Values;
    double r;
    std::vector<double> energyParamDerivs; 
    std::vector<std::vector<double> >& atomComputedValues;
};

/**
 * This function is called to create an instance of an appropriate subclass for the current CPU.
 */
CpuCustomNonbondedForce* createCpuCustomNonbondedForce(ThreadPool& threads, const CpuNeighborList& neighbors);

} // namespace OpenMM

#endif // OPENMM_CPU_CUSTOM_NONBONDED_FORCE_H__
