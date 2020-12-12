
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

#ifndef OPENMM_CPU_CUSTOM_NONBONDED_FORCE_H__
#define OPENMM_CPU_CUSTOM_NONBONDED_FORCE_H__

#include "AlignedArray.h"
#include "CpuNeighborList.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
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

       CpuCustomNonbondedForce(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpression,
                               const std::vector<std::string>& parameterNames, const std::vector<std::set<int> >& exclusions,
                               const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions, ThreadPool& threads);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~CpuCustomNonbondedForce();

      /**---------------------------------------------------------------------------------------

         Set the force to use a cutoff.

         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use

         --------------------------------------------------------------------------------------- */

      void setUseCutoff(double distance, const CpuNeighborList& neighbors);

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
private:
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
    const std::vector<std::set<int> > exclusions;
    std::vector<ThreadData*> threadData;
    std::vector<std::string> paramNames;
    std::vector<std::pair<int, int> > groupInteractions;
    std::vector<double> threadEnergy;
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
     * @param boxSize          the inverse size of the periodic box
     */
    void calculateOneIxn(int atom1, int atom2, ThreadData& data, float* forces, double& totalEnergy, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Compute the displacement and squared distance between two points, optionally using
     * periodic boundary conditions.
     */
    void getDeltaR(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const;
};

class CpuCustomNonbondedForce::ThreadData {
public:
    ThreadData(const Lepton::CompiledExpression& energyExpression, const Lepton::CompiledExpression& forceExpression, const std::vector<std::string>& parameterNames,
            const std::vector<Lepton::CompiledExpression> energyParamDerivExpressions);
    Lepton::CompiledExpression energyExpression;
    Lepton::CompiledExpression forceExpression;
    std::vector<Lepton::CompiledExpression> energyParamDerivExpressions;
    CompiledExpressionSet expressionSet;
    std::vector<double> particleParam;
    double r;
    std::vector<double> energyParamDerivs; 
};

} // namespace OpenMM

#endif // OPENMM_CPU_CUSTOM_NONBONDED_FORCE_H__
