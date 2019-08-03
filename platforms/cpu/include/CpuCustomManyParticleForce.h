
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

#ifndef OPENMM_CPU_CUSTOM_MANY_PARTICLE_FORCE_H__
#define OPENMM_CPU_CUSTOM_MANY_PARTICLE_FORCE_H__

#include "ReferenceForce.h"
#include "ReferenceBondIxn.h"
#include "CpuNeighborList.h"
#include "openmm/CustomManyParticleForce.h"
#include "openmm/internal/CompiledExpressionSet.h"
#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/vectorize.h"
#include "lepton/CompiledExpression.h"
#include "lepton/ParsedExpression.h"
#include <atomic>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace OpenMM {

class CpuCustomManyParticleForce {
private:

    class ParticleTermInfo;
    class DistanceTermInfo;
    class AngleTermInfo;
    class DihedralTermInfo;
    class ThreadData;
    int numParticles, numParticlesPerSet, numPerParticleParameters, numTypes;
    bool useCutoff, usePeriodic, triclinic, centralParticleMode;
    double cutoffDistance;
    float recipBoxSize[3];
    Vec3 periodicBoxVectors[3];
    AlignedArray<fvec4> periodicBoxVec4;
    CpuNeighborList* neighborList;
    ThreadPool& threads;
    std::vector<std::set<int> > exclusions;
    std::vector<int> particleTypes;
    std::vector<int> orderIndex;
    std::vector<std::vector<int> > particleOrder;
    std::vector<std::vector<int> > particleNeighbors;
    std::vector<ThreadData*> threadData;
    // The following variables are used to make information accessible to the individual threads.
    float* posq;
    std::vector<double>* particleParameters;        
    const std::map<std::string, double>* globalParameters;
    std::vector<AlignedArray<float> >* threadForce;
    bool includeForces, includeEnergy;
    std::atomic<int> atomicCounter;

    /**
     * This routine contains the code executed by each thread.
     */
    void threadComputeForce(ThreadPool& threads, int threadIndex);

    /**
     * This is called recursively to loop over all possible combination of a set of particles and evaluate the
     * interaction for each one.
     */
    void loopOverInteractions(std::vector<int>& availableParticles, std::vector<int>& particleSet, int loopIndex, int startIndex,
                              std::vector<double>* particleParameters, float* forces, ThreadData& data, const fvec4& boxSize, const fvec4& invBoxSize);

    /**---------------------------------------------------------------------------------------

       Calculate custom interaction for one set of particles

       @param particleSet        the indices of the particles
       @param posq               atom coordinates in float format
       @param particleParameters particle parameter values (particleParameters[particleIndex][parameterIndex])
       @param forces             force array (forces added)
       @param totalEnergy        total energy

       --------------------------------------------------------------------------------------- */

    /**
     * Calculate the interaction for one set of particles
     * 
     * @param particleSet        the indices of the particles
     * @param particleParameters particle parameter values (particleParameters[particleIndex][parameterIndex])
     * @param data               information and workspace for the current thread
     * @param boxSize            the size of the periodic box
     * @param invBoxSize         the inverse size of the periodic box
     */
    void calculateOneIxn(std::vector<int>& particleSet, std::vector<double>* particleParameters, float* forces, ThreadData& data, const fvec4& boxSize, const fvec4& invBoxSize);

    /**
     * Compute the displacement and squared distance between two points, optionally using
     * periodic boundary conditions.
     */
    void computeDelta(const fvec4& posI, const fvec4& posJ, fvec4& deltaR, float& r2, const fvec4& boxSize, const fvec4& invBoxSize) const;
    
    static float computeAngle(const fvec4& vi, const fvec4& vj, float v2i, float v2j, float sign);
    
    static float getDihedralAngleBetweenThreeVectors(const fvec4& v1, const fvec4& v2, const fvec4& v3, fvec4& cross1, fvec4& cross2, const fvec4& signVector);

public:
    /**
     * Create a new CpuCustomManyParticleForce.
     *
     * @param force      the CustomManyParticleForce to create it for
     * @param threads    the thread pool to use
     */
    CpuCustomManyParticleForce(const OpenMM::CustomManyParticleForce& force, ThreadPool& threads);

    ~CpuCustomManyParticleForce();

    /**
     * Set the force to use a cutoff.
     * 
     * @param distance   the cutoff distance
     */
    void setUseCutoff(double distance);

    /**
     * Set the force to use periodic boundary conditions.  This requires that a cutoff has
     * already been set, and the smallest side of the periodic box is at least twice the cutoff
     * distance.
     * 
     * @param periodicBoxVectors    the vectors defining the periodic box
     */
    void setPeriodic(Vec3* periodicBoxVectors);

    /**
     * Calculate the interaction.
     * 
     * @param posq               atom coordinates in float format
     * @param particleParameters particle parameter values (particleParameters[particleIndex][parameterIndex])
     * @param globalParameters   the values of global parameters
     * @param threadForce        the collection of arrays for each thread to add forces to
     * @param includeForce       whether to compute forces
     * @param includeEnergy      whether to compute energy
     * @param energy             the total energy is added to this
     */
    void calculateIxn(AlignedArray<float>& posq, std::vector<std::vector<double> >& particleParameters, const std::map<std::string, double>& globalParameters,
                      std::vector<AlignedArray<float> >& threadForce, bool includeForces, bool includeEnergy, double& energy);
};

class CpuCustomManyParticleForce::ParticleTermInfo {
public:
    std::string name;
    int atom, component, variableIndex;
    Lepton::CompiledExpression forceExpression;
    ParticleTermInfo(const std::string& name, int atom, int component, const Lepton::CompiledExpression& forceExpression, ThreadData& data);
};

class CpuCustomManyParticleForce::DistanceTermInfo {
public:
    std::string name;
    int p1, p2, variableIndex;
    Lepton::CompiledExpression forceExpression;
    int delta;
    float deltaSign;
    DistanceTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::CompiledExpression& forceExpression, ThreadData& data);
};

class CpuCustomManyParticleForce::AngleTermInfo {
public:
    std::string name;
    int p1, p2, p3, variableIndex;
    Lepton::CompiledExpression forceExpression;
    int delta1, delta2;
    float delta1Sign, delta2Sign;
    AngleTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::CompiledExpression& forceExpression, ThreadData& data);
};

class CpuCustomManyParticleForce::DihedralTermInfo {
public:
    std::string name;
    int p1, p2, p3, p4, variableIndex;
    Lepton::CompiledExpression forceExpression;
    int delta1, delta2, delta3;
    DihedralTermInfo(const std::string& name, const std::vector<int>& atoms, const Lepton::CompiledExpression& forceExpression, ThreadData& data);
};

class CpuCustomManyParticleForce::ThreadData {
public:
    CompiledExpressionSet expressionSet;
    Lepton::CompiledExpression energyExpression;
    std::vector<std::vector<int> > particleParamIndices;
    std::vector<int> permutedParticles;
    std::vector<std::pair<int, int> > deltaPairs;
    std::vector<ParticleTermInfo> particleTerms;
    std::vector<DistanceTermInfo> distanceTerms;
    std::vector<AngleTermInfo> angleTerms;
    std::vector<DihedralTermInfo> dihedralTerms;
    AlignedArray<fvec4> delta, cross1, cross2;
    std::vector<float> normDelta;
    std::vector<float> norm2Delta;
    AlignedArray<fvec4> f;
    double energy;
    ThreadData(const CustomManyParticleForce& force, Lepton::ParsedExpression& energyExpr,
            std::map<std::string, std::vector<int> >& distances, std::map<std::string, std::vector<int> >& angles, std::map<std::string, std::vector<int> >& dihedrals);
    /**
     * Request a pair of particles whose distance or displacement vector is needed in the computation.
     */
    void requestDeltaPair(int p1, int p2, int& pairIndex, float& pairSign, bool allowReversed);
};

} // namespace OpenMM

#endif // OPENMM_CPU_CUSTOM_MANY_PARTICLE_FORCE_H__
