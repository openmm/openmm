#ifndef OPENMM_CPULCPOFORCE_H_
#define OPENMM_CPULCPOFORCE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2025 Stanford University and the Authors.           *
 * Authors: Evan Pretti                                                       *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include <map>
#include <vector>

#include "CpuNeighborList.h"
#include "openmm/internal/vectorize.h"

namespace OpenMM {

/**
 * Performs energy, force, and charge derivative calculations for the CPU kernel
 * for LCPOForce.
 */
class CpuLCPOForce {
private:
    /**
     * Holds data from a thread's request to record a neighbor.
     */
    struct NeighborInfo {
        int i, j;
        fvec4 ijData, jiData;
    };

    /**
     * Keeps track of neighbors for CpuLCPOForce.
     */
    class Neighbors {
    private:
        int numParticles;
        int maxNumNeighbors;
        std::vector<int> numNeighbors;
        std::vector<int> indices;
        std::vector<fvec4> data;

    public:
        /**
         * Creates an empty CpuLCPOForce::Neighbors.
         *
         * @param numParticles         the number of particles to track (this will be numActiveParticles of the parent CpuLCPOForce)
         * @param maxNumNeighborsInit  an initial guess for the maximum number of neighbors per particle
         */
        Neighbors(int numParticles, int maxNumNeighborsInit = 0);

        /**
         * Clears all neighbor data.
         */
        void clear();

        /**
         * Records that two particles are neighbors of each other.
         *
         * @param i       the index of the first particle
         * @param j       the index of the second particle
         * @param ijData  the data to record for particle j as a neighbor of particle i
         * @param jiData  the data to record for particle i as a neighbor of particle j
         */
        void insert(int i, int j, fvec4 ijData, fvec4 jiData);

        /**
         * Retrieves the neighbors of a particle.
         *
         * @param i              the index of the particle
         * @param iNumNeighbors  the number of neighbors of the particle
         * @param iIndices       pointer to the start of the indices of the neighbors
         * @param iData          pointer to the start of the data associated with the neighbors
         */
        void getNeighbors(int i, int& iNumNeighbors, const int*& iIndices, const fvec4*& iData) const;

        /**
         * Tests whether a particle is a neighbor of another.
         *
         * @param i  the index of the first particle
         * @param j  the index of the second particle
         * @return   whether or not the particles are neighbors
         */
        bool isNeighbor(int i, int j) const;

    private:
        void insert(int i, int j, fvec4 ijData);
        void resize(int newMaxNumNeighbors);
    };

public:
    /**
     * Creates a CpuLCPOForce.
     *
     * @param threads      thread pool for computations
     * @param indices      active particle indices for each particle (-1 if not an active particle)
     * @param particles    system particle indices for each active particle
     * @param parameters   parameters (radius, P2, P3, and P4) for each active particle
     * @param usePeriodic  whether or not to use periodic boundary conditions
     */
    CpuLCPOForce(ThreadPool& threads, const std::vector<int>& indices, const std::vector<int>& particles, const std::vector<fvec4>& parameters, bool usePeriodic);

    ~CpuLCPOForce();

    /**
     * Indicates that the radii of some particles may have changed and that the
     * maximum required cutoff distance should be updated.
     */
    void updateRadii();

    /**
     * Computes energies and forces.
     *
     * @param boxVectors     periodic box vectors
     * @param posq           particle positions
     * @param threadForce    thread force accumulators
     * @param includeForces  whether or not to calculate interaction forces
     * @param includeEnergy  whether or not to calculate interaction energies
     * @param energy         energy accumulator
     */
    void execute(Vec3* boxVectors, AlignedArray<float>& posq, std::vector<AlignedArray<float> >& threadForce, bool includeForces, bool includeEnergy, double& energy);

private:
    /**
     * Thread worker for computing energies and forces.
     */
    void threadExecute(ThreadPool& threads, int threadIndex);

    /**
     * Helper for processing a block of the neighbor list.
     */
    template<bool USE_PERIODIC>
    void processNeighborListBlock(int blockIndex, std::vector<NeighborInfo>& threadNeighborInfo);

private:
    static const int RadiusIndex = 0;
    static const int P2Index = 1;
    static const int P3Index = 2;
    static const int P4Index = 3;

    int numParticles;
    int numActiveParticles;

    ThreadPool& threads;
    const std::vector<int>& activeParticles;
    const std::vector<int>& activeParticlesInv;
    const std::vector<fvec4>& parameters;
    bool usePeriodic;

    CpuNeighborList* neighborList;
    std::vector<std::set<int> > exclusions;
    Neighbors neighbors;
    float cutoff;
    Vec3* boxVectors;
    float recipBoxSize[3];

    float* posq;
    std::vector<double> threadEnergy;
    std::vector<AlignedArray<float> >* threadForce;
    std::vector<std::vector<NeighborInfo> > threadNeighbors;
    std::atomic<int> atomicCounter;
};

} // namespace OpenMM

#endif // OPENMM_CPULCPOFORCE_H_