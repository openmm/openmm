#ifndef OPENMM_OPENCLINTEGRATIONUTILITIES_H_
#define OPENMM_OPENCLINTEGRATIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/System.h"
#include "OpenCLContext.h"
#include "openmm/internal/windowsExport.h"
#include <iosfwd>

namespace OpenMM {

/**
 * This class implements features that are used by many different integrators, including
 * common workspace arrays, random number generation, and enforcing constraints.
 */

class OPENMM_EXPORT OpenCLIntegrationUtilities {
public:
    OpenCLIntegrationUtilities(OpenCLContext& context, const System& system);
    ~OpenCLIntegrationUtilities();
    /**
     * Get the array which contains position deltas.
     */
    OpenCLArray<mm_float4>& getPosDelta() {
        return *posDelta;
    }
    /**
     * Get the array which contains random values.  Each element is a float4, whose components
     * are independent, normally distributed random numbers with mean 0 and variance 1.
     */
    OpenCLArray<mm_float4>& getRandom() {
        return *random;
    }
    /**
     * Get the array which contains the current step size.
     */
    OpenCLArray<mm_float2>& getStepSize() {
        return *stepSize;
    }
    /**
     * Apply constraints to the atom positions.
     *
     * @param tol             the constraint tolerance
     */
    void applyConstraints(double tol);
    /**
     * Apply constraints to the atom velocities.
     *
     * @param tol             the constraint tolerance
     */
    void applyVelocityConstraints(double tol);
    /**
     * Initialize the random number generator.
     */
    void initRandomNumberGenerator(unsigned int randomNumberSeed);
    /**
     * Ensure that sufficient random numbers are available in the array, and generate new ones if not.
     *
     * @param numValues     the number of random float4's that will be required
     * @return the index in the array at which to start reading
     */
    int prepareRandomNumbers(int numValues);
    /**
     * Compute the positions of virtual sites.
     */
    void computeVirtualSites();
    /**
     * Distribute forces from virtual sites to the atoms they are based on.
     */
    void distributeForcesFromVirtualSites();
    /**
     * Create a checkpoint recording the current state of the random number generator.
     * 
     * @param stream    an output stream the checkpoint data should be written to
     */
    void createCheckpoint(std::ostream& stream);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     * 
     * @param stream    an input stream the checkpoint data should be read from
     */
    void loadCheckpoint(std::istream& stream);
private:
    void applyConstraints(bool constrainVelocities, double tol);
    OpenCLContext& context;
    cl::Kernel settlePosKernel, settleVelKernel;
    cl::Kernel shakePosKernel, shakeVelKernel;
    cl::Kernel ccmaDirectionsKernel;
    cl::Kernel ccmaPosForceKernel, ccmaVelForceKernel;
    cl::Kernel ccmaMultiplyKernel;
    cl::Kernel ccmaPosUpdateKernel, ccmaVelUpdateKernel;
    cl::Kernel vsitePositionKernel, vsiteForceKernel;
    cl::Kernel randomKernel;
    OpenCLArray<mm_float4>* posDelta;
    OpenCLArray<mm_int4>* settleAtoms;
    OpenCLArray<mm_float2>* settleParams;
    OpenCLArray<mm_int4>* shakeAtoms;
    OpenCLArray<mm_float4>* shakeParams;
    OpenCLArray<mm_float4>* random;
    OpenCLArray<mm_int4>* randomSeed;
    OpenCLArray<mm_float2>* stepSize;
    OpenCLArray<mm_int2>* ccmaAtoms;
    OpenCLArray<mm_float4>* ccmaDistance;
    OpenCLArray<cl_float>* ccmaReducedMass;
    OpenCLArray<cl_int>* ccmaAtomConstraints;
    OpenCLArray<cl_int>* ccmaNumAtomConstraints;
    OpenCLArray<cl_int>* ccmaConstraintMatrixColumn;
    OpenCLArray<cl_float>* ccmaConstraintMatrixValue;
    OpenCLArray<cl_float>* ccmaDelta1;
    OpenCLArray<cl_float>* ccmaDelta2;
    OpenCLArray<cl_int>* ccmaConverged;
    cl::Buffer* ccmaConvergedBuffer;
    cl_int* ccmaConvergedMemory;
    OpenCLArray<mm_int4>* vsite2AvgAtoms;
    OpenCLArray<mm_float2>* vsite2AvgWeights;
    OpenCLArray<mm_int4>* vsite3AvgAtoms;
    OpenCLArray<mm_float4>* vsite3AvgWeights;
    OpenCLArray<mm_int4>* vsiteOutOfPlaneAtoms;
    OpenCLArray<mm_float4>* vsiteOutOfPlaneWeights;
    int randomPos;
    int lastSeed, numVsites;
    bool hasInitializedPosConstraintKernels, hasInitializedVelConstraintKernels;
    struct ShakeCluster;
    struct ConstraintOrderer;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLINTEGRATIONUTILITIES_H_*/
