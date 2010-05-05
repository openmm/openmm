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
     * Get the array which contains random values.
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
private:
    OpenCLContext& context;
    cl::Kernel settleKernel;
    cl::Kernel shakeKernel;
    cl::Kernel ccmaDirectionsKernel;
    cl::Kernel ccmaForceKernel;
    cl::Kernel ccmaMultiplyKernel;
    cl::Kernel ccmaUpdateKernel;
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
    int randomPos;
    int lastSeed;
    bool hasInitializedConstraintKernels;
    struct ShakeCluster;
    struct ConstraintOrderer;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLINTEGRATIONUTILITIES_H_*/
