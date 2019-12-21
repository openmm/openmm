#ifndef OPENMM_CUDAINTEGRATIONUTILITIES_H_
#define OPENMM_CUDAINTEGRATIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2019 Stanford University and the Authors.      *
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

#include "CudaArray.h"
#include "openmm/System.h"
#include "openmm/common/IntegrationUtilities.h"
#include "openmm/common/windowsExportCommon.h"
#include <cuda.h>
#include <builtin_types.h>

namespace OpenMM {

class CudaContext;

/**
 * This class implements features that are used by many different integrators, including
 * common workspace arrays, random number generation, and enforcing constraints.
 */

class OPENMM_EXPORT_COMMON CudaIntegrationUtilities : public IntegrationUtilities {
public:
    CudaIntegrationUtilities(CudaContext& context, const System& system);
    ~CudaIntegrationUtilities();
    /**
     * Get the array which contains position deltas.
     */
    CudaArray& getPosDelta();
    /**
     * Get the array which contains random values.  Each element is a float4, whose components
     * are independent, normally distributed random numbers with mean 0 and variance 1.
     */
    CudaArray& getRandom();
    /**
     * Get the array which contains the current step size.
     */
    CudaArray& getStepSize();
    /**
     * Distribute forces from virtual sites to the atoms they are based on.
     */
    void distributeForcesFromVirtualSites();
private:
    void applyConstraintsImpl(bool constrainVelocities, double tol);
    int* ccmaConvergedMemory;
    CUdeviceptr ccmaConvergedDeviceMemory;
    CUevent ccmaEvent;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAINTEGRATIONUTILITIES_H_*/
