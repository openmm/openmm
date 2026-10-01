#ifndef OPENMM_METALINTEGRATIONUTILITIES_H_
#define OPENMM_METALINTEGRATIONUTILITIES_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2009-2026 Stanford University and the Authors.      *
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

#include "MetalArray.h"
#include "openmm/System.h"
#include "openmm/common/IntegrationUtilities.h"

namespace OpenMM {

class MetalContext;

/**
 * This class implements features that are used by many different integrators, including
 * common workspace arrays, random number generation, and enforcing constraints.
 */

class MetalIntegrationUtilities : public IntegrationUtilities {
public:
    MetalIntegrationUtilities(MetalContext& context, const System& system);
    ~MetalIntegrationUtilities();
    /**
     * Distribute forces from virtual sites to the atoms they are based on.
     */
    void distributeForcesFromVirtualSites();
private:
    void applyConstraintsImpl(bool constrainVelocities, double tol);
    int* ccmaConvergedMemory;
//    CUdeviceptr ccmaConvergedDeviceMemory;
//    CUevent ccmaEvent;
};

} // namespace OpenMM

#endif /*OPENMM_METALINTEGRATIONUTILITIES_H_*/
