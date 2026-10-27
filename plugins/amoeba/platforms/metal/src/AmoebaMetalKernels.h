#ifndef AMOEBA_OPENMM_METALKERNELS_H_
#define AMOEBA_OPENMM_METALKERNELS_H_

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
 * Authors: Mark Friedrichs, Peter Eastman                                    *
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

#include "openmm/amoebaKernels.h"
#include "openmm/kernels.h"
#include "openmm/System.h"
#include "openmm/common/ComputeSort.h"
#include "MetalContext.h"
#include "MetalNonbondedUtilities.h"
#include "AmoebaCommonKernels.h"

namespace OpenMM {

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class MetalCalcAmoebaMultipoleForceKernel : public CommonCalcAmoebaMultipoleForceKernel {
public:
    MetalCalcAmoebaMultipoleForceKernel(const std::string& name, const Platform& platform, MetalContext& mc, const System& system) :
            CommonCalcAmoebaMultipoleForceKernel(name, platform, mc, system) {
    }
    /**
     * Get whether charge spreading should be done in fixed point.
     */
    bool useFixedPointChargeSpreading() const {
        return cc.getUseDoublePrecision();
    }
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_METALKERNELS_H_*/
