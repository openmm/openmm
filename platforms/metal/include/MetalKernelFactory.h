/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.                                    *
 * Source: platforms/opencl/include/OpenCLKernelFactory.h                    *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2008 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 *                                                                            *
 * Metal Platform code:                                                       *
 * Portions copyright (c) 2026 Chun-Chi Hung.                                 *
 * Authors: Chun-Chi Hung                                                     *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the               *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.       *
 * -------------------------------------------------------------------------- */

#ifndef OPENMM_METALKERNELFACTORY_H_
#define OPENMM_METALKERNELFACTORY_H_

#include "openmm/KernelFactory.h"

namespace OpenMM {

/** @brief Instantiate Common Compute kernels and the small Metal-specific adapters. */
class MetalKernelFactory : public KernelFactory {
public:
    /**
     * @brief Create a core kernel for a single Metal compute context.
     * @param name Kernel interface name.
     * @param platform Owning Metal Platform.
     * @param context Simulation Context whose device resources the kernel uses.
     * @return A new kernel owned by OpenMM.
     * @throws OpenMMException If the kernel name is not registered.
     */
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const override;
};

} // namespace OpenMM

#endif // OPENMM_METALKERNELFACTORY_H_
