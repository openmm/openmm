#ifndef AMOEBA_OPENMM_OPENCLKERNELS_H
#define AMOEBA_OPENMM_OPENCLKERNELS_H

/* -------------------------------------------------------------------------- *
 *                              OpenMMAmoeba                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2025 Stanford University and the Authors.      *
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
#include "AmoebaCommonKernels.h"
#include "OpenCLContext.h"

namespace OpenMM {

/**
 * This kernel is invoked by AmoebaMultipoleForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcAmoebaMultipoleForceKernel : public CommonCalcAmoebaMultipoleForceKernel {
public:
    OpenCLCalcAmoebaMultipoleForceKernel(const std::string& name, const Platform& platform, OpenCLContext& cl, const System& system) :
            CommonCalcAmoebaMultipoleForceKernel(name, platform, cl, system) {
    }
    /**
     * Get whether charge spreading should be done in fixed point.
     */
    bool useFixedPointChargeSpreading() const {
        return true;
    }
};


/**
 * This kernel is invoked by HippoNonbondedForce to calculate the forces acting on the system and the energy of the system.
 */
class OpenCLCalcHippoNonbondedForceKernel : public CommonCalcHippoNonbondedForceKernel {
public:
    OpenCLCalcHippoNonbondedForceKernel(const std::string& name, const Platform& platform, OpenCLContext& cl, const System& system) :
            CommonCalcHippoNonbondedForceKernel(name, platform, cl, system) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the HippoNonbondedForce this kernel will be used for
     */
    void initialize(const System& system, const HippoNonbondedForce& force);
    /**
     * Get whether charge spreading should be done in fixed point.
     */
    bool useFixedPointChargeSpreading() const {
        return true;
    }
    /**
     * Sort the atom grid indices.
     */
    void sortGridIndex();
private:
    class SortTrait : public ComputeSortImpl::SortTrait {
        int getDataSize() const {return 8;}
        int getKeySize() const {return 4;}
        const char* getDataType() const {return "int2";}
        const char* getKeyType() const {return "int";}
        const char* getMinKey() const {return "(-2147483647-1)";}
        const char* getMaxKey() const {return "2147483647";}
        const char* getMaxValue() const {return "make_int2(2147483647, 2147483647)";}
        const char* getSortKey() const {return "value.y";}
    };
    ComputeSort sort;
};

} // namespace OpenMM

#endif /*AMOEBA_OPENMM_OPENCLKERNELS_H*/
