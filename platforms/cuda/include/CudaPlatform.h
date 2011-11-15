#ifndef OPENMM_CUDAPLATFORM_H_
#define OPENMM_CUDAPLATFORM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
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

#include "openmm/Platform.h"
#include "windowsExportCuda.h"

struct _gpuContext;

namespace OpenMM {

/**
 * This Platform subclass uses CUDA implementations of the OpenMM kernels to run on NVidia GPUs.
 */

class OPENMMCUDA_EXPORT CudaPlatform : public Platform {
public:
    class PlatformData;
    CudaPlatform();
    const std::string& getName() const {
        static const std::string name = "Cuda";
        return name;
    }
    double getSpeed() const {
        return 50;
    }
    bool supportsDoublePrecision() const;
    const std::string& getPropertyValue(const Context& context, const std::string& property) const;
    void setPropertyValue(Context& context, const std::string& property, const std::string& value) const;
    void contextCreated(ContextImpl& context, const std::map<std::string, std::string>& properties) const;
    void contextDestroyed(ContextImpl& context) const;
    /**
     * This is the name of the parameter for selecting which CUDA device to use.
     */
    static const std::string& CudaDevice() {
        static const std::string key = "CudaDevice";
        return key;
    }
    /**
     * This is the name of the parameter for selecting whether CUDA should sync or spin loop while waiting for results.
     */
    static const std::string& CudaUseBlockingSync() {
        static const std::string key = "CudaUseBlockingSync";
        return key;
    }
};

class CudaPlatform::PlatformData {
public:
    OPENMMCUDA_EXPORT PlatformData(_gpuContext* gpu);
    _gpuContext* gpu;
    bool removeCM;
    bool hasBonds, hasAngles, hasPeriodicTorsions, hasRB, hasNonbonded, hasCustomNonbonded;
    int nonbondedMethod, customNonbondedMethod;
    int cmMotionFrequency;
    int stepCount, computeForceCount;
    double time, ewaldSelfEnergy, dispersionCoefficient;
    std::map<std::string, std::string> propertyValues;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPLATFORM_H_*/
