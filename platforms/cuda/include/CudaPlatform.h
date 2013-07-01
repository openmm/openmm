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
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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
#include "openmm/System.h"
#include "windowsExportCuda.h"

namespace OpenMM {
    
class CudaContext;

/**
 * This Platform subclass uses CUDA implementations of the OpenMM kernels.
 */

class OPENMM_EXPORT_CUDA CudaPlatform : public Platform {
public:
    class PlatformData;
    CudaPlatform();
    const std::string& getName() const {
        static const std::string name = "CUDA";
        return name;
    }
    double getSpeed() const;
    bool supportsDoublePrecision() const;
    const std::string& getPropertyValue(const Context& context, const std::string& property) const;
    void setPropertyValue(Context& context, const std::string& property, const std::string& value) const;
    void contextCreated(ContextImpl& context, const std::map<std::string, std::string>& properties) const;
    void contextDestroyed(ContextImpl& context) const;
    /**
     * This is the name of the parameter for selecting which CUDA device or devices to use.
     */
    static const std::string& CudaDeviceIndex() {
        static const std::string key = "CudaDeviceIndex";
        return key;
    }
    /**
     * This is the name of the parameter that reports the CUDA device or devices being used.
     */
    static const std::string& CudaDeviceName() {
        static const std::string key = "CudaDeviceName";
        return key;
    }
    /**
     * This is the name of the parameter for selecting whether CUDA should sync or spin loop while waiting for results.
     */
    static const std::string& CudaUseBlockingSync() {
        static const std::string key = "CudaUseBlockingSync";
        return key;
    }
    /**
     * This is the name of the parameter for selecting what numerical precision to use.
     */
    static const std::string& CudaPrecision() {
        static const std::string key = "CudaPrecision";
        return key;
    }
    /**
     * This is the name of the parameter for selecting whether to use the CPU based PME calculation.
     */
    static const std::string& CudaUseCpuPme() {
        static const std::string key = "CudaUseCpuPme";
        return key;
    }
    /**
     * This is the name of the parameter for specifying the path to the CUDA compiler.
     */
    static const std::string& CudaCompiler() {
        static const std::string key = "CudaCompiler";
        return key;
    }
    /**
     * This is the name of the parameter for specifying the path to the directory for creating temporary files.
     */
    static const std::string& CudaTempDirectory() {
        static const std::string key = "CudaTempDirectory";
        return key;
    }
};

class OPENMM_EXPORT_CUDA CudaPlatform::PlatformData {
public:
    PlatformData(ContextImpl* context, const System& system, const std::string& deviceIndexProperty, const std::string& blockingProperty, const std::string& precisionProperty,
            const std::string& cpuPmeProperty, const std::string& compilerProperty, const std::string& tempProperty);
    ~PlatformData();
    void initializeContexts(const System& system);
    void syncContexts();
    ContextImpl* context;
    std::vector<CudaContext*> contexts;
    std::vector<double> contextEnergy;
    bool removeCM, peerAccessSupported, useCpuPme;
    int cmMotionFrequency;
    int stepCount, computeForceCount;
    double time;
    std::map<std::string, std::string> propertyValues;
};

} // namespace OpenMM

#endif /*OPENMM_CUDAPLATFORM_H_*/
