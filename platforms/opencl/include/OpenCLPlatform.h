#ifndef OPENMM_OPENCLPLATFORM_H_
#define OPENMM_OPENCLPLATFORM_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
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
#include "openmm/internal/ThreadPool.h"
#include "windowsExportOpenCL.h"
namespace OpenMM {
    
class OpenCLContext;

/**
 * This Platform subclass uses OpenCL implementations of the OpenMM kernels.
 */

class OPENMM_EXPORT_OPENCL OpenCLPlatform : public Platform {
public:
    class PlatformData;
    OpenCLPlatform();
    const std::string& getName() const {
        static const std::string name = "OpenCL";
        return name;
    }
    double getSpeed() const;
    bool supportsDoublePrecision() const;
    static bool isPlatformSupported();
    const std::string& getPropertyValue(const Context& context, const std::string& property) const;
    void setPropertyValue(Context& context, const std::string& property, const std::string& value) const;
    void contextCreated(ContextImpl& context, const std::map<std::string, std::string>& properties) const;
    void linkedContextCreated(ContextImpl& context, ContextImpl& originalContext) const;
    void contextDestroyed(ContextImpl& context) const;
    /**
     * This is the name of the parameter for selecting which OpenCL device or devices to use.
     */
    static const std::string& OpenCLDeviceIndex() {
        static const std::string key = "DeviceIndex";
        return key;
    }
    /**
     * This is the name of the parameter that reports the OpenCL device or devices being used.
     */
    static const std::string& OpenCLDeviceName() {
        static const std::string key = "DeviceName";
        return key;
    }
    /**
     * This is the name of the parameter for selecting which OpenCL platform to use.
     */
    static const std::string& OpenCLPlatformIndex() {
        static const std::string key = "OpenCLPlatformIndex";
        return key;
    }
    /**
     * This is the name of the parameter that reports the OpenCL platform being used.
     */
    static const std::string& OpenCLPlatformName() {
        static const std::string key = "OpenCLPlatformName";
        return key;
    }
    /**
     * This is the name of the parameter for selecting what numerical precision to use.
     */
    static const std::string& OpenCLPrecision() {
        static const std::string key = "Precision";
        return key;
    }
    /**
     * This is the name of the parameter for selecting whether to use the CPU based PME calculation.
     */
    static const std::string& OpenCLUseCpuPme() {
        static const std::string key = "UseCpuPme";
        return key;
    }
    /**
     * This is the name of the parameter for selecting whether to disable use of a separate stream for PME.
     */
    static const std::string& OpenCLDisablePmeStream() {
        static const std::string key = "DisablePmeStream";
        return key;
    }
};

class OPENMM_EXPORT_OPENCL OpenCLPlatform::PlatformData {
public:
    PlatformData(const System& system, const std::string& platformPropValue, const std::string& deviceIndexProperty, const std::string& precisionProperty,
            const std::string& cpuPmeProperty, const std::string& pmeStreamProperty, int numThreads, ContextImpl* originalContext);
    ~PlatformData();
    void initializeContexts(const System& system);
    void syncContexts();
    ContextImpl* context;
    std::vector<OpenCLContext*> contexts;
    std::vector<double> contextEnergy;
    bool hasInitializedContexts, removeCM, useCpuPme, disablePmeStream;
    int cmMotionFrequency;
    int stepCount, computeForceCount;
    double time;
    std::map<std::string, std::string> propertyValues;
    ThreadPool threads;
};

} // namespace OpenMM

#endif /*OPENMM_OPENCLPLATFORM_H_*/
