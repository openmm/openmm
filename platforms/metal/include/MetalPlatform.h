/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Ported from the OpenMM OpenCL Platform.                                    *
 * Source: platforms/opencl/include/OpenCLPlatform.h                         *
 *                                                                            *
 * Original OpenCL Platform code:                                             *
 * Portions copyright (c) 2008-2026 Stanford University and the Authors.      *
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

#ifndef OPENMM_METALPLATFORM_H_
#define OPENMM_METALPLATFORM_H_

#include "openmm/Platform.h"
#include "openmm/common/windowsExportCommon.h"
#include <memory>

namespace OpenMM {

class MetalContext;

/**
 * @brief Single-GPU, single-precision Metal implementation of the core Platform.
 *
 * The host kernels use Common Compute, with the OpenCL backend as their
 * behavioral reference. CPU PME and multiple-device execution are not supported.
 */
class OPENMM_EXPORT_COMMON MetalPlatform : public Platform {
public:
    class PlatformData;
    /**
     * @brief Register factories for the core simulation kernels.
     * @param floatingAccumulators Internal default for floating-accumulator regression Contexts.
     *                             This is not a public Context property.
     */
    explicit MetalPlatform(bool floatingAccumulators=false);
    /** @return The stable Platform name, "Metal". */
    const std::string& getName() const override;
    /** @return A relative speed ranking, not a measured benchmark result. */
    double getSpeed() const override;
    /** @return False; this backend accepts only single precision. */
    bool supportsDoublePrecision() const override;
    /** @return True if the default device is a supported Apple silicon GPU. */
    static bool isPlatformSupported();
    /** @return The selected device and precision properties of a Context. */
    const std::string& getPropertyValue(const Context& context, const std::string& property) const override;
    /**
     * @brief Reject changes to properties fixed at Context creation.
     * @throws OpenMMException If the property or value differs from the current setting.
     */
    void setPropertyValue(Context& context, const std::string& property, const std::string& value) const override;
    /**
     * @brief Report the single supported default GPU when it matches all filters.
     * @param filters Optional DeviceIndex, DeviceName, Precision, and UseCpuPme values.
     * @return Zero or one device description; unsupported filter values have no match.
     */
    std::vector<std::map<std::string, std::string> > getDevices(
            const std::map<std::string, std::string>& filters={}) const override;
    /** @brief Validate properties and allocate the Context's Metal resources. */
    void contextCreated(ContextImpl& context, const std::map<std::string, std::string>& properties) const override;
    /** @brief Create an inner Context sharing the original GPU command queue. */
    void linkedContextCreated(ContextImpl& context, ContextImpl& originalContext) const override;
    /** @brief Release the resources owned by a simulation Context. */
    void contextDestroyed(ContextImpl& context) const override;
    /** @return Device-selection property; only the default device (index 0) is supported. */
    static const std::string& MetalDeviceIndex();
    /** @return The read-only name of the selected GPU. */
    static const std::string& MetalDeviceName();
    /** @return Precision-selection property; its only supported value is "single". */
    static const std::string& MetalPrecision();
    /** @return CPU PME property; its only supported value is "false". */
    static const std::string& MetalUseCpuPme();
private:
    bool floatingAccumulators; ///< Default accumulator mode; not a public Context property.
};

/** @brief Resources associated with one simulation Context, including inner Contexts. */
class OPENMM_EXPORT_COMMON MetalPlatform::PlatformData {
public:
    /** @brief Allocate the single compute context after property validation. */
    explicit PlatformData(ContextImpl& context, MetalContext* linked=nullptr, bool floatingAccumulators=false);
    /** @brief Finish and release the compute context's owned resources. */
    ~PlatformData();
    ContextImpl* context; ///< Borrowed simulation Context, valid for this object's lifetime.
    std::unique_ptr<MetalContext> computeContext; ///< Sole owned GPU compute context.
    std::map<std::string, std::string> propertyValues; ///< Resolved, immutable properties.
};

} // namespace OpenMM

#endif // OPENMM_METALPLATFORM_H_
