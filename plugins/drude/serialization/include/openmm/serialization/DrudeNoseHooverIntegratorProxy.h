#ifndef OPENMM_DRUDE_NOSE_HOOVER_INTEGRATOR_PROXY_H_
#define OPENMM_DRUDE_NOSE_HOOVER_INTEGRATOR_PROXY_H_

#include "openmm/serialization/SerializationProxy.h"
#include "openmm/internal/windowsExportDrude.h"

namespace OpenMM {

class OPENMM_EXPORT_DRUDE DrudeNoseHooverIntegratorProxy : public SerializationProxy {
public:
    DrudeNoseHooverIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_DRUDE_NOSE_HOOVER_INTEGRATOR_PROXY_H_*/
