#ifndef OPENMM_LANGEVIN_MIDDLE_INTEGRATOR_PROXY_H_
#define OPENMM_LANGEVIN_MIDDLE_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class LangevinMiddleIntegratorProxy : public SerializationProxy {
public:
    LangevinMiddleIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_LANGEVIN_MIDDLE_INTEGRATOR_PROXY_H_*/