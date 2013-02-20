#ifndef OPENMM_LANGEVIN_INTEGRATOR_PROXY_H_
#define OPENMM_LANGEVIN_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class LangevinIntegratorProxy : public SerializationProxy {
public:
    LangevinIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_LANGEVIN_INTEGRATOR_PROXY_H_*/