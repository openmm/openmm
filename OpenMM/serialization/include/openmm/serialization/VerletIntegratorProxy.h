#ifndef OPENMM_VERLET_INTEGRATOR_PROXY_H_
#define OPENMM_VERLET_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class VerletIntegratorProxy : public SerializationProxy {
public:
    VerletIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_VERLET_INTEGRATOR_PROXY_H_*/