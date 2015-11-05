#ifndef OPENMM_COMPOUND_INTEGRATOR_PROXY_H_
#define OPENMM_COMPOUND_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class CompoundIntegratorProxy : public SerializationProxy {
public:
    CompoundIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_COMPOUND_INTEGRATOR_PROXY_H_*/