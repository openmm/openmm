#ifndef OPENMM_BROWNIAN_INTEGRATOR_PROXY_H_
#define OPENMM_BROWNIAN_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

    class BrownianIntegratorProxy : public SerializationProxy {
    public:
        BrownianIntegratorProxy();
        void serialize(const void* object, SerializationNode& node) const;
        void* deserialize(const SerializationNode& node) const;
    };

}

#endif /*OPENMM_BROWNIAN_INTEGRATOR_PROXY_H_*/