#ifndef OPENMM_CUSTOM_INTEGRATOR_PROXY_H_
#define OPENMM_CUSTOM_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

    class CustomIntegratorProxy : public SerializationProxy {
    public:
        CustomIntegratorProxy();
        void serialize(const void* object, SerializationNode& node) const;
        void* deserialize(const SerializationNode& node) const;
    };

}

#endif /*OPENMM_CUSTOM_INTEGRATOR_PROXY_H_*/