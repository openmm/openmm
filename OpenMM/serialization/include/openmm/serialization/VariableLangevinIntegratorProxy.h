#ifndef OPENMM_VARIABLE_LANGEVIN_INTEGRATOR_PROXY_H_
#define OPENMM_VARIABLE_LANGEVIN_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

    class VariableLangevinIntegratorProxy : public SerializationProxy {
    public:
        VariableLangevinIntegratorProxy();
        void serialize(const void* object, SerializationNode& node) const;
        void* deserialize(const SerializationNode& node) const;
    };

}

#endif /*OPENMM_VARIABLE_LANGEVIN_INTEGRATOR_PROXY_H_*/