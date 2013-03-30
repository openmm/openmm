#ifndef OPENMM_VARIABLE_VERLET_INTEGRATOR_PROXY_H_
#define OPENMM_VARIABLE_VERLET_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

    class VariableVerletIntegratorProxy : public SerializationProxy {
    public:
        VariableVerletIntegratorProxy();
        void serialize(const void* object, SerializationNode& node) const;
        void* deserialize(const SerializationNode& node) const;
    };

}

#endif /*OPENMM_VARIABLE_VERLET_INTEGRATOR_PROXY_H_*/