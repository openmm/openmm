#ifndef OPENMM_NOSE_HOOVER_INTEGRATOR_PROXY_H_
#define OPENMM_NOSE_HOOVER_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class NoseHooverIntegratorProxy : public SerializationProxy {
public:
    NoseHooverIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_NOSE_HOOVER_INTEGRATOR_PROXY_H_*/
