#ifndef OPENMM_DPD_INTEGRATOR_PROXY_H_
#define OPENMM_DPD_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class DPDIntegratorProxy : public SerializationProxy {
public:
    DPDIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_DPD_INTEGRATOR_PROXY_H_*/