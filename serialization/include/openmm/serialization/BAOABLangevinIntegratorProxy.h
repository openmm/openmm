#ifndef OPENMM_BAOAB_LANGEVIN_INTEGRATOR_PROXY_H_
#define OPENMM_BAOAB_LANGEVIN_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class BAOABLangevinIntegratorProxy : public SerializationProxy {
public:
    BAOABLangevinIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_BAOAB_LANGEVIN_INTEGRATOR_PROXY_H_*/