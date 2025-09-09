#ifndef OPENMM_QTB_INTEGRATOR_PROXY_H_
#define OPENMM_QTB_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

class QTBIntegratorProxy : public SerializationProxy {
public:
    QTBIntegratorProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif /*OPENMM_QTB_INTEGRATOR_PROXY_H_*/