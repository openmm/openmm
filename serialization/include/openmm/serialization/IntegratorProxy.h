#ifndef OPENMM_INTEGRATOR_PROXY_H_
#define OPENMM_INTEGRATOR_PROXY_H_

#include "openmm/internal/windowsExport.h"
#include "openmm/serialization/SerializationProxy.h"

namespace OpenMM {

/**
 * This is a proxy for serializing generic Integrator objects.
 * It makes calls to the serialize/deserialize methods in
 * derived classes of Integrator.
 */

class OPENMM_EXPORT IntegratorProxy : public SerializationProxy {
public:
    IntegratorProxy();
    virtual void serialize(const void* object, SerializationNode& node) const;
    virtual void* deserialize(const SerializationNode& node) const;
};

} // namespace OpenMM

#endif /*OPENMM_INTEGRATOR_PROXY_H_*/
