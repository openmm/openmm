#ifndef OPENMM_AMOEBA_ANGLE_TORSION_FORCE_PROXY_H_
#define OPENMM_AMOEBA_ANGLE_TORSION_FORCE_PROXY_H_

#include "openmm/internal/windowsExportAmoeba.h"
#include "openmm/serialization/SerializationProxy.h"

namespace OpenMM {

class OPENMM_EXPORT_AMOEBA AmoebaAngleTorsionForceProxy : public SerializationProxy {
public:
	AmoebaAngleTorsionForceProxy();
	void serialize(const void* object, SerializationNode& node) const;
	void* deserialize(const SerializationNode& node) const;
};

} // namespace 

#endif
