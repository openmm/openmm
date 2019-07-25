#ifndef STATE_PROXY_H_
#define STATE_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM { // needs to be for friend class to work

class StateProxy : public SerializationProxy {
public:
    StateProxy();
    void serialize(const void* object, SerializationNode& node) const;
    void* deserialize(const SerializationNode& node) const;
};

}

#endif