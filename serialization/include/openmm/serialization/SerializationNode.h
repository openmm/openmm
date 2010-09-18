#ifndef OPENMM_SERIALIZATIONNODE_H_
#define OPENMM_SERIALIZATIONNODE_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/serialization/SerializationProxy.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"
#include <map>
#include <string>
#include <vector>

namespace OpenMM {

/**
 */

class OPENMM_EXPORT SerializationNode {
public:
    const std::string& getName() const;
    void setName(const std::string& name);
    const std::vector<SerializationNode>& getChildren() const;
    std::vector<SerializationNode>& getChildren();
    const SerializationNode& getChildNode(const std::string& name) const;
    SerializationNode& getChildNode(const std::string& name);
    const std::map<std::string, std::string>& getProperties() const;
    bool hasProperty(const std::string& name) const;
    const std::string& getStringProperty(const std::string& name) const;
    const std::string& getStringProperty(const std::string& name, const std::string& defaultValue) const;
    SerializationNode& setStringProperty(const std::string& name, const std::string& value);
    int getIntProperty(const std::string& name) const;
    int getIntProperty(const std::string& name, int defaultValue) const;
    SerializationNode& setIntProperty(const std::string& name, int value);
    double getDoubleProperty(const std::string& name) const;
    double getDoubleProperty(const std::string& name, double defaultValue) const;
    SerializationNode& setDoubleProperty(const std::string& name, double value);
    SerializationNode& createChildNode(const std::string& name);
    template <class T>
    SerializationNode& createChildNode(const std::string& name, const T* object) {
        const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(*object));
        SerializationNode& node = createChildNode(name);
        proxy.serialize(object, node);
        if (node.hasProperty("type"))
            throw OpenMMException(proxy.getTypeName()+" created node with reserved property 'type'");
        node.setStringProperty("type", proxy.getTypeName());
        return node;
    }
    template<class T>
    T* decodeObject() const {
        return reinterpret_cast<T*>(SerializationProxy::getProxy(getStringProperty("type")).deserialize(*this));
    }
private:
    std::string name;
    std::vector<SerializationNode> children;
    std::map<std::string, std::string> properties;
};

} // namespace OpenMM

#endif /*OPENMM_SERIALIZATIONNODE_H_*/
