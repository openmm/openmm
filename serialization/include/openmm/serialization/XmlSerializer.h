#ifndef OPENMM_XML_SERIALIZER_H_
#define OPENMM_XML_SERIALIZER_H_

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

#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/windowsExport.h"
#include <iosfwd>

class TiXmlElement;

namespace OpenMM {

/**
 */

class OPENMM_EXPORT XmlSerializer {
public:
    template <class T>
    static void serialize(const T* object, const std::string& rootName, std::ostream& stream) {
        const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(*object));
        SerializationNode node;
        node.setName(rootName);
        proxy.serialize(object, node);
        if (node.hasProperty("type"))
            throw OpenMMException(proxy.getTypeName()+" created node with reserved property 'type'");
        node.setStringProperty("type", proxy.getTypeName());
        return serialize(node, stream);
    }
    template <class T>
    static T* deserialize(std::istream& stream) {
        return reinterpret_cast<T*>(deserializeStream(stream));
    }
private:
    static void serialize(const SerializationNode& node, std::ostream& stream);
    static void* deserializeStream(std::istream& stream);
    static TiXmlElement* encodeNode(const SerializationNode& node);
    static SerializationNode* decodeNode(SerializationNode& node, const TiXmlElement& element);
};

} // namespace OpenMM

#endif /*OPENMM_XML_SERIALIZER_H_*/
