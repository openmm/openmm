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
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
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

namespace OpenMM {

/**
 * XmlSerializer is used for serializing objects as XML, and for reconstructing them again.
 */

class OPENMM_EXPORT XmlSerializer {
public:
    /**
     * Serialize an object as XML.
     *
     * @param object    the object to serialize
     * @param rootName  the name to use for the root node of the XML document
     * @param stream    an output stream to write the XML to
     */
    template <class T>
    static void serialize(const T* object, const std::string& rootName, std::ostream& stream) {
        const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(*object));
        SerializationNode node;
        node.setName(rootName);
        proxy.serialize(object, node);
        if (node.hasProperty("type"))
            throw OpenMMException(proxy.getTypeName()+" created node with reserved property 'type'");
        node.setStringProperty("type", proxy.getTypeName());
        serialize(node, stream);
    }
    /**
     * Reconstruct an object that has been serialized as XML.
     *
     * @param stream    an input stream to read the XML from
     * @return a pointer to the newly created object.  The caller assumes ownership of the object.
     */
    template <class T>
    static T* deserialize(std::istream& stream) {
        return reinterpret_cast<T*>(deserializeStream(stream));
    }
    /**
     * Clone an object by first serializing it, then deserializing it again.  This method constructs the
     * new object directly from the SerializationNodes without first converting them to XML.  This means
     * it is faster and uses less memory than making separate calls to serialize() and deserialize().
     */
    template <class T>
    static T* clone(const T& object) {
        const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(object));
        SerializationNode node;
        proxy.serialize(&object, node);
        return reinterpret_cast<T*>(proxy.deserialize(node));
    }
private:
    class StreamReader;
    static void serialize(const SerializationNode& node, std::ostream& stream);
    static void* deserializeStream(std::istream& stream);
    static void encodeNode(const SerializationNode& node, std::ostream& stream, int depth);
};

} // namespace OpenMM

#endif /*OPENMM_XML_SERIALIZER_H_*/
