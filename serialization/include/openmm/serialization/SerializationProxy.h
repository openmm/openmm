#ifndef OPENMM_SERIALIZATIONPROXY_H_
#define OPENMM_SERIALIZATIONPROXY_H_

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

#include "openmm/internal/windowsExport.h"
#include <map>
#include <string>
#include <typeinfo>

namespace OpenMM {

class SerializationNode;

/**
 * A SerializationProxy is an object that knows how to serialize and deserialize objects of a
 * particular type.  This is an abstract class.  Subclasses implement the logic for serializing
 * particular types of logic.
 *
 * A global registry maintains the list of what SerializationProxy to use for each type of
 * object.  Call registerProxy() to register the proxy for a particular type.  This is typically
 * done at application startup or by a dynamic library's initialization code.
 */

class OPENMM_EXPORT SerializationProxy {
public:
    /**
     * Create a new SerializationProxy.
     *
     * @param typeName     the name of the object type this proxy knows how to serialize.  This
     *                     name is stored in the output stream during serialization, and is used
     *                     to select a proxy during deserialization.  This typically is the class
     *                     name, although that is not a requirement.
     */
    SerializationProxy(const std::string& typeName);
    virtual ~SerializationProxy() {
    }
    /**
     * Get the name of the object type this proxy manipulates, as passed to the constructor.
     */
    const std::string& getTypeName() const;
    /**
     * Subclasses implement this method to record information about an object being serialized.
     *
     * @param object      a pointer to the object being serialized
     * @param node        all data to be serialized should be stored into this node, either directly
     *                    as properties or indirectly by adding child nodes to it
     */
    virtual void serialize(const void* object, SerializationNode& node) const = 0;
    /**
     * Reconstruct an object from its serialized data.
     *
     * @param node    a SerializationNode containing the object's description
     * @return a pointer to a new object created from the data.  The caller assumes ownership
     * of the object.
     */
    virtual void* deserialize(const SerializationNode& node) const = 0;
    /**
     * Register a SerializationProxy to be used for objects of a particular type.
     *
     * @param type    the type_info for the object type
     * @param proxy   the proxy to use for objects of the specified type
     */
    static void registerProxy(const std::type_info& type, const SerializationProxy* proxy);
    /**
     * Get the SerializationProxy to use for objects of a particular type, specified by name.
     *
     * @param typeName    the name of the object type to get a proxy for
     */
    static const SerializationProxy& getProxy(const std::string& typeName);
    /**
     * Get the SerializationProxy to use for objects of a particular type, specified by type_info.
     *
     * @param type    the type_info of the object type to get a proxy for
     */
    static const SerializationProxy& getProxy(const std::type_info& type);
private:
    std::string typeName;
    static std::map<const std::string, const SerializationProxy*>& getProxiesByType();
    static std::map<const std::string, const SerializationProxy*>& getProxiesByName();
};

} // namespace OpenMM

#endif /*OPENMM_SERIALIZATIONPROXY_H_*/
