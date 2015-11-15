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
 * A SerializationNode stores information about an object during serialization or deserialization.
 *
 * When an object is serialized, its SerializationProxy is first called to copy information about the
 * object into a SerializationNode.  That information can then be written to the output stream in the
 * desired format.
 *
 * When an object is deserialized, the input stream is read and the information is stored into a
 * SerializationNode.  The appropriate SerializationProxy is then called to reconstruct the object.
 *
 * SerializationNodes are arranged in a tree.  There will often be a one-to-one correspondence between
 * objects and SerializationNodes, but that need not always be true.  A proxy is free to create whatever
 * child nodes it wants and store information in them using whatever organization is most convenient.
 *
 * Each SerializationNode can store an arbitrary set of "properties", represented as key-value pairs.
 * The key is always a string, while the value may be a string, an int, or a double.  If a value is specified
 * using one data type and then accessed as a different data type, the node will attempt to convert the value
 * in an appropriate way.  For example, it is always reasonable to call getStringProperty() to access a
 * property as a string.  Similarly, you can use setStringProperty() to specify a property and then access it
 * using getIntProperty().  This will produce the expected result if the original value was, in fact, the
 * string representation of an int, but if the original string was non-numeric, the result is undefined.
 */

class OPENMM_EXPORT SerializationNode {
public:
    /**
     * Get the name of this SerializationNode.
     */
    const std::string& getName() const;
    /**
     * Set the name of this SerializationNode.
     *
     * @param name    the new name of the SerializationNode
     */
    void setName(const std::string& name);
    /**
     * Get a reference to this node's child nodes.
     */
    const std::vector<SerializationNode>& getChildren() const;
    /**
     * Get a reference to this node's child nodes.
     */
    std::vector<SerializationNode>& getChildren();
    /**
     * Get a reference to the child node with a particular name.  If there is no child
     * with the specified name, this throws an exception.
     *
     * @param the name of the child node to get
     */
    const SerializationNode& getChildNode(const std::string& name) const;
    /**
     * Get a reference to the child node with a particular name.  If there is no child
     * with the specified name, this throws an exception.
     *
     * @param the name of the child node to get
     */
    SerializationNode& getChildNode(const std::string& name);
    /**
     * Get a map containing all of this node's properties.
     */
    const std::map<std::string, std::string>& getProperties() const;
    /**
     * Determine whether this node has a property with a particular node.
     *
     * @param name  the name of the property to check for
     */
    bool hasProperty(const std::string& name) const;
    /**
     * Get the property with a particular name, specified as a string.  If there is no property with
     * the specified name, an exception is thrown.
     *
     * @param name   the name of the property to get
     */
    const std::string& getStringProperty(const std::string& name) const;
    /**
     * Get the property with a particular name, specified as a string.  If there is no property with
     * the specified name, a default value is returned instead.
     *
     * @param name          the name of the property to get
     * @param defaultValue  the value to return if the specified property does not exist
     */
    const std::string& getStringProperty(const std::string& name, const std::string& defaultValue) const;
    /**
     * Set the value of a property, specified as a string.
     *
     * @param name   the name of the property to set
     * @param value  the value to set for the property
     */
    SerializationNode& setStringProperty(const std::string& name, const std::string& value);
    /**
     * Get the property with a particular name, specified as an int.  If there is no property with
     * the specified name, an exception is thrown.
     *
     * @param name   the name of the property to get
     */
    int getIntProperty(const std::string& name) const;
    /**
     * Get the property with a particular name, specified as an int.  If there is no property with
     * the specified name, a default value is returned instead.
     *
     * @param name          the name of the property to get
     * @param defaultValue  the value to return if the specified property does not exist
     */
    int getIntProperty(const std::string& name, int defaultValue) const;
    /**
     * Set the value of a property, specified as an int.
     *
     * @param name   the name of the property to set
     * @param value  the value to set for the property
     */
    SerializationNode& setIntProperty(const std::string& name, int value);
    /**
     * Get the property with a particular name, specified as an bool.  If there is no property with
     * the specified name, an exception is thrown.
     *
     * @param name   the name of the property to get
     */
    bool getBoolProperty(const std::string& name) const;
    /**
     * Get the property with a particular name, specified as a bool.  If there is no property with
     * the specified name, a default value is returned instead.
     *
     * @param name          the name of the property to get
     * @param defaultValue  the value to return if the specified property does not exist
     */
    bool getBoolProperty(const std::string& name, bool defaultValue) const;
    /**
     * Set the value of a property, specified as a bool.
     *
     * @param name   the name of the property to set
     * @param value  the value to set for the property
     */
    SerializationNode& setBoolProperty(const std::string& name, bool value);
    /**
     * Get the property with a particular name, specified as a double.  If there is no property with
     * the specified name, an exception is thrown.
     *
     * @param name   the name of the property to get
     */
    double getDoubleProperty(const std::string& name) const;
    /**
     * Get the property with a particular name, specified as a double.  If there is no property with
     * the specified name, a default value is returned instead.
     *
     * @param name          the name of the property to get
     * @param defaultValue  the value to return if the specified property does not exist
     */
    double getDoubleProperty(const std::string& name, double defaultValue) const;
    /**
     * Set the value of a property, specified as a double.
     *
     * @param name   the name of the property to set
     * @param value  the value to set for the property
     */
    SerializationNode& setDoubleProperty(const std::string& name, double value);
    /**
     * Create a new child node
     *
     * @param name    the name of the new node to create
     * @return a reference to the newly created node
     */
    SerializationNode& createChildNode(const std::string& name);
    /**
     * Create a new child node by serializing an object.  A SerializationProxy is automatically
     * selected based on the object's type, then invoked to populate the newly created node.
     *
     * Note that, while this method is templatized based on the type of object being serialized,
     * the typeid() operator is used to select the proxy.  This means the template argument may
     * be a base class, and the correct proxies will still be selected for objects of different
     * subclasses.
     *
     * @param name    the name of the new node to create
     * @param object  a pointer to the object to serialize
     * @return a reference to the newly created node
     */
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
    /**
     * Reconstruct an object based on the information stored in this node.  A SerializationProxy is
     * automatically selected based on the information stored in the node, then it is invoked to
     * create the object.
     *
     * The template parameter may be either the actual type of the object, or any base class to which
     * it may be cast.
     *
     * @return a pointer to the newly created object.  The caller assumes ownership of the object.
     */
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
