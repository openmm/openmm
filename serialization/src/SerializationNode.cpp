/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit.                   *
 * See https://openmm.org/development.                                        *
 *                                                                            *
 * Portions copyright (c) 2010-2026 Stanford University and the Authors.      *
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
#include "openmm/OpenMMException.h"
#include <cstring>
#include <sstream>

using namespace OpenMM;
using namespace std;

extern "C" char* g_fmt(char*, double);
extern "C" double strtod2(const char* s00, char** se);

const string& SerializationNode::getName() const {
    return name;
}

void SerializationNode::setName(const string& name) {
    this->name = name;
}

const vector<SerializationNode>& SerializationNode::getChildren() const {
    return children;
}

vector<SerializationNode>& SerializationNode::getChildren() {
    return children;
}

const SerializationNode& SerializationNode::getChildNode(const std::string& name) const {
    for (auto& child : children)
        if (child.name == name)
            return child;
        throw OpenMMException("Unknown child '"+name+"' for node '"+getName()+"'");
}

SerializationNode& SerializationNode::getChildNode(const std::string& name) {
    for (auto& child : children)
        if (child.name == name)
            return child;
        throw OpenMMException("Unknown child '"+name+"' for node '"+getName()+"'");
}

const map<string, string> SerializationNode::getProperties() const {
    map<string, string> result;
    int start = 0;
    while (start < properties.size()) {
        string name(&properties[start]);
        start += name.size()+1;
        string value(&properties[start]);
        start += value.size()+1;
        result[name] = value;
    }
    return result;
}

const char* SerializationNode::findPropertyValue(const std::string& name, bool required) const {
    int start = 0;
    while (start < properties.size()) {
        if (strcmp(&properties[start], name.data()) == 0)
            return &properties[start+name.size()+1];

        // This isn't the property we're looking for.  Skip over the name and value.

        start += strlen(&properties[start])+1;
        start += strlen(&properties[start])+1;
    }
    if (required)
        throw OpenMMException("Unknown property '"+name+"' in node '"+getName()+"'");

    // We didn't find it, but that's ok.  Report that it wasn't found.

    return nullptr;
}

bool SerializationNode::hasProperty(const string& name) const {
    return findPropertyValue(name, false) != nullptr;
}

const string SerializationNode::getStringProperty(const string& name) const {
    return string(findPropertyValue(name, true));
}

const string SerializationNode::getStringProperty(const string& name, const string& defaultValue) const {
    const char* value = findPropertyValue(name, false);
    if (value == nullptr)
        return defaultValue;
    return string(value);
}

SerializationNode& SerializationNode::setStringProperty(const string& name, const string& value) {
    int s = properties.size()+name.size()+value.size()+2;
    properties.reserve(properties.size()+name.size()+value.size()+2);
    properties.append(name);
    properties.push_back(0);
    properties.append(value);
    properties.push_back(0);
    return *this;
}

int SerializationNode::getIntProperty(const string& name) const {
    int value;
    stringstream(findPropertyValue(name, true)) >> value;
    return value;
}

int SerializationNode::getIntProperty(const string& name, int defaultValue) const {
    const char* value = findPropertyValue(name, false);
    if (value == nullptr)
        return defaultValue;
    int result;
    stringstream(value) >> result;
    return result;
}

SerializationNode& SerializationNode::setIntProperty(const string& name, int value) {
    stringstream s;
    s << value;
    return setStringProperty(name, s.str());
}

long long SerializationNode::getLongProperty(const string& name) const {
    long long value;
    stringstream(findPropertyValue(name, true)) >> value;
    return value;
}

long long SerializationNode::getLongProperty(const string& name, long long defaultValue) const {
    const char* value = findPropertyValue(name, false);
    if (value == nullptr)
        return defaultValue;
    long long result;
    stringstream(value) >> result;
    return result;
}

SerializationNode& SerializationNode::setLongProperty(const string& name, long long value) {
    stringstream s;
    s << value;
    return setStringProperty(name, s.str());
}

bool SerializationNode::getBoolProperty(const string& name) const {
    bool value;
    stringstream(findPropertyValue(name, true)) >> value;
    return value;
}

bool SerializationNode::getBoolProperty(const string& name, bool defaultValue) const {
    const char* value = findPropertyValue(name, false);
    if (value == nullptr)
        return defaultValue;
    bool result;
    stringstream(value) >> result;
    return result;
}

SerializationNode& SerializationNode::setBoolProperty(const string& name, bool value) {
    stringstream s;
    s << value;
    return setStringProperty(name, s.str());
}

double SerializationNode::getDoubleProperty(const string& name) const {
    return strtod2(findPropertyValue(name, true), NULL);
}

double SerializationNode::getDoubleProperty(const string& name, double defaultValue) const {
    const char* value = findPropertyValue(name, false);
    if (value == nullptr)
        return defaultValue;
    return strtod2(value, NULL);
}

SerializationNode& SerializationNode::setDoubleProperty(const string& name, double value) {
    char buffer[32];
    g_fmt(buffer, value);
    return setStringProperty(name, string(buffer));
}

SerializationNode& SerializationNode::createChildNode(const std::string& name) {
    children.push_back(SerializationNode());
    children.back().setName(name);
    return children.back();
}
