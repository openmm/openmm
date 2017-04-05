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
#include "openmm/OpenMMException.h"
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

const map<string, string>& SerializationNode::getProperties() const {
    return properties;
}

bool SerializationNode::hasProperty(const string& name) const {
    return (properties.find(name) != properties.end());
}

const string& SerializationNode::getStringProperty(const string& name) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        throw OpenMMException("Unknown property '"+name+"' in node '"+getName()+"'");
    return iter->second;
}

const string& SerializationNode::getStringProperty(const string& name, const string& defaultValue) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        return defaultValue;
    return iter->second;
}

SerializationNode& SerializationNode::setStringProperty(const string& name, const string& value) {
    properties[name] = value;
    return *this;
}

int SerializationNode::getIntProperty(const string& name) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        throw OpenMMException("Unknown property '"+name+"' in node '"+getName()+"'");
    int value;
    stringstream(iter->second) >> value;
    return value;
}

int SerializationNode::getIntProperty(const string& name, int defaultValue) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        return defaultValue;
    int value;
    stringstream(iter->second) >> value;
    return value;
}

SerializationNode& SerializationNode::setIntProperty(const string& name, int value) {
    stringstream s;
    s << value;
    properties[name] = s.str();
    return *this;
}


bool SerializationNode::getBoolProperty(const string& name) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        throw OpenMMException("Unknown property '"+name+"' in node '"+getName()+"'");
    bool value;
    stringstream(iter->second) >> value;
    return value;
}

bool SerializationNode::getBoolProperty(const string& name, bool defaultValue) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        return defaultValue;
    bool value;
    stringstream(iter->second) >> value;
    return value;
}

SerializationNode& SerializationNode::setBoolProperty(const string& name, bool value) {
    stringstream s;
    s << value;
    properties[name] = s.str();
    return *this;
}

double SerializationNode::getDoubleProperty(const string& name) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        throw OpenMMException("Unknown property '"+name+"' in node '"+getName()+"'");
    return strtod2(iter->second.c_str(), NULL);
}

double SerializationNode::getDoubleProperty(const string& name, double defaultValue) const {
    map<string, string>::const_iterator iter = properties.find(name);
    if (iter == properties.end())
        return defaultValue;
    return strtod2(iter->second.c_str(), NULL);
}

SerializationNode& SerializationNode::setDoubleProperty(const string& name, double value) {
    char buffer[32];
    g_fmt(buffer, value);
    properties[name] = string(buffer);
    return *this;
}

SerializationNode& SerializationNode::createChildNode(const std::string& name) {
    children.push_back(SerializationNode());
    children.back().setName(name);
    return children.back();
}
