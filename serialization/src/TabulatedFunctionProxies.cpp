/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "openmm/serialization/TabulatedFunctionProxies.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/TabulatedFunction.h"
#include <sstream>

using namespace OpenMM;
using namespace std;

Continuous1DFunctionProxy::Continuous1DFunctionProxy() : SerializationProxy("Continuous1DFunction") {
}

void Continuous1DFunctionProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const Continuous1DFunction& function = *reinterpret_cast<const Continuous1DFunction*>(object);
    double min, max;
    vector<double> values;
    function.getFunctionParameters(values, min, max);
    node.setDoubleProperty("min", min);
    node.setDoubleProperty("max", max);
    SerializationNode& valuesNode = node.createChildNode("Values");
    for (auto v : values)
        valuesNode.createChildNode("Value").setDoubleProperty("v", v);
}

void* Continuous1DFunctionProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& valuesNode = node.getChildNode("Values");
    vector<double> values;
    for (auto& child : valuesNode.getChildren())
        values.push_back(child.getDoubleProperty("v"));
    return new Continuous1DFunction(values, node.getDoubleProperty("min"), node.getDoubleProperty("max"));
}

Continuous2DFunctionProxy::Continuous2DFunctionProxy() : SerializationProxy("Continuous2DFunction") {
}

void Continuous2DFunctionProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const Continuous2DFunction& function = *reinterpret_cast<const Continuous2DFunction*>(object);
    int xsize, ysize;
    double xmin, xmax, ymin, ymax;
    vector<double> values;
    function.getFunctionParameters(xsize, ysize, values, xmin, xmax, ymin, ymax);
    node.setDoubleProperty("xsize", xsize);
    node.setDoubleProperty("ysize", ysize);
    node.setDoubleProperty("xmin", xmin);
    node.setDoubleProperty("xmax", xmax);
    node.setDoubleProperty("ymin", ymin);
    node.setDoubleProperty("ymax", ymax);
    SerializationNode& valuesNode = node.createChildNode("Values");
    for (auto v : values)
        valuesNode.createChildNode("Value").setDoubleProperty("v", v);
}

void* Continuous2DFunctionProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& valuesNode = node.getChildNode("Values");
    vector<double> values;
    for (auto& child : valuesNode.getChildren())
        values.push_back(child.getDoubleProperty("v"));
    return new Continuous2DFunction(node.getIntProperty("xsize"), node.getIntProperty("ysize"), values,
            node.getDoubleProperty("xmin"), node.getDoubleProperty("xmax"), node.getDoubleProperty("ymin"), node.getDoubleProperty("ymax"));
}

Continuous3DFunctionProxy::Continuous3DFunctionProxy() : SerializationProxy("Continuous3DFunction") {
}

void Continuous3DFunctionProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const Continuous3DFunction& function = *reinterpret_cast<const Continuous3DFunction*>(object);
    int xsize, ysize, zsize;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    vector<double> values;
    function.getFunctionParameters(xsize, ysize, zsize, values, xmin, xmax, ymin, ymax, zmin, zmax);
    node.setDoubleProperty("xsize", xsize);
    node.setDoubleProperty("ysize", ysize);
    node.setDoubleProperty("zsize", zsize);
    node.setDoubleProperty("xmin", xmin);
    node.setDoubleProperty("xmax", xmax);
    node.setDoubleProperty("ymin", ymin);
    node.setDoubleProperty("ymax", ymax);
    node.setDoubleProperty("zmin", zmin);
    node.setDoubleProperty("zmax", zmax);
    SerializationNode& valuesNode = node.createChildNode("Values");
    for (auto v : values)
        valuesNode.createChildNode("Value").setDoubleProperty("v", v);
}

void* Continuous3DFunctionProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& valuesNode = node.getChildNode("Values");
    vector<double> values;
    for (auto& child : valuesNode.getChildren())
        values.push_back(child.getDoubleProperty("v"));
    return new Continuous3DFunction(node.getIntProperty("xsize"), node.getIntProperty("ysize"), node.getIntProperty("zsize"), values,
            node.getDoubleProperty("xmin"), node.getDoubleProperty("xmax"), node.getDoubleProperty("ymin"), node.getDoubleProperty("ymax"),
            node.getDoubleProperty("zmin"), node.getDoubleProperty("zmax"));
}

Discrete1DFunctionProxy::Discrete1DFunctionProxy() : SerializationProxy("Discrete1DFunction") {
}

void Discrete1DFunctionProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const Discrete1DFunction& function = *reinterpret_cast<const Discrete1DFunction*>(object);
    vector<double> values;
    function.getFunctionParameters(values);
    SerializationNode& valuesNode = node.createChildNode("Values");
    for (auto v : values)
        valuesNode.createChildNode("Value").setDoubleProperty("v", v);
}

void* Discrete1DFunctionProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& valuesNode = node.getChildNode("Values");
    vector<double> values;
    for (auto& child : valuesNode.getChildren())
        values.push_back(child.getDoubleProperty("v"));
    return new Discrete1DFunction(values);
}

Discrete2DFunctionProxy::Discrete2DFunctionProxy() : SerializationProxy("Discrete2DFunction") {
}

void Discrete2DFunctionProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const Discrete2DFunction& function = *reinterpret_cast<const Discrete2DFunction*>(object);
    int xsize, ysize;
    vector<double> values;
    function.getFunctionParameters(xsize, ysize, values);
    node.setDoubleProperty("xsize", xsize);
    node.setDoubleProperty("ysize", ysize);
    SerializationNode& valuesNode = node.createChildNode("Values");
    for (auto v : values)
        valuesNode.createChildNode("Value").setDoubleProperty("v", v);
}

void* Discrete2DFunctionProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& valuesNode = node.getChildNode("Values");
    vector<double> values;
    for (auto& child : valuesNode.getChildren())
        values.push_back(child.getDoubleProperty("v"));
    return new Discrete2DFunction(node.getIntProperty("xsize"), node.getIntProperty("ysize"), values);
}

Discrete3DFunctionProxy::Discrete3DFunctionProxy() : SerializationProxy("Discrete3DFunction") {
}

void Discrete3DFunctionProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const Discrete3DFunction& function = *reinterpret_cast<const Discrete3DFunction*>(object);
    int xsize, ysize, zsize;
    vector<double> values;
    function.getFunctionParameters(xsize, ysize, zsize, values);
    node.setDoubleProperty("xsize", xsize);
    node.setDoubleProperty("ysize", ysize);
    node.setDoubleProperty("zsize", zsize);
    SerializationNode& valuesNode = node.createChildNode("Values");
    for (auto v : values)
        valuesNode.createChildNode("Value").setDoubleProperty("v", v);
}

void* Discrete3DFunctionProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& valuesNode = node.getChildNode("Values");
    vector<double> values;
    for (auto& child : valuesNode.getChildren())
        values.push_back(child.getDoubleProperty("v"));
    return new Discrete3DFunction(node.getIntProperty("xsize"), node.getIntProperty("ysize"), node.getIntProperty("zsize"), values);
}
