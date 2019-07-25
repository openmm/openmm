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
#include <typeinfo>

using namespace OpenMM;
using namespace std;

map<const string, const SerializationProxy*>& SerializationProxy::getProxiesByType() {
    static map<const string, const SerializationProxy*> proxies;
    return proxies;
}

map<const string, const SerializationProxy*>& SerializationProxy::getProxiesByName() {
    static map<const string, const SerializationProxy*> proxies;
    return proxies;
}

SerializationProxy::SerializationProxy(const string& typeName) : typeName(typeName) {
}

const string& SerializationProxy::getTypeName() const {
    return typeName;
}

void SerializationProxy::registerProxy(const type_info& type, const SerializationProxy* proxy) {
    getProxiesByType()[type.name()] = proxy;
    getProxiesByName()[proxy->getTypeName()] = proxy;
}

const SerializationProxy& SerializationProxy::getProxy(const string& typeName) {
    map<const string, const SerializationProxy*>::const_iterator iter = getProxiesByName().find(typeName);
    if (iter == getProxiesByName().end())
        throw OpenMMException("There is no serialization proxy registered for type '"+string(typeName)+"'");
    return *iter->second;
}

const SerializationProxy& SerializationProxy::getProxy(const type_info& type) {
    map<const string, const SerializationProxy*>::const_iterator iter = getProxiesByType().find(type.name());
    if (iter == getProxiesByType().end())
        throw OpenMMException("There is no serialization proxy registered for type "+string(type.name()));
    return *iter->second;
}
