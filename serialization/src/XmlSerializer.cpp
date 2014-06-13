/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2014 Stanford University and the Authors.      *
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

#include "openmm/serialization/XmlSerializer.h"
#include "tinyxml.h"

using namespace OpenMM;
using namespace std;

void XmlSerializer::serialize(const SerializationNode& node, std::ostream& stream) {
    stream << "<?xml version=\"1.0\" ?>\n";
    encodeNode(node, stream, 0);
}

void XmlSerializer::encodeNode(const SerializationNode& node, std::ostream& stream, int depth) {
    for (int i = 0; i < depth; i++)
        stream << '\t';
    stream << '<' << node.getName();
    const map<string, string>& properties = node.getProperties();
    for (map<string, string>::const_iterator iter = properties.begin(); iter != properties.end(); ++iter) {
        string name, value;
        TiXmlBase::EncodeString(iter->first, &name);
        TiXmlBase::EncodeString(iter->second, &value);
        stream << ' ' << name << "=\"" << value << '\"';
    }
    const vector<SerializationNode>& children = node.getChildren();
    if (children.size() == 0)
        stream << "/>\n";
    else {
        stream << ">\n";
        for (int i = 0; i < (int) children.size(); i++)
            encodeNode(children[i], stream, depth+1);
        for (int i = 0; i < depth; i++)
            stream << '\t';
        stream << "</" << node.getName() << ">\n";
    }
}

void* XmlSerializer::deserializeStream(std::istream& stream) {
    TiXmlDocument doc;
    stream >> doc;
    SerializationNode root;
    decodeNode(root, *doc.FirstChildElement());
    const SerializationProxy& proxy = SerializationProxy::getProxy(root.getStringProperty("type"));
    return proxy.deserialize(root);
}

void XmlSerializer::decodeNode(SerializationNode& node, const TiXmlElement& element) {
    for (const TiXmlAttribute* attribute = element.FirstAttribute(); attribute != NULL; attribute = attribute->Next())
        node.setStringProperty(attribute->NameTStr(), attribute->ValueStr());
    for (const TiXmlElement* child = element.FirstChildElement(); child != NULL; child = child->NextSiblingElement()) {
        SerializationNode& childNode = node.createChildNode(child->ValueTStr());
        decodeNode(childNode, *child);
    }
}