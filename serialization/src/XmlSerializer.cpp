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

#include "openmm/serialization/XmlSerializer.h"
#include "irrXML.h"
#include <cstring>
#include <iostream>
#include <map>

using namespace OpenMM;
using namespace std;
using namespace irr;
using namespace io;

/**
 * Apply XML encoding to a string.  This is adapted from TinyXML (written by Lee Thomason).
 */
static void encodeString(const string& str, string* outString) {
    static map<char, string> entities;
    static bool hasInitialized = false;
    if (!hasInitialized) {
        hasInitialized = true;
	entities['&'] = "&amp;";
	entities['<'] = "&lt;";
	entities['>'] = "&gt;";
	entities['\"'] = "&quot;";
	entities['\''] = "&apos;";
    }

    int i=0;

    while (i<(int)str.length()) {
        unsigned char c = (unsigned char) str[i];

        if (c == '&' 
             && i < ((int)str.length() - 2)
             && str[i+1] == '#'
             && str[i+2] == 'x') {
            // Hexadecimal character reference.
            // Pass through unchanged.
            // &#xA9;	-- copyright symbol, for example.
            //
            // The -1 is a bug fix from Rob Laveaux. It keeps
            // an overflow from happening if there is no ';'.
            // There are actually 2 ways to exit this loop -
            // while fails (error case) and break (semicolon found).
            // However, there is no mechanism (currently) for
            // this function to return an error.
            while (i<(int)str.length()-1) {
                outString->append(str.c_str() + i, 1);
                ++i;
                if (str[i] == ';')
                    break;
            }
        }
        else if (entities.find(c) != entities.end()) {
            outString->append(entities[c]);
            ++i;
        }
        else if (c < 32) {
            // Easy pass at non-alpha/numeric/symbol
            // Below 32 is symbolic.
            char buf[ 32 ];

            sprintf(buf, "&#x%02X;", (unsigned) (c & 0xff));

            //*ME:	warning C4267: convert 'size_t' to 'int'
            //*ME:	Int-Cast to make compiler happy ...
            outString->append(buf, (int)strlen(buf));
            ++i;
        }
        else {
            //char realc = (char) c;
            //outString->append(&realc, 1);
            *outString += (char) c;	// somewhat more efficient function call.
            ++i;
        }
    }
}

void XmlSerializer::serialize(const SerializationNode& node, std::ostream& stream) {
    stream << "<?xml version=\"1.0\" ?>\n";
    encodeNode(node, stream, 0);
}

void XmlSerializer::encodeNode(const SerializationNode& node, std::ostream& stream, int depth) {
    for (int i = 0; i < depth; i++)
        stream << '\t';
    stream << '<' << node.getName();
    for (auto& prop : node.getProperties()) {
        string name, value;
        encodeString(prop.first, &name);
        encodeString(prop.second, &value);
        stream << ' ' << name << "=\"" << value << '\"';
    }
    const vector<SerializationNode>& children = node.getChildren();
    if (children.size() == 0)
        stream << "/>\n";
    else {
        stream << ">\n";
        for (auto& child : children)
            encodeNode(child, stream, depth+1);
        for (int i = 0; i < depth; i++)
            stream << '\t';
        stream << "</" << node.getName() << ">\n";
    }
}

/**
 * Adapter class to let irrXML read a C++ stream.
 */
class XmlSerializer::StreamReader : public IFileReadCallBack {
public:
    StreamReader(std::istream& stream) : stream(stream) {
        stream.seekg(0, ios_base::end);
        size = stream.tellg();
        stream.seekg(0);
    }
    int read(void* buffer, int sizeToRead) {
        stream.read((char*) buffer, sizeToRead);
        return stream.gcount();
    }
    int getSize() {
        return size;
    }
private:
    std::istream& stream;
    int size;
};

/**
 * Process an XML node, storing its content into a SerializationNode.
 */
static void decodeNode(SerializationNode& node, IrrXMLReader& xml) {
    for (int i = 0; i < xml.getAttributeCount(); i++)
        node.setStringProperty(xml.getAttributeName(i), xml.getAttributeValue(i));
    if (xml.isEmptyElement())
        return;
    while (xml.read()) {
        switch (xml.getNodeType()) {
            case EXN_ELEMENT:
            {
                SerializationNode& childNode = node.createChildNode(xml.getNodeName());
                decodeNode(childNode, xml);
                break;
            }
            case EXN_ELEMENT_END:
                return;
        }
    }
}

void* XmlSerializer::deserializeStream(std::istream& stream) {
    SerializationNode root;
    StreamReader reader(stream);
    IrrXMLReader* xml = createIrrXMLReader(&reader);
    
    // Find the root node in the file.
    
    while (xml->read() && xml->getNodeType() != EXN_ELEMENT)
        ;
    decodeNode(root, *xml);
    delete xml;
    
    // Process the SerializationNodes.
    
    const SerializationProxy& proxy = SerializationProxy::getProxy(root.getStringProperty("type"));
    return proxy.deserialize(root);
}
