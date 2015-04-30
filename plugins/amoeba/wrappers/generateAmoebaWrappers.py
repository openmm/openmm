from __future__ import print_function
import sys, os
import time
import getopt
import re
import xml.etree.ElementTree as etree

def trimToSingleSpace(text):
    if text is None or len(text) == 0:
        return ""
    t = text.strip()
    if len(t) == 0:
        return t
    if text[0].isspace():
        t = " %s" % t
    if text[-1].isspace():
        t = "%s " % t
    return t

def getNodeText(node):
    if node.text is not None:
        s = node.text
    else:
        s = ""
    for n in node:
        if n.tag == "para":
            s = "%s%s\n\n" % (s, getNodeText(n))
        elif n.tag == "ref":
            s = "%s%s" % (s, getNodeText(n))
        if n.tail is not None:
            s = "%s%s" % (s, n.tail)
    return s

def getText(subNodePath, node):
    s = ""
    for n in node.findall(subNodePath):
        s = "%s%s" % (s, trimToSingleSpace(getNodeText(n)))
        if n.tag == "para":
            s = "%s\n\n" % s
    return s.strip()

def convertOpenMMPrefix(name):
    return name.replace('OpenMM::', 'OpenMM_')

OPENMM_RE_PATTERN=re.compile("(.*)OpenMM:[a-zA-Z:]*:(.*)")
def stripOpenMMPrefix(name, rePattern=OPENMM_RE_PATTERN):
    try:
        m=rePattern.search(name)
        rValue = "%s%s" % m.group(1,2)
        rValue.strip()
        return rValue
    except:
        return name

def findNodes(parent, path, **args):
    nodes = []
    for node in parent.findall(path):
        match = True
        for arg in args:
            if arg not in node.attrib or node.attrib[arg] != args[arg]:
                match = False
        if match:
            nodes.append(node)
    return nodes

class WrapperGenerator:
    """This is the parent class of generators for various API wrapper files.  It defines functions common to all of them."""
    
    def __init__(self, inputDirname, output):
        self.skipClasses = []
        self.skipMethods = []
        self.hideClasses = ['Kernel', 'KernelImpl', 'KernelFactory', 'ContextImpl', 'SerializationNode', 'SerializationProxy']
        self.nodeByID={}

        # Read all the XML files and merge them into a single document.
        self.doc = etree.ElementTree(etree.Element('root'))
        for file in os.listdir(inputDirname):
            if file.lower().endswith('xml'):
                root = etree.parse(os.path.join(inputDirname, file)).getroot()
                for node in root:
                    self.doc.getroot().append(node)

        self.out = output

        self.typesByShortName = {}
        self._orderedClassNodes = self.buildOrderedClassNodes()

    def getNodeByID(self, id):
        if id not in self.nodeByID:
            for node in findNodes(self.doc.getroot(), "compounddef", id=id):
                self.nodeByID[id] = node
        return self.nodeByID[id]

    def buildOrderedClassNodes(self):
        orderedClassNodes=[]
        for node in findNodes(self.doc.getroot(), "compounddef", kind="class", prot="public"):
            self.findBaseNodes(node, orderedClassNodes)
        return orderedClassNodes

    def findBaseNodes(self, node, excludedClassNodes=[]):
        if node in excludedClassNodes:
            return
        if node.attrib['prot'] == 'private':
            return
        nodeName = getText("compoundname", node)
        if nodeName in self.skipClasses:
            return
        for baseNodePnt in findNodes(node, "basecompoundref", prot="public"):
            if "refid" in baseNodePnt.attrib:
                baseNodeID = baseNodePnt.attrib["refid"]
                baseNode = self.getNodeByID(baseNodeID)
                self.findBaseNodes(baseNode, excludedClassNodes)
        excludedClassNodes.append(node)

    def getClassMethods(self, classNode):
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        methodList = []
        for section in findNodes(classNode, "sectiondef", kind="public-static-func")+findNodes(classNode, "sectiondef", kind="public-func"):
            for memberNode in findNodes(section, "memberdef", kind="function", prot="public"):
                methodDefinition = getText("definition", memberNode)
                shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
                methodName = shortMethodDefinition.split()[-1]
                if className+'::'+methodName in self.skipMethods:
                    continue
                methodList.append(memberNode)
        return methodList
    
    def shouldHideType(self, typeName):
        if typeName.startswith('const '):
            typeName = typeName[6:].strip()
        if typeName.endswith('&') or typeName.endswith('*'):
            typeName = typeName[:-1].strip()
        return typeName in self.hideClasses
    
    def shouldHideMethod(self, methodNode):
        paramList = findNodes(methodNode, 'param')
        returnType = self.getType(getText("type", methodNode))
        if self.shouldHideType(returnType):
            return True
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if self.shouldHideType(type):
                return True
        return False

class CHeaderGenerator(WrapperGenerator):
    """This class generates the header file for the C API wrappers."""
    
    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
        self.typeTranslations = {'bool': 'OpenMM_Boolean',
                                 'Vec3': 'OpenMM_Vec3',
                                 'Context': 'OpenMM_Context',
                                 'std::string': 'char*',
                                 'const std::string &': 'const char*',
                                 'std::vector< std::string >': 'OpenMM_StringArray',
                                 'std::vector< Vec3 >': 'OpenMM_Vec3Array',
                                 'std::vector< std::pair< int, int > >': 'OpenMM_BondArray',
                                 'std::map< std::string, double >': 'OpenMM_ParameterArray',
                                 'std::map< std::string, std::string >': 'OpenMM_PropertyArray',
                                 'std::vector< double >': 'OpenMM_DoubleArray',
                                 'std::vector< int >': 'OpenMM_IntArray',
                                 'std::set< int >': 'OpenMM_IntSet',
                                 'std::vector< std::vector< int > >': 'OpenMM_2D_IntArray',
                                 'std::vector< std::vector< std::vector< double > > >': 'OpenMM_3D_DoubleArray'}
    
    def writeGlobalConstants(self):
        self.out.write("/* Global Constants */\n\n")
        node = next((x for x in findNodes(self.doc.getroot(), "compounddef", kind="namespace") if x.findtext("compoundname") == "OpenMM"))
        for section in findNodes(node, "sectiondef", kind="var"):
            for memberNode in findNodes(section, "memberdef", kind="variable", mutable="no", prot="public", static="yes"):
                vDef = convertOpenMMPrefix(getText("definition", memberNode))
                iDef = getText("initializer", memberNode)
                if iDef.startswith("="):
                    iDef = iDef[1:]
                self.out.write("static %s = %s;\n" % (vDef, iDef))

    def writeTypeDeclarations(self):
        self.out.write("\n/* Type Declarations */\n\n")
        for classNode in self._orderedClassNodes:
            className = getText("compoundname", classNode)
            shortName = stripOpenMMPrefix(className)
            typeName = convertOpenMMPrefix(className)
            self.out.write("typedef struct %s_struct %s;\n" % (typeName, typeName))
            self.typesByShortName[shortName] = typeName

    def writeClasses(self):
        for classNode in self._orderedClassNodes:
            className = stripOpenMMPrefix(getText("compoundname", classNode))
            self.out.write("\n/* %s */\n" % className)
            self.writeEnumerations(classNode)
            self.writeMethods(classNode)
        self.out.write("\n")

    def writeEnumerations(self, classNode):
        enumNodes = []
        for section in findNodes(classNode, "sectiondef", kind="public-type"):
            for node in findNodes(section, "memberdef", kind="enum", prot="public"):
                enumNodes.append(node)
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        for enumNode in enumNodes:
            enumName = getText("name", enumNode)
            enumTypeName = "%s_%s" % (typeName, enumName)
            self.out.write("typedef enum {\n  ")
            argSep=""
            for valueNode in findNodes(enumNode, "enumvalue", prot="public"):
                vName = convertOpenMMPrefix(getText("name", valueNode))
                vInit = getText("initializer", valueNode)
                if vInit.startswith("="):
                    vInit = vInit[1:].strip()
                self.out.write("%s%s_%s = %s" % (argSep, typeName, vName, vInit))
                argSep=", "
            self.out.write("\n} %s;\n" % enumTypeName)
            self.typesByShortName[enumName] = enumTypeName
        if len(enumNodes)>0: self.out.write("\n")

    def writeMethods(self, classNode):
        methodList = self.getClassMethods(classNode)
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        destructorName = '~'+shortClassName

        if not ('abstract' in classNode.attrib and classNode.attrib['abstract'] == 'yes'):
            # Write constructors
            numConstructors = 0
            for methodNode in methodList:
                methodDefinition = getText("definition", methodNode)
                shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
                methodName = shortMethodDefinition.split()[-1]
                if methodName == shortClassName:
                    if self.shouldHideMethod(methodNode):
                        continue
                    numConstructors += 1
                    if numConstructors == 1:
                        suffix = ""
                    else:
                        suffix = "_%d" % numConstructors
                    self.out.write("extern OPENMM_EXPORT_AMOEBA %s* %s_create%s(" % (typeName, typeName, suffix))
                    self.writeArguments(methodNode, False)
                    self.out.write(");\n")
    
        # Write destructor
        self.out.write("extern OPENMM_EXPORT_AMOEBA void %s_destroy(%s* target);\n" % (typeName, typeName))

        # Record method names for future reference.
        methodNames = {}
        for methodNode in methodList:
            methodDefinition = getText("definition", methodNode)
            shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
            methodNames[methodNode] = shortMethodDefinition.split()[-1]
        
        # Write other methods
        for methodNode in methodList:
            methodName = methodNames[methodNode]
            if methodName in (shortClassName, destructorName):
                continue
            if self.shouldHideMethod(methodNode):
                continue
            isConstMethod = (methodNode.attrib['const'] == 'yes')
            if isConstMethod and any(methodNames[m] == methodName and m.attrib['const'] == 'no' for m in methodList):
                # There are two identical methods that differ only in whether they are const.  Skip the const one.
                continue
            returnType = self.getType(getText("type", methodNode))
            self.out.write("extern OPENMM_EXPORT_AMOEBA %s %s_%s(" % (returnType, typeName, methodName))
            isInstanceMethod = (methodNode.attrib['static'] != 'yes')
            if isInstanceMethod:
                if isConstMethod:
                    self.out.write('const ')
                self.out.write("%s* target" % typeName)
            self.writeArguments(methodNode, isInstanceMethod)
            self.out.write(");\n")
    
    def writeArguments(self, methodNode, initialSeparator):
        paramList = findNodes(methodNode, 'param')
        if initialSeparator:
            separator = ", "
        else:
            separator = ""
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'void':
                continue
            type = self.getType(type)
            name = getText('declname', node)
            self.out.write("%s%s %s" % (separator, type, name))
            separator = ", "
    
    def getType(self, type):
        if type in self.typeTranslations:
            return self.typeTranslations[type]
        if type in self.typesByShortName:
            return self.typesByShortName[type]
        if type.startswith('const '):
            return 'const '+self.getType(type[6:].strip())
        if type.endswith('&') or type.endswith('*'):
            return self.getType(type[:-1].strip())+'*'
        return type

    def writeOutput(self):
        print("""
#ifndef AMOEBA_OPENMM_CWRAPPER_H_
#define AMOEBA_OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT_AMOEBA
#define OPENMM_EXPORT_AMOEBA
#endif""", file=self.out)
        self.writeGlobalConstants()
        self.writeTypeDeclarations()
        print("""
typedef struct OpenMM_2D_IntArray_struct OpenMM_2D_IntArray;
typedef struct OpenMM_3D_DoubleArray_struct OpenMM_3D_DoubleArray;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_3D_DoubleArray */
OPENMM_EXPORT_AMOEBA OpenMM_3D_DoubleArray* OpenMM_3D_DoubleArray_create(int size1, int size2, int size3);
OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray* array, int index1, int index2, OpenMM_DoubleArray* values);
OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray* array);""", file=self.out)

        self.writeClasses()

        print("""
#if defined(__cplusplus)
}
#endif

#endif /*AMOEBA_OPENMM_CWRAPPER_H_*/""", file=self.out)


class CSourceGenerator(WrapperGenerator):
    """This class generates the source file for the C API wrappers."""

    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
        self.typeTranslations = {'bool': 'OpenMM_Boolean',
                                 'Vec3': 'OpenMM_Vec3',
                                 'Context': 'OpenMM_Context',
                                 'std::string': 'char*',
                                 'const std::string &': 'const char*',
                                 'std::vector< std::string >': 'OpenMM_StringArray',
                                 'std::vector< Vec3 >': 'OpenMM_Vec3Array',
                                 'std::vector< std::pair< int, int > >': 'OpenMM_BondArray',
                                 'std::map< std::string, double >': 'OpenMM_ParameterArray',
                                 'std::map< std::string, std::string >': 'OpenMM_PropertyArray',
                                 'std::vector< double >': 'OpenMM_DoubleArray',
                                 'std::vector< int >': 'OpenMM_IntArray',
                                 'std::set< int >': 'OpenMM_IntSet',
                                 'std::vector< std::vector< int > >': 'OpenMM_2D_IntArray',
                                 'std::vector< std::vector< std::vector< double > > >': 'OpenMM_3D_DoubleArray'}
        self.inverseTranslations = dict((self.typeTranslations[key], key) for key in self.typeTranslations)
        self.classesByShortName = {}
        self.enumerationTypes = {}
        self.findTypes()
    
    def findTypes(self):
        for classNode in self._orderedClassNodes:
            className = getText("compoundname", classNode)
            shortName = stripOpenMMPrefix(className)
            typeName = convertOpenMMPrefix(className)
            self.typesByShortName[shortName] = typeName
            self.classesByShortName[shortName] = className

    def findEnumerations(self, classNode):
        enumNodes = []
        for section in findNodes(classNode, "sectiondef", kind="public-type"):
            for node in findNodes(section, "memberdef", kind="enum", prot="public"):
                enumNodes.append(node)
        className = getText("compoundname", classNode)
        typeName = convertOpenMMPrefix(className)
        for enumNode in enumNodes:
            enumName = getText("name", enumNode)
            enumTypeName = "%s_%s" % (typeName, enumName)
            enumClassName = "%s::%s" % (className, enumName)
            self.typesByShortName[enumName] = enumTypeName
            self.classesByShortName[enumName] = enumClassName
            self.enumerationTypes[enumClassName] = enumTypeName

    def writeClasses(self):
        for classNode in self._orderedClassNodes:
            className = stripOpenMMPrefix(getText("compoundname", classNode))
            self.out.write("\n/* OpenMM::%s */\n" % className)
            self.findEnumerations(classNode)
            self.writeMethods(classNode)
        self.out.write("\n")

    def writeMethods(self, classNode):
        methodList = self.getClassMethods(classNode)
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        destructorName = '~'+shortClassName

        if not ('abstract' in classNode.attrib and classNode.attrib['abstract'] == 'yes'):
            # Write constructors
            numConstructors = 0
            for methodNode in methodList:
                methodDefinition = getText("definition", methodNode)
                shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
                methodName = shortMethodDefinition.split()[-1]
                if methodName == shortClassName:
                    if self.shouldHideMethod(methodNode):
                        continue
                    numConstructors += 1
                    if numConstructors == 1:
                        suffix = ""
                    else:
                        suffix = "_%d" % numConstructors
                    self.out.write("OPENMM_EXPORT_AMOEBA %s* %s_create%s(" % (typeName, typeName, suffix))
                    self.writeArguments(methodNode, False)
                    self.out.write(") {\n")
                    self.out.write("    return reinterpret_cast<%s*>(new %s(" % (typeName, className))
                    self.writeInvocationArguments(methodNode, False)
                    self.out.write("));\n")
                    self.out.write("}\n")
    
        # Write destructor
        self.out.write("OPENMM_EXPORT_AMOEBA void %s_destroy(%s* target) {\n" % (typeName, typeName))
        self.out.write("    delete reinterpret_cast<%s*>(target);\n" % className)
        self.out.write("}\n")

        # Record method names for future reference.
        methodNames = {}
        for methodNode in methodList:
            methodDefinition = getText("definition", methodNode)
            shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
            methodNames[methodNode] = shortMethodDefinition.split()[-1]
        
        # Write other methods
        for methodNode in methodList:
            methodName = methodNames[methodNode]
            if methodName in (shortClassName, destructorName):
                continue
            if self.shouldHideMethod(methodNode):
                continue
            isConstMethod = (methodNode.attrib['const'] == 'yes')
            if isConstMethod and any(methodNames[m] == methodName and m.attrib['const'] == 'no' for m in methodList):
                # There are two identical methods that differ only in whether they are const.  Skip the const one.
                continue
            methodType = getText("type", methodNode)
            returnType = self.getType(methodType)
            if methodType in self.classesByShortName:
                methodType = self.classesByShortName[methodType]
            self.out.write("OPENMM_EXPORT_AMOEBA %s %s_%s(" % (returnType, typeName, methodName))
            isInstanceMethod = (methodNode.attrib['static'] != 'yes')
            if isInstanceMethod:
                if isConstMethod:
                    self.out.write('const ')
                self.out.write("%s* target" % typeName)
            self.writeArguments(methodNode, isInstanceMethod)
            self.out.write(") {\n")
            self.out.write("    ")
            if returnType != 'void':
                if methodType.endswith('&'):
                    # Convert references to pointers
                    self.out.write('%s* result = &' % methodType[:-1].strip())
                else:
                    self.out.write('%s result = ' % methodType)
            if isInstanceMethod:
                self.out.write('reinterpret_cast<')
                if isConstMethod:
                    self.out.write('const ')
                self.out.write('%s*>(target)->' % className)
            else:
                self.out.write('%s::' % className)
            self.out.write('%s(' % methodName)
            self.writeInvocationArguments(methodNode, False)
            self.out.write(');\n')
            if returnType != 'void':
                self.out.write('    return %s;\n' % self.wrapValue(methodType, 'result'))
            self.out.write("}\n")
    
    def writeArguments(self, methodNode, initialSeparator):
        paramList = findNodes(methodNode, 'param')
        if initialSeparator:
            separator = ", "
        else:
            separator = ""
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'void':
                continue
            type = self.getType(type)
            name = getText('declname', node)
            self.out.write("%s%s %s" % (separator, type, name))
            separator = ", "
    
    def writeInvocationArguments(self, methodNode, initialSeparator):
        paramList = findNodes(methodNode, 'param')
        if initialSeparator:
            separator = ", "
        else:
            separator = ""
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'void':
                continue
            name = getText('declname', node)
            if self.getType(type) != type:
                name = self.unwrapValue(type, name)
            self.out.write("%s%s" % (separator, name))
            separator = ", "
    
    def getType(self, type):
        if type in self.typeTranslations:
            return self.typeTranslations[type]
        if type in self.typesByShortName:
            return self.typesByShortName[type]
        if type.startswith('const '):
            return 'const '+self.getType(type[6:].strip())
        if type.endswith('&') or type.endswith('*'):
            return self.getType(type[:-1].strip())+'*'
        return type
    
    def wrapValue(self, type, value):
        if type == 'bool':
            return '(%s ? OpenMM_True : OpenMM_False)' % value
        if type == 'std::string':
            return '%s.c_str()' % value
        if type == 'const std::string &':
            return '%s->c_str()' % value
        if type in self.enumerationTypes:
            return 'static_cast<%s>(%s)' % (self.enumerationTypes[type], value)
        wrappedType = self.getType(type)
        if wrappedType == type:
            return value;
        if type.endswith('*') or type.endswith('&'):
            return 'reinterpret_cast<%s>(%s)' % (wrappedType, value)
        return 'static_cast<%s>(%s)' % (wrappedType, value)
    
    def unwrapValue(self, type, value):
        if type.endswith('&'):
            unwrappedType = type[:-1].strip()
            if unwrappedType in self.classesByShortName:
                unwrappedType  = self.classesByShortName[unwrappedType]
            if unwrappedType == 'const std::string':
                return 'std::string(%s)' % value
            return '*'+self.unwrapValue(unwrappedType+'*', value)
        if type in self.classesByShortName:
            return 'static_cast<%s>(%s)' % (self.classesByShortName[type], value)
        if type == 'bool':
            return value
        return 'reinterpret_cast<%s>(%s)' % (type, value)

    def writeOutput(self):
        print("""
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "../../../wrappers/OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" {

/* OpenMM_2D_IntArray */
OPENMM_EXPORT_AMOEBA OpenMM_2D_IntArray* OpenMM_2D_IntArray_create(int size) {
    return reinterpret_cast<OpenMM_2D_IntArray*>(new vector<vector<int> >(size));
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_destroy(OpenMM_2D_IntArray* array) {
    delete reinterpret_cast<vector<vector<int> >*>(array);
}
OPENMM_EXPORT_AMOEBA int OpenMM_2D_IntArray_getSize(const OpenMM_2D_IntArray* array) {
    return reinterpret_cast<const vector<vector<int> >*>(array)->size();
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_resize(OpenMM_2D_IntArray* array, int size) {
    reinterpret_cast<vector<vector<int> >*>(array)->resize(size);
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_append(OpenMM_2D_IntArray* array, int index1, int value) {
    vector<vector<int> >* array2DInt = reinterpret_cast<vector<vector<int> >*>(array);
    if (array2DInt->size() <= index1) {
        array2DInt->resize(index1+1);
    }
    (*array2DInt)[index1].push_back(value);
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_set(OpenMM_2D_IntArray* array, int index1, int index2, int value) {
    vector<vector<int> >* array2DInt = reinterpret_cast<vector<vector<int> >*>(array);
    if (array2DInt->size() <= index1) {
        array2DInt->resize(index1+1);
    }
    if (array2DInt[index1].size() <= index2) {
        array2DInt[index1].resize(index2+1);
    }
    (*array2DInt)[index1][index2] = value;
}
OPENMM_EXPORT_AMOEBA void OpenMM_2D_IntArray_get(const OpenMM_2D_IntArray* array, int index1, int index2, int* value) {
    const vector<vector<int> >* array2DInt = reinterpret_cast<const vector<vector<int> >*>(array);
    if (array2DInt->size() <= index1)
        throw OpenMMException("OpenMM_2D_IntArray_get: first index out of range.");

    if ((*array2DInt)[index1].size() <= index2)
        throw OpenMMException("OpenMM_2D_IntArray_get: second index out of range.");
    *value = (*array2DInt)[index1][index2];
}

/* OpenMM_3D_DoubleArray */
OPENMM_EXPORT_AMOEBA OpenMM_3D_DoubleArray* OpenMM_3D_DoubleArray_create(int size1, int size2, int size3) {
    int ii, jj;  
    std::vector< std::vector< std::vector<double> > >* v3D_Array = new std::vector<std::vector<std::vector<double> > >(size1);

    for (ii = 0; ii < size1; ii++) {
        (*v3D_Array)[ii].resize(size2);
        for (jj = 0; jj < size2; jj++) {
           (*v3D_Array)[ii][jj].resize(size3);
        }    
    }    
    return reinterpret_cast<OpenMM_3D_DoubleArray*>(v3D_Array);
}

OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_set(OpenMM_3D_DoubleArray* array, int index1, int index2, OpenMM_DoubleArray* values) {
    unsigned int ii;
    std::vector< std::vector< std::vector<double> > >* v3D_Array = reinterpret_cast<std::vector<std::vector<std::vector<double> > >*>(array);
    std::vector<double> * value_array                            = reinterpret_cast<std::vector<double> *>(values);
    for (ii = 0; ii < (*value_array).size(); ii++) {
        (*v3D_Array)[index1][index2][ii] = (*value_array)[ii];
    }    
}

OPENMM_EXPORT_AMOEBA void OpenMM_3D_DoubleArray_destroy(OpenMM_3D_DoubleArray* array) {
    delete reinterpret_cast<std::vector<std::vector<std::vector<double> > >*>(array);
}""", file=self.out)

        self.writeClasses()
        print("}\n", file=self.out)

class FortranHeaderGenerator(WrapperGenerator):
    """This class generates the header file for the Fortran API wrappers."""
    
    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
        self.typeTranslations = {'int': 'integer*4',
                                 'bool': 'integer*4',
                                 'double': 'real*8',
                                 'char *': 'character(*)',
                                 'const char *': 'character(*)',
                                 'std::string': 'character(*)',
                                 'const std::string &': 'character(*)',
                                 'std::vector< std::string >': 'type (OpenMM_StringArray)',
                                 'std::vector< Vec3 >': 'type (OpenMM_Vec3Array)',
                                 'std::vector< std::pair< int, int > >': 'type (OpenMM_BondArray)',
                                 'std::map< std::string, double >': 'type (OpenMM_ParameterArray)',
                                 'std::map< std::string, std::string >': 'type (OpenMM_PropertyArray)',
                                 'std::vector< double >': 'type (OpenMM_DoubleArray)',
                                 'std::vector< int >': 'type (OpenMM_IntArray)',
                                 'std::set< int >': 'type (OpenMM_IntSet)'}
        self.enumerationTypes = set()
    
    def writeGlobalConstants(self):
        self.out.write("    ! Global Constants\n\n")
        node = next((x for x in findNodes(self.doc.getroot(), "compounddef", kind="namespace") if x.findtext("compoundname") == "OpenMM"))
        for section in findNodes(node, "sectiondef", kind="var"):
            for memberNode in findNodes(section, "memberdef", kind="variable", mutable="no", prot="public", static="yes"):
                vDef = convertOpenMMPrefix(getText("name", memberNode))
                iDef = getText("initializer", memberNode)
                if iDef.startswith("="):
                    iDef = iDef[1:]
                self.out.write("    real*8, parameter :: %s = %s\n" % (vDef, iDef))

    def writeTypeDeclarations(self):
        self.out.write("\n    ! Type Declarations\n")
        for classNode in self._orderedClassNodes:
            className = getText("compoundname", classNode)
            shortName = stripOpenMMPrefix(className)
            typeName = convertOpenMMPrefix(className)
            self.out.write("\n    type %s\n" % typeName)
            self.out.write("        integer*8 :: handle = 0\n")
            self.out.write("    end type\n")
            self.typesByShortName[shortName] = typeName

    def writeClasses(self):
        for classNode in self._orderedClassNodes:
            className = getText("compoundname", classNode)
            self.out.write("\n        ! %s\n" % className)
            self.writeMethods(classNode)
        self.out.write("\n")

    def writeEnumerations(self, classNode):
        enumNodes = []
        for section in findNodes(classNode, "sectiondef", kind="public-type"):
            for node in findNodes(section, "memberdef", kind="enum", prot="public"):
                enumNodes.append(node)
        className = getText("compoundname", classNode)
        typeName = convertOpenMMPrefix(className)
        for enumNode in enumNodes:
            for valueNode in findNodes(enumNode, "enumvalue", prot="public"):
                vName = convertOpenMMPrefix(getText("name", valueNode))
                vInit = getText("initializer", valueNode)
                if vInit.startswith("="):
                    vInit = vInit[1:].strip()
                self.out.write("    integer*4, parameter :: %s_%s = %s\n" % (typeName, vName, vInit))
            enumName = getText("name", enumNode)
            enumTypeName = "%s_%s" % (typeName, enumName)
            self.typesByShortName[enumName] = enumTypeName
            self.enumerationTypes.add(enumName)
        if len(enumNodes)>0: self.out.write("\n")

    def writeMethods(self, classNode):
        methodList = self.getClassMethods(classNode)
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        destructorName = '~'+shortClassName

        if not ('abstract' in classNode.attrib and classNode.attrib['abstract'] == 'yes'):
            # Write constructors
            numConstructors = 0
            for methodNode in methodList:
                methodDefinition = getText("definition", methodNode)
                shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
                methodName = shortMethodDefinition.split()[-1]
                if methodName == shortClassName:
                    if self.shouldHideMethod(methodNode):
                        continue
                    numConstructors += 1
                    if numConstructors == 1:
                        suffix = ""
                    else:
                        suffix = "_%d" % numConstructors
                    self.out.write("        subroutine %s_create%s(result" % (typeName, suffix))
                    self.writeArguments(methodNode, True)
                    self.out.write(")\n")
                    self.out.write("            use OpenMM_Types; implicit none\n")
                    self.out.write("            type (%s) result\n" % typeName)
                    self.declareArguments(methodNode)
                    self.out.write("        end subroutine\n")
    
        # Write destructor
        self.out.write("        subroutine %s_destroy(destroy)\n" % typeName)
        self.out.write("            use OpenMM_Types; implicit none\n")
        self.out.write("            type (%s) destroy\n" % typeName)
        self.out.write("        end subroutine\n")

        # Record method names for future reference.
        methodNames = {}
        for methodNode in methodList:
            methodDefinition = getText("definition", methodNode)
            shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
            methodNames[methodNode] = shortMethodDefinition.split()[-1]
        
        # Write other methods
        for methodNode in methodList:
            methodName = methodNames[methodNode]
            if methodName in (shortClassName, destructorName):
                continue
            if self.shouldHideMethod(methodNode):
                continue
            isConstMethod = (methodNode.attrib['const'] == 'yes')
            if isConstMethod and any(methodNames[m] == methodName and m.attrib['const'] == 'no' for m in methodList):
                # There are two identical methods that differ only in whether they are const.  Skip the const one.
                continue
            returnType = self.getType(getText("type", methodNode))
            hasReturnValue = (returnType in ('integer*4', 'real*8'))
            hasReturnArg = not (hasReturnValue or returnType == 'void')
            functionName = "%s_%s" % (typeName, methodName)
            if hasReturnValue:
                self.out.write("        function ")
            else:
                self.out.write("        subroutine ")
            self.out.write("%s(" % functionName)
            isInstanceMethod = (methodNode.attrib['static'] != 'yes')
            if isInstanceMethod:
                self.out.write("target")
            numArgs = self.writeArguments(methodNode, isInstanceMethod)
            if hasReturnArg:
                if isInstanceMethod or numArgs > 0:
                    self.out.write(", ")
                self.out.write("result")
            self.out.write(")\n")
            self.out.write("            use OpenMM_Types; implicit none\n")
            self.out.write("            type (%s) target\n" % typeName)
            self.declareArguments(methodNode)
            if hasReturnValue:
                self.declareOneArgument(returnType, functionName)
            if hasReturnArg:
                self.declareOneArgument(returnType, 'result')
            if hasReturnValue:
                self.out.write("        end function\n")
            else:
                self.out.write("        end subroutine\n")
    
    def writeArguments(self, methodNode, initialSeparator):
        paramList = findNodes(methodNode, 'param')
        if initialSeparator:
            separator = ", "
        else:
            separator = ""
        numArgs = 0
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'void':
                continue
            name = getText('declname', node)
            self.out.write("%s%s" % (separator, name))
            separator = ", &\n"
            numArgs += 1
        return numArgs
    
    def declareOneArgument(self, type, name):
        if type == 'void':
            return
        type = self.getType(type)
        if type == 'Vec3':
            self.out.write("            real*8 %s(3)\n" % name)
        else:
            self.out.write("            %s %s\n" % (type, name))
    
    def declareArguments(self, methodNode):
        paramList = findNodes(methodNode, 'param')
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            name = getText('declname', node)
            self.declareOneArgument(type, name)
    
    def getType(self, type):
        if type in self.typeTranslations:
            return self.typeTranslations[type]
        if type in self.typesByShortName:
            return 'type (%s)' % self.typesByShortName[type]
        if type.startswith('const '):
            return self.getType(type[6:].strip())
        if type.endswith('&') or type.endswith('*'):
            return self.getType(type[:-1].strip())
        return type

    def writeOutput(self):
        print("""
MODULE OpenMM_Types
    implicit none
""", file=self.out)
        self.writeGlobalConstants()
        self.writeTypeDeclarations()
        print("""
    ! Enumerations

    integer*4, parameter :: OpenMM_False = 0
    integer*4, parameter :: OpenMM_True = 1""", file=self.out)

        for classNode in self._orderedClassNodes:
            self.writeEnumerations(classNode)
        print("""
END MODULE OpenMM_Types

MODULE OpenMM
    use OpenMM_Types; implicit none
    interface
""", file=self.out)
        
        self.writeClasses()
        
        print("""
    end interface
END MODULE OpenMM""", file=self.out)


class FortranSourceGenerator(WrapperGenerator):
    """This class generates the source file for the Fortran API wrappers."""

    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
        self.typeTranslations = {'bool': 'OpenMM_Boolean',
                                 'Vec3': 'OpenMM_Vec3',
                                 'Context': 'OpenMM_Context',
                                 'std::string': 'char*',
                                 'const std::string &': 'const char*',
                                 'std::vector< std::string >': 'OpenMM_StringArray',
                                 'std::vector< Vec3 >': 'OpenMM_Vec3Array',
                                 'std::vector< std::pair< int, int > >': 'OpenMM_BondArray',
                                 'std::map< std::string, double >': 'OpenMM_ParameterArray',
                                 'std::map< std::string, std::string >': 'OpenMM_PropertyArray',
                                 'std::vector< double >': 'OpenMM_DoubleArray',
                                 'std::vector< int >': 'OpenMM_IntArray',
                                 'std::set< int >': 'OpenMM_IntSet',
                                 'std::vector< std::vector< int > >': 'OpenMM_2D_IntArray',
                                 'std::vector< std::vector< std::vector< double > > >': 'OpenMM_3D_DoubleArray'}
        self.inverseTranslations = dict((self.typeTranslations[key], key) for key in self.typeTranslations)
        self.classesByShortName = {}
        self.enumerationTypes = {}
        self.findTypes()
    
    def findTypes(self):
        for classNode in self._orderedClassNodes:
            className = getText("compoundname", classNode)
            shortName = stripOpenMMPrefix(className)
            typeName = convertOpenMMPrefix(className)
            self.typesByShortName[shortName] = typeName
            self.classesByShortName[shortName] = className

    def findEnumerations(self, classNode):
        enumNodes = []
        for section in findNodes(classNode, "sectiondef", kind="public-type"):
            for node in findNodes(section, "memberdef", kind="enum", prot="public"):
                enumNodes.append(node)
        className = getText("compoundname", classNode)
        typeName = convertOpenMMPrefix(className)
        for enumNode in enumNodes:
            enumName = getText("name", enumNode)
            enumTypeName = "%s_%s" % (typeName, enumName)
            enumClassName = "%s::%s" % (className, enumName)
            self.typesByShortName[enumName] = enumTypeName
            self.classesByShortName[enumName] = enumClassName
            self.enumerationTypes[enumClassName] = enumTypeName

    def writeClasses(self):
        for classNode in self._orderedClassNodes:
            className = stripOpenMMPrefix(getText("compoundname", classNode))
            self.out.write("\n/* OpenMM::%s */\n" % className)
            self.findEnumerations(classNode)
            self.writeMethods(classNode)
        self.out.write("\n")

    def writeMethods(self, classNode):
        methodList = self.getClassMethods(classNode)
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        destructorName = '~'+shortClassName

        if not ('abstract' in classNode.attrib and classNode.attrib['abstract'] == 'yes'):
            # Write constructors
            numConstructors = 0
            for methodNode in methodList:
                methodDefinition = getText("definition", methodNode)
                shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
                methodName = shortMethodDefinition.split()[-1]
                if methodName == shortClassName:
                    if self.shouldHideMethod(methodNode):
                        continue
                    numConstructors += 1
                    if numConstructors == 1:
                        suffix = ""
                    else:
                        suffix = "_%d" % numConstructors
                    functionName = "%s_create%s" % (typeName, suffix)
                    self.writeOneConstructor(classNode, methodNode, functionName, functionName.lower()+'_')
                    self.writeOneConstructor(classNode, methodNode, functionName, functionName.upper())
    
        # Write destructor
        functionName = "%s_destroy" % typeName
        self.writeOneDestructor(typeName, functionName.lower()+'_')
        self.writeOneDestructor(typeName, functionName.upper())

        # Record method names for future reference.
        methodNames = {}
        for methodNode in methodList:
            methodDefinition = getText("definition", methodNode)
            shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
            methodNames[methodNode] = shortMethodDefinition.split()[-1]
        
        # Write other methods
        for methodNode in methodList:
            methodName = methodNames[methodNode]
            if methodName in (shortClassName, destructorName):
                continue
            if '~' in methodName:
                print('***', methodName, destructorName)
            if self.shouldHideMethod(methodNode):
                continue
            isConstMethod = (methodNode.attrib['const'] == 'yes')
            if isConstMethod and any(methodNames[m] == methodName and m.attrib['const'] == 'no' for m in methodList):
                # There are two identical methods that differ only in whether they are const.  Skip the const one.
                continue
            functionName = "%s_%s" % (typeName, methodName)
            self.writeOneMethod(classNode, methodNode, functionName, functionName.lower()+'_')
            self.writeOneMethod(classNode, methodNode, functionName, functionName.upper())
    
    def writeOneConstructor(self, classNode, methodNode, functionName, wrapperFunctionName):
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        self.out.write("OPENMM_EXPORT_AMOEBA void %s(%s*& result" % (wrapperFunctionName, typeName))
        self.writeArguments(methodNode, True)
        self.out.write(") {\n")
        self.out.write("    result = %s(" % functionName)
        self.writeInvocationArguments(methodNode, False)
        self.out.write(");\n")
        self.out.write("}\n")
    
    def writeOneDestructor(self, typeName, wrapperFunctionName):
        self.out.write("OPENMM_EXPORT_AMOEBA void %s(%s*& destroy) {\n" % (wrapperFunctionName, typeName))
        self.out.write("    %s_destroy(destroy);\n" % typeName)
        self.out.write("    destroy = 0;\n")
        self.out.write("}\n")
    
    def writeOneMethod(self, classNode, methodNode, methodName, wrapperFunctionName):
        className = getText("compoundname", classNode)
        typeName = convertOpenMMPrefix(className)


        isConstMethod = (methodNode.attrib['const'] == 'yes')
        methodType = getText("type", methodNode)
        returnType = self.getType(methodType)
        hasReturnValue = (returnType in ('int', 'bool', 'double'))
        hasReturnArg = not (hasReturnValue or returnType == 'void')
        if methodType in self.classesByShortName:
            methodType = self.classesByShortName[methodType]
        self.out.write("OPENMM_EXPORT_AMOEBA ")
        if hasReturnValue:
            self.out.write(returnType)
        else:
            self.out.write('void')
        self.out.write(" %s(" % wrapperFunctionName)
        isInstanceMethod = (methodNode.attrib['static'] != 'yes')
        if isInstanceMethod:
            if isConstMethod:
                self.out.write('const ')
            self.out.write("%s*& target" % typeName)
        returnArg = None
        if hasReturnArg:
            if returnType == 'const char*':
                # We need a non-const buffer to copy the result into
                returnArg = 'char* result'
            else:
                returnArg = "%s& result" % returnType
        numArgs = self.writeArguments(methodNode, isInstanceMethod, returnArg)
        if hasReturnArg and returnType == 'const char*':
            self.out.write(", int result_length")
        self.out.write(") {\n")
        self.out.write("    ")
        if hasReturnValue:
            self.out.write("return ")
        if hasReturnArg:
            if returnType == 'const char*':
                self.out.write("const char* result_chars = ")
            else:
                self.out.write("result = ")
        self.out.write("%s(" % methodName)
        if isInstanceMethod:
            self.out.write("target")
        self.writeInvocationArguments(methodNode, isInstanceMethod)
        self.out.write(');\n')
        if hasReturnArg and returnType == 'const char*':
            self.out.write("    copyAndPadString(result, result_chars, result_length);\n")
        self.out.write("}\n")
    
    def writeArguments(self, methodNode, initialSeparator, extraArg=None):
        paramList = findNodes(methodNode, 'param')
        if initialSeparator:
            separator = ", "
        else:
            separator = ""
        numArgs = 0
        
        # Write the arguments.
        
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'void':
                continue
            type = self.getType(type)
            if self.isHandleType(type):
                type = type+'&'
            elif type[-1] not in ('&', '*'):
                type = type+' const&'
            name = getText('declname', node)
            self.out.write("%s%s %s" % (separator, type, name))
            separator = ", "
            numArgs += 1
        
        # If an extra argument is needed for the return value, write it.
        
        if extraArg is not None:
            self.out.write("%s%s" % (separator, extraArg))
            separator = ", "
            numArgs += 1
        
        # Write length arguments for strings.
        
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'const std::string &':
                name = getText('declname', node)
                self.out.write(", int %s_length" % name)
                numArgs += 1
        return numArgs
    
    def writeInvocationArguments(self, methodNode, initialSeparator):
        paramList = findNodes(methodNode, 'param')
        if initialSeparator:
            separator = ", "
        else:
            separator = ""
        for node in paramList:
            try:
                type = getText('type', node)
            except IndexError:
                type = getText('type/ref', node)
            if type == 'void':
                continue
            name = getText('declname', node)
            if type == 'const std::string &':
                name = 'makeString(%s, %s_length).c_str()' % (name, name)
            self.out.write("%s%s" % (separator, name))
            separator = ", "
    
    def getType(self, type):
        if type in self.typeTranslations:
            return self.typeTranslations[type]
        if type in self.typesByShortName:
            return self.typesByShortName[type]
        if type.startswith('const '):
            return 'const '+self.getType(type[6:].strip())
        if type.endswith('&') or type.endswith('*'):
            return self.getType(type[:-1].strip())+'*'
        return type
    
    def isHandleType(self, type):
        if type == 'OpenMM_Vec3':
            return False
        if type.endswith('*') or type.endswith('&'):
            return self.isHandleType(type[:-1].strip())
        if type.startswith('const '):
            return self.isHandleType(type[6:].strip())
        if type.startswith('OpenMM_'):
            return True;
        return False

    def writeOutput(self):
        print("""
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "../../../wrappers/OpenMMCWrapper.h"
#include "AmoebaOpenMMCWrapper.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

/* Utilities for dealing with Fortran's blank-padded strings. */
static void copyAndPadString(char* dest, const char* source, int length) {
    bool reachedEnd = false;
    for (int i = 0; i < length; i++) {
        if (source[i] == 0)
            reachedEnd = true;
        dest[i] = (reachedEnd ? ' ' : source[i]);
    }
}

static string makeString(const char* fsrc, int length) {
    while (length && fsrc[length-1]==' ')
        --length;
    return string(fsrc, length);
}

extern "C" {
""", file=self.out)

        self.writeClasses()
        print("}", file=self.out)

inputDirname = sys.argv[1]
builder = CHeaderGenerator(inputDirname, open(os.path.join(sys.argv[2], 'AmoebaOpenMMCWrapper.h'), 'w'))
builder.writeOutput()
builder = CSourceGenerator(inputDirname, open(os.path.join(sys.argv[2], 'AmoebaOpenMMCWrapper.cpp'), 'w'))
builder.writeOutput()
builder = FortranHeaderGenerator(inputDirname, open(os.path.join(sys.argv[2], 'AmoebaOpenMMFortranModule.f90'), 'w'))
builder.writeOutput()
builder = FortranSourceGenerator(inputDirname, open(os.path.join(sys.argv[2], 'AmoebaOpenMMFortranWrapper.cpp'), 'w'))
builder.writeOutput()
