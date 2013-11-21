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
        self.skipClasses = ['OpenMM::Vec3', 'OpenMM::XmlSerializer', 'OpenMM::Kernel', 'OpenMM::KernelImpl', 'OpenMM::KernelFactory', 'OpenMM::ContextImpl', 'OpenMM::SerializationNode', 'OpenMM::SerializationProxy']
        self.skipMethods = ['OpenMM::Context::getState', 'OpenMM::Platform::loadPluginsFromDirectory', 'OpenMM::Context::createCheckpoint', 'OpenMM::Context::loadCheckpoint']
        self.hideClasses = ['Kernel', 'KernelImpl', 'KernelFactory', 'ContextImpl', 'SerializationNode', 'SerializationProxy']
        self.typeTranslations = {'bool': 'OpenMM_Boolean',
                                 'Vec3': 'OpenMM_Vec3',
                                 'std::string': 'char*',
                                 'const std::string &': 'const char*',
                                 'std::vector< std::string >': 'OpenMM_StringArray',
                                 'std::vector< Vec3 >': 'OpenMM_Vec3Array',
                                 'std::vector< std::pair< int, int > >': 'OpenMM_BondArray',
                                 'std::map< std::string, double >': 'OpenMM_ParameterArray',
                                 'std::map< std::string, std::string >': 'OpenMM_PropertyArray',
                                 'std::vector< double >': 'OpenMM_DoubleArray',
                                 'std::vector< int >': 'OpenMM_IntArray',
                                 'std::set< int >': 'OpenMM_IntSet'}
        self.inverseTranslations = dict((self.typeTranslations[key], key) for key in self.typeTranslations)
        self.nodeByID={}

        # Read all the XML files and merge them into a single document.
        self.doc = etree.ElementTree(etree.Element('root'))
        for file in os.listdir(inputDirname):
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
                    self.out.write("extern OPENMM_EXPORT %s* %s_create%s(" % (typeName, typeName, suffix))
                    self.writeArguments(methodNode, False)
                    self.out.write(");\n")
    
            # Write destructor
            self.out.write("extern OPENMM_EXPORT void %s_destroy(%s* target);\n" % (typeName, typeName))

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
            self.out.write("extern OPENMM_EXPORT %s %s_%s(" % (returnType, typeName, methodName))
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
        print >>self.out, """
#ifndef OPENMM_CWRAPPER_H_
#define OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT
#define OPENMM_EXPORT
#endif
"""
        self.writeGlobalConstants()
        self.writeTypeDeclarations()
        print >>self.out, """
typedef struct OpenMM_Vec3Array_struct OpenMM_Vec3Array;
typedef struct OpenMM_StringArray_struct OpenMM_StringArray;
typedef struct OpenMM_BondArray_struct OpenMM_BondArray;
typedef struct OpenMM_ParameterArray_struct OpenMM_ParameterArray;
typedef struct OpenMM_PropertyArray_struct OpenMM_PropertyArray;
typedef struct OpenMM_DoubleArray_struct OpenMM_DoubleArray;
typedef struct OpenMM_IntArray_struct OpenMM_IntArray;
typedef struct OpenMM_IntSet_struct OpenMM_IntSet;
typedef struct {double x, y, z;} OpenMM_Vec3;

typedef enum {OpenMM_False = 0, OpenMM_True = 1} OpenMM_Boolean;

#if defined(__cplusplus)
extern "C" {
#endif

/* OpenMM_Vec3 */
extern OPENMM_EXPORT OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale);

/* OpenMM_Vec3Array */
extern OPENMM_EXPORT OpenMM_Vec3Array* OpenMM_Vec3Array_create(int size);
extern OPENMM_EXPORT void OpenMM_Vec3Array_destroy(OpenMM_Vec3Array* array);
extern OPENMM_EXPORT int OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array* array);
extern OPENMM_EXPORT void OpenMM_Vec3Array_resize(OpenMM_Vec3Array* array, int size);
extern OPENMM_EXPORT void OpenMM_Vec3Array_append(OpenMM_Vec3Array* array, const OpenMM_Vec3 vec);
extern OPENMM_EXPORT void OpenMM_Vec3Array_set(OpenMM_Vec3Array* array, int index, const OpenMM_Vec3 vec);
extern OPENMM_EXPORT const OpenMM_Vec3* OpenMM_Vec3Array_get(const OpenMM_Vec3Array* array, int index);

/* OpenMM_StringArray */
extern OPENMM_EXPORT OpenMM_StringArray* OpenMM_StringArray_create(int size);
extern OPENMM_EXPORT void OpenMM_StringArray_destroy(OpenMM_StringArray* array);
extern OPENMM_EXPORT int OpenMM_StringArray_getSize(const OpenMM_StringArray* array);
extern OPENMM_EXPORT void OpenMM_StringArray_resize(OpenMM_StringArray* array, int size);
extern OPENMM_EXPORT void OpenMM_StringArray_append(OpenMM_StringArray* array, const char* string);
extern OPENMM_EXPORT void OpenMM_StringArray_set(OpenMM_StringArray* array, int index, const char* string);
extern OPENMM_EXPORT const char* OpenMM_StringArray_get(const OpenMM_StringArray* array, int index);

/* OpenMM_BondArray */
extern OPENMM_EXPORT OpenMM_BondArray* OpenMM_BondArray_create(int size);
extern OPENMM_EXPORT void OpenMM_BondArray_destroy(OpenMM_BondArray* array);
extern OPENMM_EXPORT int OpenMM_BondArray_getSize(const OpenMM_BondArray* array);
extern OPENMM_EXPORT void OpenMM_BondArray_resize(OpenMM_BondArray* array, int size);
extern OPENMM_EXPORT void OpenMM_BondArray_append(OpenMM_BondArray* array, int particle1, int particle2);
extern OPENMM_EXPORT void OpenMM_BondArray_set(OpenMM_BondArray* array, int index, int particle1, int particle2);
extern OPENMM_EXPORT void OpenMM_BondArray_get(const OpenMM_BondArray* array, int index, int* particle1, int* particle2);

/* OpenMM_ParameterArray */
extern OPENMM_EXPORT int OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray* array);
extern OPENMM_EXPORT double OpenMM_ParameterArray_get(const OpenMM_ParameterArray* array, const char* name);

/* OpenMM_PropertyArray */
extern OPENMM_EXPORT int OpenMM_PropertyArray_getSize(const OpenMM_PropertyArray* array);
extern OPENMM_EXPORT const char* OpenMM_PropertyArray_get(const OpenMM_PropertyArray* array, const char* name);"""

        for type in ('double', 'int'):
            name = 'OpenMM_%sArray' % type.capitalize()
            values = {'type':type, 'name':name}
            print >>self.out, """
/* %(name)s */
extern OPENMM_EXPORT %(name)s* %(name)s_create(int size);
extern OPENMM_EXPORT void %(name)s_destroy(%(name)s* array);
extern OPENMM_EXPORT int %(name)s_getSize(const %(name)s* array);
extern OPENMM_EXPORT void %(name)s_resize(%(name)s* array, int size);
extern OPENMM_EXPORT void %(name)s_append(%(name)s* array, %(type)s value);
extern OPENMM_EXPORT void %(name)s_set(%(name)s* array, int index, %(type)s value);
extern OPENMM_EXPORT %(type)s %(name)s_get(const %(name)s* array, int index);""" % values

        for type in ('int',):
            name = 'OpenMM_%sSet' % type.capitalize()
            values = {'type':type, 'name':name}
            print >>self.out, """
/* %(name)s */
extern OPENMM_EXPORT %(name)s* %(name)s_create();
extern OPENMM_EXPORT void %(name)s_destroy(%(name)s* set);
extern OPENMM_EXPORT int %(name)s_getSize(const %(name)s* set);
extern OPENMM_EXPORT void %(name)s_insert(%(name)s* set, %(type)s value);""" % values

        print >>self.out, """
/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
extern OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types, int enforcePeriodicBox);
extern OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory);"""

        self.writeClasses()

        print >>self.out, """
#if defined(__cplusplus)
}
#endif

#endif /*OPENMM_CWRAPPER_H_*/"""


class CSourceGenerator(WrapperGenerator):
    """This class generates the source file for the C API wrappers."""

    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
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
                    self.out.write("OPENMM_EXPORT %s* %s_create%s(" % (typeName, typeName, suffix))
                    self.writeArguments(methodNode, False)
                    self.out.write(") {\n")
                    self.out.write("    return reinterpret_cast<%s*>(new %s(" % (typeName, className))
                    self.writeInvocationArguments(methodNode, False)
                    self.out.write("));\n")
                    self.out.write("}\n")
    
            # Write destructor
            self.out.write("OPENMM_EXPORT void %s_destroy(%s* target) {\n" % (typeName, typeName))
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
            self.out.write("OPENMM_EXPORT %s %s_%s(" % (returnType, typeName, methodName))
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
            return '*'+self.unwrapValue(unwrappedType+'*', value)
        if type in self.classesByShortName:
            return 'static_cast<%s>(%s)' % (self.classesByShortName[type], value)
        if type == 'bool':
            return value
        return 'reinterpret_cast<%s>(%s)' % (type, value)

    def writeOutput(self):
        print >>self.out, """
#include "OpenMM.h"
#include "OpenMMCWrapper.h"
#include <cstring>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" {

/* OpenMM_Vec3 */
OPENMM_EXPORT OpenMM_Vec3 OpenMM_Vec3_scale(const OpenMM_Vec3 vec, double scale) {
    OpenMM_Vec3 result = {vec.x*scale, vec.y*scale, vec.z*scale};
    return result;
}

/* OpenMM_Vec3Array */
OPENMM_EXPORT OpenMM_Vec3Array* OpenMM_Vec3Array_create(int size) {
    return reinterpret_cast<OpenMM_Vec3Array*>(new vector<Vec3>(size));
}
OPENMM_EXPORT void OpenMM_Vec3Array_destroy(OpenMM_Vec3Array* array) {
    delete reinterpret_cast<vector<Vec3>*>(array);
}
OPENMM_EXPORT int OpenMM_Vec3Array_getSize(const OpenMM_Vec3Array* array) {
    return reinterpret_cast<const vector<Vec3>*>(array)->size();
}
OPENMM_EXPORT void OpenMM_Vec3Array_resize(OpenMM_Vec3Array* array, int size) {
    reinterpret_cast<vector<Vec3>*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_Vec3Array_append(OpenMM_Vec3Array* array, const OpenMM_Vec3 vec) {
    reinterpret_cast<vector<Vec3>*>(array)->push_back(Vec3(vec.x, vec.y, vec.z));
}
OPENMM_EXPORT void OpenMM_Vec3Array_set(OpenMM_Vec3Array* array, int index, const OpenMM_Vec3 vec) {
    (*reinterpret_cast<vector<Vec3>*>(array))[index] = Vec3(vec.x, vec.y, vec.z);
}
OPENMM_EXPORT const OpenMM_Vec3* OpenMM_Vec3Array_get(const OpenMM_Vec3Array* array, int index) {
    return reinterpret_cast<const OpenMM_Vec3*>((&(*reinterpret_cast<const vector<Vec3>*>(array))[index]));
}

/* OpenMM_StringArray */
OPENMM_EXPORT OpenMM_StringArray* OpenMM_StringArray_create(int size) {
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(size));
}
OPENMM_EXPORT void OpenMM_StringArray_destroy(OpenMM_StringArray* array) {
    delete reinterpret_cast<vector<string>*>(array);
}
OPENMM_EXPORT int OpenMM_StringArray_getSize(const OpenMM_StringArray* array) {
    return reinterpret_cast<const vector<string>*>(array)->size();
}
OPENMM_EXPORT void OpenMM_StringArray_resize(OpenMM_StringArray* array, int size) {
    reinterpret_cast<vector<string>*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_StringArray_append(OpenMM_StringArray* array, const char* str) {
    reinterpret_cast<vector<string>*>(array)->push_back(string(str));
}
OPENMM_EXPORT void OpenMM_StringArray_set(OpenMM_StringArray* array, int index, const char* str) {
    (*reinterpret_cast<vector<string>*>(array))[index] = string(str);
}
OPENMM_EXPORT const char* OpenMM_StringArray_get(const OpenMM_StringArray* array, int index) {
    return (*reinterpret_cast<const vector<string>*>(array))[index].c_str();
}

/* OpenMM_BondArray */
OPENMM_EXPORT OpenMM_BondArray* OpenMM_BondArray_create(int size) {
    return reinterpret_cast<OpenMM_BondArray*>(new vector<pair<int, int> >(size));
}
OPENMM_EXPORT void OpenMM_BondArray_destroy(OpenMM_BondArray* array) {
    delete reinterpret_cast<vector<pair<int, int> >*>(array);
}
OPENMM_EXPORT int OpenMM_BondArray_getSize(const OpenMM_BondArray* array) {
    return reinterpret_cast<const vector<pair<int, int> >*>(array)->size();
}
OPENMM_EXPORT void OpenMM_BondArray_resize(OpenMM_BondArray* array, int size) {
    reinterpret_cast<vector<pair<int, int> >*>(array)->resize(size);
}
OPENMM_EXPORT void OpenMM_BondArray_append(OpenMM_BondArray* array, int particle1, int particle2) {
    reinterpret_cast<vector<pair<int, int> >*>(array)->push_back(pair<int, int>(particle1, particle2));
}
OPENMM_EXPORT void OpenMM_BondArray_set(OpenMM_BondArray* array, int index, int particle1, int particle2) {
    (*reinterpret_cast<vector<pair<int, int> >*>(array))[index] = pair<int, int>(particle1, particle2);
}
OPENMM_EXPORT void OpenMM_BondArray_get(const OpenMM_BondArray* array, int index, int* particle1, int* particle2) {
    pair<int, int> particles = (*reinterpret_cast<const vector<pair<int, int> >*>(array))[index];
    *particle1 = particles.first;
    *particle2 = particles.second;
}

/* OpenMM_ParameterArray */
OPENMM_EXPORT int OpenMM_ParameterArray_getSize(const OpenMM_ParameterArray* array) {
    return reinterpret_cast<const map<string, double>*>(array)->size();
}
OPENMM_EXPORT double OpenMM_ParameterArray_get(const OpenMM_ParameterArray* array, const char* name) {
    const map<string, double>* params = reinterpret_cast<const map<string, double>*>(array);
    const map<string, double>::const_iterator iter = params->find(string(name));
    if (iter == params->end())
        throw OpenMMException("OpenMM_ParameterArray_get: No such parameter");
    return iter->second;
}

/* OpenMM_PropertyArray */
OPENMM_EXPORT int OpenMM_PropertyArray_getSize(const OpenMM_PropertyArray* array) {
    return reinterpret_cast<const map<string, double>*>(array)->size();
}
OPENMM_EXPORT const char* OpenMM_PropertyArray_get(const OpenMM_PropertyArray* array, const char* name) {
    const map<string, string>* params = reinterpret_cast<const map<string, string>*>(array);
    const map<string, string>::const_iterator iter = params->find(string(name));
    if (iter == params->end())
        throw OpenMMException("OpenMM_PropertyArray_get: No such property");
    return iter->second.c_str();
}"""

        for type in ('double', 'int'):
            name = 'OpenMM_%sArray' % type.capitalize()
            values = {'type':type, 'name':name}
            print >>self.out, """
/* %(name)s */
OPENMM_EXPORT %(name)s* %(name)s_create(int size) {
    return reinterpret_cast<%(name)s*>(new vector<%(type)s>(size));
}
OPENMM_EXPORT void %(name)s_destroy(%(name)s* array) {
    delete reinterpret_cast<vector<%(type)s>*>(array);
}
OPENMM_EXPORT int %(name)s_getSize(const %(name)s* array) {
    return reinterpret_cast<const vector<%(type)s>*>(array)->size();
}
OPENMM_EXPORT void %(name)s_resize(%(name)s* array, int size) {
    reinterpret_cast<vector<%(type)s>*>(array)->resize(size);
}
OPENMM_EXPORT void %(name)s_append(%(name)s* array, %(type)s value) {
    reinterpret_cast<vector<%(type)s>*>(array)->push_back(value);
}
OPENMM_EXPORT void %(name)s_set(%(name)s* array, int index, %(type)s value) {
    (*reinterpret_cast<vector<%(type)s>*>(array))[index] = value;
}
OPENMM_EXPORT %(type)s %(name)s_get(const %(name)s* array, int index) {
    return (*reinterpret_cast<const vector<%(type)s>*>(array))[index];
}""" % values

        for type in ('int',):
            name = 'OpenMM_%sSet' % type.capitalize()
            values = {'type':type, 'name':name}
            print >>self.out, """
/* %(name)s */
OPENMM_EXPORT %(name)s* %(name)s_create() {
    return reinterpret_cast<%(name)s*>(new set<%(type)s>());
}
OPENMM_EXPORT void %(name)s_destroy(%(name)s* s) {
    delete reinterpret_cast<set<%(type)s>*>(s);
}
OPENMM_EXPORT int %(name)s_getSize(const %(name)s* s) {
    return reinterpret_cast<const set<%(type)s>*>(s)->size();
}
OPENMM_EXPORT void %(name)s_insert(%(name)s* s, %(type)s value) {
    reinterpret_cast<set<%(type)s>*>(s)->insert(value);
}""" % values

        self.writeClasses()
        print >>self.out, "}\n"

#inputDirname = '/Users/peastman/workspace/openmm/bin-release/wrappers/doxygen/xml'
inputDirname = sys.argv[1]
builder = CHeaderGenerator(inputDirname, open(os.path.join(sys.argv[2], 'OpenMMCWrapper.h'), 'w'))
builder.writeOutput()
builder = CSourceGenerator(inputDirname, open(os.path.join(sys.argv[2], 'OpenMMCWrapper.cpp'), 'w'))
builder.writeOutput()
