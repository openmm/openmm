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
        self.nodeByID={}

        # Read all the XML files and merge them into a single document.
        self.doc = etree.ElementTree(etree.Element('root'))
        for file in os.listdir(inputDirname):
            root = etree.parse(os.path.join(inputDirname, file)).getroot()
            for node in root:
                self.doc.getroot().append(node)

        self.fOut = output

        self._enumByClassname={}
        self.typesByShortName = {}
        self._orderedClassNodes = self._buildOrderedClassNodes()

    def _getNodeByID(self, id):
        if id not in self.nodeByID:
            for node in findNodes(self.doc.getroot(), "compounddef", id=id):
                self.nodeByID[id] = node
        return self.nodeByID[id]

    def _buildOrderedClassNodes(self):
        orderedClassNodes=[]
        for node in findNodes(self.doc.getroot(), "compounddef", kind="class", prot="public"):
            self._findBaseNodes(node, orderedClassNodes)
        return orderedClassNodes

    def _findBaseNodes(self, node, excludedClassNodes=[]):
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
                baseNode = self._getNodeByID(baseNodeID)
                self._findBaseNodes(baseNode, excludedClassNodes)
        excludedClassNodes.append(node)

    def writeGlobalConstants(self):
        self.fOut.write("/* Global Constants */\n\n")
        node = next((x for x in findNodes(self.doc.getroot(), "compounddef", kind="namespace") if x.findtext("compoundname") == "OpenMM"))
        for section in findNodes(node, "sectiondef", kind="var"):
            for memberNode in findNodes(section, "memberdef", kind="variable", mutable="no", prot="public", static="yes"):
                vDef = convertOpenMMPrefix(getText("definition", memberNode))
                iDef = getText("initializer", memberNode)
                if iDef.startswith("="):
                    iDef = iDef[1:]
                self.fOut.write("static %s = %s;\n" % (vDef, iDef))

    def writeTypeDeclarations(self):
        self.fOut.write("\n/* Type Declarations */\n\n")
        for classNode in self._orderedClassNodes:
            className = getText("compoundname", classNode)
            shortName = stripOpenMMPrefix(className)
            typeName = convertOpenMMPrefix(className)
            self.fOut.write("typedef struct %s_struct %s;\n" % (typeName, typeName))
            self.typesByShortName[shortName] = typeName

    def writeClasses(self):
        for classNode in self._orderedClassNodes:
            className = stripOpenMMPrefix(getText("compoundname", classNode))
            self.fOut.write("\n/* %s */\n" % className)
            self.writeEnumerations(classNode)
            self.writeMethods(classNode)
        self.fOut.write("\n")

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
            self.fOut.write("typedef enum {\n  ")
            argSep=""
            for valueNode in findNodes(enumNode, "enumvalue", prot="public"):
                vName = convertOpenMMPrefix(getText("name", valueNode))
                vInit = getText("initializer", valueNode)
                if vInit.startswith("="):
                    vInit = vInit[1:].strip()
                self.fOut.write("%s%s_%s = %s" % (argSep, typeName, vName, vInit))
                argSep=", "
            self.fOut.write("\n} %s;\n" % enumTypeName)
            self.typesByShortName[enumName] = enumTypeName
        if len(enumNodes)>0: self.fOut.write("\n")

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
                    self.fOut.write("extern OPENMM_EXPORT %s* %s_create%s(" % (typeName, typeName, suffix))
                    self.writeArguments(methodNode, False)
                    self.fOut.write(");\n")
    
            # Write destructor
            self.fOut.write("extern OPENMM_EXPORT void %s_destroy(%s* target);\n" % (typeName, typeName))

        # Write other methods
        for methodNode in methodList:
            methodDefinition = getText("definition", methodNode)
            shortMethodDefinition = stripOpenMMPrefix(methodDefinition)
            methodName = shortMethodDefinition.split()[-1]
            if methodName in (shortClassName, destructorName):
                continue
            if self.shouldHideMethod(methodNode):
                continue
            returnType = self.getType(getText("type", methodNode))
            self.fOut.write("extern OPENMM_EXPORT %s %s_%s(" % (returnType, typeName, methodName))
            isInstanceMethod = (methodNode.attrib['static'] != 'yes')
            if isInstanceMethod:
                if methodNode.attrib['const'] == 'yes':
                    self.fOut.write('const ')
                self.fOut.write("%s* target" % typeName)
            self.writeArguments(methodNode, isInstanceMethod)
            self.fOut.write(");\n")
    
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
            self.fOut.write("%s%s %s" % (separator, type, name))
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


inputDirname = '/Users/peastman/workspace/openmm/bin-release/wrappers/python/src/swig_doxygen/doxygen/xml'
out = sys.stdout
builder = WrapperGenerator(inputDirname, out)
print >>out, """
#ifndef OPENMM_CWRAPPER_H_
#define OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT
#define OPENMM_EXPORT
#endif
"""
builder.writeGlobalConstants()
builder.writeTypeDeclarations()
print >>out, """
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
    print >>out, """
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
    print >>out, """
/* %(name)s */
extern OPENMM_EXPORT %(name)s* %(name)s_create();
extern OPENMM_EXPORT void %(name)s_destroy(%(name)s* set);
extern OPENMM_EXPORT int %(name)s_getSize(const %(name)s* set);
extern OPENMM_EXPORT void %(name)s_insert(%(name)s* set, %(type)s value);""" % values

print >>out, """
/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
extern OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types, int enforcePeriodicBox);
extern OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory);"""

builder.writeClasses()

print >>out, """
#if defined(__cplusplus)
}
#endif

#endif /*OPENMM_CWRAPPER_H_*/"""