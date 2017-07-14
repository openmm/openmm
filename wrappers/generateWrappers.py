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

OPENMM_RE_PATTERN=re.compile("(.*)OpenMM:[a-zA-Z0-9_:]*:(.*)")
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
        self.skipMethods = ['State OpenMM::Context::getState',
                            'void OpenMM::Context::createCheckpoint',
                            'void OpenMM::Context::loadCheckpoint',
                            'const std::vector<std::vector<int> >& OpenMM::Context::getMolecules',
                            'static std::vector<std::string> OpenMM::Platform::getPluginLoadFailures',
                            'static std::vector<std::string> OpenMM::Platform::loadPluginsFromDirectory',
                            'Vec3 OpenMM::LocalCoordinatesSite::getOriginWeights',
                            'Vec3 OpenMM::LocalCoordinatesSite::getXWeights',
                            'Vec3 OpenMM::LocalCoordinatesSite::getYWeights']
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
                if methodDefinition in self.skipMethods:
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
        isAbstract = any('virt' in method.attrib and method.attrib['virt'] == 'pure-virtual' for method in classNode.getiterator('memberdef'))

        if not isAbstract:
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
        print("""
#ifndef OPENMM_CWRAPPER_H_
#define OPENMM_CWRAPPER_H_

#ifndef OPENMM_EXPORT
#define OPENMM_EXPORT
#endif
""", file=self.out)
        self.writeGlobalConstants()
        self.writeTypeDeclarations()
        print("""
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
extern OPENMM_EXPORT const char* OpenMM_PropertyArray_get(const OpenMM_PropertyArray* array, const char* name);""", file=self.out)

        for type in ('double', 'int'):
            name = 'OpenMM_%sArray' % type.capitalize()
            values = {'type':type, 'name':name}
            print("""
/* %(name)s */
extern OPENMM_EXPORT %(name)s* %(name)s_create(int size);
extern OPENMM_EXPORT void %(name)s_destroy(%(name)s* array);
extern OPENMM_EXPORT int %(name)s_getSize(const %(name)s* array);
extern OPENMM_EXPORT void %(name)s_resize(%(name)s* array, int size);
extern OPENMM_EXPORT void %(name)s_append(%(name)s* array, %(type)s value);
extern OPENMM_EXPORT void %(name)s_set(%(name)s* array, int index, %(type)s value);
extern OPENMM_EXPORT %(type)s %(name)s_get(const %(name)s* array, int index);""" % values, file=self.out)

        for type in ('int',):
            name = 'OpenMM_%sSet' % type.capitalize()
            values = {'type':type, 'name':name}
            print("""
/* %(name)s */
extern OPENMM_EXPORT %(name)s* %(name)s_create();
extern OPENMM_EXPORT void %(name)s_destroy(%(name)s* set);
extern OPENMM_EXPORT int %(name)s_getSize(const %(name)s* set);
extern OPENMM_EXPORT void %(name)s_insert(%(name)s* set, %(type)s value);""" % values, file=self.out)

        print("""
/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
extern OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types, int enforcePeriodicBox);
extern OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState_2(const OpenMM_Context* target, int types, int enforcePeriodicBox, int groups);
extern OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory);
extern OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_getPluginLoadFailures();
extern OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeSystem(const OpenMM_System* system);
extern OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeState(const OpenMM_State* state);
extern OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeIntegrator(const OpenMM_Integrator* integrator);
extern OPENMM_EXPORT OpenMM_System* OpenMM_XmlSerializer_deserializeSystem(const char* xml);
extern OPENMM_EXPORT OpenMM_State* OpenMM_XmlSerializer_deserializeState(const char* xml);
extern OPENMM_EXPORT OpenMM_Integrator* OpenMM_XmlSerializer_deserializeIntegrator(const char* xml);""", file=self.out)

        self.writeClasses()

        print("""
#if defined(__cplusplus)
}
#endif

#endif /*OPENMM_CWRAPPER_H_*/""", file=self.out)


class CSourceGenerator(WrapperGenerator):
    """This class generates the source file for the C API wrappers."""

    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
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
        isAbstract = any('virt' in method.attrib and method.attrib['virt'] == 'pure-virtual' for method in classNode.getiterator('memberdef'))

        if not isAbstract:
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
#include "OpenMMCWrapper.h"
#include <cstdlib>
#include <cstring>
#include <sstream>
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
}""", file=self.out)

        for type in ('double', 'int'):
            name = 'OpenMM_%sArray' % type.capitalize()
            values = {'type':type, 'name':name}
            print("""
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
}""" % values, file=self.out)

        for type in ('int',):
            name = 'OpenMM_%sSet' % type.capitalize()
            values = {'type':type, 'name':name}
            print("""
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
}""" % values, file=self.out)

        print("""
/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState(const OpenMM_Context* target, int types, int enforcePeriodicBox) {
    State result = reinterpret_cast<const Context*>(target)->getState(types, enforcePeriodicBox);
    return reinterpret_cast<OpenMM_State*>(new State(result));
}
OPENMM_EXPORT OpenMM_State* OpenMM_Context_getState_2(const OpenMM_Context* target, int types, int enforcePeriodicBox, int groups) {
    State result = reinterpret_cast<const Context*>(target)->getState(types, enforcePeriodicBox, groups);
    return reinterpret_cast<OpenMM_State*>(new State(result));
}
OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_loadPluginsFromDirectory(const char* directory) {
    vector<string> result = Platform::loadPluginsFromDirectory(string(directory));
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(result));
}
OPENMM_EXPORT OpenMM_StringArray* OpenMM_Platform_getPluginLoadFailures() {
    vector<string> result = Platform::getPluginLoadFailures();
    return reinterpret_cast<OpenMM_StringArray*>(new vector<string>(result));
}
static char* createStringFromStream(stringstream& stream) {
    int length = stream.str().size();
    char* result = (char*) malloc(length+1);
    stream.str().copy(result, length);
    result[length] = 0;
    return result;
}
OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeSystem(const OpenMM_System* system) {
    stringstream stream;
    OpenMM::XmlSerializer::serialize<OpenMM::System>(reinterpret_cast<const OpenMM::System*>(system), "System", stream);
    return createStringFromStream(stream);
}
OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeState(const OpenMM_State* state) {
    stringstream stream;
    OpenMM::XmlSerializer::serialize<OpenMM::State>(reinterpret_cast<const OpenMM::State*>(state), "State", stream);
    return createStringFromStream(stream);
}
OPENMM_EXPORT char* OpenMM_XmlSerializer_serializeIntegrator(const OpenMM_Integrator* integrator) {
    stringstream stream;
    OpenMM::XmlSerializer::serialize<OpenMM::Integrator>(reinterpret_cast<const OpenMM::Integrator*>(integrator), "Integrator", stream);
    return createStringFromStream(stream);
}
OPENMM_EXPORT OpenMM_System* OpenMM_XmlSerializer_deserializeSystem(const char* xml) {
    string input(xml);
    stringstream stream(input);
    return reinterpret_cast<OpenMM_System*>(OpenMM::XmlSerializer::deserialize<OpenMM::System>(stream));
}
OPENMM_EXPORT OpenMM_State* OpenMM_XmlSerializer_deserializeState(const char* xml) {
    string input(xml);
    stringstream stream(input);
    return reinterpret_cast<OpenMM_State*>(OpenMM::XmlSerializer::deserialize<OpenMM::State>(stream));
}
OPENMM_EXPORT OpenMM_Integrator* OpenMM_XmlSerializer_deserializeIntegrator(const char* xml) {
    string input(xml);
    stringstream stream(input);
    return reinterpret_cast<OpenMM_Integrator*>(OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(stream));
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
                self.out.write("    real*8, parameter :: OpenMM_%s = %s\n" % (vDef, iDef))

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
        isAbstract = any('virt' in method.attrib and method.attrib['virt'] == 'pure-virtual' for method in classNode.getiterator('memberdef'))

        if not isAbstract:
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
            functionName = functionName[:63]
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
        if type in self.enumerationTypes:
            return 'integer*4'
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
    type OpenMM_Vec3Array
        integer*8 :: handle = 0
    end type

    type OpenMM_StringArray
        integer*8 :: handle = 0
    end type

    type OpenMM_BondArray
        integer*8 :: handle = 0
    end type

    type OpenMM_ParameterArray
        integer*8 :: handle = 0
    end type

    type OpenMM_PropertyArray
        integer*8 :: handle = 0
    end type

    type OpenMM_DoubleArray
        integer*8 :: handle = 0
    end type

    type OpenMM_IntArray
        integer*8 :: handle = 0
    end type

    type OpenMM_IntSet
        integer*8 :: handle = 0
    end type

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

        ! OpenMM_Vec3
        subroutine OpenMM_Vec3_scale(vec, scale, result)
            use OpenMM_Types; implicit none
            real*8 vec(3)
            real*8 scale
            real*8 result(3)
        end subroutine

        ! OpenMM_Vec3Array
        subroutine OpenMM_Vec3Array_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_Vec3Array) result
        end subroutine
        subroutine OpenMM_Vec3Array_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) destroy
        end subroutine
        function OpenMM_Vec3Array_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 OpenMM_Vec3Array_getSize
        end function
        subroutine OpenMM_Vec3Array_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 size
        end subroutine
        subroutine OpenMM_Vec3Array_append(target, vec)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            real*8 vec(3)
        end subroutine
        subroutine OpenMM_Vec3Array_set(target, index, vec)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 index
            real*8 vec(3)
        end subroutine
        subroutine OpenMM_Vec3Array_get(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Vec3Array) target
            integer*4 index
            real*8 result(3)
        end subroutine

        ! OpenMM_StringArray
        subroutine OpenMM_StringArray_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_StringArray) result
        end subroutine
        subroutine OpenMM_StringArray_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) destroy
        end subroutine
        function OpenMM_StringArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 OpenMM_StringArray_getSize
        end function
        subroutine OpenMM_StringArray_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 size
        end subroutine
        subroutine OpenMM_StringArray_append(target, str)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            character(*) str
        end subroutine
        subroutine OpenMM_StringArray_set(target, index, str)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 index
            character(*) str
        end subroutine
        subroutine OpenMM_StringArray_get(target, index, result)
            use OpenMM_Types; implicit none
            type (OpenMM_StringArray) target
            integer*4 index
            character(*) result
        end subroutine

        ! OpenMM_BondArray
        subroutine OpenMM_BondArray_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (OpenMM_BondArray) result
        end subroutine
        subroutine OpenMM_BondArray_destroy(destroy)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) destroy
        end subroutine
        function OpenMM_BondArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 OpenMM_BondArray_getSize
        end function
        subroutine OpenMM_BondArray_resize(target, size)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 size
        end subroutine
        subroutine OpenMM_BondArray_append(target, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 particle1
            integer*4 particle2
        end subroutine
        subroutine OpenMM_BondArray_set(target, index, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
        end subroutine
        subroutine OpenMM_BondArray_get(target, index, particle1, particle2)
            use OpenMM_Types; implicit none
            type (OpenMM_BondArray) target
            integer*4 index
            integer*4 particle1
            integer*4 particle2
        end subroutine

        ! OpenMM_ParameterArray
        function OpenMM_ParameterArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            integer*4 OpenMM_ParameterArray_getSize
        end function
        subroutine OpenMM_ParameterArray_get(target, name, result)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            character(*) name
            character(*) result
        end subroutine

        ! OpenMM_PropertyArray
        function OpenMM_PropertyArray_getSize(target)
            use OpenMM_Types; implicit none
            type (OpenMM_ParameterArray) target
            integer*4 OpenMM_PropertyArray_getSize
        end function
        subroutine OpenMM_PropertyArray_get(target, name, result)
            use OpenMM_Types; implicit none
            type (OpenMM_PropertyArray) target
            character(*) name
            character(*) result
        end subroutine""", file=self.out)

        arrayTypes = {'OpenMM_DoubleArray':'real*8', 'OpenMM_IntArray':'integer*4'}
        for name in arrayTypes:
            values = {'type':arrayTypes[name], 'name':name}
            print("""
        ! %(name)s
        subroutine %(name)s_create(result, size)
            use OpenMM_Types; implicit none
            integer*4 size
            type (%(name)s) result
        end subroutine
        subroutine %(name)s_destroy(destroy)
            use OpenMM_Types; implicit none
            type (%(name)s) destroy
        end subroutine
        function %(name)s_getSize(target)
            use OpenMM_Types; implicit none
            type (%(name)s) target
            integer*4 %(name)s_getSize
        end function
        subroutine %(name)s_resize(target, size)
            use OpenMM_Types; implicit none
            type (%(name)s) target
            integer*4 size
        end subroutine
        subroutine %(name)s_append(target, value)
            use OpenMM_Types; implicit none
            type (%(name)s) target
            %(type)s value
        end subroutine
        subroutine %(name)s_set(target, index, value)
            use OpenMM_Types; implicit none
            type (%(name)s) target
            integer*4 index
            %(type)s value
        end subroutine
        subroutine %(name)s_get(target, index, result)
            use OpenMM_Types; implicit none
            type (%(name)s) target
            integer*4 index
            %(type)s result
        end subroutine""" % values, file=self.out)
        
        print("""
        ! These methods need to be handled specially, since their C++ APIs cannot be directly translated to Fortran.
        ! Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself.
        subroutine OpenMM_Context_getState(target, types, enforcePeriodicBox, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            integer*4 types
            integer*4 enforcePeriodicBox
            type(OpenMM_State) result
        end subroutine
        subroutine OpenMM_Context_getState_2(target, types, enforcePeriodicBox, groups, result)
            use OpenMM_Types; implicit none
            type (OpenMM_Context) target
            integer*4 types
            integer*4 enforcePeriodicBox
            integer*4 groups
            type(OpenMM_State) result
        end subroutine
        subroutine OpenMM_Platform_loadPluginsFromDirectory(directory, result)
            use OpenMM_Types; implicit none
            character(*) directory
            type(OpenMM_StringArray) result
        end subroutine
        subroutine OpenMM_Platform_getPluginLoadFailures(result)
            use OpenMM_Types; implicit none
            type(OpenMM_StringArray) result
        end subroutine
        subroutine OpenMM_XmlSerializer_serializeSystemToC(system, result, result_length)
            use iso_c_binding; use OpenMM_Types; implicit none
            type(OpenMM_System), intent(in) :: system
            type(c_ptr), intent(out) :: result
            integer, intent(out) :: result_length
        end subroutine
        subroutine OpenMM_XmlSerializer_serializeStateToC(state, result, result_length)
            use iso_c_binding; use OpenMM_Types; implicit none
            type(OpenMM_State), intent(in) :: state
            type(c_ptr), intent(out) :: result
            integer, intent(out) :: result_length
        end subroutine
        subroutine OpenMM_XmlSerializer_serializeIntegratorToC(integrator, result, result_length)
            use iso_c_binding; use OpenMM_Types; implicit none
            type(OpenMM_Integrator), intent(in) :: integrator
            type(c_ptr), intent(out) :: result
            integer, intent(out) :: result_length
        end subroutine
        subroutine OpenMM_XmlSerializer_deserializeSystem(xml, result)
            use OpenMM_Types; implicit none
            character(*) xml
            type(OpenMM_System) result
        end subroutine
        subroutine OpenMM_XmlSerializer_deserializeState(xml, result)
            use OpenMM_Types; implicit none
            character(*) xml
            type(OpenMM_State) result
        end subroutine
        subroutine OpenMM_XmlSerializer_deserializeIntegrator(xml, result)
            use OpenMM_Types; implicit none
            character(*) xml
            type(OpenMM_Integrator) result
        end subroutine""", file=self.out)
        
        self.writeClasses()
        
        print("""
    end interface

contains

    subroutine OpenMM_XmlSerializer_serializeSystem(system, result)
        use iso_c_binding, only: c_ptr, c_int, c_char, c_f_pointer
        type(OpenMM_System), intent(in) :: system
        character(len=1), allocatable, dimension(:), intent(out) :: result

        character(kind=c_char), pointer, dimension(:) :: fstr
        type(c_ptr) :: cstr
        integer :: i
        integer(kind=c_int) :: result_length

        call OpenMM_XmlSerializer_serializeSystemToC(system, cstr, result_length)
        call c_f_pointer(cstr, fstr, [ result_length ])
        allocate(character(len=1) :: result(result_length))
        do i=1,result_length
           result(i) = fstr(i)
        end do
    end subroutine

    subroutine OpenMM_XmlSerializer_serializeState(state, result)
        use iso_c_binding, only: c_ptr, c_int, c_char, c_f_pointer
        type(OpenMM_State), intent(in) :: state
        character(len=1), allocatable, dimension(:), intent(out) :: result

        character(kind=c_char), pointer, dimension(:) :: fstr
        type(c_ptr) :: cstr
        integer :: i
        integer(kind=c_int) :: result_length

        call OpenMM_XmlSerializer_serializeStateToC(state, cstr, result_length)
        call c_f_pointer(cstr, fstr, [ result_length ])
        allocate(character(len=1) :: result(result_length))
        do i=1,result_length
           result(i) = fstr(i)
        end do
    end subroutine

    subroutine OpenMM_XmlSerializer_serializeIntegrator(integrator, result)
        use iso_c_binding, only: c_ptr, c_int, c_char, c_f_pointer
        type(OpenMM_Integrator), intent(in) :: integrator
        character(len=1), allocatable, dimension(:), intent(out) :: result

        character(kind=c_char), pointer, dimension(:) :: fstr
        type(c_ptr) :: cstr
        integer :: i
        integer(kind=c_int) :: result_length

        call OpenMM_XmlSerializer_serializeIntegratorToC(integrator, cstr, result_length)
        call c_f_pointer(cstr, fstr, [ result_length ])
        allocate(character(len=1) :: result(result_length))
        do i=1,result_length
           result(i) = fstr(i)
        end do
    end subroutine

END MODULE OpenMM""", file=self.out)


class FortranSourceGenerator(WrapperGenerator):
    """This class generates the source file for the Fortran API wrappers."""

    def __init__(self, inputDirname, output):
        WrapperGenerator.__init__(self, inputDirname, output)
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
        isAbstract = any('virt' in method.attrib and method.attrib['virt'] == 'pure-virtual' for method in classNode.getiterator('memberdef'))

        if not isAbstract:
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
            truncatedName = functionName[:63]
            self.writeOneMethod(classNode, methodNode, functionName, truncatedName.lower()+'_')
            self.writeOneMethod(classNode, methodNode, functionName, truncatedName.upper())
    
    def writeOneConstructor(self, classNode, methodNode, functionName, wrapperFunctionName):
        className = getText("compoundname", classNode)
        shortClassName = stripOpenMMPrefix(className)
        typeName = convertOpenMMPrefix(className)
        self.out.write("OPENMM_EXPORT void %s(%s*& result" % (wrapperFunctionName, typeName))
        self.writeArguments(methodNode, True)
        self.out.write(") {\n")
        self.out.write("    result = %s(" % functionName)
        self.writeInvocationArguments(methodNode, False)
        self.out.write(");\n")
        self.out.write("}\n")
    
    def writeOneDestructor(self, typeName, wrapperFunctionName):
        self.out.write("OPENMM_EXPORT void %s(%s*& destroy) {\n" % (wrapperFunctionName, typeName))
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
        self.out.write("OPENMM_EXPORT ")
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
#include "OpenMMCWrapper.h"
#include <cstring>
#include <vector>
#include <cstdlib>

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

static void convertStringToChars(char* source, char*& cstr, int& length) {
	length = strlen(source);
	cstr = new char[length+1];
	strcpy(cstr, source);
    free(source);
}

extern "C" {

/* OpenMM_Vec3 */
OPENMM_EXPORT void openmm_vec3_scale_(const OpenMM_Vec3& vec, double const& scale, OpenMM_Vec3& result) {
    result = OpenMM_Vec3_scale(vec, scale);
}
OPENMM_EXPORT void OPENMM_VEC3_SCALE(const OpenMM_Vec3& vec, double const& scale, OpenMM_Vec3& result) {
    result = OpenMM_Vec3_scale(vec, scale);
}

/* OpenMM_Vec3Array */
OPENMM_EXPORT void openmm_vec3array_create_(OpenMM_Vec3Array*& result, const int& size) {
    result = OpenMM_Vec3Array_create(size);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_CREATE(OpenMM_Vec3Array*& result, const int& size) {
    result = OpenMM_Vec3Array_create(size);
}
OPENMM_EXPORT void openmm_vec3array_destroy_(OpenMM_Vec3Array*& array) {
    OpenMM_Vec3Array_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_DESTROY(OpenMM_Vec3Array*& array) {
    OpenMM_Vec3Array_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_vec3array_getsize_(const OpenMM_Vec3Array* const& array) {
    return OpenMM_Vec3Array_getSize(array);
}
OPENMM_EXPORT int OPENMM_VEC3ARRAY_GETSIZE(const OpenMM_Vec3Array* const& array) {
    return OpenMM_Vec3Array_getSize(array);
}
OPENMM_EXPORT void openmm_vec3array_resize_(OpenMM_Vec3Array* const& array, const int& size) {
    OpenMM_Vec3Array_resize(array, size);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_RESIZE(OpenMM_Vec3Array* const& array, const int& size) {
    OpenMM_Vec3Array_resize(array, size);
}
OPENMM_EXPORT void openmm_vec3array_append_(OpenMM_Vec3Array* const& array, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_append(array, vec);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_APPEND(OpenMM_Vec3Array* const& array, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_append(array, vec);
}
OPENMM_EXPORT void openmm_vec3array_set_(OpenMM_Vec3Array* const& array, const int& index, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_set(array, index-1, vec);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_SET(OpenMM_Vec3Array* const& array, const int& index, const OpenMM_Vec3& vec) {
    OpenMM_Vec3Array_set(array, index-1, vec);
}
OPENMM_EXPORT void openmm_vec3array_get_(const OpenMM_Vec3Array* const& array, const int& index, OpenMM_Vec3& result) {
    result = *OpenMM_Vec3Array_get(array, index-1);
}
OPENMM_EXPORT void OPENMM_VEC3ARRAY_GET(const OpenMM_Vec3Array* const& array, const int& index, OpenMM_Vec3& result) {
    result = *OpenMM_Vec3Array_get(array, index-1);
}

/* OpenMM_StringArray */
OPENMM_EXPORT void openmm_stringarray_create_(OpenMM_StringArray*& result, const int& size) {
    result = OpenMM_StringArray_create(size);
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_CREATE(OpenMM_StringArray*& result, const int& size) {
    result = OpenMM_StringArray_create(size);
}
OPENMM_EXPORT void openmm_stringarray_destroy_(OpenMM_StringArray*& array) {
    OpenMM_StringArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_DESTROY(OpenMM_StringArray*& array) {
    OpenMM_StringArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_stringarray_getsize_(const OpenMM_StringArray* const& array) {
    return OpenMM_StringArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_STRINGARRAY_GETSIZE(const OpenMM_StringArray* const& array) {
    return OpenMM_StringArray_getSize(array);
}
OPENMM_EXPORT void openmm_stringarray_resize_(OpenMM_StringArray* const& array, const int& size) {
    OpenMM_StringArray_resize(array, size);
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_RESIZE(OpenMM_StringArray* const& array, const int& size) {
    OpenMM_StringArray_resize(array, size);
}
OPENMM_EXPORT void openmm_stringarray_append_(OpenMM_StringArray* const& array, const char* str, int length) {
    OpenMM_StringArray_append(array, makeString(str, length).c_str());
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_APPEND(OpenMM_StringArray* const& array, const char* str, int length) {
    OpenMM_StringArray_append(array, makeString(str, length).c_str());
}
OPENMM_EXPORT void openmm_stringarray_set_(OpenMM_StringArray* const& array, const int& index, const char* str, int length) {
  OpenMM_StringArray_set(array, index-1, makeString(str, length).c_str());
  }
OPENMM_EXPORT void OPENMM_STRINGARRAY_SET(OpenMM_StringArray* const& array, const int& index, const char* str, int length) {
  OpenMM_StringArray_set(array, index-1, makeString(str, length).c_str());
}
OPENMM_EXPORT void openmm_stringarray_get_(const OpenMM_StringArray* const& array, const int& index, char* result, int length) {
    const char* str = OpenMM_StringArray_get(array, index-1);
    copyAndPadString(result, str, length);
}
OPENMM_EXPORT void OPENMM_STRINGARRAY_GET(const OpenMM_StringArray* const& array, const int& index, char* result, int length) {
    const char* str = OpenMM_StringArray_get(array, index-1);
    copyAndPadString(result, str, length);
}

/* OpenMM_BondArray */
OPENMM_EXPORT void openmm_bondarray_create_(OpenMM_BondArray*& result, const int& size) {
    result = OpenMM_BondArray_create(size);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_CREATE(OpenMM_BondArray*& result, const int& size) {
    result = OpenMM_BondArray_create(size);
}
OPENMM_EXPORT void openmm_bondarray_destroy_(OpenMM_BondArray*& array) {
    OpenMM_BondArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT void OPENMM_BONDARRAY_DESTROY(OpenMM_BondArray*& array) {
    OpenMM_BondArray_destroy(array);
    array = 0;
}
OPENMM_EXPORT int openmm_bondarray_getsize_(const OpenMM_BondArray* const& array) {
    return OpenMM_BondArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_BONDARRAY_GETSIZE(const OpenMM_BondArray* const& array) {
    return OpenMM_BondArray_getSize(array);
}
OPENMM_EXPORT void openmm_bondarray_resize_(OpenMM_BondArray* const& array, const int& size) {
    OpenMM_BondArray_resize(array, size);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_RESIZE(OpenMM_BondArray* const& array, const int& size) {
    OpenMM_BondArray_resize(array, size);
}
OPENMM_EXPORT void openmm_bondarray_append_(OpenMM_BondArray* const& array, const int& particle1, const int& particle2) {
    OpenMM_BondArray_append(array, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_APPEND(OpenMM_BondArray* const& array, const int& particle1, const int& particle2) {
    OpenMM_BondArray_append(array, particle1, particle2);
}
OPENMM_EXPORT void openmm_bondarray_set_(OpenMM_BondArray* const& array, const int& index, const int& particle1, const int& particle2) {
    OpenMM_BondArray_set(array, index-1, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_SET(OpenMM_BondArray* const& array, const int& index, const int& particle1, const int& particle2) {
    OpenMM_BondArray_set(array, index-1, particle1, particle2);
}
OPENMM_EXPORT void openmm_bondarray_get_(const OpenMM_BondArray* const& array, const int& index, int* particle1, int* particle2) {
    OpenMM_BondArray_get(array, index-1, particle1, particle2);
}
OPENMM_EXPORT void OPENMM_BONDARRAY_GET(const OpenMM_BondArray* const& array, const int& index, int* particle1, int* particle2) {
    OpenMM_BondArray_get(array, index-1, particle1, particle2);
}

/* OpenMM_ParameterArray */
OPENMM_EXPORT int openmm_parameterarray_getsize_(const OpenMM_ParameterArray* const& array) {
    return OpenMM_ParameterArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_PARAMETERARRAY_GETSIZE(const OpenMM_ParameterArray* const& array) {
    return OpenMM_ParameterArray_getSize(array);
}
OPENMM_EXPORT double openmm_parameterarray_get_(const OpenMM_ParameterArray* const& array, const char* name, int length) {
    return OpenMM_ParameterArray_get(array, makeString(name, length).c_str());
}
OPENMM_EXPORT double OPENMM_PARAMETERARRAY_GET(const OpenMM_ParameterArray* const& array, const char* name, int length) {
    return OpenMM_ParameterArray_get(array, makeString(name, length).c_str());
}

/* OpenMM_PropertyArray */
OPENMM_EXPORT int openmm_propertyarray_getsize_(const OpenMM_PropertyArray* const& array) {
    return OpenMM_PropertyArray_getSize(array);
}
OPENMM_EXPORT int OPENMM_PROPERTYARRAY_GETSIZE(const OpenMM_PropertyArray* const& array) {
    return OpenMM_PropertyArray_getSize(array);
}
OPENMM_EXPORT const char* openmm_propertyarray_get_(const OpenMM_PropertyArray* const& array, const char* name, int length) {
    return OpenMM_PropertyArray_get(array, makeString(name, length).c_str());
}
OPENMM_EXPORT const char* OPENMM_PROPERTYARRAY_GET(const OpenMM_PropertyArray* const& array, const char* name, int length) {
    return OpenMM_PropertyArray_get(array, makeString(name, length).c_str());
}""", file=self.out)

        for type in ('double', 'int'):
            name = 'OpenMM_%sArray' % type.capitalize()
            values = {'type':type, 'name':name, 'name_lower':name.lower(), 'name_upper':name.upper()}
            print("""
/* %(name)s */
OPENMM_EXPORT void %(name_lower)s_create_(%(name)s*& result, const int& size) {
    result = %(name)s_create(size);
}
OPENMM_EXPORT void %(name_upper)s_CREATE(%(name)s*& result, const int& size) {
    result = %(name)s_create(size);
}
OPENMM_EXPORT void %(name_lower)s_destroy_(%(name)s*& array) {
    %(name)s_destroy(array);
    array = 0;
}
OPENMM_EXPORT void %(name_upper)s_DESTROY(%(name)s*& array) {
    %(name)s_destroy(array);
    array = 0;
}
OPENMM_EXPORT int %(name_lower)s_getsize_(const %(name)s* const& array) {
    return %(name)s_getSize(array);
}
OPENMM_EXPORT int %(name_upper)s_GETSIZE(const %(name)s* const& array) {
    return %(name)s_getSize(array);
}
OPENMM_EXPORT void %(name_lower)s_resize_(%(name)s* const& array, const int& size) {
    %(name)s_resize(array, size);
}
OPENMM_EXPORT void %(name_upper)s_RESIZE(%(name)s* const& array, const int& size) {
    %(name)s_resize(array, size);
}
OPENMM_EXPORT void %(name_lower)s_append_(%(name)s* const& array, const %(type)s& value) {
    %(name)s_append(array, value);
}
OPENMM_EXPORT void %(name_upper)s_APPEND(%(name)s* const& array, const %(type)s& value) {
    %(name)s_append(array, value);
}
OPENMM_EXPORT void %(name_lower)s_set_(%(name)s* const& array, const int& index, const %(type)s& value) {
    %(name)s_set(array, index-1, value);
}
OPENMM_EXPORT void %(name_upper)s_SET(%(name)s* const& array, const int& index, const %(type)s& value) {
    %(name)s_set(array, index-1, value);
}
OPENMM_EXPORT void %(name_lower)s_get_(const %(name)s* const& array, const int& index, %(type)s& result) {
    result = %(name)s_get(array, index-1);
}
OPENMM_EXPORT void %(name_upper)s_GET(const %(name)s* const& array, const int& index, %(type)s& result) {
    result = %(name)s_get(array, index-1);
}""" % values, file=self.out)

        for type in ('int', ):
            name = 'OpenMM_%sSet' % type.capitalize()
            values = {'type':type, 'name':name, 'name_lower':name.lower(), 'name_upper':name.upper()}
            print("""
/* %(name)s */
OPENMM_EXPORT void %(name_lower)s_create_(%(name)s*& result) {
    result = %(name)s_create();
}
OPENMM_EXPORT void %(name_upper)s_CREATE(%(name)s*& result) {
    result = %(name)s_create();
}
OPENMM_EXPORT void %(name_lower)s_destroy_(%(name)s*& array) {
    %(name)s_destroy(array);
    array = 0;
}
OPENMM_EXPORT void %(name_upper)s_DESTROY(%(name)s*& array) {
    %(name)s_destroy(array);
    array = 0;
}
OPENMM_EXPORT int %(name_lower)s_getsize_(const %(name)s* const& array) {
    return %(name)s_getSize(array);
}
OPENMM_EXPORT int %(name_upper)s_GETSIZE(const %(name)s* const& array) {
    return %(name)s_getSize(array);
}
OPENMM_EXPORT void %(name_lower)s_insert_(%(name)s* const& array, const %(type)s& value) {
    %(name)s_insert(array, value);
}
OPENMM_EXPORT void %(name_upper)s_INSERT(%(name)s* const& array, const %(type)s& value) {
    %(name)s_insert(array, value);
}""" % values, file=self.out)

        print("""
/* These methods need to be handled specially, since their C++ APIs cannot be directly translated to C.
   Unlike the C++ versions, the return value is allocated on the heap, and you must delete it yourself. */
OPENMM_EXPORT void openmm_context_getstate_(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, OpenMM_State*& result) {
    result = OpenMM_Context_getState(target, types, enforcePeriodicBox);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETSTATE(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, OpenMM_State*& result) {
    result = OpenMM_Context_getState(target, types, enforcePeriodicBox);
}
OPENMM_EXPORT void openmm_context_getstate_2_(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, int const& groups, OpenMM_State*& result) {
    result = OpenMM_Context_getState_2(target, types, enforcePeriodicBox, groups);
}
OPENMM_EXPORT void OPENMM_CONTEXT_GETSTATE_2(const OpenMM_Context*& target, int const& types, int const& enforcePeriodicBox, int const& groups, OpenMM_State*& result) {
    result = OpenMM_Context_getState_2(target, types, enforcePeriodicBox, groups);
}
OPENMM_EXPORT void openmm_platform_loadpluginsfromdirectory_(const char* directory, OpenMM_StringArray*& result, int length) {
    result = OpenMM_Platform_loadPluginsFromDirectory(makeString(directory, length).c_str());
}
OPENMM_EXPORT void OPENMM_PLATFORM_LOADPLUGINSFROMDIRECTORY(const char* directory, OpenMM_StringArray*& result, int length) {
    result = OpenMM_Platform_loadPluginsFromDirectory(makeString(directory, length).c_str());
}
OPENMM_EXPORT void openmm_platform_getpluginloadfailures_(OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPluginLoadFailures();
}
OPENMM_EXPORT void OPENMM_PLATFORM_GETPLUGINLOADFAILURES(OpenMM_StringArray*& result) {
    result = OpenMM_Platform_getPluginLoadFailures();
}
OPENMM_EXPORT void openmm_xmlserializer_serializesystemtoc_(OpenMM_System*& system, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeSystem(system), result, result_length);
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_SERIALIZESYSTEMTOC(OpenMM_System*& system, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeSystem(system), result, result_length);
}
OPENMM_EXPORT void openmm_xmlserializer_serializestatetoc_(OpenMM_State*& state, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeState(state), result, result_length);
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_SERIALIZESTATETOC(OpenMM_State*& state, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeState(state), result, result_length);
}
OPENMM_EXPORT void openmm_xmlserializer_serializeintegratortoc_(OpenMM_Integrator*& integrator, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeIntegrator(integrator), result, result_length);
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_SERIALIZEINTEGRATORTOC(OpenMM_Integrator*& integrator, char*& result, int& result_length) {
    convertStringToChars(OpenMM_XmlSerializer_serializeIntegrator(integrator), result, result_length);
}
OPENMM_EXPORT void openmm_xmlserializer_deserializesystem_(const char* xml, OpenMM_System*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeSystem(makeString(xml, length).c_str());
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_DESERIALIZESYSTEM(const char* xml, OpenMM_System*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeSystem(makeString(xml, length).c_str());
}
OPENMM_EXPORT void openmm_xmlserializer_deserializestate_(const char* xml, OpenMM_State*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeState(makeString(xml, length).c_str());
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_DESERIALIZESTATE(const char* xml, OpenMM_State*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeState(makeString(xml, length).c_str());
}
OPENMM_EXPORT void openmm_xmlserializer_deserializeintegrator_(const char* xml, OpenMM_Integrator*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeIntegrator(makeString(xml, length).c_str());
}
OPENMM_EXPORT void OPENMM_XMLSERIALIZER_DESERIALIZEINTEGRATOR(const char* xml, OpenMM_Integrator*& result, int length) {
    result = OpenMM_XmlSerializer_deserializeIntegrator(makeString(xml, length).c_str());
}""", file=self.out)

        self.writeClasses()
        print("}", file=self.out)

inputDirname = sys.argv[1]
builder = CHeaderGenerator(inputDirname, open(os.path.join(sys.argv[2], 'OpenMMCWrapper.h'), 'w'))
builder.writeOutput()
builder = CSourceGenerator(inputDirname, open(os.path.join(sys.argv[2], 'OpenMMCWrapper.cpp'), 'w'))
builder.writeOutput()
builder = FortranHeaderGenerator(inputDirname, open(os.path.join(sys.argv[2], 'OpenMMFortranModule.f90'), 'w'))
builder.writeOutput()
builder = FortranSourceGenerator(inputDirname, open(os.path.join(sys.argv[2], 'OpenMMFortranWrapper.cpp'), 'w'))
builder.writeOutput()
