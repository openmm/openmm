#!/usr/bin/env python
"""Build swig imput file from xml encoded header files (see gccxml)."""
__author__ = "Randall J. Radmer"
__version__ = "1.0"


import sys, os
import time
import getopt
import re
import xml.etree.ElementTree as etree
from distutils.version import LooseVersion
import copy

try:
    from html.parser import HTMLParser
except ImportError:
    # python 2
    from HTMLParser import HTMLParser

INDENT = "   "
docTags = {'emphasis':'i', 'bold':'b', 'itemizedlist':'ul', 'listitem':'li', 'preformatted':'pre', 'computeroutput':'tt', 'subscript':'sub'}

def striphtmltags(s):
    """Strip a couple html tags used inside docstrings in the C++ source
    to produce something more easily read as plain text.
    """
    class ConvertLists(HTMLParser):
        def reset(self):
            HTMLParser.reset(self)
            self.out = []

        def handle_starttag(self, tag, attrs):
            if tag == 'li':
                self.out.append('\n - ')
        def handle_data(self, data):
            self.out.append(data.strip())

    convertlists = ConvertLists()

    def replace_ul_tags(m):
        a, b = m.span()
        sub = s[a:b]

        convertlists.reset()
        convertlists.feed(sub)
        return '\n%s\n\n' % ''.join(convertlists.out)

    s = s.replace('<i>', '_').replace('</i>', '_')
    s = s.replace('<b>', '*').replace('</b>', '*')

    s = re.sub('\s*(<ul>.*?</ul>\s*)', replace_ul_tags, s, flags=re.MULTILINE | re.DOTALL)
    return s

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
        elif n.tag == "xrefsect":
            title = n.find("xreftitle")
            description = n.find("xrefdescription")
            if title is not None and description is not None and getNodeText(title).lower() == "deprecated":
                s = "%s\n@deprecated %s\n\n" % (s, getNodeText(description))
        else:
            if n.tag in docTags:
                tag = docTags[n.tag]
                s = "%s<%s>%s</%s>" % (s, tag, getNodeText(n), tag)
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

OPENMM_RE_PATTERN=re.compile("(.*)OpenMM:[a-zA-Z:]*:(.*)")
def stripOpenmmPrefix(name, rePattern=OPENMM_RE_PATTERN):
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

def getClassMethodList(classNode, skipMethods):
    className = getText("compoundname", classNode)
    shortClassName=stripOpenmmPrefix(className)
    methodList=[]
    for section in findNodes(classNode, "sectiondef", kind="public-static-func")+findNodes(classNode, "sectiondef", kind="public-func"):
        for memberNode in findNodes(section, "memberdef", kind="function", prot="public"):
            methDefinition = getText("definition", memberNode)
            shortMethDefinition=stripOpenmmPrefix(methDefinition)
            methName=shortMethDefinition.split()[-1]
            if (shortClassName, methName) in skipMethods: continue
            numParams=len(findNodes(memberNode, 'param'))
            if (shortClassName, methName, numParams) in skipMethods: continue
            for catchString in ['Factory', 'Impl', 'Info', 'Kernel']:
                if shortClassName.endswith(catchString):
                    sys.stderr.write("Warning: Including class %s\n" %
                                     shortClassName)
                    continue

            if (shortClassName, methName) in skipMethods: continue

            # set template info

            templateType = getText("templateparamlist/param/type", memberNode)
            templateName = getText("templateparamlist/param/declname", memberNode)

            methodList.append( (shortClassName,
                                memberNode,
                                shortMethDefinition,
                                methName,
                                shortClassName==methName,
                                '~'+shortClassName==methName, templateType, templateName ) )
    return methodList


def docstringTypemap(cpptype):
    """Translate a C++ type to Python for inclusion in the Python docstrings.
    This doesn't need to be perfectly accurate -- it's not used for generating
    the actual swig wrapper code. It's only used for generating the docstrings.
    """
    pytype = cpptype
    if pytype.startswith('const '):
        pytype = pytype[6:]
    if pytype.startswith('std::'):
        pytype = pytype[5:]
    pytype = pytype.strip('&')
    return pytype.strip()


class SwigInputBuilder:
    def __init__(self,
                 inputDirname,
                 configFilename,
                 outputFilename=None,
                 docstringFilename=None,
                 pythonprependFilename=None,
                 pythonappendFilename=None,
                 skipAdditionalMethods=[],
                 SWIG_VERSION='3.0.2'):
        self.nodeByID={}
        self.SWIG_COMPACT_ARGUMENTS = LooseVersion(SWIG_VERSION) < LooseVersion('3.0.5')

        self.configModule = __import__(os.path.splitext(configFilename)[0])

        self.skipMethods=self.configModule.SKIP_METHODS[:]
        for skipMethod in skipAdditionalMethods:
            items=skipMethod.split('::')
            if len(items)==3:
                items[2]=int(items[2])
            self.skipMethods.append(tuple(items))

        # Read all the XML files and merge them into a single document.
        self.doc = etree.ElementTree(etree.Element('root'))
        for file in os.listdir(inputDirname):
            if file.lower().endswith('xml'):
                root = etree.parse(os.path.join(inputDirname, file)).getroot()
                for node in root:
                    self.doc.getroot().append(node)

        if outputFilename:
            self.fOut = open(outputFilename, 'w')
        else:
            self.fOut = sys.stdout

        if docstringFilename:
            self.fOutDocstring = open(docstringFilename, 'w')
        else:
            self.fOutDocstring = None

        if pythonprependFilename:
            self.fOutPythonprepend = open(pythonprependFilename, 'w')
        else:
            self.fOutPythonprepend = None

        if pythonappendFilename:
            self.fOutPythonappend = open(pythonappendFilename, 'w')
        else:
            self.fOutPythonappend = None

        self._enumByClassname={}

        self._orderedClassNodes=self._buildOrderedClassNodes()

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
        if node in excludedClassNodes: return
        nodeName = getText("compoundname", node)
        if (nodeName.split("::")[-1],) in self.skipMethods:
            return
        for baseNodePnt in findNodes(node, "basecompoundref", prot="public"):
            if "refid" in baseNodePnt.attrib:
                baseNodeID = baseNodePnt.attrib["refid"]
                baseNode = self._getNodeByID(baseNodeID)
                self._findBaseNodes(baseNode, excludedClassNodes)
        excludedClassNodes.append(node)


    def writeFactories(self):
        self.fOut.write("\n/* Declare factories */\n\n")
        forceSubclassList = []
        integratorSubclassList = []
        tabulatedFunctionSubclassList = []
        for classNode in findNodes(self.doc.getroot(), "compounddef", kind="class", prot="public"):
            className = getText("compoundname", classNode)
            shortClassName=stripOpenmmPrefix(className)
            if (className.split("::")[-1],) in self.skipMethods:
                continue
            for baseNodePnt in findNodes(classNode, "basecompoundref", prot="public"):
                if "refid" in baseNodePnt.attrib:
                    baseNodeID=baseNodePnt.attrib["refid"]
                    baseNode=self._getNodeByID(baseNodeID)
                    baseName = getText("compoundname", baseNode)
                    if baseName == 'OpenMM::Force':
                        forceSubclassList.append(shortClassName)
                    elif baseName == 'OpenMM::Integrator':
                        integratorSubclassList.append(shortClassName)
                    elif baseName == 'OpenMM::TabulatedFunction':
                        tabulatedFunctionSubclassList.append(shortClassName)

        self.fOut.write("%factory(OpenMM::Force& OpenMM::System::getForce")
        for name in sorted(forceSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Force* OpenMM::Force::__copy__")
        for name in sorted(forceSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Force* OpenMM_XmlSerializer__deserializeForce")
        for name in sorted(forceSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Force& OpenMM::CustomCVForce::getCollectiveVariable")
        for name in sorted(forceSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Integrator* OpenMM::Integrator::__copy__")
        for name in sorted(integratorSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Integrator* OpenMM_XmlSerializer__deserializeIntegrator")
        for name in sorted(integratorSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Integrator& OpenMM::Context::getIntegrator")
        for name in sorted(integratorSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::Integrator& OpenMM::CompoundIntegrator::getIntegrator")
        for name in sorted(integratorSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::TabulatedFunction* OpenMM::TabulatedFunction::__copy__")
        for name in sorted(tabulatedFunctionSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::TabulatedFunction* OpenMM_XmlSerializer__deserializeTabulatedFunction")
        for name in sorted(tabulatedFunctionSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

        for classNode in self._orderedClassNodes:
            methodList=getClassMethodList(classNode, self.skipMethods)
            for items in methodList:
                (shortClassName, memberNode,
                 shortMethDefinition, methName,
                 isConstructors, isDestructor, templateType, templateName) = items
                if shortMethDefinition == 'TabulatedFunction& getTabulatedFunction':
                    self.fOut.write("%factory(OpenMM::TabulatedFunction& OpenMM::")
                    self.fOut.write("%s::%s" % (shortClassName, methName))
                    for name in sorted(tabulatedFunctionSubclassList):
                        self.fOut.write(",\n         OpenMM::%s" % name)
                    self.fOut.write(");\n\n")

        self.fOut.write("%factory(OpenMM::VirtualSite& OpenMM::System::getVirtualSite, OpenMM::TwoParticleAverageSite, OpenMM::ThreeParticleAverageSite, OpenMM::OutOfPlaneSite, OpenMM::LocalCoordinatesSite);\n\n")
        self.fOut.write("\n")

    def writeGlobalConstants(self):
        self.fOut.write("/* Global Constants */\n\n")
        node = next((x for x in findNodes(self.doc.getroot(), "compounddef", kind="namespace") if x.findtext("compoundname") == "OpenMM"))
        for section in findNodes(node, "sectiondef", kind="var"):
            for memberNode in findNodes(section, "memberdef", kind="variable", mutable="no", prot="public", static="yes"):
                vDef = stripOpenmmPrefix(getText("definition", memberNode))
                iDef = getText("initializer", memberNode)
                if iDef.startswith("="):
                    iDef = iDef[1:]
                self.fOut.write("static %s = %s;\n" % (vDef, iDef))
        self.fOut.write("\n")

    def writeForwardDeclarations(self):
        self.fOut.write("\n/* Forward Declarations */\n\n")

        for classNode in self._orderedClassNodes:
            hasConstructor=False
            methodList=getClassMethodList(classNode, self.skipMethods)
            for items in methodList:
                (shortClassName, memberNode,
                 shortMethDefinition, methName,
                 isConstructors, isDestructor, templateType, templateName) = items
                if isConstructors:
                    hasConstructor=True

            className = stripOpenmmPrefix(getText("compoundname", classNode))
            # If has a constructor then tell swig tell to make a copy method
            if hasConstructor:
                self.fOut.write("%%copyctor %s ;\n" % className)
            self.fOut.write("class %s ;\n" % className)
        self.fOut.write("\n")

    def writeClassDeclarations(self):
        self.fOut.write("\n/* Class Declarations */\n\n")
        for classNode in self._orderedClassNodes:
            className = stripOpenmmPrefix(getText("compoundname", classNode))
            if self.fOutDocstring:
                dNode = classNode.find('detaileddescription')
                if dNode is not None:
                    docstring = getNodeText(dNode).strip().replace('"', '\\"')
                    docstring = striphtmltags(docstring)
                    self.fOutDocstring.write('%%feature("docstring") %s "%s";\n' % (className, docstring))
            self.fOut.write("class %s" % className)
            if className in self.configModule.MISSING_BASE_CLASSES:
                self.fOut.write(" : public %s" %
                                self.configModule.MISSING_BASE_CLASSES[className])

            for baseNodePnt in findNodes(classNode, "basecompoundref", prot="public"):
                if "refid" in baseNodePnt.attrib:
                    baseName = stripOpenmmPrefix(getText(".", baseNodePnt))
                    self.fOut.write(" : public %s" % baseName)
            self.fOut.write(" {\n")
            self.fOut.write("public:\n")
            self.writeEnumerations(classNode)
            self.writeMethods(classNode)
            self.fOut.write("};\n\n")
        self.fOut.write("\n")

    def writeEnumerations(self, classNode):
        enumNodes = []
        for section in findNodes(classNode, "sectiondef", kind="public-type"):
            for node in findNodes(section, "memberdef", kind="enum", prot="public"):
                enumNodes.append(node)
        for enumNode in enumNodes:
            className = getText("compoundname", classNode)
            shortClassName=stripOpenmmPrefix(className)
            enumName = getText("name", enumNode)
            try:
                self._enumByClassname[shortClassName].append(enumName)
            except KeyError:
                self._enumByClassname[shortClassName]=[enumName]
            self.fOut.write("%senum %s {" % (INDENT, enumName))
            argSep="\n"
            for valueNode in findNodes(enumNode, "enumvalue", prot="public"):
                vName = getText("name", valueNode)
                vInit = getText("initializer", valueNode)
                if vInit.startswith("="):
                    vInit = vInit[1:]
                self.fOut.write("%s%s%s = %s" % (argSep, 2*INDENT, vName, vInit))
                argSep=",\n"
            self.fOut.write("\n%s};\n" % INDENT)
        if len(enumNodes)>0: self.fOut.write("\n")

    def writeMethods(self, classNode):
        methodList=getClassMethodList(classNode, self.skipMethods)

        #write only Constructors
        for items in methodList:
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName) = items
            if isConstructors:
                mArgsstring = getText("argsstring", memberNode)
                try:
                    pExceptions = " %s" % getText('exceptions', memberNode)
                except IndexError:
                    pExceptions = ""
                self.fOut.write("%s%s%s%s;\n" % (INDENT, shortMethDefinition,
                                                 mArgsstring, pExceptions))
        #write only Destructors
        for items in methodList:
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName) = items
            if isDestructor:
                mArgsstring = getText("argsstring", memberNode)
                try:
                    pExceptions = " %s" % getText('exceptions', memberNode)
                except IndexError:
                    pExceptions = ""
                self.fOut.write("%s%s%s%s;\n" % (INDENT, shortMethDefinition,
                                                 mArgsstring, pExceptions))

        #write only non Constructor and Destructor methods and python mods
        self.fOut.write("\n")
        methodsWithOutputArgs = set()
        for items in methodList:
            clearOutput=""
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName) = items
            if isConstructors or isDestructor: continue

            key = (shortClassName, methName)
            if key in self.configModule.DOC_STRINGS:
                self.fOut.write('%%feature("autodoc", "%s") %s;\n' %
                                (self.configModule.DOC_STRINGS[key], methName))

            paramList=findNodes(memberNode, 'param')
            for pNode in paramList:
                try:
                    pType = getText('type', pNode)
                except IndexError:
                    pType = getText('type/ref', pNode)
                pName = getText('declname', pNode)
                key = (shortClassName, methName, pName)
                if pType.find('&')>=0 and \
                   'const' not in pType.split():
                    if key not in self.configModule.NO_OUTPUT_ARGS:
                        eType = pType.split()[0]
                        if shortClassName in self._enumByClassname and \
                           eType in self._enumByClassname[shortClassName]:
                            simpleType = re.sub(eType, 'int', pType)
                        else:
                            simpleType = pType
                        self.fOut.write("%s%%apply %s OUTPUT { %s %s };\n" %
                                        (INDENT, simpleType, pType, pName))
                        clearOutput = "%s%s%%clear %s %s;\n" \
                                     % (clearOutput, INDENT, pType, pName)
                        methodsWithOutputArgs.add((shortClassName, methName))

            mArgsstring = getText("argsstring", memberNode)
            try:
                pExceptions = " %s" % getText('exceptions', memberNode)
            except IndexError:
                pExceptions = ""
            if memberNode.attrib["virt"].strip()!='non-virtual':
                if 'virtual' not in shortMethDefinition.split():
                    shortMethDefinition="virtual %s" % shortMethDefinition
            if( len(templateType) > 0 ):
                self.fOut.write("%stemplate<%s %s> %s%s%s;\n" % (INDENT, templateType, templateName, shortMethDefinition, mArgsstring, pExceptions))
            else:
                self.fOut.write("%s%s%s%s;\n" % (INDENT, shortMethDefinition, mArgsstring, pExceptions))
            if clearOutput:
                self.fOut.write(clearOutput)

        #write python mod blocks
        for items in methodList:
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName) = items
            paramList = findNodes(memberNode, 'param')

            # write pythonprepend blocks
            if isConstructors:
                mArgsstring = '' # specifying args to constructors seems to prevent append and prepend from working
            else:
                mArgsstring = getText("argsstring", memberNode)
            if self.fOutPythonprepend and \
               len(paramList) and \
               mArgsstring.find('=0') < 0:
                text = '''
%pythonprepend OpenMM::{shortClassName}::{methName}{mArgsstring} %{{{{{{0}}
%}}}}'''.format(shortClassName=shortClassName, methName=methName, mArgsstring=mArgsstring)
                textInside = ''
                key = (shortClassName, methName)
                for argNum in self.configModule.STEAL_OWNERSHIP.get(key, []):
                    if self.SWIG_COMPACT_ARGUMENTS:
                        argName = 'args[%s]' % argNum
                    else:
                        argName = getText('declname', paramList[argNum])

                    textInside += '''
    if not {argName}.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)'''.format(argName=argName)

                # Convert input arguments to the proper units, if specified.
                if key not in methodsWithOutputArgs:
                    if key in self.configModule.UNITS:
                        argUnits=self.configModule.UNITS[key][1]
                    elif ("*", methName) in self.configModule.UNITS:
                        argUnits=self.configModule.UNITS[("*", methName)][1]
                    else:
                        argUnits = ()
                    if len(argUnits) > 0 and (self.SWIG_COMPACT_ARGUMENTS or isConstructors):
                        textInside += '''
    args = list(args)'''
                    for i, units in enumerate(argUnits):
                        if units is not None:
                            if self.SWIG_COMPACT_ARGUMENTS or isConstructors:
                                argName = 'args[%s]' % i
                            else:
                                argName = getText('declname', paramList[i])
                            textInside += '''
    if unit.is_quantity({argName}):
        {argName} = {argName}.value_in_unit({units})'''.format(argName=argName, units=units)

                for argNum in self.configModule.REQUIRE_ORDERED_SET.get(key, []):
                    if self.SWIG_COMPACT_ARGUMENTS:
                        argName = 'args[%s]' % argNum
                    else:
                        argName = getText('declname', paramList[argNum])

                    textInside += '''
    {argName} = list({argName})'''.format(argName=argName)
                if textInside:
                    self.fOutPythonprepend.write(text.format(textInside))

            # write pythonappend blocks
            if self.fOutPythonappend \
               and mArgsstring.find('=0') < 0:
                key = (shortClassName, methName)
                #print "key %s %s \n" % (shortClassName, methName)
                addText=''
                returnType = getText("type", memberNode)

                if key in self.configModule.UNITS:
                    valueUnits=self.configModule.UNITS[key]
                elif ("*", methName) in self.configModule.UNITS:
                    valueUnits=self.configModule.UNITS[("*", methName)]
                elif methName.startswith('get') and returnType not in ('void', 'int', 'bool', 'std::string', 'const std::string &'):
                    s = 'do not know how to add units to %s %s::%s' \
                        % (returnType, shortClassName, methName)
                    raise Exception(s)
                else:
                    valueUnits=[None, ()]

                index=0
                if valueUnits[0] is not None:
                    sys.stdout.write("%s.%s() returns %s\n" %
                                     (shortClassName, methName, valueUnits[0]))
                    if len(valueUnits[1])>0:
                        addText = "%s%sval[%d]=unit.Quantity(val[%d], %s)\n" \
                                 % (addText, INDENT,
                                    index, index,
                                    valueUnits[0])
                        index+=1
                    else:
                        addText = "%s%sval=unit.Quantity(val, %s)\n" \
                                 % (addText, INDENT, valueUnits[0])

                for vUnit in valueUnits[1]:
                    if vUnit is not None and key in methodsWithOutputArgs:
                        addText = "%s%sval[%s]=unit.Quantity(val[%s], %s)\n" \
                                     % (addText, INDENT, index, index, vUnit)
                    index+=1

                if key in self.configModule.STEAL_OWNERSHIP:
                    for argNum in self.configModule.STEAL_OWNERSHIP[key]:
                        if self.SWIG_COMPACT_ARGUMENTS:
                            argName = 'args[%s]' % argNum
                        else:
                            argName = getText('declname', paramList[argNum])
                        addText = "%s%s%s.thisown=0\n" \
                                % (addText, INDENT, argName)

                if addText:
                    self.fOutPythonappend.write("%pythonappend")
                    self.fOutPythonappend.write(" OpenMM::%s::%s(" % key)
                    sepChar=''
                    outputIndex=0
                    for pNode in paramList:
                        try:
                            pType = getText('type', pNode)
                        except IndexError:
                            pType = getText('type/ref', pNode)
                        pName = getText('declname', pNode)
                        self.fOutPythonappend.write("%s%s %s" % (sepChar, pType, pName))
                        sepChar=', '

                        if pType.find('&')>=0 and \
                          'const' not in pType.split() and \
                          key not in self.configModule.NO_OUTPUT_ARGS and \
                          len(valueUnits[1])>0:
                            try:
                                unitType=valueUnits[1][outputIndex]
                            except IndexError:
                                s = "missing unit type for %s.%s() arg named %s" \
                                   % (shortClassName, methName, pName)
                                raise Exception(s)
                            sys.stdout.write("%s.%s() returns %s as %s\n" %
                                             (shortClassName, methName,
                                              pName, unitType))
                            outputIndex+=1
                    if memberNode.attrib["const"]=="yes":
                        constString=' const'
                    else:
                        constString=''
                    self.fOutPythonappend.write(")%s %%{\n" % constString)
                    self.fOutPythonappend.write(addText)
                    self.fOutPythonappend.write("%}\n\n")


        #print "Done python mod blocks\n"
        #write Docstring info
        for items in methodList:
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName ) = items

            if self.fOutDocstring:
                signatureParams = findNodes(memberNode, 'param')
                assert len(findNodes(memberNode, 'detaileddescription')) == 1
                dNode = findNodes(memberNode, 'detaileddescription')[0]

                try:
                    description=getText('para', dNode)
                    description.strip()
                except IndexError:
                    description = ''
                params = findNodes(dNode, 'para/parameterlist/parameteritem')

                paramString = ['Parameters', '----------']
                returnString = ['Returns', '-------']

                if len(params) > 0:
                    if len(signatureParams) != len(params):
                        raise ValueError('docstring in %s.%s does not match the signature' % (shortClassName, methName))

                    for pNode, pSignatureNode in zip(params, signatureParams):
                        parameterNameNode = findNodes(pNode, 'parameternamelist/parametername')[0]
                        argDoc = getText('parameterdescription/para', pNode)
                        argName = getNodeText(parameterNameNode)
                        argType = docstringTypemap(getText('type', pSignatureNode))

                        isOutput = parameterNameNode.get('direction') == 'out'
                        if isOutput:
                            returnString.extend(['%s : %s' % (argName, argType), '    %s' % argDoc])
                        else:
                            paramString.extend(['%s : %s' % (argName, argType), '    %s' % argDoc])


                returnSection = findNodes(dNode, 'para/simplesect')
                if len(returnSection) > 0:
                    returnNode = returnSection[0]
                    if returnNode.get('kind') == 'return':
                        argType = getNodeText(findNodes(memberNode, 'type')[0])
                        argType = docstringTypemap(argType)
                        returnString.extend([argType, '    %s' % getNodeText(returnNode).strip()])

                dString = '\n'.join(
                    ([description] + [''] if len(description) > 0 else []) +
                    (paramString + [''] if len(paramString) > 2 else []) +
                    (returnString if len(returnString) > 2 else [])).strip()
                if dString:
                    dString = re.sub(r'([^\\])"', r'\g<1>\"', dString)
                    s = '%%feature("docstring") OpenMM::%s::%s "%s";' \
                       % (shortClassName, methName, dString)
                    self.fOutDocstring.write("%s\n" % s)

                self.fOutDocstring.write("\n\n")


    def writeSwigFile(self):
        self.fOut.write("/* Swig input file,\n")
        self.fOut.write("%sgenerated by %s on %s\n*/\n\n\n"
                        % (INDENT, sys.argv[0], time.asctime()))
        self.fOut.write("\nnamespace OpenMM {\n\n")
        self.writeFactories()
        self.writeGlobalConstants()
        self.writeForwardDeclarations()
        self.writeClassDeclarations()
        self.fOut.write("\n} // namespace OpenMM\n\n")


def parseCommandLine():
    opts, args_proper = getopt.getopt(sys.argv[1:], 'hi:c:o:d:a:z:s:v:')
    inputDirname = None
    configFilename = None
    outputFilename = ""
    docstringFilename = ""
    pythonprependFilename = ""
    pythonappendFilename = ""
    skipAdditionalMethods = []
    swigVersion = '3.0.2'
    for option, parameter in opts:
        if option=='-h': usageError()
        if option=='-i': inputDirname = parameter
        if option=='-c': configFilename=parameter
        if option=='-o': outputFilename = parameter
        if option=='-d': docstringFilename = parameter
        if option=='-a': pythonprependFilename=parameter
        if option=='-z': pythonappendFilename=parameter
        if option=='-s': skipAdditionalMethods.append(parameter)
        if option=='-v': swigVersion = parameter
    if not inputDirname: usageError()
    if not configFilename: usageError()
    return (args_proper, inputDirname, configFilename, outputFilename,
            docstringFilename, pythonprependFilename, pythonappendFilename,
            skipAdditionalMethods, swigVersion)

def main():
    (args_proper, inputDirname, configFilename, outputFilename,
     docstringFilename, pythonprependFilename, pythonappendFilename,
     skipAdditionalMethods, swigVersion) = parseCommandLine()
    sBuilder = SwigInputBuilder(inputDirname, configFilename, outputFilename,
                                docstringFilename, pythonprependFilename,
                                pythonappendFilename, skipAdditionalMethods,
                                swigVersion)
    #print "Calling writeSwigFile\n"
    sBuilder.writeSwigFile()
    #print "Done writeSwigFile\n"

    return


def usageError():
    sys.stdout.write('usage: %s -i xmlHeadersFilename \\\n' \
         % os.path.basename(sys.argv[0]))
    sys.stdout.write('       %s -c inputConfigFilename\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-o swigInputDirname]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-d docstringFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-a pythonprependFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-z pythonappendFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-s skippedClasses]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-v swigVersion]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.exit(1)

if __name__=='__main__':
    main()


