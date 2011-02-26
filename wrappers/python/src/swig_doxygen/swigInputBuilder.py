#!/bin/env python
#
#
 
"""Build swig imput file from xml encoded header files (see gccxml)."""
__author__ = "Randall J. Radmer"
__version__ = "1.0"
  
 
import sys, os
import time
import getopt
import re
import xml.dom.minidom as minidom
import xpath

#

INDENT = "   ";


def getText(subNodePath, node):
    xPath="%s/text() | %s/ref/text()" % (subNodePath, subNodePath)
    subNodes = xpath.find(xPath, node)
    s=""
    for n in xpath.find(xPath, node):
        s="%s%s" % (s, n.data)
    return s.strip()

OPENMM_RE_PATTERN=re.compile("(.*)OpenMM:[a-zA-Z:]*:(.*)")
def stripOpenmmPrefix(name, rePattern=OPENMM_RE_PATTERN):
    m=rePattern.search(name)
    rValue = "%s%s" % m.group(1,2)
    rValue.strip()
    return rValue


def getClassMethodList(classNode, skipMethods):
    className = getText("compoundname", classNode)
    shortClassName=stripOpenmmPrefix(className)
    xPath1 = "sectiondef[@kind='public-static-func']/memberdef[@kind='function' and @prot='public']"
    xPath2 = "sectiondef[@kind='public-func']/memberdef[@kind='function' and @prot='public']"
    methodList=[]
    for memberNode in xpath.find(xPath1, classNode) \
                     +xpath.find(xPath2, classNode):
        methDefinition = getText("definition", memberNode)
        shortMethDefinition=stripOpenmmPrefix(methDefinition)
        methName=shortMethDefinition.split()[-1]
        if (shortClassName, methName) in skipMethods: continue
        numParams=len(xpath.find('param', memberNode))
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


class SwigInputBuilder:
    def __init__(self,
                 inputFilename,
                 configFilename,
                 outputFilename=None,
                 docstringFilename=None,
                 pythonprependFilename=None,
                 pythonappendFilename=None,
                 skipAdditionalMethods=[]):
        self.nodeByID={}

        self.configModule = __import__(os.path.splitext(configFilename)[0])

        self.skipMethods=self.configModule.SKIP_METHODS[:]
        for skipMethod in skipAdditionalMethods:
            items=skipMethod.split('::')
            if len(items)==3:
                items[2]=int(items[2])
            self.skipMethods.append(tuple(items))

        self.doc = minidom.parse(inputFilename)

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
            xPath = "/doxygen/compounddef[@id='%s']" % id
            self.nodeByID[id] = xpath.find(xPath, self.doc)[0]
        return self.nodeByID[id]

    def _buildOrderedClassNodes(self):
        orderedClassNodes=[]
        xPath = "/doxygen/compounddef[@kind='class' and @prot='public']"
        for node in xpath.find(xPath, self.doc):
            self._findBaseNodes(node, orderedClassNodes)
        return orderedClassNodes

    def _findBaseNodes(self, node, excludedClassNodes=[]):
        if node in excludedClassNodes: return
        nodeName = getText("compoundname", node)
        if (nodeName.split("::")[-1],) in self.skipMethods:
            return
        xPath = "basecompoundref[@prot='public']"
        for baseNodePnt in xpath.find(xPath, node):
            baseNodeID=baseNodePnt.getAttribute("refid")
            baseNode=self._getNodeByID(baseNodeID)
            self._findBaseNodes(baseNode, excludedClassNodes)
        excludedClassNodes.append(node)


    def writeForceSubclassList(self):
        self.fOut.write("\n/* Force subclasses */\n\n")
        forceSubclassList=[]
        xPath = "/doxygen/compounddef[@kind='class' and @prot='public']"
        for classNode in xpath.find(xPath, self.doc):
            className = getText("compoundname", classNode)
            shortClassName=stripOpenmmPrefix(className)
            if (className.split("::")[-1],) in self.skipMethods:
                continue
            #print className
            #print classNode.toxml()
            xPath = "basecompoundref[@prot='public']"
            for baseNodePnt in xpath.find(xPath, classNode):
                baseNodeID=baseNodePnt.getAttribute("refid")
                baseNode=self._getNodeByID(baseNodeID)
                baseName = getText("compoundname", baseNode)
                if baseName == 'OpenMM::Force':
                    forceSubclassList.append(shortClassName)
        self.fOut.write("%factory(OpenMM::Force& OpenMM::System::getForce")
        for name in sorted(forceSubclassList):
            self.fOut.write(",\n         OpenMM::%s" % name)
        self.fOut.write(");\n\n")

#        for classNode in self._orderedClassNodes:
#            className = stripOpenmmPrefix(getText("compoundname", classNode))
#            self.fOut.write("class %s ;\n" % className)
        self.fOut.write("\n")

    def writeGlobalConstants(self):
        self.fOut.write("/* Global Constants */\n\n")
        xPath = "/doxygen/compounddef[@kind='namespace' and compoundname='OpenMM']"
        node = xpath.find(xPath, self.doc)[0]
        xPath="sectiondef[@kind='var']/memberdef[@kind='variable' and @mutable='no' and @prot='public' and @static='yes']"
        for memberNode in xpath.find(xPath, node):
            vDef = stripOpenmmPrefix(getText("definition", memberNode))
            iDef = getText("initializer", memberNode)
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
            self.fOut.write("class %s" % className)
            if className in self.configModule.MISSING_BASE_CLASSES:
                self.fOut.write(" : public %s" %
                                self.configModule.MISSING_BASE_CLASSES[className])

            xPath = "basecompoundref[@prot='public']"
            for baseNodePnt in xpath.find(xPath, classNode):
                baseName = stripOpenmmPrefix(getText(".", baseNodePnt))
                self.fOut.write(" : public %s" % baseName)
            self.fOut.write(" {\n")
            self.fOut.write("public:\n")
            self.writeEnumerations(classNode)
            self.writeMethods(classNode)
            self.fOut.write("};\n\n")
        self.fOut.write("\n")

    def writeEnumerations(self, classNode):
        xPath = "sectiondef[@kind='public-type']/memberdef[@kind='enum' and @prot='public']"
        enumNodes=xpath.find(xPath, classNode)
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
            for valueNode in xpath.find("enumvalue[@prot='public']", enumNode):
                vName = getText("name", valueNode)
                vInit = getText("initializer", valueNode)
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

            paramList=xpath.find('param', memberNode)
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

            mArgsstring = getText("argsstring", memberNode)
            try:
                pExceptions = " %s" % getText('exceptions', memberNode)
            except IndexError:
                pExceptions = ""
            if memberNode.getAttribute("virt").strip()!='non-virtual':
                if 'virtual' not in shortMethDefinition.split():
                    shortMethDefinition="virtual %s" % shortMethDefinition
            if( len(templateType) > 0 ):
                self.fOut.write("%stemplate<%s %s> %s%s%s;\n" % (INDENT, templateType, templateName, shortMethDefinition, mArgsstring, pExceptions))
            else:
                self.fOut.write("%s%s%s%s;\n" % (INDENT, shortMethDefinition, mArgsstring, pExceptions))
            if clearOutput:
                self.fOut.write(clearOutput)

        isXmlSerializer = 0
        isSystem        = 0
        for items in methodList:
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName) = items
            if( shortClassName == 'XmlSerializer' ):isXmlSerializer = 1
            if( shortClassName == 'System' ):isSystem               = 1

        if( isXmlSerializer == 1 ):
            extendString =     "%extend {\n"
            extendString +=    "       static std::string serializeSystem( const OpenMM::System *object ){\n"
            extendString +=    "           std::stringstream ss;\n"
            extendString +=    "           XmlSerializer::serialize<OpenMM::System>( object, \"System\", ss );\n"
            extendString +=    "           return ss.str();\n"
            extendString +=    "       }\n"
            extendString +=    "\n"
            extendString +=    "       static OpenMM::System* deserializeSystem( const char* inputString ){\n"
            extendString +=    "           std::stringstream ss;\n"
            extendString +=    "           ss << inputString;\n" 
            extendString +=    "           return XmlSerializer::deserialize<OpenMM::System>( ss );\n"
            extendString +=    "       }\n"
            extendString +=    "    };\n"   

            self.fOut.write("%s%s\n" % (INDENT, extendString))

        #write python mod blocks
        for items in methodList:
            (shortClassName, memberNode,
             shortMethDefinition, methName,
             isConstructors, isDestructor, templateType, templateName) = items
            paramList=xpath.find('param', memberNode)

            #write pythonprepend blocks
            mArgsstring = getText("argsstring", memberNode)
            if self.fOutPythonprepend and \
               len(paramList) and \
               mArgsstring.find('=0')<0:
                key=(shortClassName, methName)
                if key in self.configModule.STEAL_OWNERSHIP:
                    for argNum in self.configModule.STEAL_OWNERSHIP[key]:
                        self.fOutPythonprepend.write("%pythonprepend")
                        self.fOutPythonprepend.write(" OpenMM::%s::%s%s %%{\n"
                                                     % (shortClassName,
                                                        methName,
                                                        mArgsstring))
                        self.fOutPythonprepend.write(
                                         "%sif not args[%s].thisown:\n"
                                         % (INDENT, argNum))
                        s = 's = "the %s object does not own its'
                        s = '%s corresponding OpenMM object" \\' % s
                        self.fOutPythonprepend.write("%s   %s\n" % (INDENT, s))

                        s = '   %% args[%s].__class__.__name__' % argNum
                        self.fOutPythonprepend.write("%s   %s\n" % (INDENT, s))

                        s = "raise Exception(s)"
                        self.fOutPythonprepend.write("%s   %s\n" % (INDENT, s))

                        self.fOutPythonprepend.write("%}\n\n")

            #write pythonappend blocks
            if self.fOutPythonappend \
               and mArgsstring.find('=0')<0:
                key=(shortClassName, methName)
                #print "key %s %s \n" % (shortClassName, methName)
                addText=''

                if key in self.configModule.UNITS:
                    valueUnits=self.configModule.UNITS[key]
                elif ("*", methName) in self.configModule.UNITS:
                    valueUnits=self.configModule.UNITS[("*", methName)]
                elif methName.startswith('get'):
                    s = 'do not know how to add units to %s::%s' \
                        % (shortClassName, methName)
                    raise Exception(s)
                else:
                    valueUnits=[None, ()]

                index=0
                if valueUnits[0]:
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
                        if vUnit!=None:
                            addText = "%s%sval[%s]=unit.Quantity(val[%s], %s)\n" \
                                     % (addText, INDENT, index, index, vUnit)
                        index+=1

                if key in self.configModule.STEAL_OWNERSHIP:
                    for argNum in self.configModule.STEAL_OWNERSHIP[key]:
                        addText = "%s%sargs[%s].thisown=0\n" \
                                % (addText, INDENT, argNum)

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
                    if memberNode.getAttribute("const")=="yes":
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
                for dNode in xpath.find('detaileddescription', memberNode):
                    dString=""
                    try:
                        description=getText('para', dNode)
                        description.strip()
                        if description:
                            dString=description
                    except IndexError:
                        pass
                    for pNode in xpath.find('para/parameterlist/parameteritem', dNode):
                        argName = getText('parameternamelist/parametername', pNode)
                        argDoc = getText('parameterdescription/para', pNode)
                        dString="%s\n   %s -- %s" % (dString, argName, argDoc)
                        dString.strip()
                    if dString:
                        dString=re.sub(r'([^\\])"', r'\g<1>\"', dString)
                        s = '%%feature("docstring") OpenMM::%s::%s "%s";' \
                           % (shortClassName, methName, dString)
                        self.fOutDocstring.write("%s\n" % s)
                self.fOutDocstring.write("\n\n")
        #print "Done write Docstring info\n"



    def writeSwigFile(self):
        self.fOut.write("/* Swig input file,\n")
        self.fOut.write("%sgenerated by %s on %s\n*/\n\n\n"
                        % (INDENT, sys.argv[0], time.asctime()))
        self.fOut.write("\nnamespace OpenMM {\n\n")
        self.writeForceSubclassList()
        self.writeGlobalConstants()
        self.writeForwardDeclarations()
        self.writeClassDeclarations()
        self.fOut.write("\n} // namespace OpenMM\n\n")


def parseCommandLine():
    opts, args_proper = getopt.getopt(sys.argv[1:], 'hi:c:o:d:a:z:s:')
    inputFilename = None
    configFilename = None
    outputFilename = ""
    docstringFilename = ""
    pythonprependFilename = ""
    pythonappendFilename = ""
    skipAdditionalMethods = []
    for option, parameter in opts:
        if option=='-h': usageError()
        if option=='-i': inputFilename = parameter
        if option=='-c': configFilename=parameter
        if option=='-o': outputFilename = parameter
        if option=='-d': docstringFilename = parameter
        if option=='-a': pythonprependFilename=parameter
        if option=='-z': pythonappendFilename=parameter
        if option=='-s': skipAdditionalMethods.append(parameter)
    if not inputFilename: usageError()
    if not configFilename: usageError()
    return (args_proper, inputFilename, configFilename, outputFilename,
            docstringFilename,
            pythonprependFilename, pythonappendFilename, skipAdditionalMethods)

def main():
    (args_proper, inputFilename, configFilename, outputFilename,
     docstringFilename, pythonprependFilename, pythonappendFilename,
     skipAdditionalMethods) = parseCommandLine()
    sBuilder = SwigInputBuilder(inputFilename, configFilename, outputFilename,
                                docstringFilename, pythonprependFilename,
                                pythonappendFilename, skipAdditionalMethods)
    #print "Calling writeSwigFile\n"
    sBuilder.writeSwigFile()
    #print "Done writeSwigFile\n"

    return


def usageError():
    sys.stdout.write('usage: %s -i xmlHeadersFilename \\\n' \
         % os.path.basename(sys.argv[0]))
    sys.stdout.write('       %s -c inputConfigFilename\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-o swigInputFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-d docstringFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-a pythonprependFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-z pythonappendFilename]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.stdout.write('       %s[-s skippedClasses]\n' \
         % (' '*len(os.path.basename(sys.argv[0]))))
    sys.exit(1)

if __name__=='__main__':
    main()


