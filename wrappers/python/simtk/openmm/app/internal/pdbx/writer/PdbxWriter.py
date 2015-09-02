##
# File:  PdbxWriter.py
# Date:  2011-10-09 Jdw  Adapted from PdbxParser.py
#
# Updates:
#    5-Apr-2011 jdw  Using the double quote format preference
#   23-Oct-2012 jdw  update path details and reorganize.
#
###
"""
Classes for writing data and dictionary containers in PDBx/mmCIF format.

"""
from __future__ import absolute_import
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.01"


from simtk.openmm.app.internal.pdbx.reader.PdbxContainers import *

class PdbxError(Exception):
    """ Class for catch general errors 
    """
    pass

class PdbxWriter(object):
    """Write PDBx data files or dictionaries using the input container
       or container list.
    """
    def __init__(self,ofh=sys.stdout):
        self.__ofh=ofh
        self.__containerList=[]
        self.__MAXIMUM_LINE_LENGTH = 2048
        self.__SPACING = 2
        self.__INDENT_DEFINITION = 3
        self.__indentSpace = " " * self.__INDENT_DEFINITION        
        self.__doDefinitionIndent=False
        # Maximum number of rows checked for value length and format
        self.__rowPartition=None

    def setRowPartition(self,numRows):
        ''' Maximum number of rows checked for value length and format
        '''
        self.__rowPartition=numRows
        
    def write(self, containerList):
        self.__containerList=containerList
        for  container in self.__containerList:
            self.writeContainer(container)

    def writeContainer(self,container):
        indS=" " *  self.__INDENT_DEFINITION
        if isinstance(container, DefinitionContainer):
            self.__write("save_%s\n" % container.getName())
            self.__doDefinitionIndent=True
            self.__write(indS+"#\n")            
        elif isinstance(container, DataContainer):
            if (container.getGlobal()):
                self.__write("global_\n")
                self.__doDefinitionIndent=False
                self.__write("\n")            
            else:
                self.__write("data_%s\n" % container.getName())
                self.__doDefinitionIndent=False
                self.__write("#\n")                            
        
        for nm in container.getObjNameList():
            obj=container.getObj(nm)
            objL=obj.getRowList()

            # Skip empty objects
            if len(objL) == 0:
                continue

            # Item - value formattting 
            elif len(objL) == 1:
                self.__writeItemValueFormat(obj)

            # Table formatting - 
            elif len(objL) > 1 and len(obj.getAttributeList()) > 0:
                self.__writeTableFormat(obj)
            else:
                raise PdbxError()

            if self.__doDefinitionIndent:
                self.__write(indS+"#")
            else:
                self.__write("#")

        # Add a trailing saveframe reserved word 
        if isinstance(container, DefinitionContainer):
            self.__write("\nsave_\n")
        self.__write("#\n")            

    def __write(self, st):
        self.__ofh.write(st)

    def __writeItemValueFormat(self, myCategory):

        # Compute the maximum item name length within this category - 
        attributeNameLengthMax  = 0
        for attributeName in myCategory.getAttributeList():
            attributeNameLengthMax = max(attributeNameLengthMax, len(attributeName))
        itemNameLengthMax = self.__SPACING + len(myCategory.getName()) + attributeNameLengthMax + 2
        #
        lineList=[]
        lineList.append("#\n")        
        for attributeName,iPos in myCategory.getAttributeListWithOrder():        
            if  self.__doDefinitionIndent:
                #        - add indent --                
                lineList.append(self.__indentSpace)
                
            itemName = "_%s.%s" % (myCategory.getName(), attributeName)                
            lineList.append(itemName.ljust(itemNameLengthMax))
            
            lineList.append(myCategory.getValueFormatted(attributeName,0))            
            lineList.append("\n")
            
        self.__write("".join(lineList))

    def __writeTableFormat(self, myCategory):
        
        # Write the declaration of the loop_
        #
        lineList=[]
        lineList.append('#\n')
        if  self.__doDefinitionIndent:            
            lineList.append(self.__indentSpace)
        lineList.append("loop_")
        for attributeName in myCategory.getAttributeList():
            lineList.append('\n')            
            if  self.__doDefinitionIndent:            
                lineList.append(self.__indentSpace)
            itemName = "_%s.%s" % (myCategory.getName(), attributeName)
            lineList.append(itemName)
        self.__write("".join(lineList))

        #
        # Write the data in tabular format - 
        #
        #print myCategory.getName()
        #print myCategory.getAttributeList()

        #    For speed make the following evaluation on a portion of the table
        if self.__rowPartition is not None:
            numSteps=max(1,myCategory.getRowCount()/self.__rowPartition)
        else:
            numSteps=1
            
        formatTypeList,dataTypeList=myCategory.getFormatTypeList(steps=numSteps)
        maxLengthList=myCategory.getAttributeValueMaxLengthList(steps=numSteps)
        spacing   = " " * self.__SPACING
        #

        #print formatTypeList
        #print dataTypeList
        #print maxLengthList
        #
        for iRow in range(myCategory.getRowCount()):
            lineList     = []            
            lineList.append('\n')
            if  self.__doDefinitionIndent:
                lineList.append(self.__indentSpace + "  ")

            for iAt in range(myCategory.getAttributeCount()):
                formatType = formatTypeList[iAt]
                maxLength  = maxLengthList[iAt]

                if (formatType == 'FT_UNQUOTED_STRING' or formatType == 'FT_NULL_VALUE'):
                    val=myCategory.getValueFormattedByIndex(iAt,iRow)
                    lineList.append(val.ljust(maxLength))

                elif formatType == 'FT_NUMBER':
                    val=myCategory.getValueFormattedByIndex(iAt,iRow)
                    lineList.append(val.rjust(maxLength))
                    
                elif formatType == 'FT_QUOTED_STRING':
                    val=myCategory.getValueFormattedByIndex(iAt,iRow)

                    lineList.append(val.ljust(maxLength+2))

                elif formatType == "FT_MULTI_LINE_STRING":
                    val=myCategory.getValueFormattedByIndex(iAt,iRow)
                    lineList.append(val)
                    
                lineList.append(spacing)                    

            self.__write("".join(lineList))
        self.__write("\n")
