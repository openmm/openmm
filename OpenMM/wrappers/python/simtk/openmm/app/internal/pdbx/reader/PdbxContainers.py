##
#
# File:     PdbxContainers.py
# Original: 02-Feb-2009   jdw
#
# Update:
#   23-Mar-2011   jdw Added method to rename attributes in category containers.
#   05-Apr-2011   jdw Change cif writer to select double quoting as preferred 
#                     quoting style where possible.
#   16-Jan-2012   jdw Create base class for DataCategory class
#   22-Mar-2012   jdw when append attributes to existing categories update
#                     existing rows with placeholder null values.
#    2-Sep-2012   jdw add option to avoid embedded quoting that might
#                     confuse simple parsers.
#   28-Jun-2013   jdw export remove method
#   29-Jun-2013   jdw export remove row method
##
"""

A collection of container classes supporting the PDBx/mmCIF storage model.

A base container class is defined which supports common features of
data and definition containers.   PDBx data files are organized in
sections called data blocks which are mapped to data containers.
PDBx dictionaries contain definition sections and data sections
which are mapped to definition and data containes respectively.

Data in both PDBx data files and dictionaries are organized in
data categories. In the PDBx syntax individual items or data
identified by labels of the form '_categoryName.attributeName'.
The terms category and attribute in PDBx jargon are analogous
table and column in relational data model, or class and attribute
in an object oriented data model.

The DataCategory class provides base storage container for instance
data and definition meta data.

"""
from __future__ import absolute_import

__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.01"

import re,sys,traceback

class CifName(object):
    ''' Class of utilities for CIF-style data names -
    '''
    def __init__(self):
        pass

    @staticmethod
    def categoryPart(name):
        tname=""
        if name.startswith("_"):
            tname=name[1:]
        else:
            tname=name

        i = tname.find(".")
        if i == -1:
            return tname
        else:
            return tname[:i]

    @staticmethod
    def attributePart(name):
        i = name.find(".")
        if i == -1:
            return None
        else:
            return name[i+1:]


class ContainerBase(object):
    ''' Container base class for data and definition objects.
    '''
    def __init__(self,name):
        # The enclosing scope of the data container (e.g. data_/save_)
        self.__name = name
        # List of category names within this container -  
        self.__objNameList=[]
        # dictionary of DataCategory objects keyed by category name.
        self.__objCatalog={}
        self.__type=None

    def getType(self):
        return self.__type

    def setType(self,type):
        self.__type=type
    
    def getName(self):
        return self.__name

    def setName(self,name):
        self.__name=name

    def exists(self,name):
        if name in self.__objCatalog:
            return True
        else:
            return False
        
    def getObj(self,name):
        if name in self.__objCatalog:
            return self.__objCatalog[name]
        else:
            return None

    def getObjNameList(self):
        return self.__objNameList
        
    def append(self,obj):
        """ Add the input object to the current object catalog. An existing object
            of the same name will be overwritten.
        """
        if obj.getName() is not None:
            if obj.getName() not in self.__objCatalog:
                # self.__objNameList is keeping track of object order here -- 
                self.__objNameList.append(obj.getName())
            self.__objCatalog[obj.getName()]=obj

    def replace(self,obj):
        """ Replace an existing object with the input object
        """
        if ((obj.getName() is not None) and (obj.getName() in self.__objCatalog) ):
            self.__objCatalog[obj.getName()]=obj
        

    def printIt(self,fh=sys.stdout,type="brief"):
        fh.write("+ %s container: %30s contains %4d categories\n" %
                 (self.getType(),self.getName(),len(self.__objNameList)))
        for nm in self.__objNameList:
            fh.write("--------------------------------------------\n")
            fh.write("Data category: %s\n" % nm)
            if type == 'brief':
                self.__objCatalog[nm].printIt(fh)
            else:
                self.__objCatalog[nm].dumpIt(fh)


    def rename(self,curName,newName):
        """ Change the name of an object in place -
        """
        try:
            i=self.__objNameList.index(curName)
            self.__objNameList[i]=newName
            self.__objCatalog[newName]=self.__objCatalog[curName]
            self.__objCatalog[newName].setName(newName)
            return True
        except:
            return False

    def remove(self,curName):
        """ Revmove object by name.  Return True on success or False otherwise.
        """
        try:
            if curName in self.__objCatalog:
                del self.__objCatalog[curName]
                i=self.__objNameList.index(curName)
                del self.__objNameList[i]
                return True
            else:
                return False
        except:
            pass

        return False


class DefinitionContainer(ContainerBase):
    def __init__(self,name):
        super(DefinitionContainer,self).__init__(name)        
        self.setType('definition')

    def isCategory(self):
        if self.exists('category'):
            return True
        return False

    def isAttribute(self):
        if self.exists('item'):
            return True
        return False
    

    def printIt(self,fh=sys.stdout,type="brief"):
        fh.write("Definition container: %30s contains %4d categories\n" %
                 (self.getName(),len(self.getObjNameList())))
        if self.isCategory():
            fh.write("Definition type: category\n")
        elif self.isAttribute():
            fh.write("Definition type: item\n")
        else:
            fh.write("Definition type: undefined\n")            

        for nm in self.getObjNameList():
            fh.write("--------------------------------------------\n")
            fh.write("Definition category: %s\n" % nm)
            if type == 'brief':
                self.getObj(nm).printIt(fh)
            else:
                self.getObj(nm).dumpId(fh)



class DataContainer(ContainerBase):
    ''' Container class for DataCategory objects.
    '''
    def __init__(self,name):
        super(DataContainer,self).__init__(name)
        self.setType('data')
        self.__globalFlag=False
        
    def invokeDataBlockMethod(self,type,method,db):
        self.__currentRow = 1
        exec(method.getInline())

    def setGlobal(self):
        self.__globalFlag=True

    def getGlobal(self):
        return  self.__globalFlag

        

class DataCategoryBase(object):
    """ Base object definition for a data category -
    """
    def __init__(self,name,attributeNameList=None,rowList=None):
        self._name = name
        #
        if rowList is not None:
            self._rowList=rowList
        else:
            self._rowList=[]

        if attributeNameList is not None:
            self._attributeNameList=attributeNameList
        else:
            self._attributeNameList=[]
        #
        # Derived class data -
        #
        self._catalog={}
        self._numAttributes=0
        #
        self.__setup()

    def __setup(self):
        self._numAttributes = len(self._attributeNameList)
        self._catalog={}
        for attributeName in self._attributeNameList:
            attributeNameLC  = attributeName.lower()
            self._catalog[attributeNameLC] = attributeName
    #
    def setRowList(self,rowList):
        self._rowList=rowList

    def setAttributeNameList(self,attributeNameList):
        self._attributeNameList=attributeNameList
        self.__setup()

    def setName(self,name):
        self._name=name

    def get(self):
        return (self._name,self._attributeNameList,self._rowList)

    
        
class DataCategory(DataCategoryBase):
    """  Methods for creating, accessing, and formatting PDBx cif data categories.  
    """
    def __init__(self,name,attributeNameList=None,rowList=None):
        super(DataCategory,self).__init__(name,attributeNameList,rowList)        
        #
        self.__lfh = sys.stdout
        
        self.__currentRowIndex=0
        self.__currentAttribute=None
        #
        self.__avoidEmbeddedQuoting=False
        #
        # --------------------------------------------------------------------
        # any whitespace 
        self.__wsRe=re.compile(r"\s")
        self.__wsAndQuotesRe=re.compile(r"[\s'\"]")        
        # any newline or carriage control
        self.__nlRe=re.compile(r"[\n\r]")
        #
        # single quote 
        self.__sqRe=re.compile(r"[']")
        #
        self.__sqWsRe=re.compile(r"('\s)|(\s')")
        
        # double quote 
        self.__dqRe=re.compile(r'["]')
        self.__dqWsRe=re.compile(r'("\s)|(\s")')
        #
        self.__intRe=re.compile(r'^[0-9]+$')
        self.__floatRe=re.compile(r'^-?(([0-9]+)[.]?|([0-9]*[.][0-9]+))([(][0-9]+[)])?([eE][+-]?[0-9]+)?$')
        #
        self.__dataTypeList=['DT_NULL_VALUE','DT_INTEGER','DT_FLOAT','DT_UNQUOTED_STRING','DT_ITEM_NAME',
                             'DT_DOUBLE_QUOTED_STRING','DT_SINGLE_QUOTED_STRING','DT_MULTI_LINE_STRING']
        self.__formatTypeList=['FT_NULL_VALUE','FT_NUMBER','FT_NUMBER','FT_UNQUOTED_STRING',
                               'FT_QUOTED_STRING','FT_QUOTED_STRING','FT_QUOTED_STRING','FT_MULTI_LINE_STRING']
        #


    def __getitem__(self, x):
        """  Implements list-type functionality - 
             Implements op[x] for some special cases -
                x=integer - returns the row in category (normal list behavior)
                x=string  - returns the value of attribute 'x' in first row.
        """
        if isinstance(x, int):
            #return self._rowList.__getitem__(x)
            return self._rowList[x]

        elif isinstance(x, str):
            try:
                #return self._rowList[0][x]
                ii=self.getAttributeIndex(x)
                return self._rowList[0][ii]
            except (IndexError, KeyError):
                raise KeyError
        raise TypeError(x)
        
    
    def getCurrentAttribute(self):
        return self.__currentAttribute


    def getRowIndex(self):
        return self.__currentRowIndex

    def getRowList(self):
        return self._rowList

    def getRowCount(self):
        return (len(self._rowList))

    def getRow(self,index):
        try:
            return self._rowList[index]
        except:
            return []

    def removeRow(self,index):
        try:
            if ((index >= 0) and (index < len(self._rowList))):
                del self._rowList[index] 
                if self.__currentRowIndex >= len(self._rowList):
                    self.__currentRowIndex = len(self._rowList) -1
                return True
            else:
                pass
        except:
            pass

        return False

    def getFullRow(self,index):
        """ Return a full row based on the length of the the attribute list.
        """
        try:
            if (len(self._rowList[index]) < self._numAttributes):
                for ii in range( self._numAttributes-len(self._rowList[index])):
                    self._rowList[index].append('?')
            return self._rowList[index]
        except:
            return ['?' for ii in range(self._numAttributes)]

    def getName(self):
        return self._name
    
    def getAttributeList(self):
        return self._attributeNameList

    def getAttributeCount(self):
        return len(self._attributeNameList)

    def getAttributeListWithOrder(self):
        oL=[]
        for ii,att in enumerate(self._attributeNameList):
            oL.append((att,ii))
        return oL

    def getAttributeIndex(self,attributeName):
        try:
            return self._attributeNameList.index(attributeName)
        except:
            return -1

    def hasAttribute(self,attributeName):
        return attributeName in self._attributeNameList
    
    def getIndex(self,attributeName):
        try:
            return self._attributeNameList.index(attributeName)
        except:
            return -1

    def getItemNameList(self):
        itemNameList=[]
        for att in self._attributeNameList:
            itemNameList.append("_"+self._name+"."+att)
        return itemNameList
    
    def append(self,row):
        #self.__lfh.write("PdbxContainer(append) category %s row %r\n"  % (self._name,row))
        self._rowList.append(row)

    def appendAttribute(self,attributeName):
        attributeNameLC  = attributeName.lower()
        if attributeNameLC  in self._catalog:
            i = self._attributeNameList.index(self._catalog[attributeNameLC])
            self._attributeNameList[i] = attributeName
            self._catalog[attributeNameLC] = attributeName
            #self.__lfh.write("Appending existing attribute %s\n" % attributeName)
        else:
            #self.__lfh.write("Appending existing attribute %s\n" % attributeName)            
            self._attributeNameList.append(attributeName)
            self._catalog[attributeNameLC] = attributeName
            #
        self._numAttributes = len(self._attributeNameList)


    def appendAttributeExtendRows(self,attributeName):
        attributeNameLC  = attributeName.lower()
        if attributeNameLC  in self._catalog:
            i = self._attributeNameList.index(self._catalog[attributeNameLC])
            self._attributeNameList[i] = attributeName
            self._catalog[attributeNameLC] = attributeName
            self.__lfh.write("Appending existing attribute %s\n" % attributeName)
        else:
            self._attributeNameList.append(attributeName)
            self._catalog[attributeNameLC] = attributeName
            # add a placeholder to any existing rows for the new attribute.
            if (len(self._rowList) > 0):
                for row in self._rowList:
                    row.append("?")
            #
        self._numAttributes = len(self._attributeNameList)

        

    def getValue(self,attributeName=None,rowIndex=None):
        if attributeName is None:
            attribute = self.__currentAttribute
        else:
            attribute = attributeName
        if rowIndex is None:
            rowI = self.__currentRowIndex
        else:
            rowI =rowIndex
            
        if isinstance(attribute, str) and isinstance(rowI,int):
            try:
                return self._rowList[rowI][self._attributeNameList.index(attribute)]
            except (IndexError):
                raise IndexError        
        raise IndexError(attribute)

    def setValue(self,value,attributeName=None,rowIndex=None):
        if attributeName is None:
            attribute=self.__currentAttribute
        else:
            attribute=attributeName

        if rowIndex is None:
            rowI = self.__currentRowIndex
        else:
            rowI = rowIndex

        if isinstance(attribute, str) and isinstance(rowI,int):
            try:
                # if row index is out of range - add the rows -
                for ii in range(rowI+1 - len(self._rowList)):
                    self._rowList.append(self.__emptyRow())                
                # self._rowList[rowI][attribute]=value
                ll=len(self._rowList[rowI])
                ind=self._attributeNameList.index(attribute)

                # extend the list if needed - 
                if ( ind >= ll):
                    self._rowList[rowI].extend([None for ii in xrange(2*ind -ll)])    
                self._rowList[rowI][ind]=value                    
            except (IndexError):
                self.__lfh.write("DataCategory(setvalue) index error category %s attribute %s index %d value %r\n" %
                                 (self._name,attribute,rowI,value))
                traceback.print_exc(file=self.__lfh)                                
                #raise IndexError
            except (ValueError):
                self.__lfh.write("DataCategory(setvalue) value error category %s attribute %s index %d value %r\n" %
                                 (self._name,attribute,rowI,value))                
                traceback.print_exc(file=self.__lfh)                
                #raise ValueError

    def __emptyRow(self):
        return [None for ii in range(len(self._attributeNameList))]
        
    def replaceValue(self,oldValue,newValue,attributeName):
        numReplace=0
        if attributeName not in self._attributeNameList:
            return numReplace
        ind=self._attributeNameList.index(attributeName)
        for row in self._rowList:
            if row[ind] == oldValue:
                row[ind]=newValue
                numReplace += 1
        return numReplace

    def replaceSubstring(self,oldValue,newValue,attributeName):
        ok=False
        if attributeName not in self._attributeNameList:            
            return ok
        ind=self._attributeNameList.index(attributeName)
        for row in self._rowList:
            val=row[ind]
            row[ind]=val.replace(oldValue,newValue)
            if val != row[ind]:
                ok=True
        return ok
        
    def invokeAttributeMethod(self,attributeName,type,method,db):
        self.__currentRowIndex = 0
        self.__currentAttribute=attributeName
        self.appendAttribute(attributeName)
        currentRowIndex=self.__currentRowIndex
        #
        ind=self._attributeNameList.index(attributeName)
        if len(self._rowList) == 0:
            row=[None for ii in xrange(len(self._attributeNameList)*2)]
            row[ind]=None
            self._rowList.append(row)
            
        for row in self._rowList:
            ll = len(row)
            if (ind >= ll):
                row.extend([None for ii in xrange(2*ind-ll)])
                row[ind]=None
            exec(method.getInline())
            self.__currentRowIndex+=1
            currentRowIndex=self.__currentRowIndex

    def invokeCategoryMethod(self,type,method,db):
        self.__currentRowIndex = 0
        exec(method.getInline())

    def getAttributeLengthMaximumList(self):
        mList=[0 for i in len(self._attributeNameList)]
        for row in self._rowList:
            for indx,val in enumerate(row):
                mList[indx] = max(mList[indx],len(val))
        return mList
            
    def renameAttribute(self,curAttributeName,newAttributeName):
        """ Change the name of an attribute in place -
        """
        try:
            i=self._attributeNameList.index(curAttributeName)
            self._attributeNameList[i]=newAttributeName
            del self._catalog[curAttributeName.lower()]            
            self._catalog[newAttributeName.lower()]=newAttributeName
            return True
        except:
            return False
        
    def printIt(self,fh=sys.stdout):
        fh.write("--------------------------------------------\n")
        fh.write("  Category: %s attribute list length: %d\n" %
                 (self._name,len(self._attributeNameList)))
        for at in self._attributeNameList:
            fh.write("  Category: %s attribute: %s\n" % (self._name,at))
            
        fh.write("  Row value list length: %d\n" % len(self._rowList))
        #
        for row in self._rowList[:2]:
            #
            if len(row) == len(self._attributeNameList):
                for ii,v in enumerate(row):
                    fh.write("        %30s: %s ...\n" % (self._attributeNameList[ii],str(v)[:30]))
            else:
                fh.write("+WARNING - %s data length %d attribute name length %s mismatched\n" %
                         (self._name,len(row),len(self._attributeNameList)))

    def dumpIt(self,fh=sys.stdout):
        fh.write("--------------------------------------------\n")
        fh.write("  Category: %s attribute list length: %d\n" %
                 (self._name,len(self._attributeNameList)))
        for at in self._attributeNameList:
            fh.write("  Category: %s attribute: %s\n" % (self._name,at))
            
        fh.write("  Value list length: %d\n" % len(self._rowList))
        for row in self._rowList:
            for ii,v in enumerate(row):
                fh.write("        %30s: %s\n" % (self._attributeNameList[ii],v))


    def __formatPdbx(self, inp):
        """ Format input data following PDBx quoting rules - 
        """
        try:
            if (inp is None):
                return ("?",'DT_NULL_VALUE')

            # pure numerical values are returned as unquoted strings
            if (isinstance(inp,int) or self.__intRe.search(str(inp))):
                return ( [str(inp)],'DT_INTEGER')

            if (isinstance(inp,float) or self.__floatRe.search(str(inp))):
                return ([str(inp)],'DT_FLOAT')

            # null value handling -

            if (inp == "." or inp == "?"):
                return ([inp],'DT_NULL_VALUE')

            if (inp == ""):
                return (["."],'DT_NULL_VALUE')

            # Contains white space or quotes ?
            if not self.__wsAndQuotesRe.search(inp):
                if inp.startswith("_"):
                    return (self.__doubleQuotedList(inp),'DT_ITEM_NAME')
                else:
                    return ([str(inp)],'DT_UNQUOTED_STRING')
            else:
                if self.__nlRe.search(inp):
                    return (self.__semiColonQuotedList(inp),'DT_MULTI_LINE_STRING')
                else:
                    if (self.__avoidEmbeddedQuoting):                    
                        # change priority to choose double quoting where possible.
                        if not self.__dqRe.search(inp) and not self.__sqWsRe.search(inp):
                            return (self.__doubleQuotedList(inp),'DT_DOUBLE_QUOTED_STRING')                                        
                        elif not self.__sqRe.search(inp) and not self.__dqWsRe.search(inp):
                            return (self.__singleQuotedList(inp),'DT_SINGLE_QUOTED_STRING')                                        
                        else:
                            return (self.__semiColonQuotedList(inp),'DT_MULTI_LINE_STRING')
                    else:
                        # change priority to choose double quoting where possible.
                        if not self.__dqRe.search(inp):
                            return (self.__doubleQuotedList(inp),'DT_DOUBLE_QUOTED_STRING')                                        
                        elif not self.__sqRe.search(inp):
                            return (self.__singleQuotedList(inp),'DT_SINGLE_QUOTED_STRING')                                        
                        else:
                            return (self.__semiColonQuotedList(inp),'DT_MULTI_LINE_STRING')
                        
                        
        except:
            traceback.print_exc(file=self.__lfh)                            

    def __dataTypePdbx(self, inp):
        """ Detect the PDBx data type - 
        """
        if (inp is None):
            return ('DT_NULL_VALUE')
        
        # pure numerical values are returned as unquoted strings
        if isinstance(inp,int) or self.__intRe.search(str(inp)):
            return ('DT_INTEGER')

        if isinstance(inp,float) or self.__floatRe.search(str(inp)):
            return ('DT_FLOAT')

        # null value handling -

        if (inp == "." or inp == "?"):
            return ('DT_NULL_VALUE')

        if (inp == ""):
            return ('DT_NULL_VALUE')

        # Contains white space or quotes ?
        if not self.__wsAndQuotesRe.search(inp):
            if inp.startswith("_"):
                return ('DT_ITEM_NAME')
            else:
                return ('DT_UNQUOTED_STRING')
        else:
            if self.__nlRe.search(inp):
                return ('DT_MULTI_LINE_STRING')
            else:
                if (self.__avoidEmbeddedQuoting):
                    if not self.__sqRe.search(inp) and not self.__dqWsRe.search(inp):
                        return ('DT_DOUBLE_QUOTED_STRING')                                        
                    elif not self.__dqRe.search(inp) and not self.__sqWsRe.search(inp):
                        return ('DT_SINGLE_QUOTED_STRING')                                        
                    else:
                        return ('DT_MULTI_LINE_STRING')
                else:
                    if not self.__sqRe.search(inp):
                        return ('DT_DOUBLE_QUOTED_STRING')                                        
                    elif not self.__dqRe.search(inp):
                        return ('DT_SINGLE_QUOTED_STRING')                                        
                    else:
                        return ('DT_MULTI_LINE_STRING')                    

    def __singleQuotedList(self,inp):
        l=[]
        l.append("'")
        l.append(inp)
        l.append("'")        
        return(l)

    def __doubleQuotedList(self,inp):
        l=[]
        l.append('"')
        l.append(inp)
        l.append('"')        
        return(l)
    
    def __semiColonQuotedList(self,inp):
        l=[]
        l.append("\n")            
        if inp[-1] == '\n':
            l.append(";")
            l.append(inp)
            l.append(";")
            l.append("\n")                        
        else:
            l.append(";")
            l.append(inp)
            l.append("\n")                        
            l.append(";")
            l.append("\n")                                    

        return(l)

    def getValueFormatted(self,attributeName=None,rowIndex=None):
        if attributeName is None:
            attribute=self.__currentAttribute
        else:
            attribute=attributeName

        if rowIndex is None:
            rowI = self.__currentRowIndex
        else:
            rowI = rowIndex
            
        if isinstance(attribute, str) and isinstance(rowI,int):
            try:
                list,type=self.__formatPdbx(self._rowList[rowI][self._attributeNameList.index(attribute)])
                return "".join(list)
            except (IndexError):
                self.__lfh.write("attributeName %s rowI %r rowdata %r\n" % (attributeName,rowI,self._rowList[rowI]))
                raise IndexError        
        raise TypeError(attribute)
        

    def getValueFormattedByIndex(self,attributeIndex,rowIndex):
        try:
            list,type=self.__formatPdbx(self._rowList[rowIndex][attributeIndex])
            return "".join(list)
        except (IndexError):
            raise IndexError        

    def getAttributeValueMaxLengthList(self,steps=1):
        mList=[0 for i in range(len(self._attributeNameList))]
        for row in self._rowList[::steps]:
            for indx in range(len(self._attributeNameList)):
                val=row[indx]                
                mList[indx] = max(mList[indx],len(str(val)))
        return mList

    def getFormatTypeList(self,steps=1):
        try:
            curDataTypeList=['DT_NULL_VALUE' for i in range(len(self._attributeNameList))]
            for row in self._rowList[::steps]:
                for indx in range(len(self._attributeNameList)):
                    val=row[indx]
                    # print "index ",indx," val ",val
                    dType=self.__dataTypePdbx(val)
                    dIndx=self.__dataTypeList.index(dType)
                    # print "d type", dType, " d type index ",dIndx
                
                    cType=curDataTypeList[indx]
                    cIndx=self.__dataTypeList.index(cType)           
                    cIndx= max(cIndx,dIndx)
                    curDataTypeList[indx]=self.__dataTypeList[cIndx]

            # Map the format types to the data types
            curFormatTypeList=[]
            for dt in curDataTypeList:
                ii=self.__dataTypeList.index(dt)
                curFormatTypeList.append(self.__formatTypeList[ii])
        except:
            self.__lfh.write("PdbxDataCategory(getFormatTypeList) ++Index error at index %d in row %r\n" % (indx,row))

        return curFormatTypeList,curDataTypeList

    def getFormatTypeListX(self):
        curDataTypeList=['DT_NULL_VALUE' for i in range(len(self._attributeNameList))]
        for row in self._rowList:
            for indx in range(len(self._attributeNameList)):
                val=row[indx]
                #print "index ",indx," val ",val
                dType=self.__dataTypePdbx(val)
                dIndx=self.__dataTypeList.index(dType)
                #print "d type", dType, " d type index ",dIndx
                
                cType=curDataTypeList[indx]
                cIndx=self.__dataTypeList.index(cType)           
                cIndx= max(cIndx,dIndx)
                curDataTypeList[indx]=self.__dataTypeList[cIndx]

        # Map the format types to the data types
        curFormatTypeList=[]
        for dt in curDataTypeList:
            ii=self.__dataTypeList.index(dt)
            curFormatTypeList.append(self.__formatTypeList[ii])
        return curFormatTypeList,curDataTypeList
    

