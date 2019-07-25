##
# File:    PdbxReadWriteTests.py
# Author:  jdw
# Date:    9-Oct-2011
# Version: 0.001
#
# Updated:
#         24-Oct-2012 jdw update path details and reorganize.
#
##
"""  Various tests caess for PDBx/mmCIF data file and dictionary reader and writer. 
"""
from __future__ import absolute_import

__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.01"

import sys, unittest, traceback
import sys, time, os, os.path, shutil

from simtk.openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from simtk.openmm.app.internal.pdbx.writer.PdbxWriter import PdbxWriter
from simtk.openmm.app.internal.pdbx.reader.PdbxContainers import *


class PdbxReadWriteTests(unittest.TestCase):
    def setUp(self):
        self.lfh=sys.stdout
        self.verbose=False
        self.pathPdbxDataFile     = "../tests/1kip.cif"
        self.pathOutputFile       = "testOutputDataFile.cif"

    def tearDown(self):
        pass


    def testSimpleInitialization(self):
        """Test case -  Simple initialization of a data category and data block
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            fn="test-simple.cif"
            attributeNameList=['aOne','aTwo','aThree','aFour','aFive','aSix','aSeven','aEight','aNine','aTen']
            rowList=[[1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10],
                     [1,2,3,4,5,6,7,8,9,10] 
                     ]
            nameCat='myCategory'
            #
            #
            curContainer=DataContainer("myblock")
            aCat=DataCategory(nameCat,attributeNameList,rowList)
            aCat.printIt()
            curContainer.append(aCat)
            curContainer.printIt()
            #
            myContainerList=[]
            myContainerList.append(curContainer)
            ofh = open(fn, "w")        
            pdbxW=PdbxWriter(ofh)
            pdbxW.write(myContainerList)
            ofh.close()

            myContainerList=[]            
            ifh = open(fn, "r")
            pRd=PdbxReader(ifh)
            pRd.read(myContainerList)
            ifh.close()
            for container in myContainerList:
                for objName in container.getObjNameList():
                    name,aList,rList=container.getObj(objName).get()
                    self.lfh.write("Recovered data category  %s\n" % name)
                    self.lfh.write("Attribute list           %r\n" % repr(aList))
                    self.lfh.write("Row list                 %r\n" % repr(rList))                                        
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()
            
        
    def testWriteDataFile(self): 
        """Test case -  write data file 
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            myDataList=[]
            ofh = open("test-output.cif", "w")
            curContainer=DataContainer("myblock")
            aCat=DataCategory("pdbx_seqtool_mapping_ref")
            aCat.appendAttribute("ordinal")
            aCat.appendAttribute("entity_id")
            aCat.appendAttribute("auth_mon_id")
            aCat.appendAttribute("auth_mon_num")
            aCat.appendAttribute("pdb_chain_id")
            aCat.appendAttribute("ref_mon_id")
            aCat.appendAttribute("ref_mon_num")                        
            aCat.append([1,2,3,4,5,6,7])
            aCat.append([1,2,3,4,5,6,7])
            aCat.append([1,2,3,4,5,6,7])
            aCat.append([1,2,3,4,5,6,7])
            aCat.append([7,6,5,4,3,2,1])
            aCat.printIt()            
            curContainer.append(aCat)
            curContainer.printIt()
            #
            myDataList.append(curContainer)
            pdbxW=PdbxWriter(ofh)
            pdbxW.write(myDataList)
            ofh.close()
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def testUpdateDataFile(self): 
        """Test case -  update data file 
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            # Create a initial data file --
            #
            myDataList=[]

            curContainer=DataContainer("myblock")
            aCat=DataCategory("pdbx_seqtool_mapping_ref")
            aCat.appendAttribute("ordinal")
            aCat.appendAttribute("entity_id")
            aCat.appendAttribute("auth_mon_id")
            aCat.appendAttribute("auth_mon_num")
            aCat.appendAttribute("pdb_chain_id")
            aCat.appendAttribute("ref_mon_id")
            aCat.appendAttribute("ref_mon_num")                        
            aCat.append([9,2,3,4,5,6,7])
            aCat.append([10,2,3,4,5,6,7])
            aCat.append([11,2,3,4,5,6,7])
            aCat.append([12,2,3,4,5,6,7])
            
            #self.lfh.write("Assigned data category state-----------------\n")            
            #aCat.dumpIt(fh=self.lfh)

            curContainer.append(aCat)
            myDataList.append(curContainer)
            ofh = open("test-output-1.cif", "w")            
            pdbxW=PdbxWriter(ofh)
            pdbxW.write(myDataList)
            ofh.close()
            #
            #
            # Read and update the data -
            # 
            myDataList=[]
            ifh = open("test-output-1.cif", "r")
            pRd=PdbxReader(ifh)
            pRd.read(myDataList)
            ifh.close()
            #
            myBlock=myDataList[0]
            myBlock.printIt()
            myCat=myBlock.getObj('pdbx_seqtool_mapping_ref')
            myCat.printIt()
            for iRow in xrange(0,myCat.getRowCount()):
                myCat.setValue('some value', 'ref_mon_id',iRow)
                myCat.setValue(100, 'ref_mon_num',iRow)
            ofh = open("test-output-2.cif", "w")            
            pdbxW=PdbxWriter(ofh)
            pdbxW.write(myDataList)
            ofh.close()
            #
            
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def testReadDataFile(self): 
        """Test case -  read data file 
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            myDataList=[]
            ifh = open(self.pathPdbxDataFile, "r")
            pRd=PdbxReader(ifh)
            pRd.read(myDataList)
            ifh.close()            
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()

    def testReadWriteDataFile(self): 
        """Test case -  data file read write test
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            myDataList=[]            
            ifh = open(self.pathPdbxDataFile, "r")
            pRd=PdbxReader(ifh)
            pRd.read(myDataList)
            ifh.close()            
            
            ofh = open(self.pathOutputFile, "w")
            pWr=PdbxWriter(ofh)
            pWr.write(myDataList)        
            ofh.close()
        except:
            traceback.print_exc(file=self.lfh)
            self.fail()


def simpleSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxReadWriteTests("testSimpleInitialization"))
    suiteSelect.addTest(PdbxReadWriteTests("testUpdateDataFile"))        
    suiteSelect.addTest(PdbxReadWriteTests("testReadWriteDataFile"))    
    return suiteSelect


if __name__ == '__main__':
    #
    mySuite=simpleSuite()      
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #

