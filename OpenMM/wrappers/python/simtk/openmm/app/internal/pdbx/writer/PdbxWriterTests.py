##
# File:    PdbxWriterTests.py
# Author:  jdw
# Date:    3-November-2009
# Version: 0.001
#
# Update:
#  5-Apr-2011 jdw   Using the double quote format preference
# 24-Oct-2012 jdw   Update path and examples.
##
"""
Test implementing PDBx/mmCIF write and formatting operations.

"""
from __future__ import absolute_import
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.01"



import sys, unittest, traceback
import sys, time, os, os.path, shutil

from simtk.openmm.app.internal.pdbx.reader.PdbxReader  import PdbxReader
from simtk.openmm.app.internal.pdbx.writer.PdbxWriter  import PdbxWriter
from simtk.openmm.app.internal.pdbx.reader.PdbxContainers import *

class PdbxWriterTests(unittest.TestCase):
    def setUp(self):
        self.lfh=sys.stderr
        self.verbose=False
        self.pathPdbxDataFile     ="../tests/1kip.cif"
        self.pathOutputFile       ="testOutputDataFile.cif"

    def tearDown(self):
        pass

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
            aCat.append((1,2,3,4,5,6,7))
            aCat.append((1,2,3,4,5,6,7))
            aCat.append((1,2,3,4,5,6,7))
            aCat.append((1,2,3,4,5,6,7))
            curContainer.append(aCat)
            myDataList.append(curContainer)
            pdbxW=PdbxWriter(ofh)
            pdbxW.write(myDataList)
            ofh.close()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testUpdateDataFile(self): 
        """Test case -  write data file 
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            # Create a initial data file --
            #
            myDataList=[]
            ofh = open("test-output-1.cif", "w")
            curContainer=DataContainer("myblock")
            aCat=DataCategory("pdbx_seqtool_mapping_ref")
            aCat.appendAttribute("ordinal")
            aCat.appendAttribute("entity_id")
            aCat.appendAttribute("auth_mon_id")
            aCat.appendAttribute("auth_mon_num")
            aCat.appendAttribute("pdb_chain_id")
            aCat.appendAttribute("ref_mon_id")
            aCat.appendAttribute("ref_mon_num")                        
            aCat.append((1,2,3,4,5,6,7))
            aCat.append((1,2,3,4,5,6,7))
            aCat.append((1,2,3,4,5,6,7))
            aCat.append((1,2,3,4,5,6,7))
            curContainer.append(aCat)
            myDataList.append(curContainer)
            pdbxW=PdbxWriter(ofh)
            pdbxW.write(myDataList)
            ofh.close()
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
            
        except:
            traceback.print_exc(file=sys.stderr)
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
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testReadWriteDataFile(self): 
        """Test case -  data file read write test
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
            
            ofh = open(self.pathOutputFile, "w")
            pWr=PdbxWriter(ofh)
            pWr.write(myDataList)        
            ofh.close()
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

def suite():
    return unittest.makeSuite(PdbxWriterTests,'test')

if __name__ == '__main__':
    unittest.main()
