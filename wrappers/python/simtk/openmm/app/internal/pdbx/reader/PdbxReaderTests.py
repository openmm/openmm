##
# File:    PdbxReaderTests.py
# Author:  jdw
# Date:    9-Jan-2012
# Version: 0.001
#
# Update:
#  27-Sep-2012  jdw add test case for reading PDBx structure factor file 
#
##
"""
Test cases for reading PDBx/mmCIF data files PdbxReader class -

"""
from __future__ import absolute_import
import sys, unittest, traceback
import sys, time, os, os.path, shutil

from simtk.openmm.app.internal.pdbx.reader.PdbxReader import PdbxReader
from simtk.openmm.app.internal.pdbx.reader.PdbxContainers import *

class PdbxReaderTests(unittest.TestCase):
    def setUp(self):
        self.lfh=sys.stderr
        self.verbose=False
        self.pathPdbxDataFile     ="../tests/1kip.cif"
        self.pathBigPdbxDataFile  ="../tests/1ffk.cif"
        self.pathSFDataFile       ="../tests/1kip-sf.cif"

    def tearDown(self):
        pass

    def testReadSmallDataFile(self): 
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

    def testReadBigDataFile(self): 
        """Test case -  read data file 
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            myDataList=[]
            ifh = open(self.pathBigPdbxDataFile, "r")
            pRd=PdbxReader(ifh)
            pRd.read(myDataList)
            ifh.close()            
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()

    def testReadSFDataFile(self): 
        """Test case -  read PDB structure factor data  file and compute statistics on f/sig(f).
        """
        self.lfh.write("\nStarting %s %s\n" % (self.__class__.__name__,
                                               sys._getframe().f_code.co_name))
        try:
            #
            myContainerList=[]
            ifh = open(self.pathSFDataFile, "r")
            pRd=PdbxReader(ifh)
            pRd.read(myContainerList)
            c0=myContainerList[0]
            #
            catObj=c0.getObj("refln")
            if catObj is None:
                return false
        
            nRows=catObj.getRowCount()
            #
            # Get column name index.
            #
            itDict={}
            itNameList=catObj.getItemNameList()
            for idxIt,itName in enumerate(itNameList):
                itDict[str(itName).lower()]=idxIt
                #
            idf=itDict['_refln.f_meas_au']
            idsigf=itDict['_refln.f_meas_sigma_au']
            minR=100
            maxR=-1
            sumR=0
            icount=0
            for row in  catObj.getRowList():
                try:
                    f=float(row[idf])
                    sigf=float(row[idsigf])
                    ratio=sigf/f
                    #self.lfh.write(" %f %f %f\n" % (f,sigf,ratio))
                    maxR=max(maxR,ratio)
                    minR=min(minR,ratio)
                    sumR+=ratio
                    icount+=1
                except:
                    continue
            
            ifh.close()
            self.lfh.write("f/sig(f) min %f max %f avg %f count %d\n" % (minR, maxR, sumR/icount,icount))
        except:
            traceback.print_exc(file=sys.stderr)
            self.fail()


def simpleSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(PdbxReaderTests("testReadBigDataFile"))
    suiteSelect.addTest(PdbxReaderTests("testReadSmallDataFile"))    
    suiteSelect.addTest(PdbxReaderTests("testReadSFDataFile"))
    return suiteSelect

if __name__ == '__main__':
    mySuite=simpleSuite()    
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
