import unittest
import tempfile
from openmm import app
import openmm as mm
from openmm import unit
from random import random
import os
import struct


def _read_dcd_header(file):
    with open(file, "r+b") as f:
        f.seek(8, os.SEEK_SET)
        modelCount = struct.unpack("<i", f.read(4))[0]
        f.seek(20, os.SEEK_SET)
        currStep = struct.unpack("<i", f.read(4))[0]
        return modelCount, currStep

class TestDCDFile(unittest.TestCase):
    def test_dcd(self):
        """ Test the DCD file """
        fname = tempfile.mktemp(suffix='.dcd')
        pdbfile = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        natom = len(list(pdbfile.topology.atoms()))
        with open(fname, 'wb') as f:
            dcd = app.DCDFile(f, pdbfile.topology, 0.001)
            for i in range(5):
                dcd.writeModel([mm.Vec3(random(), random(), random()) for j in range(natom)]*unit.angstroms)
        os.remove(fname)
    
    def testLongTrajectory(self):
        """Test writing a trajectory that has more than 2^31 steps."""
        fname = tempfile.mktemp(suffix='.dcd')
        pdbfile = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        natom = len(list(pdbfile.topology.atoms()))
        with open(fname, 'wb') as f:
            dcd = app.DCDFile(f, pdbfile.topology, 0.001, interval=1000000000)
            for i in range(5):
                dcd.writeModel([mm.Vec3(random(), random(), random()) for j in range(natom)]*unit.angstroms)
        os.remove(fname)
    
    def testAppend(self):
        """Test appending to an existing trajectory."""
        fname = tempfile.mktemp(suffix='.dcd')
        pdb = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)
        
        # Create a simulation and write some frames to a DCD file.
        
        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        dcd = app.DCDReporter(fname, 2)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        self.assertEqual(5, dcd._dcd._modelCount)
        del simulation
        del dcd
        len1 = os.stat(fname).st_size
        modelCount, currStep = _read_dcd_header(fname)
        self.assertEqual(5, modelCount)
        self.assertEqual(10, currStep)

        # Create a new simulation and have it append some more frames.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        simulation.currentStep = 10
        dcd = app.DCDReporter(fname, 2, append=True)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        self.assertEqual(10, dcd._dcd._modelCount)
        len2 = os.stat(fname).st_size
        self.assertTrue(len2-len1 > 3*4*5*system.getNumParticles())
        del simulation
        del dcd
        modelCount, currStep = _read_dcd_header(fname)
        self.assertEqual(10, modelCount)
        self.assertEqual(20, currStep)
        os.remove(fname)

    def testAtomSubset(self):
        """Test writing a DCD file containing a subset of atoms"""
        fname = tempfile.mktemp(suffix='.dcd')
        pdb = app.PDBFile('systems/alanine-dipeptide-explicit.pdb')
        ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)

        # Create a simulation and write some frames to a DCD file.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        atomSubset = [atom.index for atom in next(pdb.topology.chains()).atoms()]
        dcd = app.DCDReporter(fname, 2, atomSubset=atomSubset)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        self.assertEqual(5, dcd._dcd._modelCount)
        del simulation
        del dcd
        len1 = os.stat(fname).st_size
        modelCount, currStep = _read_dcd_header(fname)
        self.assertEqual(5, modelCount)
        self.assertEqual(10, currStep)

        # Create a new simulation and have it append some more frames.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        simulation.currentStep = 10
        dcd = app.DCDReporter(fname, 2, append=True, atomSubset=atomSubset)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        self.assertEqual(10, dcd._dcd._modelCount)
        len2 = os.stat(fname).st_size
        self.assertTrue(len2-len1 > 3*4*5*len(atomSubset))
        del simulation
        del dcd
        modelCount, currStep = _read_dcd_header(fname)
        self.assertEqual(10, modelCount)
        self.assertEqual(20, currStep)
        os.remove(fname)

    def testAppendAtomCountMismatch(self):
        """Test that appending to a DCD file with a different number of atoms raises an error."""
        fname = tempfile.mktemp(suffix='.dcd')
        pdb = app.PDBFile('systems/alanine-dipeptide-explicit.pdb')
        ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)

        # Create a simulation and write some frames to a DCD file.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        atomSubset1 = [atom.index for chain, _ in zip(pdb.topology.chains(), range(1)) for atom in chain.atoms()]
        atomSubset2 = [atom.index for chain, _ in zip(pdb.topology.chains(), range(2)) for atom in chain.atoms()]
        dcd = app.DCDReporter(fname, 2, atomSubset=atomSubset1)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        del simulation
        del dcd

        # Create a new simulation and have it append some more frames.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        simulation.currentStep = 10
        dcd = app.DCDReporter(fname, 2, append=True, atomSubset=atomSubset2)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        with self.assertRaises(ValueError):
            simulation.step(10)

    def testAppendLongCommentBlock(self):
        """Test appending to an existing trajectory with a long comment block."""
        fname = tempfile.mktemp(suffix='.dcd')
        pdb = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology)

        # Create a simulation and write some frames to a DCD file.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        dcd = app.DCDReporter(fname, 2)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        self.assertEqual(5, dcd._dcd._modelCount)
        del simulation
        del dcd
        len1 = os.stat(fname).st_size

        # Some software writes more than 2 80-byte "comment lines" to a DCD
        # file.  Modify the DCD to simulate this and ensure we can append.

        commentLines = 10
        with open(fname, "rb") as dcdFile:
            dcdHeader = dcdFile.read(92)
            dcdCommentsLength = struct.unpack("<i", dcdFile.read(4))[0]
            dcdFile.read(dcdCommentsLength + 4)
            dcdContents = dcdFile.read()
        with open(fname, "wb") as dcdFile:
            dcdFile.write(dcdHeader)
            dcdNewCommentsLength = 4 + 80 * commentLines
            dcdFile.write(struct.pack("<2i", dcdNewCommentsLength, commentLines))
            dcdFile.write(bytes(80 * commentLines))
            dcdFile.write(struct.pack("<i", dcdNewCommentsLength))
            dcdFile.write(dcdContents)

        # Create a new simulation and have it append some more frames.

        integrator = mm.VerletIntegrator(0.001*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator, mm.Platform.getPlatform('Reference'))
        simulation.currentStep = 10
        dcd = app.DCDReporter(fname, 2, append=True)
        simulation.reporters.append(dcd)
        simulation.context.setPositions(pdb.positions)
        simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
        simulation.step(10)
        self.assertEqual(10, dcd._dcd._modelCount)
        len2 = os.stat(fname).st_size
        self.assertEqual(len2 - len1, dcdNewCommentsLength - dcdCommentsLength + 5 * 3 * 4 * (system.getNumParticles() + 2))
        del simulation
        del dcd
        modelCount, currStep = _read_dcd_header(fname)
        self.assertEqual(10, modelCount)
        self.assertEqual(20, currStep)
        os.remove(fname)


if __name__ == '__main__':
    unittest.main()
