__author__ = "Raul P. Pelaez"
import sys
import unittest
from openmm import unit
import tempfile
from openmm import app
from random import random
import openmm as mm

if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO


class TestXtcFile(unittest.TestCase):
    def test_xtc(self):
        """Test the XTC file"""
        with tempfile.NamedTemporaryFile() as temp:
            fname = temp.name
            pdbfile = app.PDBFile("systems/alanine-dipeptide-implicit.pdb")
            natom = len(list(pdbfile.topology.atoms()))
            xtc = app.XTCFile(fname, pdbfile.topology, 0.001)
            for i in range(5):
                xtc.writeModel(
                    [mm.Vec3(random(), random(), random()) for j in range(natom)]
                    * unit.angstroms
                )

    def testLongTrajectory(self):
        """Test writing a trajectory that has more than 2^31 steps."""
        with tempfile.NamedTemporaryFile() as temp:
            fname = temp.name
            pdbfile = app.PDBFile("systems/alanine-dipeptide-implicit.pdb")
            natom = len(list(pdbfile.topology.atoms()))
            xtc = app.XTCFile(fname, pdbfile.topology, 0.001, interval=1000000000)
            for i in range(5):
                xtc.writeModel(
                    [mm.Vec3(random(), random(), random()) for j in range(natom)]
                    * unit.angstroms
                )

    def testAppend(self):
        """Test appending to an existing trajectory."""
        with tempfile.NamedTemporaryFile() as temp:
            fname = temp.name
            pdb = app.PDBFile("systems/alanine-dipeptide-implicit.pdb")
            ff = app.ForceField("amber99sb.xml", "tip3p.xml")
            system = ff.createSystem(pdb.topology)

            # Create a simulation and write some frames to a XTC file.

            integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
            simulation = app.Simulation(
                pdb.topology,
                system,
                integrator,
                mm.Platform.getPlatformByName("Reference"),
            )
            xtc = app.XTCReporter(fname, 2)
            simulation.reporters.append(xtc)
            simulation.context.setPositions(pdb.positions)
            simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
            simulation.step(10)
            self.assertEqual(5, xtc._xtc._modelCount)
            self.assertEqual(5, xtc._xtc._getNumFrames())
            del simulation
            del xtc

            # Create a new simulation and have it append some more frames.

            integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
            simulation = app.Simulation(
                pdb.topology,
                system,
                integrator,
                mm.Platform.getPlatformByName("Reference"),
            )
            xtc = app.XTCReporter(fname, 2, append=True)
            simulation.reporters.append(xtc)
            simulation.context.setPositions(pdb.positions)
            simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)
            simulation.step(10)
            self.assertEqual(10, xtc._xtc._modelCount)
            self.assertEqual(10, xtc._xtc._getNumFrames())
            del simulation
            del xtc


if __name__ == "__main__":
    unittest.main()
