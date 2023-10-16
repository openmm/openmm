__author__ = "Raul P. Pelaez"
import sys
import os
import unittest
from openmm import unit
import tempfile
from openmm import app
from random import random
import openmm as mm
import numpy as np
from openmm.app.internal.xtc_utils import read_xtc

class TestXtcFile(unittest.TestCase):
    def test_xtc_triclinic(self):
        """Test the XTC file by writing a trajectory and reading it back. Using a triclinic box"""
        with tempfile.TemporaryDirectory() as temp:
            fname = os.path.join(temp, 'traj.xtc')
            pdbfile = app.PDBFile("systems/alanine-dipeptide-implicit.pdb")
            # Set some arbitrary size for the unit cell so that a box is included in the trajectory
            pdbfile.topology.setUnitCellDimensions([10, 10, 10])
            natom = len(list(pdbfile.topology.atoms()))
            nframes = 20
            xtc = app.XTCFile(fname, pdbfile.topology, 0.001)
            coords = []
            box = []
            for i in range(nframes):
                coords.append(
                    [mm.Vec3(random(), random(), random()) for j in range(natom)]
                    * unit.nanometers
                )
                box_i = (
                    mm.Vec3(random(), random(), random()) * unit.nanometers,
                    mm.Vec3(random(), random(), random()) * unit.nanometers,
                    mm.Vec3(random(), random(), random()) * unit.nanometers,
                )
                box.append(np.array([[vec.x, vec.y, vec.z] for vec in box_i]))
                xtc.writeModel(coords[i], periodicBoxVectors=box_i)
            # The  XTCFile class  does not  provide a  way to  read the
            # trajectory back, but the underlying XTC library does
            coords_read, box_read, time, step = read_xtc(fname.encode("utf-8"))
            self.assertEqual(coords_read.shape, (natom, 3, nframes))
            self.assertEqual(box_read.shape, (3, 3, nframes))
            self.assertEqual(len(time), nframes)
            self.assertEqual(len(step), nframes)
            coords = np.array(
                [c.value_in_unit(unit.nanometers) for c in coords]
            ).transpose(1, 2, 0)
            self.assertTrue(np.allclose(coords_read, coords, atol=1e-3))
            box = np.array(box).transpose(1, 2, 0)
            self.assertTrue(np.allclose(box_read, box, atol=1e-3))
            self.assertTrue(
                np.allclose(time, np.arange(0, nframes) * 0.001, atol=1e-5)
            )
            self.assertTrue(np.allclose(step, np.arange(0, nframes), atol=1e-5))

    def test_xtc_cubic(self):
        """Test the XTC file by writing a trajectory and reading it back. Using a cubic box"""
        with tempfile.TemporaryDirectory() as temp:
            fname = os.path.join(temp, 'traj.xtc')
            pdbfile = app.PDBFile("systems/alanine-dipeptide-implicit.pdb")
            # Set some arbitrary size for the unit cell so that a box is included in the trajectory
            pdbfile.topology.setUnitCellDimensions([10, 10, 10])
            natom = len(list(pdbfile.topology.atoms()))
            nframes=20
            xtc = app.XTCFile(fname, pdbfile.topology, 0.001)
            coords = []
            box = []
            for i in range(nframes):
                coords.append(
                    [mm.Vec3(random(), random(), random()) for j in range(natom)]
                    * unit.nanometers
                )
                box_i = mm.Vec3(random(), random(), random()) * unit.nanometers
                box.append(np.diag(box_i[:3].value_in_unit(unit.nanometers)))
                xtc.writeModel(coords[i], unitCellDimensions=box_i)
            #The  XTCFile class  does not  provide a  way to  read the
            #trajectory back, but the underlying XTC library does
            coords_read, box_read, time, step = read_xtc(fname.encode("utf-8"))
            self.assertEqual(coords_read.shape, (natom,3,nframes))
            self.assertEqual(box_read.shape, (3,3,nframes))
            self.assertEqual(len(time), nframes)
            self.assertEqual(len(step), nframes)
            coords = np.array(
                [c.value_in_unit(unit.nanometers) for c in coords]
            ).transpose(1, 2, 0)
            self.assertTrue(np.allclose(coords_read, coords, atol=1e-3))
            box = np.array(box).transpose(1, 2, 0)
            self.assertTrue(np.allclose(box_read, box, atol=1e-3))
            self.assertTrue(
                np.allclose(time, np.arange(0, nframes) * 0.001, atol=1e-5)
            )
            self.assertTrue(np.allclose(step, np.arange(0, nframes), atol=1e-5))

    def test_xtc_box_from_topology(self):
        """Test the XTC file by writing a trajectory and reading it back. Letting the box be set from the topology"""
        with tempfile.TemporaryDirectory() as temp:
            fname = os.path.join(temp, 'traj.xtc')
            pdbfile = app.PDBFile("systems/alanine-dipeptide-implicit.pdb")
            # Set some arbitrary size for the unit cell so that a box is included in the trajectory
            unitCell = mm.Vec3(random(), random(), random()) * unit.nanometers
            pdbfile.topology.setUnitCellDimensions(unitCell)
            natom = len(list(pdbfile.topology.atoms()))
            nframes = 20
            xtc = app.XTCFile(fname, pdbfile.topology, 0.001)
            coords = []
            for i in range(nframes):
                coords.append(
                    [mm.Vec3(random(), random(), random()) for j in range(natom)]
                    * unit.nanometers
                )
                xtc.writeModel(coords[i])
            # The  XTCFile class  does not  provide a  way to  read the
            # trajectory back, but the underlying XTC library does
            coords_read, box_read, time, step = read_xtc(fname.encode("utf-8"))
            self.assertEqual(coords_read.shape, (natom, 3, nframes))
            self.assertEqual(box_read.shape, (3, 3, nframes))
            self.assertEqual(len(time), nframes)
            self.assertEqual(len(step), nframes)
            coords = np.array(
                [c.value_in_unit(unit.nanometers) for c in coords]
            ).transpose(1, 2, 0)
            self.assertTrue(np.allclose(coords_read, coords, atol=1e-3))
            boxVectors = (
                mm.Vec3(unitCell[0].value_in_unit(unit.nanometers), 0, 0),
                mm.Vec3(0, unitCell[1].value_in_unit(unit.nanometers), 0),
                mm.Vec3(0, 0, unitCell[2].value_in_unit(unit.nanometers)),
            )
            boxVectors = np.array(
                [[vec.x, vec.y, vec.z] for vec in boxVectors], dtype=np.float32
            )
            box = np.array(np.tile(boxVectors, (nframes, 1, 1))).transpose(1, 2, 0)
            self.assertTrue(np.allclose(box_read, box, atol=1e-3))
            self.assertTrue(
                np.allclose(time, np.arange(0, nframes) * 0.001, atol=1e-5)
            )
            self.assertTrue(np.allclose(step, np.arange(0, nframes), atol=1e-5))

    def testLongTrajectory(self):
        """Test writing a trajectory that has more than 2^31 steps."""
        with tempfile.TemporaryDirectory() as temp:
            fname = os.path.join(temp, 'traj.xtc')
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
        with tempfile.TemporaryDirectory() as temp:
            fname = os.path.join(temp, 'traj.xtc')
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
