__author__ = "Raul P. Pelaez"
import unittest
import tempfile
from openmm import app
import openmm as mm
from openmm import unit
from openmm.unit.unit_math import norm
from openmm.unit import Quantity, nanometers
import os
import gc
import numpy as np


class TestXTCReporter(unittest.TestCase):
    def setUp(self):
        self.pdb = app.PDBFile("systems/alanine-dipeptide-explicit.pdb")
        self.forcefield = app.ForceField("amber99sbildn.xml", "tip3p.xml")
        self.system = self.forcefield.createSystem(
            self.pdb.topology,
            nonbondedMethod=app.CutoffNonPeriodic,
            constraints=app.HBonds,
        )

    def test_OneStep(self):
        """Test all atoms are written for a single step"""

        with tempfile.TemporaryDirectory() as tempdir:
            filename_xtc = os.path.join(tempdir, "temptraj.xtc")
            simulation = app.Simulation(
                self.pdb.topology,
                self.system,
                mm.LangevinMiddleIntegrator(
                    300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
                ),
            )
            simulation.context.setPositions(self.pdb.positions)

            simulation.reporters.append(app.XTCReporter(filename_xtc, 1))
            simulation.step(1)
            # clear reporters to ensure PDBReporter calls writeFooter and file.close
            simulation.reporters.clear()
            checkxtc = app.XTCFile(filename_xtc)
            # check the positions are correct
            validPositions = (
                simulation.context.getState(getPositions=True)
                .getPositions(True)
                .value_in_unit(unit.nanometer)
            )
            coords, _, _, _ = checkxtc.read()
            self.assertEqual(coords[:, :, 0].shape, validPositions.shape)
            self.assertTrue(np.allclose(coords[:, :, 0], validPositions, atol=1e-3))


    def test_SeveralSteps(self):
        """Test that the reporter can write several steps"""

        with tempfile.TemporaryDirectory() as tempdir:
            filename_xtc = os.path.join(tempdir, "temptraj.xtc")
            simulation = app.Simulation(
                self.pdb.topology,
                self.system,
                mm.LangevinMiddleIntegrator(
                    300 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds
                ),
            )
            simulation.context.setPositions(self.pdb.positions)

            simulation.reporters.append(app.XTCReporter(filename_xtc, 1))
            simulation.step(10)
            # clear reporters to ensure PDBReporter calls writeFooter and file.close
            simulation.reporters.clear()
            checkxtc = app.XTCFile(filename_xtc)
            # check the positions are correct
            validPositions = (
                simulation.context.getState(getPositions=True)
                .getPositions(True)
                .value_in_unit(unit.nanometer)
            )
            coords, _, _, _ = checkxtc.read()
            self.assertEqual(coords[:, :, -1].shape, validPositions.shape)
            self.assertTrue(np.allclose(coords[:, :, -1], validPositions, atol=1e-3))


if __name__ == "__main__":
    unittest.main()
