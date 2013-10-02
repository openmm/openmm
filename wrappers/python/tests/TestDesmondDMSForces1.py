"""
This suite tests the forces produced by OpenMM when loading a system from
dms.createSystem() or loading the same system using ff.createSystem() with
the same topology/positions. This only works as a test for forcefields that
exist in both viparr (for building the dms file) and OpenMM's XML format.
"""

import os
import sys
import numpy as np
import unittest
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem


class TestDesmondDMSForces1(unittest.TestCase):
    def setUp(self):
        """Set up the tests by loading the input files."""

        # alanine dipeptide with explicit water
        path = os.path.join(os.path.dirname(__file__), 'systems/alanine-dipeptide-explicit-amber99SBILDN-tip3p.dms')
        self.dms = DesmondDMSFile(path)
        
    def testForces(self):
        system_dms = self.dms.createSystem()
        context_dms = Context(system_dms, VerletIntegrator(0.0))
        context_dms.setPositions(self.dms.positions)
        forces_dms = np.array(context_dms.getState(getForces=True).getForces(asNumpy=True))
    
        forcefield_mm = ForceField('amber99sbildn.xml', 'tip3p.xml')
        system_mm = forcefield_mm.createSystem(self.dms.topology, constraints=HBonds)
        context_mm = Context(system_mm, VerletIntegrator(0.0))
        context_mm.setPositions(self.dms.positions)
        forces_mm = np.array(context_mm.getState(getForces=True).getForces(asNumpy=True))
        
        # Agreement within 0.1 kJ/mol/nm
        np.testing.assert_array_almost_equal(forces_dms, forces_mm, decimal=1)
