import sys
import unittest
import openmm as mm
from openmm import app
from openmm import unit


class TestOpcWater(unittest.TestCase):

    def test_OpcEnergy(self):
        pdb = app.PDBFile('../openmm/app/data/opcbox.pdb')
        topology, positions = pdb.topology, pdb.positions
        self.assertEqual(len(positions), 864)
        forcefield = app.ForceField('opc.xml')
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=0.7*unit.nanometer,
            constraints=app.HBonds,
            rigidWater=True,
        )

        box_vectors = pdb.topology.getPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(*box_vectors)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 2.0/unit.picoseconds, 2.0*unit.femtoseconds)
        simulation = app.Simulation(topology, system, integrator)
        context = simulation.context
        context.setPositions(positions)
        context.setPeriodicBoxVectors(*box_vectors)

        energy_amber = -2647.6233 # kcal/mol
        energy_tolerance = 1.0

        state = context.getState(getEnergy=True)
        energy1 = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        # -2647.2222697324237
        self.assertTrue(abs(energy1 - energy_amber) < energy_tolerance)

        context.applyConstraints(1e-12)
        state = context.getState(getEnergy=True)
        state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        energy2 = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        # -2647.441600693312
        self.assertTrue(abs(energy1 - energy_amber) < energy_tolerance)
        self.assertTrue(abs(energy1 - energy2) < energy_tolerance)


    def test_Opc3Energy(self):
        pdb = app.PDBFile('../openmm/app/data/opc3box.pdb')
        topology, positions = pdb.topology, pdb.positions
        self.assertEqual(len(positions), 648)
        forcefield = app.ForceField('opc3.xml')
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=0.7*unit.nanometer,
            constraints=app.HBonds,
            rigidWater=True,
        )

        box_vectors = pdb.topology.getPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(*box_vectors)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 2.0/unit.picoseconds, 2.0*unit.femtoseconds)
        simulation = app.Simulation(topology, system, integrator)
        context = simulation.context
        context.setPositions(positions)
        context.setPeriodicBoxVectors(*box_vectors)

        energy_amber = -2532.1414 # kcal/mol
        energy_tolerance = 1.0

        state = context.getState(getEnergy=True)
        energy1 = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        # -2532.4862082354407
        self.assertTrue(abs(energy1 - energy_amber) < energy_tolerance)


if __name__ == '__main__':
    unittest.main()
