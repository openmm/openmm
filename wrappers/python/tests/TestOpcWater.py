import sys
import unittest
import openmm as mm
from openmm import app
from openmm import unit

class TestOpcWater(unittest.TestCase):

    def test_Energy(self):
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

        box_vectors = [
            [1.8864844000000003, 0.       , 0.       ],
            [0.       , 1.8478108, 0.       ],
            [0.       , 0.       , 1.9006241],
        ]
        system.setDefaultPeriodicBoxVectors(*box_vectors)
        integrator = mm.LangevinIntegrator(300*unit.kelvin, 2.0/unit.picoseconds, 2.0*unit.femtoseconds)
        simulation = app.Simulation(topology, system, integrator)
        context = simulation.context
        context.setPositions(positions)
        context.setPeriodicBoxVectors(*box_vectors)

        energy_amber = -2647.6233
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
        

if __name__ == '__main__':
    unittest.main()
