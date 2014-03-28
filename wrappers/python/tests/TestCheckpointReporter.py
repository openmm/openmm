import unittest
import tempfile
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit


class TestCheckpointReporter(unittest.TestCase):
    def setUp(self):
        pdb = app.PDBFile('systems/alanine-dipeptide-implicit.pdb')
        forcefield = app.ForceField('amber99sbildn.xml')
        system = forcefield.createSystem(pdb.topology, 
            nonbondedMethod=app.CutoffNonPeriodic, nonbondedCutoff=1.0*unit.nanometers, 
            constraints=app.HBonds)
        self.simulation = app.Simulation(pdb.topology, system, mm.VerletIntegrator(0.002*unit.picoseconds))
        self.simulation.context.setPositions(pdb.positions)
        
    def test_1(self):
        file = tempfile.NamedTemporaryFile()
        self.simulation.reporters.append(app.CheckpointReporter(file, 1))
        self.simulation.step(1)
        
        with open(file.name, 'rb') as f:
            self.simulation.context.loadCheckpoint(f.read())

        file.close()

if __name__ == '__main__':
    unittest.main()
