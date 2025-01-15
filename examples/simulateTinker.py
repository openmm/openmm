from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

tinker = TinkerFiles(xyz='amoeba_solvated_phenol.xyz', key=['amoeba_phenol.prm', 'amoebabio18.prm'])
system = tinker.createSystem(
    polarization="mutual",
    mutualInducedTargetEpsilon=1e-5,
    nonbondedMethod=PME,
    useDispersionCorrection=True,
)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(tinker.topology, system, integrator)
simulation.context.setPositions(tinker.getPositions())
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
