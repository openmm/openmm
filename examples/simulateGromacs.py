from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile('input.gro')
top = GromacsTopFile('input.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
