from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

gro = GromacsGroFile('input.gro')
top = GromacsTopFile('input.top', periodicBoxVectors=gro.getPeriodicBoxVectors())
system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(top.topology, system, integrator)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
