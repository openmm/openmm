from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

tinker = TinkerFiles('amoeba_solvated_phenol.xyz', ['amoeba_phenol.prm', 'amoebabio18.prm'])
system = tinker.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.7*nanometer, vdwCutoff=0.9*nanometer)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(tinker.topology, system, integrator)
simulation.context.setPositions(tinker.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
