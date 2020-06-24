from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import mdtraj as mdt
import argparse

prmtop = AmberPrmtopFile('ala2.prmtop')
inpcrd = AmberInpcrdFile('ala2.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
platform = Platform.getPlatformByName('OpenCL')
integrator = LangevinIntegrator(0*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
print('Minimizing')
simulation.minimizeEnergy()
print('Heating from 0 to 300K')
simulation.reporters.append(StateDataReporter(stdout, 100, step=True,
        potentialEnergy=True, temperature=True))
# Initialize the simulation XTC writer from some random state, and then replace the xyz coordinates
simulation.reporters.append(XTCReporter('traj.xtc', 1000, True))
for temp in np.arange(0.0,301.0):
    integrator.setTemperature(temp*kelvin)
    simulation.step(100)
n_frames = 50
for i in range(n_frames):
    simulation.step(1000)
