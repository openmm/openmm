# Import OpenMM modules.
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# Load the Amber format parameters and topology files.
prmtop = AmberPrmtopFile('input.prmtop')
inpcrd = AmberInpcrdFile('input.inpcrd', loadBoxVectors=True)
# Create a System object using the parameters defined in the prmtop file.
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
# Create a Langevin integrator with specified temperature, collision rate, and timestep.
temperature = 300*kelvin
collision_rate = 1/picosecond
timestep = 0.002*picoseconds
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
# Create a Simulation from the topology, system, and integrator.
simulation = Simulation(prmtop.topology, system, integrator)
# Set the initial atomic positions and box vectors for the simulation from the inpcrd file.
simulation.context.setPositions(inpcrd.positions)
simulation.context.setBoxVectors(inpcrd.getBoxVectors)
# Minimize the energy prior to simulation.
simulation.minimizeEnergy()
# Add a few reporters to generate output during the simulation.
report_interval = 1000
simulation.reporters.append(PDBReporter('output.pdb', report_interval))
simulation.reporters.append(StateDataReporter(stdout, report_interval, step=True,
                                              potentialEnergy=True, temperature=True))
# Run the simulation for a specified number of timesteps.
simulation.step(10000)
