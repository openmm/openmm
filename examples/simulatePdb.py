# Import OpenMM modules.
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# Load topology and positions from the PDB file.
pdb = PDBFile('input.pdb')
# Load the forcefield parameters for AMBER99SB and TIP3P water.
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
# Create a System object from the topology defined in the PDB file.
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer, constraints=HBonds)
# Create a Langevin integrator with specified temperature, collision rate, and timestep.
temperature = 300*kelvin
collision_rate = 1/picosecond
timestep = 0.002*picoseconds
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
# Create a Simulation from the topology, system, and integrator.
simulation = Simulation(pdb.topology, system, integrator)
# Set the initial atomic positions for the simulation from the PDB file.
simulation.context.setPositions(pdb.positions)
# Minimize the energy prior to simulation.
simulation.minimizeEnergy()
# Add a few reporters to generate output during the simulation.
report_interval = 1000
simulation.reporters.append(PDBReporter('output.pdb', report_interval))
simulation.reporters.append(StateDataReporter(stdout, report_interval, step=True,
                                              potentialEnergy=True, temperature=True))
# Run the simulation for a specified number of timesteps.
simulation.step(10000)
