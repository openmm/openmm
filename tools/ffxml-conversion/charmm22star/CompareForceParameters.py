from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk import unit as units
import Comparison

protein_system = "hp35"

cms_filename = "./%s/desmond/out.dms"%protein_system
pdb_filename = "./%s/native.pdb"%protein_system

pdb = PDBFile(pdb_filename)

forcefield = app.ForceField('./charmm22star.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,nonbondedCutoff=1.1*units.nanometer, constraints=None)

integrator = VerletIntegrator(0.0001*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)

simulation.context.setPositions(pdb.positions)

comp = Comparison.Comparison(cms_filename,pdb,system)
comp.compare()
