from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk import unit as units
import force_groups
import numpy as np
import sys

forcefield = app.ForceField('./charmm22star.xml')
eu = units.kilocalories_per_mole

protein_system = "hp35"

pdb_filename = "./%s/native.pdb"%protein_system
pdb = PDBFile(pdb_filename)

system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff,constraints=None)

integrator = VerletIntegrator(0.000003*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
forces = force_groups.create_force_groups(simulation)
simulation.context.setPositions(pdb.positions)

total_energy,energies = force_groups.get_energies(simulation)

s = simulation.context.getState(getEnergy=True)

dms_energies = np.loadtxt("./%s/desmond/egrp"%protein_system,dtype='str',skiprows=1)
dms_angle = dms_energies[3,2].astype('float')*eu
dms_cmap = dms_energies[4,2].astype('float')*eu
dms_dihedral = dms_energies[5,2].astype('float')*eu
dms_improper = dms_energies[8,2].astype('float')*eu
dms_stretch = dms_energies[11,2].astype('float')*eu
dms_total = dms_energies[12,2].astype('float')*eu

dms_far_excl = dms_energies[6,2].astype('float')*eu
dms_far_terms = dms_energies[7,2].astype('float')*eu

dms_pair_vdw = dms_energies[10,2].astype('float')*eu
dms_pair_elec = dms_energies[9,2].astype('float')*eu

#I couldn't get agreement with Desmond, but Gromacs looks OK so we're using that
if protein_system == "2JOF":
	gromacs_nb_terms = np.array([ 2.15805e+02  ,  4.39477e+03  , -3.42246e+02 ,   0.00000e+00  , -5.13550e+03])
elif protein_system == "hp35":
	gromacs_nb_terms = np.array([ 5.25292e+02  ,  8.01089e+03 ,  -9.15275e+02 ,   0.00000e+00  , -1.03375e+04])

dms_nb = sum(gromacs_nb_terms)*units.kilojoules_per_mole

print("************************************")
print("OpenMM Energies:")
print("Bond + Urey-Bradley: ",energies[0]+energies[4])
print("Angle : ",energies[1])
print("Dihedral : ",energies[2])
print("Improper : ",energies[3])
print("CMAP : ",energies[6])
print("Nonbonded : ",energies[5])


print("************************************")
print("Energy Comparison:")
print("Bond + Urey-Bradley Difference: ",energies[0]+energies[4] - dms_stretch)
print("Angle Difference: ",energies[1] - dms_angle)
print("Dihedral Difference: ",energies[2] - dms_dihedral)
print("Improper Difference: ",energies[3] - dms_improper)
print("CMAP Difference: ",energies[6] - dms_cmap)
print("Nonbonded Difference: ",energies[5] - dms_nb)

