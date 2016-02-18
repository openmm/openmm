from simtk.openmm import *
from simtk.openmm.app import *
import simtk.unit as u
from parmed.charmm import CharmmPsfFile, CharmmParameterSet
import parmed.openmm as openmm

platform = Platform.getPlatformByName('Reference')

ff = ForceField('../ffxml/charmm36.xml')
pdb = PDBFile('methanol_solvated.pdb')
# # #pdb2 = app.PDBFile('../../../../../devtools/forcefields/charmm/tests/ala_ala_ala.pdb')
system_omm = ff.createSystem(pdb.topology)
integrator1 = LangevinIntegrator(300*unit.kelvin, 1.0/u.picoseconds, 1.0*u.femtosecond)
sim = Simulation(pdb.topology, system_omm, integrator1, platform)
sim.context.setPositions(pdb.positions)
state = sim.context.getState(getEnergy=True)
parm_omm = openmm.load_topology(pdb.topology, system_omm, xyz=pdb.positions)

param = CharmmParameterSet('../charmm/toppar/top_all36_cgenff.rtf',
                            '../charmm/toppar/par_all36_cgenff.prm',
                            '../charmm/toppar/toppar_water_ions.str')
structure = CharmmPsfFile('methanol_solvated.psf')
structure.load_parameters(param)
structure.positions = pdb.positions
system_parmed = structure.createSystem(param)

integrator2 = LangevinIntegrator(300*unit.kelvin, 1.0/u.picoseconds, 1.0*u.femtosecond)
sim2 = Simulation(pdb.topology, system_parmed, integrator2, platform)
sim2.context.setPositions(pdb.positions)
state2 = sim2.context.getState(getEnergy=True)

print('potential energy from openmm %s' % state.getPotentialEnergy())
print ('nonbonded exceptions from openmm and charmm36.xml %s' % system_omm.getForces()[-2].getNumExceptions())
print('energy decompostion from openmm with charmm36.xml %s' % openmm.energy_decomposition_system(parm_omm, system_omm))
print('potential energy from parmed system %s' % state2.getPotentialEnergy())
print ('nonbonded exceptions from parmed openmm system %s' % system_parmed.getForces()[-2].getNumExceptions())
print('energy decomposition from parmed openmm system %s' % openmm.energy_decomposition_system(structure, system_parmed))

# def test_ffconversion():
#     psf = CharmmPsfFile('../ala_ala_ala.psf')
#     params = CharmmParameterSet('../charmm/toppar/top_all36_prot.rtf', '../charmm/toppar/par_all36_prot.prm')
#     pdb = PDBFile('../ala_ala_ala.pdb')
#
#     system1 = psf.createSystem(params=params)
#     integrator = VerletIntegrator(1*femtosecond)
#     platform = Platform.getPlatformByName('Reference')
#
#     con = Context(system1, integrator, platform )
#     con.setPositions(pdb.positions)
#
#     energy1 = con.getState(getEnergy=True).getPotentialEnergy()
#
#     ff = ForceField('../ffxml/charmm36_prot.xml')
#     system2 = ff.createSystem(pdb.topology)
#
#     con2 = Context(system2, VerletIntegrator(1*femtosecond), platform)
#     con2.setPositions(pdb.positions)
#
#     energy2 = con2.getState(getEnergy=True).getPotentialEnergy()
