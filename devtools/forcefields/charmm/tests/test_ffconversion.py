from simtk.openmm import *
from simtk.openmm.app import *
import simtk.unit as u
from parmed.charmm import CharmmPsfFile, CharmmParameterSet
import parmed.openmm as openmm

platform = Platform.getPlatformByName('Reference')

def compare_energies_decomposition(ff, param, pdb, structure):
    system_omm = ff.createSystem(pdb.topology)
    system_parmed = structure.createSystem(param)
    structure.positions = pdb.positions
    structure_omm = openmm.load_topology(pdb.topology, system_omm, xyz=pdb.positions)
    parmed_energies = openmm.energy_decomposition_system(structure, system_parmed)
    omm_energies = openmm.energy_decomposition_system(structure_omm, system_omm)

    print('pamed energies: %s' % parmed_energies)
    print('omm eneriges %s' % omm_energies)

def compare_energies(ff, param, pdb, structure, plat):
    system_omm = ff.createSystem(pdb.topology)
    system_parmed = structure.createSystem(param)
    structure.positions = pdb.positions
    structure_omm = openmm.load_topology(pdb.topology, system_omm, xyz=pdb.positions)

    con_omm = Context(system_omm, VerletIntegrator(2*u.femtoseconds), plat)
    con_omm.setPositions(pdb.positions)
    state_omm = con_omm.getState(getEnergy=True)
    ene_omm = state_omm.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)

    con_charmm = Context(system_parmed, VerletIntegrator(2*u.femtoseconds), plat)
    con_charmm.setPositions(pdb.positions)
    state_charmm = con_charmm.getState(getEnergy=True)
    ene_charmm = state_charmm.getPotentialEnergy().value_in_unit(u.kilocalories_per_mole)

    print('omm total energy')
    print(ene_omm)
    print('charmm total energy')
    print(ene_charmm)




ff = ForceField('../ffxml/charmm36.xml')
pdb = PDBFile('N_N_dimethylaniline.pdb')
structure = CharmmPsfFile('N_N_dimethylaniline.psf')
param = CharmmParameterSet('../toppar/top_all36_cgenff.rtf',
                           '../toppar/par_all36_cgenff.prm',
                           '../toppar/toppar_water_ions.str')

compare_energies_decomposition(ff, param, pdb, structure)
compare_energies(ff, param, pdb, structure, platform)

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
