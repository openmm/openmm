from __future__ import print_function
import parmed
from parmed.utils.six.moves import StringIO, zip
import simtk.openmm.app as app
from simtk.unit import *
import os
import re
from numpy import testing
import tempfile
import yaml
from distutils.spawn import find_executable
import warnings
warnings.filterwarnings('error')

_loadoffre = re.compile(r'loadoff (\S*)', re.I)

def convert(filename, ignore=None, reference=''):
    basename = os.path.basename(filename)
    if not os.path.exists('ffxml/'):
        os.mkdir('ffxml')
    ffxml_name = 'ffxml/' + '.'.join((basename.split('.')[1:] + ['xml']))
    provenance = {}
    provenance['Reference'] = reference
    print('Preparing %s for conversion...' % basename)
    with open(filename) as f:
        lines = map(lambda line:
                    line if '#' not in line else line[:line.index('#')], f)
    if ignore is not None:
        new_lines = []
        for line in lines:
            if _loadoffre.findall(line) and _loadoffre.findall(line)[0] in ignore:
                continue
            else:
                new_lines.append(line)
    else:
        new_lines = lines
    leaprc = StringIO(''.join(new_lines))
    print('Converting to ffxml...')
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    params.write(ffxml_name, provenance=provenance, write_unused=False)
    print('Ffxml successfully written!')
    return ffxml_name

def validate_protein(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    ala3_top = tempfile.mkstemp()
    ala3_crd = tempfile.mkstemp()
    villin_top = tempfile.mkstemp()
    villin_crd = tempfile.mkstemp()
    leap_script_ala3_file = tempfile.mkstemp()
    leap_script_villin_file = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    leap_script_ala3_string = """source %s
x = loadPdb files/ala3.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
    leap_script_villin_string = """source %s
x = loadPdb files/villin.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, villin_top[1], villin_crd[1])
    os.write(leap_script_ala3_file[0], leap_script_ala3_string)
    os.write(leap_script_villin_file[0], leap_script_villin_string)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_ala3_file[1])
    if os.path.getsize(ala3_top[1]) == 0 or os.path.getsize(ala3_crd[1]) == 0:
        raise RuntimeError('Ala_ala_ala LEaP fail for %s' % leaprc_name)
    os.system('tleap -f %s' % leap_script_villin_file[1])
    if os.path.getsize(villin_top[1]) == 0 or os.path.getsize(villin_crd[1]) == 0:
        raise RuntimeError('Villin headpiece LEaP fail for %s' % leaprc_name)

    # Test ala3
    # AMBER
    print('Calculating ala_ala_ala energies...')
    parm_amber = parmed.load_file(ala3_top[1], ala3_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    ala3_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    ala3_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)
    # Test villin headpiece
    # AMBER
    print('Calculating villin headpiece energies...')
    parm_amber = parmed.load_file(villin_top[1], villin_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    villin_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    villin_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (ala3_top, ala3_crd, villin_top, villin_crd, leap_script_ala3_file,
             leap_script_villin_file):
        os.unlink(f[1])

    print('Asserting ala_ala_ala energies...')
    counter = 0
    for i, j in zip(ala3_amber_energies, ala3_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('Ala_ala_ala energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance
            testing.assert_allclose(j[1], i[1], rtol=1e-2,
            err_msg=('Ala_ala_ala energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('Ala_ala_ala energy validation successful!')

    print('Asserting villin headpiece energies...')
    counter = 0
    for i, j in zip(villin_amber_energies, villin_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('Villin energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance
            testing.assert_allclose(j[1], i[1], rtol=1e-2,
            err_msg=('Villin energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('Villin headpiece energy validation successful!')
    print('Done!')

if __name__ == '__main__':
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    if os.getenv('AMBERHOME'):
        AMBERHOME = os.getenv('AMBERHOME')
    else:
        tleap_path = find_executable('tleap')
        AMBERHOME = os.path.split(tleap_path)[0]
        AMBERHOME = os.path.join(AMBERHOME, '../')
    data = yaml.load(open('files/master.yaml'))
    ignore = {'solvents.lib', 'atomic_ions.lib', 'ions94.lib', 'ions91.lib'}
    for entry in data:
        leaprc_name = entry['Source']
        leaprc_reference = entry['Reference']
        print('Converting %s to ffxml...' % leaprc_name)
        filename = os.path.join(AMBERHOME, 'dat/leap/cmd', leaprc_name)
        ffxml_name = convert(filename, ignore=ignore, reference=leaprc_reference)
        print('Validating the conversion...')
        validate_protein(ffxml_name, leaprc_name)
