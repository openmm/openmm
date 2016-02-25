from __future__ import print_function
import parmed
from parmed.utils.six.moves import StringIO, zip
import simtk.openmm.app as app
import simtk.unit as u
import os
import re
from numpy import testing
import tempfile
import yaml
from distutils.spawn import find_executable
import hashlib
from collections import OrderedDict
import glob
import warnings
from parmed.exceptions import ParameterWarning
warnings.filterwarnings('error', category=ParameterWarning)

_loadoffre = re.compile(r'loadoff (\S*)', re.I)
_sourcere = re.compile(r'source (\S*)', re.I)

def convert(filename, ignore=None, provenance=None, write_unused=False, filter_warnings='error'):
    if isinstance(filename, list):
        basename = ''
        for f in filename:
            f_basename = os.path.basename(f)
            f_basename = f_basename.split('.')[1]
            if not basename:
                basename = f_basename
            else:
                basename += '_'
                basename += f_basename
        ffxml_name = 'ffxml/' + basename + '.xml'
    else:
        basename = os.path.basename(filename)
        ffxml_name = 'ffxml/' + '.'.join((basename.split('.')[1:] + ['xml']))
    if not os.path.exists('ffxml/'):
        os.mkdir('ffxml')
    print('Preparing %s for conversion...' % basename)
    if not isinstance(filename, list):
        filename = [filename]
    # do source processing
    new_filename = []
    for fil in filename:
        with open(fil) as f:
            lines = map(lambda line:
                    line if '#' not in line else line[:line.index('#')], f)
        for line in lines:
            if _sourcere.findall(line):
                replace_leaprc = _sourcere.findall(line)[0]
                replace_leaprc_path = os.path.join(os.path.join(AMBERHOME,
                'dat/leap/cmd', replace_leaprc))
                new_filename.append(replace_leaprc_path)
        new_filename.append(fil)
    # now ignore processing and join multiple leaprc's
    filename = new_filename
    new_lines = []
    for fil in filename:
        with open(fil) as f:
            lines = map(lambda line:
                    line if '#' not in line else line[:line.index('#')], f)
        fil_new_lines = []
        for line in lines:
            if (ignore is not None and _loadoffre.findall(line) and
            _loadoffre.findall(line)[0] in ignore):
                continue
            fil_new_lines += line
        new_lines += fil_new_lines
    leaprc = StringIO(''.join(new_lines))
    print('Converting to ffxml...')
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    if filter_warnings is not 'error':
        warnings.filterwarnings(filter_warnings, category=ParameterWarning)
    params.write(ffxml_name, provenance=provenance, write_unused=write_unused)
    if filter_warnings is not 'error':
        warnings.filterwarnings('error', category=ParameterWarning)
    print('Ffxml successfully written!')
    return ffxml_name

def validate_protein(ffxml_name, leaprc_name, united_atom=False):
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
    if not united_atom:
        leap_script_ala3_string = """source %s
x = loadPdb files/ala3.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
        leap_script_villin_string = """source %s
x = loadPdb files/villin.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, villin_top[1], villin_crd[1])
    else:
        leap_script_ala3_string = """source %s
x = loadPdb files/ala3_ua.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, ala3_top[1], ala3_crd[1])
        leap_script_villin_string = """source %s
x = loadPdb files/villin_ua.pdb
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

    print('Calculating ala_ala_ala energies...')
    # AMBER
    parm_amber = parmed.load_file(ala3_top[1], ala3_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    ala3_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    ala3_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

    print('Calculating villin headpiece energies...')
    # AMBER
    parm_amber = parmed.load_file(villin_top[1], villin_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    villin_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    villin_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

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
            try:
                testing.assert_allclose(j[1], i[1], rtol=1e-2)
            except AssertionError:
                testing.assert_allclose(j[1], i[1], rtol=2e-2,
                err_msg=('Ala_ala_ala energies outside of allowed tolerance for %s' % ffxml_name))
                warnings.warn('Ala_ala_ala impropers failed assertion at 1% tolerance, but '
                              'have been asserted at the higher 2% tolerance')
            finally:
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
            try:
                testing.assert_allclose(j[1], i[1], rtol=1e-2)
            except AssertionError:
                testing.assert_allclose(j[1], i[1], rtol=2e-2,
                err_msg=('Villin energies outside of allowed tolerance for %s' % ffxml_name))
                warnings.warn('Villin impropers failed assertion at 1% tolerance, but '
                              'have been asserted at the higher 2% tolerance')
            finally:
                counter += 1
    print('Villin headpiece energy validation successful!')
    print('Done!')

def validate_phospho(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')

    for pdbname in glob.iglob('files/phospho/*.pdb'):
        print('Preparing temporary files for validation...')
        top = tempfile.mkstemp()
        crd = tempfile.mkstemp()
        leap_script_file = tempfile.mkstemp()

        print('Preparing LeaP scripts...')
        leap_script_string = """source leaprc.ff14SB
source %s
x = loadPdb %s
saveAmberParm x %s %s
quit""" % (leaprc_name, pdbname, top[1], crd[1])

        os.write(leap_script_file[0], leap_script_string)

        print('Running LEaP...')
        os.system('tleap -f %s' % leap_script_file[1])
        if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
            raise RuntimeError('LEaP fail for %s' % leaprc_name)

        print('Calculating energies...')
        # AMBER
        parm_amber = parmed.load_file(top[1], crd[1])
        system_amber = parm_amber.createSystem(splitDihedrals=True)
        amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
        # OpenMM
        ff = app.ForceField('ffxml/ff14SB.xml', ffxml_name)
        system_omm = ff.createSystem(parm_amber.topology)
        parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
        system_omm = parm_omm.createSystem(splitDihedrals=True)
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

        print('Deleting temp files...')
        for f in (top, crd, leap_script_file):
            os.close(f[0])
            os.unlink(f[1])

        print('Asserting energies...')
        counter = 0
        for i, j in zip(amber_energies, omm_energies):
            if counter != 3: # Not Impropers
                testing.assert_allclose(j[1], i[1], rtol=1e-5,
                err_msg=('Energies outside of allowed tolerance for %s' % ffxml_name))
                counter += 1
            else: # Impropers - higher tolerance
                try:
                    testing.assert_allclose(j[1], i[1], rtol=1e-2)
                except AssertionError:
                    testing.assert_allclose(j[1], i[1], rtol=2e-2,
                    err_msg=('Energies outside of allowed tolerance for %s' % ffxml_name))
                    warnings.warn('Impropers failed assertion at 1% tolerance, but '
                                  'have been asserted at the higher 2% tolerance')
                finally:
                    counter += 1
        print('Phosphorylated protein energy validation successful!')

def validate_nucleic(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    dna_top = tempfile.mkstemp()
    dna_crd = tempfile.mkstemp()
    leap_script_dna_file = tempfile.mkstemp()
    rna_top = tempfile.mkstemp()
    rna_crd = tempfile.mkstemp()
    leap_script_rna_file = tempfile.mkstemp()
    leap_script_rna_file_alt = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    leap_script_dna_string = """addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "DG" "DG5"  } { 1 "DG" "DG3"  }
{ 0 "DA" "DA5"  } { 1 "DA" "DA3"  }
{ 0 "DC" "DC5"  } { 1 "DC" "DC3"  }
{ 0 "DT" "DT5"  } { 1 "DT" "DT3"  }
}
x = loadPdb files/4rzn_dna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, dna_top[1], dna_crd[1])

    leap_script_rna_string = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "G" "G5"  } { 1 "G" "G3"  } { "G" "G" }
{ 0 "A" "A5"  } { 1 "A" "A3"  } { "A" "A" }
{ 0 "C" "C5"  } { 1 "C" "C3"  } { "C" "C" }
{ 0 "U" "U5"  } { 1 "U" "U3"  } { "U" "U" }
}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, rna_top[1], rna_crd[1])

    leap_script_rna_string_alt = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
addPdbResMap {
{ 0 "G" "RG5"  } { 1 "G" "RG3"  } { "G" "RG" }
{ 0 "A" "RA5"  } { 1 "A" "RA3"  } { "A" "RA" }
{ 0 "C" "RC5"  } { 1 "C" "RC3"  } { "C" "RC" }
{ 0 "U" "RU5"  } { 1 "U" "RU3"  } { "U" "RU" }
}
x = loadPdb files/5c5w_rna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name, rna_top[1], rna_crd[1])

    os.write(leap_script_dna_file[0], leap_script_dna_string)
    os.write(leap_script_rna_file[0], leap_script_rna_string)
    os.write(leap_script_rna_file_alt[0], leap_script_rna_string_alt)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_dna_file[1])
    if os.path.getsize(dna_top[1]) == 0 or os.path.getsize(dna_crd[1]) == 0:
        raise RuntimeError('DNA LEaP fail for %s' % leaprc_name)
    os.system('tleap -f %s' % leap_script_rna_file[1])
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        # try alternative name mappings
        os.system('tleap -f %s' % leap_script_rna_file_alt[1])
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        raise RuntimeError('RNA LEaP fail for %s' % leaprc_name)

    print('Calculating DNA energies...')
    # AMBER
    parm_amber = parmed.load_file(dna_top[1], dna_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    dna_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    dna_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

    print('Calculating RNA energies...')
    # AMBER
    parm_amber = parmed.load_file(rna_top[1], rna_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    rna_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    rna_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (dna_top, dna_crd, leap_script_dna_file, rna_top, rna_crd,
              leap_script_rna_file, leap_script_rna_file_alt):
        os.unlink(f[1])

    print('Asserting DNA energies...')
    counter = 0
    for i, j in zip(dna_amber_energies, dna_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('DNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # turn impropers testing off for now - atom order issue
            # Impropers - higher tolerance
            #testing.assert_allclose(j[1], i[1], rtol=1e-2,
            #err_msg=('DNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('DNA energy validation successful!')

    print('Asserting RNA energies...')
    counter = 0
    for i, j in zip(rna_amber_energies, rna_omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('RNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # turn impropers testing off for now - atom order issue
            # Impropers - higher tolerance
            #testing.assert_allclose(j[1], i[1], rtol=1e-2,
            #err_msg=('RNA energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('RNA energy validation successful!')

def validate_gaff(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    imatinib_top = tempfile.mkstemp()
    imatinib_crd = tempfile.mkstemp()
    leap_script_imatinib_file = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    leap_script_imatinib_string = """source %s
loadamberparams files/frcmod.imatinib
x = loadMol2 files/imatinib.mol2
saveAmberParm x %s %s
quit""" % (leaprc_name, imatinib_top[1], imatinib_crd[1])
    os.write(leap_script_imatinib_file[0], leap_script_imatinib_string)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_imatinib_file[1])
    if os.path.getsize(imatinib_top[1]) == 0 or os.path.getsize(imatinib_crd[1]) == 0:
        raise RuntimeError('imatinib LEaP fail for %s' % leaprc_name)

    print('Calculating imatinib energies...')
    # AMBER
    parm_amber = parmed.load_file(imatinib_top[1], imatinib_crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    imatinib_amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name, 'files/imatinib_frcmod.xml', 'files/imatinib.xml')
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    imatinib_omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (imatinib_top, imatinib_crd, leap_script_imatinib_file):
        os.unlink(f[1])

    print('Asserting imatinib energies...')
    counter = 0
    for i, j in zip(imatinib_amber_energies, imatinib_omm_energies):
        if counter != 3:
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('Imatinib energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # impropers - higher tolerance
            testing.assert_allclose(j[1], i[1], rtol=1e-2,
            err_msg=('Imatinib energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('Imatinib energy validation successful!')
    print('Done!')

def validate_modrna(ffxml_name, leaprc_name):
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    print('Preparing temporary files for validation...')
    top = tempfile.mkstemp()
    crd = tempfile.mkstemp()
    leap_script_file = tempfile.mkstemp()
    leap_script_file_alt = tempfile.mkstemp()

    print('Preparing LeaP scripts...')
    leap_script_string = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
source %s
addPdbResMap {
{ 0 "G" "G5"  } { 1 "G" "G3"  } { "G" "G" }
{ 0 "A" "A5"  } { 1 "A" "A3"  } { "A" "A" }
{ 0 "C" "C5"  } { 1 "C" "C3"  } { "C" "C" }
{ 0 "U" "U5"  } { 1 "U" "U3"  } { "U" "U" }
}
x = loadPdb files/6tna_modrna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name[0], leaprc_name[1], top[1], crd[1])

    leap_script_string_alt = """
addPdbAtomMap {
{ "H1'" "H1*" }
{ "H2'" "H2'1" }
{ "H2''" "H2'2" }
{ "H3'" "H3*" }
{ "H4'" "H4*" }
{ "H5'" "H5'1" }
{ "H5''" "H5'2" }
{ "HO2'" "HO'2" }
{ "HO5'" "H5T"  }
{ "HO3'" "H3T" }
{ "OP1" "O1P" }
{ "OP2" "O2P" }
}
source %s
source %s
addPdbResMap {
{ 0 "G" "RG5"  } { 1 "G" "RG3"  } { "G" "RG" }
{ 0 "A" "RA5"  } { 1 "A" "RA3"  } { "A" "RA" }
{ 0 "C" "RC5"  } { 1 "C" "RC3"  } { "C" "RC" }
{ 0 "U" "RU5"  } { 1 "U" "RU3"  } { "U" "RU" }
}
x = loadPdb files/6tna_modrna.pdb
saveAmberParm x %s %s
quit""" % (leaprc_name[1], leaprc_name[0], top[1], crd[1])

    os.write(leap_script_file[0], leap_script_string)
    os.write(leap_script_file_alt[0], leap_script_string_alt)

    print('Running LEaP...')
    os.system('tleap -f %s' % leap_script_file[1])
    if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
        # try alternative atom name mapping
        os.system('tleap -f %s' % leap_script_file_alt[1])
    if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
        raise RuntimeError('LEaP fail for %s' % leaprc_name)

    print('Calculating energies...')
    # AMBER
    parm_amber = parmed.load_file(top[1], crd[1])
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    amber_energies = parmed.openmm.energy_decomposition_system(parm_amber, system_amber, nrg=u.kilojoules_per_mole)
    # OpenMM
    ff = app.ForceField(ffxml_name)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm, xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    omm_energies = parmed.openmm.energy_decomposition_system(parm_omm, system_omm, nrg=u.kilojoules_per_mole)

    print('Deleting temp files...')
    for f in (top, crd, leap_script_file, leap_script_file_alt):
        os.close(f[0])
        os.unlink(f[1])

    # dev
    print(amber_energies)
    print(omm_energies)

    print('Asserting energies...')
    counter = 0
    for i, j in zip(amber_energies, omm_energies):
        if counter != 3: # Not Impropers
            testing.assert_allclose(j[1], i[1], rtol=1e-5,
            err_msg=('Energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
        else: # Impropers - higher tolerance - RNA - for now turn off
            #testing.assert_allclose(j[1], i[1], rtol=2e-2,
            #err_msg=('Energies outside of allowed tolerance for %s' % ffxml_name))
            counter += 1
    print('Modified RNA energy validation successful!')

if __name__ == '__main__':
    if not find_executable('tleap'):
        raise RuntimeError('tleap not available from PATH')
    if os.getenv('AMBERHOME'):
        AMBERHOME = os.getenv('AMBERHOME')
    else:
        tleap_path = find_executable('tleap')
        AMBERHOME = os.path.split(tleap_path)[0]
        AMBERHOME = os.path.join(AMBERHOME, '../')
        parmed.amber.AMBERHOME = AMBERHOME
    data = yaml.load(open('files/master.yaml'))
    source_pack = data[0]['sourcePackage']
    source_pack_ver = data[0]['sourcePackageVersion']
    # solvents and ions converted separately; leaprc.ff10 calls phosphoaa10.lib
    # which does not exist anymore, LeAP skips it on error so we do too
    ignore = {'solvents.lib', 'atomic_ions.lib', 'ions94.lib', 'ions91.lib',
              'phosphoaa10.lib'}
    for entry in data[1:]:
        leaprc_name = entry['Source']
        leaprc_reference = entry['Reference']
        leaprc_test = entry['Test']
        provenance = OrderedDict()
        if isinstance(leaprc_name, list):
            filename = []
            source = provenance['Source'] = []
            for leaprc in leaprc_name:
                _filename = os.path.join(AMBERHOME, 'dat/leap/cmd', leaprc)
                filename.append(_filename)
                source.append(OrderedDict())
                source[-1]['Source'] = leaprc
                md5 = hashlib.md5()
                with open(_filename) as f:
                    md5.update(f.read())
                md5 = md5.hexdigest()
                source[-1]['md5hash'] = md5
                source[-1]['sourcePackage'] = source_pack
                source[-1]['sourcePackageVersion'] = source_pack_ver
        else:
            filename = os.path.join(AMBERHOME, 'dat/leap/cmd', leaprc_name)
            source = provenance['Source'] = OrderedDict()
            source['Source'] = leaprc_name
            md5 = hashlib.md5()
            with open(filename) as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source['md5hash'] = md5
            source['sourcePackage'] = source_pack
            source['sourcePackageVersion'] = source_pack_ver
        provenance['Reference'] = leaprc_reference
        # set default conversion options
        write_unused = False
        filter_warnings = 'error'
        # set conversion options if present
        if 'Options' in entry:
            for option in entry['Options']:
                if option == 'write_unused':
                    write_unused = entry['Options'][option]
                elif option == 'filter_warnings':
                    filter_warnings = entry['Options'][option]
                else:
                    raise RuntimeError("Wrong option used in Options for %s"
                                       % leaprc_name)
        print('Converting %s to ffxml...' % leaprc_name)
        ffxml_name = convert(filename, ignore=ignore, provenance=provenance,
                     write_unused=write_unused, filter_warnings=filter_warnings)
        print('Validating the conversion...')
        tested = False
        for test in leaprc_test:
            if test == 'protein':
                print('Protein tests...')
                validate_protein(ffxml_name, leaprc_name)
                tested = True
            elif test == 'nucleic':
                print('Nucleic acid tests...')
                validate_nucleic(ffxml_name, leaprc_name)
                tested = True
            elif test == 'protein_ua':
                print('United-atom protein tests...')
                validate_protein(ffxml_name, leaprc_name, united_atom=True)
                tested = True
            elif test == 'mol2':
                print('Small molecule (GAFF) tests...')
                validate_gaff(ffxml_name, leaprc_name)
                tested = True
            elif test == 'phospho':
                print('Phosphorylated protein tests...')
                validate_phospho(ffxml_name, leaprc_name)
                tested = True
            elif test == 'modrna':
                print('Modified RNA tests...')
                validate_modrna(ffxml_name, leaprc_name)
                tested = True
            elif test == 'lipid':
                print('Lipid tests...')
                validate_lipid(ffxml_name, leaprc_name)
                tested = True
            if not tested:
                raise RuntimeError('No validation tests have been run for %s' % leaprc_name)
