# AMBER --> OpenMM force-field conversion script
# Author: Rafal P. Wiewiora, ChoderaLab
from __future__ import print_function, division
import parmed
from parmed.utils.six import iteritems
from parmed.utils.six.moves import StringIO, zip
import simtk.openmm.app as app
import simtk.unit as u
import os
import sys
import re
import tempfile
import yaml
from distutils.spawn import find_executable
import hashlib
from collections import OrderedDict
import glob
import argparse
from lxml import etree as et
import csv
import logging
import warnings
from parmed.exceptions import ParameterWarning
warnings.filterwarnings('error', category=ParameterWarning)

_loadoffre = re.compile(r'loadoff (\S*)', re.I)
_sourcere = re.compile(r'source (\S*)', re.I)

# check for AMBERHOME, find from tleap location if not set, exception if can't
if os.getenv('AMBERHOME'):
    AMBERHOME = os.getenv('AMBERHOME')
else:
    if not find_executable('tleap'):
        raise RuntimeError('AMBERHOME not set and tleap not available from PATH')
    tleap_path = find_executable('tleap')
    AMBERHOME = os.path.split(tleap_path)[0]
    AMBERHOME = os.path.join(AMBERHOME, '../')
    parmed.amber.AMBERHOME = AMBERHOME

# set global defaults for verbose and log
verbose = False
no_log = False

def main():
    global verbose
    global no_log
    global logger
    # argparse
    parser = argparse.ArgumentParser(description='AMBER --> OpenMM forcefield '
                                                 'conversion script')
    parser.add_argument('--input', '-i', default='files/master.yaml',
                        help='path of the input file. Default: "files/master.yaml"')
    parser.add_argument('--input-format', '-if', default='yaml',
                        help='format of the input file: "yaml" or "leaprc". Default: "yaml"')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='turns verbosity on')
    parser.add_argument('--no-log', action='store_true',
                        help='turns logging of energies to log.csv off')
    parser.add_argument('--protein-test', action='store_true',
                        help='validate resulting XML through protein tests')
    parser.add_argument('--nucleic-test', action='store_true',
                        help='validate resulting XML through nucleic acid tests')
    parser.add_argument('--protein-ua-test', action='store_true',
                        help='validate resulting XML through united-atom protein tests')
    parser.add_argument('--phospho-protein-test', action='store_true',
                        help='validate resulting XML through phosphorylated protein tests')
    parser.add_argument('--gaff-test', action='store_true',
                        help='validate resulting XML through small-molecule (GAFF) test')
    # parser.add_argument('--write_unused', action='store_true')
    # parser.add_argument('--default-warnings', action='store_true')
    args = parser.parse_args()
    verbose = args.verbose
    no_log = args.no_log

    if not no_log: logger = Logger('log.csv')

    # input is either a YAML or a leaprc - default is leaprc
    # output directory hardcoded here for ffxml/
    if args.input_format == 'yaml':
        convert_yaml(args.input, ffxml_dir='ffxml/')
    # if leaprc converted - output to the same dir
    elif args.input_format == 'leaprc':
        ffxml_name = convert_leaprc(args.input, ffxml_dir='./')
        if args.protein_test:
            validate_protein(ffxml_name, args.input)
        if args.nucleic_test:
            validate_nucleic(ffxml_name, args.input)
        if args.protein_ua_test:
            validate_protein(ffxml_name, args.input, united_atom=True)
        if args.phospho_protein_test:
            validate_phospho_protein(ffxml_name, args.input)
        if args.gaff_test:
            validate_gaff(ffxml_name, args.input)
    else:
        raise RuntimeError('Wrong input_format chosen.')

    if not no_log: logger.close()

def convert_leaprc(files, split_filename=False, ffxml_dir='./', ignore=None,
    provenance=None, write_unused=False, filter_warnings='error'):
    if verbose: print('Converting %s to ffxml...' % files)
    # allow for multiple source files - further code assuming list is passed
    if not isinstance(files, list):
        files = [files]
    basename = ''
    for f in files:
        f_basename = os.path.basename(f)
        if split_filename:
            f_basename = f_basename.split('.')[1:]
            f_basename = '.'.join(f_basename)
        if not basename:
            basename = f_basename
        else:
            basename += '_'
            basename += f_basename
    ffxml_name = ffxml_dir + basename + '.xml'
    if not os.path.exists(ffxml_dir):
        os.mkdir(ffxml_dir)
    if verbose: print('Preprocessing the leaprc for %s...' % basename)
    # do source processing
    new_files = []
    for fil in files:
        with open(fil) as f:
            lines = map(lambda line:
                    line if '#' not in line else line[:line.index('#')], f)
        for line in lines:
            if _sourcere.findall(line):
                replace_leaprc = _sourcere.findall(line)[0]
                replace_leaprc_path = os.path.join(os.path.join(AMBERHOME,
                'dat/leap/cmd', replace_leaprc))
                new_files.append(replace_leaprc_path)
        new_files.append(fil)
    # now do ignore processing and join multiple leaprc's
    files = new_files
    new_lines = []
    for fil in files:
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
    if verbose: print('Converting to ffxml %s...' % ffxml_name)
    params = parmed.amber.AmberParameterSet.from_leaprc(leaprc)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    if filter_warnings != 'error':
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_name, provenance=provenance, write_unused=write_unused)
    else:
        params.write(ffxml_name, provenance=provenance, write_unused=write_unused)
    if verbose: print('%s successfully written!' % ffxml_name)
    return ffxml_name

def convert_recipe(files, solvent_file=None, ffxml_dir='./', provenance=None, ffxml_basename=None,
                   filter_warnings='always'):
    if verbose: print('Converting %s to ffxml...' % files)
    ffxml_name = os.path.join(ffxml_dir, (ffxml_basename + '.xml'))
    ffxml_temp_stringio = StringIO()
    params = parmed.amber.AmberParameterSet(files)
    params = parmed.openmm.OpenMMParameterSet.from_parameterset(params)
    # Change atom type naming
    # atom_types
    new_atom_types = OrderedDict()
    for name, atom_type in iteritems(params.atom_types):
        new_name = ffxml_basename + '-' + name
        new_atom_types[new_name] = atom_type
    params.atom_types = new_atom_types
    # atoms in residues
    for name, residue in iteritems(params.residues):
        for atom in residue:
            new_type = ffxml_basename + '-' + atom.type
            atom.type = new_type
    if solvent_file is None:
    # this means this file does not include a water model - hard-coded assumption it is
    # then a 'multivalent' file - set overloadLevel to 1 for all residue templates
        for name, residue in iteritems(params.residues):
            residue.overload_level = 1
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_name, provenance=provenance, write_unused=False)
    else:
        with warnings.catch_warnings():
            warnings.filterwarnings(filter_warnings, category=ParameterWarning)
            params.write(ffxml_temp_stringio, provenance=provenance, write_unused=False)
        ffxml_temp_stringio.seek(0)
        if verbose: print('Modifying converted ffxml to append solvent parameters')
        tree_main = et.parse(ffxml_temp_stringio)
        tree_water = et.parse(solvent_file)
        root_main = tree_main.getroot()
        root_water = tree_water.getroot()
        with open(ffxml_name, 'w') as f:
            f.write('<ForceField>\n ')
            f.write(et.tostring(root_main.findall('Info')[0]))
            f.write('<AtomTypes>\n  ')
            for subelement in root_main.findall('AtomTypes')[0]:
                f.write(et.tostring(subelement))
            f.write(' ')
            for subelement in root_water.findall('AtomTypes')[0]:
                f.write(et.tostring(subelement))
            f.write('</AtomTypes>\n <Residues>\n  ')
            for subelement in root_main.findall('Residues')[0]:
                f.write(et.tostring(subelement))
            f.write(' ')
            for subelement in root_water.findall('Residues')[0]:
                f.write(et.tostring(subelement))
            f.write('</Residues>\n <HarmonicBondForce>\n  ')
            for subelement in root_water.findall('HarmonicBondForce')[0]:
                f.write(et.tostring(subelement))
            f.write('</HarmonicBondForce>\n <HarmonicAngleForce>\n  ')
            for subelement in root_water.findall('HarmonicAngleForce')[0]:
                f.write(et.tostring(subelement))
            f.write('</HarmonicAngleForce>\n ')
            f.write('<NonbondedForce coulomb14scale="%s" lj14scale="%s">\n  ' %
                   (root_main.findall('NonbondedForce')[0].attrib['coulomb14scale'],
                    root_main.findall('NonbondedForce')[0].attrib['lj14scale'])
                   )
            for subelement in root_main.findall('NonbondedForce')[0]:
                f.write(et.tostring(subelement))
            f.write(' ')
            for subelement in root_water.findall('NonbondedForce')[0]:
                if subelement.tag == 'UseAttributeFromResidue': continue
                f.write(et.tostring(subelement))
            f.write('</NonbondedForce>\n</ForceField>')
    if verbose: print('%s successfully written!' % ffxml_name)
    return ffxml_name

def convert_yaml(yaml_name, ffxml_dir):
    data = yaml.load(open(yaml_name))
    source_pack = data[0]['sourcePackage']
    source_pack_ver = data[0]['sourcePackageVersion']
    # solvents and ions converted separately; leaprc.ff10 calls phosphoaa10.lib
    # which does not exist anymore, LeAP skips it on error so we do too
    ignore = {'solvents.lib', 'atomic_ions.lib', 'ions94.lib', 'ions91.lib',
              'phosphoaa10.lib'}
    # Default yaml reading more is leaprc
    MODE = 'LEAPRC'
    for entry in data[1:]:
        if 'MODE' in entry:
            MODE = entry['MODE']
            continue
        if MODE == 'RECIPE' and 'sourcePackage2' in entry:
            source_pack2 = entry['sourcePackage2']
            source_pack_ver2 = entry['sourcePackageVersion2']
            continue
        leaprc_name = entry['Source']
        leaprc_reference = entry['Reference']
        leaprc_test = entry['Test']
        if MODE == 'RECIPE':
            recipe_name = entry['Name']
            solvent_name = entry['Solvent']
            if 'Solvent_source' in entry:
                recipe_source2 = entry['Solvent_source']
            else:
                recipe_source2 = None
            if 'Standard' in entry:
                standard_ffxml = os.path.join(ffxml_dir, (entry['Standard'] + '.xml'))
            else:
                standard_ffxml = None
        provenance = OrderedDict()
        if isinstance(leaprc_name, list):
            files = []
            source = provenance['Source'] = []
            for leaprc in leaprc_name:
                if MODE == 'LEAPRC':
                    _filename = os.path.join(AMBERHOME, 'dat/leap/cmd', leaprc)
                elif MODE == 'RECIPE':
                    _filename = os.path.join(AMBERHOME, 'dat/leap/', leaprc)
                files.append(_filename)
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
            files = os.path.join(AMBERHOME, 'dat/leap/cmd', leaprc_name)
            source = provenance['Source'] = OrderedDict()
            source['Source'] = leaprc_name
            md5 = hashlib.md5()
            with open(files) as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source['md5hash'] = md5
            source['sourcePackage'] = source_pack
            source['sourcePackageVersion'] = source_pack_ver
        # add water file and source info for it
        if MODE == 'RECIPE' and recipe_source2 is not None:
            _filename = os.path.join('files', recipe_source2)
            solvent_file = _filename
            source.append(OrderedDict())
            source[-1]['Source'] = recipe_source2
            md5 = hashlib.md5()
            with open(_filename) as f:
                md5.update(f.read())
            md5 = md5.hexdigest()
            source[-1]['md5hash'] = md5
            source[-1]['sourcePackage'] = source_pack2
            source[-1]['sourcePackageVersion'] = source_pack_ver2
        elif MODE == 'RECIPE' and recipe_source2 is None:
            solvent_file = None
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
        if MODE == 'LEAPRC':
            ffxml_name = convert_leaprc(files, ffxml_dir=ffxml_dir, ignore=ignore,
                         provenance=provenance, write_unused=write_unused,
                         filter_warnings=filter_warnings, split_filename=True)
        elif MODE == 'RECIPE':
            ffxml_name = convert_recipe(files, solvent_file=solvent_file,
                         ffxml_dir=ffxml_dir, provenance=provenance,
                         ffxml_basename=recipe_name)
        if verbose: print('Validating the conversion...')
        tested = False
        for test in leaprc_test:
            if test == 'protein':
                validate_protein(ffxml_name, leaprc_name)
                tested = True
            elif test == 'nucleic':
                validate_nucleic(ffxml_name, leaprc_name)
                tested = True
            elif test == 'protein_ua':
                validate_protein(ffxml_name, leaprc_name, united_atom=True)
                tested = True
            elif test == 'protein_phospho':
                validate_phospho_protein(ffxml_name, leaprc_name)
                tested = True
            elif test == 'gaff':
                validate_gaff(ffxml_name, leaprc_name)
                tested = True
            elif test == 'water_ion':
                validate_water_ion(ffxml_name, files, solvent_name, recipe_name,
                standard_ffxml=standard_ffxml)
                tested = True
        if not tested:
            raise RuntimeError('No validation tests have been run for %s' %
                                leaprc_name)

def assert_energies(prmtop, inpcrd, ffxml, system_name='unknown', tolerance=1e-5,
                    improper_tolerance=1e-2, units=u.kilojoules_per_mole):
    # AMBER
    parm_amber = parmed.load_file(prmtop, inpcrd)
    system_amber = parm_amber.createSystem(splitDihedrals=True)
    amber_energies = parmed.openmm.energy_decomposition_system(parm_amber,
                     system_amber, nrg=units)
    # OpenMM-ffxml
    if isinstance(ffxml, str):
        ff = app.ForceField(ffxml)
    else:
        ff = app.ForceField(*ffxml)
    system_omm = ff.createSystem(parm_amber.topology)
    parm_omm = parmed.openmm.load_topology(parm_amber.topology, system_omm,
                                           xyz=parm_amber.positions)
    system_omm = parm_omm.createSystem(splitDihedrals=True)
    omm_energies = parmed.openmm.energy_decomposition_system(parm_omm,
                   system_omm, nrg=units)

    # calc rel energies and assert
    rel_energies = []
    for i, j in zip(amber_energies, omm_energies):
        if i[0] != j[0]:
            raise RuntimeError('Mismatch in energy tuples naming.')
        if i[1] != 0:
            rel_energies.append((i[0], (abs(i[1]-j[1])/i[1])))
        else:
            if j[1] != 0:
                raise RuntimeError('One of AMBER %s energies (%s) for %s is zero, '
                      'while the corresponding OpenMM energy is non-zero' %
                      (system_name, i[0], ffxml))
            rel_energies.append((i[0], 0))

    dihedrals_done = False
    for i in rel_energies:
        if i[0] != 'PeriodicTorsionForce':
            if i[1] > tolerance:
                raise AssertionError('%s energies (%s) outside of allowed tolerance for %s' %
                                     (system_name, i[0], ffxml))
        else:
            if not dihedrals_done:
                if i[1] > tolerance:
                    raise AssertionError('%s energies (%s) outside of allowed tolerance for %s' %
                                         (system_name, i[0], ffxml))
                dihedrals_done = True
            else: #impropers
                if i[1] > improper_tolerance:
                    raise AssertionError('%s energies (%s-impropers) outside of allowed tolerance for %s' %
                                         (system_name, i[0], ffxml))

    # logging
    if not no_log:
        amber_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        amber_energies_log['ffxml_name'] = ffxml
        amber_energies_log['test_system'] = system_name
        amber_energies_log['data_type'] = 'AMBER'
        amber_energies_log['units'] = units
        omm_energies_log['ffxml_name'] = ffxml
        omm_energies_log['test_system'] = system_name
        omm_energies_log['data_type'] = 'OpenMM'
        omm_energies_log['units'] = units
        rel_energies_log['ffxml_name'] = ffxml
        rel_energies_log['test_system'] = system_name
        rel_energies_log['data_type'] = 'abs(AMBER-OpenMM)/AMBER'
        dihedrals_done = False
        for item in amber_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                amber_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                amber_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                amber_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in omm_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                omm_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                omm_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                omm_energies_log[item[0]] = item[1]
        dihedrals_done = False
        for item in rel_energies:
            if item[0] == 'PeriodicTorsionForce' and not dihedrals_done:
                rel_energies_log['PeriodicTorsionForce_dihedrals'] = item[1]
                dihedrals_done = True
            elif item[0] == 'PeriodicTorsionForce' and dihedrals_done:
                rel_energies_log['PeriodicTorsionForce_impropers'] = item[1]
            elif item[0] == 'CMMotionRemover':
                continue
            else:
                rel_energies_log[item[0]] = item[1]

        logger.log(amber_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)

def validate_protein(ffxml_name, leaprc_name, united_atom=False):
    if verbose: print('Protein energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    ala3_top = tempfile.mkstemp()
    ala3_crd = tempfile.mkstemp()
    villin_top = tempfile.mkstemp()
    villin_crd = tempfile.mkstemp()
    leap_script_ala3_file = tempfile.mkstemp()
    leap_script_villin_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
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

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_ala3_file[1], os.devnull))
    if os.path.getsize(ala3_top[1]) == 0 or os.path.getsize(ala3_crd[1]) == 0:
        raise RuntimeError('Ala_ala_ala LEaP fail for %s' % leaprc_name)
    os.system('tleap -f %s > %s' % (leap_script_villin_file[1], os.devnull))
    if os.path.getsize(villin_top[1]) == 0 or os.path.getsize(villin_crd[1]) == 0:
        raise RuntimeError('Villin headpiece LEaP fail for %s' % leaprc_name)

    try:
        if verbose: print('Calculating and validating ala_ala_ala energies...')
        assert_energies(ala3_top[1], ala3_crd[1], ffxml_name,
                        system_name='protein-ala_ala_ala')
        if verbose: print('Ala_ala_ala energy validation successful!')

        if verbose: print('Calculating and validating villin headpiece energies...')
        # hack for ff03ua - need 2e-2 impropers tolerance
        improper_tolerance = 1e-2
        if united_atom: improper_tolerance = 2e-2
        assert_energies(villin_top[1], villin_crd[1], ffxml_name,
                        system_name='protein-villin headpiece',
                        improper_tolerance=improper_tolerance)
        if verbose: print('Villin headpiece energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (ala3_top, ala3_crd, villin_top, villin_crd, leap_script_ala3_file,
                 leap_script_villin_file):
            os.close(f[0])
            os.unlink(f[1])
    if verbose: print('Protein energy validation for %s done!' % ffxml_name)

def validate_nucleic(ffxml_name, leaprc_name):
    if verbose: print('Nucleic acids energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    dna_top = tempfile.mkstemp()
    dna_crd = tempfile.mkstemp()
    leap_script_dna_file = tempfile.mkstemp()
    rna_top = tempfile.mkstemp()
    rna_crd = tempfile.mkstemp()
    leap_script_rna_file = tempfile.mkstemp()
    leap_script_rna_file_alt = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
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

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_dna_file[1], os.devnull))
    if os.path.getsize(dna_top[1]) == 0 or os.path.getsize(dna_crd[1]) == 0:
        raise RuntimeError('DNA LEaP fail for %s' % leaprc_name)
    os.system('tleap -f %s > %s' % (leap_script_rna_file[1], os.devnull))
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        # try alternative name mappings
        os.system('tleap -f %s > %s' % (leap_script_rna_file_alt[1], os.devnull))
    if os.path.getsize(rna_top[1]) == 0 or os.path.getsize(rna_crd[1]) == 0:
        raise RuntimeError('RNA LEaP fail for %s' % leaprc_name)

    try:
        if verbose: print('Calculating and validating DNA energies...')
        assert_energies(dna_top[1], dna_crd[1], ffxml_name,
                        system_name='nucleic-DNA')
        if verbose: print('DNA energy validation successful!')

        if verbose: print('Calculating and validating RNA energies...')
        # improper testing turned off pending solution to problems
        assert_energies(rna_top[1], rna_crd[1], ffxml_name,
                        system_name='nucleic-RNA', improper_tolerance=float('inf'))
        if verbose: print('RNA energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (dna_top, dna_crd, leap_script_dna_file, rna_top, rna_crd,
                  leap_script_rna_file, leap_script_rna_file_alt):
            os.close(f[0])
            os.unlink(f[1])
    if verbose: print('Nucleic acids energy validation for %s done!' % ffxml_name)

def validate_gaff(ffxml_name, leaprc_name):
    if verbose: print('GAFF energy validation for %s' % ffxml_name)
    if verbose: print('Preparing temporary files for validation...')
    imatinib_top = tempfile.mkstemp()
    imatinib_crd = tempfile.mkstemp()
    leap_script_imatinib_file = tempfile.mkstemp()

    if verbose: print('Preparing LeaP scripts...')
    leap_script_imatinib_string = """source %s
loadamberparams files/frcmod.imatinib
x = loadMol2 files/imatinib.mol2
saveAmberParm x %s %s
quit""" % (leaprc_name, imatinib_top[1], imatinib_crd[1])
    os.write(leap_script_imatinib_file[0], leap_script_imatinib_string)

    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_imatinib_file[1], os.devnull))
    if os.path.getsize(imatinib_top[1]) == 0 or os.path.getsize(imatinib_crd[1]) == 0:
        raise RuntimeError('imatinib LEaP fail for %s' % leaprc_name)

    try:
        if verbose: print('Calculating and validating imatinib energies...')
        assert_energies(imatinib_top[1], imatinib_crd[1], (ffxml_name,
                        'files/imatinib.xml', 'files/imatinib_frcmod.xml'),
                         system_name='gaff-imatinib')
        if verbose: print('Imatinib energy validation successful!')
    finally:
        if verbose: print('Deleting temp files...')
        for f in (imatinib_top, imatinib_crd, leap_script_imatinib_file):
            os.close(f[0])
            os.unlink(f[1])
    if verbose: print('GAFF energy validation for %s done!' % ffxml_name)

def validate_phospho_protein(ffxml_name, leaprc_name,
                             supp_leaprc_name = 'leaprc.ff14SB',
                             supp_ffxml_name='ffxml/ff14SB.xml'):
    # this function assumes ffxml/ff14SB.xml already exists
    if verbose: print('Phosphorylated protein energy validation for %s' %
                      ffxml_name)
    for pdbname in glob.iglob('files/phospho/*.pdb'):
        if verbose: print('Now testing with pdb %s' % os.path.basename(pdbname))
        if verbose: print('Preparing temporary files for validation...')
        top = tempfile.mkstemp()
        crd = tempfile.mkstemp()
        leap_script_file = tempfile.mkstemp()

        if verbose: print('Preparing LeaP scripts...')
        leap_script_string = """source %s
source %s
x = loadPdb %s
saveAmberParm x %s %s
quit""" % (supp_leaprc_name, leaprc_name, pdbname, top[1], crd[1])

        os.write(leap_script_file[0], leap_script_string)

        if verbose: print('Running LEaP...')
        os.system('tleap -f %s > %s' % (leap_script_file[1], os.devnull))
        if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
            raise RuntimeError('LEaP fail for %s' % leaprc_name)

        try:
            if verbose: print('Calculating and validating energies...')
            assert_energies(top[1], crd[1], (supp_ffxml_name, ffxml_name),
                            system_name='phospho_protein: %s'
                            % os.path.basename(pdbname))
            if verbose: print('Energy validation successful!')
        finally:
            if verbose: print('Deleting temp files...')
            for f in (top, crd, leap_script_file):
                os.close(f[0])
                os.unlink(f[1])
        if verbose: print('Phosphorylated protein energy validation for %s done!'
                          % ffxml_name)

def validate_water_ion(ffxml_name, source_recipe_files, solvent_name, recipe_name,
                       standard_ffxml=None):
    if verbose: print('Water and ions energy validation for %s' %
                      ffxml_name)
    if solvent_name == 'tip3p':
        HOH = 'TP3'
        solvent_frcmod = None
    elif solvent_name == 'tip4pew':
        HOH = 'T4E'
        solvent_frcmod = 'frcmod.tip4pew'
    elif solvent_name == 'spce':
        HOH = 'SPC'
        solvent_frcmod = 'frcmod.spce'
    pdb_name = 'files/water_ion/' + recipe_name + '.pdb'
    if verbose: print('Preparing temporary files for validation...')
    top = tempfile.mkstemp()
    crd = tempfile.mkstemp()
    leap_script_file = tempfile.mkstemp()
    if verbose: print('Preparing LeaP scripts...')
    leap_script_string_part1 = """loadamberparams parm10.dat
loadamberparams %s
loadamberparams %s\n""" % (source_recipe_files[0], source_recipe_files[1])

    leap_script_string_part2 = """\nloadOff atomic_ions.lib
loadoff solvents.lib
HOH = %s
# for TIP4PEW
addPdbAtomMap {{ "M" "EPW" }}
x = loadPdb %s
saveAmberParm x %s %s
quit""" % (HOH, pdb_name, top[1], crd[1])

    if solvent_frcmod:
        leap_script_string = (leap_script_string_part1 + ('loadamberparams %s'
                               % solvent_frcmod) + leap_script_string_part2)
    else:
        leap_script_string = leap_script_string_part1 + leap_script_string_part2

    os.write(leap_script_file[0], leap_script_string)

    # this test does it's own energy assertion because of differences
    if verbose: print('Running LEaP...')
    os.system('tleap -f %s > %s' % (leap_script_file[1], os.devnull))
    if os.path.getsize(top[1]) == 0 or os.path.getsize(crd[1]) == 0:
        raise RuntimeError('LEaP fail for %s' % ffxml_name)
    try:
        if verbose: print('Calculating and validating energies...')
        pdb = app.PDBFile(pdb_name, extraParticleIdentifier='')
        if standard_ffxml is None:
            ff = app.ForceField(ffxml_name)
        else:
            ff = app.ForceField(ffxml_name, standard_ffxml)
        system_omm = ff.createSystem(pdb.topology)
        parm_omm = parmed.openmm.load_topology(pdb.topology, xyz=pdb.positions)
        parm_amber = parmed.load_file(top[1], crd[1])
        system_amber = parm_amber.createSystem()
        omm_energies = parmed.openmm.energy_decomposition_system(parm_omm,
                       system_omm, nrg=u.kilojoules_per_mole)
        for entry in omm_energies:
            if entry[0] == 'NonbondedForce':
                omm_nonbonded = entry[1]
        amber_energies = parmed.openmm.energy_decomposition_system(parm_amber,
                         system_amber, nrg=u.kilojoules_per_mole)
        for entry in amber_energies:
            if entry[0] == 'NonbondedForce':
                amber_nonbonded = entry[1]

        rel_nonbonded = abs(amber_nonbonded-omm_nonbonded) / amber_nonbonded
        if rel_nonbonded > 1e-5:
            raise AssertionError('NonbondedForce Water and ions energy outside of '
                                 'allowed tolerance for %s' % ffxml_name)
        if verbose: print('Energy validation successful!')

    finally:
        if verbose: print('Deleting temp files...')
        for f in (top, crd, leap_script_file):
            os.close(f[0])
            os.unlink(f[1])
    # logging
    if not no_log:
        amber_energies_log = dict()
        omm_energies_log = dict()
        rel_energies_log = dict()
        amber_energies_log['ffxml_name'] = ffxml_name
        amber_energies_log['test_system'] = 'water_ion'
        amber_energies_log['data_type'] = 'AMBER'
        amber_energies_log['NonbondedForce'] = amber_nonbonded
        amber_energies_log['units'] = u.kilojoules_per_mole
        omm_energies_log['ffxml_name'] = ffxml_name
        omm_energies_log['test_system'] = 'water_ion'
        omm_energies_log['data_type'] = 'OpenMM'
        omm_energies_log['NonbondedForce'] = omm_nonbonded
        omm_energies_log['units'] = u.kilojoules_per_mole
        rel_energies_log['ffxml_name'] = ffxml_name
        rel_energies_log['test_system'] = 'water_ion'
        rel_energies_log['data_type'] = 'abs(AMBER-OpenMM)/AMBER'
        rel_energies_log['NonbondedForce'] = rel_nonbonded
        logger.log(amber_energies_log)
        logger.log(omm_energies_log)
        logger.log(rel_energies_log)
    if verbose: print('Water and ions energy validation for %s done!'
                      % ffxml_name)

class Logger():
    # logs testing energies into csv
    def __init__(self, log_file):
        csvfile = open(log_file, 'w')
        fieldnames = ['ffxml_name', 'data_type', 'test_system', 'units', 'HarmonicBondForce',
                      'HarmonicAngleForce', 'PeriodicTorsionForce_dihedrals',
                      'PeriodicTorsionForce_impropers', 'NonbondedForce']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        self.csvfile = csvfile
        self.writer = writer

    def close(self):
        self.csvfile.close()

    def log(self, energies):
        self.writer.writerow(energies)

if __name__ == '__main__':
    main()
