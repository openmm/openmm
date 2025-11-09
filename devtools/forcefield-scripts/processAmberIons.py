import re
import os
import sys
import glob
import datetime
import subprocess

# from openmm import unit

kilocalorie_to_kilojoule =  4.184 # 1.0*unit.kilocalorie/unit.kilojoule
angstrom_to_nanometer = 0.1 # 1.0*unit.angstrom/unit.nanometer


def md5sum(file):
    result = subprocess.check_output(['md5sum', file]).decode('utf-8').split()[0]
    return result


def amber_frcmod_ions_to_openmm_xml(
    frcmod_file, 
    atomic_ions_file, 
    xml_file, 
    setname=None, 
    include_extra_comments=True,
):
    with open(atomic_ions_file) as fh:
        lines = fh.readlines()

    ionname_to_residuename = {}
    ionname_to_charge = {}
    for i, line in enumerate(lines):
        if 'unit.atoms table' in line:
            name, type_, typex, resx, flags, seq, elmnt, chg = lines[i+1].split()
            name = name.replace('"', '')
            type_ = type_.replace('"', '')
            # print(name, type_, typex, resx, flags, seq, elmnt, chg)
            ionname_to_residuename[type_] = name
            ionname_to_charge[type_] = float(chg)
            
    ionname_to_residuename['Na+'] = 'NA' # this is a duplicate in atomic_ions.lib
    ionname_to_residuename['Cl-'] = 'CL' # this is a duplicate in atomic_ions.lib
    ionname_to_residuename['Cu2+'] = 'CU' # consistent with previous residue names
    ionname_to_residuename['Cu+'] = 'CU1'

    with open(frcmod_file) as fh:
        lines = fh.readlines()
        lines = [x.strip() for x in lines]
        
    mass_data = []
    nonbon_data = []
    comments_data = []
    do_mass = False
    do_nonbon = False

    for i, line in enumerate(lines):
        if line == 'MASS':
            do_mass = True
            do_nonbon = False
        elif line == 'NONBON':
            do_nonbon = True
            do_mass = False
        elif line == '':
            do_mass = False
            do_nonbon = False
        else:
            if do_mass:
                m = re.match(r'^(\S+)\s+(\S+)(.*)', line)
                if m:
                    ionname, mass, extra = m.groups()
                    mass_data.append((ionname, mass, extra))
                else:
                    raise ValueError(f'Could not parse line {i}: {line}')
            elif do_nonbon:
                m = re.match(r'^(\S+)\s+(\S+)\s+(\S+)(.*)', line)
                if m:
                    ionname, sigma, epsilon, extra = m.groups()
                    nonbon_data.append((ionname, float(sigma), float(epsilon), extra))
                else:
                    raise ValueError(f'Could not parse line {i}: {line}')
            else:
                comments_data.append(line)

    if not setname:
        setname = os.path.basename(frcmod_file)
        setname = setname.replace('frcmod.', '')

    date_generated = datetime.datetime.utcnow().isoformat()
    ambertools_version = re.sub(r'.*/ambertools-([^/]*)/.*', r'\1', frcmod_file)
    frcmod_file_short = re.sub(r'.*/dat/leap/', '', frcmod_file)
    atomic_ions_file_short = re.sub(r'.*/dat/leap/', '', atomic_ions_file)

    output_lines = []

    output_lines.append('<ForceField>')

    output_lines.append('  <Info>')
    for comment in comments_data:
        output_lines.append(f'    <!-- {comment} -->')
    output_lines.append(f'    <DateGenerated>{date_generated}</DateGenerated>')
    output_lines.append(f'    <Source Source="{frcmod_file_short}" md5hash="{md5sum(frcmod_file)}" sourcePackage="AmberTools" sourcePackageVersion="{ambertools_version}">{frcmod_file_short}</Source>')
    output_lines.append(f'    <Source Source="{atomic_ions_file_short}" md5hash="{md5sum(atomic_ions_file)}" sourcePackage="AmberTools" sourcePackageVersion="{ambertools_version}">{atomic_ions_file_short}</Source>')
    output_lines.append('  </Info>')

    output_lines.append('  <AtomTypes>')
    for ionname, mass, extra in mass_data:
        element = re.sub(r'\d*[\-\+]*', '', ionname)
        line = f'    <Type class="{setname}-{ionname}" element="{element}" mass="{mass}" name="{setname}-{ionname}"/>'
        if extra and include_extra_comments:
            line += f'<!-- {extra} -->'
        output_lines.append(line)
    output_lines.append('  </AtomTypes>')

    output_lines.append('  <Residues>')
    # NOTE this is a hack to sort residues in a similar way to the older converted xml files, for easier diff
    # NOTE it would be better to keep them in the order in which ions are listed
    for ionname in sorted([x[0] for x in mass_data], key=str.casefold):
        if not ionname in ionname_to_residuename:
            print(f'Warning: ion name {ionname} not found in list of atomic ions')
            continue
        residuename = ionname_to_residuename[ionname]
        charge = ionname_to_charge[ionname]
        output_lines.append(f'    <Residue name="{residuename}">')
        output_lines.append(f'      <Atom charge="{charge:.1f}" name="{residuename}" type="{setname}-{ionname}"/>')
        output_lines.append(f'    </Residue>')
    output_lines.append('  </Residues>')

    output_lines.append('  <NonbondedForce coulomb14scale="0.8333333333333334" lj14scale="0.5">')
    output_lines.append('    <UseAttributeFromResidue name="charge"/>')
    for ionname, sigma, epsilon, extra in nonbon_data:
        epsilon_converted = epsilon*kilocalorie_to_kilojoule
        sigma_converted = (sigma * 2 / (2**(1/6)))*angstrom_to_nanometer
        line = f'    <Atom epsilon="{epsilon_converted}" sigma="{sigma_converted}" type="{setname}-{ionname}"/>'
        if extra and include_extra_comments:
            line += f'<!-- {extra} -->'
        output_lines.append(line)
    output_lines.append('  </NonbondedForce>')
    output_lines.append('</ForceField>')
    
    with open(xml_file, 'w') as fh:
        for line in output_lines:
            fh.write(line+'\n')
            
    print(f'Converted {len(mass_data)} ions from {frcmod_file} to {xml_file}')

if len(sys.argv) < 2:
    print('processAmberIons.py script to convert all Amber/AmberTools frcmod.ions* files to equivalent openmm xml files')
    print('Usage: python ./devtools/forcefield-scripts/processAmberIons.py /path/to/ambertools')
    sys.exit()

ambertools_path = sys.argv[1]
ambertools_dat_leap_path = os.path.join(ambertools_path, 'dat/leap/')
frcmod_files = glob.glob(ambertools_dat_leap_path + 'parm/frcmod.ions*')
atomic_ions_file = ambertools_dat_leap_path + 'lib/atomic_ions.lib'

output_path = 'wrappers/python/openmm/app/data/ions'
os.makedirs(output_path, exist_ok=True)

for frcmod_file in frcmod_files:
    if '1264' in frcmod_file:
        print(f'Skipped {frcmod_file} because it has 12-6-4 parameters, only 12-6 are supported')
        continue
    setname = os.path.basename(frcmod_file)
    setname = setname.replace('frcmod.', '')
    xml_file = f'{output_path}/{setname}.xml'
    amber_frcmod_ions_to_openmm_xml(
        frcmod_file, 
        atomic_ions_file, 
        xml_file,
        setname=setname)
