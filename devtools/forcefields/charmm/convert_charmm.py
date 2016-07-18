from parmed.charmm import CharmmParameterSet, CharmmPsfFile
from parmed import openmm
import glob
import yaml
from collections import OrderedDict
import hashlib
import os
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u

data = yaml.safe_load(open('charmm36.yaml'))
source_pack = data[0]['sourcePackage']
source_pack_ver = data[0]['sourcePackageVersion']

for entry in data[1:]:
    charmm_references = entry['References']
    source_files = entry['Source']


# files that should be excluded from conversion.
exclude_files = set(source_files['exclude'])

# charmm36 main top and par files
charmm_files = source_files['include']

# add stream files
for files in source_files['stream']:
    charmm_files.extend(glob.glob(files))

# exclude files from conversion
charmm_files = set(charmm_files) - exclude_files

provenance = OrderedDict()
source = provenance['Source'] = []
for fi in charmm_files:
    source.append(OrderedDict())
    source[-1]['Source'] = fi
    md5 = hashlib.md5()
    with open(fi) as f:
        md5.update(f.read())
    md5 = md5.hexdigest()
    source[-1]['md5hash'] = md5
    source[-1]['sourcePackage'] = source_pack
    source[-1]['sourcePackageVersion'] = source_pack_ver

references = provenance['Reference'] = []
for ff in charmm_references:
    for cite in charmm_references[ff]:
        references.append(OrderedDict())
        if type(cite) is dict:
            for key in cite.keys():
                citation = cite[key]
                references[-1]['Reference'] = citation
                references[-1]['forcefield'] = ff
                references[-1]['type'] = key
        else:
            citation = cite
            references[-1]['Reference'] = citation
            references[-1]['forcefield'] = ff


#generate recommended combination for charmm36
params = CharmmParameterSet(*charmm_files)
params_omm = openmm.OpenMMParameterSet.from_parameterset(params)
params_omm.write('ffxml/charmm36.xml', provenance=provenance)


