from parmed.charmm import CharmmParameterSet
from parmed import openmm
import glob
import yaml

data = yaml.load(open('charmm36.yaml'))
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
files_to_convert = set(charmm_files) - exclude_files


# # top and par files for Charmm36
# charmm_files = glob.glob('charmm/toppar/*_all36*')
# # stream files
# charmm_files.extend(glob.glob('charmm/toppar/*.str'))
# charmm_files.extend((glob.glob('charmm/toppar/stream/prot/*.str')))
# charmm_files.extend(glob.glob('charmm/toppar/stream/carb/*.str'))
# charmm_files.extend(set(glob.glob('charmm/toppar/stream/lipid/*.str')) - exclude_files)
# charmm_files.extend(set(glob.glob('charmm/toppar/stream/na/*.str')) - exclude_files)
# charmm_files.extend(glob.glob('charmm/toppar/stream/misc/*.str'))

#generate recommended combination for charmm36
params = CharmmParameterSet(*files_to_convert)
params_omm = openmm.OpenMMParameterSet.from_parameterset(params)
params_omm.write('ffxml/charmm36.xml')

# for file in (set(glob.glob('charmm/toppar/stream/na/*.str')) - exclude_files):
#     print(file)
#     charmm_files.extend([file])
#     #print(charmm_files)
#     params = CharmmParameterSet(*charmm_files)
#     params_omm = openmm.OpenMMParameterSet.from_parameterset(params)
#     params_omm.write('test.xml')