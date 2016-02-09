from parmed.charmm import CharmmParameterSet
from parmed import openmm
import glob
import yaml

# top and par files for Charmm36 and Charmm22
charmm_files = glob.glob('charmm/toppar/*_all3*')
# stream files
charmm_files.extend(glob.glob('charmm/toppar/*.str'))
# charmm_files.extend((glob.glob('charmm/toppar/stream/prot/*.str')))
# charmm_files.extend(glob.glob('charmm/toppar/stream/carb/*.str'))
# charmm_files.extend(glob.glob('charmm/toppar/stream/lipid/*.str'))
# charmm_files.extend(glob.glob('charmm/toppar/stream/na/*.str'))
# charmm_files.extend(glob.glob('charmm/toppar/stream/misc/*.str'))

# convert top and par files to ffxml
#for f in charmm_files:
#    param = CharmmParameterSet(f)
#    param_omm = openmm.OpenMMParameterSet.from_parameterset(param)
#    param_omm.write('ffxml/%s.xml' % f[14:-4])

# generate recommended combination for charmm36
params = CharmmParameterSet(*charmm_files)
#charmm_files.remove('charmm/toppar/stream/lipid/toppar_all36_lipid_list.str')
params_omm = openmm.OpenMMParameterSet.from_parameterset(params)
params_omm.write('ffxml/charmm36.xml')
