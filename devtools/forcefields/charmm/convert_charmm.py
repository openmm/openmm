from parmed.charmm import CharmmParameterSet
from parmed import openmm
import glob

# top and par files for Charmm36 and Charmm22
charmm_files = glob.glob('charmm/toppar/*_all3*')
charmm_files.extend(glob.glob('charmm/toppar/*_all2*'))

# stream files
stream_files = glob.glob('charmm/toppar/stream/prot/*.str')
stream_files.extend(glob.glob('charmm/toppar/stream/carb/*.str'))
stream_files.extend(glob.glob('charmm/toppar/stream/lipid/*.str'))
stream_files.extend(glob.glob('charmm/toppar/stream/na/*.str'))
stream_files.extend(glob.glob('charmm/toppar/stream/misc/*.str'))

# convert top and par files to ffxml
for file in charmm_files:
    param = CharmmParameterSet(file)
    param_omm = openmm.OpenMMParameterSet.from_parameterset(param)
    param_omm.write('ffxml/%s.xml' % file[14:-4])

# convert stream files to ffxml
for file in stream_files:
    param = CharmmParameterSet(file)
    param_omm = openmm.OpenMMParameterSet.from_parameterset(param)
    param_omm.write('ffxml/%s.xml' % file[14:-4])