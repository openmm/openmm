import yaml
from parmed.charmm import CharmmParameterSet
from parmed import openmm

with open('charmm36.yaml', 'r') as f:
    yaml_content = yaml.load(f)

for xmlfile in yaml_content:
    if 'load' in yaml_content[xmlfile]:
        params = CharmmParameterSet.load_set(*yaml_content[xmlfile]['load'])
        params_omm = openmm.OpenMMParameterSet.from_parameterset(params)
        params_omm.write('ffxml/%s' % xmlfile)