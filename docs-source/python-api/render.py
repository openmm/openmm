from os.path import dirname, join, splitext
from glob import glob
import inspect

import jinja2
import simtk.openmm as mm
from simtk.openmm import app

def fullname(klass):
    return klass.__module__ + '.' + klass.__name__


#        'integrators': [],
#        'library_extras': [],
#         'forces': [],

def library_template_variables():
    data = {
        'integrators': [],
        'library_extras': [],
        'forces': [],
    }

    mm_klasses = inspect.getmembers(mm, predicate=inspect.isclass)

    # gather all Force subclasses
    for name, klass in mm_klasses:
        if issubclass(klass, mm.Force):
            data['forces'].append(fullname(klass))

    # gather all Integrator subclasses
    for _, klass in mm_klasses:
        if issubclass(klass, mm.Integrator):
            data['integrators'].append(fullname(klass))

    # gather all extra subclasses in simtk.openmm.openmm
    exclude = ['simtk.openmm.openmm.Platform', 'simtk.openmm.openmm.Context',
              'simtk.openmm.openmm.System', 'simtk.openmm.openmm.State']
    exclude.extend(data['forces'])
    exclude.extend(data['integrators'])
    exclude.extend([
        'simtk.openmm.openmm.SwigPyIterator',
        'simtk.openmm.openmm.OpenMMException'])

    for _, klass in mm_klasses:
        full = fullname(klass)
        if full not in exclude and not klass.__name__[0].islower():
            data['library_extras'].append(full)

    return data


def app_template_variables():
    data = {
        'reporters': [],
        'fileclasses': [],
        'app_extras': [],
    }

    app_klasses = inspect.getmembers(app, predicate=inspect.isclass)

    # gather all Reporters
    for name, klass in app_klasses:
        if name.endswith('Reporter'):
            data['reporters'].append(fullname(klass))

    # gather all classes with "File" in the name
    for name, klass in app_klasses:
        if 'File' in name:
            data['fileclasses'].append(fullname(klass))

    # gather all extra subclasses in simtk.openmm.app
    exclude = ['simtk.openmm.app.topology.Topology',
               'simtk.openmm.app.modeller.Modeller',
               'simtk.openmm.app.forcefield.ForceField',
               'simtk.openmm.app.simulation.Simulation']
    exclude.extend(data['reporters'])
    exclude.extend(data['fileclasses'])

    for _, klass in app_klasses:
        full = fullname(klass)
        if full not in exclude and not klass.__name__[0].islower():
            data['app_extras'].append(full)

    return data


def main():
    here = dirname(__file__)
    templateLoader = jinja2.FileSystemLoader(here)
    templateEnv = jinja2.Environment(loader=templateLoader)
    data = library_template_variables()
    data.update(app_template_variables())

    for template_fn in glob(join(here, '*.jinja2')):
        output_fn = splitext(template_fn)[0]
        print('Rendering %s to %s...' % (template_fn, output_fn))

        template = templateEnv.get_template(template_fn)
        output_text = template.render(data)
        with open(output_fn, 'w') as f:
            f.write(output_text)


if __name__ == '__main__':
    main()
