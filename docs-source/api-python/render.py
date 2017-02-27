"""
The function of this script is to render the Jinja2 templates in the current
directory into input files for sphinx. It introspects the OpenMM Python module
to find all of the classes and formats them for inclusion into the templates.
"""
from os.path import dirname, join, splitext, basename
from glob import glob
import inspect

import jinja2
import simtk.openmm
import simtk.openmm.app



def fullname(klass):
    return klass.__module__ + '.' + klass.__name__


def library_template_variables():
    """Create the data structure available to the Jinja2 renderer when
    filling in the templates.

    This function extracts all of classes in ``simtk.openmm.openmm`` and returns
    a dictionary with them grouped into three lists, the integrators, the forces,
    and the remainder (library_extras).

    A couple core classes are skipped, because they're included manually in the
    template.
    """
    data = {
        'integrators': [],
        'forces': [],
        'library_extras': [],
        'units': [],
    }

    mm_klasses = inspect.getmembers(simtk.openmm, predicate=inspect.isclass)


    # gather all Force subclasses
    for name, klass in mm_klasses:
        if issubclass(klass, simtk.openmm.openmm.Force):
            data['forces'].append(fullname(klass))

    # gather all Integrator subclasses
    for _, klass in mm_klasses:
        if issubclass(klass, simtk.openmm.openmm.Integrator):
            data['integrators'].append(fullname(klass))

    # gather all extra subclasses in simtk.openmm.openmm

    # core classes that are already included in library.rst.jinja2
    exclude = ['simtk.openmm.openmm.Platform', 'simtk.openmm.openmm.Context',
              'simtk.openmm.openmm.System', 'simtk.openmm.openmm.State']
    exclude.extend(data['forces'])
    exclude.extend(data['integrators'])

    # these classes are useless and not worth documenting.
    exclude.extend([
        'simtk.openmm.openmm.SwigPyIterator',
        'simtk.openmm.openmm.OpenMMException'])

    for _, klass in mm_klasses:
        full = fullname(klass)
        if full not in exclude and not klass.__name__[0].islower():
            data['library_extras'].append(full)

    # gather units related classes
    unit_klasses = inspect.getmembers(simtk.unit, predicate=inspect.isclass)
    for name, klass in unit_klasses:
        data['units'].append(fullname(klass))

    return data


def app_template_variables():
    """Create the data structure available to the Jinja2 renderer when
    filling in the templates.

    This function extracts all of classes in ``simtk.openmm.app`` and returns
    a dictionary with them grouped into three lists, the reporters, the
    classes with the word "File" in the name, and the remainder.

    Four classes are skipped (see exclude), because they're included manually
    in the template.
    """
    data = {
        'reporters': [],
        'fileclasses': [],
        'app_extras': [],
    }
    app_klasses = inspect.getmembers(simtk.openmm.app, predicate=inspect.isclass)

    # gather all Reporters
    for name, klass in app_klasses:
        if name.endswith('Reporter'):
            data['reporters'].append(fullname(klass))

    # gather all classes with "File" in the name
    for name, klass in app_klasses:
        if 'File' in name or 'CharmmParameterSet' in name:
            data['fileclasses'].append(fullname(klass))

    # gather all extra subclasses in simtk.openmm.app
    exclude = ['simtk.openmm.app.topology.Topology',
               'simtk.openmm.app.topology.Chain',
               'simtk.openmm.app.topology.Residue',
               'simtk.openmm.app.topology.Atom',
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

    for template_fn in map(basename, glob(join(here, '*.jinja2'))):
        output_fn = splitext(template_fn)[0]
        print('Rendering %s to %s...' % (template_fn, output_fn))

        template = templateEnv.get_template(template_fn)
        output_text = template.render(data)
        with open(output_fn, 'w') as f:
            f.write(output_text)


if __name__ == '__main__':
    main()
