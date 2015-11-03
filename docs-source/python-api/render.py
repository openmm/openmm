from os.path import dirname, join, splitext
from glob import glob
import inspect

import jinja2
import simtk.openmm as mm
from simtk.openmm import app

def fullname(klass):
    return klass.__module__ + '.' + klass.__name__


def template_variables():
    data = {
        'reporters': [],
        'forces': [],
        'integrators': [],
        'library_extras': []
    }


    app_klasses = inspect.getmembers(app, predicate=inspect.isclass)
    mm_klasses = inspect.getmembers(mm, predicate=inspect.isclass)

    # gather all Reporters
    for name, klass in app_klasses:
        if name.endswith('Reporter'):
            data['reporters'].append(fullname(klass))

    # gather all Force subclasses
    for name, klass in mm_klasses:
        if issubclass(klass, mm.Force):
            data['forces'].append(fullname(klass))

    # gather all Integrator subclasses
    for _, klass in mm_klasses:
        if issubclass(klass, mm.Integrator):
            data['integrators'].append(fullname(klass))

    # gather all extra subclasses
    for _, klass in mm_klasses:
        full = fullname(klass)

        if full in data['forces']:
            continue
        if full in data['integrators']:
            continue
        if full in ('simtk.openmm.openmm.Platform', 'simtk.openmm.openmm.Context',
                    'simtk.openmm.openmm.System'):
            continue
        if klass.__name__[0].islower():
            continue
        if klass.__name__ in ['SwigPyIterator', 'OpenMMException']:
            continue

        data['library_extras'].append(full)

    return data


def main():
    here = dirname(__file__)
    templateLoader = jinja2.FileSystemLoader(here)
    templateEnv = jinja2.Environment(loader=templateLoader)
    data = template_variables()

    for template_fn in glob(join(here, '*.jinja2')):
        output_fn = splitext(template_fn)[0]
        print('Rendering %s to %s...' % (template_fn, output_fn))

        template = templateEnv.get_template(template_fn)
        output_text = template.render(data)
        with open(output_fn, 'w') as f:
            f.write(output_text)


if __name__ == '__main__':
    main()
