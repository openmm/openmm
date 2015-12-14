from __future__ import print_function
import os
import sys
from functools import reduce
from os.path import basename, dirname, join, splitext
from glob import glob

import jinja2
import lxml.etree as ET


def load_doxygen_xml(doxygen_xml):
    files = [os.path.join(doxygen_xml, f)
             for f in os.listdir(doxygen_xml)
             if f.lower().endswith('.xml') and not f.startswith('._')]
    if len(files) == 0:
        raise err

    document = ET.ElementTree(ET.Element('root')).getroot()
    for file in files:
        root = ET.parse(file).getroot()
        for node in root:
            document.append(node)
    return document


def subclasses(root, parent):
    parent_el = root.xpath('.//compounddef/compoundname[text()="%s"]/..' % parent)
    if len(parent_el) == 1:
        parent_id = parent_el[0].get('id')
    else:
        raise ValueError("Can't find %s" % parent)

    xp_query = ('.//compounddef/basecompoundref[@refid="%s"]'
                '/../compoundname') % parent_id
    return [parent] + [n.text.strip() for n in root.xpath(xp_query)]


def allclasses(root):
    xp_query = './/compounddef[@kind="class" and @prot="public"]/compoundname'
    return [e.text for e in root.xpath(xp_query)]


def template_data(root):
    data = {
        'core': ('OpenMM::System', 'OpenMM::Context', 'OpenMM::State', 'OpenMM::Platform'),
        'forces': sorted(subclasses(root, 'OpenMM::Force')),
        'integrators': sorted(subclasses(root, 'OpenMM::Integrator')),

    }
    data['extras'] = sorted(set(allclasses(root)) -
                            reduce(set.union, map(set, data.values())))
    return data


def main():
    if len(sys.argv) == 1:
        print('usage: %s <doxygen_xml_path>' % sys.argv[0], file=sys.stderr)
        exit(1)
    doxygen_xml_path = sys.argv[1]

    root = load_doxygen_xml(doxygen_xml_path)
    data = template_data(root)

    here = dirname(__file__)
    templateLoader = jinja2.FileSystemLoader(here)
    templateEnv = jinja2.Environment(loader=templateLoader)

    for template_fn in map(basename, glob(join(here, '*.jinja2'))):
        output_fn = splitext(template_fn)[0]
        print('Rendering %s to %s...' % (template_fn, output_fn))

        template = templateEnv.get_template(template_fn)
        output_text = template.render(data)
        with open(output_fn, 'w') as f:
            f.write(output_text)


if __name__ == '__main__':
    main()
