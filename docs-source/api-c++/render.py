from __future__ import print_function

import os
import re
import sys
from functools import reduce
from glob import glob
from os.path import basename, dirname, join, splitext
from subprocess import run

import jinja2
import lxml.etree as ET


def load_doxygen_xml(doxygen_xml):
    files = [
        os.path.join(doxygen_xml, f)
        for f in os.listdir(doxygen_xml)
        if f.lower().endswith(".xml") and not f.startswith("._")
    ]
    if len(files) == 0:
        raise err

    document = ET.ElementTree(ET.Element("root")).getroot()
    for file in files:
        root = ET.parse(file).getroot()
        for node in root:
            document.append(node)
    return document


def subclasses(root, parent):
    parent_el = root.xpath('.//compounddef/compoundname[text()="%s"]/..' % parent)
    if len(parent_el) == 1:
        parent_id = parent_el[0].get("id")
    else:
        raise ValueError("Can't find %s" % parent)

    xp_query = (
        './/compounddef/basecompoundref[@refid="%s"]/../compoundname'
    ) % parent_id
    return [parent] + [n.text.strip() for n in root.xpath(xp_query)]


def allclasses(root):
    xp_query = './/compounddef[@kind="class" and @prot="public"]/compoundname'
    return [e.text for e in root.xpath(xp_query)]


def brief_description(root, compoundname):
    xp_query = (
        f'.//compounddef/compoundname[text()="{compoundname}"]/../briefdescription'
    )
    nodes = ["".join(node.itertext()) for node in root.xpath(xp_query)]
    return "".join(nodes).strip()


def detailed_description(root, compoundname):
    xp_query = (
        f'.//compounddef/compoundname[text()="{compoundname}"]/../detaileddescription'
    )
    nodes = ["".join(node.itertext()) for node in root.xpath(xp_query)]
    return "".join(nodes).strip()


def summarize(text):
    # If there's a blank line, then we can assume the first sentence /
    # paragraph has ended, so anything after shouldn't be part of the
    # summary
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if not line.strip():
            lines = lines[:i]
            break

    # Try to find the "first sentence", which may span multiple lines
    m = re.search(r"^([A-Z].*?\.)(?:\s|$)", " ".join(lines).strip())
    if m:
        return m.group(1).strip()
    elif lines:
        return lines[0].strip()
    else:
        return ""


def class_summaries(root):
    outdict = {}
    for compoundname in allclasses(root):
        brief_description_text = brief_description(root, compoundname)
        if brief_description_text:
            outdict[compoundname] = brief_description_text
        else:
            outdict[compoundname] = summarize(detailed_description(root, compoundname))

    return outdict


def template_data(root):
    data = {
        "core": (
            "OpenMM::System",
            "OpenMM::Context",
            "OpenMM::State",
            "OpenMM::Platform",
        ),
        "forces": sorted(subclasses(root, "OpenMM::Force")),
        "integrators": sorted(subclasses(root, "OpenMM::Integrator")),
    }
    data["extras"] = sorted(
        set(allclasses(root)) - reduce(set.union, map(set, data.values()))
    )

    data["class_summaries"] = class_summaries(root)
    return data


def main():
    if len(sys.argv) == 1:
        print("usage: %s <doxygen_xml_path>" % sys.argv[0], file=sys.stderr)
        exit(1)
    doxygen_xml_path = sys.argv[1]

    root = load_doxygen_xml(doxygen_xml_path)
    data = template_data(root)

    here = dirname(__file__)
    templateLoader = jinja2.FileSystemLoader(here)
    templateEnv = jinja2.Environment(loader=templateLoader)

    for template_fn in map(basename, glob(join(here, "*.jinja2"))):
        output_fn = splitext(template_fn)[0]
        print("Rendering %s to %s..." % (template_fn, output_fn))

        template = templateEnv.get_template(template_fn)
        output_text = template.render(data)
        with open(output_fn, "w") as f:
            f.write(output_text)

    run(
        [
            "python3",
            "breathe-apidoc.py",
            "--output-dir=generated",
            "--generate=class",
            "--members",
            "--force",
            "--brief-titles",
            "--rename-output",
            "--flat-output",
            "--quiet",
            doxygen_xml_path,
        ],
        check=True,
    )


if __name__ == "__main__":
    main()
