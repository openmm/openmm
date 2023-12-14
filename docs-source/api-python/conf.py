# -*- coding: utf-8 -*-

import os
import sys

import openmm.version

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "process-docstring",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "inherited-members": True,
    "member-order": "bysource",
}

source_suffix = ".rst"
master_doc = "index"

project = u"OpenMM Python API"
copyright = u"2015, Stanford University and the Authors"

version = openmm.version.short_version
release = openmm.version.full_version

exclude_patterns = ["_build", "_templates"]
html_static_path = ["_static"]
templates_path = ["_templates"]

pygments_style = "sphinx"

html_theme = "alabaster"
html_theme_options = {
    "github_button": False,
    "github_user": "openmm",
    "github_repo": "openmm",
    "logo_name": True,
    "logo": "logo.png",
    "extra_nav_links": [
        {
            "title": "OpenMM.org",
            "uri": "https://openmm.org",
            "relative": False,
        },
        {
            "title": "User's Manual",
            "uri": "../userguide/",
            "relative": True,
        },
        {
            "title": "Developer Guide",
            "uri": "../developerguide/",
            "relative": True,
        },
        {
            "title": "C++ API reference",
            "uri": "../api-c++/",
            "relative": True,
        },
        {
            "title": "Cookbook & Tutorials",
            "uri": "https://openmm.github.io/openmm-cookbook/",
            "relative": False,
        },
        {
            "title": "GitHub",
            "uri": "https://github.com/openmm",
            "relative": False,
        },
    ],
    "show_relbar_bottom": True,
}
html_sidebars = {
    "**": [
        "about.html",
        "searchbox.html",
        "navigation.html",
    ]
}

# Napoleon settings
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
