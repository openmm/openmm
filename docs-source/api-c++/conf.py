import os
import sys

extensions = [
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinxcontrib.lunrsearch",
    "sphinxcontrib.autodoc_doxygen",
]

autosummary_generate = True
autodoc_member_order = "bysource"

source_suffix = ".rst"
master_doc = "index"

project = u"OpenMM C++ API"
copyright = u"2015, Stanford University and the Authors"

version = "@OPENMM_MAJOR_VERSION@.@OPENMM_MINOR_VERSION@"
release = "@OPENMM_MAJOR_VERSION@.@OPENMM_MINOR_VERSION@"

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
    "extra_nav_links": {
        "OpenMM.org": "https://openmm.org",
        "Developer's Guide": "http://docs.openmm.org/latest/developerguide/",
        "User's Guide": "http://docs.openmm.org/latest/userguide/",
        "Python API reference": "http://docs.openmm.org/latest/api-python/",
        "GitHub": "https://github.com/openmm",
    },
    "show_relbar_bottom": True,
}
html_sidebars = {
    "**": [
        "about.html",
        "lunrsearch.html",
        "navigation.html",
    ]
}

doxygen_xml = "doxygen/xml"
