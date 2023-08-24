import os
import sys

extensions = ["sphinx.ext.mathjax", "breathe"]

autosummary_generate = False
autodoc_member_order = "bysource"

breathe_projects = {
    "api-c++": "doxygen/xml",
}
breathe_default_project = "api-c++"

# Tell sphinx what the primary language being documented is.
primary_domain = "cpp"

# Tell sphinx what the pygments highlight language should be.
highlight_language = "cpp"

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
            "title": "Python API reference",
            "uri": "../api-python/",
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

doxygen_xml = "doxygen/xml"
