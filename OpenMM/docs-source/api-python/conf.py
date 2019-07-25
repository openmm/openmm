# -*- coding: utf-8 -*-

import sys
import os
import simtk.openmm.version

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.ifconfig', 'sphinx.ext.autosummary',
              'sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'process-docstring',
              'sphinxcontrib.lunrsearch']

autosummary_generate = True
autodoc_default_flags = ['members', 'inherited-members']
autodoc_member_order = 'bysource'

source_suffix = '.rst'
master_doc = 'index'

project = u'OpenMM'
copyright = u'2015, Stanford University and the Authors'

version = simtk.openmm.version.short_version
release = simtk.openmm.version.full_version

exclude_patterns = ['_build', '_templates']
html_static_path = ['_static']
templates_path = ['_templates']

pygments_style = 'sphinx'

html_theme = "alabaster"
html_theme_options = {
    'description': "High performance molecular simulation on GPUs",
    'github_button': False,
    # 'github_user': 'pandegroup',
    # 'github_repo': 'openmm',
    'logo_name': False,
    'logo': 'logo.png',
}
html_sidebars = {
    '**': [
        'about.html',
        'searchbox.html',
        'navigation.html',
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
