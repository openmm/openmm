# -*- coding: utf-8 -*-

import sys
import os
import simtk.openmm.version

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.ifconfig', 'sphinx.ext.autosummary',
              'sphinx.ext.autodoc', 'numpydoc', 'sphinx.ext.intersphinx',]

autosummary_generate = True
autodoc_default_flags = ['members', 'inherited-members']

_python_doc_base = 'http://docs.python.org/2.7'
intersphinx_mapping = {
    _python_doc_base: None,
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'OpenMM'
copyright = u'2015, Stanford University and the Authors'

# The short X.Y version.
version = simtk.openmm.version.short_version
# The full version, including alpha/beta/rc tags.
release = simtk.openmm.version.full_version

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build']


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = "alabaster"
html_theme_options = {
    'logo': 'OPENMMLogo.png',
    'description': 'A high performance GPU molecular simulation toolkit',
    'github_user': 'pandegroup',
    'github_repo': 'openmm',
    'travis_button': True,
}


# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = []

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'searchbox.html',
    ]
}


autodoc_member_order = 'bysource'

# stackoverflow.com/questions/12206334
numpydoc_show_class_members = False
