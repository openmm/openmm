# -*- coding: utf-8 -*-

import sys
import os
import simtk.openmm.version

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.ifconfig', 'sphinx.ext.autosummary',
              'sphinx.ext.autodoc', 'numpydoc',]

autosummary_generate = True
autodoc_default_flags = ['members', 'inherited-members']

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
html_theme = "sphinx_rtd_theme"
html_logo = 'logo.png'
#html_favicon = None

html_static_path = ['_static']

autodoc_member_order = 'bysource'

# stackoverflow.com/questions/12206334
numpydoc_show_class_members = False
