
from __future__ import print_function
import re
import sys
import textwrap
from inspect import cleandoc
from numpydoc.docscrape import NumpyDocString


# Doxygen does a bad job of generating documentation based on docstrings.
# This script is run as a filter
# on each file, and converts the docstrings into Doxygen style comments
# so we get better documentation.

input = open(sys.argv[1])

while True:
    line = input.readline()
    if len(line) == 0:
        break
    stripped = line.lstrip()
    if stripped.startswith('def') or stripped.startswith('class'):
        prefix = line[:len(line)-len(stripped)]
        split = stripped.split()
        if split[0] == 'class' and split[1][0].islower():
            # Classes that start with a lowercase letter were defined by SWIG.  We want to hide them.
            print("%s## @private" % prefix, end='')
        if split[1][0] == '_' and split[1][1] != '_':
            # Names starting with a single _ are assumed to be private.
            print("%s## @private" % prefix, end='')

        # We're at the start of a class or function definition.  Find all lines that contain the declaration.

        declaration = line
        while len(line) > 0 and line.find(':') == -1:
            line = input.readline()
            declaration += line

        # Now look for a docstring.
        docstringlines = []
        line = input.readline()
        l = line.strip()
        if l.startswith('"""') or l.startswith("'''"):
            if l.count('"""') == 2 or l.count("'''") == 2:
                docstringlines.append(l[3:-3])
            else:
                docstringlines.append(l[3:])
                line = input.readline()
                l = line.strip()
                while l.find('"""') == -1 and l.find("'''") == -1:
                    docstringlines.append(line.rstrip())
                    line = input.readline()
                    l = line.strip()
                if line.rstrip().endswith('"""') or line.rstrip().endswith("'''"):
                    # last line
                    docstringlines.append(line.rstrip()[:-3])

        docstring = NumpyDocString(cleandoc('\n'.join(docstringlines)))
        # print(docstringlines)

        # Print out the docstring in Doxygen syntax, followed by the declaration.
        for line in docstring['Summary']:
            print('{prefix}## {line}'.format(prefix=prefix, line=line))
        if len(docstring['Extended Summary']) > 0:
            print('{prefix}##'.format(prefix=prefix))
            for line in docstring['Extended Summary']:
                print('{prefix}## {line}'.format(prefix=prefix, line=line.strip()))
        print('{prefix}##'.format(prefix=prefix))
        for name, type, descr in docstring['Parameters']:
            print('{prefix}## @param {name} ({type}) {descr}'.format(prefix=prefix, type=type, name=name, descr=' '.join(descr)))
        for name, type, descr in docstring['Returns']:
            if type == '':
                type = name
            print('{prefix}## @return ({type}) {descr}'.format(prefix=prefix, type=type, name=name, descr=''.join(descr)))

        print(declaration, end='')
        if len(docstringlines) == 0:
            print(line, end='')
    else:
        print(line, end='')
