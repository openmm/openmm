from __future__ import print_function
from numpydoc.docscrape import NumpyDocString
import sys
import textwrap
from inspect import cleandoc


# Doxygen does a bad job of generating documentation based on docstrings.  This script is run as a filter
# on each file, and converts the docstrings into Doxygen style comments so we get better documentation.

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
        stripped = line.lstrip()
        if stripped.startswith('"""') or stripped.startswith("'''"):
            line = stripped[3:]
            readingParameters = False
            if line.find('"""') != -1 or line.find("'''") != -1:
                docstringlines.append(line.rstrip()[:-3])
            else:
                while line.find('"""') == -1 and line.find("'''") == -1:
                    docstringlines.append(line)
                    line = input.readline()

        docstring = NumpyDocString(cleandoc(''.join(docstringlines)))

        # Print out the docstring in Doxygen syntax, followed by the declaration.
        for line in docstring['Summary']:
            print('{prefix}## {line}'.format(prefix=prefix, line=line))
        if len(docstring['Extended Summary']) > 0:
            print('{prefix}##'.format(prefix=prefix))
            for line in docstring['Extended Summary']:
                print('{prefix}## {line}'.format(prefix=prefix, line=line.strip()))
        print('{prefix}##'.format(prefix=prefix))
        for name, type, descr in docstring['Parameters']:
            print('{prefix}## @param {name} ({type}) {descr}'.format(prefix=prefix, type=type, name=name, descr=''.join(descr)))
        for name, type, descr in docstring['Returns']:
            if type == '':
                type = name
            print('{prefix}## @return ({type}) {descr}'.format(prefix=prefix, type=type, name=name, descr=''.join(descr)))

        print(declaration)
        if len(docstringlines) == 0:
            print(line, end='')
    else:
        print(line, end='')
