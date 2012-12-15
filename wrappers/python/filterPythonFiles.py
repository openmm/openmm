import sys

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
            print "%s## @private" % prefix
        if split[1][0] == '_' and split[1][1] != '_':
            # Names starting with a single _ are assumed to be private.
            print "%s## @private" % prefix
        
        # We're at the start of a class or function definition.  Find all lines that contain the declaration.
        
        declaration = line
        while len(line) > 0 and line.find(':') == -1:
            line = input.readline()
            declaration += line
        
        # Now look for a docstring.
        
        docstrings = []
        line = input.readline()
        stripped = line.lstrip()
        if stripped.startswith('"""'):
            line = stripped[3:]
            readingParameters = False
            while line.find('"""') == -1:
                docstrings.append(line)
                line = input.readline()
                if line.strip() == 'Parameters:':
                    readingParameters = True
                    line = input.readline()
                stripped = line.lstrip()
                if readingParameters and stripped.startswith('- '):
                    line = "@param %s" % stripped[2:]
                elif stripped.startswith('Returns:'):
                    line = "@return %s" % stripped[8:]
            line = line[:line.find('"""')]
            docstrings.append(line)            
        
        # Print out the docstring in Doxygen syntax, followed by the declaration.
        
        for s in docstrings:
            print "%s##%s" % (prefix, s.strip())
        print declaration
        if len(docstrings) == 0:
            print line
    else:
        print line