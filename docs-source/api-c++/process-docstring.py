import os
import re
import sys

def replace_code(match):
    lines = match.group(1).split('\n')
    for i, line in enumerate(list(lines)):
        if line.startswith('*'):
            lines[i] = line[1:]
    return f"\n* .. code-block:: c++\n*    {'\n*    '.join(lines)}"


def process_file(file):
    """Edit docstrings in an XML file to remove the Python versions of code examples,
    and put the correct Sphinx directives around the C++ versions."""
    with open(file) as input:
        content = input.read()
    processed = re.sub(r'&lt;c\+\+&gt;((.|\n)*?)&lt;/c\+\+&gt;', replace_code, content)
    processed = re.sub(r'&lt;python&gt;((.|\n)*?)&lt;/python&gt;', '', processed)
    print(file, processed != content)
    if processed != content:
        with open(file, 'w') as output:
            output.write(processed)


dir = sys.argv[1]
for file in os.listdir(dir):
    process_file(os.path.join(dir, file))
