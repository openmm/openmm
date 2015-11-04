import re

def process_docstring(app, what, name, obj, options, lines):
    """This hook edits the docstrings to replace "<tt><pre>" html tags
    with sphinx directives.
    """
    def repl(m):
        s = m.group(1)
        if not s.startswith(linesep):
            s = linesep + s
        newline = '.. code-block:: c++' + linesep
        return newline + '    ' + s.replace(linesep, linesep + '    ')

    if name == 'simtk.openmm.openmm.CustomTorsionForce':
         print(lines)

    linesep = '|LINEBREAK|'
    joined = linesep.join(lines)

    joined = re.sub(r'<tt><pre>((|LINEBREAK|)?.*?)</pre></tt>', repl, joined)
    joined = re.sub(r'<tt>(.*?)</tt>', repl, joined)
    lines[:] = [(l if not l.isspace() else '') for l in joined.split(linesep)]


def setup(app):
    app.connect('autodoc-process-docstring', process_docstring)


def test():
    lines    = ['Hello World', '<tt><pre>', 'contents', '</pre></tt>', '', '<tt>contents2</tt>']
    linesRef = ['Hello World', '.. code-block:: c++', '', '    contents', '', '', '.. code-block:: c++', '', '    contents2']
    process_docstring(None, None, None, None, None, lines)
    assert lines == linesRef

if __name__ == '__main__':
    test()
