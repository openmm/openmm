import re

def count_leading_whitespace(s):
    count = 0
    for c in s:
        if c.isspace():
            count += 1
        else:
            break
    return count

def process_docstring(app, what, name, obj, options, lines):
    """This hook edits the docstrings to replace "<tt><pre>" html tags,
    Breathe-style verbatim embed blocks, and deprecated markers with
    sphinx directives.
    """
    def repl(m):
        s = m.group(1)
        if not s.startswith(linesep):
            s = linesep + s
        newline = '.. code-block:: c++' + linesep
        return newline + '    ' + s.replace(linesep, linesep + '    ')
    def repl2(m):
        s = m.group(1)
        if not s.startswith(linesep):
            s = linesep + s
        newline = '|LINEBREAK|.. admonition::|LINEBREAK|   Deprecated' + linesep
        return newline + '    ' + s.replace(linesep, linesep + '    ')
    def repl3(m):
        s = m.group(1)
        return '*' + s + '*'
    def repl4(m):
        s = m.group(1)
        if s.startswith("embed:rst"):
            lines = s.split(linesep)[1:]

            # Match behaviour of Breathe
            if s.startswith("embed:rst:leading-asterisk"):
                lines = [l.replace("*", " ", 1) for l in lines]
            elif s.startswith("embed:rst:leading-slashes"):
                lines = [l.replace("///", "   ", 1) for l in lines]

            # Strip leading whitespace to match first line
            to_strip = count_leading_whitespace(lines[0])
            lines = [l[to_strip:] for l in lines]

            return linesep.join(lines)
        else:
            s = m.group(1)
            if not s.startswith(linesep):
                s = linesep + s
            newline = '.. verbatim::' + linesep
            return newline + '    ' + s.replace(linesep, linesep + '    ')
    def replace_code(m):
        s = m.group(1)
        if not s.startswith(linesep):
            s = linesep + s
        newline = linesep + '.. code-block:: python' + linesep
        return newline + '    ' + s.replace(linesep, linesep + '    ')
    def replace_subscript(m):
        """ Replace subscript tags. """
        return r'\ :sub:`{}`\ '.format(m.group(1))
    def replace_superscript(m):
        """ Replace superscript tags. """
        return r'\ :sup:`{}`\ '.format(m.group(1))

    linesep = '|LINEBREAK|'
    joined = linesep.join(lines)

    joined = re.sub(r'<tt><pre>((|LINEBREAK|)?.*?)</pre></tt>', repl, joined)
    joined = re.sub(r'<tt>(.*?)</tt>', repl, joined)
    joined = re.sub(r'@deprecated(.*?\|LINEBREAK\|)', repl2, joined, flags=re.IGNORECASE)
    joined = re.sub(r'<i>(.*?)</i>', repl3, joined)
    joined = re.sub(r'<verbatim>(.*?)</verbatim>', repl4, joined)
    joined = re.sub(r'<sub>(.*?)</sub>', replace_subscript, joined)
    joined = re.sub(r'<sup>(.*?)</sup>', replace_superscript, joined)
    joined = re.sub(r'<c\+\+>(.*?)</c\+\+>', '', joined)
    joined = re.sub(r'<python>(.*?)</python>', replace_code, joined)

    lines[:] = [(l if not l.isspace() else '') for l in joined.split(linesep)]


substitutions = {'double':'float', 'long long':'int', 'string':'str',
                 'pairii':'tuple[int, int]',
                 'vectord':'tuple[float, ...]',
                 'vectorvectorvectord':'tuple[tuple[tuple[float, ...], ...], ...]',
                 'vectori':'tuple[int, ...]',
                 'vectorvectori':'tuple[tuple[int, ...], ...]',
                 'vectorpairii':'tuple[tuple[int, int], ...]',
                 'vectorstring':'tuple[str, ...]',
                 'mapstringstring':'Mapping[str, str]',
                 'mapstringdouble':'Mapping[str, float]',
                 'mapii':'Mapping[int, int]',
                 'seti':'set[int]'
                }

def convert_type(type):
    if type in substitutions:
        type = substitutions[type]
    if '<' in type:
        match = re.match(r'vector<(.*?),.*>', type)
        if match is not None:
            type = f'tuple[{convert_type(match[1])}]'
        match = re.match(r'map<(.*?),(.*?),.*>', type)
        if match is not None:
            type = f'Mapping[{convert_type(match[1])}, {convert_type(match[2])}]'
    return type

def process_signature(app, what, name, obj, options, signature, return_annotation):
    if return_annotation is not None:
        # Convert C++ types to Python types
        if return_annotation.startswith('std::'):
            return_annotation = return_annotation[5:]
        if return_annotation.endswith(' const &'):
            return_annotation = return_annotation[:-8]
        return_annotation = convert_type(return_annotation)
    return (signature, return_annotation)


def setup(app):
    app.connect('autodoc-process-docstring', process_docstring)
    app.connect('autodoc-process-signature', process_signature)


def test():
    lines    = ['Hello World', '<tt><pre>', 'contents', '</pre></tt>', '', '<tt>contents2</tt>', 'r<sub>1</sub>', 'r<sup>1</sup>']
    linesRef = ['Hello World', '.. code-block:: c++', '', '    contents', '', '', '.. code-block:: c++', '', '    contents2', 'r\\ :sub:`1`\\ ', 'r\\ :sup:`1`\\ ']
    process_docstring(None, None, None, None, None, lines)
    assert lines == linesRef, (lines, linesRef)

if __name__ == '__main__':
    test()
