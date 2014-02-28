from docutils.parsers.rst import Directive
from docutils.nodes import compound, raw

class SamepageDirective(Directive):
    
    has_content = True

    def run(self):
        prefix = raw('', '\\par\\begin{samepage}', format='latex')
        suffix = raw('', '\\end{samepage}\\par', format='latex')
        text = '\n'.join(self.content)
        content_node = compound(rawsource=text)
        self.state.nested_parse(self.content, self.content_offset, content_node)
        return [prefix, content_node, suffix]

def setup(app):
    app.add_directive('samepage', SamepageDirective)
