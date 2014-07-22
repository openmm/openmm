from docutils.parsers.rst import Directive
from docutils.nodes import compound, raw

class CaptionDirective(Directive):
    
    has_content = True

    def run(self):
        latexPrefix = raw('', '{\\centering', format='latex')
        latexSuffix = raw('', '\\par}\\bigskip', format='latex')
        text = '\n'.join(self.content)
        content_node = compound(rawsource=text)
        self.state.nested_parse(self.content, self.content_offset, content_node)
        content_node.attributes['classes'].append('caption')
        return [latexPrefix, content_node, latexSuffix]

def setup(app):
    app.add_directive('caption', CaptionDirective)
