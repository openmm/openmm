"""
Changes section references to be the section number
instead of the title of the section.
"""

from docutils import nodes
import sphinx.domains.std

class CustomStandardDomain(sphinx.domains.std.StandardDomain):

    def __init__(self, env):
        env.settings['footnote_references'] = 'superscript'
        sphinx.domains.std.StandardDomain.__init__(self, env)

    def resolve_xref(self, env, fromdocname, builder,
                     typ, target, node, contnode):
        res = super(CustomStandardDomain, self).resolve_xref(env, fromdocname, builder,
                                                            typ, target, node, contnode)
        
        if res is None:
            return res
        
        if typ == 'ref' and not node['refexplicit']:
            docname, labelid, sectname = self.data['labels'].get(target, ('','',''))
            res['refdocname'] = docname
        
        return res

def doctree_resolved(app, doctree, docname):
    secnums = app.builder.env.toc_secnumbers
    for node in doctree.traverse(nodes.reference):
        if 'refdocname' in node:
            refdocname = node['refdocname']
            if refdocname in secnums:
                secnum = secnums[refdocname]
                emphnode = node.children[0]
                textnode = emphnode.children[0]
                
                toclist = app.builder.env.tocs[refdocname]
                anchorname = None
                for refnode in toclist.traverse(nodes.reference):
                    if refnode.astext() == textnode.astext():
                        anchorname = refnode['anchorname']
                if anchorname is None:
                    continue
                linktext = '.'.join(map(str, secnum[anchorname]))
                node.replace(emphnode, nodes.Text(linktext))

def setup(app):
    app.override_domain(CustomStandardDomain)
    app.connect('doctree-resolved', doctree_resolved)

