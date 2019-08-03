from docutils.parsers.rst import roles
from docutils.nodes import Text, reference, section, label
from sphinx.roles import XRefRole

class autonumber(label):
    pass

class autonumber_ref(reference):
    pass

def autonumber_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    return ([autonumber(text=text)], [])

def doctree_resolved(app, doctree, docname):
    index = {};
    refTable = {}
    if app.config.autonumber_by_chapter:
        # Record the number of each chapter

        env = app.builder.env
        sectionNumbers = {}
        for doc in env.toc_secnumbers:
            sections = env.toc_secnumbers[doc]
            for sectionId in sections:
                sectionNumbers[sectionId[1:]] = sections[sectionId]
        lastChapter = -1

    # Assign numbers to all the autonumbered objects.

    for node in doctree.traverse(autonumber):
        category = node.astext().split(',')[0]
        if category in index:
            nextNumber = index[category]+1
        else:
            nextNumber = 1
        if app.config.autonumber_by_chapter:
            parent = node.parent
            chapter = None
            while chapter is None:
                if isinstance(parent, section):
                    chapter = parent
                parent = parent.parent
            chapter = sectionNumbers[chapter.attributes['ids'][0]][0]
            if chapter != lastChapter:
                index = {}
                nextNumber = 1
            newNode = Text('%s %d-%d' % (category, chapter, nextNumber))
            lastChapter = chapter
        else:
            newNode = Text('%s %d' % (category, nextNumber))
        index[category] = nextNumber
        refTable[node.astext()] = newNode
        node.parent.replace(node, newNode)

    # Replace references with the name of the referenced object

    for ref_info in doctree.traverse(autonumber_ref):
        target = ref_info['reftarget']
        if target not in refTable:
            raise ValueError('Unknown target for autonumber reference: '+target)
        ref_info.replace_self(Text(refTable[target].astext()))

def setup(app):
    app.add_config_value('autonumber_by_chapter', True, False)
    roles.register_local_role('autonumber', autonumber_role)
    app.add_node(autonumber)
    app.add_node(autonumber_ref)
    app.add_role('autonumref', XRefRole(nodeclass=autonumber_ref))
    app.connect('doctree-resolved', doctree_resolved)

