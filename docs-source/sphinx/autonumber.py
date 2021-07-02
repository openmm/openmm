from docutils.nodes import Text, label, reference, section
from docutils.parsers.rst import roles
from sphinx.roles import XRefRole


class autonumber(label):
    pass


class autonumber_ref(reference):
    pass


def autonumber_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    return ([autonumber(text=text)], [])


def doctree_resolved(app, doctree, docname):
    index = {}
    ref_table = {}
    if app.config.autonumber_by_chapter:
        # Record the number of each chapter
        env = app.builder.env
        section_numbers = {}
        for doc in env.toc_secnumbers:
            sections = env.toc_secnumbers[doc]
            for section_id in sections:
                section_numbers[section_id[1:]] = sections[section_id]
        last_chapter = -1

    # Assign numbers to all the autonumbered objects.

    for node in doctree.traverse(autonumber):
        category = node.astext().split(",")[0]
        if category in index:
            next_number = index[category] + 1
        else:
            next_number = 1
        if app.config.autonumber_by_chapter:
            parent = node.parent
            chapter = None
            while chapter is None:
                if isinstance(parent, section):
                    chapter = parent
                parent = parent.parent
            # If assigning a number fails, ignore app.config.autonumber_by_chapter
            try:
                chapter = section_numbers[chapter.attributes["ids"][0]][0]
            except Exception:
                index[category] = next_number
                new_node = Text("%s %d" % (category, next_number))
            else:
                if chapter != last_chapter:
                    index = {}
                    next_number = 1
                new_node = Text("%s %d-%d" % (category, chapter, next_number))
                last_chapter = chapter
        else:
            new_node = Text("%s %d" % (category, next_number))
        index[category] = next_number
        ref_table[node.astext()] = new_node
        node.parent.replace(node, new_node)

    # Replace references with the name of the referenced object

    for ref_info in doctree.traverse(autonumber_ref):
        target = ref_info["reftarget"]
        if target not in ref_table:
            raise ValueError("Unknown target for autonumber reference: " + target)
        ref_info.replace_self(Text(ref_table[target].astext()))


def setup(app):
    app.add_config_value("autonumber_by_chapter", True, False)
    roles.register_local_role("autonumber", autonumber_role)
    app.add_node(autonumber)
    app.add_node(autonumber_ref)
    app.add_role("autonumref", XRefRole(nodeclass=autonumber_ref))
    app.connect("doctree-resolved", doctree_resolved)
