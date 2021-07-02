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
    autonum_depth = app.config.autonumber_by_depth
    if autonum_depth:
        # Record the number of each chapter
        env = app.builder.env
        section_numbers = {}
        for doc in env.toc_secnumbers:
            sections = env.toc_secnumbers[doc]
            for section_id in sections:
                section_numbers[section_id[1:]] = sections[section_id]
        last_chapter = None

    # Assign numbers to all the autonumbered objects.

    for node in doctree.traverse(autonumber):
        category = node.astext().split(",")[0]
        if category in index:
            next_number = index[category] + 1
        else:
            next_number = 1
        if autonum_depth:
            parent = node.parent
            chapter = None
            while chapter is None:
                if isinstance(parent, section):
                    chapter = parent
                parent = parent.parent
            # If assigning a number fails, ignore autonumber_by_depth config option
            try:
                chapter = section_numbers[chapter.attributes["ids"][0]][:autonum_depth]
            except Exception:
                index[category] = next_number
                new_node = Text("%s %d" % (category, next_number))
            else:
                if chapter != last_chapter:
                    index = {}
                    next_number = 1
                chapter_fmt = ".".join(str(i) for i in chapter)
                new_node = Text("%s %s-%d" % (category, chapter_fmt, next_number))
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
    app.add_config_value("autonumber_by_depth", 1, "env")
    roles.register_local_role("autonumber", autonumber_role)
    app.add_node(autonumber)
    app.add_node(autonumber_ref)
    app.add_role("autonumref", XRefRole(nodeclass=autonumber_ref))
    app.connect("doctree-resolved", doctree_resolved)
