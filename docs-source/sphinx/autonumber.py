from pathlib import Path

from docutils.nodes import Text, label, reference, section
from docutils.parsers.rst import roles
from sphinx.roles import XRefRole


class autonumber(label):
    pass


class autonumber_ref(reference):
    pass


def autonumber_role(name, rawtext, text, lineno, inliner, options={}, content=[]):
    return ([autonumber(text=text)], [])


def warn(*args, **kwargs):
    """Issue a warning/notification to the user. Same signature as `print`"""
    print("\nAutonumber: ", *args, **kwargs)


def chapter_numbers_by_section(env):
    """Collect a mapping from a section's reference key to its position in the TOC

    The reference key is of the form `f'{path-to-doc}:{section-id}'`.
    The position in the TOC is useful eg. to label sections with numbers.
    It is returned as a tuple of ints"""
    # Record the number of each chapter
    section_numbers = {}
    for doc in env.toc_secnumbers:
        sections = env.toc_secnumbers[doc]
        for section_id in sections:
            # Include the src to disambiguate duplicates
            src = str(Path(env.srcdir)/doc)
            key = f"{src}:{section_id[1:]}"
            if key in section_numbers:
                warn(f"{section_id} is duplicated in {src}")
            section_numbers[key] = sections[section_id]
    return section_numbers


def get_chapter(node, depth, section_numbers):
    """Get the numerical position of the chapter in which node resides

    args:
        node:
            A docutils node whose chapter we want the number of
        depth:
            How many levels deep into the toctree is a "chapter"
        section_numbers:
            The output of chapter_numbers_by_section(env)
    """
    parent = node.parent
    chapter = None
    while chapter is None:
        if isinstance(parent, section):
            chapter = parent
        parent = parent.parent
    src = str(Path(chapter.source).with_suffix(""))
    chapter_id = chapter.attributes["ids"][0]
    key = src + ":" + chapter_id
    try:
        chapter = section_numbers[key][:depth]
    except KeyError:
        # The above will fail if the section is at the top of a file;
        # There doesn't seem to be a way to get the top section label
        # in chapter_numbers_by_section, so we'll just assume that if
        # the above fails, we're looking for a section with no label:
        key = src + ":"
        warn(f"Assuming {repr(chapter_id)} is a top level section")
        chapter = section_numbers[key][:depth]
    return ".".join(str(i) for i in chapter)


def proc_autonumber(app, doctree, docname):
    index = {}
    ref_table = {}
    autonum_depth = app.config.autonumber_by_depth
    if autonum_depth:
        section_numbers = chapter_numbers_by_section(app.builder.env)
        last_chapter = None

    # Assign numbers to all the autonumbered objects.

    for node in doctree.traverse(autonumber):
        category = node.astext().split(",")[0]
        if category in index:
            next_number = index[category] + 1
        else:
            next_number = 1
        if autonum_depth:
            chapter = get_chapter(node, autonum_depth, section_numbers)
            if chapter != last_chapter:
                index = {}
                next_number = 1
            new_node = Text(f"{category} {chapter}-{next_number}")
            last_chapter = chapter
        else:
            new_node = Text(f"{category} {next_number}")
        index[category] = next_number
        ref_table[node.astext()] = new_node
        node.parent.replace(node, new_node)

    # Replace references with the name of the referenced object

    for ref_info in doctree.traverse(autonumber_ref):
        target = ref_info["reftarget"]
        if target not in ref_table:
            raise ValueError(f"Unknown target for autonumber reference: {target}")
        ref_info.replace_self(Text(ref_table[target].astext()))


def setup(app):
    app.add_config_value("autonumber_by_depth", 1, "env")
    roles.register_local_role("autonumber", autonumber_role)
    app.add_node(autonumber)
    app.add_node(autonumber_ref)
    app.add_role("autonumref", XRefRole(nodeclass=autonumber_ref))
    app.connect("doctree-resolved", proc_autonumber)
