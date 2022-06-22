""" continuous.py

A Sphinx extension to number chapters continuously across toctrees

Derived from sphinx-multitoc-numbering with thanks
https://github.com/executablebooks/sphinx-multitoc-numbering"""

from typing import Dict, Iterable, List, Set, Tuple, cast

from docutils import nodes
from docutils.nodes import Element
from sphinx import addnodes
from sphinx.environment import BuildEnvironment
from sphinx.environment.collectors.toctree import TocTreeCollector
from sphinx.locale import __
from sphinx.util import logging, url_re

logger = logging.getLogger(__name__)


def find_numbered_toctree_nodes(
    env, iterable: Iterable[addnodes.toctree]
) -> Tuple[Set[str], List[addnodes.toctree]]:
    """Recursively walk the toctree, recording docnames and numbered nodes"""
    toctree_nodes = []
    assigned = set([])
    for node in iterable:
        if node["numbered"]:
            toctree_nodes.append(node)
        else:
            for _, ref in node["entries"]:
                assigned.add(ref)
                doctree = env.get_doctree(ref)
                inner_toctree_nodes = doctree.traverse(addnodes.toctree)
                # RECURSION
                # Base case: All nodes in inner_toctree_nodes are numbered,
                #            or inner_toctree_nodes is empty
                # This'll happen eventually as long as the toctree is acyclic,
                # as we are guaranteed to get a level deeper with every call.
                # Sphinx already requires an acyclic toctree, so we're fine.
                inner_assigned, inner_toctree_nodes = find_numbered_toctree_nodes(
                    env, inner_toctree_nodes
                )
                assigned.update(inner_assigned)
                toctree_nodes.extend(inner_toctree_nodes)
    return assigned, toctree_nodes


def get_toctree_nodes(env: BuildEnvironment) -> Tuple[Set[str], List[addnodes.toctree]]:
    """Get all numbered toctrees, in the order they appear in the document

    Walks the entire toctree, starting with the index, and records numbered
    toctrees as it finds them."""
    ## Get the toctrees in the correct order
    # Sphinx toctrees always start at index
    index_doctree = env.get_doctree("index")
    index_toctree_nodes = index_doctree.traverse(addnodes.toctree)

    assigned, toctree_nodes = find_numbered_toctree_nodes(env, index_toctree_nodes)
    assigned.add("index")

    if assigned != env.numbered_toctrees:
        logger.warning(
            "Couldn't number some toctrees: {env.numbered_toctrees - assigned}",
            type="toc",
        )
    return assigned, toctree_nodes


def assign_section_numbers(self, env: BuildEnvironment) -> List[str]:
    """Assign a section number to each heading under a numbered toctree."""
    # a list of all docnames whose section numbers changed
    rewrite_needed = []

    old_secnumbers = env.toc_secnumbers
    env.toc_secnumbers = {}
    self.last_chapter_number = 0

    assigned, toctree_nodes = get_toctree_nodes(env)

    def _walk_toc(
        node: Element, secnums: Dict, depth: int, titlenode: nodes.title = None
    ) -> None:
        # titlenode is the title of the document, it will get assigned a
        # secnumber too, so that it shows up in next/prev/parent rellinks
        for subnode in node.children:
            if isinstance(subnode, nodes.bullet_list):
                numstack.append(0)
                _walk_toc(subnode, secnums, depth - 1, titlenode)
                numstack.pop()
                titlenode = None
            elif isinstance(subnode, nodes.list_item):
                _walk_toc(subnode, secnums, depth, titlenode)
                titlenode = None
            elif isinstance(subnode, addnodes.only):
                # at this stage we don't know yet which sections are going
                # to be included; just include all of them, even if it leads
                # to gaps in the numbering
                _walk_toc(subnode, secnums, depth, titlenode)
                titlenode = None
            elif isinstance(subnode, addnodes.compact_paragraph):
                numstack[-1] += 1
                reference = cast(nodes.reference, subnode[0])

                # if a new chapter is encountered increment the chapter number
                if len(numstack) == 1:
                    self.last_chapter_number += 1
                if depth > 0:
                    number = list(numstack)
                    secnums[reference["anchorname"]] = tuple(numstack)
                else:
                    number = None
                    secnums[reference["anchorname"]] = None
                reference["secnumber"] = number
                if titlenode:
                    titlenode["secnumber"] = number
                    titlenode = None
            elif isinstance(subnode, addnodes.toctree):
                _walk_toctree(subnode, depth)

    def _walk_toctree(toctreenode: addnodes.toctree, depth: int) -> None:
        if depth == 0:
            return
        for _, ref in toctreenode["entries"]:
            if url_re.match(ref) or ref == "self":
                # don't mess with those
                continue
            elif ref in assigned:
                logger.warning(
                    __(
                        "%s is already assigned section numbers (nested numbered toctree?)"
                    ),
                    ref,
                    location=toctreenode,
                    type="toc",
                    subtype="secnum",
                )
            elif ref in env.tocs:
                secnums = {}  # type: Dict[str, Tuple[int, ...]]
                env.toc_secnumbers[ref] = secnums
                assigned.add(ref)
                _walk_toc(env.tocs[ref], secnums, depth, env.titles.get(ref))
                if secnums != old_secnumbers.get(ref):
                    rewrite_needed.append(ref)

    for toctreenode in toctree_nodes:
        depth = toctreenode.get("numbered", 0)
        if depth:
            # every numbered toctree continues the numbering
            numstack = [self.last_chapter_number]
            _walk_toctree(toctreenode, depth)

    return rewrite_needed


def setup(app):
    TocTreeCollector.assign_section_numbers = assign_section_numbers
