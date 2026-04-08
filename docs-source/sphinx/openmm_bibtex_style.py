"""Pybtex formatting style extending ``unsrt`` without DOI hyperlinks.

Bibliography entries still show ``doi:10.xxxx/...`` for readers, but omit
``https://doi.org/...`` anchors. Automated link checkers often receive HTTP 403
from doi.org when using a non-browser user agent, which breaks CI.
"""

from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.template import field, join, optional


class Style(UnsrtStyle):
    """Same as :class:`pybtex.style.formatting.unsrt.Style` but plain-text DOIs."""

    def format_doi(self, e):
        return optional [join ["doi:", field("doi", raw=True)]]
