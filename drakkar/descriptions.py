"""What makes a description valid, shared by dataset scripts."""

from __future__ import annotations

from typing import LiteralString


def valid(alias: LiteralString | None = None) -> LiteralString:
    """The WHERE conditions selecting valid PPI descriptions of the `dataset` view.

    A description is valid when it is curated, neither accession is obsolete
    and it is the current revision. `alias` is the view's alias in the query,
    for scripts that join it. The interactome is not part of it: restrict it
    with `type = 'hh'` or `type = 'vh'`.
    """
    prefix: LiteralString = f"{alias}." if alias else ""
    return f"""
    {prefix}state = 'curated'
    AND {prefix}is_obsolete1 IS FALSE
    AND {prefix}is_obsolete2 IS FALSE
    AND {prefix}deleted_at IS NULL
"""
