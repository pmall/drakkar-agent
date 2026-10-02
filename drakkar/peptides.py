"""Peptides of a mapping, shared by dataset scripts.

A peptide is a mapping entry of 5-20 aa inclusive, a short binding region as
opposed to a domain-scale mapping. Entries are read through
`drakkar.mappings.occurrences`; this module settles the length bounds and the
sequences an entry is matched against.
"""

from __future__ import annotations

from collections.abc import Iterator

import psycopg

from drakkar.mappings import Entry, Occurrence, occurrences

PEPTIDE_MIN, PEPTIDE_MAX = 5, 20


def load_sequences(
    cur: psycopg.Cursor, protein_ids: set[int]
) -> dict[int, dict[str, str]]:
    """The `sequences` of each protein row: canonical accession and isoforms."""
    cur.execute(
        "SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (sorted(protein_ids),)
    )
    return dict(cur.fetchall())


def sources(
    protein: dict[str, str], accession: str, start: int, stop: int
) -> dict[str, str]:
    """The sequences a mapping on this protein can have been curated on.

    A protein curated full length -- every human protein, and a viral protein
    that is not a mature one -- brings its canonical sequence and its
    isoforms, since a peptide may have been located on an isoform. A mature
    viral protein is a region of a polyprotein and has no isoforms of its
    own, so its mature region is the only sequence to match against.
    """
    canonical = protein[accession]
    if start == 1 and stop == len(canonical):
        return protein
    return {accession: canonical[start - 1 : stop]}


def peptide_occurrences(
    mapping: list[Entry], on: dict[str, str]
) -> Iterator[Occurrence]:
    """Every occurrence of a peptide of the mapping, on the sequences `on`."""
    return occurrences(mapping, on, PEPTIDE_MIN, PEPTIDE_MAX)
