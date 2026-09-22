"""Reading `mapping1` / `mapping2`, shared by dataset scripts.

A mapping is a list of entries, each a sequence curators reported as
sufficient to bind the partner, with the positions they recorded it at:

    [{"sequence": "...",
      "isoforms": [{"accession": "...",
                    "occurrences": [{"start": 1, "stop": 9, "identity": 100}]}]}]

`occurrences` is the whole module: it normalises that (coordinates and
identities are stored as JSON strings about a third of the time), flattens
it to one row per occurrence, and corrects each position against the protein
sequence it was curated on.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Literal, NamedTuple, TypedDict


class RecordedOccurrence(TypedDict):
    """One recorded position, as stored -- numbers or the strings of them."""

    start: int | str
    stop: int | str
    identity: int | float | str


class Isoform(TypedDict):
    """The sequence an entry was recorded against, and where on it."""

    accession: str
    occurrences: list[RecordedOccurrence]


class Entry(TypedDict):
    """One curated sequence and the sequences it was recorded on."""

    sequence: str
    isoforms: list[Isoform]


type How = Literal[
    "recorded",
    "relocated",
    "searched",
    "not verified",
    "recorded despite mismatch",
]


class Occurrence(NamedTuple):
    """One mapping occurrence, placed on a protein sequence.

    `start` and `stop` are 1-based on the sequence given for `accession`.

    `how` says where the position comes from:

    - `recorded` -- the recorded position reads back as the sequence.
    - `relocated` -- it does not read back, but the sequence is on the
      protein, so the match nearest it is taken. A position curated in
      another frame lands here: the polyprotein rather than the mature
      protein, or one isoform rather than another.
    - `searched` -- the entry recorded no position, or recorded one against
      an accession with no sequence given, so the sequence was searched for.
      Every match is an occurrence: a sequence does occur twice on the same
      protein, and there is nothing to disambiguate against.
    - `not verified` -- curated below identity 100, so it was aligned rather
      than matched; searching for it is meaningless and the recorded
      position is kept as is.
    - `recorded despite mismatch` -- the position does not read back and the
      sequence is nowhere on the protein, typically because the curated
      sequence includes a residue the protein does not have. The recorded
      position is kept: it still says where the sequence was meant to sit.

    An occurrence that can be placed nowhere -- a recorded position that
    does not fit the protein, or a sequence searched for and never found --
    yields no row at all.
    """

    sequence: str
    accession: str
    start: int
    stop: int
    how: How


def occurrences(
    mapping: list[Entry],
    sequences: dict[str, str],
    min_length: int,
    max_length: int,
) -> Iterator[Occurrence]:
    """Yield every occurrence of the mapping, corrected, in document order.

    `sequences` maps accession to protein sequence: the canonical accession
    and its isoforms for a protein curated full length, or the single
    accession of a mature viral protein, whose sequence is then the mature
    region `accession2[start2:stop2]` rather than the polyprotein. Recorded
    coordinates are 1-based on the sequence of the accession they were
    recorded against.

    Entries whose sequence is shorter than `min_length` or longer than
    `max_length` are skipped.
    """
    for entry in mapping:
        curated = entry["sequence"]
        if not min_length <= len(curated) <= max_length:
            continue
        if not entry["isoforms"]:
            yield from _search(curated, sequences)
            continue
        for isoform in entry["isoforms"]:
            accession = isoform["accession"]
            for occurrence in isoform["occurrences"]:
                sequence = sequences.get(accession)
                if sequence is None:
                    yield from _search(curated, sequences)
                    continue
                start = int(occurrence["start"])
                identity = float(occurrence["identity"])
                for position, how in _correct(sequence, curated, start, identity):
                    yield _row(curated, accession, position, how)


def _search(curated: str, sequences: dict[str, str]) -> Iterator[Occurrence]:
    """Yield an occurrence per exact match of `curated` in `sequences`."""
    for accession, sequence in sequences.items():
        for position in _positions(sequence, curated):
            yield _row(curated, accession, position, "searched")


def _correct(
    sequence: str, curated: str, start: int, identity: float
) -> list[tuple[int, How]]:
    """Where `curated` really sits on `sequence`, given its recorded start.

    An occurrence curated below identity 100 was aligned rather than
    matched, so searching for it is meaningless and its position is kept.
    """
    exact = identity == 100
    fits = 1 <= start and start + len(curated) - 1 <= len(sequence)
    if fits and sequence[start - 1 : start - 1 + len(curated)] == curated:
        return [(start, "recorded")]
    if exact:
        hits = _positions(sequence, curated)
        if hits:
            return [(min(hits, key=lambda hit: abs(hit - start)), "relocated")]
    if not fits:
        return []
    if not exact:
        return [(start, "not verified")]
    return [(start, "recorded despite mismatch")]


def _positions(sequence: str, curated: str) -> list[int]:
    """Every 1-based position where `curated` occurs in `sequence`."""
    return [i + 1 for i in range(len(sequence)) if sequence.startswith(curated, i)]


def _row(curated: str, accession: str, position: int, how: How) -> Occurrence:
    return Occurrence(curated, accession, position, position + len(curated) - 1, how)
