"""Peptides reported as binding a given list of human proteins.

Reads human UniProt accessions from standard input -- so the list can come
from any command producing them -- and exports one row per distinct
(peptide, target, source protein, position):

  sequence, target_accession, target_name,
  source_type, source_ncbi_taxon_id, source_taxon, source_name,
  source_accession, source_start, source_stop,
  peptide_start, peptide_stop, source_sequence

A peptide is a mapping entry of 5-20 aa inclusive located on the target's
*partner* in a valid description: another human protein (`hh`, the target on
either side) or a viral protein of any virus (`vh`). Entries are read
through `drakkar.mappings.occurrences`, which places each one on the
protein's sequences, so an entry that fits nowhere on its protein yields no
peptide at all.

The source is the protein the peptide lies on, identified by the
`(accession, start, stop)` triple of the description's partner side, and
`source_sequence` is its sequence -- the mature region for a viral
polyprotein. A peptide is on an isoform only when it is absent from the
canonical sequence (see `occurrences`); the isoform is then the source, as
`(isoform, 1, length)`.
`peptide_start` / `peptide_stop` are 1-based on `source_sequence`.
`source_type` is `h` or `v`.

`source_ncbi_taxon_id` and `source_taxon` are the source's taxon as curated,
strain-level for a virus and 9606 / Homo sapiens for a human protein.
`source_name` is the name the curation team gave a viral protein (`name2`,
stable per mature region) and the UniProt gene name of a human one.

Human names -- `target_name`, and `source_name` of a human source -- come
from the database, never from the input: the protein's distinct `name1`
values, `", "`-separated should there be several, or its `hh` `name2` values
for a protein curated only as side 2.

Input: the first whitespace-separated field of each line, kept when it reads
as a UniProt accession, so a header row, blank lines and trailing columns
need no stripping -- the lines dropped are counted on stderr. Accessions are
deduplicated, and one with no valid description yields no row.

Valid PPI filter applied (curated, both accessions live, current revision).

Run: cut -f3 data/HH-Ferroptosis_UniProt.tsv | uv run python
     scripts/human_target_peptides.py --output data/ferroptosis_peptides.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from collections.abc import Iterable
from pathlib import Path

import polars as pl

from drakkar.db import connect
from drakkar.mappings import occurrences
from drakkar.runs import log_run

PEPTIDE_MIN, PEPTIDE_MAX = 5, 20

ACCESSION = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)

# Spelled out so a target list with no peptide at all still writes a header.
SCHEMA: dict[str, pl.DataType] = {
    "sequence": pl.String(),
    "target_accession": pl.String(),
    "target_name": pl.String(),
    "source_type": pl.String(),
    "source_ncbi_taxon_id": pl.Int64(),
    "source_taxon": pl.String(),
    "source_name": pl.String(),
    "source_accession": pl.String(),
    "source_start": pl.Int64(),
    "source_stop": pl.Int64(),
    "peptide_start": pl.Int64(),
    "peptide_stop": pl.Int64(),
    "source_sequence": pl.String(),
}

VALID = """
    d.state = 'curated'
    AND d.is_obsolete1 IS FALSE
    AND d.is_obsolete2 IS FALSE
    AND d.deleted_at IS NULL
"""

# A human protein's name, from the side it is curated on: `name1` is the
# UniProt gene name, `name2` the curator-chosen one, used only for a human
# protein that never appears as side 1.
NAMES_QUERY = f"""
WITH named AS (
    SELECT 1 AS side, d.accession1 AS accession, d.name1 AS name
    FROM dataset d
    WHERE {VALID} AND d.accession1 = ANY(%(accessions)s)
    UNION
    SELECT 2 AS side, d.accession2 AS accession, d.name2 AS name
    FROM dataset d
    WHERE {VALID} AND d.type = 'hh' AND d.accession2 = ANY(%(accessions)s)
)
SELECT accession, string_agg(DISTINCT name, ', ' ORDER BY name) AS name
FROM named
GROUP BY accession, side
ORDER BY accession, side
"""

# Descriptions with a mapping involving a target, from both interactomes and
# any virus. Which side is the target and which the peptide source is sorted
# out in peptides().
MAPPINGS_QUERY = f"""
SELECT
    protein1_id, type1, accession1, start1, stop1,
    name1, ncbi_taxon_id1, taxon1, mapping1,
    protein2_id, type2, accession2, start2, stop2,
    name2, ncbi_taxon_id2, taxon2, mapping2
FROM dataset d
WHERE {VALID}
  AND (json_array_length(d.mapping1) > 0 OR json_array_length(d.mapping2) > 0)
  AND (
    d.accession1 = ANY(%(targets)s)
    OR (d.type = 'hh' AND d.accession2 = ANY(%(targets)s))
  )
"""


def accessions(lines: Iterable[str]) -> tuple[list[str], int]:
    """The accessions read off the input lines, and the number of lines dropped.

    One accession per line, in the first whitespace-separated field, kept in
    input order without duplicates. A line whose first field does not read as
    a UniProt accession -- a header, a comment, a blank line -- is dropped
    and counted.
    """
    kept: list[str] = []
    dropped = 0
    for line in lines:
        fields = line.split()
        if fields and ACCESSION.fullmatch(fields[0]):
            kept.append(fields[0])
        elif fields:
            dropped += 1
    return list(dict.fromkeys(kept)), dropped


def names(humans: list[str]) -> dict[str, str]:
    """The name of every human protein the database holds a valid description for."""
    with connect() as conn, conn.cursor() as cur:
        cur.execute(NAMES_QUERY.encode(), {"accessions": humans})
        rows = cur.fetchall()
    # Ordered by side, so side 1 wins where a protein is curated on both.
    found: dict[str, str] = {}
    for accession, name in rows:
        found.setdefault(accession, name)
    return found


def _sources(
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


def peptides(targets: dict[str, str]) -> pl.DataFrame:
    with connect() as conn, conn.cursor() as cur:
        cur.execute(MAPPINGS_QUERY.encode(), {"targets": list(targets)})
        rows = cur.fetchall()
        ids = sorted({row[0] for row in rows} | {row[9] for row in rows})
        cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
        sequences = dict(cur.fetchall())

    # A human source carries its UniProt gene name rather than the name the
    # curators gave it as side 2, so its name is resolved like a target's.
    # Side 1 is always human, side 2 only in `hh`.
    humans = {row[2] for row in rows} | {row[11] for row in rows if row[10] == "h"}
    human_names = names(sorted(humans))

    records: list[dict[str, str | int]] = []
    for row in rows:
        side1, side2 = row[:9], row[9:]
        for target, source in ((side1, side2), (side2, side1)):
            target_accession = target[2]
            if target_accession not in targets:
                continue
            protein_id, type_, accession, start, stop = source[:5]
            name, ncbi_taxon_id, taxon, mapping = source[5:]
            protein = sequences[protein_id]
            sources = _sources(protein, accession, start, stop)
            for occurrence in occurrences(mapping, sources, PEPTIDE_MIN, PEPTIDE_MAX):
                # An isoform source is full length: its triple is
                # (isoform, 1, length).
                source_sequence = sources[occurrence.accession]
                on_isoform = occurrence.accession != accession
                records.append(
                    {
                        "sequence": occurrence.sequence,
                        "target_accession": target_accession,
                        "target_name": targets[target_accession],
                        "source_type": type_,
                        "source_ncbi_taxon_id": ncbi_taxon_id,
                        "source_taxon": taxon,
                        "source_name": (
                            human_names[accession] if type_ == "h" else name
                        ),
                        "source_accession": occurrence.accession,
                        "source_start": 1 if on_isoform else start,
                        "source_stop": len(source_sequence) if on_isoform else stop,
                        "peptide_start": occurrence.start,
                        "peptide_stop": occurrence.stop,
                        "source_sequence": source_sequence,
                    }
                )
    # The same peptide at the same place, reported by several descriptions, is
    # one row.
    return (
        pl.DataFrame(records, schema=SCHEMA)
        .unique()
        # Grouped by target, then by source: `h` sorts before `v`, so a
        # target's human sources come first, and each taxon's sources are
        # together. The remaining keys only make the order total.
        .sort(
            "target_accession",
            "source_type",
            "source_ncbi_taxon_id",
            "source_name",
            "source_accession",
            "source_start",
            "sequence",
            "peptide_start",
        )
    )


def main(output: Path) -> None:
    targets, dropped = accessions(sys.stdin)
    if dropped:
        print(f"{dropped} input line(s) carried no accession", file=sys.stderr)
    if not targets:
        raise SystemExit("no accession on standard input")

    found = names(targets)
    peps = peptides(found)
    output.parent.mkdir(parents=True, exist_ok=True)
    peps.write_csv(output, separator="\t")

    source_types = peps["source_type"]
    log_run(
        "scripts/human_target_peptides.py",
        str(output),
        {
            "accessions read": len(targets),
            "targets in the database": len(found),
            "targets with a peptide": peps["target_accession"].n_unique(),
            "peptides": peps.height,
            "distinct peptide sequences": peps["sequence"].n_unique(),
            "peptides from human sources": int((source_types == "h").sum()),
            "peptides from viral sources": int((source_types == "v").sum()),
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export the peptides binding the human proteins read on stdin."
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="TSV file the peptides go to"
    )
    args = parser.parse_args()
    main(args.output)
