"""Human targets of a viral taxon, and the peptides binding them.

Takes an NCBI taxon id and exports two TSV files into `data/`, named after
the taxon's scientific name (e.g. `ebolavirus_human_targets.tsv`).

`<taxon>_human_targets.tsv` -- one row per human protein interacting with a
viral protein under the taxon, i.e. `accession1` of the valid `vh`
descriptions whose `left_value2` falls in the taxon's nested-set bounds, so
every strain below it counts:

  accession, name, descriptions, interactions, viral_proteins

`name` is the target's distinct `name1` values, `", "`-separated should
there be several. `descriptions` is the number of those descriptions
involving the target, `interactions` the number of distinct viral proteins
among them -- an interaction is a protein pair, and the target is one side
of each. `viral_proteins` lists those viral proteins as
`accession:start-stop` triples, `", "`-separated.

`<taxon>_target_peptides.tsv` -- the peptides reported as binding one of
those targets, one row per distinct (peptide, target, source protein,
position):

  sequence, target_accession, target_name,
  source_type, source_accession, source_start, source_stop,
  peptide_start, peptide_stop, source_sequence

A peptide is a mapping entry of 5-20 aa inclusive located on the target's
*partner* in a valid description: another human protein (`hh`, the target
on either side) or a viral protein of any virus (`vh`). Entries are
read through `drakkar.mappings.occurrences`, which places each one on the
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

Valid PPI filter applied (curated, both accessions live, current revision).

Run: uv run python scripts/taxon_human_target_peptides.py --ncbi-taxon-id 186536
"""

from __future__ import annotations

import argparse
import re

import polars as pl

from drakkar.db import ROOT, connect, fetch
from drakkar.mappings import occurrences
from drakkar.runs import log_run

PEPTIDE_MIN, PEPTIDE_MAX = 5, 20

VALID = """
    d.state = 'curated'
    AND d.is_obsolete1 IS FALSE
    AND d.is_obsolete2 IS FALSE
    AND d.deleted_at IS NULL
"""

TAXON_QUERY = """
SELECT t.left_value, t.right_value, n.name
FROM taxon t
JOIN taxon_name n ON n.taxon_id = t.taxon_id AND n.name_class = 'scientific name'
WHERE t.ncbi_taxon_id = %(ncbi_taxon_id)s
"""

TARGETS_QUERY = f"""
SELECT
    d.accession1 AS accession,
    string_agg(DISTINCT d.name1, ', ' ORDER BY d.name1) AS name,
    count(*) AS descriptions,
    count(DISTINCT (d.accession2, d.start2, d.stop2)) AS interactions,
    string_agg(
        DISTINCT d.accession2 || ':' || d.start2 || '-' || d.stop2, ', '
    ) AS viral_proteins
FROM dataset d
WHERE {VALID}
  AND d.type = 'vh'
  AND d.left_value2 BETWEEN %(left)s AND %(right)s
GROUP BY d.accession1
ORDER BY descriptions DESC, accession
"""

# Descriptions with a mapping involving a target, from both interactomes and
# any virus. Which side is the target and which the peptide source is sorted
# out in peptides().
MAPPINGS_QUERY = f"""
SELECT
    protein1_id, type1, accession1, start1, stop1, mapping1,
    protein2_id, type2, accession2, start2, stop2, mapping2
FROM dataset d
WHERE {VALID}
  AND (json_array_length(d.mapping1) > 0 OR json_array_length(d.mapping2) > 0)
  AND (d.accession1 = ANY(%(targets)s) OR d.accession2 = ANY(%(targets)s))
"""


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
        cur.execute(MAPPINGS_QUERY, {"targets": list(targets)})
        rows = cur.fetchall()
        ids = sorted({row[0] for row in rows} | {row[6] for row in rows})
        cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
        sequences = dict(cur.fetchall())

    records = []
    for row in rows:
        side1, side2 = row[:6], row[6:]
        for target, source in ((side1, side2), (side2, side1)):
            target_accession = target[2]
            if target_accession not in targets:
                continue
            protein_id, type_, accession, start, stop, mapping = source
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
        pl.DataFrame(records)
        .unique()
        .sort(
            "target_accession",
            "sequence",
            "source_accession",
            "source_start",
            "peptide_start",
        )
    )


def main(ncbi_taxon_id: int) -> None:
    taxon = fetch(TAXON_QUERY, {"ncbi_taxon_id": ncbi_taxon_id})
    if taxon.is_empty():
        raise SystemExit(f"no taxon with NCBI id {ncbi_taxon_id}")
    left, right, name = taxon.row(0)
    slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
    targets_tsv = ROOT / "data" / f"{slug}_human_targets.tsv"
    peptides_tsv = ROOT / "data" / f"{slug}_target_peptides.tsv"
    targets_tsv.parent.mkdir(parents=True, exist_ok=True)

    targets = fetch(TARGETS_QUERY, {"left": left, "right": right})
    if targets.is_empty():
        raise SystemExit(f"no valid vh description under {name} ({ncbi_taxon_id})")
    targets.write_csv(targets_tsv, separator="\t")

    names = dict(zip(targets["accession"], targets["name"], strict=True))
    peps = peptides(names)
    peps.write_csv(peptides_tsv, separator="\t")

    source_types = peps["source_type"]
    log_run(
        "scripts/taxon_human_target_peptides.py",
        f"{targets_tsv.relative_to(ROOT)} + {peptides_tsv.relative_to(ROOT)}",
        {
            "taxon": f"{name} (NCBI {ncbi_taxon_id})",
            "human targets": targets.height,
            "target descriptions": int(targets["descriptions"].sum()),
            "target interactions": int(targets["interactions"].sum()),
            "peptides": peps.height,
            "distinct peptide sequences": peps["sequence"].n_unique(),
            "peptides from human sources": int((source_types == "h").sum()),
            "peptides from viral sources": int((source_types == "v").sum()),
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export the human targets of a viral taxon and their peptides."
    )
    parser.add_argument(
        "--ncbi-taxon-id", type=int, required=True, help="NCBI taxon id of the virus"
    )
    main(parser.parse_args().ncbi_taxon_id)
