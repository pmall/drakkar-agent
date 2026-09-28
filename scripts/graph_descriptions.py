"""All valid PPI descriptions with their viral proteins and peptides.

Exported as four TSV files into the folder given on the command line, a
graph of the curated interactome: descriptions are the edges, viral proteins
and peptides hang off them by protein triple / `stable_id`.

`descriptions_hh.tsv` / `descriptions_vh.tsv` -- one row per valid
description of each interactome, i.e. the MV grain *(publication x detection
method x protein 1 region x protein 2 region)*, with the same columns:

  stable_id, pmid, psimi_id,
  accession1, start1, stop1, name1, ncbi_taxon_id1,
  accession2, start2, stop2, name2, ncbi_taxon_id2

Each protein is identified by its `(accession, start, stop)` triple -- never
by accession alone, since one viral polyprotein accession yields many mature
proteins at different coordinates. Human proteins are always full length, so
their triple is always `(accession, 1, length)`. `ncbi_taxon_idN` is the
protein's NCBI taxon: 9606 for a human protein, the strain-level taxon for a
viral one.

`viral_proteins.tsv` -- one row per distinct viral protein of
`descriptions_vh.tsv`, mature or full length:

  accession, start, stop, name, ncbi_taxon_id, sequence

keyed on the `(accession2, start2, stop2)` triple, with the curator-chosen
`name2`, and `sequence` the `start`-`stop` region of the accession's
canonical sequence -- the mature protein itself for a polyprotein region.

`peptides.tsv` -- the short binding regions curators reported, both sides
together, one row per distinct peptide of a description's protein:

  stable_id, sequence, source_type, source_accession, source_start, source_stop

A peptide is a `mapping1` / `mapping2` entry whose sequence is 5-20 aa
inclusive (a one-time exception to the 4-20 team convention). Entries are
read through `drakkar.mappings.occurrences`, which places each one on the
sequence it was curated on, so an entry that fits nowhere on its protein
yields no peptide at all. The position it was placed at is *not* exported:
the peptide is attached to the protein it belongs to, by that protein's
identifier triple, which is the description's `(accessionN, startN, stopN)`
and joins back to the descriptions files. `source_type` is that protein's
`typeN` -- `h` for a human protein, `v` for a viral one, so side 1 is always
`h`. Occurrences that differ only by position therefore collapse to one
row. The file covers both interactomes.

Valid PPI filter applied (curated, both accessions live, current revision).

Run: uv run python scripts/graph_descriptions.py data/graph
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Literal

import polars as pl

from drakkar.db import connect, fetch, stream
from drakkar.mappings import occurrences
from drakkar.runs import log_run

type Interactome = Literal["hh", "vh"]

PEPTIDE_MIN, PEPTIDE_MAX = 5, 20

VALID = """
    d.state = 'curated'
    AND d.is_obsolete1 IS FALSE
    AND d.is_obsolete2 IS FALSE
    AND d.deleted_at IS NULL
"""

DESCRIPTIONS_QUERY = f"""
SELECT
    stable_id,
    pmid,
    psimi_id,
    accession1,
    start1,
    stop1,
    name1,
    ncbi_taxon_id1,
    accession2,
    start2,
    stop2,
    name2,
    ncbi_taxon_id2
FROM dataset d
WHERE {VALID}
  AND d.type = '{{interactome}}'
ORDER BY pmid, accession1, accession2, start2, stop2, psimi_id, stable_id
"""

# name2 and the taxon are stable per triple, so DISTINCT yields one row per
# viral protein; main() checks it does.
VIRAL_PROTEINS_QUERY = f"""
SELECT DISTINCT
    d.accession2 AS accession,
    d.start2 AS start,
    d.stop2 AS stop,
    d.name2 AS name,
    d.ncbi_taxon_id2 AS ncbi_taxon_id,
    substring(p.sequences->>d.accession2 FROM d.start2 FOR d.stop2 - d.start2 + 1)
        AS sequence
FROM dataset d
JOIN proteins p ON p.id = d.protein2_id
WHERE {VALID}
  AND d.type = 'vh'
ORDER BY accession, start, stop
"""

# Descriptions reporting a binding region on either side -- a few percent of
# the view, so the peptide pass never holds the whole thing.
MAPPINGS_QUERY = f"""
SELECT
    stable_id,
    protein1_id, type1, accession1, start1, stop1, mapping1,
    protein2_id, type2, accession2, start2, stop2, mapping2
FROM dataset d
WHERE {VALID}
  AND (json_array_length(mapping1) > 0 OR json_array_length(mapping2) > 0)
ORDER BY stable_id
"""


def descriptions_query(interactome: Interactome) -> str:
    # COPY takes no bind parameters, so the interactome is formatted in; it is
    # an Interactome literal, never outside input.
    return DESCRIPTIONS_QUERY.format(interactome=interactome)


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


def peptides() -> pl.DataFrame:
    with connect() as conn, conn.cursor() as cur:
        cur.execute(MAPPINGS_QUERY)
        rows = cur.fetchall()
        ids = sorted({row[1] for row in rows} | {row[7] for row in rows})
        cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
        sequences = dict(cur.fetchall())

    records = []
    for stable_id, *columns in rows:
        side1, side2 = columns[:6], columns[6:]
        for protein_id, type_, accession, start, stop, mapping in (side1, side2):
            sources = _sources(sequences[protein_id], accession, start, stop)
            for occurrence in occurrences(mapping, sources, PEPTIDE_MIN, PEPTIDE_MAX):
                records.append(
                    {
                        "stable_id": stable_id,
                        "sequence": occurrence.sequence,
                        "source_type": type_,
                        "source_accession": accession,
                        "source_start": start,
                        "source_stop": stop,
                    }
                )
    # A peptide occurring twice on the same protein is one peptide of it once
    # its position is dropped.
    return pl.DataFrame(records).unique(maintain_order=True)


def viral_proteins() -> pl.DataFrame:
    proteins = fetch(VIRAL_PROTEINS_QUERY)
    triples = proteins.select("accession", "start", "stop").n_unique()
    if triples != proteins.height:
        raise RuntimeError(
            f"{proteins.height - triples} viral protein triple(s) with more than "
            "one name, taxon or sequence"
        )
    return proteins


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    hh_descriptions_tsv = out / "descriptions_hh.tsv"
    vh_descriptions_tsv = out / "descriptions_vh.tsv"
    viral_proteins_tsv = out / "viral_proteins.tsv"
    peptides_tsv = out / "peptides.tsv"

    # 400k rows: streamed to disk rather than materialised as a DataFrame.
    stream(descriptions_query("hh"), hh_descriptions_tsv)
    stream(descriptions_query("vh"), vh_descriptions_tsv)

    viral = viral_proteins()
    viral.write_csv(viral_proteins_tsv, separator="\t")

    peps = peptides()
    peps.write_csv(peptides_tsv, separator="\t")

    log_run(
        "scripts/graph_descriptions.py",
        f"{hh_descriptions_tsv}, {vh_descriptions_tsv}, "
        f"{viral_proteins_tsv}, {peptides_tsv}",
        {"viral proteins": viral.height, "peptides": peps.height},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export the drakkar PPI graph.")
    parser.add_argument("out", type=Path, help="folder the four TSV files go to")
    main(parser.parse_args().out)
