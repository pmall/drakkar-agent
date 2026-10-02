"""All valid PPI descriptions with their viral proteins and peptides.

Exported as four TSV files into a folder of `data/<database>/`, a
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
inclusive, the team convention. Entries are
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
Curation sometimes enters the same *(publication x detection method x protein
pair)* under several `stable_id`s, an `hh` pair in either order; only the
earliest one (`created_at`, then `stable_id`) is kept, in the descriptions and
in the peptides alike.

Run: uv run python scripts/graph_descriptions.py --output graph
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Literal

import polars as pl

from drakkar.db import connect, dataset_path, fetch, stream
from drakkar.descriptions import valid
from drakkar.peptides import load_sequences, peptide_occurrences, sources
from drakkar.stats import show

type Interactome = Literal["hh", "vh"]

# The two protein triples of a description, ordered: an `hh` pair may be entered
# either way round, so (A, B) and (B, A) are one pair.
PAIR = """
        LEAST(ROW(accession1, start1, stop1), ROW(accession2, start2, stop2)),
        GREATEST(ROW(accession1, start1, stop1), ROW(accession2, start2, stop2))"""

# The valid descriptions, one per (pmid, method, protein pair): the first
# curated, `stable_id` breaking a tie on `created_at`.
FIRST_CURATED = f"""(
    SELECT DISTINCT ON (type, pmid, psimi_id,{PAIR}) *
    FROM dataset d
    WHERE {valid("d")}
    ORDER BY type, pmid, psimi_id,{PAIR},
        created_at, stable_id
)"""

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
FROM {FIRST_CURATED} d
WHERE d.type = '{{interactome}}'
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
WHERE {valid("d")}
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
FROM {FIRST_CURATED} d
WHERE (json_array_length(mapping1) > 0 OR json_array_length(mapping2) > 0)
ORDER BY stable_id
"""


def descriptions_query(interactome: Interactome) -> str:
    # COPY takes no bind parameters, so the interactome is formatted in; it is
    # an Interactome literal, never outside input.
    return DESCRIPTIONS_QUERY.format(interactome=interactome)


def peptides() -> pl.DataFrame:
    with connect() as conn, conn.cursor() as cur:
        cur.execute(MAPPINGS_QUERY)
        rows = cur.fetchall()
        sequences = load_sequences(
            cur, {row[1] for row in rows} | {row[7] for row in rows}
        )

    records = []
    for stable_id, *columns in rows:
        side1, side2 = columns[:6], columns[6:]
        for protein_id, type_, accession, start, stop, mapping in (side1, side2):
            on = sources(sequences[protein_id], accession, start, stop)
            for occurrence in peptide_occurrences(mapping, on):
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


def main(folder: Path) -> None:
    out = dataset_path(folder)
    out.mkdir(exist_ok=True)
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

    show(
        {"viral proteins": viral.height, "peptides": peps.height},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export the drakkar PPI graph.")
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="folder the four TSV files go to, in data/<database>/",
    )
    main(parser.parse_args().output)
