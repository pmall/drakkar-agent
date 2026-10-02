"""Human-interactor degree of each strain-abstracted viral protein.

A viral protein is `(virus, name2)`: the NCBI species of `taxon2` and the
curated `name2`, strains collapsed, so "Nucleoprotein" curated on 40
Influenza A strains is one protein. Its degree is the number of distinct
human interactors (`accession1`) across the valid `vh` descriptions of all
those strains:

  virus, name2, n_human

A `(viral protein, human protein)` pair is supported by the descriptions of
every strain collapsed into it: its evidence is the number of distinct
publications (`pmid`) among them, and it is binary when at least one of them
used a binary detection method (`drakkar.binary_methods`). `--only-binary`
keeps the binary pairs, and `--min-publications` the pairs reaching that many
publications; with both, a pair must satisfy both. With neither, every pair
counts. The degree counts the kept pairs only, and a viral protein with none
is dropped.

Rows are sorted by degree, highest first.

Valid PPI filter applied, vh only.

Run: uv run python scripts/viral_protein_degree.py --output viral_protein_degree.tsv
     [--only-binary] [--min-publications 2]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import polars as pl

from drakkar.binary_methods import BINARY_PSIMI_IDS
from drakkar.db import dataset_path, fetch
from drakkar.descriptions import valid
from drakkar.stats import show
from drakkar.taxonomy import roll_up

DESCRIPTIONS_QUERY = f"""
SELECT DISTINCT name2, taxon2, left_value2, accession1, pmid, psimi_id
FROM dataset d
WHERE {valid("d")} AND d.type = 'vh'
"""


def main(only_binary: bool, min_publications: int | None, filename: Path) -> None:
    rolled = roll_up(fetch(DESCRIPTIONS_QUERY))

    pairs = rolled.group_by("virus", "name2", "accession1").agg(
        publications=pl.col("pmid").n_unique(),
        is_binary=pl.col("psimi_id").is_in(BINARY_PSIMI_IDS).any(),
    )
    if only_binary:
        pairs = pairs.filter("is_binary")
    if min_publications is not None:
        pairs = pairs.filter(pl.col("publications") >= min_publications)

    degree = (
        pairs.group_by("virus", "name2")
        .agg(n_human=pl.len())
        .sort(["n_human", "virus", "name2"], descending=[True, False, False])
    )
    degree.write_csv(dataset_path(filename), separator="\t")

    show(
        {
            "only binary": "yes" if only_binary else "no",
            "minimum publications": (
                "none" if min_publications is None else min_publications
            ),
            "viral proteins": degree.height,
            "interactions": pairs.height,
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export the human-interactor degree of each viral protein."
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="TSV file name the degrees go to, in data/<database>/",
    )
    parser.add_argument(
        "--only-binary",
        action="store_true",
        help="keep the pairs supported by a binary detection method",
    )
    parser.add_argument(
        "--min-publications",
        type=int,
        help="keep the pairs supported by at least this many publications",
    )
    args = parser.parse_args()
    main(args.only_binary, args.min_publications, args.output)
