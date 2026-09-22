"""Unique viral proteins per virus, abstracting away the strain.

``taxon2`` is strain-level (individual SARS-CoV-2 / Influenza A / HBV isolates).
We roll each viral protein up the NCBI nested-set tree to a fixed rank and count
distinct ``(taxon at that rank, curated name2)`` pairs -- so "Nucleoprotein"
curated on 40 Influenza A strains counts once.

The rollup is done in Python: one ~3k-row scan of the MV for the distinct viral
proteins, one pull of the taxon table, then an interval match. Avoids a
range-join against the unindexed 420k-row MV (that query never returns).

Reports the count at species / genus / family rank so the right granularity is
visible. Strain-level rows with no ancestor at a given rank keep their own
``taxon2`` string as the bucket.

Valid PPI filter applied, vh only.

Note: species rank is NCBI species, so SARS-CoV-1 and SARS-CoV-2 both collapse
into "Severe acute respiratory syndrome-related coronavirus".

Exports the species-rank table to data/viral_proteins_by_virus.tsv.
"""

import polars as pl

from drakkar.db import fetch

VIRAL = """
select distinct name2, taxon2, left_value2
from dataset
where state = 'curated'
  and is_obsolete1 is false
  and is_obsolete2 is false
  and deleted_at is null
  and type2 = 'v'
"""

TAXON = """
select t.node_rank, t.left_value, t.right_value, n.name
from taxon t
join taxon_name n on n.taxon_id = t.taxon_id and n.name_class = 'scientific name'
where t.node_rank in ('species', 'genus', 'family')
"""


def rollup(viral: pl.DataFrame, nodes: pl.DataFrame, rank: str) -> pl.DataFrame:
    r = nodes.filter(pl.col("node_rank") == rank).select(
        "left_value", "right_value", virus=pl.col("name")
    )
    matched = viral.join_where(
        r,
        pl.col("left_value2") >= pl.col("left_value"),
        pl.col("left_value2") <= pl.col("right_value"),
    ).select("name2", "taxon2", "virus")
    # rows with no ancestor at this rank: fall back to their strain string
    missing = viral.join(
        matched.select("name2", "taxon2").unique(),
        on=["name2", "taxon2"],
        how="anti",
    ).with_columns(virus=pl.col("taxon2"))
    return pl.concat([matched, missing.select("name2", "taxon2", "virus")])


def main() -> None:
    viral = fetch(VIRAL)
    nodes = fetch(TAXON)

    for rank in ("species", "genus", "family"):
        rolled = rollup(viral, nodes, rank)
        n_taxa = rolled["virus"].n_unique()
        n_pairs = rolled.select("virus", "name2").unique().height
        print(f"{rank:8s}: {n_taxa:4d} taxa, {n_pairs:5d} unique (taxon, name) entries")

    species = (
        rollup(viral, nodes, "species")
        .select("virus", "name2")
        .unique()
        .group_by("virus")
        .len()
        .rename({"len": "n_proteins"})
        .sort(["n_proteins", "virus"], descending=[True, False])
    )
    species.write_csv("data/viral_proteins_by_virus.tsv", separator="\t")
    print()
    print(species.head(25))


if __name__ == "__main__":
    main()
