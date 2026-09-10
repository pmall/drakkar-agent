"""Human-interactor degree of each strain-abstracted viral protein.

Builds on scripts/viral_proteins_by_virus.py: a viral protein is
``(NCBI species of taxon2, curated name2)`` -- strains collapsed. For each such
viral protein we count distinct human interactors (``accession1``), then look at
how those degrees are distributed: a few hub proteins vs. an even spread.

One ~55k-row scan of the MV for the distinct (viral protein, human partner)
pairs, one pull of the taxon table, rollup in polars.

Valid PPI filter applied, vh only.

Exports data/viral_protein_human_degree.tsv.
"""

import polars as pl

from drakkar.db import fetch

PAIRS = """
select distinct name2, taxon2, left_value2, accession1
from dataset
where state = 'curated'
  and is_obsolete1 is false
  and is_obsolete2 is false
  and deleted_at is null
  and type2 = 'v'
"""

SPECIES = """
select t.left_value, t.right_value, n.name as virus
from taxon t
join taxon_name n on n.taxon_id = t.taxon_id and n.name_class = 'scientific name'
where t.node_rank = 'species'
"""


def main() -> None:
    pairs = fetch(PAIRS)
    species = fetch(SPECIES)

    matched = pairs.join_where(
        species,
        pl.col("left_value2") >= pl.col("left_value"),
        pl.col("left_value2") <= pl.col("right_value"),
    ).select("name2", "accession1", "virus")
    missing = pairs.join(
        matched.select("name2", pl.col("virus")).unique(),
        on="name2",
        how="anti",
    )
    # fall back to the strain string when taxon2 has no species-rank ancestor
    missing = missing.with_columns(virus=pl.col("taxon2")).select(
        "name2", "accession1", "virus"
    )
    rolled = pl.concat([matched, missing]).unique()

    degree = (
        rolled.group_by("virus", "name2")
        .agg(n_human=pl.col("accession1").n_unique())
        .sort("n_human", descending=True)
    )
    degree.write_csv("data/viral_protein_human_degree.tsv", separator="\t")

    n_prot = degree.height
    total_pairs = degree["n_human"].sum()
    d = degree["n_human"]

    print(f"strain-abstracted viral proteins : {n_prot}")
    print(f"unique (viral, human) pairs       : {total_pairs}")
    print()
    print("degree summary:")
    print(degree["n_human"].describe())
    print()
    for q in (0.5, 0.75, 0.9, 0.95, 0.99):
        print(f"  p{int(q * 100):>2} = {d.quantile(q):.0f}")
    print()

    # concentration: what share of all pairs sits in the top-k viral proteins
    cum = degree.with_columns(cum=pl.col("n_human").cum_sum())
    for frac in (0.01, 0.05, 0.10, 0.25):
        k = max(1, round(n_prot * frac))
        share = cum["cum"][k - 1] / total_pairs
        print(f"  top {frac:>4.0%} ({k:>4} proteins) hold {share:6.1%} of all pairs")
    print()

    # histogram of degree
    breaks = [2, 3, 6, 11, 26, 51, 101, 251]
    labels = ["1", "2", "3-5", "6-10", "11-25", "26-50", "51-100", "101-250", "250+"]
    hist = (
        degree.with_columns(
            bucket=pl.col("n_human").cut(breaks, labels=labels, left_closed=True)
        )
        .group_by("bucket")
        .agg(n_proteins=pl.len(), pairs=pl.col("n_human").sum())
        .sort("bucket")
    )
    print(hist)
    print()
    print("top 20 viral proteins by human degree:")
    print(degree.head(20))


if __name__ == "__main__":
    main()
