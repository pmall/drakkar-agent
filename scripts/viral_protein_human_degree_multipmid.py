"""Human-interactor degree of strain-abstracted viral proteins -- multi-publication only.

Variant of scripts/viral_protein_human_degree.py. Same grain (a viral protein is
``(NCBI species of taxon2, curated name2)``, strains collapsed), but an
interaction is kept only when the *(viral protein, human protein)* pair is
described by **more than one distinct PubMed id** (counted after strain
abstraction: two single-pmid curations on two strains of the same virus count as
two publications for the collapsed pair).

Writes the degree table and regenerates the markdown report.

Valid PPI filter applied, vh only.

Exports data/viral_protein_human_degree_multipmid.tsv and regenerates the
markdown report.
"""

import polars as pl

from drakkar.db import fetch

PAIRS = """
select distinct name2, taxon2, left_value2, accession1, pmid
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

TSV = "data/viral_protein_human_degree_multipmid.tsv"
MD = "reports/viral_protein_human_degree_multipmid.md"


def roll_to_virus(pairs: pl.DataFrame, species: pl.DataFrame) -> pl.DataFrame:
    """Attach the NCBI species name to each row; fall back to taxon2 when none."""
    p = pairs.with_row_index("_rid")
    matched = p.join_where(
        species,
        pl.col("left_value2") >= pl.col("left_value"),
        pl.col("left_value2") <= pl.col("right_value"),
    ).select("_rid", "virus")
    return (
        p.join(matched, on="_rid", how="left")
        .with_columns(virus=pl.col("virus").fill_null(pl.col("taxon2")))
        .drop("_rid")
    )


def degree_report(degree: pl.DataFrame, title: str, intro: str) -> str:
    n_prot = degree.height
    total = int(degree["n_human"].sum())
    d = degree["n_human"]

    qs = {q: int(d.quantile(q)) for q in (0.5, 0.75, 0.9, 0.95, 0.99)}
    cum = degree.with_columns(cum=pl.col("n_human").cum_sum())
    conc = {}
    for frac in (0.01, 0.05, 0.10, 0.25):
        k = max(1, round(n_prot * frac))
        conc[frac] = (k, cum["cum"][k - 1] / total)

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

    lines = [
        f"# {title}",
        "",
        "**Date:** 2026-09-04  ",
        "**Source:** drakkar `dataset` MV, valid PPI filter, `vh` interactions only  ",
        "**Script:** `scripts/viral_protein_human_degree_multipmid.py`  ",
        f"**Export:** `{TSV}`",
        "",
        intro,
        "",
        "## Result",
        "",
        f"- **{n_prot}** strain-abstracted viral proteins",
        f"- **{total:,}** unique *(viral protein, human protein)* pairs",
        "",
        "### Degree summary",
        "",
        "| statistic | value |",
        "|-----------|------:|",
        f"| mean   | {d.mean():.1f} |",
        f"| median | {int(d.median())} |",
        f"| p75    | {qs[0.75]} |",
        f"| p90    | {qs[0.9]} |",
        f"| p95    | {qs[0.95]} |",
        f"| p99    | {qs[0.99]} |",
        f"| max    | {int(d.max()):,} |",
        "",
        "### Concentration",
        "",
        "| top viral proteins | count | share of all pairs |",
        "|--------------------|------:|-------------------:|",
    ]
    for frac, (k, share) in conc.items():
        lines.append(f"| top {frac:.0%} | {k} | {share:.1%} |")
    lines += [
        "",
        "### Degree distribution",
        "",
        "| interactors | # viral proteins | pairs contributed |",
        "|-------------|-----------------:|------------------:|",
    ]
    for row in hist.iter_rows(named=True):
        lines.append(f"| {row['bucket']} | {row['n_proteins']:,} | {row['pairs']:,} |")
    lines += [
        "",
        "## Top 20 viral proteins by human degree",
        "",
        "| virus | protein | human interactors |",
        "|-------|---------|------------------:|",
    ]
    for row in degree.head(20).iter_rows(named=True):
        lines.append(f"| {row['virus']} | {row['name2']} | {row['n_human']:,} |")

    pv = (
        degree.group_by("virus")
        .agg(n_prot=pl.len(), pairs=pl.col("n_human").sum())
        .sort("pairs", descending=True)
        .head(10)
    )
    lines += [
        "",
        "## Top viruses by total interactions",
        "",
        "| virus | viral proteins | pairs |",
        "|-------|---------------:|------:|",
    ]
    for row in pv.iter_rows(named=True):
        lines.append(f"| {row['virus']} | {row['n_prot']} | {row['pairs']:,} |")

    lines += [
        "",
        "## Caveat",
        "",
        (
            "Requiring two independent publications removes most single-paper "
            "high-throughput screens, so this view is closer to a "
            "reproducibility-weighted interactome -- but it is still shaped by which "
            "interactions attract repeat study (well-known oncoproteins, druggable "
            "targets) rather than by biology alone."
        ),
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    pairs = fetch(PAIRS)
    species = fetch(SPECIES)

    rolled = roll_to_virus(pairs, species)

    # keep only (virus, name2, human) pairs backed by >1 distinct pmid
    multi = (
        rolled.group_by("virus", "name2", "accession1")
        .agg(n_pmid=pl.col("pmid").n_unique())
        .filter(pl.col("n_pmid") > 1)
    )

    degree = (
        multi.group_by("virus", "name2")
        .agg(n_human=pl.len())
        .sort("n_human", descending=True)
    )
    degree.write_csv(TSV, separator="\t")

    report = degree_report(
        degree,
        "Human-interactor degree of viral proteins -- multi-publication interactions",
        (
            "Only *(viral protein, human protein)* pairs described by **more than "
            "one distinct PubMed id** are counted (pmids pooled across strains of "
            "the same virus). This drops single-paper interactions, including most "
            "one-off high-throughput screens."
        ),
    )
    with open(MD, "w") as fh:
        fh.write(report)

    print(f"pairs kept (>1 pmid): {int(degree['n_human'].sum()):,}")
    print(f"viral proteins      : {degree.height}")
    print(degree.head(15))


if __name__ == "__main__":
    main()
