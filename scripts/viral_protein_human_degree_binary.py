"""Human-interactor degree of strain-abstracted viral proteins -- binary methods only.

Variant of scripts/viral_protein_human_degree.py. Same grain (a viral protein is
``(NCBI species of taxon2, curated name2)``, strains collapsed), but a
*(viral protein, human protein)* pair is kept only when at least one supporting
description used a **binary** detection method -- one that reports a direct 1:1
physical contact, as opposed to co-complex / affinity capture, proximity
labelling, cross-linking, enzymatic-reaction or nucleic-acid methods.

The allowlist of binary methods is `drakkar.binary_methods`: stable PSI-MI
terms (`BINARY_PSIMI_IDS`), matched on `dataset.psimi_id` -- never on the
internal `methods.id`, which is not stable across database versions. See that
module for the grouped allowlist and the rationale for what it excludes.

Valid PPI filter applied, vh only.

Exports data/viral_protein_human_degree_binary.tsv and regenerates the
markdown report.
"""

from datetime import UTC, datetime
from typing import Any

import polars as pl

from drakkar.binary_methods import BINARY_PSIMI_IDS
from drakkar.db import fetch
from drakkar.runs import log_run

PAIRS = """
select distinct name2, taxon2, left_value2, accession1, psimi_id
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

TSV = "data/viral_protein_human_degree_binary.tsv"
MD = "reports/viral_protein_human_degree_binary.md"


def roll_to_virus(pairs: pl.DataFrame, species: pl.DataFrame) -> pl.DataFrame:
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


def summarise(degree: pl.DataFrame) -> dict[str, Any]:
    """The degree distribution -- what the report and the run log both quote."""
    n_human = pl.col("n_human")
    summary = degree.select(
        proteins=pl.len(),
        pairs=n_human.sum(),
        mean=n_human.mean(),
        median=n_human.median(),
        max=n_human.max(),
        p75=n_human.quantile(0.75),
        p90=n_human.quantile(0.9),
        p95=n_human.quantile(0.95),
        p99=n_human.quantile(0.99),
    ).row(0, named=True)

    # the table is sorted, so the top k proteins are its first k rows
    cum = degree.select(n_human.cum_sum())["n_human"]
    summary["top"] = {}
    for frac in (0.01, 0.05, 0.10, 0.25):
        k = max(1, round(summary["proteins"] * frac))
        summary["top"][frac] = (k, cum[k - 1] / summary["pairs"])
    return summary


def degree_report(
    degree: pl.DataFrame, summary: dict[str, Any], title: str, intro: str
) -> str:
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
        f"**Date:** {datetime.now(UTC).date().isoformat()}  ",
        "**Source:** drakkar `dataset` MV, valid PPI filter, `vh` interactions only  ",
        "**Script:** `scripts/viral_protein_human_degree_binary.py`  ",
        f"**Export:** `{TSV}`",
        "",
        intro,
        "",
        "## Result",
        "",
        f"- **{summary['proteins']}** strain-abstracted viral proteins",
        f"- **{summary['pairs']:,}** unique *(viral protein, human protein)* pairs",
        "",
        "### Degree summary",
        "",
        "| statistic | value |",
        "|-----------|------:|",
        f"| mean   | {summary['mean']:.1f} |",
        f"| median | {int(summary['median'])} |",
        f"| p75    | {int(summary['p75'])} |",
        f"| p90    | {int(summary['p90'])} |",
        f"| p95    | {int(summary['p95'])} |",
        f"| p99    | {int(summary['p99'])} |",
        f"| max    | {summary['max']:,} |",
        "",
        "### Concentration",
        "",
        "| top viral proteins | count | share of all pairs |",
        "|--------------------|------:|-------------------:|",
    ]
    for frac, (k, share) in summary["top"].items():
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
        .sort(["pairs", "virus"], descending=[True, False])
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
        "## Binary detection methods used",
        "",
        (
            "See `drakkar/binary_methods.py` for the full rationale. "
            "Allowlist = two-hybrid & reporter-recruitment, "
            "protein-fragment complementation, resonance energy transfer, direct "
            "biophysical binding of purified partners, structural methods, and "
            "immobilised-bait / display binding. Excluded: AP-MS, co-IP, pull-down, "
            "affinity chromatography, TAP, proximity labelling / BioID, proximity "
            "ligation, cross-linking, EM, co-sedimentation / co-migration, "
            "scattering, EMSA, antibody array, western blot, protein three-hybrid, "
            "enzymatic-reaction methods."
        ),
        "",
        "## Caveat",
        "",
        (
            "Binary assays are themselves biased: two-hybrid over-represents pairs "
            "that fold and localise in yeast, and large binary screens (Y2H, GPCA) "
            "still contribute the hubs. This is a direct-contact view of the "
            "curated data, not an unbiased interactome."
        ),
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    pairs = fetch(PAIRS)
    species = fetch(SPECIES)

    rolled = roll_to_virus(pairs, species)

    binary_pairs = (
        rolled.filter(pl.col("psimi_id").is_in(BINARY_PSIMI_IDS))
        .select("virus", "name2", "accession1")
        .unique()
    )

    degree = (
        binary_pairs.group_by("virus", "name2")
        .agg(n_human=pl.len())
        .sort(["n_human", "virus", "name2"], descending=[True, False, False])
    )
    degree.write_csv(TSV, separator="\t")

    summary = summarise(degree)
    report = degree_report(
        degree,
        summary,
        "Human-interactor degree of viral proteins -- binary-method interactions",
        (
            "Only *(viral protein, human protein)* pairs with at least one "
            "supporting description from a **binary** detection method are counted "
            "-- methods that report a direct 1:1 physical contact. Co-complex / "
            "affinity-capture, proximity-labelling, cross-linking and "
            "enzymatic-reaction evidence is excluded."
        ),
    )
    with open(MD, "w") as fh:
        fh.write(report)

    print(f"binary-supported pairs: {summary['pairs']:,}")
    print(f"viral proteins        : {summary['proteins']}")
    print(degree.head(15))

    _, top1 = summary["top"][0.01]
    _, top10 = summary["top"][0.10]
    log_run(
        "scripts/viral_protein_human_degree_binary.py",
        f"{TSV} + {MD}",
        {
            "binary-supported pairs": summary["pairs"],
            "viral proteins": summary["proteins"],
            "degree": (
                f"mean {summary['mean']:.1f}, median {int(summary['median'])}, "
                f"max {summary['max']:,}, p90 = {int(summary['p90'])}, "
                f"p99 = {int(summary['p99'])}"
            ),
            "concentration": (
                f"top 1% hold {top1:.1%} of pairs, top 10% hold {top10:.1%}"
            ),
        },
    )


if __name__ == "__main__":
    main()
