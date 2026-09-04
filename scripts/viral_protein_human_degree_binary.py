"""Human-interactor degree of strain-abstracted viral proteins -- binary methods only.

Variant of scripts/viral_protein_human_degree.py. Same grain (a viral protein is
``(NCBI species of taxon2, curated name2)``, strains collapsed), but a
*(viral protein, human protein)* pair is kept only when at least one supporting
description used a **binary** detection method -- one that reports a direct 1:1
physical contact, as opposed to co-complex / affinity capture, proximity
labelling, cross-linking, enzymatic-reaction or nucleic-acid methods.

BINARY_METHOD_IDS below is the curated allowlist (drakkar ``methods.id``). It was
built from the PSI-MI detection-method branch, restricted to methods actually
present in valid vh rows, and covers:

  * two-hybrid and reporter-recruitment assays (Y2H and derivatives, MAPPIT,
    GAL4-VP16, lambda-repressor 2H)
  * protein-fragment complementation (PCA, BiFC, split-ubiquitin, split-luc,
    GPCA, beta-gal/beta-lac/adenylate-cyclase complementation, LUMIER, BRET)
  * resonance energy transfer / homogeneous proximity (FRET, HTRF, AlphaScreen)
  * direct biophysical binding of purified partners (SPR, ITC, BLI, MST, FP,
    FCS, fluorescence spectroscopy, thermal shift, DSC, QCM, biosensor, filter
    binding, competition binding, solid-phase / scintillation proximity, CD)
  * structural methods implying atomic contact (X-ray, NMR, HDX-MS, disulfide
    bond)
  * binding to immobilised bait / display selection (protein & peptide arrays,
    phage / T7 / lambda / yeast display, far-western, ELISA, sandwich immunoassay)

Deliberately EXCLUDED as non-binary: mass spectrometry of complexes, all
co-immunoprecipitation / pull-down / affinity chromatography / TAP, proximity
labelling & BioID, proximity ligation assay, cross-linking studies, EM,
co-sedimentation / co-migration / gel filtration / native PAGE, light /
X-ray / neutron scattering, DLS, EPR, EMSA / mobility shift, antibody array,
western blot, protein three-hybrid, and enzymatic-reaction methods
(phosphorylation, ubiquitination, cleavage, ...).

Borderline calls (protein array, far-western, ELISA, HDX-MS, disulfide bond,
CD, DSC) are included; dropping them does not change the qualitative picture.

Writes the degree table and regenerates the markdown report.

Valid PPI filter applied, vh only.

Run 2026-09-04:
  binary-supported pairs        : 13,626
  strain-abstracted viral proteins with >=1 such pair : 950
  degree: mean 14.3, median 2, max 430, p90 = 27, p99 = 271
  top  1% of viral proteins (10)  hold 25.4% of all pairs
  top 10% (95) hold 75.2%
  exports: data/viral_protein_human_degree_binary.tsv
           reports/viral_protein_human_degree_binary.md
"""

import polars as pl

from drakkar.db import fetch

# drakkar methods.id -- see module docstring for the rationale / grouping
BINARY_METHOD_IDS = (
    # two-hybrid & reporter recruitment
    18,
    254,
    255,
    861,
    1103,
    1259,
    532,
    138,
    459,
    # protein-fragment complementation
    78,
    564,
    99,
    952,
    786,
    1214,
    1263,
    15,
    12,
    11,
    533,
    13,
    # resonance energy transfer / homogeneous proximity
    49,
    357,
    656,
    # direct biophysical binding of purified partners
    95,
    672,
    55,
    718,
    995,
    47,
    46,
    17,
    45,
    983,
    1059,
    1240,
    717,
    44,
    262,
    643,
    87,
    16,
    # structural
    101,
    66,
    852,
    579,
    693,
    265,
    # immobilised-bait binding / display selection
    77,
    69,
    72,
    43,
    56,
    96,
    102,
    42,
    267,
    499,
)

PAIRS = """
select distinct name2, taxon2, left_value2, accession1, method_id
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
        "**Script:** `scripts/viral_protein_human_degree_binary.py`  ",
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
        "## Binary detection methods used",
        "",
        (
            "See `scripts/viral_protein_human_degree_binary.py` docstring for the "
            "full rationale. Allowlist = two-hybrid & reporter-recruitment, "
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
        rolled.filter(pl.col("method_id").is_in(BINARY_METHOD_IDS))
        .select("virus", "name2", "accession1")
        .unique()
    )

    degree = (
        binary_pairs.group_by("virus", "name2")
        .agg(n_human=pl.len())
        .sort("n_human", descending=True)
    )
    degree.write_csv(TSV, separator="\t")

    report = degree_report(
        degree,
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

    print(f"binary-supported pairs: {int(degree['n_human'].sum()):,}")
    print(f"viral proteins        : {degree.height}")
    print(degree.head(15))


if __name__ == "__main__":
    main()
