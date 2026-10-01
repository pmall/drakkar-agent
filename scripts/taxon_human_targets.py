"""Human targets of a viral taxon.

Takes an NCBI taxon id and exports one row per human protein interacting
with a viral protein under that taxon, i.e. `accession1` of the valid `vh`
descriptions whose `left_value2` falls in the taxon's nested-set bounds, so
every strain below it counts:

  accession, name, descriptions, interactions, viral_proteins

`name` is the target's distinct `name1` values, `", "`-separated should
there be several. `descriptions` is the number of those descriptions
involving the target, `interactions` the number of distinct viral proteins
among them -- an interaction is a protein pair, and the target is one side
of each. `viral_proteins` lists those viral proteins as
`accession:start-stop` triples, `", "`-separated.

Evidence: an interaction (human protein, viral protein) is supported by its
descriptions; its evidence is the number of distinct publications (`pmid`)
and of distinct detection methods (`psimi_id`) among them.
`--min-publications` and `--min-methods` keep the interactions reaching
either threshold -- an OR, so `--min-publications 2 --min-methods 2` is at
least two publications or at least two methods. With only one of them, that
one applies alone; with neither, every interaction is kept. A target with no
kept interaction is dropped, and `descriptions`, `interactions` and
`viral_proteins` count the kept interactions only.

The accession column feeds `scripts/human_target_peptides.py`, which exports
the peptides binding these targets.

Valid PPI filter applied (curated, both accessions live, current revision).

Run: uv run python scripts/taxon_human_targets.py --ncbi-taxon-id 186536
     --output data/ebolavirus_targets.tsv [--min-publications 2 --min-methods 2]
"""

from __future__ import annotations

import argparse
from pathlib import Path

from drakkar.db import fetch
from drakkar.runs import log_run

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

# Interactions with enough evidence: {evidence} is the OR of the thresholds
# given, or TRUE when there are none.
TARGETS_QUERY = f"""
WITH supported AS (
    SELECT d.accession1, d.accession2, d.start2, d.stop2
    FROM dataset d
    WHERE {VALID}
      AND d.type = 'vh'
      AND d.left_value2 BETWEEN %(left)s AND %(right)s
    GROUP BY d.accession1, d.accession2, d.start2, d.stop2
    HAVING {{evidence}}
)
SELECT
    d.accession1 AS accession,
    string_agg(DISTINCT d.name1, ', ' ORDER BY d.name1) AS name,
    count(*) AS descriptions,
    count(DISTINCT (d.accession2, d.start2, d.stop2)) AS interactions,
    string_agg(
        DISTINCT d.accession2 || ':' || d.start2 || '-' || d.stop2, ', '
    ) AS viral_proteins
FROM dataset d
JOIN supported s
  ON s.accession1 = d.accession1
 AND s.accession2 = d.accession2
 AND s.start2 = d.start2
 AND s.stop2 = d.stop2
WHERE {VALID}
  AND d.type = 'vh'
  AND d.left_value2 BETWEEN %(left)s AND %(right)s
GROUP BY d.accession1
ORDER BY descriptions DESC, accession
"""


def _evidence(min_publications: int | None, min_methods: int | None) -> str:
    """The HAVING condition on an interaction's evidence: thresholds OR-ed."""
    conditions = []
    if min_publications is not None:
        conditions.append(f"count(DISTINCT d.pmid) >= {min_publications}")
    if min_methods is not None:
        conditions.append(f"count(DISTINCT d.psimi_id) >= {min_methods}")
    return " OR ".join(conditions) or "TRUE"


def _threshold(value: int | None) -> str:
    return "none" if value is None else str(value)


def main(
    ncbi_taxon_id: int,
    min_publications: int | None,
    min_methods: int | None,
    output: Path,
) -> None:
    taxon = fetch(TAXON_QUERY, {"ncbi_taxon_id": ncbi_taxon_id})
    if taxon.is_empty():
        raise SystemExit(f"no taxon with NCBI id {ncbi_taxon_id}")
    left, right, name = taxon.row(0)

    query = TARGETS_QUERY.format(evidence=_evidence(min_publications, min_methods))
    targets = fetch(query, {"left": left, "right": right})
    if targets.is_empty():
        raise SystemExit(f"no valid vh interaction under {name} ({ncbi_taxon_id})")

    output.parent.mkdir(parents=True, exist_ok=True)
    targets.write_csv(output, separator="\t")

    log_run(
        "scripts/taxon_human_targets.py",
        str(output),
        {
            "taxon": f"{name} (NCBI {ncbi_taxon_id})",
            "minimum publications": _threshold(min_publications),
            "minimum methods": _threshold(min_methods),
            "human targets": targets.height,
            "target descriptions": int(targets["descriptions"].sum()),
            "target interactions": int(targets["interactions"].sum()),
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Export the human targets of a viral taxon."
    )
    parser.add_argument(
        "--ncbi-taxon-id", type=int, required=True, help="NCBI taxon id of the virus"
    )
    parser.add_argument(
        "--output", type=Path, required=True, help="TSV file the targets go to"
    )
    parser.add_argument(
        "--min-publications",
        type=int,
        help="keep interactions supported by at least this many publications",
    )
    parser.add_argument(
        "--min-methods",
        type=int,
        help="keep interactions supported by at least this many detection methods"
        " (OR-ed with --min-publications)",
    )
    args = parser.parse_args()
    main(
        args.ncbi_taxon_id,
        args.min_publications,
        args.min_methods,
        args.output,
    )
