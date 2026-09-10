"""All human/viral interactions with human and mature viral protein sequences.

One row = one vh interaction, i.e. one (human accession, viral polyprotein
accession + start2/stop2) pair. All supporting descriptions (valid PPI filter)
are collapsed into the row: distinct pmids and PSI-MI ids are aggregated as
comma-separated cells, sequences are attached, and a binary flag is computed.

Columns:
  accession1, name1, human_seq, human_len,
  accession2, name2, start2, stop2, viral_seq, viral_len,
  len_sum_squared = (human_len + viral_len)^2,
  taxon2 (strain-level scientific name), ncbi_taxon_id2,
  n_descriptions, n_pmid, pmids (", "-separated, numeric order),
  psimi_ids (", "-separated, lexicographic order),
  method_ids (", "-separated drakkar methods.id, numeric order),
  is_binary (TRUE when at least one supporting description used a binary
  detection method).

Sequences come from the exact protein rows referenced by the descriptions
(proteins.id = dataset.protein1_id / protein2_id): human_seq is the canonical
full-length sequence, viral_seq is the mature-protein substring of the
polyprotein (start2..stop2). A pair always points to a single protein-row
combo, so sequences are unambiguous per row.

Binary methods: is_binary is TRUE when at least one supporting description
has a psimi_id in drakkar.binary_methods.BINARY_PSIMI_IDS (stable PSI-MI
terms -- matching is never done on internal methods.id). See that module for
the full grouped rationale and use drakkar.binary_methods.format_binary_methods
to render the list; excluded as non-binary are AP-MS / mass spec of complexes,
co-IP, pull-down, affinity chromatography, TAP, proximity labelling / BioID,
proximity ligation, cross-linking, EM, co-sedimentation / co-migration,
scattering, EMSA, antibody array, western blot, protein three-hybrid and
enzymatic-reaction methods.

Valid PPI filter applied, vh only.
"""

from drakkar.binary_methods import BINARY_PSIMI_IDS
from drakkar.db import export
from drakkar.runs import log_run

VALID = """
    d.state = 'curated'
    AND d.is_obsolete1 IS FALSE
    AND d.is_obsolete2 IS FALSE
    AND d.deleted_at IS NULL
    AND d.type = 'vh'
"""

QUERY = f"""
SELECT
    g.accession1,
    g.name1,
    g.human_seq,
    g.human_len,
    g.accession2,
    g.name2,
    g.start2,
    g.stop2,
    g.viral_seq,
    g.viral_len,
    (g.human_len + g.viral_len) * (g.human_len + g.viral_len)
        AS len_sum_squared,
    g.taxon2,
    g.ncbi_taxon_id2,
    g.n_descriptions,
    g.n_pmid,
    g.pmids,
    g.psimi_ids,
    g.method_ids,
    g.is_binary
FROM (
    SELECT
        d.accession1,
        max(d.name1) AS name1,
        max(p1.sequences ->> d.accession1) AS human_seq,
        max(length(p1.sequences ->> d.accession1)) AS human_len,
        d.accession2,
        max(d.name2) AS name2,
        d.start2,
        d.stop2,
        max(substring(
            p2.sequences ->> d.accession2 FROM d.start2 FOR d.stop2 - d.start2 + 1
        )) AS viral_seq,
        max(d.stop2 - d.start2 + 1) AS viral_len,
        max(d.taxon2) AS taxon2,
        max(d.ncbi_taxon_id2) AS ncbi_taxon_id2,
        count(*) AS n_descriptions,
        count(DISTINCT d.pmid) AS n_pmid,
        array_to_string(array_agg(DISTINCT d.pmid ORDER BY d.pmid), ', ')
            AS pmids,
        array_to_string(
            array_agg(DISTINCT d.psimi_id ORDER BY d.psimi_id), ', '
        ) AS psimi_ids,
        array_to_string(
            array_agg(DISTINCT d.method_id ORDER BY d.method_id), ', '
        ) AS method_ids,
        bool_or(d.psimi_id IN ({", ".join(f"'{p}'" for p in BINARY_PSIMI_IDS)}))
            AS is_binary
    FROM dataset d
    JOIN proteins p1 ON p1.id = d.protein1_id
    JOIN proteins p2 ON p2.id = d.protein2_id
    WHERE {VALID}
    GROUP BY d.accession1, d.accession2, d.start2, d.stop2
) g
ORDER BY g.accession1, g.accession2, g.start2, g.stop2
"""

TSV = "data/vh_interactions_with_sequences.tsv"


def main() -> None:
    df = export(QUERY, TSV)
    n_bin = int(df["is_binary"].sum())
    print(f"interactions        : {df.height:,}")
    print(f"binary interactions : {n_bin:,}")
    print(df.head(5))
    log_run(
        "scripts/vh_interactions_with_sequences.py",
        TSV,
        {"interactions": df.height, "binary interactions": n_bin},
    )


if __name__ == "__main__":
    main()
