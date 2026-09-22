"""Summary counts for viral peptides (viral-side mappings of 4-20 aa).

A *viral peptide* is a `mapping2` entry whose `sequence` is 4-20 aa inclusive
(team convention), reported on the viral side of a valid `vh` description as a
region sufficient to bind the human partner.

Entries are read through `drakkar.mappings.occurrences`, which places each one
on the sequence it was curated on: the mature protein for a mature viral
protein, the canonical sequence and its isoforms for one curated full length.
A peptide that sits nowhere on its protein is therefore not counted --
`scripts/mapping_corrections.py` audits those.

Reports, over valid `vh` descriptions:
  - unique peptide sequences
  - unique source viral (mature) proteins, keyed on (accession2, start2, stop2)
  - unique human targets (accession1)
  - unique (peptide sequence, human target) associations

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL AND type = 'vh'.

Run: uv run python scripts/viral_peptide_summary.py
Writes the per-sequence table to data/viral_peptide_summary.tsv.
"""

from collections import defaultdict

import polars as pl

from drakkar.db import connect
from drakkar.mappings import occurrences
from drakkar.runs import log_run

type ViralProtein = tuple[str, int, int]  # accession2, start2, stop2

OUTPUT = "data/viral_peptide_summary.tsv"

PEPTIDE_MIN, PEPTIDE_MAX = 4, 20

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
    AND type = 'vh'
    AND json_array_length(mapping2) > 0
"""

with connect() as conn, conn.cursor() as cur:
    cur.execute(f"""
        SELECT accession1, protein2_id, accession2, start2, stop2, mapping2
        FROM dataset
        WHERE {VALID}
    """)
    rows = cur.fetchall()

    ids = sorted({r[1] for r in rows})
    cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
    sequences = dict(cur.fetchall())

viral_proteins_of: dict[str, set[ViralProtein]] = defaultdict(set)
human_targets_of: dict[str, set[str]] = defaultdict(set)

for acc1, protein2_id, acc2, start2, stop2, mapping2 in rows:
    protein = sequences[protein2_id]
    canonical = protein[acc2]

    # A mature protein is a region of a polyprotein and has no isoforms of its
    # own; a viral protein curated full length may have some, and a peptide may
    # be curated on one.
    full_length = start2 == 1 and stop2 == len(canonical)
    sources = protein if full_length else {acc2: canonical[start2 - 1 : stop2]}

    for occ in occurrences(mapping2, sources, PEPTIDE_MIN, PEPTIDE_MAX):
        viral_proteins_of[occ.sequence].add((acc2, start2, stop2))
        human_targets_of[occ.sequence].add(acc1)

# per-sequence table
table = pl.DataFrame(
    {
        "sequence": list(viral_proteins_of),
        "n_source_viral_proteins": [len(v) for v in viral_proteins_of.values()],
        "n_human_targets": [len(human_targets_of[s]) for s in viral_proteins_of],
    }
).sort(["n_human_targets", "sequence"], descending=[True, False])
table.write_csv(OUTPUT, separator="\t")

n_viral_proteins = len({vp for vps in viral_proteins_of.values() for vp in vps})
n_human_targets = len({acc for accs in human_targets_of.values() for acc in accs})
n_pep_target = sum(len(targets) for targets in human_targets_of.values())

print(f"valid vh descriptions with a mapping:     {len(rows):,}")
print(f"unique viral peptide sequences:           {len(viral_proteins_of):,}")
print(f"unique source viral (mature) proteins:    {n_viral_proteins:,}")
print(f"unique human targets:                     {n_human_targets:,}")
print(f"unique (peptide, human target) assoc:     {n_pep_target:,}")

log_run(
    "scripts/viral_peptide_summary.py",
    OUTPUT,
    {
        "valid vh descriptions with a mapping": len(rows),
        "unique viral peptide sequences": len(viral_proteins_of),
        "unique source viral (mature) proteins": n_viral_proteins,
        "unique human targets": n_human_targets,
        "unique (peptide, human target) associations": n_pep_target,
    },
)
