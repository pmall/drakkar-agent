"""Summary counts for viral peptides (viral-side mappings of 4-20 aa).

A *viral peptide* is a `mapping2` entry whose `sequence` is 4-20 aa inclusive
(team convention), reported on the viral side of a valid `vh` description as a
region sufficient to bind the human partner.

Reports, over valid `vh` descriptions:
  - unique peptide sequences
  - unique source viral (mature) proteins, keyed on (accession2, start2, stop2)
  - unique human targets (accession1)
  - unique (peptide sequence, human target) associations

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL.

Run: uv run python scripts/viral_peptide_summary.py
Writes the per-sequence table to data/viral_peptide_summary.tsv.
"""

import polars as pl

from drakkar.db import connect

PEPTIDE_MIN, PEPTIDE_MAX = 4, 20

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
    AND type = 'vh'
"""

with connect() as conn, conn.cursor() as cur:
    cur.execute(f"""
        SELECT accession1, accession2, start2, stop2, mapping2
        FROM dataset
        WHERE {VALID}
    """)
    rows = cur.fetchall()

sequences = set()
viral_proteins = set()
human_targets = set()
pep_target = set()
seq_rows = {}  # sequence -> [set of viral proteins, set of human targets]

for acc1, acc2, start2, stop2, mapping2 in rows:
    vp = (acc2, start2, stop2)
    for entry in mapping2 or []:
        seq = entry["sequence"]
        if not PEPTIDE_MIN <= len(seq) <= PEPTIDE_MAX:
            continue
        sequences.add(seq)
        viral_proteins.add(vp)
        human_targets.add(acc1)
        pep_target.add((seq, acc1))
        d = seq_rows.setdefault(seq, [set(), set()])
        d[0].add(vp)
        d[1].add(acc1)

# per-sequence table
table = pl.DataFrame(
    {
        "sequence": list(seq_rows),
        "n_source_viral_proteins": [len(v[0]) for v in seq_rows.values()],
        "n_human_targets": [len(v[1]) for v in seq_rows.values()],
    }
).sort("n_human_targets", descending=True)
table.write_csv("data/viral_peptide_summary.tsv", separator="\t")

print(f"valid vh descriptions scanned:            {len(rows):,}")
print(f"unique viral peptide sequences:           {len(sequences):,}")
print(f"unique source viral (mature) proteins:    {len(viral_proteins):,}")
print(f"unique human targets:                     {len(human_targets):,}")
print(f"unique (peptide, human target) assoc:     {len(pep_target):,}")
