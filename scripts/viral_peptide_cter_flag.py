"""Viral peptides (viral-side mappings of 5-20 aa), flagged for C-terminal position.

A *viral peptide* is a `mapping2` entry whose `sequence` is 4-20 aa inclusive
by team convention; this dataset uses a one-time exception of 5-20 aa
inclusive. Peptides are reported on the viral side of a valid `vh`
description as a region sufficient to bind the human partner. This dataset
locates each peptide occurrence on its mature viral protein and flags
whether it ends exactly at that protein's C-terminus (`is_cter`).

A handful of `mapping2` entries carry a `sequence` but an empty `isoforms`
list -- the peptide was curated without ever being positioned (a data gap,
confirmed in the base `descriptions` table, not a view artifact). Those have
no recorded coordinate to trust or fall back to, so every exact match of the
sequence on the protein's sequences is taken as an occurrence.

Peptides are matched against every sequence of their protein: the mature
protein alone when it is one, since a mature protein is a region of a
polyprotein and has no isoforms of its own, and the canonical sequence plus
its isoforms for a viral protein curated full length. A peptide curated on an
isoform is therefore positioned on that isoform.

One row per distinct peptide sequence (a sequence can occur on more than one
mature protein, or more than once on the same one). Columns:

    sequence      -- peptide as curated
    length        -- peptide length
    is_cter       -- True if the peptide ends exactly at the C-terminus of
                     the sequence it sits on, in at least one occurrence
    accession2,
    name2         -- one occurrence's viral protein (polyprotein accession +
                     curator name)
    source        -- the sequence that occurrence sits on: `accession2` for
                     the mature or canonical protein, an isoform accession
                     otherwise
    pep_start,
    pep_stop,
    source_length -- that occurrence's peptide position on `source`, and the
                     length of `source`

The example occurrence is C-terminal when the peptide has any C-terminal
occurrence, and sits on the mature or canonical sequence rather than on an
isoform whenever one does.

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL AND type = 'vh'.

Run: uv run python scripts/viral_peptide_cter_flag.py
Writes data/viral_peptide_cter_flag.tsv.
"""

from collections import Counter

import polars as pl

from drakkar.db import connect
from drakkar.mappings import occurrences as mapping_occurrences
from drakkar.runs import log_run

PEPTIDE_MIN, PEPTIDE_MAX = 5, 20

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
    AND type = 'vh'
"""

with connect() as conn, conn.cursor() as cur:
    cur.execute(f"""
        SELECT protein2_id, accession2, name2, start2, stop2, mapping2
        FROM dataset
        WHERE {VALID}
    """)
    rows = cur.fetchall()

    ids = sorted({r[0] for r in rows})
    cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
    sequences = dict(cur.fetchall())

print(f"{len(rows)} valid vh descriptions, {len(sequences)} viral protein versions")


records = []
stats = Counter()

for protein2_id, acc2, name2, start2, stop2, mapping2 in rows:
    if not mapping2:
        continue
    protein = sequences[protein2_id]
    canonical = protein[acc2]
    mature = canonical[start2 - 1 : stop2]

    # A mature protein is a region of a polyprotein and has no isoforms of its
    # own, so it is the only sequence to match against. A viral protein curated
    # full length may have isoforms, and a peptide may be curated on one.
    full_length = start2 == 1 and stop2 == len(canonical)
    sources = protein if full_length else {acc2: mature}

    for occ in mapping_occurrences(mapping2, sources, PEPTIDE_MIN, PEPTIDE_MAX):
        stats[occ.how] += 1
        source_length = len(sources[occ.accession])
        records.append(
            {
                "sequence": occ.sequence,
                "length": len(occ.sequence),
                "is_cter": occ.stop == source_length,
                "accession2": acc2,
                "name2": name2,
                "source": occ.accession,
                "is_isoform": occ.accession != acc2,
                "pep_start": occ.start,
                "pep_stop": occ.stop,
                "source_length": source_length,
            }
        )

occurrences_df = pl.DataFrame(records)

# The example occurrence: a C-terminal one when the peptide has any, and the
# canonical or mature sequence over an isoform either way. group_by is
# unordered and .first() follows the input order, so both sorts are total --
# otherwise the example, and the row order of the file, differ from one run
# to the next.
example_cols = [
    "accession2",
    "name2",
    "source",
    "pep_start",
    "pep_stop",
    "source_length",
]
example = (
    occurrences_df.sort(
        ["is_cter", "is_isoform", "accession2", "name2", "pep_start"],
        descending=[True, False, False, False, False],
    )
    .group_by("sequence")
    .agg(pl.col(c).first() for c in example_cols)
)

df = (
    occurrences_df.group_by("sequence")
    .agg(pl.col("length").first(), pl.col("is_cter").any())
    .join(example, on="sequence")
    .sort(["length", "sequence"])
)
df.write_csv("data/viral_peptide_cter_flag.tsv", separator="\t")

print(f"peptide occurrences: {dict(stats)}")
print(f"{len(occurrences_df)} occurrences, {len(df)} distinct sequences")
print("-> data/viral_peptide_cter_flag.tsv")
print(f"is_cter = True: {df['is_cter'].sum()}")

log_run(
    "scripts/viral_peptide_cter_flag.py",
    "data/viral_peptide_cter_flag.tsv",
    {
        "valid vh descriptions scanned": len(rows),
        "peptide occurrences resolved": len(occurrences_df),
        "distinct peptide sequences": len(df),
        "is_cter true": int(df["is_cter"].sum()),
    },
)
