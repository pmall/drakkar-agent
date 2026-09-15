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
no recorded coordinate to trust or fall back to, so every exact substring
match of the sequence on the mature protein is taken as an occurrence.

One row per distinct peptide sequence (a sequence can occur on more than one
mature protein, or more than once on the same one). Columns:

    sequence      -- peptide as curated
    length        -- peptide length
    is_cter       -- True if the peptide ends exactly at the C-terminus of
                     the mature protein in at least one occurrence
    accession2,
    name2         -- one occurrence's mature viral protein (polyprotein
                     accession + curator name) -- a C-terminal occurrence
                     when is_cter is True, an arbitrary one otherwise
    pep_start,
    pep_stop,
    mature_length -- that occurrence's peptide position and protein length

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL AND type = 'vh'.

Run: uv run python scripts/viral_peptide_cter_flag.py
Writes data/viral_peptide_cter_flag.tsv.
"""

from collections import Counter

import polars as pl

from drakkar.db import connect
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


def peptides(mapping):
    """Yield (sequence, recorded start) for every 4-20 aa mapping occurrence.

    `recorded start` is None when the entry carries a sequence but no
    isoforms -- a curation gap where the peptide was never positioned. Those
    have no coordinate to trust or fall back to except a direct substring
    search on the mature protein.
    """
    for entry in mapping or []:
        seq = entry["sequence"]
        if not PEPTIDE_MIN <= len(seq) <= PEPTIDE_MAX:
            continue
        if not entry["isoforms"]:
            yield seq, None
            continue
        for isoform in entry["isoforms"]:
            for occ in isoform["occurrences"]:
                yield seq, occ["start"]


def locate(mature, seq, start):
    """Position(s) of `seq` in `mature`, 1-based.

    `start` is the curated coordinate, normally already relative to the
    mature protein. A few descriptions record it in the polyprotein frame
    instead, so verify it and fall back to an exact search (peptide
    occurrences here are all identity=100, so an exact search is valid).

    `start` is None when no coordinate was curated at all (see `peptides`):
    every substring match is then returned, since there is nothing to
    disambiguate against.

    A handful of descriptions (all Q99IB8 Core) record `start` in the
    polyprotein frame rather than the mature-protein frame, one residue off
    from the mature protein's own start -- there is no substring match to
    find at all, since the recorded sequence includes a residue the mature
    protein doesn't have. There, the curated position is trusted as-is
    rather than dropped, since it still says where on the mature protein the
    peptide was meant to sit.
    """
    hits = [i + 1 for i in range(len(mature)) if mature.startswith(seq, i)]
    if start is None:
        return hits, "no coordinate, searched"
    if mature[start - 1 : start - 1 + len(seq)] == seq:
        return [start], "recorded"
    if hits:
        return [min(hits, key=lambda h: abs(h - start))], "relocated"
    return [start], "recorded despite mismatch"


records = []
stats = Counter()
unresolved = []

for protein2_id, acc2, name2, start2, stop2, mapping2 in rows:
    occurrences = list(peptides(mapping2))
    if not occurrences:
        continue
    mature = sequences[protein2_id][acc2][start2 - 1 : stop2]
    mature_length = len(mature)

    for seq, recorded in occurrences:
        positions, how = locate(mature, seq, recorded)
        stats[how] += 1
        if not positions:
            unresolved.append((acc2, name2, start2, stop2, seq, recorded))
            continue
        for pos in positions:
            pep_stop = pos + len(seq) - 1
            records.append(
                {
                    "sequence": seq,
                    "length": len(seq),
                    "is_cter": pep_stop == mature_length,
                    "accession2": acc2,
                    "name2": name2,
                    "pep_start": pos,
                    "pep_stop": pep_stop,
                    "mature_length": mature_length,
                }
            )

occurrences_df = pl.DataFrame(records)

example_cols = ["accession2", "name2", "pep_start", "pep_stop", "mature_length"]
example = (
    occurrences_df.sort("is_cter", descending=True)
    .group_by("sequence")
    .agg(pl.col(c).first() for c in example_cols)
)

df = (
    occurrences_df.group_by("sequence")
    .agg(pl.col("length").first(), pl.col("is_cter").any())
    .join(example, on="sequence")
    .sort("length", descending=False)
)
df.write_csv("data/viral_peptide_cter_flag.tsv", separator="\t")

print(f"peptide occurrences: {dict(stats)}")
print(f"{len(occurrences_df)} occurrences, {len(df)} distinct sequences")
print("-> data/viral_peptide_cter_flag.tsv")
print(f"is_cter = True: {df['is_cter'].sum()}")
for u in unresolved:
    print("  unresolved:", u)

log_run(
    "scripts/viral_peptide_cter_flag.py",
    "data/viral_peptide_cter_flag.tsv",
    {
        "valid vh descriptions scanned": len(rows),
        "peptide occurrences resolved": len(occurrences_df),
        "distinct peptide sequences": len(df),
        "unresolved": len(unresolved),
        "is_cter true": int(df["is_cter"].sum()),
    },
)
