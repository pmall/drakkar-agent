"""FASTA of viral mature proteins carrying a peptide (mapping of 4-20 aa).

Built for multiple alignment of the mature proteins, then extraction of the
conserved residues under the peptides.

One record per peptide: the sequence is the whole mature protein, repeated when a
protein carries several peptides. Header:

    >accession|start2-stop2|name2|pepstart-pepstop|taxon

Peptide coordinates are 1-based inclusive on the emitted sequence.

Run: uv run python scripts/viral_peptide_mature_fasta.py

Produced 2026-09-02 -> data/viral_mature_with_peptides.fasta:
    116,015 valid vh descriptions
    2,467 peptide occurrences: 2,446 at the recorded coordinate, 19 relocated,
    2 skipped (see `locate` below)
    1,264 records over 625 mature proteins, 571 distinct sequences
"""

from collections import Counter
from pathlib import Path

from drakkar.db import connect

OUT = Path("data/viral_mature_with_peptides.fasta")
PEPTIDE_MIN, PEPTIDE_MAX = 4, 20
WRAP = 60

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
    AND type = 'vh'
"""

with connect() as conn, conn.cursor() as cur:
    cur.execute(f"""
        SELECT protein2_id, accession2, name2, start2, stop2, taxon2, mapping2
        FROM dataset
        WHERE {VALID}
    """)
    rows = cur.fetchall()

    ids = sorted({r[0] for r in rows})
    cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
    sequences = dict(cur.fetchall())

print(f"{len(rows)} valid vh descriptions, {len(sequences)} viral protein versions")


def peptides(mapping):
    """Yield (sequence, recorded start) for every 4-20 aa mapping occurrence."""
    for entry in mapping:
        seq = entry["sequence"]
        if not PEPTIDE_MIN <= len(seq) <= PEPTIDE_MAX:
            continue
        for isoform in entry["isoforms"]:
            for occ in isoform["occurrences"]:
                yield seq, occ["start"]


def locate(mature, seq, start):
    """Position of `seq` in `mature`, 1-based, or None.

    `start` is the curated coordinate, normally already relative to the mature
    protein. A few descriptions record it in the polyprotein frame instead (all
    Q99IB8 Core, start2=2), and one is a plain off-by-one (P06935 Core), so verify
    it and fall back to an exact search. Viral peptide occurrences are all
    identity=100, which makes the search exact.
    """
    if mature[start - 1 : start - 1 + len(seq)] == seq:
        return start, "recorded"
    hits = [i + 1 for i in range(len(mature)) if mature.startswith(seq, i)]
    if not hits:
        return None, "not found"
    return min(hits, key=lambda h: abs(h - start)), "relocated"


records = {}  # header -> mature sequence
stats = Counter()
unresolved = []

for protein_id, acc, name, start, stop, taxon, mapping in rows:
    occurrences = list(peptides(mapping))
    if not occurrences:
        continue
    mature = sequences[protein_id][acc][start - 1 : stop]

    for seq, recorded in occurrences:
        pos, how = locate(mature, seq, recorded)
        stats[how] += 1
        if pos is None:
            unresolved.append((acc, name, start, stop, seq, recorded))
            continue
        header = f"{acc}|{start}-{stop}|{name}|{pos}-{pos + len(seq) - 1}|{taxon}"
        records[header] = mature

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w") as fh:
    for header in sorted(records):
        seq = records[header]
        fh.write(f">{header}\n")
        fh.writelines(seq[i : i + WRAP] + "\n" for i in range(0, len(seq), WRAP))

print(f"peptide occurrences: {dict(stats)}")
print(
    f"{len(records)} records, {len(set(records.values()))} distinct sequences -> {OUT}"
)
for u in unresolved:
    print("  unresolved:", u)
