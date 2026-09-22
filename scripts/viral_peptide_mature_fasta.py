"""FASTA of viral mature proteins carrying a peptide (mapping of 4-20 aa).

Built for multiple alignment of the mature proteins, then extraction of the
conserved residues under the peptides.

One record per peptide: the sequence is the whole mature protein, repeated when a
protein carries several peptides. Header:

    >accession|start2-stop2|name2|pepstart-pepstop|taxon

Peptide coordinates are 1-based inclusive on the emitted sequence. Entries are
read through `drakkar.mappings.occurrences` against the mature protein alone --
the sequence the record carries, so the coordinates are in its frame. An entry
recorded on an isoform is therefore searched for on the mature protein, and one
that is nowhere on it yields no record.

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL AND type = 'vh'.

Run: uv run python scripts/viral_peptide_mature_fasta.py
Writes data/viral_mature_with_peptides.fasta.
"""

from collections import Counter
from pathlib import Path

from drakkar.db import connect
from drakkar.mappings import How, occurrences
from drakkar.runs import log_run

OUT = Path("data/viral_mature_with_peptides.fasta")
PEPTIDE_MIN, PEPTIDE_MAX = 4, 20
WRAP = 60

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
        SELECT protein2_id, accession2, name2, start2, stop2, taxon2, mapping2
        FROM dataset
        WHERE {VALID}
    """)
    rows = cur.fetchall()

    ids = sorted({r[0] for r in rows})
    cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
    sequences = dict(cur.fetchall())

print(f"{len(rows)} valid vh descriptions with a mapping, {len(sequences)} proteins")

records: dict[str, str] = {}  # header -> mature sequence
mature_proteins: set[tuple[str, int, int]] = set()
stats: Counter[How] = Counter()

for protein_id, acc, name, start, stop, taxon, mapping in rows:
    mature = sequences[protein_id][acc][start - 1 : stop]
    for occ in occurrences(mapping, {acc: mature}, PEPTIDE_MIN, PEPTIDE_MAX):
        stats[occ.how] += 1
        header = f"{acc}|{start}-{stop}|{name}|{occ.start}-{occ.stop}|{taxon}"
        records[header] = mature
        mature_proteins.add((acc, start, stop))

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

log_run(
    "scripts/viral_peptide_mature_fasta.py",
    str(OUT),
    {
        "valid vh descriptions with a mapping": len(rows),
        "peptide occurrences": sum(stats.values()),
        "records": len(records),
        "mature proteins": len(mature_proteins),
        "distinct sequences": len(set(records.values())),
    }
    | {f"occurrences {how}": n for how, n in sorted(stats.items())},
)
