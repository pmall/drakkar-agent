"""Mappings that require correction, on valid descriptions (hh and vh).

Audits every `mapping1` / `mapping2` entry for malformation and coordinate
problems. A valid mapping entry is an object `{sequence, isoforms}` with a
non-empty `sequence`, at least one isoform `{accession, occurrences}`, and every
isoform at least one occurrence `{start, stop, identity}`. The isoform accession
must be the protein's canonical accession or one of its isoforms present in
`proteins.sequences`, and never a non-canonical isoform on a mature protein.

Occurrence coordinates are checked against the protein as curated (the mature
protein `startN..stopN` for the canonical accession, the whole isoform
sequence otherwise). An occurrence recorded with identity 100 must match the
mapping sequence exactly; one below 100 must reach at least the recorded
identity minus one point against the region (ungapped or gapped match count).
`X` residues, numbers stored as strings, and isoforms or repeat positions left
out by the curator are allowed and not reported.

One row per problem. Columns:

    stable_id, type, side  -- description and mapping side (1 or 2)
    accession, name,
    start, stop            -- protein (mature protein on the viral side)
    is_peptide             -- mapping sequence is 4-20 aa inclusive
    problem                -- empty_isoforms, empty_occurrences,
                              duplicate_occurrence, repeated_isoform_block,
                              malformed_entry, invalid_isoform,
                              coordinates_mismatch, coordinates_out_of_bounds
    sequence               -- mapping sequence
    isoform                -- isoform accession concerned (when any)
    recorded               -- recorded occurrence(s)
    finding                -- what is actually there / where the sequence is

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL.

Run: uv run python scripts/mapping_corrections.py
Writes data/mapping_corrections.tsv.
"""

import json
import re
from collections import Counter
from difflib import SequenceMatcher

import polars as pl

from drakkar.db import connect
from drakkar.runs import log_run

OUTPUT = "data/mapping_corrections.tsv"

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
"""

with connect() as conn, conn.cursor() as cur:
    cur.execute(f"""
        SELECT stable_id, type,
               protein1_id, accession1, name1, start1, stop1, mapping1::text,
               protein2_id, accession2, name2, start2, stop2, mapping2::text
        FROM dataset
        WHERE {VALID}
          AND (json_array_length(mapping1) > 0 OR json_array_length(mapping2) > 0)
    """)
    rows = cur.fetchall()
    ids = sorted({r[2] for r in rows} | {r[8] for r in rows})
    cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
    sequences = dict(cur.fetchall())

print(f"{len(rows)} valid descriptions with at least one mapping")


def number(value):
    """Numeric value of an occurrence field (strings tolerated), else None."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return value
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def aligned_identity(region, seq):
    """Identity (%) of `seq` against `region`, matches over the longer length.

    Residue-by-residue for equal lengths; when that falls short (or lengths
    differ) a gapped match count from difflib, which never overestimates.
    """
    longest = max(len(region), len(seq))
    ungapped = 0
    if len(region) == len(seq):
        ungapped = sum(a == b for a, b in zip(region, seq, strict=True))
    matcher = SequenceMatcher(None, region, seq, autojunk=False)
    gapped = sum(block.size for block in matcher.get_matching_blocks())
    return 100 * max(ungapped, gapped) / longest


def positions(frame, seq):
    return [i + 1 for i in range(len(frame)) if frame.startswith(seq, i)]


def span(starts, length):
    return ", ".join(f"{s}-{s + length - 1}" for s in starts)


records = []
stats = Counter()


def audit(stable_id, type_, side, protein_id, acc, name, start, stop, raw):
    base = {
        "stable_id": stable_id.strip(),
        "type": type_.strip(),
        "side": side,
        "accession": acc,
        "name": name,
        "start": start,
        "stop": stop,
    }

    def report(problem, seq, isoform="", recorded="", finding=""):
        records.append(
            base
            | {
                "is_peptide": isinstance(seq, str) and 4 <= len(seq) <= 20,
                "problem": problem,
                "sequence": seq if isinstance(seq, str) else json.dumps(seq),
                "isoform": isoform,
                "recorded": recorded,
                "finding": finding,
            }
        )

    mapping = json.loads(raw)
    if not isinstance(mapping, list):
        report("malformed_entry", None, finding="mapping is not an array")
        return
    protein = sequences[protein_id]
    canonical = protein[acc]
    mature = canonical[start - 1 : stop]
    full_length = start == 1 and stop == len(canonical)

    for entry in mapping:
        stats["entries"] += 1
        if not isinstance(entry, dict) or set(entry) != {"sequence", "isoforms"}:
            report(
                "malformed_entry",
                entry,
                finding="entry keys are not sequence, isoforms",
            )
            continue
        seq, isoforms = entry["sequence"], entry["isoforms"]
        if not isinstance(seq, str) or not seq:
            report("malformed_entry", seq, finding="sequence missing or empty")
            continue
        if not isinstance(isoforms, list) or not isoforms:
            poly = positions(canonical, seq)
            if not full_length and poly:
                finding = (
                    f"not on the mature protein; found at {span(poly, len(seq))} "
                    f"on {acc}"
                )
            elif poly:
                finding = f"found at {span(poly, len(seq))} on {acc}"
            else:
                finding = f"sequence not found on {acc}"
            report("empty_isoforms", seq, finding=finding)
            continue

        blocks = Counter(
            iso.get("accession") for iso in isoforms if isinstance(iso, dict)
        )
        for iso_acc, n in blocks.items():
            if n > 1:
                report(
                    "repeated_isoform_block",
                    seq,
                    iso_acc,
                    finding=f"{iso_acc} appears in {n} isoform blocks",
                )

        for iso in isoforms:
            if not isinstance(iso, dict) or set(iso) != {"accession", "occurrences"}:
                report(
                    "malformed_entry",
                    seq,
                    finding="isoform keys are not accession, occurrences",
                )
                continue
            iso_acc, occurrences = iso["accession"], iso["occurrences"]
            if (
                not isinstance(iso_acc, str)
                or not re.fullmatch(re.escape(acc) + r"(-\d+)?", iso_acc)
                or iso_acc not in protein
            ):
                report(
                    "invalid_isoform",
                    seq,
                    str(iso_acc),
                    finding="not an isoform of the protein",
                )
                continue
            if iso_acc != acc and not full_length:
                report(
                    "invalid_isoform",
                    seq,
                    iso_acc,
                    finding="non-canonical isoform on a mature protein",
                )
                continue
            if not isinstance(occurrences, list) or not occurrences:
                report("empty_occurrences", seq, iso_acc)
                continue
            frame = mature if iso_acc == acc else protein[iso_acc]

            seen = Counter()
            for occ in occurrences:
                stats["occurrences"] += 1
                if not isinstance(occ, dict) or set(occ) != {
                    "start",
                    "stop",
                    "identity",
                }:
                    report(
                        "malformed_entry",
                        seq,
                        iso_acc,
                        json.dumps(occ),
                        "occurrence keys are not start, stop, identity",
                    )
                    continue
                o_start, o_stop, identity = (
                    number(occ["start"]),
                    number(occ["stop"]),
                    number(occ["identity"]),
                )
                if (
                    None in (o_start, o_stop, identity)
                    or o_start != int(o_start)
                    or o_stop != int(o_stop)
                    or not 0 < identity <= 100
                ):
                    report(
                        "malformed_entry",
                        seq,
                        iso_acc,
                        json.dumps(occ),
                        "non-numeric or out-of-range value",
                    )
                    continue
                o_start, o_stop = int(o_start), int(o_stop)
                recorded = f"{o_start}-{o_stop} ({identity:g}%)"
                seen[o_start, o_stop] += 1
                if seen[o_start, o_stop] == 2:
                    report(
                        "duplicate_occurrence",
                        seq,
                        iso_acc,
                        recorded,
                        "same occurrence listed more than once",
                    )
                    continue

                hits = positions(frame, seq)
                where = (
                    f"sequence found at {span(hits, len(seq))}"
                    if hits
                    else "sequence not found on the protein"
                )
                if not 1 <= o_start <= o_stop or o_stop > len(frame):
                    report(
                        "coordinates_out_of_bounds",
                        seq,
                        iso_acc,
                        recorded,
                        f"protein length {len(frame)}; {where}",
                    )
                    continue
                region = frame[o_start - 1 : o_stop]
                if identity == 100:
                    if region != seq:
                        report(
                            "coordinates_mismatch",
                            seq,
                            iso_acc,
                            recorded,
                            f"region reads {region}; {where}",
                        )
                    continue
                actual = aligned_identity(region, seq)
                if actual < identity - 1:
                    report(
                        "coordinates_mismatch",
                        seq,
                        iso_acc,
                        recorded,
                        f"region reads {region}; at most {actual:.1f}% identity",
                    )


for r in rows:
    audit(r[0], r[1], 1, *r[2:8])
    audit(r[0], r[1], 2, *r[8:14])

df = pl.DataFrame(
    records,
    schema={
        "stable_id": pl.String,
        "type": pl.String,
        "side": pl.Int64,
        "accession": pl.String,
        "name": pl.String,
        "start": pl.Int64,
        "stop": pl.Int64,
        "is_peptide": pl.Boolean,
        "problem": pl.String,
        "sequence": pl.String,
        "isoform": pl.String,
        "recorded": pl.String,
        "finding": pl.String,
    },
).sort("problem", "accession", "start", "stable_id", "side")
df.write_csv(OUTPUT, separator="\t")

print(f"{stats['entries']} mapping entries, {stats['occurrences']} occurrences audited")
print(
    df.group_by("problem").agg(
        pl.len(), pl.col("stable_id").n_unique().alias("descriptions")
    )
)
print(f"-> {OUTPUT}")

log_run(
    "scripts/mapping_corrections.py",
    OUTPUT,
    {
        "valid descriptions with mappings": len(rows),
        "mapping entries audited": stats["entries"],
        "occurrences audited": stats["occurrences"],
        "problems": len(df),
        "descriptions requiring correction": df["stable_id"].n_unique(),
        "peptide problems": int(df["is_peptide"].sum()),
    }
    | {f"problem {k}": v for k, v in sorted(Counter(df["problem"].to_list()).items())},
)
