"""Integrity problems on valid descriptions (hh and vh).

Audits what curators enter, looking for descriptions that are valid by the
usual filter yet malformed, so they can be sent back for correction. Three
families of check, reported side by side.

`description` -- the description as a whole:

    Descriptions absent from the `dataset` view. The view reaches the
    taxonomy through inner joins, so a description whose protein carries an
    NCBI taxon the `taxon` table does not have silently disappears from it,
    and from every dataset built on it. Found by reading
    `descriptions`/`associations` and taking what the view is missing.

    Live descriptions sharing a publication, a detection method and both
    protein regions. That is the grain of a description, so they are one
    description several times over -- typically proteins curated as separate
    entries that a later UniProt release merged into one. They need a
    curation pass to keep one and delete the rest.

`protein` -- the proteins as curated, `(accession, name, start, stop)`:

    A human protein is curated full length, so `startN` must be 1 and
    `stopN` the length of `accessionN` in the UniProt release the
    description was curated against, and `nameN` must read as UniProt's gene
    name does at that release. A mature viral protein is a region of its
    polyprotein, so `start2..stop2` must fit on it, under a curator-chosen
    `name2`: datasets key on `(accession2, name2)` rather than the
    coordinate triple, which holds only while one name means one region of
    one polyprotein, both ways round.

`mapping` -- every `mapping1` / `mapping2` entry, checked for malformation
and coordinate problems. A valid entry is an object `{sequence, isoforms}`
with a non-empty `sequence`, at least one isoform `{accession, occurrences}`,
and every isoform at least one occurrence `{start, stop, identity}`. The
isoform accession must be the protein's canonical accession or one of its
isoforms present in `proteins.sequences`, and never a non-canonical isoform
on a mature protein.

    Occurrence coordinates are checked against the protein as curated (the
    mature protein `startN..stopN` for the canonical accession, the whole
    isoform sequence otherwise). An occurrence recorded with identity 100
    must match the mapping sequence exactly; one below 100 must reach at
    least the recorded identity minus one point against the region (ungapped
    or gapped match count). `X` residues, numbers stored as strings, and
    isoforms or repeat positions left out by the curator are allowed and not
    reported.

`EXPLANATIONS` says in one line what each problem means; the markdown report
carries them next to the descriptions concerned.

One row per problem. Columns:

    stable_id, type          -- description and interactome
    check                    -- description, protein or mapping
    problem                  -- the check that failed
    side                     -- protein side concerned (1 or 2), when one is
    accession, name,
    start, stop              -- protein as curated (mature protein on the
                                viral side)
    is_peptide               -- mapping sequence is 4-20 aa inclusive
    sequence                 -- mapping sequence
    isoform                  -- isoform accession concerned
    recorded                 -- what the curator recorded
    finding                  -- what is actually there

Valid PPI filter: state = 'curated' AND is_obsolete1/2 IS FALSE AND
deleted_at IS NULL. The descriptions absent from the view are read from
`descriptions`/`associations` instead, where that filter cannot apply.

Run: uv run python scripts/dataset_integrity.py
Writes data/dataset_integrity.tsv and regenerates the markdown report.
"""

import json
import re
from collections import Counter
from datetime import UTC, datetime
from difflib import SequenceMatcher
from pathlib import Path
from typing import Literal, NamedTuple

import polars as pl

from drakkar.db import connect
from drakkar.runs import log_run

type Check = Literal["description", "protein", "mapping"]

type Problem = Literal[
    # description
    "not_in_dataset_view",
    "duplicated",
    # protein
    "not_full_length",
    "region_off_sequence",
    "sequence_missing",
    "name_not_uniprot",
    "name_at_other_region",
    "region_with_other_name",
    # mapping
    "no_isoform",
    "no_occurrence",
    "duplicated_occurrence",
    "repeated_isoform",
    "malformed_entry",
    "invalid_isoform",
    "position_mismatch",
    "position_off_sequence",
]

# What each problem means, for the curators reading the report.
EXPLANATIONS: dict[str, str] = {
    "not_in_dataset_view": (
        "The description is live and curated, but absent from the `dataset` "
        "view -- an NCBI taxon the `taxon` table does not have drops it from "
        "every dataset built on the view."
    ),
    "duplicated": (
        "Several live descriptions share a publication, a detection method "
        "and both protein regions, which is the grain of a description. "
        "Usually proteins curated as separate entries that a later UniProt "
        "release merged into one: keep one description, delete the rest."
    ),
    "not_full_length": (
        "A human protein is curated full length, so its region must start at "
        "1 and stop at the length of its canonical sequence in the UniProt "
        "release it was curated against."
    ),
    "region_off_sequence": (
        "The mature viral protein region does not fit on its polyprotein."
    ),
    "sequence_missing": (
        "No sequence is stored for the accession at that release, so nothing "
        "can be checked against it."
    ),
    "name_not_uniprot": (
        "The name of a human protein is its UniProt gene name, so it must "
        "read as UniProt has it at that release -- a UniProt rename the "
        "description was never brought in line with lands here."
    ),
    "name_at_other_region": (
        "The same mature-protein name is used at other coordinates of the "
        "same polyprotein, so the name no longer identifies one protein."
    ),
    "region_with_other_name": (
        "The same polyprotein region is curated under another name, so the "
        "mature protein has two identities."
    ),
    "no_isoform": "The mapping entry records no isoform at all.",
    "no_occurrence": "The isoform records no occurrence at all.",
    "duplicated_occurrence": "The same occurrence is listed more than once.",
    "repeated_isoform": "The same accession opens several isoform blocks.",
    "malformed_entry": ("Keys or values that are not what a mapping entry calls for."),
    "invalid_isoform": (
        "An accession that is not an isoform of the protein, or a "
        "non-canonical isoform on a mature protein, which the curation "
        "interface does not allow."
    ),
    "position_mismatch": (
        "The recorded position does not read back as the curated sequence, "
        "to the recorded identity."
    ),
    "position_off_sequence": (
        "The recorded position does not fit the protein it was curated on."
    ),
}

TSV = "data/dataset_integrity.tsv"
MD = "reports/dataset_integrity.md"

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
"""


class Finding(NamedTuple):
    """One problem found on one description."""

    stable_id: str
    type: str
    check: Check
    problem: Problem
    side: int | None = None
    accession: str | None = None
    name: str | None = None
    start: int | None = None
    stop: int | None = None
    is_peptide: bool | None = None
    sequence: str | None = None
    isoform: str | None = None
    recorded: str = ""
    finding: str = ""


findings: list[Finding] = []
stats: Counter[str] = Counter()


# -- descriptions the view drops --------------------------------------------

# `dataset` reaches the taxonomy through inner joins, so a description whose
# protein carries a taxon absent from `taxon` never appears in it. Read the
# base tables and take what the view is missing.
MISSING = """
    SELECT d.stable_id,
           CASE WHEN p2.type = 'v' THEN 'vh' ELSE 'hh' END,
           p1.accession, d.name1, p1.ncbi_taxon_id, t1.taxon_id IS NULL,
           p2.accession, d.name2, p2.ncbi_taxon_id, t2.taxon_id IS NULL
    FROM descriptions d
    JOIN associations a ON a.id = d.association_id
    JOIN proteins p1 ON p1.id = d.protein1_id
    JOIN proteins p2 ON p2.id = d.protein2_id
    LEFT JOIN taxon t1 ON t1.ncbi_taxon_id = p1.ncbi_taxon_id
    LEFT JOIN taxon t2 ON t2.ncbi_taxon_id = p2.ncbi_taxon_id
    WHERE d.deleted_at IS NULL
      AND a.state = 'curated'
      AND NOT EXISTS (
          SELECT 1 FROM dataset v
          WHERE v.stable_id = d.stable_id AND v.deleted_at IS NULL
      )
    ORDER BY d.stable_id
"""


def missing_from_dataset(
    stable_id: str,
    type_: str,
    accession1: str,
    name1: str,
    taxon1: int,
    taxon1_unknown: bool,
    accession2: str,
    name2: str,
    taxon2: int,
    taxon2_unknown: bool,
) -> Finding:
    """One live curated description the dataset view does not carry."""
    unknown = (
        (1, accession1, name1, taxon1)
        if taxon1_unknown
        else (2, accession2, name2, taxon2)
        if taxon2_unknown
        else None
    )
    if unknown is None:
        return Finding(
            stable_id=stable_id.strip(),
            type=type_,
            check="description",
            problem="not_in_dataset_view",
            recorded=f"{accession1} x {accession2}",
            finding="live curated description absent from the dataset view",
        )
    side, accession, name, taxon = unknown
    return Finding(
        stable_id=stable_id.strip(),
        type=type_,
        check="description",
        problem="not_in_dataset_view",
        side=side,
        accession=accession,
        name=name,
        recorded=f"{accession1} x {accession2}",
        finding=(
            f"NCBI taxon {taxon} is absent from taxon, so the dataset view "
            f"drops the description"
        ),
    )


# -- descriptions curated twice ---------------------------------------------

# A description is (publication x method x region 1 x region 2), so two live
# ones sharing all four are the same description twice over.
DUPLICATES = f"""
    SELECT stable_id, type, pmid, psimi_id,
           accession1, accession2, name2, start2, stop2,
           copies, mappings_differ
    FROM (
        SELECT stable_id, type, pmid, psimi_id,
               accession1, accession2, name2, start2, stop2,
               count(*) OVER grain AS copies,
               min(mapping1::text || mapping2::text) OVER grain
                   <> max(mapping1::text || mapping2::text) OVER grain
                   AS mappings_differ
        FROM dataset
        WHERE {VALID}
        WINDOW grain AS (
            PARTITION BY pmid, psimi_id,
                         accession1, start1, stop1,
                         accession2, start2, stop2
        )
    ) grains
    WHERE copies > 1
    ORDER BY pmid, accession1, accession2, start2, stable_id
"""


def duplicate_description(
    stable_id: str,
    type_: str,
    pmid: int,
    psimi_id: str,
    accession1: str,
    accession2: str,
    name2: str,
    start2: int,
    stop2: int,
    copies: int,
    mappings_differ: bool,
) -> Finding:
    """One of several live descriptions sharing a grain."""
    differ = "their mappings differ" if mappings_differ else "mappings included"
    return Finding(
        stable_id=stable_id.strip(),
        type=type_,
        check="description",
        problem="duplicated",
        side=2,
        accession=accession2,
        name=name2,
        start=start2,
        stop=stop2,
        recorded=f"pmid {pmid}, {psimi_id}, {accession1} x {accession2}",
        finding=f"{copies} live descriptions of the same thing, {differ}",
    )


# -- proteins as curated ----------------------------------------------------

# A human protein spans its whole canonical sequence and carries its UniProt
# gene name; a mature viral protein is a region of its polyprotein, under a
# name the curator chose.
PROTEINS = f"""
    SELECT stable_id, type, side, accession, name, uniprot_name,
           region_start, region_stop, seq_length, version
    FROM (
        SELECT d.stable_id, d.type, 1 AS side,
               d.accession1 AS accession, d.name1 AS name,
               p.name AS uniprot_name,
               d.start1 AS region_start, d.stop1 AS region_stop,
               length(p.sequences ->> d.accession1) AS seq_length,
               p.version AS version
        FROM dataset d
        JOIN proteins p ON p.id = d.protein1_id
        WHERE {VALID}
        UNION ALL
        SELECT d.stable_id, d.type, 2,
               d.accession2, d.name2, p.name,
               d.start2, d.stop2,
               length(p.sequences ->> d.accession2),
               p.version
        FROM dataset d
        JOIN proteins p ON p.id = d.protein2_id
        WHERE {VALID}
    ) curated
    WHERE seq_length IS NULL
       OR region_start < 1
       OR region_stop < region_start
       OR region_stop > seq_length
       OR ((side = 1 OR type = 'hh')
           AND (region_start <> 1
                OR region_stop <> seq_length
                OR name <> uniprot_name))
    ORDER BY accession, region_start, stable_id, side
"""


def protein_problems(
    stable_id: str,
    type_: str,
    side: int,
    accession: str,
    name: str,
    uniprot_name: str,
    start: int,
    stop: int,
    seq_length: int | None,
    version: str,
) -> list[Finding]:
    """The region and naming problems of one protein as curated."""

    def problem(problem: Problem, recorded: str, what: str) -> Finding:
        return Finding(
            stable_id=stable_id.strip(),
            type=type_,
            check="protein",
            problem=problem,
            side=side,
            accession=accession,
            name=name,
            start=start,
            stop=stop,
            recorded=recorded,
            finding=what,
        )

    if seq_length is None:
        return [
            problem(
                "sequence_missing",
                f"{start}-{stop}",
                f"{accession} has no sequence at release {version}",
            )
        ]
    human = side == 1 or type_ == "hh"
    length = f"{accession} is {seq_length} aa at release {version}"
    reported: list[Finding] = []
    if human and (start != 1 or stop != seq_length):
        reported.append(problem("not_full_length", f"{start}-{stop}", length))
    elif not human and not 1 <= start <= stop <= seq_length:
        reported.append(problem("region_off_sequence", f"{start}-{stop}", length))
    if human and name != uniprot_name:
        reported.append(
            problem(
                "name_not_uniprot",
                name,
                f"{accession} is {uniprot_name} at release {version}",
            )
        )
    return reported


# -- mature viral protein names ---------------------------------------------

# A mature protein is (accession2, start2, stop2), and `name2` is expected to
# name exactly that -- datasets key on the pair rather than the triple.
VIRAL_NAMES = f"""
    SELECT stable_id, type, accession2, name2, start2, stop2,
           name_spans, span_names
    FROM (
        SELECT stable_id, type, accession2, name2, start2, stop2,
               min(start2) OVER named <> max(start2) OVER named
                   OR min(stop2) OVER named <> max(stop2) OVER named
                   AS name_spans,
               min(name2) OVER span <> max(name2) OVER span AS span_names
        FROM dataset
        WHERE {VALID} AND type = 'vh'
        WINDOW named AS (PARTITION BY accession2, name2),
               span AS (PARTITION BY accession2, start2, stop2)
    ) mature
    WHERE name_spans OR span_names
    ORDER BY accession2, start2, name2, stable_id
"""


def viral_name(
    stable_id: str,
    type_: str,
    accession: str,
    name: str,
    start: int,
    stop: int,
    name_spans: bool,
    span_names: bool,
) -> list[Finding]:
    """The naming problems of one mature protein, one or both of them."""

    def problem(problem: Problem, what: str) -> Finding:
        return Finding(
            stable_id=stable_id.strip(),
            type=type_,
            check="protein",
            problem=problem,
            side=2,
            accession=accession,
            name=name,
            start=start,
            stop=stop,
            recorded=f"{name} {start}-{stop}",
            finding=what,
        )

    reported: list[Finding] = []
    if name_spans:
        reported.append(
            problem(
                "name_at_other_region",
                f"{name} also names another region of {accession}",
            )
        )
    if span_names:
        reported.append(
            problem(
                "region_with_other_name",
                f"{start}-{stop} of {accession} is also named otherwise",
            )
        )
    return reported


# -- mappings ---------------------------------------------------------------

MAPPINGS = f"""
    SELECT stable_id, type,
           protein1_id, accession1, name1, start1, stop1, mapping1::text,
           protein2_id, accession2, name2, start2, stop2, mapping2::text
    FROM dataset
    WHERE {VALID}
      AND (json_array_length(mapping1) > 0 OR json_array_length(mapping2) > 0)
"""


def _number(value: object) -> float:
    """Numeric value of an occurrence field, stored as a string about a third
    of the time. Raises when it is not a number -- `float` does the raising
    for a string that does not read as one."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise TypeError(f"not a number: {value!r}")
    return float(value)


def coordinate(value: object) -> int:
    """A recorded coordinate. Raises ValueError when it is not a whole number."""
    number = _number(value)
    if number != int(number):
        raise ValueError(f"not a whole number: {value!r}")
    return int(number)


def percentage(value: object) -> float:
    """A recorded identity. Raises ValueError when it is out of range."""
    number = _number(value)
    if not 0 < number <= 100:
        raise ValueError(f"not a percentage: {value!r}")
    return number


def aligned_identity(read: str, seq: str) -> float:
    """Identity (%) of `seq` against `read`, matches over the longer length.

    Residue-by-residue for equal lengths; when that falls short (or lengths
    differ) a gapped match count from difflib, which never overestimates.
    """
    longest = max(len(read), len(seq))
    ungapped = 0
    if len(read) == len(seq):
        ungapped = sum(a == b for a, b in zip(read, seq, strict=True))
    matcher = SequenceMatcher(None, read, seq, autojunk=False)
    gapped = sum(block.size for block in matcher.get_matching_blocks())
    return 100 * max(ungapped, gapped) / longest


def positions(frame: str, seq: str) -> list[int]:
    return [i + 1 for i in range(len(frame)) if frame.startswith(seq, i)]


def span(starts: list[int], length: int) -> str:
    return ", ".join(f"{s}-{s + length - 1}" for s in starts)


def audit(
    sequences: dict[int, dict[str, str]],
    stable_id: str,
    type_: str,
    side: int,
    protein_id: int,
    acc: str,
    name: str,
    start: int,
    stop: int,
    raw: str,
) -> None:
    """Report every problem in one side's mapping."""

    def report(
        problem: Problem,
        seq: object,
        isoform: str | None = "",
        recorded: str = "",
        finding: str = "",
    ) -> None:
        findings.append(
            Finding(
                stable_id=stable_id.strip(),
                type=type_,
                check="mapping",
                problem=problem,
                side=side,
                accession=acc,
                name=name,
                start=start,
                stop=stop,
                is_peptide=isinstance(seq, str) and 4 <= len(seq) <= 20,
                sequence=seq if isinstance(seq, str) else json.dumps(seq),
                isoform=isoform,
                recorded=recorded,
                finding=finding,
            )
        )

    mapping = json.loads(raw)
    if not isinstance(mapping, list):
        report("malformed_entry", None, finding="mapping is not an array")
        return
    protein = sequences[protein_id]
    if acc not in protein:
        return  # reported by the region check, and nothing to audit against
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
            report("no_isoform", seq, finding=finding)
            continue

        blocks = Counter(
            iso.get("accession") for iso in isoforms if isinstance(iso, dict)
        )
        for iso_acc, n in blocks.items():
            if n > 1:
                report(
                    "repeated_isoform",
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
                report("no_occurrence", seq, iso_acc)
                continue
            frame = mature if iso_acc == acc else protein[iso_acc]

            seen: Counter[tuple[int, int]] = Counter()
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
                try:
                    o_start = coordinate(occ["start"])
                    o_stop = coordinate(occ["stop"])
                    identity = percentage(occ["identity"])
                except TypeError, ValueError:
                    report(
                        "malformed_entry",
                        seq,
                        iso_acc,
                        json.dumps(occ),
                        "non-numeric or out-of-range value",
                    )
                    continue
                recorded = f"{o_start}-{o_stop} ({identity:g}%)"
                seen[o_start, o_stop] += 1
                if seen[o_start, o_stop] == 2:
                    report(
                        "duplicated_occurrence",
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
                        "position_off_sequence",
                        seq,
                        iso_acc,
                        recorded,
                        f"protein length {len(frame)}; {where}",
                    )
                    continue
                span_read = frame[o_start - 1 : o_stop]
                if identity == 100:
                    if span_read != seq:
                        report(
                            "position_mismatch",
                            seq,
                            iso_acc,
                            recorded,
                            f"region reads {span_read}; {where}",
                        )
                    continue
                actual = aligned_identity(span_read, seq)
                if actual < identity - 1:
                    report(
                        "position_mismatch",
                        seq,
                        iso_acc,
                        recorded,
                        f"region reads {span_read}; at most {actual:.1f}% identity",
                    )


# -- run every check --------------------------------------------------------

# Every valid description goes through the description and protein checks;
# the mapping audit only sees those carrying a mapping, the common case
# being none at all.
SCOPE = f"SELECT count(*) FROM dataset WHERE {VALID}"

with connect() as conn, conn.cursor() as cur:
    cur.execute(SCOPE)
    audited = cur.fetchone()
    if audited is None:
        raise ValueError("counting the valid descriptions returned nothing")
    stats["valid descriptions"] = audited[0]

    cur.execute(MISSING)
    findings += [missing_from_dataset(*row) for row in cur.fetchall()]

    cur.execute(DUPLICATES)
    findings += [duplicate_description(*row) for row in cur.fetchall()]

    cur.execute(PROTEINS)
    findings += [f for row in cur.fetchall() for f in protein_problems(*row)]

    cur.execute(VIRAL_NAMES)
    findings += [f for row in cur.fetchall() for f in viral_name(*row)]

    cur.execute(MAPPINGS)
    mappings = cur.fetchall()
    ids = sorted({r[2] for r in mappings} | {r[8] for r in mappings})
    cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
    sequences = dict(cur.fetchall())

stats["descriptions with a mapping"] = len(mappings)
for r in mappings:
    audit(sequences, r[0], r[1], 1, *r[2:8])
    audit(sequences, r[0], r[1], 2, *r[8:14])

df = pl.DataFrame(
    findings,
    schema={
        "stable_id": pl.String,
        "type": pl.String,
        "check": pl.String,
        "problem": pl.String,
        "side": pl.Int64,
        "accession": pl.String,
        "name": pl.String,
        "start": pl.Int64,
        "stop": pl.Int64,
        "is_peptide": pl.Boolean,
        "sequence": pl.String,
        "isoform": pl.String,
        "recorded": pl.String,
        "finding": pl.String,
    },
    orient="row",
).sort("check", "problem", "accession", "start", "stable_id", "side", "sequence")
df.write_csv(TSV, separator="\t")


# -- markdown report --------------------------------------------------------

# The export is one row per problem, for triage in a spreadsheet; the report
# is the same findings grouped, so a curation pass reads one line per thing
# to fix rather than one per description it shows up in.
GROUPED = [
    "side",
    "accession",
    "name",
    "start",
    "stop",
    "sequence",
    "isoform",
    "recorded",
    "finding",
]

LISTED = 8  # stable_ids spelled out per row before the rest are counted

CELL = 60  # curated text is cut here -- a whole sequence widens a table out


def trimmed(rows: pl.DataFrame) -> pl.DataFrame:
    """Cut the curated text columns to `CELL`."""
    return rows.with_columns(
        pl.when(pl.col(c).str.len_chars() > CELL)
        .then(pl.col(c).str.slice(0, CELL - 3) + "...")
        .otherwise(pl.col(c))
        for c in rows.columns
        if c in GROUPED and rows[c].dtype == pl.String
    )


def markdown_table(rows: pl.DataFrame) -> str:
    """Render a frame as a markdown table, numbers right-aligned."""
    align = ["---:" if rows[c].dtype.is_numeric() else "---" for c in rows.columns]
    header = f"| {' | '.join(rows.columns)} |\n| {' | '.join(align)} |\n"
    body = "".join(
        "| " + " | ".join("" if v is None else str(v) for v in row) + " |\n"
        for row in rows.iter_rows()
    )
    return header + body


def count(n: int, thing: str) -> str:
    """`n` of something, pluralised."""
    return f"{n:,} {thing}" if n == 1 else f"{n:,} {thing}s"


def section(check: str, problem: str, rows: pl.DataFrame) -> str:
    """One problem: what it means, and every finding of it, grouped."""
    grouped = (
        rows.group_by(GROUPED, maintain_order=True)
        .agg(pl.col("stable_id"))
        .with_columns(descriptions=pl.col("stable_id").list.len())
        .with_columns(
            stable_ids=pl.when(pl.col("descriptions") > LISTED)
            .then(
                pl.col("stable_id").list.head(LISTED).list.join(", ")
                + pl.format(" (+{} more)", pl.col("descriptions") - LISTED)
            )
            .otherwise(pl.col("stable_id").list.join(", "))
        )
        .drop("stable_id")
    )
    # A column no finding of this problem fills says nothing here.
    grouped = grouped.drop(
        c for c in GROUPED if grouped[c].null_count() == len(grouped)
    )
    descriptions = rows["stable_id"].n_unique()
    return (
        f"## {check} -- {problem}\n\n"
        f"{EXPLANATIONS[problem]}\n\n"
        f"{count(len(rows), 'problem')} over {count(descriptions, 'description')}.\n\n"
        f"{markdown_table(trimmed(grouped))}\n"
    )


def report(df: pl.DataFrame) -> str:
    """The whole audit as markdown, a section per problem found."""
    summary = (
        df.group_by("check", "problem")
        .agg(problems=pl.len(), descriptions=pl.col("stable_id").n_unique())
        .sort("check", "problem")
    )
    head = (
        "# Integrity problems on valid descriptions\n\n"
        f"**Date:** {datetime.now(UTC).date().isoformat()}  \n"
        "**Source:** drakkar `dataset` MV, valid PPI filter, `hh` and `vh`  \n"
        "**Script:** `scripts/dataset_integrity.py`  \n"
        f"**Export:** `{TSV}`\n\n"
        "Descriptions that pass the valid PPI filter yet are malformed, for a "
        "curation pass. Each section says what the problem means and lists "
        "every finding of it, one line per thing to fix. Long cells are cut "
        "here; the export carries them whole.\n\n"
        f"Audited {stats['valid descriptions']:,} valid descriptions, "
        f"{stats['descriptions with a mapping']:,} of them carrying a "
        "mapping, plus the live curated descriptions the view leaves out.\n\n"
        "## Summary\n\n"
        f"{markdown_table(summary)}\n"
    )
    sections = [
        section(check, problem, df.filter(pl.col("problem") == problem))
        for check, problem in summary.select("check", "problem").iter_rows()
    ]
    return head + "".join(sections)


Path(MD).parent.mkdir(parents=True, exist_ok=True)
Path(MD).write_text(report(df))

print(
    f"{stats['valid descriptions']:,} valid descriptions audited, "
    f"{stats['descriptions with a mapping']:,} of them carrying a mapping: "
    f"{stats['entries']:,} entries, {stats['occurrences']:,} occurrences"
)
print(df.group_by("check", "problem").agg(pl.len()).sort("check", "problem"))
print(f"-> {TSV} + {MD}")

log_run(
    "scripts/dataset_integrity.py",
    f"{TSV} + {MD}",
    {
        "valid descriptions audited": stats["valid descriptions"],
        "of them carrying a mapping": stats["descriptions with a mapping"],
        "mapping entries audited": stats["entries"],
        "mapping occurrences audited": stats["occurrences"],
        "problems": len(df),
        "descriptions requiring correction": df["stable_id"].n_unique(),
        "mapping problems on peptides": int(df["is_peptide"].sum()),
    }
    # A label no longer says which family it belongs to, so the log does.
    | {
        f"problem {check}/{problem}": n
        for (check, problem), n in sorted(
            Counter(zip(df["check"], df["problem"], strict=True)).items()
        )
    },
)
