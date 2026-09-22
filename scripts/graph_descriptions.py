"""All valid PPI descriptions with their proteins, publications and peptides.

Exported as three TSV files into the folder given on the command line, a
graph of the curated interactome: descriptions are the edges, publications
and peptides hang off them by `stable_id` / `pmid`.

`descriptions.tsv` -- one row per valid description, i.e. the MV grain
*(publication x detection method x protein 1 region x protein 2 region)*.
Both interactomes are included; `type` tells them apart (`hh` / `vh`).

  stable_id, type, pmid,
  psimi_id, method,
  accession1, start1, stop1, name1, description1,
  accession2, start2, stop2, name2, description2

Each protein is identified by its `(accession, start, stop)` triple -- never
by accession alone, since one viral polyprotein accession yields many mature
proteins at different coordinates. Human proteins (side 1) are always full
length, so their triple is always `(accession1, 1, length)`. Sequences are
deliberately not exported; recover them from `proteins.sequences`.

`publications.tsv` -- one row per distinct pmid referenced above, so the
article metadata is not repeated on every description:

  pmid, title, year, journal, authors, abstract

Parsed out of the PubMed record stored in `publications.metadata`: `authors`
is `"LastName Initials"` per author (collective names kept as is), `"; "`
-separated in author order; `year` falls back to the leading year of
`MedlineDate` when `PubDate` carries no `Year`; a structured abstract's parts
are joined with a space; missing fields come out empty.

`peptides.tsv` -- the short binding regions curators reported, both sides
together, one row per distinct peptide of a description's protein:

  stable_id, sequence, source_type, source_accession, source_start, source_stop

A peptide is a `mapping1` / `mapping2` entry whose sequence is 5-20 aa
inclusive (a one-time exception to the 4-20 team convention). Entries are
read through `drakkar.mappings.occurrences`, which places each one on the
sequence it was curated on, so an entry that fits nowhere on its protein
yields no peptide at all. The position it was placed at is *not* exported:
the peptide is attached to the protein it belongs to, by that protein's
identifier triple, which is the description's `(accessionN, startN, stopN)`
and joins back to `descriptions.tsv`. `source_type` is that protein's
`typeN` -- `h` for a human protein, `v` for a viral one, so side 1 is always
`h`. Occurrences that differ only by position therefore collapse to one
row.

Valid PPI filter applied (curated, both accessions live, current revision).

Run: uv run python scripts/graph_descriptions.py data/graph
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import polars as pl

from drakkar.db import connect, fetch, stream
from drakkar.mappings import occurrences
from drakkar.runs import log_run

PEPTIDE_MIN, PEPTIDE_MAX = 5, 20

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
"""

DESCRIPTIONS_QUERY = f"""
SELECT
    stable_id,
    type,
    pmid,
    psimi_id,
    method,
    accession1,
    start1,
    stop1,
    name1,
    description1,
    accession2,
    start2,
    stop2,
    name2,
    description2
FROM dataset
WHERE {VALID}
ORDER BY pmid, accession1, accession2, start2, stop2, psimi_id
"""

# Only the Article subtree is pulled back: the full PubMed record carries
# MeSH, chemicals and history we do not export.
#
# GROUP BY, not SELECT DISTINCT: pmid is the publications primary key, so
# grouping on it alone is enough to select the metadata, and the subtree is
# then extracted once per publication. Under DISTINCT the extracted ~5 kB of
# JSON is part of the dedup key instead, so all 400k joined rows carry it
# through the sort.
PUBLICATIONS_QUERY = f"""
SELECT
    p.pmid,
    (p.metadata -> 'PubmedArticle' -> 'MedlineCitation' -> 'Article')::text
        AS article
FROM dataset d
JOIN publications p ON p.pmid = d.pmid
WHERE {VALID}
GROUP BY p.pmid
ORDER BY p.pmid
"""

# Descriptions reporting a binding region on either side -- a few percent of
# the view, so the peptide pass never holds the whole thing.
MAPPINGS_QUERY = f"""
SELECT
    stable_id,
    protein1_id, type1, accession1, start1, stop1, mapping1,
    protein2_id, type2, accession2, start2, stop2, mapping2
FROM dataset
WHERE {VALID}
  AND (json_array_length(mapping1) > 0 OR json_array_length(mapping2) > 0)
ORDER BY stable_id
"""

YEAR_RE = re.compile(r"\d{4}")


def _text(node: Any) -> str:
    """Flatten a PubMed JSON text node to a plain string.

    A node is a string, a list of nodes (structured abstract, inline markup)
    or an object whose non-attribute values hold the text.
    """
    if node is None:
        return ""
    if isinstance(node, str):
        return node.strip()
    if isinstance(node, list):
        return " ".join(filter(None, (_text(item) for item in node)))
    if isinstance(node, dict):
        return " ".join(
            filter(
                None,
                (_text(v) for k, v in node.items() if k != "@attributes"),
            )
        )
    return str(node)


def _author(author: Any) -> str:
    if not isinstance(author, dict):
        return _text(author)
    collective = _text(author.get("CollectiveName"))
    if collective:
        return collective
    initials = _text(author.get("Initials")) or _text(author.get("ForeName"))
    parts = (_text(author.get("LastName")), initials, _text(author.get("Suffix")))
    return " ".join(part for part in parts if part)


def _authors(article: dict) -> str:
    authors = (article.get("AuthorList") or {}).get("Author")
    if authors is None:
        return ""
    if not isinstance(authors, list):
        authors = [authors]
    return "; ".join(filter(None, (_author(a) for a in authors)))


def _year(article: dict) -> str:
    journal = article.get("Journal") or {}
    pubdate = (journal.get("JournalIssue") or {}).get("PubDate") or {}
    year = _text(pubdate.get("Year"))
    if year:
        return year
    match = YEAR_RE.search(_text(pubdate.get("MedlineDate")))
    return match.group(0) if match else ""


def _publication(pmid: int, article: dict) -> dict[str, Any]:
    journal = article.get("Journal") or {}
    return {
        "pmid": pmid,
        "title": _text(article.get("ArticleTitle")),
        "year": _year(article),
        "journal": _text(journal.get("Title")),
        "authors": _authors(article),
        "abstract": _text((article.get("Abstract") or {}).get("AbstractText")),
    }


def publications() -> pl.DataFrame:
    rows = fetch(PUBLICATIONS_QUERY)
    missing = rows.filter(pl.col("article").is_null())["pmid"].to_list()
    if missing:
        # Book chapters and unpopulated records have no Article subtree; fail
        # loudly rather than exporting rows with every field blank.
        raise RuntimeError(
            f"{len(missing)} pmid(s) without a PubmedArticle record, e.g. {missing[:5]}"
        )
    return pl.DataFrame(
        [_publication(pmid, json.loads(article)) for pmid, article in rows.iter_rows()]
    )


def _sources(
    protein: dict[str, str], accession: str, start: int, stop: int
) -> dict[str, str]:
    """The sequences a mapping on this protein can have been curated on.

    A protein curated full length -- every human protein, and a viral protein
    that is not a mature one -- brings its canonical sequence and its
    isoforms, since a peptide may have been located on an isoform. A mature
    viral protein is a region of a polyprotein and has no isoforms of its
    own, so its mature region is the only sequence to match against.
    """
    canonical = protein[accession]
    if start == 1 and stop == len(canonical):
        return protein
    return {accession: canonical[start - 1 : stop]}


def peptides() -> pl.DataFrame:
    with connect() as conn, conn.cursor() as cur:
        cur.execute(MAPPINGS_QUERY)
        rows = cur.fetchall()
        ids = sorted({row[1] for row in rows} | {row[7] for row in rows})
        cur.execute("SELECT id, sequences FROM proteins WHERE id = ANY(%s)", (ids,))
        sequences = dict(cur.fetchall())

    records = []
    for stable_id, *columns in rows:
        side1, side2 = columns[:6], columns[6:]
        for protein_id, type_, accession, start, stop, mapping in (side1, side2):
            sources = _sources(sequences[protein_id], accession, start, stop)
            for occurrence in occurrences(mapping, sources, PEPTIDE_MIN, PEPTIDE_MAX):
                records.append(
                    {
                        "stable_id": stable_id,
                        "sequence": occurrence.sequence,
                        "source_type": type_,
                        "source_accession": accession,
                        "source_start": start,
                        "source_stop": stop,
                    }
                )
    # A peptide occurring twice on the same protein is one peptide of it once
    # its position is dropped.
    return pl.DataFrame(records).unique(maintain_order=True)


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    descriptions_tsv = out / "descriptions.tsv"
    publications_tsv = out / "publications.tsv"
    peptides_tsv = out / "peptides.tsv"

    # 400k rows: streamed to disk rather than materialised as a DataFrame.
    stream(DESCRIPTIONS_QUERY, descriptions_tsv)

    pubs = publications()
    pubs.write_csv(publications_tsv, separator="\t")

    peps = peptides()
    peps.write_csv(peptides_tsv, separator="\t")

    log_run(
        "scripts/graph_descriptions.py",
        f"{descriptions_tsv}, {publications_tsv}, {peptides_tsv}",
        {"publications": pubs.height, "peptides": peps.height},
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export the drakkar PPI graph.")
    parser.add_argument("out", type=Path, help="folder the three TSV files go to")
    main(parser.parse_args().out)
