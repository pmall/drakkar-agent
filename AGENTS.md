# drakkar-agent

Assistant for querying the local **drakkar** protein–protein interaction database and
producing datasets. Curated human–human (`hh`) and human–viral (`vh`) PPIs.

## Setup

Credentials live in `.env` (gitignored; see `.env.example`). Never print the password.

Layout:

```
drakkar/db.py  connection + query helpers (loads .env via python-dotenv)
scripts/       one dataset script per file
data/          exports (gitignored)
```

Python env is managed with **uv**:

```python
from drakkar.db import fetch, export

df = fetch("select ... from dataset where ...")  # -> polars DataFrame
export("select ...", "data/my_dataset.tsv")  # .csv / .tsv / .parquet
```

`drakkar` is installed into the venv as a package, so scripts run from anywhere:

```bash
uv run python scripts/my_dataset.py
```

`scripts/example_dataset.py` is a working template. **All exports go to `data/`.**
`psql` is fine for quick counts and exploration.

`data/` is gitignored, so a script is the only record of how its dataset was built.
Treat `scripts/` as an append-only provenance log: commit every dataset script, one
commit per dataset, and record the run date and the counts it produced in the
docstring. Nothing imports them, so they are never edited to keep them running —
a new dataset is a new script.

Before considering a script finished, run it, then format and lint the repo:

```bash
uv run ruff format .
uv run ruff check --fix .
```

## Data model

Curation source tables: `runs → associations → descriptions`, plus `methods`,
`proteins` (+ `proteins_versions`), `publications`, `taxon`/`taxon_name`, `keywords`.

Everything needed for datasets is denormalized in the **`dataset` materialized view**
(423,652 rows). Query that; drop to the base tables only for things it doesn't carry
(sequences, publication metadata, UniProt features).

Schema varies between database versions — verify a table exists before depending on it
rather than assuming this file is exhaustive.

### Grain

One row per **description** = *(publication × detection method × protein 1 region ×
protein 2 region)*. `stable_id` identifies a description across its revisions;
`deleted_at IS NULL` selects the current revision.

A **run** is a batch of publications the curation team processed together
(`run`/`run_id`, 8 `hh` runs and 60 `vh` runs).

### Valid PPI description — always apply

```sql
WHERE state = 'curated'
  AND is_obsolete1 IS FALSE
  AND is_obsolete2 IS FALSE
  AND deleted_at IS NULL
```

That is **401,705 rows**: 285,690 `hh` + 116,015 `vh`, over 15,888 human proteins,
14,798 publications, 209 detection methods. Use this filter for every dataset unless
explicitly told otherwise; if a request needs it relaxed, say so.

- `state`: `curated` / `selected` / `discarded` (`pending` never reaches the view).
- `is_obsoleteN`: the accession no longer exists in current UniProt
  (`original_versionN` = release at curation time, `current_versionN` = current release,
  NULL when obsolete).

### Sides

Protein 1 is **always human** (`type1 = 'h'`, taxon 9606). Protein 2 is human (`hh`)
or viral (`vh`) — check `type` / `type2`.

| | side 1 (human) | side 2 |
|---|---|---|
| accession | `accession1` | `accession2` |
| name | `name1` — UniProt gene name | `name2` — **curator-chosen** |
| region | `start1`/`stop1` (always full length, `start1 = 1`) | `start2`/`stop2` |
| taxon | 9606, `taxon1 = 'Homo sapiens'` | `ncbi_taxon_id2`, `taxon2`, `left_value2`/`right_value2` |

### Mature viral proteins — important

Viral proteins are **mature proteins**: subsequences of a polyprotein. A viral protein's
identity is therefore the triple

```
(accession2, start2, stop2)   -- polyprotein accession + coordinates on it
```

with `name2` the name assigned by the curation team. **Never key viral proteins on
`accession2` alone** — one polyprotein yields many mature proteins (e.g. `P0DTD1`
SARS-CoV-2 ORF1ab → 15 mature proteins). Across valid rows: 3,483 distinct mature
proteins over 2,984 accessions. `name2` is stable per triple and not reused across
different coordinates, so `(accession2, name2)` also identifies a mature protein
in practice — but prefer the coordinate triple.

To recover a mature protein sequence, substring the polyprotein from
`proteins.sequences`:

```sql
substring(p.sequences->>d.accession2 FROM d.start2 FOR d.stop2 - d.start2 + 1)
```

Human proteins are full-length, so `(accession1)` identifies them.

### Mappings — binding regions

`mapping1` / `mapping2` (JSON) are the subsequences curators reported as **sufficient
to bind the partner**. Empty `[]` when the publication reported no such region: among
valid rows only 9,512 have `mapping1`, 13,423 have `mapping2`, 4,514 have both.

```json
[{"sequence": "<peptide>",
  "isoforms": [{"accession": "<isoform acc>",
                "occurrences": [{"start": 1203, "stop": 1264, "identity": 100}]}]}]
```

- `sequence` is the curated peptide as reported.
- `isoforms[].accession` keys into `proteins.sequences`, a JSON object holding the
  canonical accession plus `-2`, `-3`, … isoforms.
  - **Human side:** frequently a non-canonical isoform (`P12345-2`) — the peptide was
    located on that isoform rather than the canonical sequence.
  - **Viral side:** for a *mature* protein (`start2 > 1`) it is **always** the canonical
    accession. **The curation interface does not allow a mature protein to be defined on
    an isoform** — mature proteins are always coordinates on the canonical sequence. Zero
    exceptions in the whole table. Non-canonical viral isoform references therefore occur
    only for full-length viral proteins (`start2 = 1`): 12 cases (Tat, E1A, ICP22, LT,
    HBZ, UL37, VSV M, HEV ORF2).
  - A polyprotein entry may still *have* isoforms in UniProt (12 of the 323 accessions
    used with `start2 > 1` do — FMDV, PRRSV, FIPV frameshift/alternative-ORF products);
    they are simply never used as a mature-protein frame.
- `identity` is a percentage; it is **not always 100** — the peptide is aligned, not
  necessarily exact. Filter on it when exactness matters.
- **Occurrence coordinates are 1-based relative to the protein as curated, i.e. to
  `startN`, not to the polyprotein.** For mature viral proteins the absolute position
  on the polyprotein is `start2 + occurrence.start - 1` (verified against
  `proteins.sequences`). On the human side `start1 = 1`, so the two frames coincide.

#### Peptides

By team convention a **peptide** is a mapping `sequence` of **4–20 aa inclusive** — a
short binding region, as opposed to a domain-scale mapping. Nothing in the schema flags
it; derive it by length:

```sql
WHERE length(e->>'sequence') BETWEEN 4 AND 20
```

Across all mappings that is 4,604 occurrences / 4,386 distinct
`(stable_id, side, sequence)` out of 26,777 mapping sequences.

### Taxonomy

`taxon`/`taxon_name` hold the full NCBI taxonomy as a **nested set**. `taxon2` is
strain-level (SARS-CoV-2, individual Zika/Influenza A strains, HHV-4/5 strains…), so
rolling up to species/family/order requires the nested-set bounds, not string matching:

```sql
JOIN taxon t ON d.left_value2 BETWEEN t.left_value AND t.right_value
JOIN taxon_name n ON n.taxon_id = t.taxon_id AND n.name_class = 'scientific name'
WHERE n.name = 'Orthomyxoviridae'
```

`taxon2` is the NCBI scientific name of `ncbi_taxon_id2`.

## Gotchas

- The `dataset` MV has **no indexes** — every query is a ~420k-row seq scan. Fine for
  one pass; think before self-joining it.
- One live description is absent from the MV (its `ncbi_taxon_id` is missing from
  `taxon`, dropped by an inner join).
- `character(N)` columns (`stable_id`, `type`, `type1`, `original_versionN`) are
  blank-padded; compare with care.
- The MV is materialized — `REFRESH MATERIALIZED VIEW dataset;` after base-table changes.
