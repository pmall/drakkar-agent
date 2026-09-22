# drakkar-agent

Assistant for querying the local **drakkar** protein–protein interaction database and
producing datasets. Curated human–human (`hh`) and human–viral (`vh`) PPIs.

## Repo

```
drakkar/db.py              connection + query helpers (loads .env via python-dotenv)
drakkar/mappings.py        reading mapping1/mapping2 occurrences
drakkar/binary_methods.py  shared binary detection-method definition (PSI-MI ids)
drakkar/runs.py            append-only run log helper
scripts/                   one dataset script per file
data/                      exports (gitignored)
RUNS.md                    committed log of dataset runs (date + counts per run)
```

Credentials live in `.env` (gitignored; see `.env.example`). Never print the password.

`drakkar` is installed into the venv as a package, so scripts run from anywhere:

```python
from drakkar.db import export, fetch, stream

df = fetch("select ... from dataset where ...")  # -> polars DataFrame
export("select ...", "data/my_dataset.tsv")  # .csv / .tsv / .parquet
stream("select ...", "data/big_dataset.tsv")  # COPY, never held in memory
```

```bash
uv run python scripts/my_dataset.py
```

`psql` is fine for quick counts and exploration.

## Datasets

One script per dataset, **all exports to `data/`**; `scripts/example_dataset.py` is a
working template.

`data/` is gitignored, so the script is the only record of how its dataset was built.
Treat `scripts/` as an append-only provenance log: commit every dataset script, one
commit per dataset. Nothing imports them, so they are never edited to keep them
running — a new dataset is a new script.

Docstrings describe the dataset, never a specific run: run dates, database name and
counts go to the committed `RUNS.md` via `drakkar.runs.log_run`.

A script is finished once it has been run and **Verification** below passes.

## Writing code

`uv` is the package manager and Python is 3.14. Use `uv add` / `uv sync` / `uv run`,
never bare `pip` or `python` — bare `python` resolves to a pyenv version that is not
installed and errors out.

Type everything — parameters, returns, attributes — and keep annotations precise
rather than convenient: a literal or an enum over a bare `str`, a concrete collection
type over `Any`.

`None` and `Optional` are statements about the domain, not escape hatches. Use them
only where absence is genuinely possible in the data, never to satisfy the type
checker, silence an error, or stand in for a value that is not initialized yet. An
awkward type usually means the shape of the code is wrong.

## Verification

At the end of every coding session, in order, fixing what they surface rather than
working around it:

1. Delete dead code, and imports the change introduced or left orphaned.
2. `uv run ruff format .`
3. `uv run ruff check --fix .`
4. `uv run pyright`
5. `uv run python -m compileall -q drakkar scripts` — add other code roots as they
   appear.
6. Exercise the entrypoints the change touched. There are no unit tests: this is
   query and export code, verified by running it.
7. Re-read prose the change added — docs, docstrings, comments — with fresh eyes, and
   make it match the voice of what is already there.
8. Read the diff and confirm it holds only intentional changes.

## Data model

Curation source tables: `runs → associations → descriptions`, plus `methods`,
`proteins` (+ `proteins_versions`), `publications`, `taxon`/`taxon_name`, `keywords`.

**The database is read-only.** Only ever `SELECT` from it: no `INSERT` / `UPDATE` /
`DELETE`, no DDL, no `REFRESH`. The curation team owns the data; a new version of it
arrives as a new database.

Everything needed for datasets is denormalized in the **`dataset` materialized view**.
Query that; drop to the base tables only for things it doesn't carry (sequences,
publication metadata, UniProt features).

Schema varies between database versions — verify a table exists before depending on it
rather than assuming this file is exhaustive.

### Grain

One row per **description** = *(publication × detection method × protein 1 region ×
protein 2 region)*. `stable_id` identifies a description across its revisions;
`deleted_at IS NULL` selects the current revision.

A **run** is a batch of publications the curation team processed together
(`run`/`run_id`).

### Valid PPI description — always apply

```sql
WHERE state = 'curated'
  AND is_obsolete1 IS FALSE
  AND is_obsolete2 IS FALSE
  AND deleted_at IS NULL
```

Use this filter for every dataset unless explicitly told otherwise; if a request needs
it relaxed, say so. Restrict to one interactome with `AND type = 'hh'` or
`AND type = 'vh'` (equivalent to filtering on `type1`/`type2`, since side 1 is always
human — prefer the `type` form).

- `state`: `curated` / `selected` / `discarded` (`pending` never reaches the view).
- `is_obsoleteN`: the accession no longer exists in current UniProt
  (`original_versionN` = release at curation time, `current_versionN` = current
  release, NULL when obsolete).

### Sides

Protein 1 is **always human** (`type1 = 'h'`, taxon 9606). Protein 2 is human (`hh`)
or viral (`vh`) — check `type` / `type2`.

| | side 1 (human) | side 2 |
|---|---|---|
| accession | `accession1` | `accession2` |
| name | `name1` — UniProt gene name | `name2` — **curator-chosen** |
| region | `start1`/`stop1` (always full length, `start1 = 1`) | `start2`/`stop2` |
| taxon | 9606, `taxon1 = 'Homo sapiens'` | `ncbi_taxon_id2`, `taxon2`, `left_value2`/`right_value2` |

Human proteins are full length, so `accession1` alone identifies one.

### Mature viral proteins

Viral proteins are **mature proteins**: subsequences of a polyprotein. A viral
protein's identity is therefore the triple

```
(accession2, start2, stop2)   -- polyprotein accession + coordinates on it
```

with `name2` the name assigned by the curation team. **Never key viral proteins on
`accession2` alone** — one polyprotein yields many mature proteins (e.g. `P0DTD1`
SARS-CoV-2 ORF1ab → 15 mature proteins). `name2` is stable per triple and not reused
across different coordinates, so `(accession2, name2)` also identifies a mature
protein in practice — but prefer the coordinate triple.

To recover a mature protein sequence, substring the polyprotein from
`proteins.sequences`:

```sql
substring(p.sequences->>d.accession2 FROM d.start2 FOR d.stop2 - d.start2 + 1)
```

### Mappings — binding regions

`mapping1` / `mapping2` (JSON) are the subsequences curators reported as **sufficient
to bind the partner**. Empty `[]` when the publication reported no such region (the
common case).

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
  - **Viral side:** for a *mature* protein (`start2 > 1`) it is **always** the
    canonical accession — the curation interface does not allow a mature protein to be
    defined on an isoform. Non-canonical viral isoform references therefore occur only
    for full-length viral proteins (`start2 = 1`), e.g. Tat, E1A, ICP22, LT, HBZ,
    UL37, VSV M, HEV ORF2. A polyprotein entry may still *have* isoforms in UniProt
    (FMDV, PRRSV, FIPV frameshift/alternative-ORF products); they are simply never
    used as a mature-protein frame.
- `identity` is a percentage; it is **not always 100** — the peptide is aligned, not
  necessarily exact. Filter on it when exactness matters.
- **Occurrence coordinates are 1-based relative to the protein as curated, i.e. to
  `startN`, not to the polyprotein.** For mature viral proteins the absolute position
  on the polyprotein is `start2 + occurrence.start - 1`. On the human side
  `start1 = 1`, so the two frames coincide.

Read them through `drakkar.mappings.occurrences` rather than by hand: it normalizes
the JSON (coordinates and identities are stored as strings about a third of the time),
flattens it to one row per occurrence, and corrects each position against the sequence
it was curated on.

#### Peptides

By team convention a **peptide** is a mapping `sequence` of **4–20 aa inclusive** — a
short binding region, as opposed to a domain-scale mapping. Nothing in the schema
flags it; derive it by length:

```sql
WHERE length(e->>'sequence') BETWEEN 4 AND 20
```

### Taxonomy

`taxon`/`taxon_name` hold the full NCBI taxonomy as a **nested set**. `taxon2` is the
NCBI scientific name of `ncbi_taxon_id2`, and it is strain-level (SARS-CoV-2,
individual Zika/Influenza A strains, HHV-4/5 strains…), so rolling up to
species/family/order requires the nested-set bounds, not string matching:

```sql
JOIN taxon t ON d.left_value2 BETWEEN t.left_value AND t.right_value
JOIN taxon_name n ON n.taxon_id = t.taxon_id AND n.name_class = 'scientific name'
WHERE n.name = 'Orthomyxoviridae'
```

## Gotchas

- Join `proteins` through `dataset.protein1_id` / `protein2_id` (`proteins.id`), never
  on `accession` — the table holds several rows per accession (UniProt versions) and
  an accession join duplicates dataset rows.
- Internal ids (`methods.id`, `proteins.id`, …) are join keys only. Never hardcode
  them as meaningful values: match detection methods on the stable `psimi_id` (see
  `drakkar/binary_methods.py`), and resolve proteins through the dataset's
  `proteinN_id` links.
- The `dataset` MV has **no indexes** — every query is a ~420k-row seq scan. Fine for
  one pass; think before self-joining it.
- One live description is absent from the MV (its `ncbi_taxon_id` is missing from
  `taxon`, dropped by an inner join).
