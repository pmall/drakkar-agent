# drakkar-agent

Assistant for querying the local **drakkar** protein–protein interaction database and
producing datasets. Curated human–human (`hh`) and human–viral (`vh`) PPIs.

## Repo

```
drakkar/db.py              connection + query helpers (loads .env via python-dotenv)
drakkar/descriptions.py    the valid-description filter, `valid()`
drakkar/mappings.py        reading mapping1/mapping2 occurrences
drakkar/peptides.py        peptide length bounds (5-20) and the sequences a peptide is read on
drakkar/taxonomy.py        rolling strain-level taxa up to species/genus/family
drakkar/binary_methods.py  shared binary detection-method definition (PSI-MI ids)
drakkar/stats.py           prints a dataset's counts, for its log entry
scripts/                   one dataset script per file
data/<database>/           exports, one folder per source database (gitignored)
```

Credentials live in `.env` (gitignored; see `.env.example`). Never print the password.

`drakkar` is installed into the venv as a package, so scripts run from anywhere:

```python
from drakkar.db import dataset_path, export, fetch, stream

df = fetch("select ... from dataset where ...")  # -> polars DataFrame
export("select ...", dataset_path("my_dataset.tsv"))  # .csv / .tsv / .parquet
stream("select ...", dataset_path("big_dataset.tsv"))  # COPY, never held in memory
```

```bash
uv run python scripts/my_dataset.py
```

`psql` is fine for quick counts and exploration.

## Datasets

**Every export goes to `data/<database>/`**, the database being the one in `.env`
(`drakkar.db.dataset_path` builds the path; `report_path` does the same under
`reports/`). A dataset is a function of its source database, so each database
version has its own folder and versions never overwrite each other; the same command
on the same database gives the same file. Check **Generic scripts** below before writing a
script; `scripts/example_dataset.py` is a template for the rare dataset none of them
covers.

`data/` is gitignored, so the script is the only record of how its dataset was built.
A one-off script is an append-only provenance record: nothing imports it, so it is
never edited to keep it running. The generic scripts are the exception — they are
maintained, on top of the shared `drakkar/` modules, so that one script serves every
dataset of its kind.

### Generic scripts

Parameterized scripts that cover the common requests. Use one, with its options,
rather than writing a new script; run them from the repo root with `uv run python`.
Each prints its counts. A table that is a few manipulations of a generic
script's output (a flag, a count per group, a filter, a join) is made from that
output ad hoc, in a polars one-liner, never as a new script.

| Script | Use it for | Input → output |
|---|---|---|
| `taxon_human_targets.py` | the human proteins a virus interacts with | `--ncbi-taxon-id` (covers every strain below it), `--output`, optional `--min-publications` / `--min-methods` (OR-ed) → `accession, name, descriptions, interactions, viral_proteins` |
| `human_target_peptides.py` | the peptides binding a list of human proteins | accessions on stdin (first field of each line), `--output` → one row per peptide, target, source protein and position, with the source and its sequence |
| `viral_protein_degree.py` | how many human interactors each viral protein has, strains collapsed to species | `--output`, optional `--only-binary`, `--min-publications` → `virus, name2, n_human` |
| `vh_interactions_with_sequences.py` | every `vh` interaction with sequences, pmids, methods and a binary flag | none → `vh_interactions_with_sequences.tsv` |
| `graph_descriptions.py` | the whole curated interactome: descriptions, viral proteins, peptides | `--output` folder → four TSVs |
| `dataset_integrity.py` | auditing malformed curated descriptions | none → `dataset_integrity.tsv` and a report |

Recipes:

- **Peptides of a virus**: `taxon_human_targets.py`, then its accession column into
  `human_target_peptides.py`. Peptides *on* viral proteins are the rows with
  `source_type = v`; the viral peptides of every virus come from running the
  targets script on taxon 10239 (Viruses).
- **Peptides of a gene list**: its accessions straight into `human_target_peptides.py`.
- **Peptide on a C-terminus**: `peptide_stop == len(source_sequence)`.
- **Per-sequence counts** (sources, targets): group the peptides output by `sequence`.
- **Viral proteins per virus**: count the rows of the degree output per `virus`.

Docstrings describe the dataset, never a specific run: run dates, database name and
counts do not go there. A script prints its counts (`drakkar.stats.show`) and logs
nothing; the dataset's entry is written by the agent, see **Dataset log**.

A script is finished once it has been run and **Verification** below passes.

### Conventions

Defaults for every dataset script; follow them without asking.

- **Protein identity.** A human protein is its accession alone; a viral protein is its
  `(accession, start, stop)` triple. Group and count on those only. Names are labels:
  never part of a key; aggregate them with `string_agg(DISTINCT name, ', ')` when
  grouping.
- **Taxon input.** A virus is given as an NCBI taxon id and covers every taxon below
  it (nested set, see **Taxonomy**). It selects the viral side it is asked for and
  nothing else: when the request says "all sources" or "any", every other step spans
  both interactomes and every virus.
- **Peptides.** 5–20 aa inclusive (see **Peptides**), identity below 100 kept, read
  through `drakkar.mappings.occurrences`. A peptide binding a protein is a
  mapping on that protein's *partner* side of the description, never on the protein
  itself; in `hh` the protein may be on either side.
- **Peptide columns** depend on the script: the distinct sequences alone, or with the
  source, the target, or both, optionally with coordinates on the source. A source is
  its `(accession, start, stop)` triple, coordinates are on the source sequence, and
  a peptide on the canonical sequence is reported there even if curated on an isoform
  (`occurrences` does this). If the request doesn't make the columns clear, ask.
- **Output files.** Sequences go in the last columns, since they make a file
  unreadable otherwise. File names are the dataset name with an underscore suffix
  (`descriptions_hh.tsv`), parameterized ones named after the parameter
  (`ebolavirus_human_targets.tsv`).
- **CLI.** Script parameters are named options (`--ncbi-taxon-id`), not positionals.
  `--output` is a file name (or folder) inside `data/<database>/`, never a path.

### Dataset log

Each `data/<database>/` folder holds one `DATASETS.md`, uncommitted like the data
beside it, that explains how every dataset of that database was made. A script is only
part of that: a dataset may come from several scripts piped together (the accessions
of one fed to another, a list read from a file), from options that make the file's
name no clue to its content, or from manipulations of a script's output. Scripts log
nothing, so the agent writes the entry when a dataset is delivered, headed
`## <file or folder> -- <date>`:

- what the dataset is, in terms of the data;
- the exact commands, in order, with their inputs and any ad hoc step after;
- the counts the scripts printed.

The same data re-laid-out (columns, sort order) keeps its entry; a change in the
counts is a new entry. A dataset file with no entry is unexplained: add the entry or
delete the file.

### Working with the user

- Ask only about choices that change the output, phrased in terms of the data ("peptides
  on viral proteins of any virus, or only of this taxon?"), and apply the answer
  literally.
- Report what the output is. Checks that turned out fine and change nothing are not
  reported.
- Before calling an output row odd, read how `drakkar.mappings` documents the case —
  most are handled by design.
- Never edit `drakkar/` modules without the user's consent.
- **Never run `git commit`.** Finished work is left as changes in
  the working tree; the user commits, and writes the message. Nothing in this file is
  permission to commit.
- Commit messages carry no co-authorship or attribution line.

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

Runs > descriptions > interactions.

A **run** is a batch of publications the curation team processed together
(`run`/`run_id`).

One row per **description** = *(publication × detection method × protein 1 region ×
protein 2 region)*. `stable_id` identifies a description across its revisions;
`deleted_at IS NULL` selects the current revision.

An **interaction** is a protein pair: *(human accession, human accession)*, unordered,
in `hh`, *(human accession, viral triple)* in `vh` (see **Mature viral proteins**).
Many descriptions — publications, methods — support one interaction. "How many
interactions" counts distinct pairs, never description rows; a protein's interactions
are its distinct partners.

A list of interactions spanning both `hh` and `vh` always gives the partner as the
three columns `accession, start, stop`, human partners included: their `start` /
`stop` are either NULL or `(1, length)`.

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
it was curated on. An occurrence curated on an isoform is reported on the canonical
sequence whenever its sequence is there; it stays on the isoform only otherwise.

Odd-looking occurrences are usually documented cases, flagged by `how` — e.g. a
curated peptide carrying an initiator Met its mature protein lacks is kept at its
recorded position as `recorded despite mismatch`.

#### Peptides

A **peptide** is a mapping `sequence` of **5–20 aa inclusive** — a short binding
region, as opposed to a domain-scale mapping. Nothing in the schema flags it; derive
it by length:

```sql
WHERE length(e->>'sequence') BETWEEN 5 AND 20
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
