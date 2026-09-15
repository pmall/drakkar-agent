# Dataset runs

Append-only log of dataset exports. Each run appends one section below via
`drakkar.runs.log_run` -- never edit past entries.

## 2026-09-02 -- scripts/viral_peptide_mature_fasta.py

- output: `data/viral_mature_with_peptides.fasta`
- database: drakkar-other-version
- valid vh descriptions: 116,015
- peptide occurrences: 2,467 (2,446 recorded, 19 relocated, 2 skipped)
- records: 1,264 over 625 mature proteins, 571 distinct sequences

## 2026-09-04 -- scripts/viral_proteins_by_virus.py

- output: `data/viral_proteins_by_virus.tsv`
- database: drakkar-other-version
- species: 354 viruses, 1,734 unique (virus, name) entries
- genus: 148 genera, 1,258 unique (genus, name) entries
- family: 51 families, 977 unique (family, name) entries

## 2026-09-04 -- scripts/viral_protein_human_degree.py

- output: `data/viral_protein_human_degree.tsv`
- database: drakkar-other-version
- strain-abstracted viral proteins: 1,734
- unique (viral, human) pairs: 90,683
- degree: mean 52.3, median 4, max 2,817, p90 = 126, p99 = 814
- concentration: top 1% hold 24.5% of pairs, top 10% hold 78.2% (hub-dominated)

## 2026-09-04 -- scripts/viral_protein_human_degree_multipmid.py

- output: `data/viral_protein_human_degree_multipmid.tsv + reports/viral_protein_human_degree_multipmid.md`
- database: drakkar-other-version
- pairs with >1 pmid: 7,055
- viral proteins: 460
- degree: mean 15.3, median 2, max 411, p90 = 34, p99 = 240
- concentration: top 1% hold 24.0% of pairs, top 10% hold 72.0%

## 2026-09-04 -- scripts/viral_protein_human_degree_binary.py

- output: `data/viral_protein_human_degree_binary.tsv + reports/viral_protein_human_degree_binary.md`
- database: drakkar-other-version
- binary-supported pairs: 13,626
- viral proteins: 950
- degree: mean 14.3, median 2, max 430, p90 = 27, p99 = 271
- concentration: top 1% hold 25.4% of pairs, top 10% hold 75.2%

## 2026-09-04 -- scripts/viral_peptide_summary.py

- output: `data/viral_peptide_summary.tsv`
- database: drakkar-other-version
- valid vh descriptions scanned: 116,015
- unique viral peptide sequences: 1,111
- unique source viral (mature) proteins: 628
- unique human targets: 684
- unique (peptide, human target) associations: 2,011

## 2026-09-10 -- scripts/vh_interactions_with_sequences.py

- output: `data/vh_interactions_with_sequences.tsv`
- database: drakkar_2026_09
- interactions: 103,758
- binary interactions: 15,656

## 2026-09-15 -- scripts/viral_peptide_cter_flag.py

- output: `data/viral_peptide_cter_flag.tsv`
- database: drakkar_2026_09
- valid vh descriptions scanned: 117,431
- peptide occurrences resolved: 2,866
- distinct peptide sequences: 1,359
- unresolved: 16
- is_cter true: 227
