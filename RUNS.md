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

## 2026-09-15 -- scripts/mapping_corrections.py

- output: `data/mapping_corrections.tsv`
- database: drakkar_2026_09
- valid descriptions with mappings: 18,825
- mapping entries audited: 25,925
- occurrences audited: 44,487
- problems: 75
- descriptions requiring correction: 58
- peptide problems: 45
- problem coordinates_mismatch: 36
- problem coordinates_out_of_bounds: 1
- problem duplicate_occurrence: 9
- problem empty_isoforms: 16
- problem repeated_isoform_block: 13

## 2026-09-22 -- scripts/graph_descriptions.py

- output: `data/graph-2026-09-09/descriptions.tsv, data/graph-2026-09-09/publications.tsv, data/graph-2026-09-09/peptides.tsv`
- database: drakkar_2026_09
- publications: 14,929
- peptides: 4,593

## 2026-09-22 -- scripts/viral_peptide_mature_fasta.py

- output: `data/viral_mature_with_peptides.fasta`
- database: drakkar_2026_09
- valid vh descriptions with a mapping: 8,513
- peptide occurrences: 3,003
- records: 1,609
- mature proteins: 764
- distinct sequences: 708
- occurrences recorded: 2,949
- occurrences recorded despite mismatch: 2
- occurrences relocated: 19
- occurrences searched: 33

## 2026-09-22 -- scripts/viral_peptide_summary.py

- output: `data/viral_peptide_summary.tsv`
- database: drakkar_2026_09
- valid vh descriptions with a mapping: 8,513
- unique viral peptide sequences: 1,430
- unique source viral (mature) proteins: 764
- unique human targets: 696
- unique (peptide, human target) associations: 2,471

## 2026-09-22 -- scripts/viral_protein_human_degree_binary.py

- output: `data/viral_protein_human_degree_binary.tsv + reports/viral_protein_human_degree_binary.md`
- database: drakkar_2026_09
- binary-supported pairs: 13,876
- viral proteins: 1,051
- degree: mean 13.2, median 2, max 430, p90 = 23, p99 = 230
- concentration: top 1% hold 26.6% of pairs, top 10% hold 75.9%

## 2026-09-22 -- scripts/dataset_integrity.py

- output: `data/dataset_integrity.tsv + reports/dataset_integrity.md`
- database: drakkar_2026_09
- valid descriptions audited: 403,121
- of them carrying a mapping: 18,825
- mapping entries audited: 25,925
- mapping occurrences audited: 44,487
- problems: 914
- descriptions requiring correction: 888
- mapping problems on peptides: 45
- problem description/duplicated: 216
- problem description/not_in_dataset_view: 1
- problem mapping/duplicated_occurrence: 9
- problem mapping/no_isoform: 16
- problem mapping/position_mismatch: 36
- problem mapping/position_off_sequence: 1
- problem mapping/repeated_isoform: 13
- problem protein/name_not_uniprot: 2
- problem protein/not_full_length: 620
