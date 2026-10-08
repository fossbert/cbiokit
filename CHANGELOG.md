# Changelog

## 0.4.0

### Added
- `drivers_only=False` for `driver_event_matrix`, `driver_gene_matrix` and `alteration_matrix`:
  all alterations of the chosen types, not only those labelled '(driver)' (variants of unknown
  significance included, needed for locus-specific variant analyses). Gene-level columns of
  `driver_gene_matrix` are then named `<GENE>_ANY`. On TCGA STAD + EAC: 1116 instead of 649
  mutation events.

## 0.3.0

### Added
- `fusion_events`, `fusion_matrix`: gene fusions, which the exports never label as drivers. Events
  are parsed to 5'/3' partners ('A-B', 'A-B fusion', reported under one or both partners) and
  deduplicated; `groups` merges variants into one row (e.g. `{"CLDN18-ARHGAP6/26":
  r"CLDN18-ARHGAP(6|26)"}`); absence is only called where all queried partners were profiled for
  structural variants. Non-fusion structural variants (intragenic, deletions within a
  transcript) are excluded. On TCGA STAD + EAC: 10 CLDN18-ARHGAP6/26 fusions, identical to the
  former notebook routine.

## 0.2.0

### Added
- `concat_exports`: combine exports of several studies; genes missing from an export count as
  not profiled for its samples.
- `alteration_matrix`: binary driver matrix rows x samples (gene or `GENE:EVENT` level, chosen
  alteration types, samples or patients as columns, NaN = not profiled). Checked on TCGA
  STAD + EAC exports: identical to the former notebook routine (0 of 23,636 values differ).

## 0.1.0

- Extracted from `chordutils.cbioportal` 0.3.0 (`read_alteration_export`, `driver_event_matrix`,
  `driver_gene_matrix`, `panel_genes_from_export`, `AlterationExport`), unchanged in behaviour.
