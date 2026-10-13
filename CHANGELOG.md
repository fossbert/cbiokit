# Changelog

## 0.8.0

### Changed (behaviour)
- `alteration_matrix` now counts **all alterations by default** (`drivers_only=False`); before, only events labelled
  '(driver)' counted. The driver label exists only for annotated genes (e.g. RHOA deletions and VUS have none), so it
  is no neutral filter. Code that relied on the old default must pass `drivers_only=True` or use `export.drivers()`.
  `driver_event_matrix` and `driver_gene_matrix` keep their driver default (their names say so; `chordutils` uses them).

### Added
- `AlterationExport.drivers()`: copy of the export with only the driver events (samples and profiling unchanged);
  works with every matrix function. `fusion_matrix` of such an export is empty instead of an error.

## 0.7.0

### Added
- `calc_combo_muts`, `make_combo_df` (module `cbiokit.combos`): count co-occurring alterations per
  patient (all combinations of size >= 2 above `min_occurence`) and build the binary combination
  table for them. Input is patients x alterations (0/1), e.g. `alteration_matrix(...).T`. New
  options `max_size` (limit combination size, the number of combinations grows exponentially) and
  `sep` (id separator, for names containing '_').

## 0.6.0

### Added
- `variant_class`, `variant_class_matrix`: classes of mutations from the protein change in the
  export (truncating = nonsense, frameshift, splice site, stop-inserting; inframe; missense;
  other) and a binary matrix with rows such as `CDH1:truncating`, with the usual NaN = not
  profiled convention. For comparing mutation classes of a gene, e.g. with
  `pyrea.compare_phenotypes`. The classes are a heuristic on the event text; on the STAD, EAC,
  BRCA and GBM exports 6 of 2,807 distinct variants fall into `other` (start loss, stop loss,
  generic labels).

## 0.5.1

### Fixed
- Copy number: cBioPortal's `<GENE>: AMP` and `<GENE>: HOMDEL` columns are identical (both list all
  copy number events), so `types=("AMP",)` also counted deletions and `types=("HOMDEL",)` also
  amplifications. The type of a copy number event is now taken from the event ('AMP' or 'HOMDEL'),
  duplicates are dropped. Checked on STAD, EAC, BRCA and GBM exports: e.g. GBM CDKN2A had 322
  deletions counted as 'AMP'. The default `types=("MUT", "AMP", "HOMDEL")` and all MUT/fusion
  results were not affected; ERBB2 (STAD, BRCA) and EGFR (GBM) amplification counts are unchanged
  because those genes have no deletions.

## 0.5.0

### Added
- `fusion_matrix(groups=...)`: a group can be a pair of regexes `(5' partner, 3' partner)`
  (None = any gene), e.g. `(r"CLDN\d+", r"ARHGAP\d+")`; `ignore_orientation=True` also accepts
  the reverse orientation. Plain regexes work as before.
- `fusion_summary`: table of the fusions in an export (partners, number of samples/patients,
  queried genes), optionally filtered by partner regexes.

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
