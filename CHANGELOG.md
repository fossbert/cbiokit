# Changelog

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
