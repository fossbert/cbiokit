# cbiokit

Helpers for [cBioPortal](https://www.cbioportal.org) exports. First module:
driver (OncoKB / hotspot) annotations from the *Alterations across samples*
download, parsed into events, driver matrices and per gene-and-type profiling
information ("not profiled" stays missing, never wild type).

Extracted from `chordutils.cbioportal`; `chordutils` now imports from here.

## Install

```bash
pip install -e '.[test]'
```

## Usage

```python
import cbiokit as cbk

ex = cbk.read_alteration_export("alterations_across_samples.tsv")
events = cbk.driver_event_matrix(ex)      # samples x 'KRAS:G12D', 'ERBB2:AMP', ...; NaN = not profiled
genes = cbk.driver_gene_matrix(ex)        # <GENE>_DRIVER, <GENE>_SV per sample
```

See the module docstring of `cbiokit.alterations` for the export format and
the checks made on twelve real exports (MSK-CHORD, MSK-MET, GENIE BPC, TCGA, CPTAC).

Exports are derived study data: keep them out of version control when the study
licence forbids redistribution.
