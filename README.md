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

# several cohorts, rows x samples with NaN = not profiled (input for e.g. pyrea's MPS)
ex = cbk.concat_exports([cbk.read_alteration_export(f) for f in files])
muts = cbk.alteration_matrix(ex, types=("MUT",), by="PATIENT_ID", min_samples=10)   # genes x patients
amps = cbk.alteration_matrix(ex, types=("AMP",), level="gene")
variants = cbk.alteration_matrix(ex, types=("MUT",), level="event")                 # 'GENE:EVENT' rows
fusions = cbk.fusion_matrix(ex, groups={"CLDN18-ARHGAP6/26": r"CLDN18-ARHGAP(6|26)"})  # gene fusions
```

See the module docstring of `cbiokit.alterations` for the export format and
the checks made on twelve real exports (MSK-CHORD, MSK-MET, GENIE BPC, TCGA, CPTAC).

Exports are derived study data: keep them out of version control when the study
licence forbids redistribution.

## Drivers only, or all variants?

cBioPortal labels some alterations '(driver)' (OncoKB / hotspots). By default the matrices only
contain these. Pass `drivers_only=False` to include every alteration of the chosen `types`,
variants of unknown significance included:

```python
drivers = cbk.alteration_matrix(ex, types=("MUT",), level="event")                       # 649 rows (STAD+EAC)
allvars = cbk.alteration_matrix(ex, types=("MUT",), level="event", drivers_only=False)   # 1116 rows
anymut  = cbk.alteration_matrix(ex, types=("MUT",), level="gene",  drivers_only=False)   # any non-silent mutation
```

An event matrix has a 0 for a variant in every sample that does not carry it, including carriers
of *other* variants of the gene; those are not wild type. When comparing variants with wild type,
use `pyrea.locus_specific_mps`, which derives the gene state itself (see the pyrea README), and do
not drop rare variants beforehand (`min_samples` in `alteration_matrix`).

## Gene fusions

Fusions are never labelled as drivers and therefore have their own matrix:

```python
# which fusions does the export contain? (to decide what to group)
cbk.fusion_summary(ex, five=r"CLDN\d+", ignore_orientation=True)

# one row per group; a pattern is a regex on 'GENE5-GENE3' or a pair of regexes (5', 3')
fus = cbk.fusion_matrix(ex, by="PATIENT_ID", ignore_orientation=True,
                        groups={"CLDN-ARHGAP": (r"CLDN\d+", r"ARHGAP\d+"),     # any claudin with any ARHGAP
                                "CLDN18-ARHGAP6/26": r"CLDN18-ARHGAP(6|26)"})  # first matching group wins
```

A fusion that the export lists under both partners counts once per sample; 'A-B' and 'B-A' are
different fusions unless `ignore_orientation=True`; `None` in a pair means any gene.

**Only fusions listed under a queried gene are in an export.** If the download was made with
CLDN18 and ARHGAP26 only, claudin-ARHGAP fusions involving other genes cannot be found. For a
general analysis (or another tumour type) download the export with all family members in the
gene query (e.g. every CLDN* and ARHGAP*).
