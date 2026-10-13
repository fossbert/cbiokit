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

## All alterations by default, drivers on request

cBioPortal labels some alterations '(driver)' (OncoKB / hotspots). The label exists only for annotated genes and
alteration types, so it is not a neutral filter: e.g. `RHOA` deletions or variants of unknown significance carry no
label and would silently disappear. **The matrices therefore contain all alterations by default.** To restrict to the
labelled drivers, filter the export:

```python
ex = cbk.read_alteration_export("alterations_across_samples.tsv")

anymut = cbk.alteration_matrix(ex, types=("MUT",))                          # any non-silent mutation (VUS included)
variants = cbk.alteration_matrix(ex, types=("MUT",), level="event")         # 'GENE:EVENT' rows, all variants
drivers = cbk.alteration_matrix(ex.drivers(), types=("MUT",), level="event")   # only events labelled '(driver)'
```

`ex.drivers()` returns a copy of the export with only the driver events (samples and profiling unchanged, so 'not
profiled' stays missing and unannotated genes become 0); every function accepts it. `drivers_only=True` in
`alteration_matrix` is a shorthand for the same. Fusions are never labelled as drivers, so `fusion_matrix(ex.drivers())`
is empty.

The functions named `driver_event_matrix` and `driver_gene_matrix` keep their original meaning (default
`drivers_only=True`, the OncoKB matrices of `chordutils`).

An event matrix has a 0 for a variant in every sample that does not carry it, including carriers of *other* variants of
the gene; those are not wild type. When comparing variants with wild type, use `pyrea.locus_specific_mps`, which derives
the gene state itself (see the pyrea README), and do not drop rare variants beforehand (`min_samples` in
`alteration_matrix`).

## Mutation classes

```python
cbk.variant_class("S70Pfs*13")                           # 'truncating' (also nonsense, splice site)
classes = cbk.variant_class_matrix(ex, by="PATIENT_ID")   # rows 'CDH1:truncating', 'CDH1:missense', 'CDH1:inframe', ...
classes.loc[classes.index.str.startswith("CDH1:")].sum(axis=1)
```

All mutations count by default (`drivers_only=False`): missense variants are mostly unlabelled. The
classes come from the text of the protein change, so look at the events of your gene if it matters.
A sample with a mutation of another class is 0 in a class row; to compare a class with wild type,
set such samples to NaN or compare all classes together (`pyrea.compare_phenotypes`).

## Co-occurring alterations

```python
mat = cbk.alteration_matrix(ex, by="PATIENT_ID").T.fillna(0)     # patients x alterations, 0/1
combos = cbk.calc_combo_muts(mat, min_occurence=10, max_size=3)   # Series 'KRAS_TP53' -> n patients
combo_df = cbk.make_combo_df(mat, combos.index)                   # patients x combinations, 0/1
```

A combination counts for a patient if all its members are present (other alterations may be too).
`max_size` limits the combination size; without it, patients with many alterations produce
exponentially many combinations. If alteration names contain '_', pass `sep=` to
`calc_combo_muts` and the same value as `splitter=` to `make_combo_df`.

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
