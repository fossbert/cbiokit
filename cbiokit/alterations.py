"""Alterations from cBioPortal "alterations across samples" exports.

cBioPortal lets you download, for a gene query, one table per study via
*Download -> Alterations across samples*. Alterations are annotated as "driver" (by default
based on OncoKB and cancer hotspots). This module parses such an export into events, driver
matrices and per-gene/type profiling information. Variants of unknown significance, which
raw mutation files contain, are separable from drivers this way.

Export format
-------------
One row per sample with ``Study ID``, ``Sample ID``, ``Patient ID``, ``Altered`` and, per
queried gene, a summary column ``<GENE>`` plus typed columns ``<GENE>: MUT``,
``<GENE>: AMP``, ``<GENE>: HOMDEL`` and ``<GENE>: FUSION`` (some studies lack a type, e.g. no
FUSION column without structural-variant data). A cell holds ``no alteration``,
``not profiled`` or the events joined by ', ', each optionally followed by ``(driver)``,
e.g. ``R175H (driver), P322Hfs*23``.

Checked on twelve exports (MSK-CHORD, MSK-MET, GENIE BPC CRC, TCGA CRC/STAD/EAC/BRCA/GBM/LUAD/
LUSC, CPTAC COAD), 2026-09-30:

- The summary column is always the union of the typed columns, so only the typed columns are
  used here and the type of every event is known (no guessing from the text).
- ``not profiled`` differs between types: e.g. samples without mutation data but with copy
  number data show ``no alteration`` in the summary column. Profiling is therefore tracked
  per gene AND type; 'not profiled' is missing, never 'wild type'.
- Structural variants are written in many forms ('X-Y fusion', 'X-Y', 'APC-intragenic',
  'Deletion within transcript: mid-exon', ...) and are never labelled '(driver)'. They are
  reported separately (``<GENE>_SV``), not as drivers.

The exports are derived data of the respective studies; keep them out of version control
when the study licence forbids redistribution (e.g. MSK-CHORD: CC BY-NC-ND 4.0).
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Sequence, Union

import numpy as np
import pandas as pd

TYPES = ("MUT", "AMP", "HOMDEL", "FUSION")
DRIVER_TAG = " (driver)"
NO_ALTERATION = "no alteration"
NOT_PROFILED = "not profiled"
_META = {"Study ID": "STUDY_ID", "Sample ID": "SAMPLE_ID", "Patient ID": "PATIENT_ID"}


@dataclass
class AlterationExport:
    """Parsed cBioPortal export.

    Attributes
    ----------
    samples : pd.DataFrame
        SAMPLE_ID, PATIENT_ID, STUDY_ID in file order.
    events : pd.DataFrame
        One row per alteration: SAMPLE_ID, GENE, TYPE (MUT/AMP/HOMDEL/FUSION), EVENT (e.g.
        'G12D', 'AMP', 'TRAP1-CREBBP fusion'), DRIVER (bool).
    profiled : pd.DataFrame
        Boolean, index SAMPLE_ID, columns MultiIndex (GENE, TYPE): was the gene profiled for
        this alteration type in this sample. Types missing from the export are False.
    """

    samples: pd.DataFrame
    events: pd.DataFrame
    profiled: pd.DataFrame

    @property
    def genes(self) -> list:
        return list(dict.fromkeys(self.profiled.columns.get_level_values(0)))

    def __repr__(self):
        n_drv = int(self.events["DRIVER"].sum())
        return (f"AlterationExport({len(self.samples)} samples, {len(self.genes)} genes, "
                f"{len(self.events)} events, {n_drv} drivers)")


def _split_events(cell: str):
    """'R175H (driver), V157F' -> [('R175H', True), ('V157F', False)]; none for the status values."""

    if cell in (NO_ALTERATION, NOT_PROFILED) or pd.isna(cell):
        return []
    out = []
    for ev in cell.split(", "):
        driver = ev.endswith(DRIVER_TAG)
        out.append((ev[: -len(DRIVER_TAG)] if driver else ev, driver))
    return out


def read_alteration_export(path: Union[str, Path]) -> AlterationExport:
    """Read a cBioPortal 'alterations across samples' export (see module docstring).

    Parameters
    ----------
    path : str or Path
        The downloaded .tsv file.

    Returns
    -------
    AlterationExport
    """

    raw = pd.read_table(path, dtype=str, keep_default_na=False)
    missing = [c for c in _META if c not in raw.columns]
    if missing:
        raise ValueError(f"Not a cBioPortal alteration export, missing columns: {missing}")

    samples = raw[list(_META)].rename(columns=_META)[["SAMPLE_ID", "PATIENT_ID", "STUDY_ID"]]
    if samples["SAMPLE_ID"].duplicated().any():
        raise ValueError("Sample IDs are not unique (export from several studies?)")
    sids = samples["SAMPLE_ID"].to_numpy()

    typed = [c for c in raw.columns if ": " in c and c.split(": ", 1)[1] in TYPES]
    genes = list(dict.fromkeys(c.split(": ", 1)[0] for c in typed))
    if not genes:
        raise ValueError("No typed columns ('<GENE>: MUT', ...) found; re-download the export from cBioPortal")

    rows, profiled = [], {}
    for g in genes:
        for t in TYPES:
            col = f"{g}: {t}"
            if col not in raw.columns:
                profiled[(g, t)] = np.zeros(len(raw), dtype=bool)
                continue
            values = raw[col].to_numpy()
            profiled[(g, t)] = values != NOT_PROFILED
            for sid, cell in zip(sids, values):
                for ev, drv in _split_events(cell):
                    rows.append((sid, g, t, ev, drv))

    events = pd.DataFrame(rows, columns=["SAMPLE_ID", "GENE", "TYPE", "EVENT", "DRIVER"])
    prof = pd.DataFrame(profiled, index=pd.Index(sids, name="SAMPLE_ID"))
    prof.columns = pd.MultiIndex.from_tuples(prof.columns, names=["GENE", "TYPE"])
    return AlterationExport(samples.reset_index(drop=True), events, prof)


def panel_genes_from_export(export: AlterationExport, gene_panel: pd.Series,
                            alteration_type: str = "MUT") -> Dict[str, list]:
    """Genes on each sequencing panel, as observed in an export.

    A gene counts as being on a panel if it was profiled (``alteration_type``) in EVERY sample
    of that panel in the export. The MSK-CHORD download itself has no panel
    gene lists; the result can be passed as ``panel_genes`` to
    ``chordutils.covariates.genomic_features``.
    Only genes of the export's gene query are covered.

    Parameters
    ----------
    export : AlterationExport
    gene_panel : pd.Series
        Panel per sample, index SAMPLE_ID (e.g. ``data_gene_panel_matrix.txt``, column
        'mutations').
    alteration_type : str
        Profiling of which type decides.

    Returns
    -------
    dict of {panel: list of genes}
        Panels without any sample in the export are omitted.
    """

    prof = export.profiled.xs(alteration_type, axis=1, level="TYPE")
    panel = gene_panel.reindex(prof.index)
    out = {}
    for p, sub in prof.groupby(panel):
        out[p] = [g for g in prof.columns if sub[g].all()]
    return out


def driver_event_matrix(export: AlterationExport, types: Sequence[str] = ("MUT", "AMP", "HOMDEL"),
                        min_samples: int = 1, order: bool = True) -> pd.DataFrame:
    """Samples x driver events ('GENE:EVENT', e.g. 'KRAS:G12D', 'ERBB2:AMP').

    Parameters
    ----------
    export : AlterationExport
    types : sequence of {'MUT', 'AMP', 'HOMDEL', 'FUSION'}
        Alteration types included. FUSION events are never labelled as drivers in the exports
        checked so far; see ``driver_gene_matrix`` for structural variants.
    min_samples : int
        Events present in fewer samples are dropped.
    order : bool
        Order genes by their total number of driver events and, within a gene, events by
        frequency (both descending), e.g. for oncoprints.

    Returns
    -------
    pd.DataFrame
        Index SAMPLE_ID; 1 = event present, 0 = absent, missing = the gene was not profiled
        for the event's type in this sample.
    """

    ev = export.events[export.events["DRIVER"] & export.events["TYPE"].isin(types)]
    ev = ev.assign(COLUMN=ev["GENE"] + ":" + ev["EVENT"]).drop_duplicates(["SAMPLE_ID", "COLUMN"])

    counts = ev.groupby(["GENE", "COLUMN"]).size().rename("n").reset_index()
    counts = counts[counts["n"] >= min_samples]
    if order:
        gene_total = counts.groupby("GENE")["n"].sum()
        counts = counts.assign(_g=counts["GENE"].map(gene_total)) \
                       .sort_values(["_g", "GENE", "n", "COLUMN"], ascending=[False, True, False, True])

    columns = pd.Index(counts["COLUMN"].tolist())
    index = export.profiled.index
    values = np.zeros((len(index), len(columns)))
    hits = ev[ev["COLUMN"].isin(columns)]
    values[index.get_indexer(hits["SAMPLE_ID"]), columns.get_indexer(hits["COLUMN"])] = 1.0
    m = pd.DataFrame(values, index=index, columns=columns)

    event_type = ev.drop_duplicates("COLUMN").set_index("COLUMN")[["GENE", "TYPE"]]
    for c in columns:
        g, t = event_type.loc[c]
        m[c] = m[c].where(export.profiled[(g, t)])
    return m


def driver_gene_matrix(export: AlterationExport, types: Sequence[str] = ("MUT", "AMP", "HOMDEL"),
                       gene_groups: Optional[Dict[str, Sequence[str]]] = None,
                       fusions: bool = True) -> pd.DataFrame:
    """One row per sample: driver status per gene (and gene group), plus structural variants.

    Columns
    -------
    <GENE>_DRIVER
        1 if the gene has a driver event of one of ``types``; 0 if it was profiled for ALL of
        ``types`` and has none; missing otherwise (e.g. no mutation data for the sample).
    <GROUP>_DRIVER, <GROUP>_N_DRIVER
        For ``gene_groups`` (e.g. ``config.gene_groups``): 1 if any gene of the group has a
        driver, 0 if all genes are 0, else missing; N_DRIVER = number of genes with a driver.
    <GENE>_SV (with ``fusions=True``)
        1 if the export lists any structural variant for the gene (fusion, intragenic
        deletion/duplication, ...; the exports carry no driver label for these), 0 if
        profiled without one, missing if not profiled. Inspect the individual events in
        ``export.events.query("TYPE == 'FUSION'")`` before using a column.

    Returns
    -------
    pd.DataFrame
        SAMPLE_ID + the columns above. Column names do not overlap with
        ``chordutils.covariates.genomic_features`` (``_ALT``, ``_N``, ``_FUSION``), so both can
        be merged on SAMPLE_ID.
    """

    prof = export.profiled
    ev = export.events
    drv = ev[ev["DRIVER"] & ev["TYPE"].isin(types)]
    out = {}

    for g in export.genes:
        profiled_all = prof[[(g, t) for t in types]].all(axis=1)
        has = prof.index.isin(drv.loc[drv["GENE"] == g, "SAMPLE_ID"])
        out[f"{g}_DRIVER"] = np.where(has, 1.0, np.where(profiled_all, 0.0, np.nan))

    res = pd.DataFrame(out, index=prof.index)

    for grp, gl in (gene_groups or {}).items():
        cols = [f"{g}_DRIVER" for g in gl if f"{g}_DRIVER" in res.columns]
        absent = [g for g in gl if f"{g}_DRIVER" not in res.columns]
        if absent:
            raise ValueError(f"Genes of group {grp!r} not in the export: {absent}")
        sub = res[cols]
        res[f"{grp}_N_DRIVER"] = sub.sum(axis=1)
        res[f"{grp}_DRIVER"] = np.where(sub.eq(1).any(axis=1), 1.0, np.where(sub.notna().all(axis=1), 0.0, np.nan))

    if fusions:
        sv = ev[ev["TYPE"] == "FUSION"]
        for g in export.genes:
            has = prof.index.isin(sv.loc[sv["GENE"] == g, "SAMPLE_ID"])
            res[f"{g}_SV"] = np.where(has, 1.0, np.where(prof[(g, "FUSION")], 0.0, np.nan))

    return res.reset_index()


def concat_exports(exports: Sequence[AlterationExport]) -> AlterationExport:
    """Combine exports of several studies/queries into one.

    Samples are stacked (IDs must be unique across the exports). Genes missing from one
    export count as not profiled in its samples, so no sample is ever treated as wild type for
    a gene it was never queried for.

    Parameters
    ----------
    exports : sequence of AlterationExport
        e.g. one export per TCGA cohort (STAD, ESCA) read with ``read_alteration_export``.

    Returns
    -------
    AlterationExport
    """

    exports = list(exports)
    if not exports:
        raise ValueError("No exports given")
    samples = pd.concat([e.samples for e in exports], ignore_index=True)
    dup = samples.loc[samples["SAMPLE_ID"].duplicated(), "SAMPLE_ID"]
    if len(dup):
        raise ValueError(f"Sample IDs are not unique across exports, e.g. {dup.iloc[0]!r}")
    events = pd.concat([e.events for e in exports], ignore_index=True)
    profiled = pd.concat([e.profiled for e in exports], axis=0)
    profiled = profiled.fillna(False).astype(bool)
    return AlterationExport(samples, events, profiled)


def alteration_matrix(export: AlterationExport, types: Sequence[str] = ("MUT", "AMP", "HOMDEL"),
                      level: str = "gene", by: str = "SAMPLE_ID", min_samples: int = 0) -> pd.DataFrame:
    """Binary driver matrix, rows x samples (e.g. as the ``mutations`` input of ``pyrea``).

    Parameters
    ----------
    export : AlterationExport
    types : sequence of {'MUT', 'AMP', 'HOMDEL', 'FUSION'}
        Alteration types counted, e.g. ``("MUT",)`` for mutations only or ``("AMP",)``.
    level : {'gene', 'event'}
        ``'gene'``: one row per gene (1 = driver of one of ``types``; 0 only if profiled for
        ALL ``types``). ``'event'``: one row per ``GENE:EVENT`` (see ``driver_event_matrix``).
    by : {'SAMPLE_ID', 'PATIENT_ID'}
        Columns are samples, or patients (several samples of a patient are combined: 1 if any
        is 1, else 0 if any is 0, else missing).
    min_samples : int
        Rows with fewer samples carrying the alteration are dropped (default 0 keeps all genes,
        also those without any driver).

    Returns
    -------
    pd.DataFrame
        1 = present, 0 = absent, NaN = not profiled. Rows in the order of the underlying
        matrix (events: by gene frequency).
    """

    if level == "gene":
        m = driver_gene_matrix(export, types=types, fusions=False).set_index("SAMPLE_ID")
        m.columns = m.columns.str.removesuffix("_DRIVER")
    elif level == "event":
        m = driver_event_matrix(export, types=types)
    else:
        raise ValueError("level must be 'gene' or 'event'")

    m = _aggregate_by(m, export, by)
    return m.loc[m.sum(axis=1) >= min_samples]


def _aggregate_by(m: pd.DataFrame, export: AlterationExport, by: str) -> pd.DataFrame:
    """samples x rows -> rows x (samples or patients); patients: any 1, else any 0, else NaN."""

    if by == "PATIENT_ID":
        patient = export.samples.set_index("SAMPLE_ID")["PATIENT_ID"].reindex(m.index)
        m = m.groupby(patient.to_numpy(), sort=False).max()
        m.index.name = "PATIENT_ID"
    elif by != "SAMPLE_ID":
        raise ValueError("by must be 'SAMPLE_ID' or 'PATIENT_ID'")
    return m.T


_NOT_A_PARTNER = {"intragenic"}


def _parse_fusion(gene: str, event: str):
    """Event listed under ``gene`` -> (5' gene, 3' gene), or None if it is no gene fusion.

    'CLDN18-ARHGAP26 fusion' / 'ERBB2-PSMB3' under CLDN18 / ERBB2; the queried gene is one of
    the two partners and the rest of the text is the other. 'APC-intragenic' and
    non-fusion structural variants give None.
    """

    name = event[: -len(" fusion")] if event.endswith(" fusion") else event
    if name.startswith(gene + "-"):
        five, three = gene, name[len(gene) + 1:]
        other = three
    elif name.endswith("-" + gene):
        five, three = name[: -len(gene) - 1], gene
        other = five
    else:
        return None
    if not other or other in _NOT_A_PARTNER or re.search(r"\s", other):
        return None
    return five, three


def fusion_events(export: AlterationExport) -> pd.DataFrame:
    """Gene fusions of an export, one row per sample, queried gene and event.

    Columns: SAMPLE_ID, QUERIED (gene under which the export lists the event), GENE5, GENE3,
    FUSION ('GENE5-GENE3'). A fusion between two queried genes appears twice (once per
    partner); orientation is kept, 'A-B' and 'B-A' are different fusions. Structural variants
    that are no gene fusions (intragenic, deletions within a transcript) are left out.
    """

    ev = export.events[export.events["TYPE"] == "FUSION"]
    rows = []
    for sid, gene, event in zip(ev["SAMPLE_ID"], ev["GENE"], ev["EVENT"]):
        parsed = _parse_fusion(gene, event)
        if parsed:
            rows.append((sid, gene, *parsed))
    out = pd.DataFrame(rows, columns=["SAMPLE_ID", "QUERIED", "GENE5", "GENE3"])
    out["FUSION"] = out["GENE5"] + "-" + out["GENE3"]
    return out


def fusion_matrix(export: AlterationExport, groups: Optional[Dict[str, str]] = None, others: bool = True,
                  by: str = "SAMPLE_ID", min_samples: int = 0) -> pd.DataFrame:
    """Binary gene-fusion matrix, fusions x samples (same conventions as ``alteration_matrix``).

    Exports never label fusions as drivers, so every fusion is listed. Each fusion appears once
    per sample even if the export reports it under both partners.

    Parameters
    ----------
    export : AlterationExport
    groups : dict of {row name: regex}, optional
        Fusions whose 'GENE5-GENE3' fully matches the regex are merged into one row, e.g.
        ``{"CLDN18-ARHGAP6/26": r"CLDN18-ARHGAP(6|26)"}``. The first matching group wins.
        Groups without any matching fusion are omitted.
    others : bool
        Keep fusions matching no group as their own rows (default); otherwise drop them.
    by : {'SAMPLE_ID', 'PATIENT_ID'}
    min_samples : int
        Rows with fewer samples are dropped.

    Returns
    -------
    pd.DataFrame
        1 = fusion present; 0 = absent and every queried partner gene was profiled for
        structural variants; NaN otherwise. Rows: groups first, then single fusions by
        frequency (descending).
    """

    fe = fusion_events(export)
    index = export.profiled.index
    prof = export.profiled.xs("FUSION", axis=1, level="TYPE")

    def _row(fusion: str) -> Optional[str]:
        for name, pattern in (groups or {}).items():
            if re.fullmatch(pattern, fusion):
                return name
        return fusion if others else None

    fe["ROW"] = fe["FUSION"].map({f: _row(f) for f in fe["FUSION"].unique()})
    fe = fe.dropna(subset=["ROW"])

    counts = fe.drop_duplicates(["SAMPLE_ID", "ROW"]).groupby("ROW").size()
    order = [g for g in (groups or {}) if g in counts.index]
    order += counts.drop(order).sort_values(ascending=False, kind="stable").index.tolist()

    cols = {}
    for row in order:
        sub = fe[fe["ROW"] == row]
        has = index.isin(sub["SAMPLE_ID"])
        evaluable = prof[sorted(set(sub["QUERIED"]))].all(axis=1).to_numpy()
        cols[row] = np.where(has, 1.0, np.where(evaluable, 0.0, np.nan))
    m = pd.DataFrame(cols, index=index)

    m = _aggregate_by(m, export, by)
    return m.loc[m.sum(axis=1) >= min_samples]
