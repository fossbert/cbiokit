"""concat_exports and alteration_matrix."""

import numpy as np
import pandas as pd
import pytest

from cbiokit import alteration_matrix, concat_exports, read_alteration_export
from test_alterations import export  # noqa: F401  (fixture)

NA, NP = "no alteration", "not profiled"


def _write(tmp_path, name, rows, genes):
    """rows: (sample, patient, {gene: (MUT, AMP)}) -> one export file with MUT and AMP columns."""
    recs = []
    for s, p, vals in rows:
        rec = {"Study ID": name, "Sample ID": s, "Patient ID": p, "Altered": 0}
        for g in genes:
            mut, amp = vals[g]
            rec[g] = ", ".join(c for c in (mut, amp) if c not in (NA, NP)) or NA
            rec[f"{g}: MUT"], rec[f"{g}: AMP"] = mut, amp
        recs.append(rec)
    f = tmp_path / f"{name}.tsv"
    pd.DataFrame(recs).to_csv(f, sep="\t", index=False)
    return read_alteration_export(f)


@pytest.fixture
def two_exports(tmp_path):
    a = _write(tmp_path, "stad", [("A1", "PA", {"KRAS": ("G12D (driver)", NA), "ERBB2": (NA, "AMP (driver)")}),
                                  ("A2", "PB", {"KRAS": (NA, NA), "ERBB2": (NP, NA)})], ["KRAS", "ERBB2"])
    b = _write(tmp_path, "eac", [("B1", "PC", {"KRAS": ("G12V (driver)", NA), "TP53": ("R175H (driver)", NA)})],
               ["KRAS", "TP53"])
    return a, b


def test_concat_profiling_of_missing_genes(two_exports):
    ex = concat_exports(two_exports)
    assert ex.samples["SAMPLE_ID"].tolist() == ["A1", "A2", "B1"]
    assert sorted(ex.genes) == ["ERBB2", "KRAS", "TP53"]
    m = alteration_matrix(ex, types=("MUT",))
    assert m.loc["KRAS"].tolist() == [1, 0, 1]
    assert m.loc["ERBB2", ["A1", "B1"]].isna().tolist() == [False, True]  # B1: ERBB2 never queried
    assert m.loc["TP53"].isna().tolist() == [True, True, False]           # TP53 not queried in stad
    assert len(ex.events) == sum(len(e.events) for e in two_exports)


def test_concat_rejects_duplicate_samples(two_exports):
    a, _ = two_exports
    with pytest.raises(ValueError, match="not unique"):
        concat_exports([a, a])
    with pytest.raises(ValueError, match="No exports"):
        concat_exports([])


def test_gene_matrix_orientation_and_types(export):  # noqa: F811
    m = alteration_matrix(export)
    assert list(m.columns) == ["S1", "S2", "S3", "S4", "S5"]
    assert m.loc["KRAS", ["S1", "S2", "S4"]].tolist() == [1, 1, 1]
    assert m.loc["KRAS", ["S3", "S5"]].isna().all()
    assert not any(str(i).endswith("_DRIVER") for i in m.index)
    amp = alteration_matrix(export, types=("AMP",))
    assert amp.loc["ERBB2"].tolist()[:3] == [1, 0, 0]           # S3: CNA profiled -> 0
    assert amp.loc["KRAS"].isna().tolist() == [False] * 4 + [True]


def test_event_matrix_and_min_samples(export):  # noqa: F811
    m = alteration_matrix(export, level="event")
    assert m.index.tolist() == ["KRAS:G12D", "KRAS:G12V", "TP53:HOMDEL", "TP53:MUTATED", "TP53:R175H", "ERBB2:AMP"]
    assert alteration_matrix(export, level="event", min_samples=2).index.tolist() == ["KRAS:G12D"]


def test_by_patient_combines_samples(export):  # noqa: F811
    m = alteration_matrix(export, by="PATIENT_ID", types=("MUT",))
    assert list(m.columns) == ["P1", "P2", "P3", "P4"]          # S1+S2 -> P1
    assert m.loc["ERBB2", "P1"] == 0 and m.loc["TP53", "P1"] == 1  # ERBB2 driver is an AMP; S310F is no driver
    assert m.loc["KRAS", "P1"] == 1
    assert m.loc["KRAS", ["P2", "P4"]].isna().all()


def test_invalid_arguments(export):  # noqa: F811
    with pytest.raises(ValueError):
        alteration_matrix(export, level="exon")
    with pytest.raises(ValueError):
        alteration_matrix(export, by="STUDY_ID")
