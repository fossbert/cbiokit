"""fusion_events and fusion_matrix."""

import numpy as np
import pandas as pd
import pytest

from cbiokit import fusion_events, fusion_matrix, read_alteration_export
from test_alterations import export  # noqa: F401  (fixture)

NA, NP = "no alteration", "not profiled"

# sample, patient, CLDN18 FUSION, ARHGAP26 FUSION (both genes queried, FUSION column only)
ROWS = [
    ("S1", "P1", "CLDN18-ARHGAP26 fusion", "CLDN18-ARHGAP26 fusion"),            # listed under both partners
    ("S2", "P2", "CLDN18-ARHGAP6 fusion", NA),                                      # ARHGAP6 not queried
    ("S3", "P3", "ARHGAP26-CLDN18 fusion", "ARHGAP26-CLDN18 fusion"),              # reciprocal orientation
    ("S4", "P4", NA, NA),
    ("S5", "P5", NP, NA),                                                           # CLDN18 not profiled
    ("S6", "P1", "CLDN18-intragenic, LYZ-CLDN18", NA),                              # 2nd sample of P1
]


@pytest.fixture
def fusion_export(tmp_path):
    recs = []
    for s, p, c, a in ROWS:
        rec = {"Study ID": "demo", "Sample ID": s, "Patient ID": p, "Altered": 0, "CLDN18": c, "ARHGAP26": a}
        rec.update({"CLDN18: MUT": NA, "CLDN18: FUSION": c, "ARHGAP26: MUT": NA, "ARHGAP26: FUSION": a})
        recs.append(rec)
    f = tmp_path / "fusions.tsv"
    pd.DataFrame(recs).to_csv(f, sep="\t", index=False)
    return read_alteration_export(f)


def test_fusion_events_parsing(fusion_export):
    fe = fusion_events(fusion_export)
    assert set(fe["FUSION"]) == {"CLDN18-ARHGAP26", "CLDN18-ARHGAP6", "ARHGAP26-CLDN18", "LYZ-CLDN18"}
    assert not fe["FUSION"].str.contains("intragenic").any()          # not a gene fusion
    row = fe[fe["FUSION"] == "LYZ-CLDN18"].iloc[0]
    assert (row["GENE5"], row["GENE3"], row["QUERIED"]) == ("LYZ", "CLDN18", "CLDN18")


def test_groups_merge_variants_and_dedupe(fusion_export):
    m = fusion_matrix(fusion_export, groups={"CLDN18-ARHGAP6/26": r"CLDN18-ARHGAP(6|26)"})
    assert m.index.tolist()[0] == "CLDN18-ARHGAP6/26"                 # groups first
    assert m.loc["CLDN18-ARHGAP6/26"].tolist()[:4] == [1, 1, 0, 0]
    assert m.loc["ARHGAP26-CLDN18", ["S1", "S3"]].tolist() == [0, 1]  # reciprocal is its own fusion
    assert m.loc["CLDN18-ARHGAP6/26"].sum() == 2                      # S1 not double counted


def test_profiling_gives_nan(fusion_export):
    m = fusion_matrix(fusion_export, groups={"G": r"CLDN18-ARHGAP(6|26)"})
    assert np.isnan(m.loc["G", "S5"])                                 # CLDN18 FUSION not profiled
    assert m.loc["G", "S4"] == 0
    # a fusion listed under both partners needs both profiled; S5 lacks CLDN18 -> NaN, not 0
    assert np.isnan(fusion_matrix(fusion_export).loc["CLDN18-ARHGAP26", "S5"])


def test_others_min_samples_and_patient(fusion_export):
    only = fusion_matrix(fusion_export, groups={"G": r"CLDN18-ARHGAP(6|26)"}, others=False)
    assert only.index.tolist() == ["G"]
    grp = {"G": r"CLDN18-ARHGAP(6|26)"}
    assert fusion_matrix(fusion_export, groups=grp, min_samples=2).index.tolist() == ["G"]
    p = fusion_matrix(fusion_export, by="PATIENT_ID", groups={"G": r"CLDN18-ARHGAP(6|26)"})
    assert list(p.columns) == ["P1", "P2", "P3", "P4", "P5"]
    assert p.loc["G", "P1"] == 1                                      # S1 (1) + S6 (0) -> 1
    assert p.loc["LYZ-CLDN18", "P1"] == 1


def test_export_without_fusions(export):  # noqa: F811
    m = fusion_matrix(export)
    assert "KRAS-CDH1" in m.index                                     # the only gene fusion in the fixture
    assert not any("intragenic" in r or "Deletion" in r for r in m.index)
