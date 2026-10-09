"""variant_class and variant_class_matrix."""

import numpy as np
import pytest

from cbiokit import alteration_matrix, variant_class, variant_class_matrix
from test_alterations import export  # noqa: F401  (fixture)


@pytest.mark.parametrize("event, expected", [
    ("R213*", "truncating"), ("S70Pfs*13", "truncating"), ("K2Rfs*", "truncating"),
    ("X229_splice", "truncating"), ("T515_F516ins*", "truncating"), ("X1234_splice", "truncating"),
    ("E746_A750del", "inframe"), ("L220_D221delinsP", "inframe"), ("A123dup", "inframe"), ("Q61delinsHK", "inframe"),
    ("G12D", "missense"), ("R175H", "missense"),
    ("M1?", "other"), ("MUTATED", "other"), ("AMP", "other"),
])
def test_variant_class(event, expected):
    assert variant_class(event) == expected


def test_class_matrix_conventions(export):  # noqa: F811
    m = variant_class_matrix(export)                       # fixture: KRAS G12D/G12V, TP53 R175H, P322Hfs*23, MUTATED, ERBB2 S310F
    assert set(m.index) >= {"KRAS:missense", "TP53:missense", "TP53:truncating", "TP53:other", "ERBB2:missense"}
    assert m.loc["TP53:truncating"].tolist()[:2] == [1, 0]                 # S1 has P322Hfs*23 (a non-driver variant)
    assert m.loc["TP53:missense", "S1"] == 1 and m.loc["TP53:other", "S4"] == 1
    assert m.loc["ERBB2:missense", "S2"] == 1                            # S310F is no driver but counts by default
    assert m.loc["KRAS:missense"].isna().tolist() == [False, False, True, False, True]   # not profiled for MUT
    drv = variant_class_matrix(export, drivers_only=True)
    assert "TP53:truncating" not in drv.index and "ERBB2:missense" not in drv.index
    pat = variant_class_matrix(export, by="PATIENT_ID")
    assert list(pat.columns) == ["P1", "P2", "P3", "P4"]
    assert variant_class_matrix(export, min_samples=2).shape[0] < m.shape[0]
    # same samples carry a gene mutation as in the plain gene-level matrix
    gene = alteration_matrix(export, types=("MUT",), drivers_only=False)
    cls = m.loc[[i for i in m.index if i.startswith("TP53:")]].max(axis=0)
    assert (cls.fillna(-1) == gene.loc["TP53"].fillna(-1)).all()
