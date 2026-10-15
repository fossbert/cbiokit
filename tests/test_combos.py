"""calc_combo_muts and make_combo_df."""

import pandas as pd

import numpy as np
import pytest

from cbiokit import calc_combo_muts, joint_status, make_combo_df

DF = pd.DataFrame({"A": [1, 1, 1, 0], "B": [1, 1, 0, 1], "C": [1, 0, 0, 1]}, index=["p1", "p2", "p3", "p4"])


def test_counts_and_sorting():
    res = calc_combo_muts(DF, min_occurence=0)
    assert res.to_dict() == {"A_B": 2, "B_C": 2, "A_C": 1, "A_B_C": 1}
    assert list(res.values) == sorted(res.values, reverse=True)


def test_min_occurence_is_strict():
    assert calc_combo_muts(DF, min_occurence=1).to_dict() == {"A_B": 2, "B_C": 2}
    assert calc_combo_muts(DF, min_occurence=2).empty


def test_max_size_and_multiindex():
    assert "A_B_C" not in calc_combo_muts(DF, min_occurence=0, max_size=2).index
    mi = DF.copy()
    mi.columns = pd.MultiIndex.from_product([["x"], DF.columns])
    assert calc_combo_muts(mi, min_occurence=1).to_dict() == {"A_B": 2, "B_C": 2}


def test_make_combo_df():
    out = make_combo_df(DF, ["A_B", "A_B_C"])
    assert out["A_B"].tolist() == [1, 1, 0, 0]
    assert out["A_B_C"].tolist() == [1, 0, 0, 0]
    assert make_combo_df(DF, []).shape == (4, 0)


def test_custom_separator():
    d = DF.rename(columns={"A": "X_1"})
    res = calc_combo_muts(d, min_occurence=1, sep="|")
    assert "B|X_1" in res.index
    assert make_combo_df(d, res.index, splitter="|")["B|X_1"].tolist() == [1, 1, 0, 0]


def test_joint_status_four_cells_and_nan():
    a = pd.Series([1, 1, 0, 0, np.nan, 1], index=list("pqrstu"), name="KRAS")
    b = pd.Series([1, 0, 1, 0, 1, np.nan], index=list("pqrstu"), name="TP53")
    st = joint_status(a, b)
    assert list(st.cat.categories) == ["neither", "KRAS only", "TP53 only", "KRAS+TP53"] and st.cat.ordered
    assert st.loc["p"] == "KRAS+TP53" and st.loc["q"] == "KRAS only" and st.loc["r"] == "TP53 only" and st.loc["s"] == "neither"
    assert st.isna().loc[["t", "u"]].all() and st.notna().sum() == 4          # not profiled for one of them
    assert st.value_counts().to_dict() == {"neither": 1, "KRAS only": 1, "TP53 only": 1, "KRAS+TP53": 1}


def test_joint_status_names_empty_cells_and_alignment():
    a = pd.Series([1, 0, 0], index=["x", "y", "z"])
    b = pd.Series([0, 0], index=["y", "z"])                                   # sample x missing from b -> NaN
    st = joint_status(a, b, names=("a!", "b!", "ab"), neither="none")
    assert list(st.cat.categories) == ["none", "a!", "b!", "ab"]
    assert st.loc["y"] == "none" and pd.isna(st.loc["x"])
    assert st.value_counts().to_dict() == {"none": 2, "a!": 0, "b!": 0, "ab": 0}          # empty cells are kept
    unnamed = joint_status(pd.Series([1.0]), pd.Series([1.0]))
    assert unnamed.iloc[0] == "A+B"
    with pytest.raises(ValueError, match="three labels"):
        joint_status(a, b, names=("only", "two"))
    with pytest.raises(ValueError, match="must differ"):
        joint_status(a, b, names=("x", "x", "y"))


def test_joint_status_counts_match_crosstab():
    rng = np.random.default_rng(0)
    a = pd.Series(rng.integers(0, 2, 200).astype(float), name="A")
    b = pd.Series(rng.integers(0, 2, 200).astype(float), name="B")
    st = joint_status(a, b)
    tab = pd.crosstab(a, b)
    assert st.value_counts()["A+B"] == tab.loc[1, 1] and st.value_counts()["neither"] == tab.loc[0, 0]
