"""calc_combo_muts and make_combo_df."""

import pandas as pd

from cbiokit import calc_combo_muts, make_combo_df

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
