"""Co-occurring alterations: combinations of alterations per patient.

Both functions expect a binary table with one row per patient and one column per alteration
(0/1), i.e. the transpose of the matrices from :mod:`cbiokit.alterations`
(``alteration_matrix(...).T``). A column MultiIndex is flattened to its second level.
"""

from collections import defaultdict
from itertools import combinations
from typing import Optional, Sequence

import pandas as pd


def _flat_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if isinstance(df.columns, pd.MultiIndex):
        df = df.set_axis(df.columns.get_level_values(1), axis=1)
    return df


def calc_combo_muts(mut_df: pd.DataFrame, min_occurence: int = 10, max_size: Optional[int] = None,
                    sep: str = "_") -> pd.Series:
    """Count every combination of alterations (size >= 2) that occurs together in a row.

    Parameters
    ----------
    mut_df : pd.DataFrame
        Rows = patients, columns = alterations, values 0/1 (NaN counts as absent).
    min_occurence : int
        Only combinations occurring in *more than* this many rows are returned.
    max_size : int, optional
        Largest combination size; default all. Rows with many alterations have exponentially
        many combinations, so set this (e.g. 3) for dense tables.
    sep : str
        Joins the sorted alteration names into the combination id (``'KRAS_TP53'``).

    Returns
    -------
    pd.Series
        Count per combination id, sorted descending.
    """

    md = _flat_columns(mut_df)
    res = defaultdict(int)
    for _, row in md.iterrows():
        mut_list = row.index[row > 0].to_list()
        for i in range(2, (len(mut_list) if max_size is None else min(max_size, len(mut_list))) + 1):
            for c in combinations(mut_list, r=i):
                res[sep.join(sorted(map(str, c)))] += 1

    res_sorted = {k: v for k, v in sorted(res.items(), key=lambda item: item[1], reverse=True) if v > min_occurence}
    return pd.Series(res_sorted.values(), index=res_sorted.keys(), dtype=int)


def make_combo_df(mut_df: pd.DataFrame, combos: Sequence, splitter: str = "_") -> pd.DataFrame:
    """Binary table of combinations: 1 if a row has all alterations of the combination.

    Parameters
    ----------
    mut_df : pd.DataFrame
        Rows = patients, columns = alterations, values 0/1.
    combos : sequence of str
        Combination ids as returned by :func:`calc_combo_muts` (the index of its result).
    splitter : str
        Separates the alteration names in an id; use the ``sep`` of :func:`calc_combo_muts`.

    Returns
    -------
    pd.DataFrame
        Same rows as ``mut_df``, one column per combination.
    """

    md = _flat_columns(mut_df)
    comb = {}
    for c in list(combos):
        idx = c.split(splitter)
        comb[c] = (md[idx].sum(axis=1) == len(idx)).astype(int)
    if not comb:
        return pd.DataFrame(index=md.index)
    return pd.concat(comb, axis=1)
