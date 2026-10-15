"""Co-occurring alterations: combinations of alterations per patient.

:func:`calc_combo_muts` and :func:`make_combo_df` expect a binary table with one row per patient and one
column per alteration (0/1), i.e. the transpose of the matrices from :mod:`cbiokit.alterations`
(``alteration_matrix(...).T``). A column MultiIndex is flattened to its second level.
:func:`joint_status` takes two alterations as Series (e.g. rows of ``alteration_matrix``) and labels every patient
by which of the two it carries.
"""

from collections import defaultdict
from itertools import combinations
from typing import Optional, Sequence

import numpy as np
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


def joint_status(a: pd.Series, b: pd.Series, names: Optional[Sequence[str]] = None,
                 neither: str = "neither") -> pd.Series:
    """Label every sample by which of two alterations it carries (2 x 2 design).

    Parameters
    ----------
    a, b : pd.Series
        Alteration per sample: 1 (or > 0) = present, 0 = absent, NaN = not profiled. Aligned by index; a sample
        missing from one of them is NaN.
    names : sequence of 3 str, optional
        Labels for "a only", "b only" and "both". Default: built from the Series names
        (``"KRAS only"``, ``"TP53 only"``, ``"KRAS+TP53"``; ``"A"`` / ``"B"`` for unnamed Series).
    neither : str
        Label of samples with neither alteration.

    Returns
    -------
    pd.Series
        Ordered categorical with the categories ``[neither, a only, b only, both]`` (all present even if empty), NaN
        where either alteration is not profiled. Use ``.value_counts()`` or ``pd.crosstab(a, b)`` for the cell sizes.

    Examples
    --------
    >>> mat = cbk.alteration_matrix(ex, types=("MUT",), by="PATIENT_ID")
    >>> status = cbk.joint_status(mat.loc["RHOA"], mat.loc["TP53"])
    >>> status.value_counts()          # neither / RHOA only / TP53 only / RHOA+TP53
    """

    an = "A" if a.name is None else str(a.name)
    bn = "B" if b.name is None else str(b.name)
    labels = list(names) if names is not None else [f"{an} only", f"{bn} only", f"{an}+{bn}"]
    if len(labels) != 3:
        raise ValueError("names needs three labels: a only, b only, both")
    cats = [neither] + labels
    if len(set(cats)) != 4:
        raise ValueError(f"The four labels must differ: {cats}")
    idx = a.index.union(b.index, sort=False)
    x, y = a.reindex(idx).astype(float), b.reindex(idx).astype(float)
    xp, yp = x.gt(0), y.gt(0)
    out = pd.Series(np.select([xp & yp, xp & ~yp, ~xp & yp, ~xp & ~yp], [labels[2], labels[0], labels[1], neither],
                              default=""), index=idx, dtype=object)
    out[x.isna() | y.isna()] = np.nan
    return pd.Series(pd.Categorical(out, categories=cats, ordered=True), index=idx, name=f"{an}/{bn}")
