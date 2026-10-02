"""Rolling strain-level viral taxa up the NCBI taxonomy, shared by dataset scripts.

`taxon2` is strain-level, so grouping by virus needs the nested-set bounds of
the taxon table rather than string matching.
"""

from __future__ import annotations

from typing import Literal

import polars as pl

from drakkar.db import fetch

type Rank = Literal["species", "genus", "family"]


def roll_up(df: pl.DataFrame, rank: Rank = "species") -> pl.DataFrame:
    """Add a `virus` column: the name of each row's ancestor taxon at `rank`.

    `df` carries the `taxon2` and `left_value2` columns of the `dataset` view.
    A taxon with no ancestor at that rank keeps its own `taxon2` as `virus`.
    The rollup is an interval match done here: a range join against the
    unindexed view never returns.
    """
    nodes = fetch(
        """
        SELECT t.left_value, t.right_value, n.name AS virus
        FROM taxon t
        JOIN taxon_name n
          ON n.taxon_id = t.taxon_id AND n.name_class = 'scientific name'
        WHERE t.node_rank = %(rank)s
        """,
        {"rank": rank},
    )
    rows = df.with_row_index("_row")
    matched = rows.join_where(
        nodes,
        pl.col("left_value2") >= pl.col("left_value"),
        pl.col("left_value2") <= pl.col("right_value"),
    ).select("_row", "virus")
    return (
        rows.join(matched, on="_row", how="left")
        .with_columns(virus=pl.col("virus").fill_null(pl.col("taxon2")))
        .drop("_row")
    )
