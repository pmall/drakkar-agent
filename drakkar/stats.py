"""Counts of a produced dataset, printed for the person (or agent) who logs it.

Scripts log nothing: how a dataset was made is more than one script's doing, and
the agent writes it up in the database folder's `DATASETS.md`. `show` prints what
a script knows about the data it just wrote, in the form its entries take.
"""

from __future__ import annotations

from collections.abc import Mapping

from drakkar.db import database


def _format(value: str | int) -> str:
    return f"{value:,}" if isinstance(value, int) else str(value)


def show(stats: Mapping[str, str | int]) -> None:
    """Print the database name, then one `- key: value` line per count."""
    print(f"- database: {database()}")
    for key, value in stats.items():
        print(f"- {key}: {_format(value)}")
