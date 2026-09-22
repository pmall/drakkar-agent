"""Append-only run log for dataset scripts.

Run dates and counts live here (committed as ``RUNS.md``), not in script
docstrings: docstrings describe the dataset, this log records each run.
``data/`` is gitignored, so without this log a run's numbers would be lost.
"""

from __future__ import annotations

import os
from collections.abc import Mapping
from datetime import UTC, datetime
from pathlib import Path

from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]
RUNS_MD = ROOT / "RUNS.md"

_HEADER = """# Dataset runs

Append-only log of dataset exports. Each run appends one section below via
`drakkar.runs.log_run` -- never edit past entries.
"""


def _format(value: str | int) -> str:
    return f"{value:,}" if isinstance(value, int) else str(value)


def _db_name() -> str:
    load_dotenv(ROOT / ".env")
    return os.environ.get("POSTGRES_DB", "unknown")


def log_run(
    script: str,
    output: str,
    stats: Mapping[str, str | int],
    date: str | None = None,
) -> str:
    """Append a run entry to RUNS.md and return the entry text.

    ``date`` overrides today (UTC) -- only for backfilling historical runs.
    """
    day = date or datetime.now(UTC).date().isoformat()
    lines = [
        f"## {day} -- {script}",
        "",
        f"- output: `{output}`",
        f"- database: {_db_name()}",
    ]
    lines += [f"- {key}: {_format(value)}" for key, value in stats.items()]
    entry = "\n".join(lines) + "\n"
    if not RUNS_MD.exists():
        RUNS_MD.write_text(_HEADER + "\n" + entry)
    else:
        with RUNS_MD.open("a") as fh:
            fh.write("\n" + entry)
    return entry
