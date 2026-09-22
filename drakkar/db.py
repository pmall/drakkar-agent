"""Connection + query helpers for the drakkar PPI database."""

from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import psycopg
from dotenv import load_dotenv

ROOT = Path(__file__).resolve().parents[1]


def dsn() -> str:
    load_dotenv(ROOT / ".env")
    return (
        f"host={os.environ['POSTGRES_HOST']} "
        f"port={os.environ['POSTGRES_PORT']} "
        f"dbname={os.environ['POSTGRES_DB']} "
        f"user={os.environ['POSTGRES_USER']} "
        f"password={os.environ['POSTGRES_PASSWORD']}"
    )


def connect() -> psycopg.Connection:
    return psycopg.connect(dsn())


def fetch(sql: str, params: tuple | dict | None = None) -> pl.DataFrame:
    """Run a SELECT and return the rows as a polars DataFrame."""
    with connect() as conn, conn.cursor() as cur:
        cur.execute(sql, params)
        cols = [d.name for d in cur.description]
        rows = cur.fetchall()
    return pl.DataFrame(rows, schema=cols, orient="row", infer_schema_length=None)


def stream(sql: str, path: str | Path) -> None:
    """Stream a SELECT straight to .csv/.tsv via COPY, without holding rows.

    Same output as ``export`` (RFC 4180 quoting, header, empty field for
    NULL), but the result set never lands in memory on either side -- use it
    instead of ``export`` for datasets large enough to be worth not holding
    as a DataFrame.
    """
    path = Path(path)
    delimiters = {".csv": ",", ".tsv": "\t"}
    if path.suffix not in delimiters:
        raise ValueError(f"unsupported output format: {path.suffix}")
    path.parent.mkdir(parents=True, exist_ok=True)
    copy = (
        f"COPY ({sql}) TO STDOUT "
        f"WITH (FORMAT csv, DELIMITER '{delimiters[path.suffix]}', HEADER)"
    )
    with (
        connect() as conn,
        conn.cursor() as cur,
        path.open("wb") as fh,
        cur.copy(copy) as reader,
    ):
        for chunk in reader:
            fh.write(chunk)


def export(
    sql: str, path: str | Path, params: tuple | dict | None = None
) -> pl.DataFrame:
    """Run a SELECT, write it to .csv/.tsv/.parquet, and return the DataFrame."""
    df = fetch(sql, params)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    match path.suffix:
        case ".csv":
            df.write_csv(path)
        case ".tsv":
            df.write_csv(path, separator="\t")
        case ".parquet":
            df.write_parquet(path)
        case other:
            raise ValueError(f"unsupported output format: {other}")
    return df
