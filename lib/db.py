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


def export(sql: str, path: str | Path, params: tuple | dict | None = None) -> pl.DataFrame:
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
