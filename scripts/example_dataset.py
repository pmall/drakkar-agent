"""Template for a dataset script. Run: uv run python scripts/<name>.py"""

from drakkar.db import export

VALID = """
    state = 'curated'
    AND is_obsolete1 IS FALSE
    AND is_obsolete2 IS FALSE
    AND deleted_at IS NULL
"""

df = export(
    f"""
    SELECT accession1, name1, accession2, name2, start2, stop2, taxon2, pmid, method
    FROM dataset
    WHERE {VALID} AND type = 'vh'
    """,
    "data/example_vh.tsv",
)

print(df.shape)
