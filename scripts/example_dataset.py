"""Template for a dataset script. Run: uv run python scripts/<name>.py"""

from drakkar.db import dataset_path, export
from drakkar.descriptions import valid

df = export(
    f"""
    SELECT accession1, name1, accession2, name2, start2, stop2, taxon2, pmid, method
    FROM dataset
    WHERE {valid()} AND type = 'vh'
    """,
    dataset_path("example_vh.tsv"),
)

print(df.shape)
