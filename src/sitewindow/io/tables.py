"""Parquet tables carrying the schema version in the file metadata; SHA-256 manifests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from sitewindow.schema.models import SCHEMA_VERSION

METADATA_KEY = b"sitewindow_schema_version"


def write_table(
    frame: pd.DataFrame, path: Path, extra_metadata: dict[str, str] | None = None
) -> Path:
    """Write ``frame`` to Parquet with ``sitewindow_schema_version`` in the file metadata."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    meta = dict(table.schema.metadata or {})
    meta[METADATA_KEY] = SCHEMA_VERSION.encode()
    for k, v in (extra_metadata or {}).items():
        meta[k.encode()] = v.encode()
    pq.write_table(table.replace_schema_metadata(meta), path)
    return path


def read_table(path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read a Parquet table and its string metadata (schema version included)."""
    table = pq.read_table(Path(path))
    meta = {
        k.decode(): v.decode()
        for k, v in (table.schema.metadata or {}).items()
        if not k.startswith(b"pandas")
    }
    return table.to_pandas(), meta


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest(
    paths: list[Path], manifest_path: Path, extra: dict[str, Any] | None = None
) -> Path:
    """Write ``MANIFEST``-style JSON mapping relative file names to SHA-256 hashes."""
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    entries = {Path(p).name: sha256_file(Path(p)) for p in paths}
    payload: dict[str, Any] = {"files": entries}
    if extra:
        payload.update(extra)
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return manifest_path


__all__ = ["write_table", "read_table", "sha256_file", "write_manifest", "METADATA_KEY"]
