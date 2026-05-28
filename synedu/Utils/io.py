"""I/O helpers for SynEdu datasets."""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any, List

_GZIP_MAGIC = b"\x1f\x8b"


def _is_gzip(path: Path) -> bool:
    """Return True when the file starts with the gzip magic bytes."""
    try:
        with open(path, "rb") as fh:
            return fh.read(2) == _GZIP_MAGIC
    except OSError:
        return False


def load_database(path: Any) -> List[dict]:
    """
    Load a JSON or gzip-compressed JSON database and return a list of records.

    Format is detected from the file content (magic bytes), not the extension,
    so files named ``.json.gz`` that are actually plain JSON are handled correctly.

    :param path: Path to a ``.json`` or ``.json.gz`` file.
    :returns: List of record dicts.
    :rtype: list
    """
    path = Path(path)
    if _is_gzip(path):
        opener = gzip.open
    else:
        opener = open  # type: ignore[assignment]
    with opener(path, "rt", encoding="utf-8") as fh:  # type: ignore[call-overload]
        data = json.load(fh)
    if isinstance(data, list):
        return data
    if isinstance(data, dict):
        return list(data.values())
    return list(data)
