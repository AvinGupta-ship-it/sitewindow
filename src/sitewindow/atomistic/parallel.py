"""Minimal MPI awareness for workflows that may run under ``mpirun ... gpaw python``.

Only rank 0 writes files; every rank computes. When GPAW is not imported the process is
treated as a single rank.
"""

from __future__ import annotations

import sys
from typing import Any


def _world() -> Any | None:
    if "gpaw" in sys.modules or "_gpaw" in sys.modules:
        try:
            from gpaw.mpi import world

            return world
        except Exception:  # pragma: no cover - gpaw import failure
            return None
    return None


def rank() -> int:
    w = _world()
    return 0 if w is None else int(w.rank)


def size() -> int:
    w = _world()
    return 1 if w is None else int(w.size)


def is_master() -> bool:
    return rank() == 0


def barrier() -> None:
    w = _world()
    if w is not None:
        w.barrier()


def broadcast_array(values: Any) -> Any:
    """Return rank 0's copy of a float array on every rank (identity without MPI)."""
    import numpy as np

    w = _world()
    arr = np.ascontiguousarray(np.asarray(values, dtype=float))
    if w is None or w.size == 1:
        return arr
    w.broadcast(arr, 0)
    return arr


__all__ = ["rank", "size", "is_master", "barrier", "broadcast_array"]
