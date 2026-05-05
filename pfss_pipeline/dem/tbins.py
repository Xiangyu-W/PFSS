"""Resolve DEM logT range and bin edges from region_type (or custom)."""
from __future__ import annotations

import numpy as np

# Region defaults:
#   CH (coronal hole): logT 5.5 – 6.5
#   AR (active region): logT 5.5 – 6.7 (user can extend to 7.0 via region_type='custom')
_DEFAULTS = {
    "CH": (5.5, 6.5),
    "AR": (5.5, 6.7),
}

LOGT_FLOOR = 5.5  # never go below this; keeps the inversion well-conditioned


def resolve_grid(region_type: str, logT_range: list | None, bin_width: float = 0.1) -> tuple:
    """Return (T_bin_edges, n_T_bins, logT_range_resolved)."""
    if region_type in _DEFAULTS:
        lo, hi = _DEFAULTS[region_type]
    elif region_type == "custom":
        if not logT_range or len(logT_range) != 2:
            raise ValueError("region_type='custom' requires a 2-element dem.logT_range")
        lo, hi = float(logT_range[0]), float(logT_range[1])
    else:
        raise ValueError(f"unknown region_type {region_type!r}; expected 'CH'|'AR'|'custom'")

    if lo < LOGT_FLOOR:
        raise ValueError(f"logT lower bound {lo} is below floor {LOGT_FLOOR}")
    if hi <= lo:
        raise ValueError(f"logT range invalid: lo={lo} hi={hi}")

    n_bins = int(round((hi - lo) / bin_width))
    if n_bins < 1:
        raise ValueError(f"bin_width {bin_width} produces 0 bins for range [{lo}, {hi}]")

    T_bin_edges = 10 ** np.linspace(lo, hi, num=n_bins + 1)
    return T_bin_edges, n_bins, (lo, hi)
