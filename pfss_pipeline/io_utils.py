"""Shared I/O helpers across stages."""
from __future__ import annotations

import glob
import os
import re

from astropy.time import Time
from sunpy.map import Map


def find_closest_map(prep_dir: str, wl: str, target_time: Time, tolerance_minutes: float = 5):
    """Return the prepared AIA map for `wl` whose timestamp is closest to `target_time`.

    Filenames must follow `aia_prep_{wl}A_YYYYMMDD_HHMMSS.fits`.
    Raises if nothing is within `tolerance_minutes`.
    """
    candidates = sorted(glob.glob(os.path.join(prep_dir, f"aia_prep_{wl}A_*.fits")))
    assert candidates, f"No prepared map found for {wl}A in {prep_dir}"

    def parse_time(filepath: str) -> Time:
        m = re.search(r"(\d{8})_(\d{6})\.fits$", filepath)
        assert m, f"Cannot parse time from filename: {filepath}"
        date_str, time_str = m.groups()
        iso = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}T{time_str[:2]}:{time_str[2:4]}:{time_str[4:6]}"
        return Time(iso)

    times = [parse_time(c) for c in candidates]
    deltas = [(t - target_time).to("min").value for t in times]
    closest_idx = min(range(len(deltas)), key=lambda i: abs(deltas[i]))
    if abs(deltas[closest_idx]) > tolerance_minutes:
        raise FileNotFoundError(
            f"No prepared {wl}A map within {tolerance_minutes} min of {target_time.iso}; "
            f"closest = {times[closest_idx].iso} (Δ={deltas[closest_idx]:.2f} min)"
        )
    return Map(candidates[closest_idx])
