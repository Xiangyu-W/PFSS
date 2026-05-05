"""AIA Level 1 file resolution: local glob first, JSOC fallback."""
from __future__ import annotations

import glob
import logging
import os
from datetime import datetime
from functools import reduce
from pathlib import Path

from astropy.time import Time, TimeDelta
from sunpy.net import Fido, attrs as a

log = logging.getLogger(__name__)


def _parse_l1_time(f: str) -> datetime:
    ts = os.path.basename(f).split(".")[2]
    return datetime.strptime(ts, "%Y-%m-%dT%H%M%SZ")


def find_local(target_time: Time, wavelengths: list[str], data_dir: str | Path,
              tolerance_seconds: float = 120.0) -> dict[str, str]:
    """Glob all aia.lev1_euv_12s files for each wavelength; pick the one closest to target_time.

    Uses a ±tolerance_seconds window to catch exposures that straddle a minute boundary
    (e.g. raw exposed at 13:59:59 for a 14:00 target).
    """
    target_dt = datetime.strptime(target_time.iso[:19], "%Y-%m-%d %H:%M:%S")
    out: dict[str, str] = {}
    for wl in wavelengths:
        candidates = sorted(glob.glob(os.path.join(str(data_dir),
                                                  f"aia.lev1_euv_12s.*.{wl}.image_lev1.fits")))
        if not candidates:
            continue
        best = min(candidates, key=lambda f: abs((_parse_l1_time(f) - target_dt).total_seconds()))
        delta = abs((_parse_l1_time(best) - target_dt).total_seconds())
        if delta <= tolerance_seconds:
            out[wl] = best
    return out


def fetch_jsoc(target_time: Time, missing_wl: list[str], data_dir: str | Path,
               jsoc_notify: str) -> dict[str, str]:
    """Use Fido + JSOC to download the missing wavelengths."""
    if not missing_wl:
        return {}
    if not jsoc_notify:
        raise ValueError("jsoc_notify required to download from JSOC")
    log.info("downloading from JSOC: %s", missing_wl)
    primekeys = [a.jsoc.PrimeKey("WAVELNTH", wl) for wl in missing_wl]
    q = Fido.search(
        a.Time(target_time, target_time + TimeDelta(1, format="sec")),
        a.jsoc.Series("aia.lev1_euv_12s"),
        reduce(lambda x, y: x | y, primekeys),
        a.jsoc.Segment("image"),
        a.jsoc.Notify(jsoc_notify),
    )
    fetched = Fido.fetch(q, path=str(Path(data_dir) / "{file}"))
    out: dict[str, str] = {}
    for f in fetched:
        for wl in missing_wl:
            if f".{wl}.image_lev1" in str(f):
                out[wl] = str(f)
                break
    return out


def find_local_or_fetch(target_time: Time, wavelengths: list[str], data_dir: str | Path,
                       jsoc_notify: str | None) -> dict[str, str]:
    """Resolve all wavelengths to local FITS paths. Triggers JSOC fetch for any missing."""
    local = find_local(target_time, wavelengths, data_dir)
    missing = [wl for wl in wavelengths if wl not in local]
    if missing:
        log.info("missing locally: %s; querying JSOC", missing)
        local.update(fetch_jsoc(target_time, missing, data_dir, jsoc_notify or ""))

    still_missing = [wl for wl in wavelengths if wl not in local]
    if still_missing:
        raise RuntimeError(f"could not resolve AIA L1 for wavelengths: {still_missing}")
    return local
