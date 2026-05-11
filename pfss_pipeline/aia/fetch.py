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


def find_local(target_time: Time, wavelengths: list[str],
              data_dir: str | Path) -> dict[str, str]:
    """Match local AIA L1 by minute-prefix glob, closest second wins.

    Mirrors src/AIA_prep_all.ipynb cell 2: pattern is
    `aia.lev1_euv_12s.{YYYY-MM-DDTHHMM}*.{wl}.image_lev1.fits`. AIA cadence
    is 12 s, so a hit minute usually has 5 candidates; we pick the closest
    in seconds. No cross-minute fallback — if the minute is empty, JSOC
    fetch covers it.
    """
    target_dt = datetime.strptime(target_time.iso[:19], "%Y-%m-%d %H:%M:%S")
    minute_str = target_time.iso[:16].replace(":", "").replace(" ", "T")  # e.g. 2022-02-28T0502
    out: dict[str, str] = {}
    for wl in wavelengths:
        pattern = os.path.join(str(data_dir),
                              f"aia.lev1_euv_12s.{minute_str}*.{wl}.image_lev1.fits")
        candidates = sorted(glob.glob(pattern))
        if not candidates:
            continue
        if len(candidates) == 1:
            out[wl] = candidates[0]
        else:
            out[wl] = min(candidates,
                          key=lambda f: abs((_parse_l1_time(f) - target_dt).total_seconds()))
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
