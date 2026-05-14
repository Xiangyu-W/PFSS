"""AIA Level 1 file resolution: local glob first, JSOC fallback.

Two JSOC backends supported (selected by cfg.aia.fetch_backend):
  - "drms-as-is" (default): drms `protocol="as-is"` + `method="url_quick"` —
    skips JSOC's FITS-rewrite export queue. ~2-10x faster than Fido depending
    on queue state. Same Level 1 calibration; only 2 doc-class header keys
    differ from Fido (`BLD_VERS, TRECROUN` extra; `LICENSE, POLICY` absent).
  - "fido": legacy `sunpy.net.Fido` + `a.jsoc.Series` + protocol="fits".
    Falls back here automatically if as-is fails.
"""
from __future__ import annotations

import glob
import logging
import os
import re
from datetime import datetime, timedelta
from functools import reduce
from pathlib import Path

from astropy.time import Time, TimeDelta
from sunpy.net import Fido, attrs as a

log = logging.getLogger(__name__)


def _parse_l1_time(f: str) -> datetime:
    """Parse the timestamp from a local L1 filename. Handles both Fido's dashed
    convention (`2022-03-01T071647Z`) and drms's dashless (`20220301T071647Z`)."""
    ts = os.path.basename(f).split(".")[2]
    for fmt in ("%Y-%m-%dT%H%M%SZ", "%Y%m%dT%H%M%SZ"):
        try:
            return datetime.strptime(ts, fmt)
        except ValueError:
            continue
    raise ValueError(f"can't parse L1 timestamp from {f}")


def _normalize_drms_fname(path: Path) -> Path:
    """Rename drms's dashless naming to the pipeline's dashed convention so
    `find_local`'s minute-glob hits both backends uniformly.

    `aia.lev1_euv_12s.20220301T071647Z.171.image_lev1.fits`
        -> `aia.lev1_euv_12s.2022-03-01T071647Z.171.image_lev1.fits`
    """
    m = re.match(r"(aia\.lev1_euv_12s\.)(\d{4})(\d{2})(\d{2})T(\d{6})Z(.*)", path.name)
    if not m:
        return path
    new_name = f"{m.group(1)}{m.group(2)}-{m.group(3)}-{m.group(4)}T{m.group(5)}Z{m.group(6)}"
    new_path = path.parent / new_name
    if new_path != path:
        if new_path.exists():
            path.unlink()  # already normalized version present (parallel race or rerun)
        else:
            path.rename(new_path)
    return new_path


def find_local(target_time: Time, wavelengths: list[str],
              data_dir: str | Path) -> dict[str, str]:
    """Match local AIA L1 by minute-prefix glob across [-1, 0, +1] minutes;
    closest second to `target_time` wins.

    Pattern: `aia.lev1_euv_12s.{YYYY-MM-DDTHHMM}*.{wl}.image_lev1.fits`. The
    ±1-minute window catches files that fall just across a minute boundary
    from `target_time` (e.g. target 07:18:00 with a 07:17:59 file on disk).
    """
    target_dt = datetime.strptime(target_time.iso[:19], "%Y-%m-%d %H:%M:%S")
    minute_strs = [(target_dt + timedelta(minutes=d)).strftime("%Y-%m-%dT%H%M")
                   for d in (-1, 0, 1)]
    out: dict[str, str] = {}
    for wl in wavelengths:
        candidates: list[str] = []
        for ms in minute_strs:
            candidates.extend(glob.glob(os.path.join(
                str(data_dir), f"aia.lev1_euv_12s.{ms}*.{wl}.image_lev1.fits")))
        if not candidates:
            continue
        out[wl] = min(candidates,
                      key=lambda f: abs((_parse_l1_time(f) - target_dt).total_seconds()))
    return out


def fetch_drms_as_is(target_time: Time, missing_wl: list[str], data_dir: str | Path,
                     jsoc_notify: str) -> dict[str, str]:
    """JSOC `protocol="as-is"` via the drms client — bypasses the FITS-export queue.

    Submits one combined export for all wavelengths, downloads via drms,
    normalizes filenames to the pipeline's dashed convention. Picks the
    record closest to `target_time` per wavelength.
    """
    if not missing_wl:
        return {}
    if not jsoc_notify:
        raise ValueError("jsoc_notify required to download from JSOC")
    import drms  # lazy import (so non-JSOC code paths don't need it)

    # JSOC indexes aia.lev1_euv_12s by T_REC in TAI; convert target_time (UTC)
    # to TAI for the query window. AIA cadence is 12 s, so a 14 s range
    # straddling target_time guarantees 1-2 records per wavelength. We then
    # filter `req.urls` to the *closest* record per wavelength and only
    # download those indices (saves ~half the bandwidth).
    target_dt = datetime.strptime(target_time.iso[:19], "%Y-%m-%d %H:%M:%S")
    target_tai = target_time.tai
    t_a = (target_tai - TimeDelta(7, format="sec")).strftime("%Y.%m.%d_%H:%M:%S")
    t_b = (target_tai + TimeDelta(7, format="sec")).strftime("%Y.%m.%d_%H:%M:%S")
    wl_clause = ",".join(missing_wl)
    query = f"aia.lev1_euv_12s[{t_a}_TAI-{t_b}_TAI][{wl_clause}]{{image}}"
    log.info("drms as-is query: %s", query)
    c = drms.Client(email=jsoc_notify)
    req = c.export(query, protocol="as-is", method="url_quick")
    log.info("drms export status=%s n_records=%d", req.status, len(req.urls))

    # Pick the closest record per wavelength from req.urls (record format:
    # `aia.lev1_euv_12s[2022-03-01T07:17:59Z][171]{image_lev1}`).
    record_re = re.compile(r"\[(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})Z\]\[(\d+)\]")
    best: dict[str, tuple[int, float]] = {}  # wl -> (row_index, abs_dt_seconds)
    for idx, rec in enumerate(req.urls["record"].tolist()):
        m = record_re.search(rec)
        if not m:
            continue
        t_rec = datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S")
        wl = m.group(2)
        if wl not in missing_wl:
            continue
        dt_abs = abs((t_rec - target_dt).total_seconds())
        if wl not in best or dt_abs < best[wl][1]:
            best[wl] = (idx, dt_abs)
    keep_indices = sorted(v[0] for v in best.values())
    log.info("drms keeping %d of %d records (closest per wavelength)",
             len(keep_indices), len(req.urls))

    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    df = req.download(str(data_dir), index=keep_indices)

    out: dict[str, str] = {}
    for raw in df["download"].tolist():
        if not raw or not Path(raw).exists():
            continue
        p = _normalize_drms_fname(Path(raw))
        for wl in missing_wl:
            if f".{wl}.image_lev1" in p.name:
                out[wl] = str(p)
                break
    return out


def fetch_jsoc_fido(target_time: Time, missing_wl: list[str], data_dir: str | Path,
                    jsoc_notify: str) -> dict[str, str]:
    """Fallback: Fido + JSOC with `protocol="fits"` (the export queue path)."""
    if not missing_wl:
        return {}
    if not jsoc_notify:
        raise ValueError("jsoc_notify required to download from JSOC")
    log.info("Fido JSOC fits fetch: %s", missing_wl)
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


_BACKENDS = {
    "drms-as-is": fetch_drms_as_is,
    "fido": fetch_jsoc_fido,
}


def find_local_or_fetch(target_time: Time, wavelengths: list[str], data_dir: str | Path,
                       jsoc_notify: str | None, backend: str = "drms-as-is") -> dict[str, str]:
    """Resolve all wavelengths to local FITS paths. Triggers JSOC fetch for any missing.

    On `backend="drms-as-is"`, falls back to `fido` automatically if the as-is
    request raises.
    """
    if backend not in _BACKENDS:
        raise ValueError(f"unknown aia.fetch_backend={backend!r}; expected one of {list(_BACKENDS)}")

    local = find_local(target_time, wavelengths, data_dir)
    missing = [wl for wl in wavelengths if wl not in local]
    if missing:
        log.info("missing locally: %s; backend=%s", missing, backend)
        try:
            local.update(_BACKENDS[backend](target_time, missing, data_dir, jsoc_notify or ""))
        except Exception as exc:
            if backend == "drms-as-is":
                log.warning("as-is fetch failed (%s); falling back to fido", exc)
                local.update(fetch_jsoc_fido(target_time, missing, data_dir, jsoc_notify or ""))
            else:
                raise

    still_missing = [wl for wl in wavelengths if wl not in local]
    if still_missing:
        raise RuntimeError(f"could not resolve AIA L1 for wavelengths: {still_missing}")
    return local
