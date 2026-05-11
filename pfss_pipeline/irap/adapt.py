"""ADAPT GONG magnetogram fetch + load. Verbatim port from IRAP汇总.ipynb cell 14."""
from __future__ import annotations

import gzip
import logging
import re
from datetime import datetime
from io import BytesIO
from pathlib import Path

import numpy as np
import requests
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from bs4 import BeautifulSoup
from sunpy.coordinates import get_earth
from sunpy.map import Map

log = logging.getLogger(__name__)

GONG_BASE = "https://gong.nso.edu/adapt/maps/gong"

ADAPT_RE = re.compile(
    r"(adapt"
    r"(?P<Z>\d)(?P<X>\d)(?P<A>\d)(?P<B>\d)(?P<R>\d)"
    r"_(?P<CC>\d{2})(?P<E>\w)(?P<FFF>\d{3})"
    r"_(?P<Y>\d{4})(?P<M>\d{2})(?P<D>\d{2})(?P<HH>\d{2})(?P<NN>\d{2})"
    r"_(?P<T>[aifs])(?P<II>\d{2})(?P<JJ>\d{2})(?P<KK>\d{2})(?P<LL>\d{2})"
    r"(?P<G>[nfeb])(?P<Q>\d)\.fts(?:\.gz)?)"
)


def _patch_adapt_header(header):
    """Fix non-standard WCS in raw ADAPT FITS so SunPy maps Carrington correctly."""
    header["CTYPE1"] = "CRLN-CAR"
    header["CTYPE2"] = "CRLT-CAR"
    if "DATE-OBS" not in header and "MAPTIME" in header:
        header["DATE-OBS"] = header["MAPTIME"]
    if "DATE-OBS" in header:
        obs_time = Time(header["DATE-OBS"])
        earth = get_earth(obs_time)
        header["DSUN_OBS"] = earth.radius.to("m").value
        header["HGLT_OBS"] = earth.lat.deg
        header["HGLN_OBS"] = 0.0
    if "RSUN" not in header:
        header["RSUN"] = 1 * u.Rsun.to("m")


def list_adapt_files(base_url: str, year: int, lon_type: str = "0", timeout: int = 120) -> list[tuple]:
    """List Carrington-fixed ADAPT GONG files for a given year (returns (datetime, fname, url) triples)."""
    url = base_url + "/" + str(year) + "/"
    log.info("listing ADAPT: %s", url)
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")
    results = []
    for link in soup.find_all("a"):
        href = link.get("href", "")
        m = ADAPT_RE.match(href)
        if m is None or m.group("X") != lon_type:
            continue
        dt = datetime(int(m.group("Y")), int(m.group("M")), int(m.group("D")),
                      int(m.group("HH")), int(m.group("NN")))
        results.append((dt, href, url + href))
    results.sort(key=lambda x: x[0])
    log.info("found %d ADAPT maps for year %d", len(results), year)
    return results


def download_adapt(url: str, cache_dir: str | Path) -> str:
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    fname = url.split("/")[-1]
    local = cache_dir / fname
    if local.exists():
        log.info("ADAPT cached: %s", fname)
        return str(local)
    log.info("downloading ADAPT: %s", fname)
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    local.write_bytes(r.content)
    log.info("ADAPT downloaded (%.1f MB) -> %s", len(r.content) / 1e6, local)
    return str(local)


def load_adapt(filepath: str, realization: int = 0):
    """Load ADAPT FITS (or .fts.gz). Returns (data_2d, header, sunpy_map)."""
    if filepath.endswith(".gz"):
        with gzip.open(filepath, "rb") as gz:
            with fits.open(BytesIO(gz.read())) as hdul:
                data = hdul[0].data.copy().astype(np.float64)
                header = hdul[0].header.copy()
    else:
        with fits.open(filepath) as hdul:
            data = hdul[0].data.copy().astype(np.float64)
            header = hdul[0].header.copy()

    if data.ndim == 3:
        n_real = data.shape[0]
        if not (0 <= realization < n_real):
            raise IndexError(f"realization={realization} out of range; available 0..{n_real - 1}")
        log.info("%d realizations found, using #%d", n_real, realization)
        data = data[realization]
    elif data.ndim != 2:
        raise ValueError(f"unexpected ADAPT data shape: {data.shape}")

    _patch_adapt_header(header)
    adapt_map = Map((data, header))
    return data, header, adapt_map


def make_adapt_axes(data, hdr):
    """Build 1-D longitude / latitude arrays from FITS WCS keywords."""
    nlat, nlon = data.shape
    crval1 = hdr.get("CRVAL1", hdr.get("CRVAL1A"))
    cdelt1 = hdr.get("CDELT1", hdr.get("CDELT1A"))
    crpix1 = hdr.get("CRPIX1", hdr.get("CRPIX1A"))
    crval2 = hdr.get("CRVAL2", hdr.get("CRVAL2A"))
    cdelt2 = hdr.get("CDELT2", hdr.get("CDELT2A"))
    crpix2 = hdr.get("CRPIX2", hdr.get("CRPIX2A"))
    if all(v is not None for v in [crval1, cdelt1, crpix1]):
        lon = crval1 + (np.arange(nlon) + 1 - crpix1) * cdelt1
    else:
        lon = np.linspace(0, 360, nlon, endpoint=False)
    if all(v is not None for v in [crval2, cdelt2, crpix2]):
        lat = crval2 + (np.arange(nlat) + 1 - crpix2) * cdelt2
    else:
        sinlat = np.linspace(-1, 1, nlat)
        lat = np.degrees(np.arcsin(sinlat))
    return lon, lat


def fname_info(fname: str) -> dict:
    m = ADAPT_RE.match(fname)
    if not m:
        return {"raw": fname}
    src = {"0": "All", "1": "KPVT", "2": "VSM", "3": "GONG", "4": "HMI",
           "5": "FDT", "7": "SVSM+FDT", "8": "GONG+FDT", "9": "HMI+FDT"}
    evol = {"a": "assimilation", "i": "intermediate", "f": "forecast", "s": "seedmap"}
    return dict(
        source=src.get(m.group("A"), "?"),
        version=f"v{m.group('CC')}.{m.group('E')}",
        n_real=int(m.group("FFF")),
        map_time=f"{m.group('Y')}-{m.group('M')}-{m.group('D')} {m.group('HH')}:{m.group('NN')} UT",
        evol=evol.get(m.group("T"), m.group("T")),
        lag=f"{m.group('II')}d {m.group('JJ')}h {m.group('KK')}m",
    )


def _scan_cache(cache_dir: str | Path) -> list[tuple]:
    """Local scan of an ADAPT cache dir; returns (datetime, fname, local_path) triples."""
    cache_dir = Path(cache_dir)
    if not cache_dir.exists():
        return []
    out = []
    for p in cache_dir.iterdir():
        m = ADAPT_RE.match(p.name)
        if m is None:
            continue
        dt = datetime(int(m.group("Y")), int(m.group("M")), int(m.group("D")),
                      int(m.group("HH")), int(m.group("NN")))
        out.append((dt, p.name, str(p)))
    out.sort(key=lambda x: x[0])
    return out


def find_closest_adapt(base_url: str, target_dt: datetime, cache_dir: str | Path) -> tuple:
    """Return (matched_dt, local_path, fname, source).

    Always queries GONG to identify the actual closest ADAPT map to
    `target_dt` (matches the original notebook IRAP汇总.ipynb logic).
    `download_adapt` reuses the file from `cache_dir` if already present, so
    repeated calls don't re-download. Only falls back to closest cached file
    when GONG is unreachable.
    """
    try:
        files = list_adapt_files(base_url, target_dt.year)
        if target_dt.month == 1:
            files = list_adapt_files(base_url, target_dt.year - 1) + files
    except requests.RequestException as exc:
        log.warning("GONG listing failed (%s); falling back to closest cached file", exc)
        cache_hits = _scan_cache(cache_dir)
        if not cache_hits:
            raise RuntimeError("GONG unreachable and ADAPT cache is empty") from exc
        nearest = min(cache_hits, key=lambda x: abs((x[0] - target_dt).total_seconds()))
        return nearest[0], nearest[2], nearest[1], "cache"
    if not files:
        raise RuntimeError("no ADAPT GONG files found; check network access to gong.nso.edu")
    matched = min(files, key=lambda x: abs((x[0] - target_dt).total_seconds()))
    local = download_adapt(matched[2], cache_dir)
    return matched[0], local, matched[1], matched[2]
