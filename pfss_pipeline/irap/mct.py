"""IRAP MCT: URL building, Selenium-driven form trigger, ZIP fetch, parsers.

Functions ported verbatim from test_code/IRAP汇总.ipynb (cells 1, 2, 5, 8, 10),
with print() replaced by logging and one new helper `cached_zip_or_fetch`.
"""
from __future__ import annotations

import io
import logging
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import yaml
from astropy.io import fits
from astropy.table import Table

log = logging.getLogger(__name__)

BASE_URL = "https://connect-tool.irap.omp.eu"
R_SUN_M = 6.957e8
AU_M = 1.496e11
RSUN_KM = 695700.0
REQUIRED_COLS = ["curv", "r", "lon", "lat", "br", "blon", "blat"]

_MODE_RADIO_ID = {
    "SUNTIME": "reftime2-0",
    "SCTIME": "reftime2-1",
    "SUNTIMEBW": "reftime2-2",
    "SCTIMEBW": "reftime2-3",
}
_TIME_RADIO_ID = {
    "000000": "time-0",
    "060000": "time-1",
    "120000": "time-2",
    "180000": "time-3",
}
_ALL_CORONAL_IDS = ("wso", "nso", "adapt")


# ----------------------------------------------------------------------
# URL builders
# ----------------------------------------------------------------------
def build_zip_url(sc: str, coronal: str, mode: str, date_str: str, time_str: str) -> str:
    date_compact = date_str.replace("-", "")
    fname = f"{sc}_PARKER_PFSS_{mode}_{coronal}_SCIENCE_{date_compact}T{time_str}.zip"
    return f"{BASE_URL}/static/zip_files/{fname}"


def build_api_url(sc: str, coronal: str, mode: str, date_str: str, time_str: str) -> str:
    return f"{BASE_URL}/api/{sc}/{coronal}/PARKER/{mode}/{date_str}/{time_str}"


# ----------------------------------------------------------------------
# Fetch / cache
# ----------------------------------------------------------------------
def fetch_mct_zip_bytes(sc: str, coronal: str, mode: str, date_str: str, time_str: str,
                       timeout: int = 60) -> bytes:
    """Return the raw ZIP bytes by trying the static URL then the API URL."""
    for url in (build_zip_url(sc, coronal, mode, date_str, time_str),
                build_api_url(sc, coronal, mode, date_str, time_str)):
        log.info("fetching MCT ZIP: %s", url)
        resp = requests.get(url, timeout=timeout, allow_redirects=True)
        resp.raise_for_status()
        try:
            zipfile.ZipFile(io.BytesIO(resp.content)).namelist()
            return resp.content
        except zipfile.BadZipFile:
            log.warning("response from %s is not a valid ZIP, trying next URL", url)
    raise RuntimeError(f"no valid ZIP for {sc} {coronal} {mode} {date_str} {time_str}")


def unpack_zip_bytes(content: bytes) -> dict[str, bytes]:
    files: dict[str, bytes] = {}
    with zipfile.ZipFile(io.BytesIO(content)) as zf:
        for name in zf.namelist():
            files[name] = zf.read(name)
    return files


def cached_zip_or_fetch(zip_path: Path, sc: str, coronal: str, mode: str,
                       date_str: str, time_str: str,
                       chrome_binary: str = "/usr/bin/google-chrome",
                       force: bool = False) -> dict[str, bytes]:
    """ZIP cache-first lookup. Selenium triggers MCT only on cache miss + fetch failure."""
    zip_path = Path(zip_path)
    if zip_path.exists() and not force:
        log.info("using cached MCT ZIP: %s", zip_path)
        return unpack_zip_bytes(zip_path.read_bytes())

    try:
        content = fetch_mct_zip_bytes(sc, coronal, mode, date_str, time_str)
    except (requests.HTTPError, RuntimeError) as exc:
        log.info("static ZIP fetch failed (%s); triggering MCT via Selenium", exc)
        trigger_mct_and_get_urls(sc=sc, coronal=coronal, mode=mode,
                                 date=date_str, time=time_str,
                                 chrome_binary=chrome_binary)
        content = fetch_mct_zip_bytes(sc, coronal, mode, date_str, time_str)

    zip_path.parent.mkdir(parents=True, exist_ok=True)
    zip_path.write_bytes(content)
    log.info("cached MCT ZIP -> %s (%.1f MB)", zip_path, len(content) / 1e6)
    return unpack_zip_bytes(content)


def get_ascii_text(files: dict, keyword: str) -> bytes | None:
    for name, content in files.items():
        if keyword in name.lower():
            return content
    return None


# ----------------------------------------------------------------------
# Selenium trigger (only used when cache miss + static URL fails)
# ----------------------------------------------------------------------
def trigger_mct_and_get_urls(sc: str = "SOLO", coronal: str = "ADAPT",
                             mode: str = "SUNTIME", date: str = "2022-03-03",
                             time: str = "120000", timeout: int = 120,
                             chrome_binary: str = "/usr/bin/google-chrome") -> dict:
    """Headless Chrome to fill the MCT form and wait for server-side computation."""
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC

    opts = Options()
    opts.binary_location = chrome_binary
    opts.add_argument("--headless=new")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    opts.add_argument("--disable-gpu")
    opts.add_argument("--window-size=1920,1080")

    driver = webdriver.Chrome(options=opts)
    wait = WebDriverWait(driver, timeout)
    try:
        driver.get(BASE_URL)

        for cb_id in _ALL_CORONAL_IDS:
            cb = driver.find_element(By.ID, cb_id)
            if cb.get_attribute("disabled"):
                continue
            want = (cb_id == coronal.lower())
            if cb.is_selected() != want:
                cb.click()

        mode_radio = driver.find_element(By.ID, _MODE_RADIO_ID[mode])
        if not mode_radio.is_selected():
            mode_radio.click()

        date_input = driver.find_element(By.ID, "date")
        driver.execute_script(
            "arguments[0].value = arguments[1];"
            "arguments[0].dispatchEvent(new Event('change', {bubbles: true}));",
            date_input, date,
        )

        time_radio = driver.find_element(By.ID, _TIME_RADIO_ID[time])
        if not time_radio.is_selected():
            time_radio.click()

        driver.find_element(By.ID, "submit").click()

        dt = date.replace("-", "") + "T" + time
        stem_prefix = f"{sc}_PARKER_PFSS_{mode}_{coronal}"
        fits_xpath = (
            f"//a[contains(@href, '{stem_prefix}') and contains(@href, '{dt}') "
            f"and contains(@href, '_filefieldline.fits')]"
        )
        fits_elem = wait.until(EC.presence_of_element_located((By.XPATH, fits_xpath)))
        fits_url = fits_elem.get_attribute("href")
        conn_url = driver.find_element(By.XPATH,
            f"//a[contains(@href, '{stem_prefix}') and contains(@href, '{dt}') "
            f"and contains(@href, '_fileconnectivity.ascii')]").get_attribute("href")
        hcs_url = driver.find_element(By.XPATH,
            f"//a[contains(@href, '{stem_prefix}') and contains(@href, '{dt}') "
            f"and contains(@href, '_filehcs.ascii')]").get_attribute("href")
        log.info("MCT triggered: fits=%s conn=%s hcs=%s", fits_url, conn_url, hcs_url)
        return {"fieldline_fits": fits_url, "connectivity_ascii": conn_url, "hcs_ascii": hcs_url}
    finally:
        driver.quit()


# ----------------------------------------------------------------------
# Parsers (verbatim from notebook)
# ----------------------------------------------------------------------
def parse_params(content: bytes) -> dict:
    text = content.decode("utf-8", errors="ignore")
    yaml_data = yaml.safe_load(text)
    spacecraft_name = next(iter(yaml_data))
    spacecraft_data = yaml_data[spacecraft_name]
    R_SUN_KM = R_SUN_M / 1000
    AU_KM = AU_M / 1000

    def parse_position(position_list):
        if position_list is None:
            return None
        radius_km = position_list[0]
        return {
            "R_km": radius_km,
            "R_AU": radius_km / AU_KM,
            "R_Rsun": radius_km / R_SUN_KM,
            "lon_HGC_deg": position_list[1],
            "lat_HGC_deg": position_list[2],
        }

    return {
        "sc_name": spacecraft_name,
        "sc_position": parse_position(spacecraft_data.get("position_sc")),
        "source_surface_footpoints": {
            n: parse_position(p)
            for n, p in (spacecraft_data.get("position_ss") or {}).items() if p is not None
        },
        "photospheric_footpoints": {
            n: parse_position(p)
            for n, p in (spacecraft_data.get("main_connect_point") or {}).items() if p is not None
        },
        "metadata": {
            "date_insitu": str(spacecraft_data.get("date_in", "")),
            "date_surf": str(spacecraft_data.get("date_surf", "")),
            "coronal_model": spacecraft_data.get("cmodel"),
            "mag_input": spacecraft_data.get("magtype"),
            "realization_adapt": spacecraft_data.get("realization_adapt"),
            "helio_model": spacecraft_data.get("hmodel"),
            "source_surface_Rsun": spacecraft_data.get("rss"),
            "polarity": spacecraft_data.get("polarity", {}),
            "vhelio_km_s": spacecraft_data.get("vhelio", {}),
            "reliability_score": spacecraft_data.get("score"),
        },
    }


_CONN_COL_RENAME = {
    "SSW/FSW/M": "type",
    "density(%)": "prob",
    "CRLT(degrees)": "lat_CR",
    "CRLN(degrees)": "lon_CR",
    "R(m)": "R_m",
    "DIST(m)": "DIST_m",
    "HPLT(degrees)": "HPLT_deg",
    "HPLN(degrees)": "HPLN_deg",
}


def _clean_hdr_key(k: str) -> str:
    return k.replace("(m)", "_m").replace("(degrees)", "_deg")


def parse_connectivity(buf) -> tuple:
    """Parse *_fileconnectivity.ascii. Accepts a path, StringIO, or bytes."""
    if isinstance(buf, (bytes, bytearray)):
        buf = io.StringIO(buf.decode("utf-8", errors="ignore"))
    elif isinstance(buf, (str, Path)) and Path(buf).exists():
        buf = open(buf)

    comments, data_lines = [], []
    for ln in buf:
        s = ln.rstrip("\n").strip()
        if not s:
            continue
        (comments if s.startswith("#") else data_lines).append(s)

    anchor_idx = None
    for i, cl in enumerate(comments):
        if "Data is formatted as follow" in cl:
            anchor_idx = i + 1
            break
    if anchor_idx is None:
        raise ValueError("could not find '#Data is formatted as follow:' anchor")

    fmt_spec_lines = []
    col_header_line = None
    for cl in comments[anchor_idx:]:
        content = cl.lstrip("#").strip()
        if content.startswith("SSW/FSW/M") or ("SSW" in content and "FSW" in content and "density" in content):
            col_header_line = content
            break
        fmt_spec_lines.append(content)
    if col_header_line is None:
        raise ValueError("could not find '#SSW/FSW/M …' column header")

    header_dict = {}
    for fmt_line, data_line in zip(fmt_spec_lines, data_lines):
        keys = fmt_line.split()
        vals = data_line.split()
        if len(keys) == 1 and len(vals) > 1:
            header_dict[_clean_hdr_key(keys[0])] = " ".join(vals)
        else:
            for k, v in zip(keys, vals):
                header_dict[_clean_hdr_key(k)] = v
    header_df = pd.DataFrame([header_dict])

    n_meta = len(fmt_spec_lines)
    npoint_tot = int(data_lines[n_meta - 1].split()[0])
    conn_rows = data_lines[n_meta:n_meta + npoint_tot]
    if len(conn_rows) != npoint_tot:
        raise ValueError(f"expected {npoint_tot} connectivity rows, found {len(conn_rows)}")

    raw_cols = col_header_line.split()
    renamed = [_CONN_COL_RENAME.get(c, c) for c in raw_cols]
    conn_df = pd.DataFrame([row.split() for row in conn_rows], columns=renamed)
    for col in [c for c in renamed if c not in ("type", "i")]:
        conn_df[col] = conn_df[col].astype(float)
    conn_df["i"] = conn_df["i"].astype(int)
    return comments, header_df, conn_df


def parse_fieldline_fits(source) -> pd.DataFrame:
    """Parse MCT fieldline FITS. Accepts HDUList, path, or bytes."""
    if isinstance(source, (bytes, bytearray)):
        hdul = fits.open(io.BytesIO(source))
        close_after = True
    elif isinstance(source, fits.HDUList):
        hdul = source
        close_after = False
    else:
        hdul = fits.open(source)
        close_after = True

    try:
        tab = None
        for hdu in hdul:
            if not isinstance(hdu, (fits.BinTableHDU, fits.TableHDU)):
                continue
            if hdu.data is None:
                continue
            names = [c.name.lower() for c in hdu.columns]
            if all(c in names for c in REQUIRED_COLS):
                tab = Table(hdu.data).to_pandas()
                break
        if tab is None:
            raise RuntimeError(f"no HDU with columns {REQUIRED_COLS} found")
    finally:
        if close_after:
            hdul.close()

    return pd.DataFrame({
        "curv": tab["curv"].astype(float),
        "r_km": tab["r"].astype(float),
        "R_Rsun": tab["r"].astype(float) / RSUN_KM,
        "lon_HGC_deg": np.mod(tab["lon"].astype(float), 360.0),
        "lat_HGC_deg": tab["lat"].astype(float),
        "br": tab["br"].astype(float),
        "blon": tab["blon"].astype(float),
        "blat": tab["blat"].astype(float),
    })


def parse_hcs(buf) -> pd.DataFrame:
    """Parse *_filehcs.ascii. Accepts path, StringIO, or bytes."""
    if isinstance(buf, (bytes, bytearray)):
        buf = io.StringIO(buf.decode("utf-8", errors="ignore"))
    elif isinstance(buf, (str, Path)) and Path(buf).exists():
        buf = open(buf)

    rows = []
    for line in buf:
        line = line.strip()
        if not line or line.startswith("#") or line.lower().startswith("npoint"):
            continue
        p = line.split()
        if len(p) < 4:
            continue
        rows.append({
            "ipoint": int(p[0]),
            "R_m": float(p[1]),
            "lat_deg": float(p[2]),
            "lon_deg": float(p[3]),
        })
    return pd.DataFrame(rows)
