"""stage_aia_prep: AIA Level 1 → calibrated prepared FITS, one per wavelength."""
from __future__ import annotations

import glob
import logging
import os
import re
from datetime import datetime
from pathlib import Path

from astropy import units as u
from astropy.time import Time
from aiapy.calibrate.utils import get_correction_table, get_pointing_table
from sunpy.map import Map

from pfss_pipeline import io_utils, manifest as mfst
from pfss_pipeline.aia import fetch as fetch_mod
from pfss_pipeline.aia import plots as plots_mod
from pfss_pipeline.aia import prep as prep_mod

log = logging.getLogger(__name__)


def _find_existing_prep(prep_dir: Path, wl: str, target_time: Time,
                       tol_min: float) -> str | None:
    """Find a prepped FITS within `tol_min` of target_time. Returns path or None."""
    candidates = glob.glob(os.path.join(str(prep_dir), f"aia_prep_{wl}A_*.fits"))
    if not candidates:
        return None
    target_dt = datetime.strptime(target_time.iso[:19], "%Y-%m-%d %H:%M:%S")
    best, best_d = None, float("inf")
    for f in candidates:
        mt = re.search(r"(\d{8})_(\d{6})\.fits$", f)
        if not mt:
            continue
        ft = datetime.strptime(mt.group(1) + mt.group(2), "%Y%m%d%H%M%S")
        d = abs((ft - target_dt).total_seconds())
        if d < best_d:
            best_d, best = d, f
    return best if best and best_d <= tol_min * 60 else None


def run(cfg: dict, layout, force: bool = False) -> dict:
    if layout.target_time is None:
        raise RuntimeError("target_time not set; run stage_irap first")

    aia_cfg = cfg["aia"]
    wavelengths = aia_cfg["wavelengths"]
    tol_min = cfg["dem"]["match_tolerance_min"]
    layout.ensure_dirs()

    # ---- 1. Skip-if-exists: scan prep dir with tolerance ----
    existing = {wl: _find_existing_prep(layout.aia_prep_dir, wl, layout.target_time, tol_min)
                for wl in wavelengths}
    needed = [wl for wl in wavelengths if force or existing[wl] is None]

    if not needed:
        log.info("all %d wavelengths already prepped within %.1f min; skipping",
                 len(wavelengths), tol_min)
        result = {wl: existing[wl] for wl in wavelengths}
        result["images_dir"] = layout.aia_prep_images_dir
        mfst.update_stage(layout.manifest_path, "aia_prep",
                         {wl: str(existing[wl]) for wl in wavelengths} |
                         {"images_dir": str(layout.aia_prep_images_dir), "skipped": True})
        return result

    log.info("preparing %d wavelengths: %s (existing: %s)",
             len(needed), needed, [wl for wl in wavelengths if existing[wl]])

    # ---- 2. Resolve L1 files (local first; JSOC fallback) ----
    local_files = fetch_mod.find_local_or_fetch(
        layout.target_time, needed, cfg["aia_data_dir"], aia_cfg.get("jsoc_notify"),
    )

    # ---- 3. Build pointing + correction tables once ----
    sample = Map(local_files[needed[0]])
    pwh = aia_cfg["pointing_window_hours"]
    log.info("fetching pointing table (±%d h around %s)", pwh, sample.date.iso)
    pointing_table = get_pointing_table("JSOC", time_range=(sample.date - pwh * u.h,
                                                            sample.date + pwh * u.h))
    correction_table = get_correction_table()

    # ---- 4. Per-wavelength prep + comparison PNG ----
    do_psf = aia_cfg["do_psf_deconvolve"]
    new_paths: dict[str, Path] = {}
    for wl in needed:
        log.info("[%s/%s] %s Å", wavelengths.index(wl) + 1, len(wavelengths), wl)
        # Stamp the prepared file with the L1 exposure time, not target_time, so
        # different events that share an L1 file deduplicate naturally.
        ts = datetime.strptime(os.path.basename(local_files[wl]).split(".")[2],
                              "%Y-%m-%dT%H%M%SZ")
        save_path = layout.aia_prep_fits(wl, Time(ts))
        png_path = layout.aia_prep_png(wl, Time(ts))
        prep_mod.prepare_one(
            local_files[wl], pointing_table, correction_table,
            do_psf=do_psf, save_path=save_path,
        )
        new_paths[wl] = save_path
        if force or not png_path.exists():
            plots_mod.save_comparison(
                local_files[wl], save_path, png_path,
                dpi=aia_cfg["diagnostic_dpi"],
            )

    # ---- 5. Manifest ----
    final = {wl: (new_paths.get(wl) or existing[wl]) for wl in wavelengths}
    payload = {wl: str(final[wl]) for wl in wavelengths}
    payload.update({"images_dir": str(layout.aia_prep_images_dir), "do_psf": do_psf})
    mfst.update_stage(layout.manifest_path, "aia_prep", payload)
    log.info("stage_aia_prep complete; %d wavelengths in %s", len(wavelengths), layout.aia_prep_dir)

    result = dict(final)
    result["images_dir"] = layout.aia_prep_images_dir
    return result
