"""stage_aia_fetch: resolve AIA Level 1 paths (local-first, JSOC fallback).

Splits the network-bound part out of stage_aia_prep so a batch driver can
run all events' fetches serially (JSOC export hash collides under
concurrency) and then PSF-prep them in parallel.
"""
from __future__ import annotations

import logging
from pathlib import Path

from pfss_pipeline import manifest as mfst
from pfss_pipeline.aia import fetch as fetch_mod

log = logging.getLogger(__name__)


def run(cfg: dict, layout, force: bool = False) -> dict:
    if layout.target_time is None:
        raise RuntimeError("target_time not set; run stage_irap_fetch first")

    aia_cfg = cfg["aia"]
    wavelengths = aia_cfg["wavelengths"]
    layout.ensure_dirs()

    log.info("resolving %d AIA L1 wavelengths at %s",
             len(wavelengths), layout.target_time.iso)

    if force:
        # `--force` forces re-download via JSOC; drop the local-first short-circuit
        # by passing an empty data_dir to find_local (returns nothing), so all
        # wavelengths are treated as missing.
        local: dict[str, str] = {}
    else:
        local = fetch_mod.find_local(layout.target_time, wavelengths, cfg["aia_data_dir"])
    missing = [wl for wl in wavelengths if wl not in local]
    if missing:
        log.info("missing locally: %s; querying JSOC", missing)
        local.update(fetch_mod.fetch_jsoc(
            layout.target_time, missing, cfg["aia_data_dir"], aia_cfg.get("jsoc_notify") or "",
        ))

    still_missing = [wl for wl in wavelengths if wl not in local]
    if still_missing:
        raise RuntimeError(f"could not resolve AIA L1 for: {still_missing}")

    payload = {wl: str(local[wl]) for wl in wavelengths}
    payload["data_dir"] = str(cfg["aia_data_dir"])
    mfst.update_stage(layout.manifest_path, "aia_fetch", payload)
    log.info("aia_fetch complete; %d L1 files in %s", len(local), cfg["aia_data_dir"])
    return {wl: Path(p) for wl, p in local.items()}
