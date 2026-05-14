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
    backend = aia_cfg.get("fetch_backend", "drms-as-is")
    layout.ensure_dirs()

    log.info("resolving %d AIA L1 wavelengths at %s (backend=%s)",
             len(wavelengths), layout.target_time.iso, backend)

    if force:
        # `--force` forces re-download via JSOC; drop both short-circuits.
        local: dict[str, str] = {}
    else:
        # Manifest short-circuit: trust paths recorded by a previous aia_fetch
        # if every file still exists on disk. Avoids re-globbing / network when
        # the run state is already known good.
        recorded = mfst.read(layout.manifest_path).get("stages", {}).get("aia_fetch", {})
        cached = {wl: recorded[wl] for wl in wavelengths
                  if isinstance(recorded.get(wl), str) and Path(recorded[wl]).exists()}
        if len(cached) == len(wavelengths):
            log.info("all %d L1 files resolved from manifest; skipping glob/network",
                     len(wavelengths))
            return {wl: Path(p) for wl, p in cached.items()}
        local = fetch_mod.find_local(layout.target_time, wavelengths, cfg["aia_data_dir"])
    missing = [wl for wl in wavelengths if wl not in local]
    if missing:
        notify = aia_cfg.get("jsoc_notify") or ""
        log.info("missing locally: %s; backend=%s", missing, backend)
        try:
            local.update(fetch_mod._BACKENDS[backend](
                layout.target_time, missing, cfg["aia_data_dir"], notify))
        except Exception as exc:
            if backend == "drms-as-is":
                log.warning("as-is fetch failed (%s); falling back to fido", exc)
                local.update(fetch_mod.fetch_jsoc_fido(
                    layout.target_time, missing, cfg["aia_data_dir"], notify))
            else:
                raise

    still_missing = [wl for wl in wavelengths if wl not in local]
    if still_missing:
        raise RuntimeError(f"could not resolve AIA L1 for: {still_missing}")

    payload = {wl: str(local[wl]) for wl in wavelengths}
    payload["data_dir"] = str(cfg["aia_data_dir"])
    payload["backend"] = backend
    mfst.update_stage(layout.manifest_path, "aia_fetch", payload)
    log.info("aia_fetch complete; %d L1 files in %s", len(local), cfg["aia_data_dir"])
    return {wl: Path(p) for wl, p in local.items()}
