"""stage_dem: DEM inversion over the user-specified ROI."""
from __future__ import annotations

import logging

import numpy as np

from pfss_pipeline import io_utils, manifest as mfst
from pfss_pipeline.dem import core, derived, plots, tbins

log = logging.getLogger(__name__)


_CACHE_KEYS = (
    "region_type", "logT_range", "n_T_bins", "bin_width",
    "roi_bottom_left_arcsec", "roi_top_right_arcsec",
    "clip_negative", "fill_nan", "wavelengths",
)


def _current_signature(dem_cfg: dict, wavenum: list, n_T_bins: int, logT_range: tuple) -> dict:
    """Settings that affect T_mean values; if any change, the cached FITS is stale."""
    return {
        "region_type": dem_cfg["region_type"],
        "logT_range": [float(logT_range[0]), float(logT_range[1])],
        "n_T_bins": int(n_T_bins),
        "bin_width": float(dem_cfg["bin_width"]),
        "roi_bottom_left_arcsec": list(dem_cfg["roi"]["bottom_left_arcsec"]),
        "roi_top_right_arcsec": list(dem_cfg["roi"]["top_right_arcsec"]),
        "clip_negative": bool(dem_cfg["clip_negative"]),
        "fill_nan": dem_cfg["fill_nan"],
        "wavelengths": list(wavenum),
    }


def _cache_mismatch(cached: dict, current: dict) -> list:
    """Return list of keys whose stored value disagrees with the current config."""
    diffs = []
    for k in _CACHE_KEYS:
        if cached.get(k) != current.get(k):
            diffs.append((k, cached.get(k), current.get(k)))
    return diffs


def run(cfg: dict, layout, force: bool = False) -> dict:
    if layout.target_time is None:
        raise RuntimeError("target_time not set; run stage_irap first")

    dem_cfg = cfg["dem"]
    wavenum = cfg["aia"]["wavelengths"]
    layout.ensure_dirs()

    # ---- 1. Resolve T-response + bin grid ----
    T_resp_logt, T_resp_matrix, channels = core.load_t_response()
    T_bin_edges, n_T_bins, logT_range = tbins.resolve_grid(
        dem_cfg["region_type"], dem_cfg.get("logT_range"), dem_cfg["bin_width"],
    )
    log.info("region_type=%s, logT range=%s, n_T_bins=%d",
             dem_cfg["region_type"], logT_range, n_T_bins)
    current_sig = _current_signature(dem_cfg, wavenum, n_T_bins, logT_range)

    t_mean_path = layout.dem_t_mean()
    if not force and t_mean_path.exists():
        cached = mfst.read(layout.manifest_path).get("stages", {}).get("dem", {})
        diffs = _cache_mismatch(cached, current_sig)
        if not diffs:
            log.info("DEM already present at %s with matching config; skipping (use --force to redo)",
                     t_mean_path)
            payload = {
                "T_mean": str(t_mean_path),
                "T_mean_raw": str(layout.dem_t_mean_raw()),
                "T_peak": str(layout.dem_t_peak()),
                "EM": str(layout.dem_em()),
                "cube": str(layout.dem_cube()),
                "skipped": True,
            }
            mfst.update_stage(layout.manifest_path, "dem", payload)
            return {k: v for k, v in payload.items()}
        log.warning("DEM cache stale; settings changed -> recomputing. Diffs (cached -> current):")
        for k, old, new in diffs:
            log.warning("  %s: %r -> %r", k, old, new)

    # ---- 2. Find prepared AIA maps ----
    aia_maps = {
        wl: io_utils.find_closest_map(
            cfg["aia_prep_dir"], wl, layout.target_time,
            tolerance_minutes=dem_cfg["match_tolerance_min"],
        ) for wl in wavenum
    }

    # ---- 3. Submaps + errors ----
    submaps, data, err, shape = core.build_submaps(aia_maps, wavenum, dem_cfg["roi"])

    # ---- 4. dn2dem ----
    dem, dem_unc, logt_unc, chisq, dn_rec = core.run_inversion(
        data, err, T_resp_matrix, T_resp_logt, T_bin_edges,
    )

    # ---- 5. Derive maps ----
    res = derived.dem_temperature_maps(dem, T_bin_edges, clip_negative=dem_cfg["clip_negative"])
    T_mean_raw = res["T_mean"]
    T_mean = (derived.fill_nan_2d(T_mean_raw, method=dem_cfg["fill_nan"])
              if dem_cfg["fill_nan"] else T_mean_raw)
    log.info("filled %d NaN pixels in T_mean (remaining %d)",
             int(np.isnan(T_mean_raw).sum()), int(np.isnan(T_mean).sum()))

    # ---- 6. Save FITS via make_derived_map (use channel index 3 = 193 Å as WCS reference) ----
    ref_idx = wavenum.index("193") if "193" in wavenum else 3
    ref_map = submaps[ref_idx]
    fits_settings = {
        "region_type": current_sig["region_type"],
        "logT_range": current_sig["logT_range"],
        "n_T_bins": current_sig["n_T_bins"],
        "bin_width": current_sig["bin_width"],
    }
    derived.make_derived_map(T_mean, ref_map, unit_str="K", dem_settings=fits_settings
                            ).save(str(layout.dem_t_mean()), overwrite=True)
    derived.make_derived_map(T_mean_raw, ref_map, unit_str="K", dem_settings=fits_settings
                            ).save(str(layout.dem_t_mean_raw()), overwrite=True)
    derived.make_derived_map(res["T_peak"], ref_map, unit_str="K", dem_settings=fits_settings
                            ).save(str(layout.dem_t_peak()), overwrite=True)
    derived.make_derived_map(res["EM"], ref_map, unit_str="cm-5", dem_settings=fits_settings
                            ).save(str(layout.dem_em()), overwrite=True)

    # ---- 7. Save DEM cube + chisq ----
    np.savez(layout.dem_cube(),
             dem=dem, dem_uncertainty=dem_unc, logt_uncertainty=logt_unc,
             chisq=chisq, T_bin_edges=T_bin_edges,
             roi_bottom_left_arcsec=np.array(dem_cfg["roi"]["bottom_left_arcsec"], dtype=float),
             roi_top_right_arcsec=np.array(dem_cfg["roi"]["top_right_arcsec"], dtype=float),
             region_type=dem_cfg["region_type"], logT_range=np.array(logT_range))

    # ---- 8. Plots ----
    fig_dir = layout.dem_figures_dir
    plots.response_curves(T_resp_logt, T_resp_matrix, channels,
                         fig_dir / "response_curves.png", dpi=cfg["plots"]["dpi"])
    plots.reconstruction_comparison(data, dn_rec, err, wavenum,
                                   fig_dir / "reconstruction_comparison.png", dpi=cfg["plots"]["dpi"])
    plots.temperature_map_pixel(T_mean, fig_dir / "temperature_map.png", dpi=cfg["plots"]["dpi"])
    logt_centers = 0.5 * (np.log10(T_bin_edges[:-1]) + np.log10(T_bin_edges[1:]))
    plots.dem_bins_grid(dem, logt_centers, fig_dir / "dem_bins.png", dpi=cfg["plots"]["dpi"])

    # ---- 9. Manifest ----
    # Persist current_sig keys so the next run can detect a config change and
    # invalidate the cached T_mean.fits (instead of silently reusing it).
    payload = {
        "T_mean": str(layout.dem_t_mean()),
        "T_mean_raw": str(layout.dem_t_mean_raw()),
        "T_peak": str(layout.dem_t_peak()),
        "EM": str(layout.dem_em()),
        "cube": str(layout.dem_cube()),
        "figures_dir": str(fig_dir),
        "shape": list(shape),
        "chisq_median": float(np.nanmedian(chisq)),
        **current_sig,
    }
    mfst.update_stage(layout.manifest_path, "dem", payload)
    log.info("stage_dem complete; T_mean -> %s", layout.dem_t_mean())
    return {k: v for k, v in payload.items()}
