"""Config loader: defaults, validation, stage-specific gating."""
from __future__ import annotations

import copy
from pathlib import Path

import yaml

DEFAULTS: dict = {
    "irap": {
        "spacecraft": None,
        "spacecraft_time": None,
        "coronal_model": "ADAPT",
        "mode": "SUNTIME",
        "realization_adapt": 0,
        "prob_threshold_pct": 60,
        "carrington_roi": {"lon": [10, 100], "lat": [-60, 30]},
        "adapt_cache_dir": "/disk/plasma/xw2/PFSS/data/adapt_gong",
        "selenium_chrome_binary": "/usr/bin/google-chrome",
        # Manual sw_type override; required only when 'M' is absent from the
        # IRAP solarsurf footpoints (e.g. set to 'SSW' in the event YAML).
        "sw_type": None,
    },
    "results_root": None,
    "aia_data_dir": "/disk/plasma/xw2/sunpy/data",
    "aia_prep_dir": "/disk/plasma/xw2/sunpy/data/aia_prep",
    "aia": {
        "wavelengths": ["94", "131", "171", "193", "211", "335"],
        "jsoc_notify": None,
        "do_psf_deconvolve": True,
        "pointing_window_hours": 12,
        "diagnostic_dpi": 200,
    },
    "dem": {
        "roi": {"bottom_left_arcsec": None, "top_right_arcsec": None},
        "region_type": None,
        "logT_range": None,
        "bin_width": 0.1,
        "match_tolerance_min": 5,
        "gaussian_prior": {"amp": 4.0e22, "center": 6.5, "sigma": 0.15},
        "clip_negative": True,
        "fill_nan": "linear",
    },
    "extract": {
        "hull_method": "convex",
        "overlay_aia_wavelength": "193",
        # Overlay framing in Carrington (deg). If null, derive from hull bbox + pad.
        "overlay_carrington_roi": {"lon": None, "lat": None},
        "overlay_carrington_pad_deg": 15.0,
    },
    "plots": {"dpi": 300},
    "runtime": {"log_level": "INFO"},
}


def _deep_merge(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: str | Path) -> dict:
    """Load YAML config and merge with defaults. Validates pre-IRAP requirements only."""
    with open(path) as fh:
        user = yaml.safe_load(fh) or {}
    cfg = _deep_merge(DEFAULTS, user)
    _validate_pre_irap(cfg)
    return cfg


def _validate_pre_irap(cfg: dict) -> None:
    missing = []
    if not cfg["irap"]["spacecraft"]:
        missing.append("irap.spacecraft")
    if not cfg["irap"]["spacecraft_time"]:
        missing.append("irap.spacecraft_time")
    if not cfg["results_root"]:
        missing.append("results_root")
    if missing:
        raise ValueError(f"missing required config keys: {missing}")


def assert_dem_ready(cfg: dict) -> None:
    """Stages dem/extract require ROI + region_type. Called by those stages, not at load time."""
    roi = cfg["dem"]["roi"]
    bl = roi.get("bottom_left_arcsec")
    tr = roi.get("top_right_arcsec")
    if not (bl and tr and len(bl) == 2 and len(tr) == 2):
        raise ValueError(
            "dem.roi.bottom_left_arcsec/top_right_arcsec must be set "
            "after stage_irap_fetch (inspect footpoint figure first)."
        )
    rt = cfg["dem"]["region_type"]
    if rt not in ("CH", "AR", "custom"):
        raise ValueError(f"dem.region_type must be 'CH'|'AR'|'custom', got {rt!r}")
    if rt == "custom" and not cfg["dem"].get("logT_range"):
        raise ValueError("dem.region_type='custom' requires dem.logT_range to be set")
