"""DEM-derived map helpers (verbatim from test_code/DEM_note.ipynb)."""
from __future__ import annotations

import numpy as np
from scipy.interpolate import griddata
from sunpy.map import Map


def dem_temperature_maps(dem, temperatures, clip_negative: bool = True) -> dict:
    """Convert a (ny, nx, nt) DEM cube + bin edges to T_mean / T_peak / EM maps."""
    dem = np.asarray(dem, dtype=float)
    temperatures = np.asarray(temperatures, dtype=float)
    if dem.ndim != 3:
        raise ValueError(f"dem must be (ny, nx, nt), got {dem.shape}")
    ny, nx, nt = dem.shape
    if clip_negative:
        dem = np.clip(dem, 0, None)
    if temperatures.size != nt + 1:
        raise ValueError(
            f"temperatures has incompatible size: {temperatures.size}, expected nt+1={nt + 1}"
        )

    T_edges = temperatures
    dT = np.diff(T_edges)
    logT_mid = 0.5 * (np.log10(T_edges[:-1]) + np.log10(T_edges[1:]))
    T_mid = 10 ** logT_mid

    T3 = T_mid[None, None, :]
    dT3 = dT[None, None, :]
    EM = np.sum(dem * dT3, axis=2)

    numerator = np.sum(dem * T3 * dT3, axis=2)
    T_mean = np.full((ny, nx), np.nan)
    good = EM > 0
    T_mean[good] = numerator[good] / EM[good]

    logT_mean = np.full((ny, nx), np.nan)
    goodT = T_mean > 0
    logT_mean[goodT] = np.log10(T_mean[goodT])

    peak_idx = np.argmax(dem, axis=2)
    T_peak = T_mid[peak_idx].astype(float)
    logT_peak = np.log10(T_peak)
    no_signal = np.all(dem <= 0, axis=2)
    T_peak[no_signal] = np.nan
    logT_peak[no_signal] = np.nan

    return {"T_mean": T_mean, "logT_mean": logT_mean,
            "T_peak": T_peak, "logT_peak": logT_peak, "EM": EM}


def fill_nan_2d(data, method: str = "linear") -> np.ndarray:
    """Fill NaN pixels via scipy.interpolate.griddata; falls back to nearest at the hull edge."""
    filled = np.asarray(data, dtype=float).copy()
    nan_mask = np.isnan(filled)
    if not nan_mask.any():
        return filled
    ny, nx = filled.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    valid_points = np.column_stack((yy[~nan_mask], xx[~nan_mask]))
    valid_values = filled[~nan_mask]

    filled[nan_mask] = griddata(
        valid_points, valid_values,
        np.column_stack((yy[nan_mask], xx[nan_mask])), method=method,
    )
    still_nan = np.isnan(filled)
    if still_nan.any():
        filled[still_nan] = griddata(
            valid_points, valid_values,
            np.column_stack((yy[still_nan], xx[still_nan])), method="nearest",
        )
    return filled


_WCS_KEYWORDS = {
    "naxis", "naxis1", "naxis2",
    "crpix1", "crpix2", "crval1", "crval2",
    "cdelt1", "cdelt2", "ctype1", "ctype2", "cunit1", "cunit2",
    "crota2", "pc1_1", "pc1_2", "pc2_1", "pc2_2",
    "cd1_1", "cd1_2", "cd2_1", "cd2_2",
    "lonpole", "latpole",
    "date-obs", "date_obs", "date-beg", "date-end",
    "rsun_ref", "rsun_obs", "dsun_obs",
    "hgln_obs", "hglt_obs", "crln_obs", "crlt_obs",
    "solar_b0", "solar_l0", "solar_p",
    "wcsname", "wcsaxes",
}


def make_derived_map(data_2d, reference_map, unit_str: str = "",
                     dem_settings: dict | None = None):
    """Wrap a 2D array as a sunpy Map, copying only WCS keywords from `reference_map`.

    If `dem_settings` is given (region_type / logT_range / n_T_bins / bin_width),
    the values are written into the FITS header so the saved file is self-
    describing — downstream readers can recover the inversion settings without
    a sidecar npz/manifest.
    """
    ref_meta = dict(reference_map.meta)
    new_meta = {k: v for k, v in ref_meta.items() if k.lower() in _WCS_KEYWORDS}
    new_meta["bunit"] = unit_str
    if dem_settings:
        lo, hi = dem_settings["logT_range"]
        new_meta["DEMRTYPE"] = str(dem_settings["region_type"])
        new_meta["DEMLOGT1"] = float(lo)
        new_meta["DEMLOGT2"] = float(hi)
        new_meta["DEMNBINS"] = int(dem_settings["n_T_bins"])
        new_meta["DEMBINW"] = float(dem_settings["bin_width"])
    return Map(data_2d.astype(np.float64), new_meta)
