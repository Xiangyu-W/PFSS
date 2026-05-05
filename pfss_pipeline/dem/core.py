"""DEM core: T-response, ROI submaps + errors, gaussian prior, dn2dem inversion."""
from __future__ import annotations

import logging
import math

import numpy as np
import scipy.io as scio
from aiapy.calibrate import estimate_error
from aiapy.calibrate.utils import get_error_table
from astropy import units as u
from astropy.coordinates import SkyCoord
from demregpy import dn2dem
from demregpy.tresp import aia_tresp

log = logging.getLogger(__name__)


def load_t_response() -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (T_resp_logt, T_resp_matrix, channel_names)."""
    T_response = scio.readsav(aia_tresp)
    channels = [c.decode("utf-8") if isinstance(c, bytes) else str(c)
                for c in T_response["channels"]]
    T_resp_logt = np.array(T_response["logt"])
    n_channels = len(T_response["tr"][:])
    matrix = np.zeros((len(T_resp_logt), n_channels))
    for i in range(n_channels):
        matrix[:, i] = T_response["tr"][i]
    return T_resp_logt, matrix, channels


def build_submaps(aia_maps: dict, wavenum: list[str], roi: dict,
                  error_table_source: str = "SSW") -> tuple:
    """Crop each wavelength to ROI, estimate errors, return aligned data cubes.

    Returns (submaps_list, submaps_data, submaps_data_err, common_shape).
    submaps_data and submaps_data_err are (ny, nx, nf) with units stripped, exposure-normalised.
    """
    bl = SkyCoord(Tx=roi["bottom_left_arcsec"][0] * u.arcsec,
                  Ty=roi["bottom_left_arcsec"][1] * u.arcsec,
                  frame=aia_maps[wavenum[0]].coordinate_frame)
    tr = SkyCoord(Tx=roi["top_right_arcsec"][0] * u.arcsec,
                  Ty=roi["top_right_arcsec"][1] * u.arcsec,
                  frame=aia_maps[wavenum[0]].coordinate_frame)
    error_table = get_error_table(error_table_source)

    submaps, submaps_errors = [], []
    for wl in wavenum:
        sm = aia_maps[wl].submap(bl, top_right=tr)
        err = estimate_error(sm.quantity / u.pix, sm.wavelength, error_table=error_table)
        submaps.append(sm)
        submaps_errors.append(err)

    shapes = [sm.data.shape for sm in submaps]
    min_ny = min(s[0] for s in shapes)
    min_nx = min(s[1] for s in shapes)
    log.info("submap shapes %s; common (%d, %d)", shapes, min_ny, min_nx)

    nf = len(wavenum)
    data = np.zeros((min_ny, min_nx, nf))
    err = np.zeros((min_ny, min_nx, nf))
    for j, sm in enumerate(submaps):
        data[:, :, j] = sm.data[:min_ny, :min_nx] / sm.exposure_time.value
        err[:, :, j] = submaps_errors[j][:min_ny, :min_nx] / sm.exposure_time.value
    data[data < 0] = 0
    return submaps, data, err, (min_ny, min_nx)


def gaussian_prior_weight(T_resp_logt: np.ndarray, T_bin_edges: np.ndarray,
                         amp: float, center: float, sigma: float) -> np.ndarray:
    """Build the per-bin DEM-norm weight from a gaussian prior in logT."""
    root2pi = (2.0 * math.pi) ** 0.5
    prior = (amp / (root2pi * sigma)) * np.exp(-(T_resp_logt - center) ** 2 / (2 * sigma ** 2))
    logt_centers = 0.5 * (np.log10(T_bin_edges[:-1]) + np.log10(T_bin_edges[1:]))
    w = 10 ** np.interp(logt_centers, T_resp_logt, np.log10(prior))
    return w / w.max()


def run_inversion(submaps_data, submaps_data_err, T_resp_matrix, T_resp_logt, T_bin_edges,
                 dem_norm0=None) -> tuple:
    """Wrapper around demregpy.dn2dem. Returns (dem, dem_uncertainty, logt_uncertainty, chisq, dn_reconstructed)."""
    log.info("dn2dem on %dx%d pixels x %d bins ...",
             submaps_data.shape[0], submaps_data.shape[1], len(T_bin_edges) - 1)
    return dn2dem(submaps_data, submaps_data_err, T_resp_matrix, T_resp_logt, T_bin_edges)
