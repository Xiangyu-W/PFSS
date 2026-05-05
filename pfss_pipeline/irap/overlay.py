"""Plotting helpers: contour levels, AIA→Carrington reprojection, footpoint overlays."""
from __future__ import annotations

import logging

import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from sunpy.coordinates import HeliographicCarrington
from sunpy.map.header_helper import make_fitswcs_header, make_heliographic_header

log = logging.getLogger(__name__)


def make_mag_contour_levels(mag_map, n_levels: int = 6,
                            vmin_percentile: float = 80.0,
                            vmax_percentile: float = 99.5) -> u.Quantity:
    """Symmetric, log-spaced contour levels adapted to the actual |B| range."""
    abs_data = np.abs(mag_map.data[np.isfinite(mag_map.data)])
    nonzero = abs_data[abs_data > 0]
    v_min = np.percentile(nonzero, vmin_percentile)
    v_max = np.percentile(nonzero, vmax_percentile)
    pos_levels = np.logspace(np.log10(v_min), np.log10(v_max), n_levels)
    return np.concatenate((-pos_levels[::-1], pos_levels)) * mag_map.unit


def adapt_vmin_vmax(adapt_data: np.ndarray) -> tuple[float, float]:
    valid = adapt_data[adapt_data != -9999.0]
    valid = valid[np.isfinite(valid)]
    vmin = float(np.percentile(valid, 0.2))
    vmax = float(np.percentile(valid, 99.8))
    return max(vmin, -100.0), min(vmax, 100.0)


def reproject_aia_to_carrington_full(aia_map, shape: tuple[int, int] = (720, 1440),
                                     crval1: float = 180.0):
    """Reproject AIA HPC map to a full-disk Carrington raster.

    Reprojecting AIA directly onto a small ROI-scoped CAR WCS leaves limb
    artefacts (curved no-data regions). Going through the full-disk Carrington
    grid first and submapping afterwards is the clean route.
    """
    carr_header = make_heliographic_header(aia_map.date, aia_map.observer_coordinate,
                                           shape, frame="carrington")
    carr_header["CRVAL1"] = crval1
    return aia_map.reproject_to(carr_header)


def make_carrington_roi_wcs(observer, obstime,
                            lon_lo: float, lon_hi: float,
                            lat_lo: float, lat_hi: float,
                            shape: tuple[int, int]) -> WCS:
    """Build a local Carrington-CAR WCS over (lon, lat) bounds at the given output shape.

    The ROI center becomes CRVAL, and the pixel scale is derived as ROI_extent /
    shape. Output shape is set explicitly (not derived from a target pixel scale)
    to keep the hull-pixel statistics scale-invariant: we want the Carrington
    raster to preserve the source DEM's native pixel count so reprojection
    neither up- nor down-samples. Mirrors the prototype notebook's approach
    (`make_fitswcs_header((ny_in, nx_in), ..., scale=30.0/nx_in)`).
    """
    ny, nx = int(shape[0]), int(shape[1])
    scale_lon = (lon_hi - lon_lo) / nx
    scale_lat = (lat_hi - lat_lo) / ny
    center = SkyCoord(
        0.5 * (lon_lo + lon_hi) * u.deg,
        0.5 * (lat_lo + lat_hi) * u.deg,
        frame=HeliographicCarrington,
        obstime=obstime,
        observer=observer,
    )
    header = make_fitswcs_header(
        (ny, nx),
        center,
        scale=u.Quantity([scale_lon, scale_lat], u.deg / u.pix),
        projection_code="CAR",
    )
    wcs = WCS(header)
    wcs.array_shape = (ny, nx)
    return wcs


def plot_footpoints_on_adapt(adapt_map, foot_diffrot: SkyCoord, prob_values: np.ndarray,
                             hcs_diffrot: SkyCoord | None = None,
                             carrington_roi: dict | None = None,
                             title: str | None = None,
                             figsize: tuple = (16, 8)):
    """Full-disk footpoints-on-ADAPT figure (or submap if carrington_roi given)."""
    vmin, vmax = adapt_vmin_vmax(adapt_map.data)
    mp = adapt_map
    if carrington_roi:
        bl = SkyCoord(carrington_roi["lon"][0] * u.deg, carrington_roi["lat"][0] * u.deg,
                      frame=adapt_map.coordinate_frame)
        tr = SkyCoord(carrington_roi["lon"][1] * u.deg, carrington_roi["lat"][1] * u.deg,
                      frame=adapt_map.coordinate_frame)
        mp = adapt_map.submap(bottom_left=bl, top_right=tr)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(projection=mp)
    mp.plot(axes=ax, vmin=vmin, vmax=vmax, cmap="gray")
    bound = ax.axis()  # capture submap data extent
    sorted_idx = np.argsort(prob_values)
    ft = ax.scatter_coord(
        foot_diffrot[sorted_idx],
        alpha=0.9, vmin=0, vmax=float(prob_values.max()),
        c=prob_values[sorted_idx], cmap="YlOrRd",
        s=30, edgecolors="k", linewidths=0.5, zorder=6,
        label="Footpoints on solar surface",
    )
    if hcs_diffrot is not None:
        ax.plot_coord(hcs_diffrot, color="white", zorder=4, label="HCS", linestyle="--")
    ax.axis(bound)  # restore — keeps overlays from expanding the axes
    ax.legend(loc="upper right", fontsize=8, framealpha=0.7)
    fig.colorbar(ft, ax=ax, label="Probability (%)", shrink=0.8, pad=0.01)
    if title:
        ax.set_title(title)
    return fig
