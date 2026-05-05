"""Footpoint utilities: dominant filter, ConvexHull boundary, differential rotation, ECSV save/load."""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from scipy.spatial import ConvexHull
from sunpy.coordinates import HeliographicCarrington, HeliographicStonyhurst, RotatedSunFrame

log = logging.getLogger(__name__)


def select_sw_type(foot_df: pd.DataFrame) -> str:
    """Prefer 'M' if any M points exist; otherwise 'SSW'."""
    return "M" if "M" in foot_df["type"].values else "SSW"


def filter_dominant(foot_df: pd.DataFrame, sw_type: str, prob_threshold_pct: float = 60.0) -> pd.DataFrame:
    """Sort by prob desc; return rows whose CUMULATIVE probability is below threshold."""
    df = (foot_df[foot_df["type"] == sw_type]
          .sort_values("prob", ascending=False).reset_index(drop=True))
    df["cumulative_prob"] = df["prob"].cumsum()
    return df[df["cumulative_prob"].shift(1, fill_value=0) < prob_threshold_pct]


def make_skycoord_carr(df: pd.DataFrame, obstime, observer_coord) -> SkyCoord:
    """Build HeliographicCarrington SkyCoord from a DataFrame with lon_CR/lat_CR/R_m."""
    return SkyCoord(
        df["lon_CR"].values * u.deg,
        df["lat_CR"].values * u.deg,
        df["R_m"].values * u.m,
        obstime=obstime, frame=HeliographicCarrington, observer=observer_coord,
    )


def diff_rotate(coord: SkyCoord, rotated_time) -> SkyCoord:
    """Apply differential rotation from coord.obstime to rotated_time."""
    return SkyCoord(RotatedSunFrame(base=coord, rotated_time=rotated_time))


def convex_hull_polygon(coord: SkyCoord) -> tuple[np.ndarray, np.ndarray]:
    """Return (hull_lons_closed, hull_lats_closed) where the polygon is closed."""
    lons = coord.lon.deg
    lats = coord.lat.deg
    hull = ConvexHull(np.column_stack([lons, lats]))
    hull_lons = np.append(lons[hull.vertices], lons[hull.vertices[0]])
    hull_lats = np.append(lats[hull.vertices], lats[hull.vertices[0]])
    return hull_lons, hull_lats


def save_dominant_ecsv(out_path: Path, dominant_ft: SkyCoord, dominant_df: pd.DataFrame,
                      params: dict, sw_type: str, prob_threshold_pct: float) -> None:
    """Persist dominant footpoints as ECSV with full SkyCoord-rebuild metadata."""
    observer = dominant_ft.observer
    tbl = QTable({
        "lon": dominant_ft.lon.to(u.deg),
        "lat": dominant_ft.lat.to(u.deg),
        "radius": dominant_ft.radius.to(u.m),
        "prob": np.asarray(dominant_df["prob"].values, dtype=float),
        "cumulative_prob": np.asarray(dominant_df["cumulative_prob"].values, dtype=float),
        "type": np.asarray(dominant_df["type"].values, dtype=str),
    })
    tbl.meta.update({
        "sc_name": str(params.get("sc_name", "")),
        "sw_type": str(sw_type),
        "date_surf": str(params["metadata"]["date_surf"]),
        "date_insitu": str(params["metadata"].get("date_insitu", "")),
        "prob_threshold_pct": prob_threshold_pct,
        "coronal_model": params["metadata"].get("coronal_model"),
        "mag_input": params["metadata"].get("mag_input"),
        "realization_adapt": params["metadata"].get("realization_adapt"),
        "helio_model": params["metadata"].get("helio_model"),
        "source_surface_Rsun": params["metadata"].get("source_surface_Rsun"),
        "frame": "HeliographicCarrington",
        "obstime": str(dominant_ft.obstime),
        "rsun_m": float(dominant_ft.frame.rsun.to_value(u.m)),
        "observer_frame": "HeliographicStonyhurst",
        "observer_obstime": str(observer.obstime),
        "observer_lon_deg": float(observer.lon.to_value(u.deg)),
        "observer_lat_deg": float(observer.lat.to_value(u.deg)),
        "observer_radius_m": float(observer.radius.to_value(u.m)),
        "observer_rsun_m": float(observer.rsun.to_value(u.m)),
    })
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tbl.write(out_path, format="ascii.ecsv", overwrite=True)
    log.info("saved %d dominant footpoints -> %s", len(tbl), out_path)


def load_dominant_ecsv(path: Path) -> tuple[QTable, SkyCoord]:
    """Reverse of save_dominant_ecsv. Returns (table, SkyCoord rebuilt with original observer)."""
    tbl = QTable.read(path, format="ascii.ecsv")
    meta = tbl.meta
    observer = SkyCoord(
        lon=meta["observer_lon_deg"] * u.deg,
        lat=meta["observer_lat_deg"] * u.deg,
        radius=meta["observer_radius_m"] * u.m,
        frame=HeliographicStonyhurst,
        obstime=meta["observer_obstime"],
    )
    coord = SkyCoord(
        tbl["lon"], tbl["lat"], tbl["radius"],
        frame=HeliographicCarrington,
        obstime=meta["obstime"],
        observer=observer,
    )
    return tbl, coord
