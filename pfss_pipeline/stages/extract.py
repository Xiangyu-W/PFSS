"""stage_extract: clip DEM to ROI, mask by footpoint hull, save npz + summary CSV + overlay PNG.

Final deliverable. Reads:
  - stage_irap_fetch outputs (dominant_ft ECSV, ADAPT path, footpoint hull)
  - stage_dem outputs (T_mean / T_peak / EM FITS)
  - stage_aia_prep outputs (AIA 193 Å prepped for the Carrington overlay)

Writes:
  - T_mean / T_peak / EM clipped to HPC ROI as FITS
  - {SC}_T_inside_hull_*.npz with arrays values/lon/lat for pixels inside the
    diff-rotated dominant-footpoint hull
  - DEM_on_AIA_overlay PNG: AIA Carrington submap + ADAPT |B| contours +
    T_mean half-transparent + dominant footpoints + hull boundary + hull pixels
  - footpoint_temperature_summary.csv (one row per run; appended on rerun)
"""
from __future__ import annotations

import csv
import logging
from datetime import datetime, timezone

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.path import Path as MplPath
from sunpy.coordinates import RotatedSunFrame
from sunpy.map import Map

from pfss_pipeline import io_utils, manifest as mfst
from pfss_pipeline.irap import adapt as adapt_mod
from pfss_pipeline.irap import footpoints as fp_mod
from pfss_pipeline.irap import overlay

log = logging.getLogger(__name__)


SUMMARY_COLUMNS = [
    "event_id", "sc_time_utc", "target_time_utc", "spacecraft", "coronal_model", "mode",
    "region_type", "roi_blx", "roi_bly", "roi_trx", "roi_try",
    "prob_threshold_pct", "n_dominant_footpoints",
    "hull_npix_total", "hull_npix_valid", "hull_npix_nan",
    "T_mean_hull_avg",
    "completed_at",
]


def _resolve_overlay_roi(cfg: dict, hull_lons: np.ndarray, hull_lats: np.ndarray) -> tuple:
    """Return (lon_lo, lon_hi, lat_lo, lat_hi) in deg for the Carrington overlay window."""
    ov = cfg["extract"].get("overlay_carrington_roi") or {}
    lon = ov.get("lon")
    lat = ov.get("lat")
    if lon and lat:
        return float(lon[0]), float(lon[1]), float(lat[0]), float(lat[1])
    pad = float(cfg["extract"].get("overlay_carrington_pad_deg", 15.0))
    return (float(hull_lons.min() - pad), float(hull_lons.max() + pad),
            float(hull_lats.min() - pad), float(hull_lats.max() + pad))


def _append_summary_row(csv_path, row: dict) -> None:
    """Append `row` to summary CSV. If the existing header doesn't match the current
    schema, rotate it aside so the new file starts clean."""
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    if csv_path.exists():
        with open(csv_path) as fh:
            first_line = fh.readline().rstrip("\n")
        existing_cols = first_line.split(",")
        if existing_cols != SUMMARY_COLUMNS:
            archive = csv_path.with_suffix(
                f".old-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')}.csv"
            )
            csv_path.rename(archive)
            log.info("schema mismatch in summary CSV; archived old file -> %s", archive)
    write_header = not csv_path.exists()
    with open(csv_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_COLUMNS)
        if write_header:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in SUMMARY_COLUMNS})


def run(cfg: dict, layout, force: bool = False) -> dict:
    if layout.target_time is None:
        raise RuntimeError("target_time not set; run stage_irap first")

    layout.ensure_dirs()
    out = {
        "T_mean_in_roi": layout.extract_t_mean_in_roi(),
        "T_peak_in_roi": layout.extract_t_peak_in_roi(),
        "EM_in_roi": layout.extract_em_in_roi(),
        "T_in_hull_npz": layout.extract_t_in_hull_npz(),
        "overlay_png": layout.extract_overlay_png(),
        "summary_csv": layout.extract_summary_csv,
    }

    # ---- 1. Load DEM products + dominant footpoints + AIA prep + ADAPT --------
    T_mean_map = Map(str(layout.dem_t_mean()))
    T_peak_map = Map(str(layout.dem_t_peak()))
    EM_map = Map(str(layout.dem_em()))

    tbl, dominant_ft = fp_mod.load_dominant_ecsv(layout.irap_dominant_ft_path)
    log.info("loaded %d dominant footpoints from %s", len(tbl), layout.irap_dominant_ft_path)

    manifest = mfst.read(layout.manifest_path)
    adapt_path = manifest["stages"]["irap_fetch"]["adapt"]
    realization = cfg["irap"]["realization_adapt"]
    _, _, adapt_map = adapt_mod.load_adapt(adapt_path, realization=realization)

    aia_wl = cfg["extract"]["overlay_aia_wavelength"]
    aia_map = io_utils.find_closest_map(
        cfg["aia_prep_dir"], aia_wl, layout.target_time,
        tolerance_minutes=cfg["dem"]["match_tolerance_min"],
    )

    # ---- 2. ROI clip in HPC --------------------------------------------------
    roi = cfg["dem"]["roi"]
    bl = SkyCoord(Tx=roi["bottom_left_arcsec"][0] * u.arcsec,
                  Ty=roi["bottom_left_arcsec"][1] * u.arcsec,
                  frame=T_mean_map.coordinate_frame)
    tr = SkyCoord(Tx=roi["top_right_arcsec"][0] * u.arcsec,
                  Ty=roi["top_right_arcsec"][1] * u.arcsec,
                  frame=T_mean_map.coordinate_frame)
    T_mean_map.submap(bl, top_right=tr).save(str(out["T_mean_in_roi"]), overwrite=True)
    T_peak_map.submap(bl, top_right=tr).save(str(out["T_peak_in_roi"]), overwrite=True)
    EM_map.submap(bl, top_right=tr).save(str(out["EM_in_roi"]), overwrite=True)

    # ---- 3. Diff-rotate dominant footpoints to AIA epoch, build hull --------
    dominant_ft_diffrot = SkyCoord(
        RotatedSunFrame(base=dominant_ft, rotated_time=aia_map.date)
    )
    hull_lons, hull_lats = fp_mod.convex_hull_polygon(dominant_ft_diffrot)
    hull_path = MplPath(np.column_stack([hull_lons[:-1], hull_lats[:-1]]))

    # ---- 4. AIA -> full-disk Carrington (clean reproject) + submap; DEM -> local fine ROI WCS
    lon_lo, lon_hi, lat_lo, lat_hi = _resolve_overlay_roi(cfg, hull_lons, hull_lats)
    aia_carr_full = overlay.reproject_aia_to_carrington_full(aia_map)
    bl_carr = SkyCoord(lon_lo * u.deg, lat_lo * u.deg, frame=aia_carr_full.coordinate_frame)
    tr_carr = SkyCoord(lon_hi * u.deg, lat_hi * u.deg, frame=aia_carr_full.coordinate_frame)
    aia_carr_sub = aia_carr_full.submap(bottom_left=bl_carr, top_right=tr_carr)
    adapt_sub = adapt_map.submap(bottom_left=bl_carr, top_right=tr_carr)

    # Output Carrington shape = T_mean's native HPC shape: prevents up/down-
    # sampling, keeps the per-pixel area constant across runs, and makes the
    # hull statistics scale-invariant (see notebook IRAP汇总.ipynb cell that
    # builds make_fitswcs_header from `T_mean_map.data.shape`).
    roi_wcs = overlay.make_carrington_roi_wcs(
        observer=aia_map.observer_coordinate, obstime=aia_map.date,
        lon_lo=lon_lo, lon_hi=lon_hi, lat_lo=lat_lo, lat_hi=lat_hi,
        shape=T_mean_map.data.shape,
    )
    T_mean_carr = T_mean_map.reproject_to(roi_wcs)

    # ---- 5. Hull pixel mask on Carrington-reprojected T_mean ----------------
    ny_t, nx_t = T_mean_carr.data.shape
    y_idx, x_idx = np.mgrid[0:ny_t, 0:nx_t]
    world = T_mean_carr.wcs.pixel_to_world(x_idx.ravel(), y_idx.ravel())
    T_pix_lonlat = np.column_stack([world.lon.deg, world.lat.deg])
    inside_mask = hull_path.contains_points(T_pix_lonlat).reshape(ny_t, nx_t)

    # Use the interpolated T_mean (NaN-filled in the dem stage) so the hull
    # stats match the prototype notebook. Pixels still NaN here are off-ROI.
    valid_mask = inside_mask & np.isfinite(T_mean_carr.data)
    T_inside_values = T_mean_carr.data[valid_mask]
    T_inside_lon = world.lon.deg.reshape(ny_t, nx_t)[valid_mask]
    T_inside_lat = world.lat.deg.reshape(ny_t, nx_t)[valid_mask]
    hull_npix_total = int(inside_mask.sum())
    hull_npix_valid = int(valid_mask.sum())
    hull_npix_nan = hull_npix_total - hull_npix_valid
    log.info("T_mean pixels inside hull: %d (total in hull: %d, NaN excluded: %d)",
             hull_npix_valid, hull_npix_total, hull_npix_nan)

    np.savez(out["T_in_hull_npz"],
             values=T_inside_values, lon=T_inside_lon, lat=T_inside_lat)
    log.info("saved %s", out["T_in_hull_npz"])

    T_mean_hull_avg = float(np.nanmean(T_inside_values)) if T_inside_values.size else float("nan")

    # ---- 6. AIA + ADAPT contours + T_mean overlay ---------------------------
    mag_levels = overlay.make_mag_contour_levels(adapt_map)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(projection=aia_carr_sub)
    aia_carr_sub.plot(axes=ax)
    bound = ax.axis()

    cset = adapt_sub.draw_contours(levels=mag_levels, axes=ax, cmap="seismic", alpha=0.4)
    Tmap = T_mean_carr.plot(axes=ax, cmap="inferno", alpha=0.6, origin="lower")
    ax.axis(bound)

    ft = ax.scatter_coord(
        dominant_ft_diffrot, alpha=0.9, marker="X", s=20,
        vmin=0, vmax=float(tbl["prob"].max()),
        c=tbl["prob"], cmap="YlOrRd",
        edgecolors="k", linewidths=0.5, zorder=6,
        label=f"Footpoints on solar surface (prob>{cfg['irap']['prob_threshold_pct']}%)",
    )
    hull_coords = SkyCoord(hull_lons * u.deg, hull_lats * u.deg, frame=aia_carr_sub.coordinate_frame)
    ax.plot_coord(hull_coords, color="lime", linewidth=1.5, linestyle="--",
                 zorder=7, alpha=0.5, label="Dominant footpoints boundary")

    # White dots at hull-interior pixels in AIA submap pixel coords
    ny_aia, nx_aia = aia_carr_sub.data.shape
    y_aia_idx, x_aia_idx = np.mgrid[0:ny_aia, 0:nx_aia]
    aia_world = aia_carr_sub.wcs.pixel_to_world(x_aia_idx.ravel(), y_aia_idx.ravel())
    aia_inside = hull_path.contains_points(
        np.column_stack([aia_world.lon.deg, aia_world.lat.deg])
    ).reshape(ny_aia, nx_aia)
    in_y, in_x = np.where(aia_inside)
    ax.scatter(in_x, in_y, marker=".", s=1, color="white", linewidths=1, zorder=8)

    ax.legend(loc="upper right", fontsize=8, framealpha=0.7)
    fig.colorbar(cset, ax=ax, fraction=0.025, shrink=0.8, pad=0.08,
                label=f"Magnetic Field Strength [{adapt_map.unit}]",
                ticks=list(mag_levels.value) + [0])
    fig.colorbar(Tmap, ax=ax, fraction=0.025, shrink=0.8, pad=0.015, label="Temperature (K)")
    ax.set_title(f"{aia_carr_sub.date.strftime('%Y-%m-%d %H:%M:%S')}  "
                 f"(hull pixels: {hull_npix_valid}/{hull_npix_total}, NaN: {hull_npix_nan})")
    fig.savefig(out["overlay_png"], dpi=cfg["plots"]["dpi"], bbox_inches="tight")
    plt.close(fig)
    log.info("saved overlay %s", out["overlay_png"])

    # ---- 7. Append summary CSV ---------------------------------------------
    row = {
        "event_id": layout.event_id,
        "sc_time_utc": cfg["irap"]["spacecraft_time"],
        "target_time_utc": layout.target_time.iso,
        "spacecraft": cfg["irap"]["spacecraft"],
        "coronal_model": cfg["irap"]["coronal_model"],
        "mode": cfg["irap"]["mode"],
        "region_type": cfg["dem"]["region_type"],
        "roi_blx": roi["bottom_left_arcsec"][0],
        "roi_bly": roi["bottom_left_arcsec"][1],
        "roi_trx": roi["top_right_arcsec"][0],
        "roi_try": roi["top_right_arcsec"][1],
        "prob_threshold_pct": cfg["irap"]["prob_threshold_pct"],
        "n_dominant_footpoints": int(len(tbl)),
        "hull_npix_total": hull_npix_total,
        "hull_npix_valid": hull_npix_valid,
        "hull_npix_nan": hull_npix_nan,
        "T_mean_hull_avg": T_mean_hull_avg,
        "completed_at": datetime.now(timezone.utc).isoformat(),
    }
    _append_summary_row(out["summary_csv"], row)
    log.info("appended summary row -> %s (T_mean_hull_avg=%.4g K)",
             out["summary_csv"], T_mean_hull_avg)

    mfst.update_stage(layout.manifest_path, "extract",
                     {**{k: str(v) for k, v in out.items()}, **row})
    return out
