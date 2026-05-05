"""stage_irap_fetch: IRAP MCT fetch + parse + footpoint figure + derive target_time.

This stage runs first. It:
  - obtains the MCT ZIP (cache or Selenium-triggered fetch)
  - parses fileparam / connectivity / fieldline / hcs into pandas / yaml objects
  - persists each parsed product to disk (parquet + yaml)
  - downloads the ADAPT map closest to params.metadata.date_surf
  - derives target_time = ADAPT map date and writes it to the manifest
  - saves footpoints-on-ADAPT figures (full disk + submap)
  - saves dominant-footpoint ECSV for downstream stages
"""
from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml
from astropy.time import Time

from pfss_pipeline import manifest as mfst
from pfss_pipeline.irap import adapt as adapt_mod
from pfss_pipeline.irap import footpoints as fp_mod
from pfss_pipeline.irap import mct
from pfss_pipeline.irap import overlay

log = logging.getLogger(__name__)


def _date_time_strs(sct_iso: str) -> tuple[str, str]:
    t = Time(sct_iso)
    return t.strftime("%Y-%m-%d"), t.strftime("%H%M%S")


def run(cfg: dict, layout, force: bool = False) -> dict:
    ir = cfg["irap"]
    date_str, time_str = _date_time_strs(ir["spacecraft_time"])

    layout.ensure_dirs()
    irap_dir = layout.irap_dir
    fig_dir = layout.irap_figures_dir

    # ---- 1. ZIP (cache or fetch) ---------------------------------------
    zip_files = mct.cached_zip_or_fetch(
        zip_path=layout.irap_zip_path,
        sc=ir["spacecraft"], coronal=ir["coronal_model"], mode=ir["mode"],
        date_str=date_str, time_str=time_str,
        chrome_binary=ir["selenium_chrome_binary"],
        force=force,
    )

    # ---- 2. Parse all four products -----------------------------------
    fileparam_bytes = mct.get_ascii_text(zip_files, "fileparam")
    if fileparam_bytes is None:
        raise RuntimeError(f"no fileparam in MCT ZIP: {layout.irap_zip_path}")
    params = mct.parse_params(fileparam_bytes)

    fieldline_bytes = mct.get_ascii_text(zip_files, "filefieldline")
    fieldlines_df = mct.parse_fieldline_fits(fieldline_bytes) if fieldline_bytes else None

    conn_bytes = mct.get_ascii_text(zip_files, "fileconnectivity")
    if conn_bytes is None:
        raise RuntimeError(f"no fileconnectivity in MCT ZIP: {layout.irap_zip_path}")
    _, header_df, foot_solarsurf_df = mct.parse_connectivity(conn_bytes)

    hcs_bytes = mct.get_ascii_text(zip_files, "filehcs")
    hcs_df = mct.parse_hcs(hcs_bytes) if hcs_bytes else None

    # ---- 3. Persist parsed products -----------------------------------
    fileparam_path = irap_dir / "fileparam.yaml"
    fileparam_path.write_text(yaml.safe_dump(params, sort_keys=False))

    fieldlines_path = irap_dir / "fieldlines.parquet"
    if fieldlines_df is not None:
        fieldlines_df.to_parquet(fieldlines_path)

    conn_path = irap_dir / "connectivity.parquet"
    foot_solarsurf_df.to_parquet(conn_path)

    hcs_path = irap_dir / "hcs.parquet"
    if hcs_df is not None:
        hcs_df.to_parquet(hcs_path)

    # ---- 4. Locate ADAPT closest to params.metadata.date_surf ---------
    pfss_time_dt = datetime.strptime(params["metadata"]["date_surf"], "%Y-%m-%d %H:%M:%S.%f")
    matched_dt, adapt_path, adapt_fname, adapt_url = adapt_mod.find_closest_adapt(
        adapt_mod.GONG_BASE, pfss_time_dt, ir["adapt_cache_dir"],
    )
    realization = ir["realization_adapt"]
    adapt_data, adapt_hdr, adapt_map = adapt_mod.load_adapt(adapt_path, realization=realization)
    log.info("ADAPT matched %s (Δt = %+.2f h)", matched_dt,
             (pfss_time_dt - matched_dt).total_seconds() / 3600.0)

    # ---- 5. Derive target_time = ADAPT map date -----------------------
    target_time = Time(adapt_map.date)
    layout.target_time = target_time
    log.info("derived target_time = %s", target_time.iso)

    # ---- 6. Build footpoint SkyCoords ---------------------------------
    sw_type = fp_mod.select_sw_type(foot_solarsurf_df)
    type_df = foot_solarsurf_df[foot_solarsurf_df["type"] == sw_type]
    irap_foot = fp_mod.make_skycoord_carr(
        type_df, obstime=params["metadata"]["date_surf"],
        observer_coord=adapt_map.observer_coordinate,
    )
    irap_foot_diffrot = fp_mod.diff_rotate(irap_foot, rotated_time=adapt_map.date)

    irap_hcs_diffrot = None
    if hcs_df is not None and len(hcs_df):
        from astropy import units as u  # local import; avoids polluting module top
        from astropy.coordinates import SkyCoord
        from sunpy.coordinates import HeliographicCarrington
        irap_hcs = SkyCoord(
            hcs_df["lon_deg"].values * u.deg,
            hcs_df["lat_deg"].values * u.deg,
            1 * u.Rsun,
            obstime=params["metadata"]["date_surf"],
            frame=HeliographicCarrington,
            observer=adapt_map.observer_coordinate,
        )
        irap_hcs_diffrot = fp_mod.diff_rotate(irap_hcs, rotated_time=adapt_map.date)

    # ---- 7. Save figures (full + submap) ------------------------------
    fig_full = layout.irap_figures_dir / "footpoints_on_ADAPT_full.png"
    fig_sub = layout.irap_figures_dir / "footpoints_on_ADAPT_submap.png"
    if force or not fig_full.exists():
        fig = overlay.plot_footpoints_on_adapt(
            adapt_map, irap_foot_diffrot, type_df["prob"].values,
            hcs_diffrot=irap_hcs_diffrot,
            title=f"{ir['spacecraft']} IRAP {ir['coronal_model']} footpoints @ {date_str}T{time_str}",
        )
        fig.savefig(fig_full, dpi=cfg["plots"]["dpi"], bbox_inches="tight")
        plt.close(fig)
        log.info("saved %s", fig_full)
    if force or not fig_sub.exists():
        fig = overlay.plot_footpoints_on_adapt(
            adapt_map, irap_foot_diffrot, type_df["prob"].values,
            hcs_diffrot=irap_hcs_diffrot,
            carrington_roi=ir["carrington_roi"],
            title=f"{ir['spacecraft']} IRAP {ir['coronal_model']} footpoints (submap)",
        )
        fig.savefig(fig_sub, dpi=cfg["plots"]["dpi"], bbox_inches="tight")
        plt.close(fig)
        log.info("saved %s", fig_sub)

    # ---- 8. Dominant footpoints ECSV ----------------------------------
    dominant_df = fp_mod.filter_dominant(
        foot_solarsurf_df, sw_type, prob_threshold_pct=ir["prob_threshold_pct"],
    )
    dominant_ft = fp_mod.make_skycoord_carr(
        dominant_df, obstime=params["metadata"]["date_surf"],
        observer_coord=adapt_map.observer_coordinate,
    )
    fp_mod.save_dominant_ecsv(
        layout.irap_dominant_ft_path, dominant_ft, dominant_df, params,
        sw_type=sw_type, prob_threshold_pct=ir["prob_threshold_pct"],
    )

    # ---- 9. Manifest update -------------------------------------------
    paths = {
        "zip": str(layout.irap_zip_path),
        "fileparam": str(fileparam_path),
        "fieldlines": str(fieldlines_path) if fieldlines_df is not None else None,
        "connectivity": str(conn_path),
        "hcs": str(hcs_path) if hcs_df is not None else None,
        "adapt": str(adapt_path),
        "adapt_url": adapt_url,
        "adapt_fname": adapt_fname,
        "dominant_ft": str(layout.irap_dominant_ft_path),
        "figure_full": str(fig_full),
        "figure_submap": str(fig_sub),
        "target_time": target_time.iso,
        "sw_type": sw_type,
        "n_dominant": int(len(dominant_df)),
        "n_total_solarsurf": int(len(type_df)),
        "date_surf": params["metadata"]["date_surf"],
        "date_insitu": params["metadata"].get("date_insitu", ""),
    }
    mfst.update_stage(layout.manifest_path, "irap_fetch", paths)
    log.info("stage_irap_fetch complete; target_time=%s, %d dominant footpoints",
             target_time.iso, len(dominant_df))
    return {k: Path(v) if v and isinstance(v, str) and v.startswith("/") else v
            for k, v in paths.items()}
