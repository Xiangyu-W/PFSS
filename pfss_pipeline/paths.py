"""Deterministic path resolution for all pipeline artifacts."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from astropy.time import Time


def event_id(cfg: dict) -> str:
    sct = Time(cfg["irap"]["spacecraft_time"]).strftime("%Y%m%dT%H%M%S")
    return f"{cfg['irap']['spacecraft']}_{sct}_{cfg['irap']['coronal_model']}"


def stamp(t: Time) -> str:
    return t.strftime("%Y%m%d_%H%M%S")


@dataclass
class OutputLayout:
    cfg: dict
    target_time: Time | None = None  # set after irap_fetch

    @property
    def event_id(self) -> str:
        return event_id(self.cfg)

    @property
    def event_short_id(self) -> str:
        """`{spacecraft}_{spacecraft_time}` — used for ROI-helper figure filenames."""
        sct = Time(self.cfg["irap"]["spacecraft_time"]).strftime("%Y%m%dT%H%M%S")
        return f"{self.cfg['irap']['spacecraft']}_{sct}"

    @property
    def event_dir(self) -> Path:
        return Path(self.cfg["results_root"]) / self.event_id

    # ---- Cross-event ROI-helper figure aggregation ----
    @property
    def figures_for_roi_dir(self) -> Path:
        return Path(self.cfg["results_root"]) / "figures_for_roi"

    @property
    def irap_footpoints_full_path(self) -> Path:
        return self.figures_for_roi_dir / "irap" / f"{self.event_short_id}_footpoints_on_ADAPT_full.png"

    @property
    def irap_footpoints_submap_path(self) -> Path:
        return self.figures_for_roi_dir / "irap" / f"{self.event_short_id}_footpoints_on_ADAPT_submap.png"

    def aia_fulldisk_for_roi_path(self, wl: str) -> Path:
        return self.figures_for_roi_dir / "aia" / f"{self.event_short_id}_{wl}A.png"

    @property
    def manifest_path(self) -> Path:
        return self.event_dir / "run_manifest.yaml"

    # ---- AIA prep (shared cache, not under event_dir) ----
    @property
    def aia_prep_dir(self) -> Path:
        return Path(self.cfg["aia_prep_dir"])

    @property
    def aia_prep_images_dir(self) -> Path:
        return self.aia_prep_dir / "images"

    def aia_prep_fits(self, wl: str, t: Time) -> Path:
        return self.aia_prep_dir / f"aia_prep_{wl}A_{stamp(t)}.fits"

    def aia_prep_png(self, wl: str, t: Time) -> Path:
        return self.aia_prep_images_dir / f"aia_prep_comparison_{wl}A_{stamp(t)}.png"

    # ---- IRAP ----
    @property
    def irap_dir(self) -> Path:
        return self.event_dir / "irap"

    @property
    def irap_zip_dir(self) -> Path:
        return self.irap_dir / "mct_raw"

    @property
    def irap_zip_path(self) -> Path:
        ir = self.cfg["irap"]
        sct = Time(ir["spacecraft_time"]).strftime("%Y%m%dT%H%M%S")
        return self.irap_zip_dir / f"{ir['spacecraft']}_PARKER_PFSS_{ir['mode']}_{ir['coronal_model']}_{sct}.zip"

    @property
    def irap_figures_dir(self) -> Path:
        return self.irap_dir / "figures"

    @property
    def irap_dominant_ft_path(self) -> Path:
        thr = self.cfg["irap"]["prob_threshold_pct"]
        return self.irap_dir / f"dominant_ft_prob{thr}.ecsv"

    # ---- DEM ----
    @property
    def dem_dir(self) -> Path:
        return self.event_dir / "dem"

    @property
    def dem_figures_dir(self) -> Path:
        return self.dem_dir / "figures"

    def _need_t(self) -> Time:
        if self.target_time is None:
            raise RuntimeError("target_time not set; run stage_irap_fetch first or pass --target-time.")
        return self.target_time

    def dem_t_mean(self) -> Path:
        return self.dem_dir / f"T_mean_map_{stamp(self._need_t())}.fits"

    def dem_t_mean_raw(self) -> Path:
        return self.dem_dir / f"T_mean_raw_map_{stamp(self._need_t())}.fits"

    def dem_t_peak(self) -> Path:
        return self.dem_dir / f"T_peak_map_{stamp(self._need_t())}.fits"

    def dem_em(self) -> Path:
        return self.dem_dir / f"EM_map_{stamp(self._need_t())}.fits"

    def dem_cube(self) -> Path:
        return self.dem_dir / f"dem_cube_{stamp(self._need_t())}.npz"

    # ---- Extract ----
    @property
    def extract_dir(self) -> Path:
        return self.event_dir / "extract"

    def extract_t_mean_in_roi(self) -> Path:
        return self.extract_dir / f"T_mean_in_roi_{stamp(self._need_t())}.fits"

    def extract_t_in_hull_npz(self) -> Path:
        sc = self.cfg["irap"]["spacecraft"]
        return self.extract_dir / f"{sc}_T_inside_hull_{stamp(self._need_t())}.npz"

    def extract_t_peak_in_roi(self) -> Path:
        return self.extract_dir / f"T_peak_in_roi_{stamp(self._need_t())}.fits"

    def extract_em_in_roi(self) -> Path:
        return self.extract_dir / f"EM_in_roi_{stamp(self._need_t())}.fits"

    def extract_overlay_png(self) -> Path:
        return self.extract_dir / f"DEM_on_AIA_overlay_{stamp(self._need_t())}.png"

    @property
    def extract_summary_csv(self) -> Path:
        return self.extract_dir / "footpoint_temperature_summary.csv"

    @property
    def extract_summary_csv_latest(self) -> Path:
        return self.extract_dir / "footpoint_temperature_summary_latest.csv"

    def ensure_dirs(self) -> None:
        for d in [
            self.event_dir,
            self.aia_prep_dir,
            self.aia_prep_images_dir,
            self.irap_dir,
            self.irap_zip_dir,
            self.dem_dir,
            self.dem_figures_dir,
            self.extract_dir,
            self.figures_for_roi_dir / "irap",
            self.figures_for_roi_dir / "aia",
        ]:
            d.mkdir(parents=True, exist_ok=True)
