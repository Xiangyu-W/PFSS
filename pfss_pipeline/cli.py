"""argparse driver for the pipeline."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from astropy.time import Time

from pfss_pipeline import config as cfg_mod
from pfss_pipeline import manifest as mfst
from pfss_pipeline.paths import OutputLayout

STAGES = ("irap", "aia", "dem", "extract", "all-after-roi")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="pfss_pipeline")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--stage", required=True, choices=STAGES)
    p.add_argument("--force", action="store_true", help="bypass skip-if-exists")
    p.add_argument("--force-stage", choices=STAGES[:-1], help="re-run only this stage when --stage all-after-roi")
    p.add_argument("--target-time", help="ISO override for derived target_time")
    p.add_argument("--results-root", type=Path)
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--log-level", default=None, choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p


def _resolve_target_time(args, cfg: dict, layout: OutputLayout) -> Time | None:
    if args.target_time:
        return Time(args.target_time)
    m = mfst.read(layout.manifest_path)
    tt = mfst.get_target_time(m)
    return Time(tt) if tt else None


def _force_for(stage: str, args) -> bool:
    if args.force:
        return True
    if args.stage == "all-after-roi" and args.force_stage == stage:
        return True
    return False


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    cfg = cfg_mod.load_config(args.config)
    if args.results_root:
        cfg["results_root"] = str(args.results_root)
    if args.log_level:
        cfg["runtime"]["log_level"] = args.log_level
    logging.basicConfig(level=cfg["runtime"]["log_level"], format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    log = logging.getLogger("pfss_pipeline.cli")

    layout = OutputLayout(cfg)
    layout.target_time = _resolve_target_time(args, cfg, layout)

    if args.dry_run:
        return _dry_run_report(args, cfg, layout)

    layout.ensure_dirs()

    # Lazy imports so each stage's heavy deps load only when needed.
    if args.stage == "irap":
        from pfss_pipeline.stages import irap_fetch
        irap_fetch.run(cfg, layout, force=_force_for("irap", args))
    elif args.stage == "aia":
        _require_target_time(layout)
        from pfss_pipeline.stages import aia_prep
        aia_prep.run(cfg, layout, force=_force_for("aia", args))
    elif args.stage == "dem":
        _require_target_time(layout)
        cfg_mod.assert_dem_ready(cfg)
        from pfss_pipeline.stages import dem
        dem.run(cfg, layout, force=_force_for("dem", args))
    elif args.stage == "extract":
        _require_target_time(layout)
        cfg_mod.assert_dem_ready(cfg)
        from pfss_pipeline.stages import extract
        extract.run(cfg, layout, force=_force_for("extract", args))
    elif args.stage == "all-after-roi":
        cfg_mod.assert_dem_ready(cfg)
        _require_target_time(layout)
        from pfss_pipeline.stages import aia_prep, dem, extract
        aia_prep.run(cfg, layout, force=_force_for("aia", args))
        dem.run(cfg, layout, force=_force_for("dem", args))
        extract.run(cfg, layout, force=_force_for("extract", args))

    log.info("done.")
    return 0


def _require_target_time(layout: OutputLayout) -> None:
    if layout.target_time is None:
        raise SystemExit(
            "target_time not in manifest; run --stage irap first or pass --target-time."
        )


def _dry_run_report(args, cfg: dict, layout: OutputLayout) -> int:
    print(f"event_id        : {layout.event_id}")
    print(f"event_dir       : {layout.event_dir}")
    print(f"manifest        : {layout.manifest_path}")
    print(f"target_time     : {layout.target_time}")
    print(f"aia_prep_dir    : {layout.aia_prep_dir}")
    print(f"irap zip path   : {layout.irap_zip_path}")
    print(f"irap dom_ft     : {layout.irap_dominant_ft_path}")
    if layout.target_time is not None:
        print(f"dem T_mean      : {layout.dem_t_mean()}")
        print(f"extract summary : {layout.extract_summary_csv}")
    return 0
