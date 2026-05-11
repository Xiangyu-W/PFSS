"""Batch driver for pfss_pipeline over a date range at 6h cadence.

Workflow:
  1. python scripts/batch_run.py gen-configs                  # write configs/event_*.yaml
  2. python scripts/batch_run.py irap --workers 4             # fetch IRAP + ADAPT
  3. (manual) fill extract.overlay_carrington_roi in each YAML
  4. python scripts/batch_run.py aia-fetch --workers 1        # JSOC L1 download (serial)
  5. python scripts/batch_run.py aia-prep --workers 8         # PSF + register (CPU-bound)
  6. (manual) fill dem.roi / dem.region_type in each YAML
  7. python scripts/batch_run.py all-after-roi --workers 8    # DEM + extract

Use `aia` (combined) instead of `aia-fetch` + `aia-prep` if you don't need to split.

Recommended worker caps:
  irap          : <=4   Selenium/Chrome + IRAP server rate limits
  aia-fetch     : 1     JSOC export-hash collisions when concurrent
  aia-prep      : <=20  pure CPU, BLAS pinned to 1 thread/worker
  all-after-roi : <=8   DEM is CPU-bound; same BLAS pinning applies
"""
from __future__ import annotations

import argparse
import copy
import logging
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = REPO_ROOT / "configs"
TEMPLATE = CONFIG_DIR / "event_2022-03-03.yaml"

START = "2022-03-02 06:00"
END = "2022-03-07 00:00"
FREQ = "6h"

log = logging.getLogger("batch_run")


def event_times() -> list[pd.Timestamp]:
    return list(pd.date_range(START, END, freq=FREQ))


def cfg_path_for(t: pd.Timestamp) -> Path:
    return CONFIG_DIR / f"event_{t:%Y-%m-%dT%H%M}.yaml"


def gen_configs(force: bool) -> None:
    base = yaml.safe_load(TEMPLATE.read_text())
    times = event_times()
    log.info("generating %d configs from %s", len(times), TEMPLATE.name)
    for t in times:
        out = cfg_path_for(t)
        if out.exists() and not force:
            log.info("skip existing %s", out.name)
            continue
        cfg = copy.deepcopy(base)
        cfg["irap"]["spacecraft_time"] = t.strftime("%Y-%m-%dT%H:%M:%S")
        cfg["dem"]["roi"] = {"bottom_left_arcsec": None, "top_right_arcsec": None}
        cfg["extract"]["overlay_carrington_roi"] = {"lon": None, "lat": None}
        out.write_text(yaml.safe_dump(cfg, sort_keys=False))
        log.info("wrote %s", out.name)


def _resolve_results_root(cfg_path: Path) -> Path:
    cfg = yaml.safe_load(cfg_path.read_text())
    return Path(cfg.get("results_root") or REPO_ROOT / "results")


def _run_one(stage: str, idx: int, total: int, t: pd.Timestamp,
             force: bool = False) -> tuple[pd.Timestamp, int, str]:
    cfg = cfg_path_for(t)
    if not cfg.exists():
        return (t, 2, f"missing {cfg.name}")

    results_root = _resolve_results_root(cfg)
    log_dir = results_root / "_batch_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%dT%H%M%S")
    log_path = log_dir / f"{cfg.stem}_{stage}_{ts}.log"

    env = dict(os.environ)
    # Pin BLAS/OMP to 1 thread per worker so multiple parallel DEM workers
    # don't oversubscribe the CPU.
    env.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    })
    # The conda `selenium-manager` package sets SE_MANAGER_PATH via its
    # activate.d hook, which doesn't fire for non-interactive subprocesses.
    # Re-derive it from sys.executable so Selenium can find the binary.
    if "SE_MANAGER_PATH" not in env:
        candidate = Path(sys.executable).parent / "selenium-manager"
        if candidate.exists():
            env["SE_MANAGER_PATH"] = str(candidate)
    # `base` conda env on this host exports SSL_CERT_FILE pointing to a 2019
    # system bundle without modern Let's Encrypt roots, breaking JSOC HTTPS.
    # Always force certifi's bundle, which conda's DEM activate.d hook also
    # would set if conda activate had been used.
    try:
        import certifi
        env["SSL_CERT_FILE"] = certifi.where()
        env["REQUESTS_CA_BUNDLE"] = certifi.where()
    except ImportError:
        pass

    cmd = [sys.executable, "-m", "pfss_pipeline", "--config", str(cfg), "--stage", stage]
    if force:
        cmd.append("--force")
    log.info("[%d/%d] START %s -> %s", idx, total, cfg.name, log_path.name)
    with open(log_path, "w") as fh:
        fh.write(f"# cmd: {' '.join(cmd)}\n# cwd: {REPO_ROOT}\n# started: {ts}\n\n")
        fh.flush()
        rc = subprocess.run(cmd, cwd=REPO_ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
    return (t, rc, str(log_path))


def _select_times(start: int, end: int | None,
                  indices: list[int] | None = None) -> list[tuple[int, pd.Timestamp]]:
    """Return [(absolute_idx_1based, time), ...].

    If `indices` is given (1-based, any order), it takes precedence over start/end
    and the result preserves the user-supplied order with duplicates removed.
    Otherwise the result is a contiguous slice from `start` to `end` inclusive.
    """
    times = event_times()
    n = len(times)
    if indices is not None:
        seen: set[int] = set()
        ordered: list[int] = []
        for i in indices:
            if i < 1 or i > n:
                raise SystemExit(f"--indices out of range; valid: 1..{n}, got {i}")
            if i not in seen:
                seen.add(i)
                ordered.append(i)
        return [(i, times[i - 1]) for i in ordered]
    end = end if end is not None else n
    if start < 1 or start > n or end < start or end > n:
        raise SystemExit(f"--from/--to out of range; valid: 1..{n}, got from={start} to={end}")
    return [(i, times[i - 1]) for i in range(start, end + 1)]


def run_stage(stage: str, workers: int, start: int = 1,
              end: int | None = None, force: bool = False,
              indices: list[int] | None = None) -> tuple[int, list[pd.Timestamp]]:
    """Returns (rc, list_of_successful_event_times)."""
    selected = _select_times(start, end, indices)
    total = len(selected)
    sel_idx = [i for (i, _) in selected]
    log.info("running stage=%s on events %s (%d total) with %d workers (force=%s)",
             stage, sel_idx, total, workers, force)

    successes: list[pd.Timestamp] = []
    failures: list[tuple[pd.Timestamp, int, str]] = []
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_run_one, stage, abs_i, total, t, force): t
                   for (abs_i, t) in selected}
        for fut in as_completed(futures):
            t, rc, info = fut.result()
            done += 1
            if rc == 0:
                log.info("[%d/%d] OK    %s", done, total, t)
                successes.append(t)
            else:
                log.error("[%d/%d] FAIL  %s rc=%d  see %s", done, total, t, rc, info)
                failures.append((t, rc, info))

    log.info("=" * 70)
    log.info("done; %d/%d failures", len(failures), total)
    for t, rc, info in failures:
        log.error("  %s rc=%d log=%s", t, rc, info)
    return (1 if failures else 0, successes)


def _save_aia_preview(t: pd.Timestamp, wavelength: str) -> Path:
    """Make a full-disk AIA PNG with HPC arcsec axes for ROI determination."""
    import matplotlib  # lazy
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy import units as u
    from sunpy.map import Map

    cfg_path = cfg_path_for(t)
    cfg = yaml.safe_load(cfg_path.read_text())
    results_root = Path(cfg.get("results_root") or REPO_ROOT / "results")

    manifest_paths = sorted(results_root.glob(f"*{t:%Y%m%dT%H%M}*/run_manifest.yaml"))
    if not manifest_paths:
        raise FileNotFoundError(f"no run_manifest for {t}")
    mf = yaml.safe_load(manifest_paths[0].read_text())
    prep_path = mf.get("stages", {}).get("aia_prep", {}).get(wavelength)
    if not prep_path:
        raise KeyError(f"manifest has no aia_prep[{wavelength}] for {t}")

    event_dir = manifest_paths[0].parent
    out_dir = event_dir / "aia"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"aia_{wavelength}A_fulldisk_for_roi.png"

    m = Map(prep_path)
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(projection=m)
    m.plot(axes=ax, clip_interval=(1, 99.95) * u.percent)
    ax.set_title(f"{m.wavelength.value:.0f} Å  {m.date.iso[:19]}  ({event_dir.name})")
    ax.coords.grid(color="white", alpha=0.3, linestyle="--")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run_aia_with_preview(stage: str, workers: int, preview_wl: str,
                         start: int = 1, end: int | None = None,
                         force: bool = False,
                         indices: list[int] | None = None) -> int:
    rc, successes = run_stage(stage, workers, start, end, force, indices)
    log.info("=" * 70)
    log.info("generating %s Å full-disk preview PNGs for %d events", preview_wl, len(successes))
    n_ok, n_fail = 0, 0
    for t in successes:
        try:
            png = _save_aia_preview(t, preview_wl)
            log.info("preview: %s", png)
            n_ok += 1
        except Exception as e:
            log.warning("preview FAILED for %s: %s", t, e)
            n_fail += 1
    log.info("preview: %d ok, %d fail", n_ok, n_fail)
    return rc


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("gen-configs", help="write per-event YAML configs")
    g.add_argument("--force", action="store_true", help="overwrite existing configs")

    def _parse_indices(s: str) -> list[int]:
        return [int(tok) for tok in s.replace(",", " ").split() if tok]

    def _add_common(p):
        p.add_argument("--from", dest="start", type=int, default=1,
                       help="start event index, 1-based (default 1)")
        p.add_argument("--to", dest="end", type=int, default=None,
                       help="end event index, 1-based inclusive (default last)")
        p.add_argument("--indices", type=_parse_indices, default=None,
                       help="explicit 1-based event indices, comma- or space-separated "
                            "(e.g. '2,6,10,14,18'); overrides --from/--to")
        p.add_argument("--force", action="store_true",
                       help="pass --force through to pipeline (bypass skip-if-exists)")

    for name in ("irap", "aia-fetch", "all-after-roi"):
        s = sub.add_parser(name, help=f"run --stage {name} for every event")
        s.add_argument("--workers", type=int, default=1, help="concurrent workers (default 1)")
        _add_common(s)

    pp = sub.add_parser("aia-prep", help="run --stage aia-prep + save preview PNG (CPU-bound)")
    pp.add_argument("--workers", type=int, default=1, help="concurrent workers (default 1)")
    pp.add_argument("--preview-wl", default="193", help="wavelength for preview PNG (default 193)")
    _add_common(pp)

    a = sub.add_parser("aia", help="run --stage aia (fetch + prep) + save preview PNG")
    a.add_argument("--workers", type=int, default=1, help="concurrent workers (default 1)")
    a.add_argument("--preview-wl", default="193", help="wavelength for preview PNG (default 193)")
    _add_common(a)

    args = p.parse_args(argv)
    logging.basicConfig(level="INFO", format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    if args.cmd == "gen-configs":
        gen_configs(force=args.force)
        return 0
    if args.cmd in ("irap", "aia-fetch", "all-after-roi"):
        rc, _ = run_stage(args.cmd, args.workers, args.start, args.end,
                          args.force, args.indices)
        return rc
    if args.cmd in ("aia", "aia-prep"):
        return run_aia_with_preview(args.cmd, args.workers, args.preview_wl,
                                     args.start, args.end, args.force, args.indices)
    return 1


if __name__ == "__main__":
    sys.exit(main())
