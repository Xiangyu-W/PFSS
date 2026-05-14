"""Phase-A verification: compare drms `protocol="as-is"` vs Fido `protocol="fits"`.

Downloads the same AIA L1 record two ways, times each, compares FITS headers
and pixel arrays, and tries aiapy.calibrate.register on the as-is map (which
is what aia_prep does first).

Pass criteria printed at end:
  - data arrays match (allclose)
  - as-is header has the pointing keys needed by aiapy.calibrate.register
  - register() succeeds without raising

Usage:
  python scripts/verify_drms_as_is.py \
      --event configs/event_2022-03-03T0600.yaml --wl 171
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]

# Force certifi CA bundle (matches __main__.py behavior).
try:
    import certifi
    os.environ["SSL_CERT_FILE"] = certifi.where()
    os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()
except ImportError:
    pass

log = logging.getLogger("verify_drms_as_is")


def _event_id_from_cfg(cfg: dict) -> str:
    from astropy.time import Time
    sct = Time(cfg["irap"]["spacecraft_time"]).strftime("%Y%m%dT%H%M%S")
    return f"{cfg['irap']['spacecraft']}_{sct}_{cfg['irap']['coronal_model']}"


def _load_target_time(cfg_path: Path):
    from astropy.time import Time
    cfg = yaml.safe_load(cfg_path.read_text())
    results_root = Path(cfg["results_root"])
    event_id = _event_id_from_cfg(cfg)
    manifest = results_root / event_id / "run_manifest.yaml"
    if not manifest.exists():
        raise SystemExit(f"manifest missing: {manifest}; run --stage irap first")
    m = yaml.safe_load(manifest.read_text())
    tt = m["stages"]["irap_fetch"]["target_time"]
    return Time(tt), cfg


def fetch_drms_as_is(target_time, wl: str, out_dir: Path, jsoc_notify: str):
    """drms as-is + url_quick: synchronous, no export queue."""
    import drms
    c = drms.Client(email=jsoc_notify)
    # Time window: 1 minute centered on target_time, in TAI.
    # JSOC indexes aia.lev1_euv_12s by T_REC in TAI; AIA cadence is 12s so a
    # 1-minute window guarantees ~5 records — we'll grab them all and keep one.
    t0 = (target_time - 30 / 86400).strftime("%Y.%m.%d_%H:%M:%S")
    query = f"aia.lev1_euv_12s[{t0}_TAI/1m][{wl}]{{image}}"
    log.info("drms query: %s", query)
    req = c.export(query, protocol="as-is", method="url_quick")
    log.info("export status=%s n_records=%d", req.status, len(req.urls))
    out_dir.mkdir(parents=True, exist_ok=True)
    df = req.download(str(out_dir))
    # df has columns including 'download' (local path)
    paths = [p for p in df["download"].tolist() if p and Path(p).exists()]
    if not paths:
        raise RuntimeError(f"drms as-is returned no files; df={df}")
    return paths


def fetch_fido_fits(target_time, wl: str, out_dir: Path, jsoc_notify: str):
    """Current pipeline path: Fido + a.jsoc + export queue (protocol=fits)."""
    from astropy.time import TimeDelta
    from sunpy.net import Fido, attrs as a
    out_dir.mkdir(parents=True, exist_ok=True)
    q = Fido.search(
        a.Time(target_time, target_time + TimeDelta(60, format="sec")),
        a.jsoc.Series("aia.lev1_euv_12s"),
        a.jsoc.PrimeKey("WAVELNTH", wl),
        a.jsoc.Segment("image"),
        a.jsoc.Notify(jsoc_notify),
    )
    log.info("Fido search returned %d records", sum(len(r) for r in q))
    fetched = Fido.fetch(q, path=str(out_dir / "{file}"))
    paths = [str(p) for p in fetched if Path(p).exists()]
    if not paths:
        raise RuntimeError(f"Fido returned no files; fetched={fetched}")
    return paths


def _pick_closest(paths, target_time) -> str:
    """Pick the FITS file whose name's timestamp is closest to target_time.

    Handles both filename conventions:
      - dashed (Fido):  aia.lev1_euv_12s.2022-03-01T071647Z.171.image_lev1.fits
      - dashless (drms): aia.lev1_euv_12s.20220301T071647Z.171.image_lev1.fits
    """
    def fname_time(p):
        s = str(p)
        m = re.search(r"aia\.lev1_euv_12s\.(\d{4}-\d{2}-\d{2}T\d{6})Z", s)
        if m:
            return datetime.strptime(m.group(1), "%Y-%m-%dT%H%M%S")
        m = re.search(r"aia\.lev1_euv_12s\.(\d{8}T\d{6})Z", s)
        if m:
            return datetime.strptime(m.group(1), "%Y%m%dT%H%M%S")
        raise ValueError(f"can't parse timestamp from {s}")
    target_dt = datetime.strptime(target_time.iso[:19], "%Y-%m-%d %H:%M:%S")
    return str(min(paths, key=lambda p: abs((fname_time(p) - target_dt).total_seconds())))


def inspect_map(path: str, label: str) -> dict:
    from sunpy.map import Map
    m = Map(path)
    h = m.fits_header
    info = {
        "label": label,
        "path": path,
        "size_bytes": Path(path).stat().st_size,
        "date": str(m.date),
        "wavelength": str(m.wavelength),
        "dimensions": tuple(int(x.value) for x in m.dimensions),
        "exptime": h.get("EXPTIME"),
        "t_obs": h.get("T_OBS"),
        "lvl_num": h.get("LVL_NUM"),
        "x0_mp": h.get("X0_MP"),
        "y0_mp": h.get("Y0_MP"),
        "crota2": h.get("CROTA2"),
        "crpix1": h.get("CRPIX1"),
        "crpix2": h.get("CRPIX2"),
        "cdelt1": h.get("CDELT1"),
        "cdelt2": h.get("CDELT2"),
        "n_header_keys": len(h),
    }
    return info, m, h


def try_register(m, label: str) -> tuple[bool, str]:
    try:
        from aiapy.calibrate import register, update_pointing
        m_up = update_pointing(m)
        m_reg = register(m_up)
        return True, f"{label}: register OK, dims={tuple(int(x.value) for x in m_reg.dimensions)}"
    except Exception as e:
        return False, f"{label}: register FAILED: {type(e).__name__}: {e}"


def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--event", type=Path, required=True,
                  help="Path to event YAML (uses its manifest for target_time)")
    p.add_argument("--wl", default="171", help="AIA wavelength (default 171)")
    p.add_argument("--keep-tmp", action="store_true",
                  help="Don't delete the temp dirs after the run")
    p.add_argument("--skip-fido", action="store_true",
                  help="Skip the Fido side-by-side (use when JSOC queue is busy)")
    p.add_argument("--existing", type=Path, default=None,
                  help="Reuse FITS files already in this directory (skips the drms download)")
    p.add_argument("--fido-baseline", type=Path, default=None,
                  help="Compare headers against this Fido-downloaded FITS file (any time)")
    args = p.parse_args(argv)
    logging.basicConfig(level="INFO",
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    target_time, cfg = _load_target_time(args.event)
    jsoc_notify = cfg["aia"]["jsoc_notify"]
    log.info("event=%s target_time=%s wl=%s notify=%s",
             args.event.name, target_time.iso, args.wl, jsoc_notify)

    tmp_root = Path(tempfile.mkdtemp(prefix="verify_drms_"))
    tmp_a = tmp_root / "as_is"
    tmp_b = tmp_root / "fido"
    log.info("tmp_root=%s", tmp_root)

    # ---- A: drms as-is ----
    if args.existing:
        paths_a = sorted(str(p) for p in args.existing.glob("*.fits"))
        if not paths_a:
            raise SystemExit(f"no .fits files in --existing dir {args.existing}")
        dt_a = float("nan")
        log.info("AS-IS  reused %d file(s) from %s", len(paths_a), args.existing)
    else:
        t0 = time.perf_counter()
        paths_a = fetch_drms_as_is(target_time, args.wl, tmp_a, jsoc_notify)
        dt_a = time.perf_counter() - t0
        log.info("AS-IS  done in %.2fs, %d file(s)", dt_a, len(paths_a))
    path_a = _pick_closest(paths_a, target_time)

    # ---- B: Fido protocol=fits (optional) ----
    info_b = hdr_b = map_b = None
    dt_b = None
    if args.fido_baseline:
        info_b, map_b, hdr_b = inspect_map(str(args.fido_baseline), "fido-baseline")
        log.info("FIDO   reused baseline %s", args.fido_baseline)
    elif args.skip_fido:
        log.info("FIDO   skipped (--skip-fido)")
    else:
        t0 = time.perf_counter()
        paths_b = fetch_fido_fits(target_time, args.wl, tmp_b, jsoc_notify)
        dt_b = time.perf_counter() - t0
        log.info("FIDO   done in %.2fs, %d file(s)", dt_b, len(paths_b))
        path_b = _pick_closest(paths_b, target_time)
        info_b, map_b, hdr_b = inspect_map(path_b, "fido")

    # ---- Inspect ----
    info_a, map_a, hdr_a = inspect_map(path_a, "as-is")

    print("\n" + "=" * 70)
    print("HEADER + MAP SUMMARY")
    print("=" * 70)
    for info in (info_a, info_b):
        if info is None:
            continue
        print(f"\n[{info['label']}]")
        for k, v in info.items():
            print(f"  {k:>15s}: {v}")

    # Header key diff vs Fido (if available)
    allclose: object = "skipped (no Fido baseline)"
    if hdr_b is not None:
        keys_a = set(hdr_a.keys())
        keys_b = set(hdr_b.keys())
        only_a = sorted(keys_a - keys_b)
        only_b = sorted(keys_b - keys_a)
        print(f"\nheader keys: as-is={len(keys_a)}, fido={len(keys_b)}")
        if only_a:
            print(f"  keys ONLY in as-is ({len(only_a)}): {only_a}")
        if only_b:
            print(f"  keys ONLY in fido  ({len(only_b)}): {only_b}")

        same_record = info_b["label"] == "fido"   # only true when we downloaded same time
        if same_record:
            import numpy as np
            if map_a.data.shape == map_b.data.shape:
                try:
                    allclose = bool(np.allclose(map_a.data.astype("float64"),
                                                map_b.data.astype("float64"), equal_nan=True))
                except Exception as e:
                    allclose = f"compare failed: {e}"
            else:
                allclose = f"shape mismatch: {map_a.data.shape} vs {map_b.data.shape}"
            print(f"\ndata arrays allclose: {allclose}")
        else:
            print("\ndata arrays: skipped (Fido baseline is a different record)")
    else:
        print(f"\nheader keys: as-is={len(hdr_a)} (no Fido baseline to diff)")

    # Register test (the first AIA prep step)
    ok_a, msg_a = try_register(map_a, "as-is")
    print(f"\nregister() check:\n  {msg_a}")
    if map_b is not None:
        _ok_b, msg_b = try_register(map_b, "fido")
        print(f"  {msg_b}")

    # Pass criteria
    print("\n" + "=" * 70)
    print("PASS CRITERIA")
    print("=" * 70)
    crit = {
        "register_as_is":  ok_a,
        "has_exptime":     info_a["exptime"] is not None,
        "has_pointing":    info_a["x0_mp"] is not None and info_a["y0_mp"] is not None,
        "has_t_obs":       info_a["t_obs"] is not None,
    }
    if hdr_b is not None and info_b["label"] == "fido":
        crit["data_allclose"] = allclose is True
    for k, v in crit.items():
        print(f"  {'PASS' if v else 'FAIL'}  {k}")

    if dt_b is not None:
        print(f"\nTIMING:  as-is={dt_a:.2f}s   fido={dt_b:.2f}s   "
              f"speedup={dt_b / dt_a:.1f}×")
    else:
        print(f"\nTIMING:  as-is={dt_a:.2f}s   (no Fido baseline)")

    if not args.keep_tmp:
        shutil.rmtree(tmp_root, ignore_errors=True)
        log.info("cleaned up %s", tmp_root)
    else:
        log.info("kept tmp_root=%s", tmp_root)

    return 0 if all(crit.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
