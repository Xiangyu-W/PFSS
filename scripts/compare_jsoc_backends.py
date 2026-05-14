"""Side-by-side comparison of drms `as-is` vs Fido `fits` JSOC backends.

For one event (default `event_2022-03-03T0600`):
  1. Drms `as-is` fetch into <out>/drms_as_is/   (uses pfss_pipeline.aia.fetch.fetch_drms_as_is)
  2. Fido  `fits`  fetch into <out>/fido_fits/   (uses pfss_pipeline.aia.fetch.fetch_jsoc_fido)
  3. Per wavelength: compare filename T_REC, header diff, pixel data
  4. Plot 6×3 grid (as-is | fido | diff) and a per-wavelength stats panel.
  5. Write report.md

Nothing inside <out> touches the real shared `aia_data_dir` cache or the run
manifest. Both backends download fresh files into the comparison folder.

Usage:
    python scripts/compare_jsoc_backends.py \\
        --event configs/event_2022-03-03T0600.yaml \\
        [--out results/<event_id>/jsoc_backend_comparison]
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import yaml

# Force certifi CA bundle (matches __main__.py behavior).
try:
    import certifi
    os.environ["SSL_CERT_FILE"] = certifi.where()
    os.environ["REQUESTS_CA_BUNDLE"] = certifi.where()
except ImportError:
    pass

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))
log = logging.getLogger("compare_jsoc_backends")


def _event_id(cfg: dict) -> str:
    from astropy.time import Time
    sct = Time(cfg["irap"]["spacecraft_time"]).strftime("%Y%m%dT%H%M%S")
    return f"{cfg['irap']['spacecraft']}_{sct}_{cfg['irap']['coronal_model']}"


def _load_target_time(cfg_path: Path):
    from astropy.time import Time
    from pfss_pipeline import config as cfg_mod
    cfg = cfg_mod.load_config(cfg_path)  # merges DEFAULTS so cfg.aia.wavelengths exists
    results_root = Path(cfg["results_root"])
    eid = _event_id(cfg)
    manifest = results_root / eid / "run_manifest.yaml"
    m = yaml.safe_load(manifest.read_text())
    tt = m["stages"]["irap_fetch"]["target_time"]
    return Time(tt), cfg, eid


def _parse_l1_time(p) -> datetime:
    s = str(p)
    m = re.search(r"aia\.lev1_euv_12s\.(\d{4}-\d{2}-\d{2}T\d{6})Z", s)
    if m:
        return datetime.strptime(m.group(1), "%Y-%m-%dT%H%M%S")
    m = re.search(r"aia\.lev1_euv_12s\.(\d{8}T\d{6})Z", s)
    if m:
        return datetime.strptime(m.group(1), "%Y%m%dT%H%M%S")
    raise ValueError(f"can't parse timestamp from {s}")


def fetch_both(target_time, wavelengths, jsoc_notify, out_dir: Path):
    """Run both backends fresh, return (drms_paths, fido_paths, timings)."""
    from pfss_pipeline.aia import fetch as fm

    drms_dir = out_dir / "drms_as_is"
    fido_dir = out_dir / "fido_fits"
    drms_dir.mkdir(parents=True, exist_ok=True)
    fido_dir.mkdir(parents=True, exist_ok=True)

    timings = {}

    log.info("=== drms as-is ===")
    t0 = time.perf_counter()
    drms_paths = fm.fetch_drms_as_is(target_time, wavelengths, drms_dir, jsoc_notify)
    timings["drms_as_is"] = time.perf_counter() - t0
    log.info("drms as-is done in %.1fs", timings["drms_as_is"])

    log.info("=== Fido fits ===")
    t0 = time.perf_counter()
    fido_paths = fm.fetch_jsoc_fido(target_time, wavelengths, fido_dir, jsoc_notify)
    timings["fido_fits"] = time.perf_counter() - t0
    log.info("Fido fits done in %.1fs", timings["fido_fits"])

    return drms_paths, fido_paths, timings


def compare_pair(path_a: str, path_b: str):
    """Compare two AIA L1 FITS files. Returns a dict of comparison stats."""
    from astropy.io import fits
    import numpy as np

    with fits.open(path_a) as ha, fits.open(path_b) as hb:
        # Find image HDU (AIA L1 typically has the compressed image in HDU 1)
        idx_a = next(i for i, h in enumerate(ha) if h.data is not None)
        idx_b = next(i for i, h in enumerate(hb) if h.data is not None)
        hdr_a, hdr_b = ha[idx_a].header, hb[idx_b].header
        data_a, data_b = ha[idx_a].data, hb[idx_b].data

    keys_a = set(hdr_a.keys()) - {""}
    keys_b = set(hdr_b.keys()) - {""}
    only_a = sorted(keys_a - keys_b)
    only_b = sorted(keys_b - keys_a)

    # Science-key value diffs (keys present in both)
    science_keys = ["T_REC", "T_OBS", "DATE-OBS", "WAVELNTH", "EXPTIME",
                    "LVL_NUM", "X0_MP", "Y0_MP", "CROTA2",
                    "CRPIX1", "CRPIX2", "CDELT1", "CDELT2"]
    val_diffs = {}
    for k in science_keys:
        if k in hdr_a and k in hdr_b:
            va, vb = hdr_a[k], hdr_b[k]
            if isinstance(va, float) and isinstance(vb, float):
                same = (np.isnan(va) and np.isnan(vb)) or va == vb
            else:
                same = va == vb
            if not same:
                val_diffs[k] = (va, vb)

    # T_REC times
    t_rec_a = hdr_a.get("T_REC", "")
    t_rec_b = hdr_b.get("T_REC", "")
    same_record = t_rec_a == t_rec_b and t_rec_a != ""

    # Data comparison
    shape_a, shape_b = data_a.shape, data_b.shape
    same_shape = shape_a == shape_b
    pixel_stats = {}
    if same_shape:
        a64 = data_a.astype("float64")
        b64 = data_b.astype("float64")
        d = a64 - b64
        finite = np.isfinite(d)
        pixel_stats = {
            "max_abs_diff": float(np.max(np.abs(d[finite])) if finite.any() else 0),
            "mean_abs_diff": float(np.mean(np.abs(d[finite])) if finite.any() else 0),
            "n_diff_pixels": int(np.sum((d != 0) & finite)),
            "n_total_pixels": int(np.sum(finite)),
            "allclose_rtol_1e-5": bool(np.allclose(a64, b64, rtol=1e-5, equal_nan=True)),
            "byte_identical": bool(np.array_equal(data_a, data_b, equal_nan=True)),
        }

    return {
        "path_a": path_a,
        "path_b": path_b,
        "size_a_bytes": Path(path_a).stat().st_size,
        "size_b_bytes": Path(path_b).stat().st_size,
        "t_rec_a": str(t_rec_a),
        "t_rec_b": str(t_rec_b),
        "same_record": same_record,
        "fname_time_a": _parse_l1_time(path_a).isoformat(),
        "fname_time_b": _parse_l1_time(path_b).isoformat(),
        "fname_time_diff_s": (
            _parse_l1_time(path_a) - _parse_l1_time(path_b)).total_seconds(),
        "n_header_keys_a": len(hdr_a),
        "n_header_keys_b": len(hdr_b),
        "keys_only_in_a": only_a,
        "keys_only_in_b": only_b,
        "science_key_value_diffs": val_diffs,
        "shape_a": shape_a,
        "shape_b": shape_b,
        "same_shape": same_shape,
        **pixel_stats,
    }


def plot_overview(drms_paths, fido_paths, out_png: Path):
    """6 rows (one per wavelength) × 3 cols (as-is | fido | diff)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from sunpy.map import Map

    wls = sorted(drms_paths.keys(), key=int)
    n = len(wls)
    fig, axes = plt.subplots(n, 3, figsize=(13, 3.5 * n),
                             subplot_kw={"projection": None})

    for r, wl in enumerate(wls):
        m_a = Map(drms_paths[wl])
        m_b = Map(fido_paths[wl])
        d_a = np.asarray(m_a.data).astype("float64")
        d_b = np.asarray(m_b.data).astype("float64")

        # match shape if necessary (should be 4096x4096 both)
        if d_a.shape != d_b.shape:
            note = f"shape mismatch {d_a.shape} vs {d_b.shape}"
        else:
            note = ""
        diff = d_a - d_b if d_a.shape == d_b.shape else None

        vmin, vmax = np.nanpercentile(d_a, [1, 99.5])
        axes[r, 0].imshow(d_a, origin="lower", cmap="sdoaia" + wl if False else "inferno",
                          vmin=vmin, vmax=vmax)
        axes[r, 0].set_title(f"{wl} Å as-is  T_REC≈{_parse_l1_time(drms_paths[wl]).strftime('%H:%M:%S')}")
        axes[r, 0].axis("off")
        vmin, vmax = np.nanpercentile(d_b, [1, 99.5])
        axes[r, 1].imshow(d_b, origin="lower", cmap="inferno", vmin=vmin, vmax=vmax)
        axes[r, 1].set_title(f"{wl} Å fido   T_REC≈{_parse_l1_time(fido_paths[wl]).strftime('%H:%M:%S')}")
        axes[r, 1].axis("off")

        if diff is not None:
            v = float(np.nanpercentile(np.abs(diff), 99.5))
            v = max(v, 1.0)
            im2 = axes[r, 2].imshow(diff, origin="lower", cmap="RdBu_r",
                                    vmin=-v, vmax=v)
            axes[r, 2].set_title(f"{wl} Å (as-is − fido)  max|Δ|p99.5={v:.2g}")
            plt.colorbar(im2, ax=axes[r, 2], fraction=0.046, pad=0.04)
        else:
            axes[r, 2].text(0.5, 0.5, note, ha="center", va="center")
        axes[r, 2].axis("off")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_png


def plot_summary(comparisons: dict, out_png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    wls = sorted(comparisons.keys(), key=int)
    fname_diff = [comparisons[w]["fname_time_diff_s"] for w in wls]
    max_abs_diff = [comparisons[w].get("max_abs_diff", 0) for w in wls]
    n_diff = [comparisons[w].get("n_diff_pixels", 0) for w in wls]
    n_tot = [comparisons[w].get("n_total_pixels", 1) for w in wls]
    pct_diff = [100 * n_d / max(n_t, 1) for n_d, n_t in zip(n_diff, n_tot)]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].bar(wls, fname_diff)
    axes[0].set_title("T_REC time gap (as-is − fido, seconds)")
    axes[0].set_xlabel("wavelength (Å)")
    axes[0].axhline(0, color="k", lw=0.5)

    axes[1].bar(wls, max_abs_diff)
    axes[1].set_title("max |pixel diff| (DN)")
    axes[1].set_xlabel("wavelength (Å)")
    axes[1].set_yscale("symlog")

    axes[2].bar(wls, pct_diff)
    axes[2].set_title("% pixels differing (any nonzero)")
    axes[2].set_xlabel("wavelength (Å)")

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return out_png


def write_report(out_dir: Path, comparisons: dict, timings: dict, target_time):
    rep = out_dir / "report.md"
    lines = []
    lines.append("# JSOC backend comparison: drms `as-is` vs Fido `fits`\n")
    lines.append(f"- **event target_time** (UTC): `{target_time.iso}`")
    lines.append(f"- **drms as-is fetch wall-time**: `{timings['drms_as_is']:.1f} s`")
    lines.append(f"- **Fido fits fetch wall-time**: `{timings['fido_fits']:.1f} s`")
    lines.append(f"- **speedup**: `{timings['fido_fits'] / max(timings['drms_as_is'], 0.001):.2f}×`\n")

    lines.append("## Per-wavelength summary\n")
    lines.append("| WL | as-is T_REC | fido T_REC | Δt (s) | same record? | byte-identical | n_diff px | max|Δ| | mean|Δ| | only-in-as-is hdr keys | only-in-fido hdr keys |")
    lines.append("|---:|---|---|---:|:---:|:---:|---:|---:|---:|---|---|")
    for wl in sorted(comparisons.keys(), key=int):
        c = comparisons[wl]
        lines.append("| {wl} | {a} | {b} | {dt:.0f} | {same} | {bi} | {nd}/{nt} | {ma:.3g} | {me:.3g} | {oa} | {ob} |".format(
            wl=wl,
            a=c.get("t_rec_a", ""),
            b=c.get("t_rec_b", ""),
            dt=c.get("fname_time_diff_s", float("nan")),
            same="✓" if c.get("same_record") else "✗",
            bi="✓" if c.get("byte_identical") else "✗",
            nd=c.get("n_diff_pixels", "?"),
            nt=c.get("n_total_pixels", "?"),
            ma=c.get("max_abs_diff", float("nan")),
            me=c.get("mean_abs_diff", float("nan")),
            oa=", ".join(c.get("keys_only_in_a", [])) or "—",
            ob=", ".join(c.get("keys_only_in_b", [])) or "—",
        ))

    lines.append("\n## Header value diffs on science keys\n")
    for wl in sorted(comparisons.keys(), key=int):
        c = comparisons[wl]
        sd = c.get("science_key_value_diffs", {})
        if not sd:
            lines.append(f"- {wl} Å: science-key values **all match** for shared keys.")
        else:
            lines.append(f"- {wl} Å:")
            for k, (va, vb) in sd.items():
                lines.append(f"    - `{k}`: as-is=`{va}` vs fido=`{vb}`")

    lines.append("\n## Plots\n")
    lines.append("- `plots/overview.png` — per-wavelength (as-is | fido | diff)")
    lines.append("- `plots/summary.png` — per-wavelength stats bar charts")

    lines.append("\n## Full per-wavelength JSON")
    lines.append("```")
    lines.append(json.dumps(comparisons, indent=2, default=str))
    lines.append("```")

    rep.write_text("\n".join(lines) + "\n")
    return rep


def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--event", type=Path,
                  default=Path("configs/event_2022-03-03T0600.yaml"))
    p.add_argument("--out", type=Path, default=None,
                  help="Output dir (default: results/<event_id>/jsoc_backend_comparison)")
    args = p.parse_args(argv)
    logging.basicConfig(level="INFO", format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    target_time, cfg, eid = _load_target_time(args.event)
    out_dir = args.out or (REPO_ROOT / "results" / eid / "jsoc_backend_comparison")
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("event=%s target=%s out=%s", eid, target_time.iso, out_dir)

    wavelengths = cfg["aia"]["wavelengths"]
    jsoc_notify = cfg["aia"]["jsoc_notify"]

    drms_paths, fido_paths, timings = fetch_both(
        target_time, wavelengths, jsoc_notify, out_dir)

    # Per-wavelength comparisons
    comparisons = {}
    for wl in wavelengths:
        if wl not in drms_paths:
            log.warning("missing as-is %s", wl); continue
        if wl not in fido_paths:
            log.warning("missing fido %s", wl); continue
        log.info("comparing %s Å", wl)
        comparisons[wl] = compare_pair(drms_paths[wl], fido_paths[wl])

    # Plots
    log.info("rendering plots …")
    plot_overview(drms_paths, fido_paths, out_dir / "plots" / "overview.png")
    plot_summary(comparisons, out_dir / "plots" / "summary.png")

    rep = write_report(out_dir, comparisons, timings, target_time)
    log.info("report -> %s", rep)
    log.info("plots  -> %s", out_dir / "plots")
    return 0


if __name__ == "__main__":
    sys.exit(main())
