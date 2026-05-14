# PFSS Pipeline — CLI Manual

Two entry points:

| Tool | Scope | When to use |
|---|---|---|
| `python -m pfss_pipeline` | one event, one stage | drive a single event end-to-end, or rerun one stage |
| `python scripts/batch_run.py` | all events in `configs/`, one stage | drive a date-range sweep in parallel |

The batch driver just shells out to `python -m pfss_pipeline` per event, so anything you can do per-event you can do in batch.

---

## 1. Pipeline overview

Stages run in this order; each writes to `manifest.stages.<stage>` so downstream stages can pick up its outputs.

```
irap  →  aia-fetch  →  aia-prep  →  dem  →  extract
            └────── aia ──────┘
                              └────── all-after-roi ──────┘
```

| Stage | Reads | Writes | Network? |
|---|---|---|---|
| `irap` | IRAP MCT service, GONG ADAPT | IRAP zip, dominant-footpoint ECSV, ADAPT map, sets `target_time` | yes |
| `aia-fetch` | local `aia_data_dir`, JSOC if missing | L1 FITS paths in manifest | yes (only if local miss) |
| `aia-prep` | L1 paths from manifest | prepped FITS in `aia_prep_dir` | no |
| `aia` | — | runs `aia-fetch` then `aia-prep` | yes |
| `dem` | prepped AIA, DEM config | T_mean / T_peak / EM FITS | no |
| `extract` | DEM, footpoints, ADAPT, AIA | ROI FITS, hull NPZ, overlay PNG, summary CSV (history + latest) | no |
| `all-after-roi` | — | runs `aia-fetch → aia-prep → dem → extract` | yes |

**Two manual gates** between stages (the YAML can't be filled until you've looked at the output):

1. After `irap`, look at `irap/figures/*.png` and fill `extract.overlay_carrington_roi` (lon/lat box for the Carrington overlay).
2. After `aia-prep`, look at the full-disk AIA preview PNG and fill `dem.roi` (HPC arcsec bounding box) + `dem.region_type` (`CH` / `AR` / `custom`; `custom` also needs `dem.logT_range`).

---

## 2. `python -m pfss_pipeline` — per-event CLI

### Synopsis

```bash
python -m pfss_pipeline \
    --config <event.yaml> \
    --stage  <irap|aia-fetch|aia-prep|aia|dem|extract|all-after-roi> \
    [--force] [--force-stage <stage>] \
    [--target-time <ISO>] [--results-root <dir>] \
    [--log-level DEBUG|INFO|WARNING|ERROR] [--dry-run]
```

### Flags

| Flag | Required | Meaning |
|---|---|---|
| `--config` | yes | Path to event YAML. |
| `--stage` | yes | Which stage to run (see table above). |
| `--force` | no | Bypass skip-if-exists for the chosen stage. |
| `--force-stage <s>` | no | Only meaningful with `--stage all-after-roi`: forces just `<s>` to rerun while others still skip. |
| `--target-time <ISO>` | no | Override the `target_time` written by `irap`. Use to point AIA/DEM at a different timestamp (e.g. for an SC vs SS time-shift study). |
| `--results-root <dir>` | no | Overrides `results_root` from the YAML. |
| `--log-level` | no | Default comes from `runtime.log_level` (default INFO). |
| `--dry-run` | no | Print the resolved paths for this config + manifest state, then exit. No side effects. |

### Stages — details

#### `--stage irap`
- Calls IRAP MCT (Selenium → static ZIP) and downloads the ADAPT GONG map closest to `irap.spacecraft_time`.
- Writes `dominant_ft_prob<N>.ecsv` and sets `target_time = date_surf` (minute-rounded) in `run_manifest.yaml`.
- Skip-if-exists: if zip + ECSV already exist and manifest has `target_time`, the stage exits early. `--force` re-runs everything.
- Requires: `irap.spacecraft`, `irap.spacecraft_time`, `results_root` in YAML.

#### `--stage aia-fetch`
- Resolves AIA L1 paths for every wavelength in `aia.wavelengths`:
  1. **Local first**: glob `aia.lev1_euv_12s.<YYYY-MM-DDTHHMM>*.<wl>.image_lev1.fits` in `aia_data_dir`, pick the second closest to `target_time`.
  2. **JSOC fallback** for any wavelength not found locally. Needs `aia.jsoc_notify` (your email registered with JSOC).
- Writes L1 paths + `backend` field to `manifest.stages.aia_fetch`.
- Skip-if-exists: each wavelength independently. `--force` empties the local short-circuit and refetches everything from JSOC.
- **Two JSOC backends** (selected by `aia.fetch_backend` in config, default `drms-as-is`):
  - **`drms-as-is`** (default, fast): uses `drms.Client.export(protocol="as-is", method="url_quick")`. Skips JSOC's FITS-rewrite export queue; submits one combined query for all wavelengths. ~2–10× faster than `fido`. Filenames come back dashless (`20220301T071647Z`) and are auto-renamed to the pipeline's dashed convention.
  - **`fido`** (legacy): `sunpy.net.Fido` + `protocol="fits"`. Goes through the JSOC export queue (10–30 s staging per request).
  - On as-is failure (any exception), the stage **automatically falls back** to `fido`. Both backends share the same JSOC quota of "1 pending export request per user" — if you hit `[status=7] pending export requests`, wait for the prior request to complete (`drms.Client().export_from_id(<id>).wait()`) before retrying.
- **Concurrency**: still serial by default in batch_run.py (`--workers 1`) because JSOC enforces the 1-pending-per-user limit globally. as-is helps speed, not concurrency.

#### `--stage aia-prep`
- Reads L1 paths from `manifest.stages.aia_fetch` (run `aia-fetch` first or it errors).
- PSF deconvolve (`aia.do_psf_deconvolve`) → register → degradation-correct → save to `aia_prep_dir` (a shared cache, **not** under the event dir).
- Skip-if-exists: per wavelength, looks for an existing prepped FITS within `dem.match_tolerance_min` of `target_time`.
- Drift guard: after skipping, `_verify_prep_drift` re-checks that the filename time of every kept file is within tolerance of `target_time`. Catches stale manifest entries from older runs. If it fires, rerun with `--force`.

#### `--stage aia`
- Shorthand: `aia-fetch` then `aia-prep` in the same process.

#### `--stage dem`
- Requires `dem.roi.bottom_left_arcsec`, `dem.roi.top_right_arcsec`, and `dem.region_type` set in YAML.
- Reads prepped AIA, runs DEMReg, writes `T_mean`, `T_peak`, `EM` FITS to `dem/`.
- Skip-if-exists: based on the output FITS existing for the current `target_time`. `--force` rebuilds.

#### `--stage extract`
- Final stage. Combines DEM + dominant footpoints + ADAPT contours + AIA Carrington reprojection.
- Writes:
  - `T_{mean,peak}_in_roi_*.fits`, `EM_in_roi_*.fits` — DEM clipped to HPC ROI
  - `<SC>_T_inside_hull_*.npz` — temperature pixels inside the diff-rotated dominant-footpoint hull (arrays `values`, `lon`, `lat`)
  - `DEM_on_AIA_overlay_*.png` — diagnostic overlay
  - `footpoint_temperature_summary.csv` — **append-only history**, one row per rerun (audit log)
  - `footpoint_temperature_summary_latest.csv` — **upsert by key** (`event_id, target_time, spacecraft, model, mode, region_type, logT_range, roi, prob_threshold`); use this for cross-event analysis
- Schema drift: if either CSV's header doesn't match `SUMMARY_COLUMNS`, it's auto-archived as `*.old-<UTC>.csv` and a fresh file starts.

#### `--stage all-after-roi`
- Convenience: `aia-fetch → aia-prep → dem → extract`. Each sub-stage still respects its own skip-if-exists.
- `--force-stage dem` (for example) re-runs only DEM; the others still skip if their outputs exist.

### Skip-if-exists summary

| Stage | What "exists" means |
|---|---|
| `irap` | zip + ECSV present and manifest has `target_time` |
| `aia-fetch` | local L1 FITS for the target minute found |
| `aia-prep` | prepped FITS within `dem.match_tolerance_min` of `target_time` |
| `dem` | T_mean / T_peak / EM FITS for current `target_time` |
| `extract` | always runs (cheap, writes outputs unconditionally) |

`--force` bypasses skip for the chosen stage. `--force-stage` bypasses for one sub-stage of `all-after-roi`.

### Examples

Single event, full pipeline (with manual gates):

```bash
# 1. footpoints + ADAPT
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml --stage irap

# (open results/<event>/irap/figures/*.png, fill extract.overlay_carrington_roi)

# 2. AIA L1 + prep
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml --stage aia

# (open results/<event>/aia/*.png, fill dem.roi + dem.region_type)

# 3. DEM + extract
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml --stage all-after-roi
```

Other handy invocations:

```bash
# Inspect resolved paths without running anything:
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml --stage extract --dry-run

# Re-run only DEM in the final sweep, but skip aia-* and re-run extract too:
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml --stage all-after-roi --force-stage dem

# Force re-fetch AIA L1 from JSOC (ignore local cache):
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml --stage aia-fetch --force

# Run AIA/DEM/extract for a different target_time than IRAP picked:
python -m pfss_pipeline --config configs/event_2022-03-03T0600.yaml \
    --stage all-after-roi --target-time 2022-03-03T05:30:00
```

---

## 3. `python scripts/batch_run.py` — all events

Drives the per-event CLI over a date range fixed at the top of `batch_run.py`:

```python
START = "2022-03-02 06:00"
END   = "2022-03-07 00:00"
FREQ  = "6h"
```

Edit those + rerun `gen-configs` to extend the sweep.

### Subcommands

```bash
python scripts/batch_run.py gen-configs     [--force]
python scripts/batch_run.py irap            [--workers N] [--from i] [--to j] [--indices 1,3,5] [--force]
python scripts/batch_run.py aia-fetch       [--workers N] [--from i] [--to j] [--indices ...] [--force]
python scripts/batch_run.py aia-prep        [--workers N] [--preview-wl 193] [--from i] [--to j] [--indices ...] [--force]
python scripts/batch_run.py aia             [--workers N] [--preview-wl 193] [--from i] [--to j] [--indices ...] [--force]
python scripts/batch_run.py all-after-roi   [--workers N] [--from i] [--to j] [--indices ...] [--force]
```

| Subcommand | Action |
|---|---|
| `gen-configs` | Copy `configs/event_2022-03-03.yaml` as template; emit one YAML per `event_times()` slot with `irap.spacecraft_time` filled and the manual-gate fields blanked. `--force` overwrites existing configs. |
| `irap` | Run `--stage irap` for selected events. |
| `aia-fetch` | Run `--stage aia-fetch`. **Always use `--workers 1`** — JSOC's export-hash collides under concurrency. |
| `aia-prep` | Run `--stage aia-prep`, then save a full-disk preview PNG (`aia_<wl>A_fulldisk_for_roi.png`) for each successful event. |
| `aia` | `aia-fetch + aia-prep` in one go (worker = 1 to be safe), then preview PNG. |
| `all-after-roi` | Run `--stage all-after-roi`. |

### Common flags

- `--workers N` — concurrent events (ThreadPoolExecutor; each spawns a `python -m pfss_pipeline` subprocess).
- `--from i --to j` — 1-based inclusive event index slice. Indices match the order of `event_times()`.
- `--indices "2,6,10"` — explicit set; overrides `--from/--to`. Comma or space separated, order preserved.
- `--force` — passed through to the per-event CLI.

### Recommended worker caps

| Subcommand | Cap | Reason |
|---|---|---|
| `irap` | ≤6 | Selenium/Chrome + IRAP server rate limits (can't be too much)|
| `aia-fetch` | **1** | JSOC export-hash collision under concurrency |
| `aia-prep` | ~12 (variable) | Pure CPU; BLAS pinned to 1 thread per worker so workers don't oversubscribe |
| `all-after-roi` | ~12 (variable) | DEM is CPU-bound; same BLAS pinning |

### Logs

Per-event subprocess output is captured to `<results_root>/_batch_logs/<event>_<stage>_<UTC>.log`. The aggregate run only prints `OK` / `FAIL` lines; check the log file for the details of a failure.

### Recipe — full sweep

```bash
# 1. Generate per-event YAMLs
python scripts/batch_run.py gen-configs

# 2. IRAP + ADAPT
python scripts/batch_run.py irap --workers 4

# 3. (manual) fill extract.overlay_carrington_roi in each YAML.
#    irap output figures are at results/<event>/irap/figures/*.png

# 4. AIA L1 download (serial)
python scripts/batch_run.py aia-fetch --workers 1

# 5. AIA prep (parallel) + preview PNGs for ROI picking
python scripts/batch_run.py aia-prep --workers 8

# 6. (manual) fill dem.roi + dem.region_type in each YAML.
#    preview PNGs are at results/<event>/aia/aia_193A_fulldisk_for_roi.png

# 7. DEM + extract (parallel)
python scripts/batch_run.py all-after-roi --workers 8
```

Run only a subset (e.g. just the four "noon" events):

```bash
python scripts/batch_run.py all-after-roi --workers 4 --indices 3,7,11,15
```

---

## 4. Event YAML schema

Generated by `gen-configs` from the `event_2022-03-03.yaml` template. Required keys per stage:

```yaml
# Identity (required for all stages)
results_root: /disk/plasma/xw2/PFSS/results
irap:
  spacecraft: SolO                       # "SolO" | "PSP" | ...
  spacecraft_time: "2022-03-03T06:00:00" # ISO
  coronal_model: ADAPT
  mode: SUNTIME
  realization_adapt: 0
  prob_threshold_pct: 60
  carrington_roi: {lon: [10, 100], lat: [-60, 30]}
  # If IRAP returns no 'M'-type rows, set this explicitly (e.g. 'SSW').
  # If unset and 'M' is absent, irap stage raises.
  sw_type: null

# AIA — defaults usually fine
aia:
  wavelengths: ["94", "131", "171", "193", "211", "335"]
  jsoc_notify: "you@example.com"     # must be registered with JSOC
  fetch_backend: "drms-as-is"        # "drms-as-is" (fast, default) | "fido" (legacy)
  do_psf_deconvolve: true
  pointing_window_hours: 12

# DEM — region_type + logT_range gated by region_type
dem:
  roi:
    bottom_left_arcsec: [-300, -300]   # FILL after looking at AIA preview
    top_right_arcsec:  [ 300,  300]
  region_type: CH                       # "CH" | "AR" | "custom"
  logT_range: null                      # required only when region_type=custom, e.g. [5.8, 6.4]
  bin_width: 0.1
  match_tolerance_min: 5

# Extract
extract:
  overlay_aia_wavelength: "193"
  overlay_carrington_roi:               # FILL after looking at IRAP footpoint figure
    lon: [60, 130]
    lat: [-30, 30]
  overlay_carrington_pad_deg: 15.0
```

The pipeline validates in two passes:

- **At load time** (`load_config`): requires `irap.spacecraft`, `irap.spacecraft_time`, `results_root`.
- **At DEM/extract entry** (`assert_dem_ready`): requires `dem.roi.*`, `dem.region_type`. If `region_type=custom`, also `dem.logT_range`.

So you can run `irap` and `aia-*` before filling the DEM gates — they'll only error when you try `dem`/`extract`.

---

## 5. Output layout

```
<results_root>/
  <event_id>/                              # e.g. SolO_20220303T060000_ADAPT
    run_manifest.yaml                      # stage outputs + target_time
    irap/
      mct_raw/<SC>_PARKER_PFSS_<MODE>_<MODEL>_<sct>.zip
      figures/                             # footpoint / hull diagnostics
      dominant_ft_prob<N>.ecsv
    dem/
      T_mean_map_<YYYYMMDD_HHMMSS>.fits
      T_peak_map_*.fits
      EM_map_*.fits
      dem_cube_*.npz
      figures/
    extract/
      T_mean_in_roi_*.fits, T_peak_in_roi_*.fits, EM_in_roi_*.fits
      <SC>_T_inside_hull_*.npz            # values + lon + lat
      DEM_on_AIA_overlay_*.png
      footpoint_temperature_summary.csv          # append-only history
      footpoint_temperature_summary_latest.csv   # upsert-by-key
    aia/                                   # preview PNG written by batch driver
      aia_<wl>A_fulldisk_for_roi.png

  _batch_logs/                             # batch_run.py subprocess logs
    <event>_<stage>_<UTC>.log

<aia_prep_dir>/                            # shared across events (NOT under results_root)
  aia_prep_<wl>A_<YYYYMMDD_HHMMSS>.fits
  images/aia_prep_comparison_<wl>A_*.png
```

`run_manifest.yaml` is the source of truth between stages. Anything in `manifest.stages.<stage>` is what that stage wrote; downstream stages read it instead of re-globbing.

---

## 6. Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `SSL: CERTIFICATE_VERIFY_FAILED` on JSOC | Conda `base` activate hook sets `SSL_CERT_FILE=/usr/share/ssl/certs/ca-bundle.crt` (2019 bundle, missing ISRG Root X1) | `__main__.py` forces certifi already; batch_run.py also injects it into the subprocess env. If you invoke via something else, export `SSL_CERT_FILE=$(python -c 'import certifi; print(certifi.where())')`. |
| `DrmsExportError: ... 1 pending export requests ... [status=7]` | JSOC enforces a "1 pending export per user" quota; a prior request (often a queued `fido` fits export from minutes/hours ago) is still active | Wait it out, or explicitly: `conda run -n DEM python -c "import drms; drms.Client(email='YOU').export_from_id('JSOC_<id>').wait()"` with the ID from the error message. The `drms-as-is` backend triggers this less often than `fido` because as-is requests complete in seconds, but it does share the same quota. |
| `selenium-manager not found` | `selenium-manager` activate hook didn't fire for non-interactive subprocesses | batch_run.py re-derives `SE_MANAGER_PATH` from `sys.executable`. If running pfss_pipeline directly, `conda activate DEM` first or set the env var. |
| `no 'M' type in IRAP solarsurf footpoints` | IRAP returned only SSW/etc.; older code silently fell back to `SSW` | Set `irap.sw_type: SSW` (or whichever type exists) in the event YAML and rerun `irap`. |
| `aia_prep file time drift > Ns from target_time` | manifest has an aia_prep entry from a previous `target_time` (e.g. before you changed it) | rerun `--stage aia-prep --force` |
| `L1 paths missing in manifest for [...]` | `aia-prep` ran without `aia-fetch` first | `--stage aia-fetch` (or `aia`) first |
| `dem.roi.bottom_left_arcsec/top_right_arcsec must be set` | Trying `dem`/`extract` before manual ROI gate | Fill `dem.roi` + `dem.region_type` in the YAML |
| GONG listing fails | Network blip | `find_closest_adapt` falls back to the closest cached ADAPT file; if cache empty it raises. Re-run later or pre-populate `irap.adapt_cache_dir`. |
| `extract` summary CSV has duplicate rows for the same config | This is by design — the `*.csv` file is append-only history | Use `footpoint_temperature_summary_latest.csv` (upsert-by-key) for analysis |
