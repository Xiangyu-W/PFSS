# pfss_pipeline

End-to-end pipeline for Solar Orbiter / Parker Solar Probe footpoint analysis:
fetch IRAP footpoints + ADAPT magnetogram → prep AIA → DEM inversion →
DEM-on-AIA overlay with hull statistics. CLI-driven, YAML-configured,
cache-aware.

Replaces the prototype notebooks under `test_code/`.

---

## Directory layout

```
pfss_pipeline/
├── __main__.py            # `python -m pfss_pipeline` entry
├── cli.py                 # arg parsing + stage dispatch
├── config.py              # DEFAULTS, YAML deep-merge, validation
├── manifest.py            # run_manifest.yaml read/write
├── paths.py               # OutputLayout: per-event directory structure
├── io_utils.py            # find_closest_map, etc.
│
├── stages/                # one runner per pipeline stage
│   ├── irap_fetch.py      # → footpoints + ADAPT
│   ├── aia_prep.py        # → AIA L1 → PSF/register/degradation
│   ├── dem.py             # → T_mean / T_peak / EM (with cache)
│   └── extract.py         # → overlay PNG + hull stats + summary CSV
│
├── irap/                  # IRAP/ADAPT helpers (used by irap_fetch + extract)
│   ├── adapt.py           # ADAPT load + realization select
│   ├── footpoints.py      # IRAP footpoint table + convex hull
│   ├── overlay.py         # Carrington reproject helpers, contour levels
│   └── mct.py             # multi-coronal-temperature helpers
│
├── aia/                   # AIA fetch + prep + diagnostic plots
│   ├── fetch.py
│   ├── prep.py
│   └── plots.py
│
└── dem/                   # DEM inversion building blocks
    ├── tbins.py           # logT range / bin grid resolver (CH/AR/custom)
    ├── core.py            # T-response, submaps+errors, dn2dem call
    ├── derived.py         # T_mean/T_peak/EM, NaN fill, FITS wrapping
    └── plots.py           # response curves, reconstruction, T map, DEM bins
```

The package depends on `demregpy` (dn2dem), `aiapy`, `sunpy`, `astropy`,
`reproject`, `scipy`.

---

## Stage flow

```
┌───────────────┐
│ irap_fetch    │  IRAP footpoints (ECSV) + ADAPT magnetogram (FITS)
│  ─ inputs ─   │  YAML: irap.* (spacecraft, time, model, prob_threshold, …)
│  ─ outputs ─  │  irap/dominant_ft_*.ecsv, ADAPT cached path in manifest
└───────┬───────┘
        ▼
┌───────────────┐
│ aia_prep      │  Per-wavelength AIA L1 → PSF deconv → update_pointing →
│               │  register → correct_degradation
│  ─ inputs ─   │  YAML: aia.wavelengths, aia.do_psf_deconvolve, target_time
│  ─ outputs ─  │  /disk/plasma/xw2/sunpy/data/aia_prep/aia_prep_<wl>A_*.fits
└───────┬───────┘
        ▼
┌───────────────┐
│ dem           │  HPC ROI submaps + estimate_error → dn2dem → T_mean / T_peak
│               │  / EM. NaN-fill on T_mean, raw kept separately.
│  ─ inputs ─   │  YAML: dem.{roi, region_type, logT_range, bin_width,
│               │             clip_negative, fill_nan, match_tolerance_min}
│  ─ outputs ─  │  dem/{T_mean, T_mean_raw, T_peak, EM}_map_*.fits,
│               │  dem/dem_cube_*.npz, dem/figures/*.png
│  ─ cache ─    │  Skipped if manifest signature matches; auto-invalidates
│               │  on YAML change (logs the diff).
└───────┬───────┘
        ▼
┌───────────────┐
│ extract       │  Diff-rotate dominant footpoints to AIA epoch → convex hull
│               │  → AIA full-disk Carrington reproject + submap → DEM
│               │  reproject onto local Carrington ROI WCS at native shape →
│               │  hull pixel mask → DEM-on-AIA overlay PNG + summary CSV.
│  ─ inputs ─   │  YAML: extract.{overlay_aia_wavelength, overlay_carrington
│               │              _roi, overlay_carrington_pad_deg}
│  ─ outputs ─  │  extract/DEM_on_AIA_overlay_*.png,
│               │  extract/<SC>_T_inside_hull_*.npz,
│               │  extract/footpoint_temperature_summary.csv (append per run)
└───────────────┘
```

`run_manifest.yaml` lives at `<results_root>/<event_id>/run_manifest.yaml`
and records every stage's outputs + config signature for reproducibility
and cache invalidation.

---

## How to run

### Setup

1. Pick or copy a config: `configs/example.yaml` → `configs/<my_event>.yaml`
2. Edit `irap.spacecraft`, `irap.spacecraft_time`, `results_root`, etc.
3. Run from the repo root: `cd /disk/plasma/xw2/PFSS`

### Typical flow

```bash
# 1. Fetch IRAP footpoints + ADAPT (manual review of footpoint figure follows)
python -m pfss_pipeline --config configs/event_2022-03-03.yaml --stage irap

# 2. Inspect the generated footpoint figure, decide HPC ROI, fill in
#    `dem.roi.bottom_left_arcsec` / `top_right_arcsec` and `dem.region_type`
#    in the YAML. Then run everything from AIA prep onward:
python -m pfss_pipeline --config configs/event_2022-03-03.yaml --stage all-after-roi
```

### Single-stage commands

| Stage | When to run |
|---|---|
| `--stage irap` | Changed `irap.*` |
| `--stage aia` | Changed `aia.*` (e.g. wavelengths, PSF flag) |
| `--stage dem` | Changed `dem.*` |
| `--stage extract` | Changed `extract.*` |
| `--stage all-after-roi` | Anything below irap; runs aia → dem → extract |

The DEM cache auto-invalidates on YAML changes — no manual `--force` needed.
Use `--force` to ignore caches anyway, or `--force-stage <stage>` to force
just one inside `all-after-roi`:

```bash
python -m pfss_pipeline --config configs/event.yaml --stage all-after-roi --force-stage dem
```

### Useful CLI flags

```
--results-root /some/path     # override results root for this run only
--target-time 2022-03-01T14:00:00   # override IRAP-derived target time
--log-level DEBUG             # more detail in logs
--dry-run                     # print the plan, do nothing
```

### Output locations

```
<results_root>/<event_id>/
├── run_manifest.yaml
├── irap/    dominant_ft_prob*.ecsv, footpoint figures
├── aia/     (kept in /disk/plasma/xw2/sunpy/data/aia_prep, shared)
├── dem/     T_mean*.fits, T_peak*.fits, EM*.fits, dem_cube_*.npz, figures/
└── extract/ DEM_on_AIA_overlay_*.png, *_T_inside_hull_*.npz,
             footpoint_temperature_summary.csv
```

The `event_id` is auto-derived from `irap.spacecraft` + `spacecraft_time` +
`coronal_model`, e.g. `SOLO_20220303T120000_ADAPT`.

---

## Config (YAML)

Minimal example:

```yaml
irap:
  spacecraft: SOLO
  spacecraft_time: '2022-03-03T12:00:00'
  coronal_model: ADAPT
  mode: SUNTIME
  realization_adapt: 0
  prob_threshold_pct: 60
  carrington_roi: {lon: [10, 100], lat: [-60, 30]}

results_root: /disk/plasma/xw2/PFSS/results

aia:
  jsoc_notify: you@example.com    # required for JSOC requests

dem:
  roi:
    bottom_left_arcsec: [-550, -300]
    top_right_arcsec: [200, 450]
  region_type: 'CH'                # 'CH' | 'AR' | 'custom'
  # logT_range: [5.5, 6.6]         # required only if region_type='custom'

extract:
  overlay_carrington_roi:
    lon: [30, 60]                  # null/null → derive from hull bbox + pad
    lat: [-30, 0]
  # overlay_carrington_pad_deg: 15.0    # used when overlay_carrington_roi is null

plots:
  dpi: 300
runtime:
  log_level: INFO
```

Full set of defaults lives in [config.py](config.py) (`DEFAULTS` dict).
User YAML overrides defaults via deep-merge; CLI flags override both.

---

## DEM FITS metadata

The DEM-derived FITS files (`T_mean`, `T_mean_raw`, `T_peak`, `EM`) embed
the inversion settings as header keywords:

```
DEMRTYPE = 'CH'        # region_type
DEMLOGT1 = 5.5         # logT lower
DEMLOGT2 = 6.5         # logT upper
DEMNBINS = 10          # number of T bins
DEMBINW  = 0.1         # bin width in logT
```

So any reader can recover the inversion parameters from the file alone:

```python
from sunpy.map import Map
m = Map('T_mean_map_*.fits')
print(m.meta.get('DEMRTYPE'), m.meta.get('DEMLOGT1'), m.meta.get('DEMLOGT2'))
```

Older files saved before this change will return `None` for these keys.

---

## Cache invalidation (DEM stage)

`stages/dem.py` records a config signature in `run_manifest.yaml`:

```
region_type, logT_range, n_T_bins, bin_width,
roi_bottom_left_arcsec, roi_top_right_arcsec,
clip_negative, fill_nan, wavelengths
```

On a re-run, if the cached signature differs from the current config, the
stage logs the diff and recomputes. Example:

```
WARNING DEM cache stale; settings changed -> recomputing. Diffs (cached -> current):
  region_type: 'CH' -> 'custom'
  logT_range: [5.5, 6.5] -> [5.5, 6.6]
  n_T_bins: 10 -> 11
```

---

## Cross-checking with the prototype notebooks

The prototype lives in `test_code/DEM_note.ipynb` (DEM inversion) and
`test_code/IRAP汇总.ipynb` (overlay). Pipeline behaviour matches the
notebook for the 2022-03-03 SOLO event:

- Identical AIA ROI, error model, `dn2dem` call, T_mean formula
- Same Carrington reprojection strategy (AIA full-disk + submap; T_mean
  on a local CAR WCS at the DEM's native shape)
- Hull stats use the interpolated T_mean (`fill_nan_2d` linear)
- T_mean_hull_avg agrees within numerical noise

To pin a run to the notebook's exact ROI + scale, set
`extract.overlay_carrington_roi: {lon: [30, 60], lat: [-30, 0]}` in the
event YAML.
