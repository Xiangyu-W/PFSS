"""Aggregate per-event footpoint_temperature_summary_latest.csv into one CSV.

Walks `results/SOLO_*/extract/footpoint_temperature_summary_latest.csv`,
concatenates them (sorted by spacecraft time), tags each row with the
source path and the aggregation timestamp, and writes the combined CSV
to a timestamped filename so future aggregations don't collide.

Usage:
    python scripts/aggregate_summary.py
    python scripts/aggregate_summary.py --results-root /path/to/results
    python scripts/aggregate_summary.py --out my_aggregate.csv
"""
from __future__ import annotations

import argparse
import logging
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
log = logging.getLogger("aggregate_summary")


def _compact(iso_ts: str) -> str:
    """`2022-03-02T06:00:00` -> `20220302T060000`."""
    return pd.Timestamp(iso_ts).strftime("%Y%m%dT%H%M%S")


def aggregate(results_root: Path, out_path: Path | None) -> Path:
    csvs = sorted(results_root.glob("*/extract/footpoint_temperature_summary_latest.csv"))
    if not csvs:
        raise SystemExit(f"no *_latest.csv found under {results_root}")

    dfs = []
    for p in csvs:
        df = pd.read_csv(p)
        df["source_csv"] = str(p.relative_to(results_root))
        dfs.append(df)
    out = pd.concat(dfs, ignore_index=True)
    out["aggregated_at"] = datetime.now().strftime("%Y%m%dT%H%M%S")
    if "sc_time_utc" in out.columns:
        out = out.sort_values("sc_time_utc", kind="stable").reset_index(drop=True)

    if out_path is None:
        start = _compact(out["sc_time_utc"].iloc[0])
        end = _compact(out["sc_time_utc"].iloc[-1])
        out_path = results_root / f"footpoint_temperature_aggregated_{start}_to_{end}.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    log.info("aggregated %d events from %d CSVs -> %s", len(out), len(csvs), out_path)
    return out_path


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--results-root", type=Path, default=REPO_ROOT / "results",
                   help="results root (default: %(default)s)")
    p.add_argument("--out", type=Path, default=None,
                   help="output CSV path (default: <results-root>/"
                        "footpoint_temperature_aggregated_<start>_to_<end>.csv, "
                        "where start/end are the min/max sc_time_utc found)")
    args = p.parse_args(argv)
    logging.basicConfig(level="INFO",
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    aggregate(args.results_root, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
