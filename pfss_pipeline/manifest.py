"""run_manifest.yaml: cross-stage record of derived target_time + output paths."""
from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import yaml


def read(path: str | Path) -> dict:
    p = Path(path)
    if not p.exists():
        return {}
    with open(p) as fh:
        return yaml.safe_load(fh) or {}


def write(path: str | Path, data: dict) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as fh:
        yaml.safe_dump(data, fh, sort_keys=False)


def update_stage(path: str | Path, stage: str, payload: dict) -> dict:
    """Merge payload under manifest[stages][stage], stamp timestamp, persist."""
    m = read(path)
    m.setdefault("stages", {})
    payload = {**payload, "completed_at": datetime.now(timezone.utc).isoformat()}
    m["stages"][stage] = {**m["stages"].get(stage, {}), **payload}
    write(path, m)
    return m


def get_target_time(manifest: dict) -> str | None:
    try:
        return manifest["stages"]["irap_fetch"]["target_time"]
    except KeyError:
        return None
