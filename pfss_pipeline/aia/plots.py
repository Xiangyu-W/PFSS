"""Diagnostic comparison plot: Level 1 vs prepared."""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy import units as u
from sunpy.map import Map

log = logging.getLogger(__name__)


def save_comparison(level1_fits: str, prepared_fits: Path, out_path: Path,
                    dpi: int = 200) -> None:
    """1×2 panel: original Level 1 (left) vs prepared (right)."""
    aia_l1 = Map(level1_fits)
    aia_prep = Map(str(prepared_fits))

    fig = plt.figure(figsize=(16, 8))
    ax1 = fig.add_subplot(1, 2, 1, projection=aia_l1)
    aia_l1.plot(axes=ax1, clip_interval=(1, 99.99) * u.percent)
    ax1.set_title(f"{aia_l1.wavelength.value:.0f} Å — Level 1 (original)")

    ax2 = fig.add_subplot(1, 2, 2, projection=aia_prep)
    aia_prep.plot(axes=ax2, clip_interval=(1, 99.99) * u.percent)
    ax2.set_title(f"{aia_prep.wavelength.value:.0f} Å — Prepared (DN)")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("saved comparison %s", out_path)
