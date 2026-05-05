"""DEM diagnostic figures."""
from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap

log = logging.getLogger(__name__)

_RES_BOUNDS = [-3, -2, -1, -0.5, 0.5, 1, 2, 3]
_RES_COLORS = [
    "#08519c", "#3182bd", "#2ca25f",
    "#f0f0f0",
    "#fee391", "#fcae91", "#cb181d",
]


def response_curves(T_resp_logt, T_resp_matrix, channels: list[str], out_path: Path,
                    dpi: int = 150) -> None:
    fig, ax = plt.subplots(1, figsize=(6, 5))
    for i in range(T_resp_matrix.shape[1]):
        ax.plot(T_resp_logt, T_resp_matrix[:, i], label=channels[i])
    ax.legend()
    ax.set_xlabel("log T (K)")
    ax.set_ylabel("Response")
    ax.set_xticks(np.arange(4, 9, 0.5))
    ax.set_xlim(4, 8.5)
    ax.set_ylim(3e-28, 3e-24)
    ax.set_yscale("log")
    ax.set_title("AIA Temperature Response")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("saved %s", out_path)


def reconstruction_comparison(submaps_data, dn_reconstructed, submaps_data_err,
                              wavenum: list[str], out_path: Path, dpi: int = 200) -> None:
    """3-row × n-channel grid: observed / reconstructed / residual (σ)."""
    nf = len(wavenum)
    cmap = ListedColormap(_RES_COLORS)
    norm_res = BoundaryNorm(_RES_BOUNDS, cmap.N)

    fig, axes = plt.subplots(3, nf, figsize=(4 * nf, 11))
    plt.subplots_adjust(wspace=0.02, hspace=0.005)

    im_obs = im_rec = im_res = None
    for j in range(nf):
        vmax = np.percentile(submaps_data[:, :, j], 99)
        n = plt.Normalize(vmin=0, vmax=vmax)
        im_obs = axes[0, j].imshow(submaps_data[:, :, j], origin="lower", cmap="inferno", norm=n)
        axes[0, j].set_title(f"{wavenum[j]} Å observed", fontsize=10)
        im_rec = axes[1, j].imshow(dn_reconstructed[:, :, j], origin="lower", cmap="inferno", norm=n)
        axes[1, j].set_title(f"{wavenum[j]} Å reconstructed", fontsize=10)
        with np.errstate(invalid="ignore"):
            res = (submaps_data[:, :, j] - dn_reconstructed[:, :, j]) / (submaps_data_err[:, :, j] + 1e-30)
        im_res = axes[2, j].imshow(res, origin="lower", cmap=cmap, norm=norm_res)
        axes[2, j].set_title(f"{wavenum[j]} Å residual", fontsize=10)
        for r in range(3):
            axes[r, j].set_xticks([])
            axes[r, j].set_yticks([])

    axes[0, 0].set_ylabel("Observed", fontsize=12)
    axes[1, 0].set_ylabel("Reconstructed", fontsize=12)
    axes[2, 0].set_ylabel("Residual (σ)", fontsize=12)
    fig.colorbar(im_obs, ax=axes[0, :].tolist(), shrink=0.8, label="DN/s")
    fig.colorbar(im_rec, ax=axes[1, :].tolist(), shrink=0.8, label="DN/s")
    fig.colorbar(im_res, ax=axes[2, :].tolist(), shrink=0.8,
                 label="(obs − rec) / σ$_{err}$",
                 ticks=_RES_BOUNDS, spacing="proportional")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("saved %s", out_path)


def temperature_map_pixel(T_mean, out_path: Path, vmax: float = 2.1e6, dpi: int = 200) -> None:
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111)
    cax = ax.imshow(T_mean, cmap="inferno", origin="lower", vmax=vmax)
    fig.colorbar(cax).set_label("Temperature (K)")
    ax.set_title("DEM-weighted mean temperature")
    ax.set_xlabel("X (pixel)")
    ax.set_ylabel("Y (pixel)")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("saved %s", out_path)


def dem_bins_grid(dem, logt_centers, out_path: Path, dpi: int = 150) -> None:
    n = dem.shape[2]
    nrows = int(np.ceil(n / 4))
    fig = plt.figure(figsize=(16, 4 * nrows))
    for j in range(n):
        plt.subplot(nrows, 4, j + 1)
        with np.errstate(invalid="ignore"):
            plt.imshow(np.log10(dem[:, :, j] + 1e-20), "inferno",
                       vmin=17, vmax=25, origin="lower")
        plt.title(f"log T = {logt_centers[j]:.2f}")
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    log.info("saved %s", out_path)
