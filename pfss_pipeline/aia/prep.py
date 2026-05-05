"""Per-wavelength AIA calibration: PSF deconv → update_pointing → register → correct_degradation."""
from __future__ import annotations

import logging
import warnings
from pathlib import Path

import aiapy.psf
from aiapy.calibrate import correct_degradation, register, update_pointing
from astropy.utils.exceptions import AstropyUserWarning
from sunpy.map import Map

log = logging.getLogger(__name__)


def prepare_one(local_fits: str, pointing_table, correction_table,
                do_psf: bool, save_path: Path) -> Path:
    """Prepare a single AIA Level 1 file, save to `save_path`. Returns save_path."""
    aia_map = Map(local_fits)
    log.info("loaded %s (%s Å)", Path(local_fits).name, aia_map.wavelength.value)

    if do_psf:
        log.info("calculating PSF for %s Å", aia_map.wavelength.value)
        psf = aiapy.psf.calculate_psf(aia_map.wavelength)
        log.info("deconvolving %s Å", aia_map.wavelength.value)
        aia_map = aiapy.psf.deconvolve(aia_map, psf=psf)

    aia_map = update_pointing(aia_map, pointing_table=pointing_table)
    aia_map = register(aia_map)
    aia_map = correct_degradation(aia_map, correction_table=correction_table)

    aia_map.meta["history"] = (aia_map.meta.get("history", "")
                              + f" pfss_pipeline.aia.prep do_psf={do_psf};")

    save_path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyUserWarning)
        aia_map.save(str(save_path), overwrite=True)
    log.info("saved %s", save_path)
    return save_path
