"""DADA2 stage: run DADA2_Pipe.R and move the generated FASTA."""

from __future__ import annotations

import glob
import importlib.resources as pkg_resources
import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from capellini.config import CapelliniConfig

logger = logging.getLogger(__name__)


def _get_r_script_path() -> Path:
    """Locate the bundled DADA2_Pipe.R inside the installed package.

    Returns:
        Absolute path to DADA2_Pipe.R.
    """
    ref = pkg_resources.files("capellini.r_scripts").joinpath("DADA2_Pipe.R")
    return Path(str(ref))


def run_dada2(cfg: CapelliniConfig) -> Path:
    """Run the DADA2_Pipe.R script and move the output FASTA to the input folder.

    The R script is called with six arguments matching the notebook:
        Rscript DADA2_Pipe.R <bacterial_raw_fasta_folder> <dada2_folder>
            <silva_ref_path> <silva_taxmap_path> <direction> <TRUE|FALSE>

    After the run the produced ASV_sequences_{F|R|P}.fasta is moved to
    ``cfg.input_fasta_folder / cfg.bacteria_fasta_name``.

    Args:
        cfg: Populated CapelliniConfig instance.

    Returns:
        Path to the moved bacteria FASTA file.

    Raises:
        FileNotFoundError: If the FASTA is not found after the R run.
        subprocess.CalledProcessError: If Rscript exits non-zero.
    """
    logger.info("DADA2: starting R pipeline")

    r_script = _get_r_script_path()
    fasta_generation_str = "TRUE" if cfg.fasta_generation else "FALSE"

    cmd = [
        "Rscript",
        str(r_script),
        cfg.bacterial_raw_fasta_folder,
        cfg.dada2_folder,
        cfg.silva_ref_path,
        cfg.silva_taxmap_path,
        cfg.direction,
        fasta_generation_str,
    ]

    logger.info("Running DADA2 Pipeline...")
    logger.debug("Command: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=True, text=True)
    logger.debug("STDOUT:\n%s", result.stdout)
    if result.returncode != 0:
        raise subprocess.CalledProcessError(
            result.returncode,
            cmd,
            output=result.stdout,
            stderr=result.stderr,
        )

    # Move generated FASTA to input folder
    bact_path = Path(cfg.input_fasta_folder) / cfg.bacteria_fasta_name
    fasta_files = glob.glob(os.path.join(cfg.dada2_folder, "*.fasta"))

    if fasta_files:
        src = fasta_files[0]
        shutil.move(src, str(bact_path))
        logger.debug("Moved %s -> %s", src, bact_path)

    if not bact_path.is_file():
        raise FileNotFoundError(
            f"Bacteria 16S FASTA not found after DADA2 run: {bact_path}"
        )

    logger.info("DADA2 complete: %s", bact_path.name)
    return bact_path
