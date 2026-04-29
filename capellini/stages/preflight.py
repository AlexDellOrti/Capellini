"""Pre-flight stage: folder initialization and optional fresh-start cleanup."""

from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

from capellini.config import CapelliniConfig

logger = logging.getLogger(__name__)

# Bundled reference paths (never deleted by fresh_start)
_BUNDLED_16S = Path(__file__).parent.parent / "data" / "references" / "progenome16S.fasta"
_BUNDLED_SPACERS = (
    Path(__file__).parent.parent / "data" / "references" / "spacers" / "spacers_CompleteCollection.fasta"
)


def run_preflight(cfg: CapelliniConfig) -> None:
    """Create all output folders; if fresh_start=True delete previous intermediates.

    Protected files (never deleted):
        - virus FASTA in input_fasta_folder
        - bundled progenome16S.fasta (if present)
        - metadata file

    Args:
        cfg: Populated CapelliniConfig instance.
    """
    logger.info("Pre-flight: initializing folder structure")

    folders_to_manage = [
        cfg.dada2_folder,
        cfg.mmseq_folder,
        cfg.sp_folder,
        cfg.procs_folder,
    ]

    if cfg.fresh_start:
        logger.info("Fresh start: removing previous intermediates (protected files preserved)")
        for folder in folders_to_manage:
            if folder:
                shutil.rmtree(folder, ignore_errors=True)
                os.makedirs(folder, exist_ok=True)

        # Clean input fasta folder but protect critical files
        protected_names: set[str] = set()
        if cfg.virus_fasta_name:
            protected_names.add(cfg.virus_fasta_name)
        if cfg.isolate_ref_16S:
            protected_names.add("progenome16S.fasta")
        if cfg.metadata_path:
            protected_names.add(Path(cfg.metadata_path).name)

        input_folder = Path(cfg.input_fasta_folder)
        if input_folder.exists():
            for fp in input_folder.iterdir():
                if fp.is_file() and fp.name not in protected_names:
                    fp.unlink()
                    logger.debug("Removed: %s", fp)
    else:
        for folder in folders_to_manage:
            if folder:
                os.makedirs(folder, exist_ok=True)

    # Always ensure SpacePHARER subdirectories exist
    if cfg.sp_folder:
        for sub in ("spacers", "databases", "output", "tmp"):
            Path(cfg.sp_folder, sub).mkdir(parents=True, exist_ok=True)

    # Ensure protein / clustering paths exist
    if cfg.proteins_extraction_path:
        os.makedirs(cfg.proteins_extraction_path, exist_ok=True)
    if cfg.clustering_path:
        os.makedirs(cfg.clustering_path, exist_ok=True)
    if cfg.enhanced_networks_folder:
        os.makedirs(cfg.enhanced_networks_folder, exist_ok=True)
    if cfg.input_fasta_folder:
        os.makedirs(cfg.input_fasta_folder, exist_ok=True)

    logger.info("Pre-flight complete")
