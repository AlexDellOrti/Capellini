"""NCBI mapping stage: download taxonomy names and assign real NCBI taxids."""

from __future__ import annotations

import logging
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

from capellini.config import CapelliniConfig
from capellini.utils.taxonomy import (
    RANKS_FOR_TAXID,
    build_name_to_ncbi,
    assign_ncbi_taxids,
)

logger = logging.getLogger(__name__)

NCBI_TAXDMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"


def download_ncbi_names(download_path: str | Path) -> Path:
    """Download NCBI taxdmp.zip, extract names.dmp, and delete the zip.

    Skips the download if names.dmp already exists.

    Args:
        download_path: Directory where names.dmp will be saved.

    Returns:
        Path to the names.dmp file.
    """
    download_path = Path(download_path)
    download_path.mkdir(parents=True, exist_ok=True)
    names_dmp_path = download_path / "names.dmp"

    if names_dmp_path.exists():
        logger.info("names.dmp already exists — skipping download: %s", names_dmp_path)
        return names_dmp_path

    zip_path = download_path / "taxdmp.zip"
    logger.info("Downloading NCBI taxdmp.zip from %s", NCBI_TAXDMP_URL)
    urllib.request.urlretrieve(NCBI_TAXDMP_URL, str(zip_path))

    logger.debug("Extracting names.dmp from %s", zip_path)
    with zipfile.ZipFile(str(zip_path), "r") as zf:
        with zf.open("names.dmp") as src, open(str(names_dmp_path), "wb") as dst:
            dst.write(src.read())

    zip_path.unlink(missing_ok=True)
    logger.info("Saved: %s", names_dmp_path)
    return names_dmp_path


def run_ncbi_mapping(cfg: CapelliniConfig) -> pd.DataFrame:
    """Load taxonomy table, assign NCBI taxids, and return the updated DataFrame.

    Loads the DADA2-produced taxonomy_table_{F|R|P}.csv, downloads NCBI names if
    needed, looks up real NCBI taxids for each ASV (finest available rank), and adds
    NCBI_taxid and taxid_matched_rank columns.

    Args:
        cfg: Populated CapelliniConfig instance.

    Returns:
        taxonomy_table DataFrame with NCBI_taxid and taxid_matched_rank columns added.
    """
    logger.info("NCBI mapping: loading taxonomy table")

    suffix = {"forward": "F", "reverse": "R", "paired": "P"}.get(cfg.direction, "F")
    taxonomy_path = Path(cfg.dada2_folder) / f"taxonomy_table_{suffix}.csv"
    taxonomy_table = pd.read_csv(taxonomy_path)
    logger.info("Loaded taxonomy table: %s rows", len(taxonomy_table))

    names_dmp_path = download_ncbi_names(cfg.download_path)

    logger.info("Building name -> NCBI taxid mapping from names.dmp")
    name_to_ncbi = build_name_to_ncbi(str(names_dmp_path))
    logger.info("Loaded %s scientific names", len(name_to_ncbi))

    taxonomy_table = assign_ncbi_taxids(taxonomy_table, name_to_ncbi)

    logger.info("NCBI mapping complete")
    return taxonomy_table
