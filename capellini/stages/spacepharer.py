"""SpacePHARER stage: spacer extraction, DB creation, prediction, and statistics."""

from __future__ import annotations

import bz2
import logging
import os
import re
import shutil
import subprocess
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO

from capellini.config import CapelliniConfig
from capellini.utils.io import sh

logger = logging.getLogger(__name__)

# Bundled spacers collection path (Modification 2)
_BUNDLED_SPACERS = (
    Path(__file__).parent.parent / "data" / "references" / "spacers" / "spacers_CompleteCollection.fasta"
)


def check_and_install_spacepharer() -> tuple[str, str]:
    """Check that spacepharer and minced are on PATH; install via conda if missing.

    Returns:
        Tuple (spacepharer_path, minced_path).
    """
    if not (shutil.which("spacepharer") and shutil.which("minced")):
        logger.info("Installing spacepharer and minced via conda ...")
        subprocess.run(
            "conda install -y -c conda-forge -c bioconda spacepharer minced",
            shell=True,
            check=True,
        )
        logger.info("SpacePHARER and MinCED installed")
    else:
        logger.info("SpacePHARER and MinCED already installed")

    spacepharer_path = shutil.which("spacepharer") or "spacepharer"
    minced_path = shutil.which("minced") or "minced"
    return spacepharer_path, minced_path


class SpacePHARERWorkflow:
    """Wrapper around the SpacePHARER + MinCED command-line tools."""

    def __init__(self, workdir: str, spacerdir: str) -> None:
        """Initialize directories and tool paths.

        Args:
            workdir: Root SpacePHARER working directory.
            spacerdir: Directory for spacer FASTA files.
        """
        self.wd = Path(workdir)
        self.sd = Path(spacerdir)
        self.spacepharer_path = shutil.which("spacepharer") or "spacepharer"
        self.minced_path = shutil.which("minced") or "minced"
        for d in ("spacers", "databases", "output", "tmp"):
            (self.wd / d).mkdir(parents=True, exist_ok=True)
        self.sd.mkdir(parents=True, exist_ok=True)

    def extract_spacers(
        self,
        fasta_path: str | Path,
        min_n_spacers: int,
        min_length: int,
        max_length: int,
        tag: str = "spacers_CompleteCollection",
    ) -> Path:
        """Extract CRISPR spacers from a FASTA using MinCED.

        Args:
            fasta_path: Input contigs FASTA.
            min_n_spacers: Minimum number of spacers in a CRISPR array (-minNR).
            min_length: Minimum spacer repeat length (-minRL).
            max_length: Maximum spacer repeat length (-maxRL).
            tag: Output file name prefix.

        Returns:
            Path to the spacer FASTA file.
        """
        raw_txt = self.sd / "spacers" / f"{tag}.txt"
        raw_gff = self.sd / "spacers" / f"{tag}.gff"
        fa_out1 = self.sd / "spacers" / f"{tag}_spacers.fa"
        fa_out2 = self.sd / "spacers" / f"{tag}.fasta"

        (self.sd / "spacers").mkdir(parents=True, exist_ok=True)

        cmd = (
            f'"{self.minced_path}" -spacers -minNR {min_n_spacers} '
            f'-minRL {min_length} -maxRL {max_length} '
            f'"{fasta_path}" "{raw_txt}" "{raw_gff}"'
        )
        sh(cmd, "MinCED - Extracting CRISPR spacers")

        if fa_out1.exists():
            fa_out1.rename(fa_out2)
        return fa_out2

    def make_db(
        self,
        fasta_file: str | Path,
        dbname: str,
        is_spacer: bool = False,
        rev: bool = False,
    ) -> Path:
        """Create a SpacePHARER database from a FASTA file.

        Args:
            fasta_file: Input FASTA.
            dbname: Database name (placed in databases/ subdirectory).
            is_spacer: Add --extractorf-spacer 1 flag.
            rev: Add --reverse-fragments 1 flag.

        Returns:
            Path to the database.
        """
        db_path = self.wd / "databases" / dbname
        tmp_path = self.wd / "tmp"

        flag = ""
        if is_spacer:
            flag += " --extractorf-spacer 1"
        if rev:
            flag += " --reverse-fragments 1"

        cmd = (
            f'"{self.spacepharer_path}" createsetdb '
            f'"{fasta_file}" "{db_path}" "{tmp_path}"{flag}'
        )
        sh(cmd, f"Creating DB {dbname}")
        return db_path

    def predict(
        self,
        spacerDB: Path,
        viralDB: Path,
        viralctrlDB: Path,
        out: str = "phage_host_predictions.tsv",
        fdr: float = 0.05,
    ) -> Path:
        """Run SpacePHARER predictmatch.

        Args:
            spacerDB: Path to spacer database.
            viralDB: Path to viral database.
            viralctrlDB: Path to viral control database.
            out: Output TSV filename.
            fdr: False discovery rate threshold.

        Returns:
            Path to the prediction TSV.
        """
        out_path = self.wd / "output" / out
        tmp_path = self.wd / "tmp"
        cmd = (
            f'"{self.spacepharer_path}" predictmatch '
            f'"{spacerDB}" "{viralDB}" "{viralctrlDB}" "{out_path}" "{tmp_path}" '
            f"--fdr {fdr}"
        )
        sh(cmd, "Running SpacePHARER predictmatch")
        return out_path

    def quick_stats(self, tsv: str | Path) -> None:
        """Print quick interaction statistics from a prediction TSV.

        Args:
            tsv: Path to phage_host_predictions.tsv.
        """
        with open(str(tsv)) as fh:
            lines = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
        hosts = len({ln.split("\t")[0] for ln in lines})
        phages = len({ln.split("\t")[1] for ln in lines})
        print(f"Interactions: {len(lines)}  Unique Hosts: {hosts}  Unique Phages: {phages}")


def get_spacers_collection(cfg: CapelliniConfig, wf: SpacePHARERWorkflow) -> Path:
    """Return path to spacers_CompleteCollection.fasta, using the bundled file if present.

    Modification 2: checks the bundled FASTA first. If not present, downloads
    progenomes3.contigs.representatives.fasta.bz2, decompresses it, runs MinCED,
    and optionally removes the decompressed FASTA.

    Args:
        cfg: Populated CapelliniConfig instance.
        wf: Initialized SpacePHARERWorkflow instance.

    Returns:
        Path to spacers_CompleteCollection.fasta.
    """
    spacers_path = Path(cfg.download_path) / "spacers"
    spacers_collection = spacers_path / "spacers_CompleteCollection.fasta"

    # Modification 2: bundled file takes priority (unless user asked to regenerate)
    if _BUNDLED_SPACERS.exists() and not getattr(cfg, "regenerate_spacers_collection", False):
        logger.info("Bundled spacers_CompleteCollection.fasta found — using it")
        spacers_path.mkdir(parents=True, exist_ok=True)
        if not spacers_collection.exists():
            import shutil as _sh
            _sh.copy2(str(_BUNDLED_SPACERS), str(spacers_collection))
        return spacers_collection

    if spacers_collection.exists():
        logger.info("Complete spacers collection found — skipping generation")
        return spacers_collection

    # Download, decompress, extract spacers
    filename = os.path.basename(cfg.bacContigs_reference_url)
    bacContigs_reference_path = Path(cfg.download_path) / filename
    bz2_path = bacContigs_reference_path
    fasta_path = bz2_path.with_suffix("")

    if not bacContigs_reference_path.exists():
        logger.info("Downloading ProGenomes3 contigs reference (~45 GB) ...")
        urllib.request.urlretrieve(cfg.bacContigs_reference_url, str(bacContigs_reference_path))

    if bz2_path.suffix == ".bz2" and not fasta_path.exists():
        logger.info("Decompressing bz2 FASTA reference ...")
        with bz2.open(str(bz2_path), "rb") as fin, open(str(fasta_path), "wb") as fout:
            shutil.copyfileobj(fin, fout)

    spacers_path.mkdir(parents=True, exist_ok=True)
    logger.info("Extracting spacers with MinCED ...")
    wf.extract_spacers(
        fasta_path=str(fasta_path),
        min_n_spacers=cfg.min_n_spacers,
        min_length=cfg.min_length,
        max_length=cfg.max_length,
        tag="spacers_CompleteCollection",
    )

    if cfg.remove_decomp_fasta and fasta_path.exists():
        fasta_path.unlink()
        logger.info("Removed decompressed FASTA: %s", fasta_path)

    return spacers_collection


def filter_target_spacers(
    spacers_collection: Path,
    ncbi_id_target_set_int: set,
    input_fasta_folder: str | Path,
) -> Path:
    """Filter spacers_CompleteCollection.fasta to cohort-specific NCBI IDs.

    Args:
        spacers_collection: Path to the complete spacers FASTA.
        ncbi_id_target_set_int: Set of integer NCBI IDs to keep.
        input_fasta_folder: Directory where target_spacers.fasta will be written.

    Returns:
        Path to the filtered target_spacers.fasta.
    """
    target_spacers_path = Path(input_fasta_folder) / "target_spacers.fasta"
    ncbi_id_target_set_str = set(str(i) for i in ncbi_id_target_set_int)

    logger.info("Filtering spacers to %s target NCBI IDs", len(ncbi_id_target_set_str))
    with open(str(target_spacers_path), "w") as fasta_out:
        for record in SeqIO.parse(str(spacers_collection), "fasta"):
            ncbi_id = record.id.split(".")[0]
            if ncbi_id in ncbi_id_target_set_str:
                SeqIO.write(record, fasta_out, "fasta")

    logger.info("Filtered spacers written: %s", target_spacers_path)
    return target_spacers_path


def run_spacepharer(cfg: CapelliniConfig, silva_fixed: pd.DataFrame) -> None:
    """Run the full SpacePHARER stage: spacer collection, filtering, and prediction.

    Args:
        cfg: Populated CapelliniConfig instance.
        silva_fixed: Output of the MMSeqs2 stage with progenomes_taxid and GCA columns.
    """
    logger.info("SpacePHARER: starting")

    # Virus FASTA check
    virus_fasta_path = Path(cfg.input_fasta_folder) / cfg.virus_fasta_name
    if not virus_fasta_path.exists():
        raise FileNotFoundError(f"Virus FASTA not found: {virus_fasta_path}")

    check_and_install_spacepharer()

    spacers_path = Path(cfg.download_path) / "spacers"
    wf = SpacePHARERWorkflow(workdir=cfg.sp_folder, spacerdir=str(spacers_path))

    # Build target sets
    if cfg.species_level:
        ncbi_id_target_set_int = (
            set(silva_fixed["progenomes_taxid_species"].dropna().astype(int))
            | set(silva_fixed["progenomes_taxid_family"].dropna().astype(int))
        )
        gca_target_set = (
            set(silva_fixed["GCA_species"].dropna())
            | set(silva_fixed["GCA_family"].dropna())
        )
    else:
        ncbi_id_target_set_int = (
            set(silva_fixed["progenomes_taxid_genus"].dropna().astype(int))
            | set(silva_fixed["progenomes_taxid_family"].dropna().astype(int))
        )
        gca_target_set = (
            set(silva_fixed["GCA_genus"].dropna())
            | set(silva_fixed["GCA_family"].dropna())
        )

    print(f"Resolution: {'species' if cfg.species_level else 'genus'}-level")
    print(f"Total unique NCBI target IDs: {len(ncbi_id_target_set_int)}")
    print(f"Total unique GCA target IDs: {len(gca_target_set)}")

    # Get spacers collection (Modification 2)
    spacers_collection = get_spacers_collection(cfg, wf)

    # Filter to target taxa
    target_spacers_path = filter_target_spacers(
        spacers_collection, ncbi_id_target_set_int, cfg.input_fasta_folder
    )

    # Run SpacePHARER
    spacerDB = wf.make_db(target_spacers_path, "spacerDB", is_spacer=True)
    viralDB = wf.make_db(virus_fasta_path, "viralDB")
    controlDB = wf.make_db(virus_fasta_path, "viralDB_control", rev=True)

    tsv = wf.predict(spacerDB, viralDB, controlDB, fdr=cfg.fdr)
    wf.quick_stats(tsv)

    logger.info("SpacePHARER stage complete")


def compute_spacepharer_stats(
    cfg: CapelliniConfig,
    silva_fixed: pd.DataFrame,
    topBitScore_df: pd.DataFrame,
) -> None:
    """Print the full pipeline statistics summary (Sections 1-3 from notebook).

    Args:
        cfg: Populated CapelliniConfig instance.
        silva_fixed: Output DataFrame with NCBI taxid columns.
        topBitScore_df: MMSeqs2 hits DataFrame.
    """
    sp_folder = Path(cfg.sp_folder)
    virus_fasta_path = Path(cfg.input_fasta_folder) / cfg.virus_fasta_name
    spacers_collection = Path(cfg.download_path) / "spacers" / "spacers_CompleteCollection.fasta"
    reference_16S_path = Path(cfg.input_fasta_folder) / "progenome16S.fasta"

    n_viral_contigs = sum(1 for _ in SeqIO.parse(str(virus_fasta_path), "fasta"))

    n_16s_records = 0
    unique_16s_ncbi_ids: set = set()
    with open(str(reference_16S_path), "r") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            n_16s_records += 1
            ncbi_id = rec.id.split(".")[0]
            if ncbi_id.isdigit():
                unique_16s_ncbi_ids.add(ncbi_id)
    n_16s_ncbi_ids = len(unique_16s_ncbi_ids)

    n_complete_spacers = sum(1 for _ in SeqIO.parse(str(spacers_collection), "fasta")) if spacers_collection.exists() else 0

    if cfg.species_level:
        ncbi_id_target_set_int = (
            set(silva_fixed["progenomes_taxid_species"].dropna().astype(int))
            | set(silva_fixed["progenomes_taxid_family"].dropna().astype(int))
        )
    else:
        ncbi_id_target_set_int = (
            set(silva_fixed["progenomes_taxid_genus"].dropna().astype(int))
            | set(silva_fixed["progenomes_taxid_family"].dropna().astype(int))
        )

    ncbi_id_target_set_str = set(str(i) for i in ncbi_id_target_set_int)
    target_hits = topBitScore_df[topBitScore_df["NCBI ID"].isin(ncbi_id_target_set_str)]
    n_bac_genomes = target_hits["Genome Accession ID"].dropna().nunique()
    n_bac_genomes_ncbi = target_hits["NCBI ID"].dropna().nunique()

    target_spacers_path = Path(cfg.input_fasta_folder) / "target_spacers.fasta"
    n_bac_spacers = 0
    ncbi_ids_with_spacers: set = set()
    if target_spacers_path.exists():
        for record in SeqIO.parse(str(target_spacers_path), "fasta"):
            n_bac_spacers += 1
            ncbi_ids_with_spacers.add(record.id.split(".")[0])
    n_ncbi_with_spacers = len(ncbi_ids_with_spacers)

    n_asvs_total = len(silva_fixed)
    n_asvs_matched = topBitScore_df["query"].nunique()
    pct_matched = 100 * n_asvs_matched / n_asvs_total
    mean_bitscore = topBitScore_df["bitscore"].mean()
    median_bitscore = topBitScore_df["bitscore"].median()

    n_assigned_species = silva_fixed["progenomes_taxid_species"].notna().sum()
    n_assigned_genus = silva_fixed["progenomes_taxid_genus"].notna().sum()
    n_assigned_family = silva_fixed["progenomes_taxid_family"].notna().sum()
    n_assigned_any = silva_fixed[
        ["progenomes_taxid_species", "progenomes_taxid_genus", "progenomes_taxid_family"]
    ].notna().any(axis=1).sum()

    # SpacePHARER results
    pred_path = sp_folder / "output" / "phage_host_predictions.tsv"
    colnames_9 = ["spacer_id", "phage_id", "pvalue", "aln_len", "mismatch",
                  "qstart", "qend", "qprot_aln", "sprot_aln"]
    sp_df = pd.read_csv(
        str(pred_path), sep="\t", comment="#", header=None, names=colnames_9,
        engine="python", dtype=str, on_bad_lines="skip"
    ).dropna(how="all").reset_index(drop=True)
    for c in ["pvalue", "aln_len", "mismatch", "qstart", "qend"]:
        sp_df[c] = pd.to_numeric(sp_df[c], errors="coerce")
    sp_df["host_id"] = sp_df["spacer_id"].str.replace(r"_spacer_\d+$", "", regex=True)

    pairs = sp_df[["host_id", "phage_id"]].drop_duplicates().reset_index(drop=True)
    pairs = pairs.merge(
        sp_df.groupby(["host_id", "phage_id"])["spacer_id"].nunique().rename("n_spacers").reset_index(),
        on=["host_id", "phage_id"]
    )
    pairs = pairs.merge(
        sp_df.groupby(["host_id", "phage_id"])["pvalue"].min().rename("best_pval").reset_index(),
        on=["host_id", "phage_id"]
    )

    host_degree_sp = pairs.groupby("host_id")["phage_id"].nunique().rename("degree")
    virus_degree_sp = pairs.groupby("phage_id")["host_id"].nunique().rename("degree")
    n_hosts = host_degree_sp.shape[0]
    n_viruses = virus_degree_sp.shape[0]
    n_interactions = pairs.shape[0]
    n_spacer_hits = len(sp_df)
    density = n_interactions / (n_hosts * n_viruses) if n_hosts * n_viruses > 0 else 0.0

    print("═" * 64)
    print("  PIPELINE STATISTICS SUMMARY")
    print("═" * 64)
    print(f"\n  [INPUT]")
    print(f"  Total 16S ASVs                          : {n_asvs_total:>8,}")
    print(f"  Viral contigs                           : {n_viral_contigs:>8,}")
    print(f"\n  [REFERENCE FILES]")
    print(f"  16S ProGenomes ref records              : {n_16s_records:>8,}")
    print(f"  16S ProGenomes ref unique NCBI IDs      : {n_16s_ncbi_ids:>8,}")
    print(f"  ProGenomes spacers (complete collection): {n_complete_spacers:>8,}")
    print(f"\n  [TARGET SET]")
    print(f"  Bacterial genomes (unique GCAs)         : {n_bac_genomes:>8,}")
    print(f"  Bacterial genomes (unique NCBI IDs)     : {n_bac_genomes_ncbi:>8,}")
    print(f"  Bacterial spacers (cohort-specific)     : {n_bac_spacers:>8,}")
    print(f"\n  [MMSeqs2 — 16S vs ProGenomes3 16S]")
    print(f"  ASVs with ≥1 hit (bitscore ≥ {cfg.min_bitscore})      : {n_asvs_matched:>8,}  ({pct_matched:.1f}%)")
    print(f"  Mean / Median bitscore                  : {mean_bitscore:>6.1f}  / {median_bitscore:.1f}")
    print(f"\n  [NCBI TAXID ASSIGNMENT — 3-layer]")
    print(f"  ASVs assigned (any layer)               : {n_assigned_any:>8,}  ({100*n_assigned_any/n_asvs_total:.1f}%)")
    print(f"    Layer 1 — species                     : {n_assigned_species:>8,}")
    print(f"    Layer 2 — genus                       : {n_assigned_genus:>8,}")
    print(f"    Layer 3 — family                      : {n_assigned_family:>8,}")
    print(f"  NCBI IDs with ≥1 spacer                 : {n_ncbi_with_spacers:>8,}")
    print(f"\n  [SPACEPHARER NETWORK]")
    print(f"  Unique hosts (bacteria)                 : {n_hosts:>8,}")
    print(f"  Unique viruses (phages/contigs)         : {n_viruses:>8,}")
    print(f"  Unique host–phage interactions          : {n_interactions:>8,}")
    print(f"  Raw spacer-hit rows                     : {n_spacer_hits:>8,}")
    print(f"  Network density                         : {density:>8.5f}")
    print("═" * 64)


def plot_spacepharer_figures(cfg: CapelliniConfig) -> None:
    """Generate SpacePHARER network figures (Sections 4-5 from notebook).

    Only called when figures_display=True in the pipeline config.

    Args:
        cfg: Populated CapelliniConfig instance.
    """
    import matplotlib.ticker as ticker

    sp_folder = Path(cfg.sp_folder)
    FIG_DIR = sp_folder / "figures"
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.family": "Arial", "font.size": 8, "axes.titlesize": 9,
        "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "axes.linewidth": 0.8,
        "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
    })
    PALETTE = {"host": "#E07B39", "virus": "#4A90C4", "edge": "#AAAAAA"}

    def save_fig(name: str) -> None:
        for ext in ("pdf", "png"):
            plt.savefig(str(FIG_DIR / f"{name}.{ext}"))
        plt.close()

    pred_path = sp_folder / "output" / "phage_host_predictions.tsv"
    colnames_9 = ["spacer_id", "phage_id", "pvalue", "aln_len", "mismatch",
                  "qstart", "qend", "qprot_aln", "sprot_aln"]
    sp_df = pd.read_csv(
        str(pred_path), sep="\t", comment="#", header=None, names=colnames_9,
        engine="python", dtype=str, on_bad_lines="skip"
    ).dropna(how="all").reset_index(drop=True)
    for c in ["pvalue", "aln_len", "mismatch", "qstart", "qend"]:
        sp_df[c] = pd.to_numeric(sp_df[c], errors="coerce")
    sp_df["host_id"] = sp_df["spacer_id"].str.replace(r"_spacer_\d+$", "", regex=True)

    pairs = sp_df[["host_id", "phage_id"]].drop_duplicates().reset_index(drop=True)
    pairs = pairs.merge(
        sp_df.groupby(["host_id", "phage_id"])["spacer_id"].nunique().rename("n_spacers").reset_index(),
        on=["host_id", "phage_id"]
    )
    pairs = pairs.merge(
        sp_df.groupby(["host_id", "phage_id"])["pvalue"].min().rename("best_pval").reset_index(),
        on=["host_id", "phage_id"]
    )

    host_degree_sp = pairs.groupby("host_id")["phage_id"].nunique().rename("degree")
    virus_degree_sp = pairs.groupby("phage_id")["host_id"].nunique().rename("degree")

    # Figure 1: degree distributions
    fig, axes = plt.subplots(1, 2, figsize=(6, 2))
    for ax, data, label, color in [
        (axes[0], host_degree_sp, "Host degree (phages per bacterium)", PALETTE["host"]),
        (axes[1], virus_degree_sp, "Virus degree (bacteria per phage)", PALETTE["virus"]),
    ]:
        bins = np.logspace(np.log10(max(data.min(), 1)), np.log10(data.max()), 30)
        ax.hist(data, bins=bins, color=color, edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.set_xlabel(label)
        ax.set_ylabel("Count")
    fig.suptitle("Degree distributions", fontsize=9, y=1.02)
    plt.tight_layout()
    save_fig("01_degree_distributions")

    # Figure 2: spacer support
    fig, ax = plt.subplots(figsize=(1.75, 1.4))
    ax.hist(pairs["n_spacers"], bins=30, color=PALETTE["host"], edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Unique spacers per interaction")
    ax.set_ylabel("Count")
    ax.set_title("Spacer support per host–phage pair")
    save_fig("02_spacer_support")

    logger.info("SpacePHARER figures saved to: %s", FIG_DIR)
