"""Configuration dataclass for the CAPELLINI pipeline."""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class CapelliniConfig:
    """All settings for the CAPELLINI pipeline, mirroring the notebook Settings section.

    Required fields (no defaults) must be provided explicitly or via from_yaml/from_dict.
    Derived path fields are computed in __post_init__ when left as empty strings.
    """

    # ── Global ────────────────────────────────────────────────────────────────
    base: str = ""
    download_path: str = ""

    # Derived folder paths (computed from base in __post_init__ if empty)
    input_fasta_folder: str = ""
    dada2_folder: str = ""
    mmseq_folder: str = ""
    sp_folder: str = ""
    procs_folder: str = ""
    enhanced_networks_folder: str = ""

    # Reference paths
    silva_ref_path: str = ""
    silva_taxmap_path: str = ""
    full_ncbi_taxonomy_path: str = ""
    virus_fasta_name: str = ""
    metadata_path: str = ""
    bacterial_raw_fasta_folder: str = ""

    species_level: bool = False
    fresh_start: bool = False
    ref_removal: bool = True

    # Force-regenerate bundled references (off by default — bundled files are reused)
    regenerate_16S_reference: bool = False
    regenerate_spacers_collection: bool = False

    # External download URLs (fixed)
    genes_reference_url: str = (
        "http://progenomes3.embl.de/data/repGenomes/"
        "progenomes3.genes.representatives.fasta.bz2"
    )
    bacContigs_reference_url: str = (
        "http://progenomes3.embl.de/data/repGenomes/"
        "progenomes3.contigs.representatives.fasta.bz2"
    )
    protein_reference_url: str = (
        "http://progenomes3.embl.de/data/repGenomes/"
        "progenomes3.proteins.representatives.fasta.bz2"
    )

    # ── DADA2 ─────────────────────────────────────────────────────────────────
    direction: str = "forward"
    bacteria_fasta_name: str = "16S_DADA2_bacteria.fasta"
    fasta_generation: bool = True

    # ── MMSeqs2 ───────────────────────────────────────────────────────────────
    isolate_ref_16S: bool = True
    mapping_saving: bool = True
    min_bitscore: int = 50
    max_matches: int = 20
    add_taxonomy: bool = True
    extend_taxonomy: bool = True

    # ── SpacePHARER ───────────────────────────────────────────────────────────
    min_n_spacers: int = 3
    min_length: int = 23
    max_length: int = 47
    fdr: float = 0.05
    keep_spacers_collection: bool = True
    remove_decomp_fasta: bool = True

    # ── ProCs ─────────────────────────────────────────────────────────────────
    proteins_extraction_path: str = ""
    clustering_path: str = ""
    matrix_type: str = "count"
    save_single_bacgenome_collection: bool = False
    keep_coords: bool = False
    filter_1bac_1vir: bool = False
    remove_collections: bool = False
    batch_size: int = 1500

    # ── Network ───────────────────────────────────────────────────────────────
    OUTPUT_ROOT: str = ""
    OVERWRITE: bool = False
    VERBOSE: bool = True
    RUN_COMMON_ABUNDANCE: bool = True
    RUN_SHRINKAGE_CORRELATIONS: bool = True
    RUN_RAW_CRISPR_NETWORKS: bool = True
    RUN_SMOOTH_CRISPR: bool = True
    RUN_XSTAR: bool = True
    PREVALENCE: float = 0.10
    KEEP_COLUMN: str = "keep_for_analysis"
    BACTERIA_TAXONOMY_RANK: str = "target_taxids"
    BACTERIAL_RANKS: list = field(
        default_factory=lambda: ["Phylum", "Class", "Order", "Family", "progenomes_taxid_genus"]
    )
    BACTERIAL_WEIGHTS: list = field(default_factory=lambda: [1, 2, 3, 6, 8])
    CRISPR_SMOOTH_ALPHA: float = 0.95
    TRANSPOSE_RAW_CRISPR_AFTER_LOAD: bool = True
    PSEUDOCOUNT: float = 1e-6
    LAM: float = 0.5
    N_STEPS: int = 1
    PRESERVE_SCALE: bool = False
    STUDY: str = "default"

    # Network input file paths
    virus_abundance_raw: str = ""
    bacteria_otu: str = ""
    bacteria_taxonomy: str = ""
    phage_host_predictions: str = ""
    tax_bac_for_smoothing: str = ""
    tax_vir: str = ""
    viral_ranks: list = field(
        default_factory=lambda: ["lev8", "lev7", "lev6", "lev5", "lev4", "lev3", "lev2", "lev1", "lev0"]
    )
    viral_weights: list = field(default_factory=lambda: [1, 1, 2, 3, 4, 6, 8, 10, 12])
    aggregate_viral_rank: str = "lev0"

    def __post_init__(self) -> None:
        if self.base:
            if not self.input_fasta_folder:
                self.input_fasta_folder = self.base + "/Inputs/Fasta Collection"
            if not self.dada2_folder:
                self.dada2_folder = self.base + "/DADA2 output"
            if not self.mmseq_folder:
                self.mmseq_folder = self.base + "/MMSeqs2 Output"
            if not self.sp_folder:
                self.sp_folder = self.base + "/SpacePHARER output"
            if not self.procs_folder:
                self.procs_folder = self.base + "/Procs Estimations"
            if not self.enhanced_networks_folder:
                self.enhanced_networks_folder = self.base + "/Enhanced Networks"
        if self.download_path:
            if not self.full_ncbi_taxonomy_path:
                self.full_ncbi_taxonomy_path = os.path.join(self.download_path, "names.dmp")
        if self.procs_folder:
            if not self.proteins_extraction_path:
                self.proteins_extraction_path = self.procs_folder + "/Targets Proteins Extraction/"
            if not self.clustering_path:
                self.clustering_path = self.procs_folder + "/Clustering"
        if self.enhanced_networks_folder and not self.OUTPUT_ROOT:
            self.OUTPUT_ROOT = self.enhanced_networks_folder
        if self.dada2_folder:
            suffix = {"forward": "F", "reverse": "R", "paired": "P"}.get(self.direction, "F")
            if not self.bacteria_taxonomy:
                self.bacteria_taxonomy = (
                    f"{self.dada2_folder}/taxonomy_table_{suffix}_progenomeLikeNCBIIDs.csv"
                )
            if not self.bacteria_otu:
                self.bacteria_otu = f"{self.dada2_folder}/OTU_table_{suffix}.csv"
        if self.sp_folder and not self.phage_host_predictions:
            self.phage_host_predictions = f"{self.sp_folder}/output/phage_host_predictions.tsv"
        if self.dada2_folder and not self.tax_bac_for_smoothing:
            suffix = {"forward": "F", "reverse": "R", "paired": "P"}.get(self.direction, "F")
            self.tax_bac_for_smoothing = (
                f"{self.dada2_folder}/taxonomy_table_{suffix}_progenomeLikeNCBIIDs.csv"
            )

    # ── Class methods ─────────────────────────────────────────────────────────

    @classmethod
    def default(cls) -> "CapelliniConfig":
        """Return a config with all default values (paths will be empty until base is set)."""
        return cls()

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "CapelliniConfig":
        """Load config from a plain Python dict."""
        known = {f.name for f in cls.__dataclass_fields__.values()}  # type: ignore[attr-defined]
        filtered = {k: v for k, v in d.items() if k in known}
        return cls(**filtered)

    @classmethod
    def from_yaml(cls, path: str | Path) -> "CapelliniConfig":
        """Load config from a YAML file.

        Args:
            path: Path to the YAML configuration file.

        Returns:
            CapelliniConfig populated from the YAML.
        """
        with open(path, "r") as fh:
            data = yaml.safe_load(fh) or {}
        return cls.from_dict(data)

    def to_yaml(self, path: str | Path) -> None:
        """Serialize config to a YAML file.

        Args:
            path: Destination YAML file path.
        """
        import dataclasses

        data = dataclasses.asdict(self)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            yaml.dump(data, fh, default_flow_style=False, allow_unicode=True, sort_keys=False)

    def virus_fasta_path(self) -> Path:
        """Resolved absolute path to the virus FASTA file."""
        return Path(self.input_fasta_folder) / self.virus_fasta_name
