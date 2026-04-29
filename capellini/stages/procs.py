"""ProCs stage: bacterial/viral protein extraction, clustering, and PA matrix."""

from __future__ import annotations

import bz2
import logging
import os
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from capellini.config import CapelliniConfig
from capellini.utils.io import sh

logger = logging.getLogger(__name__)


def download_protein_reference(cfg: CapelliniConfig) -> Path:
    """Download the ProGenomes3 proteins FASTA if not already present.

    Args:
        cfg: Populated CapelliniConfig instance.

    Returns:
        Path to the downloaded bz2 protein reference.
    """
    filename = os.path.basename(cfg.protein_reference_url)
    protein_reference_path = Path(cfg.download_path) / filename

    if protein_reference_path.exists():
        logger.info("Protein reference found — skipping download")
    else:
        logger.info("Downloading ProGenomes3 protein reference ...")
        urllib.request.urlretrieve(cfg.protein_reference_url, str(protein_reference_path))

    return protein_reference_path


def extract_bacterial_proteins(cfg: CapelliniConfig, gca_target_set: set) -> Path:
    """Stream the bz2 protein FASTA and extract proteins for target GCAs.

    Uses a batch approach when len(gca_target_set) > cfg.batch_size to avoid
    holding all sequences in memory at once.

    Args:
        cfg: Populated CapelliniConfig instance.
        gca_target_set: Set of GCA IDs (e.g. 'GCA_000001405') to extract.

    Returns:
        Path to BacterialProteinsCollection.fasta.
    """
    protein_reference_path = download_protein_reference(cfg)
    bac_protcoll_fasta_path = Path(cfg.proteins_extraction_path) / "BacterialProteinsCollection.fasta"
    single_bacgenome_collection_path = Path(cfg.proteins_extraction_path) / "BacterialSingleGenomesCollections"

    if cfg.save_single_bacgenome_collection:
        single_bacgenome_collection_path.mkdir(parents=True, exist_ok=True)

    # Initialize output file
    with open(str(bac_protcoll_fasta_path), "w"):
        pass

    proteins_genomes: dict = {}
    protein_counter: dict = {}
    batch_counter = 0
    targets_list = gca_target_set

    entry_last_batch_checkpoint: dict = {genome_name: 0 for genome_name in targets_list}

    if len(targets_list) > cfg.batch_size:
        progress_bar = tqdm(total=cfg.batch_size, desc="Filling Progenome Genome IDs batch")
    else:
        progress_bar = tqdm(total=len(targets_list), desc="Processing Genomes assignment")

    with bz2.open(str(protein_reference_path), "rt") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            genome_name_parts = record.description.split()[0].split(".")
            genome_name_elements = genome_name_parts[2].split("_")
            genome_name = f"{genome_name_elements[0]}_{genome_name_elements[1]}"

            if genome_name in targets_list:
                seq = str(record.seq)
                if genome_name not in proteins_genomes:
                    proteins_genomes[genome_name] = set()
                    protein_counter[genome_name] = 1
                    progress_bar.update(1)
                if seq not in proteins_genomes[genome_name]:
                    proteins_genomes[genome_name].add(seq)
                    protein_counter[genome_name] += 1

            if len(targets_list) > cfg.batch_size and len(proteins_genomes.keys()) >= cfg.batch_size:
                progress_bar.close()
                batch_counter += 1
                logger.debug("Batch %s filled, saving FASTA files", batch_counter)

                _flush_batch(
                    proteins_genomes,
                    entry_last_batch_checkpoint,
                    bac_protcoll_fasta_path,
                    single_bacgenome_collection_path if cfg.save_single_bacgenome_collection else None,
                )
                proteins_genomes.clear()

                tqdm.write(f"Batch {batch_counter} completed.")
                progress_bar = tqdm(total=cfg.batch_size, desc="Filling Progenome Genome IDs batch")

    progress_bar.close()

    # Final flush
    _flush_batch(
        proteins_genomes,
        entry_last_batch_checkpoint,
        bac_protcoll_fasta_path,
        single_bacgenome_collection_path if cfg.save_single_bacgenome_collection else None,
    )

    logger.info("BacterialProteinsCollection written: %s", bac_protcoll_fasta_path)

    if cfg.ref_removal:
        try:
            os.remove(str(protein_reference_path))
            logger.info("Deleted protein reference: %s", protein_reference_path.name)
        except FileNotFoundError:
            pass

    return bac_protcoll_fasta_path


def _flush_batch(
    proteins_genomes: dict,
    entry_last_batch_checkpoint: dict,
    bac_protcoll_fasta_path: Path,
    single_genome_path=None,
) -> None:
    """Write the current batch of proteins to disk.

    Args:
        proteins_genomes: Dict genome_name -> set of protein sequences.
        entry_last_batch_checkpoint: Per-genome protein counter checkpoint.
        bac_protcoll_fasta_path: Combined output FASTA path.
        single_genome_path: Optional path for per-genome FASTAs.
    """
    for genome_accession, prot_seqs in list(proteins_genomes.items()):
        i = entry_last_batch_checkpoint.get(genome_accession, 0)
        for prot_seq in prot_seqs:
            i += 1
            entry_name = f">{genome_accession}_Protein{i}"
            if single_genome_path is not None:
                output_fasta_path = single_genome_path / f"{genome_accession}.fasta"
                with open(str(output_fasta_path), "a") as fh:
                    fh.write(f"{entry_name}\n{prot_seq}\n")
            with open(str(bac_protcoll_fasta_path), "a") as fh:
                fh.write(f"{entry_name}\n{prot_seq}\n")
        entry_last_batch_checkpoint[genome_accession] = i


def extract_viral_proteins(cfg: CapelliniConfig) -> Path:
    """Run prodigal in metagenomic mode on the virus FASTA to extract proteins.

    Args:
        cfg: Populated CapelliniConfig instance.

    Returns:
        Path to ViralProteinsCollection.fasta.
    """
    virus_fasta_path = Path(cfg.input_fasta_folder) / cfg.virus_fasta_name
    vir_protcoll_fasta_path = Path(cfg.proteins_extraction_path) / "ViralProteinsCollection.fasta"
    coords_path = Path(cfg.proteins_extraction_path) / "coords.gbk"

    cmd = (
        f"prodigal -i '{virus_fasta_path}' -o '{coords_path}' "
        f"-a '{vir_protcoll_fasta_path}' -p meta"
    )
    sh(cmd, "Running Prodigal")

    if not cfg.keep_coords and coords_path.exists():
        coords_path.unlink()
        logger.debug("Removed coords.gbk")

    logger.info("ViralProteinsCollection written: %s", vir_protcoll_fasta_path)
    return vir_protcoll_fasta_path


def combine_protein_collections(
    bac_path: Path,
    vir_path: Path,
    combined_path: Path,
) -> Path:
    """Concatenate bacterial and viral protein FASTAs into a single combined FASTA.

    Args:
        bac_path: BacterialProteinsCollection.fasta path.
        vir_path: ViralProteinsCollection.fasta path.
        combined_path: Destination CombinedProteinsCollection.fasta path.

    Returns:
        Path to the combined FASTA.
    """
    with open(str(combined_path), "w") as out_handle:
        for fasta_path in [bac_path, vir_path]:
            with open(str(fasta_path), "r") as in_handle:
                for record in SeqIO.parse(in_handle, "fasta"):
                    SeqIO.write(record, out_handle, "fasta")
    logger.info("CombinedProteinsCollection written: %s", combined_path)
    return combined_path


def run_mmseqs_clustering(combined_fasta_path: Path, clustering_path: str) -> Path:
    """Run mmseqs easy-cluster on the combined protein FASTA.

    Args:
        combined_fasta_path: CombinedProteinsCollection.fasta path.
        clustering_path: Directory for clustering outputs.

    Returns:
        Path to clusterRes (prefix; actual tsv is clusterRes_cluster.tsv).
    """
    prot_path = Path(combined_fasta_path)
    work_dir = prot_path.parent
    clustering_root = Path(clustering_path)

    cluster_prefix_rel = Path("..") / clustering_root.name / "clusterRes"
    tmp_dir_rel = Path("..") / clustering_root.name / "tmp"

    (work_dir / tmp_dir_rel).mkdir(parents=True, exist_ok=True)

    cmd = (
        f'cd "{work_dir}" && '
        f'mmseqs easy-cluster '
        f'"{prot_path.name}" '
        f'"{cluster_prefix_rel}" '
        f'"{tmp_dir_rel}"'
    )
    sh(cmd, "MMseqs2 - clustering proteins")
    return clustering_root / "clusterRes"


def build_pa_matrix(
    cluster_res_df: pd.DataFrame,
    filter_1bac_1vir: bool,
    vir_fasta: Path,
    bac_fasta: Path,
    matrix_type: str = "count",
) -> pd.DataFrame:
    """Build a presence/absence or count matrix of protein clusters per genome/virus.

    Args:
        cluster_res_df: DataFrame with Cluster and Protein columns.
        filter_1bac_1vir: If True, keep only clusters with ≥1 bacterial and ≥1 viral protein.
        vir_fasta: ViralProteinsCollection.fasta for filter_1bac_1vir logic.
        bac_fasta: BacterialProteinsCollection.fasta for filter_1bac_1vir logic.
        matrix_type: 'count' or 'binary'.

    Returns:
        Genomes/viruses x protein clusters matrix DataFrame.
    """
    if filter_1bac_1vir:
        logger.info("Applying 1-bacterial + 1-viral protein filter")
        vp_set = set(record.description for record in SeqIO.parse(str(vir_fasta), "fasta"))
        bp_set = set(record.description for record in SeqIO.parse(str(bac_fasta), "fasta"))

        cluster_res_df = cluster_res_df.copy()
        cluster_res_df["IsVirusProtein"] = cluster_res_df["Protein"].isin(vp_set)
        cluster_res_df["IsBacteriaProtein"] = cluster_res_df["Protein"].isin(bp_set)

        grouped = cluster_res_df.groupby("Cluster").agg({
            "IsVirusProtein": "any",
            "IsBacteriaProtein": "any",
        }).reset_index()

        clusters_to_keep = grouped[
            grouped["IsVirusProtein"] & grouped["IsBacteriaProtein"]
        ]["Cluster"]

        clusters_removal_count = len(grouped) - len(clusters_to_keep)
        cluster_res_df = cluster_res_df[cluster_res_df["Cluster"].isin(clusters_to_keep)]
        cluster_res_df.reset_index(drop=True, inplace=True)
        logger.info("Removed %s / %s clusters after filter", clusters_removal_count, len(grouped))

    count_p = Counter(cluster_res_df["Cluster"])
    bacterial_id: set = set()
    viral_id: set = set()
    dict_con_to_procs: dict = defaultdict(list)

    for ind in cluster_res_df.index:
        protein = cluster_res_df.loc[ind, "Protein"]
        procs = cluster_res_df.loc[ind, "Cluster"]
        if "GCA" in str(protein):
            bac_entry = protein.rsplit("_", 1)[0]
            bacterial_id.add(bac_entry)
            if count_p[procs] != 1:
                dict_con_to_procs[bac_entry].append(procs)
        else:
            viral_entry = protein.rsplit("_", 1)[0]
            viral_id.add(viral_entry)
            if count_p[procs] != 1:
                dict_con_to_procs[viral_entry].append(procs)

    s = pd.Series(dict_con_to_procs).explode()
    df = pd.crosstab(s.index, s)

    if matrix_type == "binary":
        df = (df > 0).astype(int)

    return df


def run_procs(cfg: CapelliniConfig, gca_target_set: set) -> pd.DataFrame:
    """Orchestrate the full ProCs stage.

    Steps:
        1. Extract bacterial proteins from ProGenomes3 bz2.
        2. Extract viral proteins with Prodigal.
        3. Combine into a single FASTA.
        4. Run mmseqs easy-cluster.
        5. Build PA matrix.

    Args:
        cfg: Populated CapelliniConfig instance.
        gca_target_set: Set of target GCA IDs from the MMSeqs2 stage.

    Returns:
        PA/count matrix DataFrame (genomes/viruses x protein clusters).
    """
    logger.info("ProCs: starting protein extraction and clustering")

    bac_protcoll = Path(cfg.proteins_extraction_path) / "BacterialProteinsCollection.fasta"
    vir_protcoll = Path(cfg.proteins_extraction_path) / "ViralProteinsCollection.fasta"
    combined_protcoll = Path(cfg.proteins_extraction_path) / "CombinedProteinsCollection.fasta"

    extract_bacterial_proteins(cfg, gca_target_set)
    extract_viral_proteins(cfg)
    combine_protein_collections(bac_protcoll, vir_protcoll, combined_protcoll)

    run_mmseqs_clustering(combined_protcoll, cfg.clustering_path)

    cluster_res_df = pd.read_table(
        f"{cfg.clustering_path}/clusterRes_cluster.tsv",
        header=None,
        names=["Cluster", "Protein"],
    )
    logger.info(
        "Clustering complete: %s proteins, %s clusters",
        cluster_res_df["Protein"].nunique(),
        cluster_res_df["Cluster"].nunique(),
    )

    pa_matrix = build_pa_matrix(
        cluster_res_df=cluster_res_df,
        filter_1bac_1vir=cfg.filter_1bac_1vir,
        vir_fasta=vir_protcoll,
        bac_fasta=bac_protcoll,
        matrix_type=cfg.matrix_type,
    )

    if cfg.remove_collections:
        for f in [bac_protcoll, vir_protcoll, combined_protcoll]:
            try:
                os.remove(str(f))
            except FileNotFoundError:
                pass

    logger.info("ProCs stage complete: PA matrix shape %s", pa_matrix.shape)
    return pa_matrix
