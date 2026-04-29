"""Network stage: build common-abundance, shrinkage, CRISPR, smoothed, and X* outputs."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from capellini.config import CapelliniConfig
from capellini.utils.io import read_table, write_df
from capellini.utils.taxonomy import (
    clean_df_ids,
    clean_index_ids,
    load_bacteria_taxonomy,
)
from capellini.utils.transforms import (
    double_clr_transform,
    old_style_clr_transform,
    schaefer_strimmer_corr,
)
from capellini.utils.network_utils import (
    align_abundance_from_metadata,
    build_smoothed_crispr_for_study,
    build_xstar_from_smoothed_crispr,
    build_xstar_from_smoothed_crispr_residual,
    crispr_matrix_aggregate_viruses,
    orient_W_viruses_by_bacteria,
    prepare_bacteria_genus_abundance,
    prevalence_filter_df,
    remove_disease_columns_from_virus_abundance,
    study_outdir,
)

logger = logging.getLogger(__name__)


def _load_metadata(cfg: CapelliniConfig) -> pd.DataFrame:
    return read_table(cfg.metadata_path, index_col=None)


def build_common_abundance_one(study: str, cfg: CapelliniConfig) -> dict[str, Path]:
    """Build the aligned, prevalence-filtered common abundance tables.

    Args:
        study: Study identifier.
        cfg: Populated CapelliniConfig.

    Returns:
        Mapping of artefact name to written file path.
    """
    logger.info("[%s] common abundance: loading raw inputs", study)
    V_raw = read_table(cfg.virus_abundance_raw)
    V_raw = remove_disease_columns_from_virus_abundance(V_raw)
    V_raw = clean_df_ids(V_raw)

    B_otu = read_table(cfg.bacteria_otu)
    B_tax = load_bacteria_taxonomy(cfg.bacteria_taxonomy)
    B_genus = prepare_bacteria_genus_abundance(B_otu, B_tax, rank=cfg.BACTERIA_TAXONOMY_RANK)

    meta = _load_metadata(cfg)
    V, B, meta_aligned = align_abundance_from_metadata(
        virus_abundance=V_raw,
        bacteria_abundance=B_genus,
        metadata=meta,
        keep_col=cfg.KEEP_COLUMN,
    )

    V = prevalence_filter_df(V, prevalence=cfg.PREVALENCE)
    B = prevalence_filter_df(B, prevalence=cfg.PREVALENCE)

    out_dir = study_outdir(study, "common", Path(cfg.OUTPUT_ROOT))
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "V": out_dir / "V_processed.csv",
        "B": out_dir / "B_processed.csv",
        "meta": out_dir / "metadata_aligned.csv",
    }
    write_df(V, paths["V"])
    write_df(B, paths["B"])
    write_df(meta_aligned, paths["meta"])
    logger.info("[%s] common abundance written to %s", study, out_dir)
    return paths


def build_shrinkage_one(study: str, cfg: CapelliniConfig) -> dict[str, Path]:
    """Compute shrinkage correlations on the bacteria-virus joint matrix."""
    out_dir = study_outdir(study, "shrinkage", Path(cfg.OUTPUT_ROOT))
    out_dir.mkdir(parents=True, exist_ok=True)
    common = study_outdir(study, "common", Path(cfg.OUTPUT_ROOT))
    V = read_table(common / "V_processed.csv")
    B = read_table(common / "B_processed.csv")

    V = V.loc[B.index]
    V_pref = V.add_prefix("V__")
    B_pref = B.add_prefix("B__")
    joint = pd.concat([B_pref, V_pref], axis=1)
    joint_clr = old_style_clr_transform(joint)
    corr_df, info = schaefer_strimmer_corr(joint_clr)
    bv = corr_df.loc[B_pref.columns, V_pref.columns]

    paths = {
        "corr_full": out_dir / "shrinkage_corr_full.csv.gz",
        "corr_bv": out_dir / "shrinkage_corr_BV.csv.gz",
    }
    write_df(corr_df, paths["corr_full"])
    write_df(bv, paths["corr_bv"])
    logger.info("[%s] shrinkage written (lambda=%s)", study, info)
    return paths


def build_raw_crispr_one(study: str, cfg: CapelliniConfig) -> Path:
    """Aggregate phage-host predictions into a raw CRISPR network."""
    out_dir = study_outdir(study, "crispr_raw", Path(cfg.OUTPUT_ROOT))
    out_dir.mkdir(parents=True, exist_ok=True)
    df_pred = pd.read_csv(cfg.phage_host_predictions, sep="\t", header=None, comment="#")
    vir_tax = read_table(cfg.tax_vir)
    crispr = crispr_matrix_aggregate_viruses(
        df_pred, vir_tax, vir_rank=cfg.aggregate_viral_rank
    )
    out = out_dir / "crispr_net.csv.gz"
    write_df(crispr, out)
    logger.info("[%s] raw CRISPR matrix written: %s", study, out)
    return out


def build_smooth_crispr_one(study: str, cfg: CapelliniConfig) -> dict[str, Path]:
    """Smooth the raw CRISPR matrix using bacteria/virus taxonomy kernels."""
    out_dir = study_outdir(study, "crispr_smooth", Path(cfg.OUTPUT_ROOT))
    out_dir.mkdir(parents=True, exist_ok=True)
    common = study_outdir(study, "common", Path(cfg.OUTPUT_ROOT))
    V = read_table(common / "V_processed.csv")
    B = read_table(common / "B_processed.csv")
    crispr_path = study_outdir(study, "crispr_raw", Path(cfg.OUTPUT_ROOT)) / "crispr_net.csv.gz"

    artefacts = build_smoothed_crispr_for_study(
        raw_crispr_path=str(crispr_path),
        bacteria_features=B.columns,
        virus_features=V.columns,
        tax_bac_path=cfg.tax_bac_for_smoothing,
        tax_vir_path=cfg.tax_vir,
        bacterial_ranks=cfg.BACTERIAL_RANKS,
        viral_ranks=cfg.viral_ranks,
        bacterial_weights=cfg.BACTERIAL_WEIGHTS,
        viral_weights=cfg.viral_weights,
        alpha=cfg.CRISPR_SMOOTH_ALPHA,
        transpose_after_load=cfg.TRANSPOSE_RAW_CRISPR_AFTER_LOAD,
    )
    paths = {}
    for name in ("crispr_binary", "crispr_smooth", "K_bac", "K_vir"):
        p = out_dir / f"{name}.csv.gz"
        write_df(artefacts[name], p)
        paths[name] = p
    # Also write the vir × bac orientation used by the X* stage
    p_vb = out_dir / "crispr_smooth_vir_bac.csv.gz"
    write_df(artefacts["crispr_smooth"].T, p_vb)
    paths["crispr_smooth_vir_bac"] = p_vb
    logger.info("[%s] smoothed CRISPR written to %s", study, out_dir)
    return paths


def build_xstar_one(study: str, cfg: CapelliniConfig) -> dict[str, Path]:
    """Build the CLR-transformed message-passing X* outputs."""
    out_dir = study_outdir(study, "xstar", Path(cfg.OUTPUT_ROOT))
    out_dir.mkdir(parents=True, exist_ok=True)
    common = study_outdir(study, "common", Path(cfg.OUTPUT_ROOT))
    smooth = study_outdir(study, "crispr_smooth", Path(cfg.OUTPUT_ROOT))

    V = read_table(common / "V_processed.csv")
    B = read_table(common / "B_processed.csv")
    W_vb_path = smooth / "crispr_smooth_vir_bac.csv.gz"
    if W_vb_path.exists():
        W = read_table(W_vb_path)
    else:
        # Fallback: derive vir × bac from the bac × vir file
        W = read_table(smooth / "crispr_smooth.csv.gz").T
    W = orient_W_viruses_by_bacteria(W, V, B)

    conv = build_xstar_from_smoothed_crispr(
        V_df=V,
        B_df=B,
        W_vh_smooth_df=W,
        pseudocount=cfg.PSEUDOCOUNT,
        lam=cfg.LAM,
        n_steps=cfg.N_STEPS,
        preserve_scale=cfg.PRESERVE_SCALE,
    )
    resid = build_xstar_from_smoothed_crispr_residual(
        V_df=V,
        B_df=B,
        W_vh_smooth_df=W,
        pseudocount=cfg.PSEUDOCOUNT,
        lam=cfg.LAM,
        n_steps=cfg.N_STEPS,
        preserve_scale=cfg.PRESERVE_SCALE,
    )

    paths = {
        "V_clr": out_dir / "V_clr.csv.gz",
        "B_clr": out_dir / "B_clr.csv.gz",
        "X_star_convex": out_dir / "X_star_convex.csv.gz",
        "X_star_residual": out_dir / "X_star_residual.csv.gz",
    }
    write_df(conv["V_clr"], paths["V_clr"])
    write_df(conv["B_clr"], paths["B_clr"])
    write_df(conv["X_star"], paths["X_star_convex"])
    write_df(resid["X_star"], paths["X_star_residual"])
    logger.info("[%s] X* outputs written to %s", study, out_dir)
    return paths


def run_network(cfg: CapelliniConfig) -> dict[str, dict]:
    """Run the network stage according to the RUN_* flags on the config.

    Args:
        cfg: Populated CapelliniConfig instance.

    Returns:
        Dict mapping each enabled sub-stage name to its outputs.
    """
    if not cfg.OUTPUT_ROOT:
        raise ValueError("Network stage requires cfg.OUTPUT_ROOT (or enhanced_networks_folder)")

    Path(cfg.OUTPUT_ROOT).mkdir(parents=True, exist_ok=True)
    study = cfg.STUDY or "default"
    results: dict[str, dict] = {}

    if cfg.RUN_COMMON_ABUNDANCE:
        results["common"] = build_common_abundance_one(study, cfg)
    if cfg.RUN_SHRINKAGE_CORRELATIONS:
        results["shrinkage"] = build_shrinkage_one(study, cfg)
    if cfg.RUN_RAW_CRISPR_NETWORKS:
        results["crispr_raw"] = {"crispr_net": build_raw_crispr_one(study, cfg)}
    if cfg.RUN_SMOOTH_CRISPR:
        results["crispr_smooth"] = build_smooth_crispr_one(study, cfg)
    if cfg.RUN_XSTAR:
        results["xstar"] = build_xstar_one(study, cfg)

    logger.info("Network stage complete (%d sub-stages)", len(results))
    return results
