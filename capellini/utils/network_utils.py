"""Network-level utilities: message passing, CRISPR smoothing, taxonomy kernels, abundance helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from capellini.utils.transforms import (
    closure,
    clr,
    double_clr_transform,
    row_normalize,
)
from capellini.utils.taxonomy import (
    sanitize_taxon_name,
    sanitize_index,
    clean_index_ids,
    clean_df_ids,
    parse_bool_series,
    load_bacteria_taxonomy,
    clean_bacteria_taxonomy,
    apply_custom_renames,
    rename_clostridium_sensu_stricto,
)
from capellini.utils.io import read_table, write_df


# ── Validation helpers ─────────────────────────────────────────────────────────

def _validate_inputs(V: np.ndarray, B: np.ndarray, W: np.ndarray, lam: float, n_steps: int) -> None:
    """Validate shapes and parameter ranges for message-passing functions.

    Args:
        V: Samples x viruses matrix.
        B: Samples x bacteria matrix.
        W: Viruses x bacteria CRISPR matrix.
        lam: Mixing weight.
        n_steps: Number of propagation steps.

    Raises:
        ValueError: On shape mismatch or invalid parameter values.
    """
    n, q = V.shape
    n_b, p = B.shape
    if n != n_b:
        raise ValueError(f"V and B must have the same number of samples, got {n} and {n_b}.")
    if W.shape != (q, p):
        raise ValueError(f"W_vh_smooth must have shape ({q}, {p}), got {W.shape}.")
    if lam < 0:
        raise ValueError(f"lam must be nonnegative, got {lam}.")
    if n_steps < 1:
        raise ValueError(f"n_steps must be at least 1, got {n_steps}.")


def _build_common_inputs(
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    pseudocount: float,
    eps: float,
) -> tuple:
    """Shared preprocessing for X-star builders: align samples, apply double-CLR, subset W.

    Args:
        V_df: Samples x viruses abundance.
        B_df: Samples x bacteria abundance.
        W_vh_smooth_df: Viruses x bacteria smoothed CRISPR matrix.
        pseudocount: CLR pseudocount.
        eps: Numerical stability constant.

    Returns:
        Tuple (V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub).

    Raises:
        ValueError: If no overlapping samples, viruses, or bacteria are found.
    """
    common_samples = V_df.index.intersection(B_df.index)
    if len(common_samples) == 0:
        raise ValueError("No overlapping sample IDs between V_df and B_df.")

    V0 = V_df.loc[common_samples].copy()
    B0 = B_df.loc[common_samples].copy()

    V_clr_df, B_clr_df, _ = double_clr_transform(V0, B0, pseudocount=pseudocount, eps=eps)

    viruses = V_clr_df.columns.intersection(W_vh_smooth_df.index)
    bacteria = B_clr_df.columns.intersection(W_vh_smooth_df.columns)

    if len(viruses) == 0:
        raise ValueError("No overlapping viruses between V_clr_df.columns and W_vh_smooth_df.index.")
    if len(bacteria) == 0:
        raise ValueError("No overlapping bacteria between B_clr_df.columns and W_vh_smooth_df.columns.")

    V_clr_df = V_clr_df.loc[:, viruses].copy()
    B_clr_df = B_clr_df.loc[:, bacteria].copy()
    X_clr_df = pd.concat([V_clr_df, B_clr_df], axis=1)
    W_sub = W_vh_smooth_df.loc[viruses, bacteria].copy()

    return V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub


# ── Message passing ────────────────────────────────────────────────────────────

def transform_message_passing_smoothed_crispr(
    V: np.ndarray,
    B: np.ndarray,
    W_vh_smooth: np.ndarray,
    lam: float = 0.1,
    n_steps: int = 1,
    preserve_scale: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Non-residual convex cross-domain message passing using a smoothed CRISPR matrix.

    Args:
        V: Samples x viruses matrix (CLR-transformed).
        B: Samples x bacteria matrix (CLR-transformed).
        W_vh_smooth: Viruses x bacteria smoothed CRISPR matrix.
        lam: Mixing weight (0 = no update, 1 = full cross-domain).
        n_steps: Number of propagation steps.
        preserve_scale: If True, restore original column standard deviations.

    Returns:
        Tuple (V_star, B_star).
    """
    V = np.asarray(V, dtype=float)
    B = np.asarray(B, dtype=float)
    W = np.asarray(W_vh_smooth, dtype=float)

    _validate_inputs(V, B, W, lam=lam, n_steps=n_steps)

    P_vh = row_normalize(W)    # viruses -> bacteria, shape q x p
    P_hv = row_normalize(W.T)  # bacteria -> viruses, shape p x q

    Vt = V.copy()
    Bt = B.copy()

    V_sd0 = V.std(axis=0, ddof=0)
    B_sd0 = B.std(axis=0, ddof=0)

    for _ in range(n_steps):
        V_from_hosts = Bt @ P_hv   # n x q
        B_from_virs = Vt @ P_vh    # n x p
        Vt = (1.0 - lam) * Vt + lam * V_from_hosts
        Bt = (1.0 - lam) * Bt + lam * B_from_virs

    if preserve_scale:
        V_sd = Vt.std(axis=0, ddof=0)
        B_sd = Bt.std(axis=0, ddof=0)
        V_scale = np.divide(V_sd0, V_sd, out=np.ones_like(V_sd0), where=V_sd > 1e-12)
        B_scale = np.divide(B_sd0, B_sd, out=np.ones_like(B_sd0), where=B_sd > 1e-12)
        Vt = Vt * V_scale
        Bt = Bt * B_scale

    return Vt, Bt


def transform_message_passing_smoothed_crispr_df(
    V_clr_df: pd.DataFrame,
    B_clr_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    lam: float = 0.1,
    n_steps: int = 1,
    preserve_scale: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """DataFrame wrapper for non-residual convex message passing.

    Args:
        V_clr_df: Samples x viruses CLR-transformed abundance.
        B_clr_df: Samples x bacteria CLR-transformed abundance.
        W_vh_smooth_df: Viruses x bacteria smoothed CRISPR matrix.
        lam: Mixing weight.
        n_steps: Number of propagation steps.
        preserve_scale: If True, restore original column standard deviations.

    Returns:
        Tuple (V_star_df, B_star_df, X_star_df).
    """
    viruses = V_clr_df.columns.intersection(W_vh_smooth_df.index)
    bacteria = B_clr_df.columns.intersection(W_vh_smooth_df.columns)

    if len(viruses) == 0:
        raise ValueError("No overlapping virus names between V_clr_df.columns and W_vh_smooth_df.index.")
    if len(bacteria) == 0:
        raise ValueError("No overlapping bacteria names between B_clr_df.columns and W_vh_smooth_df.columns.")

    V_sub = V_clr_df.loc[:, viruses].copy()
    B_sub = B_clr_df.loc[:, bacteria].copy()
    W_sub = W_vh_smooth_df.loc[viruses, bacteria].copy()

    if not V_sub.index.equals(B_sub.index):
        common_samples = V_sub.index.intersection(B_sub.index)
        if len(common_samples) == 0:
            raise ValueError("No overlapping sample IDs between V_clr_df and B_clr_df.")
        V_sub = V_sub.loc[common_samples]
        B_sub = B_sub.loc[common_samples]

    V_star, B_star = transform_message_passing_smoothed_crispr(
        V=V_sub.to_numpy(),
        B=B_sub.to_numpy(),
        W_vh_smooth=W_sub.to_numpy(),
        lam=lam,
        n_steps=n_steps,
        preserve_scale=preserve_scale,
    )

    V_star_df = pd.DataFrame(V_star, index=V_sub.index, columns=V_sub.columns)
    B_star_df = pd.DataFrame(B_star, index=B_sub.index, columns=B_sub.columns)
    X_star_df = pd.concat([V_star_df, B_star_df], axis=1)
    return V_star_df, B_star_df, X_star_df


def build_xstar_from_smoothed_crispr(
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    pseudocount: float = 1e-6,
    lam: float = 0.1,
    n_steps: int = 1,
    preserve_scale: bool = False,
    eps: float = 1e-12,
) -> dict[str, pd.DataFrame]:
    """Full non-residual X-star pipeline (convex message passing).

    Args:
        V_df: Samples x viruses raw abundance.
        B_df: Samples x bacteria raw abundance.
        W_vh_smooth_df: Viruses x bacteria smoothed CRISPR matrix.
        pseudocount: CLR pseudocount.
        lam: Mixing weight.
        n_steps: Number of propagation steps.
        preserve_scale: If True, restore original column standard deviations.
        eps: Numerical stability constant.

    Returns:
        Dict with keys: V_raw_aligned, B_raw_aligned, V_clr, B_clr, X_clr, V_star, B_star, X_star, W_smooth_aligned.
    """
    V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub = _build_common_inputs(
        V_df, B_df, W_vh_smooth_df, pseudocount=pseudocount, eps=eps
    )
    V_star_df, B_star_df, X_star_df = transform_message_passing_smoothed_crispr_df(
        V_clr_df, B_clr_df, W_sub, lam=lam, n_steps=n_steps, preserve_scale=preserve_scale
    )
    return {
        "V_raw_aligned": V0.loc[:, V_clr_df.columns],
        "B_raw_aligned": B0.loc[:, B_clr_df.columns],
        "V_clr": V_clr_df,
        "B_clr": B_clr_df,
        "X_clr": X_clr_df,
        "V_star": V_star_df,
        "B_star": B_star_df,
        "X_star": X_star_df,
        "W_smooth_aligned": W_sub,
    }


def residual_message_passing(
    V: np.ndarray,
    B: np.ndarray,
    W_vh_smooth: np.ndarray,
    lam: float = 0.1,
    n_steps: int = 1,
    preserve_scale: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Residual additive cross-domain message passing using a smoothed CRISPR matrix.

    Args:
        V: Samples x viruses CLR matrix.
        B: Samples x bacteria CLR matrix.
        W_vh_smooth: Viruses x bacteria smoothed CRISPR matrix.
        lam: Additive weight for cross-domain messages.
        n_steps: Number of propagation steps.
        preserve_scale: If True, restore original column standard deviations.

    Returns:
        Tuple (V_star, B_star).
    """
    V = np.asarray(V, dtype=float)
    B = np.asarray(B, dtype=float)
    W = np.asarray(W_vh_smooth, dtype=float)

    _validate_inputs(V, B, W, lam=lam, n_steps=n_steps)

    P_vh = row_normalize(W)
    P_hv = row_normalize(W.T)

    Vt = V.copy()
    Bt = B.copy()

    V_sd0 = V.std(axis=0, ddof=0)
    B_sd0 = B.std(axis=0, ddof=0)

    for _ in range(n_steps):
        Rv = Vt - Vt.mean(axis=0, keepdims=True)
        Rb = Bt - Bt.mean(axis=0, keepdims=True)
        V_from_hosts = Rb @ P_hv
        B_from_virs = Rv @ P_vh
        Vt = Vt + lam * V_from_hosts
        Bt = Bt + lam * B_from_virs

    if preserve_scale:
        V_sd = Vt.std(axis=0, ddof=0)
        B_sd = Bt.std(axis=0, ddof=0)
        V_scale = np.divide(V_sd0, V_sd, out=np.ones_like(V_sd0), where=V_sd > 1e-12)
        B_scale = np.divide(B_sd0, B_sd, out=np.ones_like(B_sd0), where=B_sd > 1e-12)
        Vt = Vt * V_scale
        Bt = Bt * B_scale

    return Vt, Bt


def residual_message_passing_df(
    V_clr_df: pd.DataFrame,
    B_clr_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    lam: float = 0.1,
    n_steps: int = 1,
    preserve_scale: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """DataFrame wrapper for residual additive message passing.

    Args:
        V_clr_df: Samples x viruses CLR-transformed abundance.
        B_clr_df: Samples x bacteria CLR-transformed abundance.
        W_vh_smooth_df: Viruses x bacteria smoothed CRISPR matrix.
        lam: Additive weight.
        n_steps: Number of propagation steps.
        preserve_scale: If True, restore original column standard deviations.

    Returns:
        Tuple (V_star_df, B_star_df, X_star_df).
    """
    viruses = V_clr_df.columns.intersection(W_vh_smooth_df.index)
    bacteria = B_clr_df.columns.intersection(W_vh_smooth_df.columns)

    if len(viruses) == 0:
        raise ValueError("No overlapping virus names between V_clr_df.columns and W_vh_smooth_df.index.")
    if len(bacteria) == 0:
        raise ValueError("No overlapping bacteria names between B_clr_df.columns and W_vh_smooth_df.columns.")

    V_sub = V_clr_df.loc[:, viruses].copy()
    B_sub = B_clr_df.loc[:, bacteria].copy()
    W_sub = W_vh_smooth_df.loc[viruses, bacteria].copy()

    if not V_sub.index.equals(B_sub.index):
        common_samples = V_sub.index.intersection(B_sub.index)
        if len(common_samples) == 0:
            raise ValueError("No overlapping sample IDs between V_clr_df and B_clr_df.")
        V_sub = V_sub.loc[common_samples]
        B_sub = B_sub.loc[common_samples]

    V_star, B_star = residual_message_passing(
        V=V_sub.to_numpy(),
        B=B_sub.to_numpy(),
        W_vh_smooth=W_sub.to_numpy(),
        lam=lam,
        n_steps=n_steps,
        preserve_scale=preserve_scale,
    )

    V_star_df = pd.DataFrame(V_star, index=V_sub.index, columns=V_sub.columns)
    B_star_df = pd.DataFrame(B_star, index=B_sub.index, columns=B_sub.columns)
    X_star_df = pd.concat([V_star_df, B_star_df], axis=1)
    return V_star_df, B_star_df, X_star_df


def build_xstar_from_smoothed_crispr_residual(
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    W_vh_smooth_df: pd.DataFrame,
    pseudocount: float = 1e-6,
    lam: float = 0.1,
    n_steps: int = 1,
    preserve_scale: bool = False,
    eps: float = 1e-12,
) -> dict[str, pd.DataFrame]:
    """Full residual additive X-star pipeline.

    Args:
        V_df: Samples x viruses raw abundance.
        B_df: Samples x bacteria raw abundance.
        W_vh_smooth_df: Viruses x bacteria smoothed CRISPR matrix.
        pseudocount: CLR pseudocount.
        lam: Additive weight.
        n_steps: Number of propagation steps.
        preserve_scale: If True, restore original column standard deviations.
        eps: Numerical stability constant.

    Returns:
        Dict with keys: V_raw_aligned, B_raw_aligned, V_clr, B_clr, X_clr, V_star, B_star, X_star, W_smooth_aligned.
    """
    V0, B0, V_clr_df, B_clr_df, X_clr_df, W_sub = _build_common_inputs(
        V_df, B_df, W_vh_smooth_df, pseudocount=pseudocount, eps=eps
    )
    V_star_df, B_star_df, X_star_df = residual_message_passing_df(
        V_clr_df, B_clr_df, W_sub, lam=lam, n_steps=n_steps, preserve_scale=preserve_scale
    )
    return {
        "V_raw_aligned": V0.loc[:, V_clr_df.columns],
        "B_raw_aligned": B0.loc[:, B_clr_df.columns],
        "V_clr": V_clr_df,
        "B_clr": B_clr_df,
        "X_clr": X_clr_df,
        "V_star": V_star_df,
        "B_star": B_star_df,
        "X_star": X_star_df,
        "W_smooth_aligned": W_sub,
    }


# ── CRISPR smoothing ───────────────────────────────────────────────────────────

def smooth_crispr_bac_vir(
    crispr_df: pd.DataFrame,
    K_bac: pd.DataFrame,
    K_vir: pd.DataFrame,
    alpha: float = 1.0,
    preserve_original: bool = True,
) -> pd.DataFrame:
    """Smooth a CRISPR matrix using taxonomy kernels.

    W_smooth = (1 - alpha) * W + alpha * K_bac @ W @ K_vir

    Args:
        crispr_df: Bacteria x viruses binary CRISPR matrix.
        K_bac: Bacteria taxonomy kernel.
        K_vir: Virus taxonomy kernel.
        alpha: Smoothing weight (1.0 = full kernel propagation).
        preserve_original: If True, blend with original; if False, use propagated only.

    Returns:
        Smoothed CRISPR DataFrame.
    """
    bac_common = crispr_df.index.intersection(K_bac.index)
    vir_common = crispr_df.columns.intersection(K_vir.index)

    if len(bac_common) == 0:
        raise ValueError("No overlapping bacteria between crispr_df.index and K_bac.index")
    if len(vir_common) == 0:
        raise ValueError("No overlapping viruses between crispr_df.columns and K_vir.index")

    W = crispr_df.loc[bac_common, vir_common].copy()
    Kb = K_bac.loc[bac_common, bac_common]
    Kv = K_vir.loc[vir_common, vir_common]

    W_prop = Kb.to_numpy() @ W.to_numpy(dtype=float) @ Kv.to_numpy()

    if preserve_original:
        W_smooth = (1 - alpha) * W.to_numpy(dtype=float) + alpha * W_prop
    else:
        W_smooth = W_prop

    return pd.DataFrame(W_smooth, index=W.index, columns=W.columns)


def build_taxonomy_kernel_from_shared_ranks(
    ids,
    tax_df: pd.DataFrame,
    ranks,
    weights=None,
    fill_value: str = "",
    normalize_rows: bool = True,
    strict: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a taxonomy kernel matrix using weighted shared rank agreement.

    K[i,j] = sum_k weights[k] * I(rank_k[i] == rank_k[j]) if rank_k[i] != fill_value.
    K[i,i] = 1.

    Args:
        ids: Ordered collection of feature IDs to include.
        tax_df: Taxonomy DataFrame indexed by the same IDs.
        ranks: Ordered list of taxonomy rank column names.
        weights: Per-rank weights (default: 1..n). Normalized to sum to 1.
        fill_value: Value treated as missing (no match score).
        normalize_rows: Row-normalize the kernel.
        strict: Raise ValueError if any IDs are missing from tax_df.

    Returns:
        Tuple (K DataFrame, aligned taxonomy DataFrame).
    """
    ids = pd.Index(ids)
    ranks = list(ranks)

    missing_ranks = [r for r in ranks if r not in tax_df.columns]
    if missing_ranks:
        raise ValueError(f"Missing taxonomy columns: {missing_ranks}")

    common = ids.intersection(tax_df.index)
    missing = ids.difference(tax_df.index)

    if len(common) == 0:
        raise ValueError("No overlap between ids and taxonomy index")

    if len(missing) > 0:
        msg = f"Missing taxonomy for {len(missing)} taxa"
        if strict:
            raise ValueError(msg)
        else:
            print(msg + " — keeping only overlapping taxa.")

    if weights is None:
        weights = np.arange(1, len(ranks) + 1, dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)

    if len(weights) != len(ranks):
        raise ValueError("weights must have same length as ranks")

    weights = weights / weights.sum()

    tax = tax_df.loc[common, ranks].copy()
    tax = tax.fillna(fill_value).astype(str)
    tax = tax.replace(
        {"nan": fill_value, "None": fill_value, "NA": fill_value,
         "N/A": fill_value, "unknown": fill_value, "unclassified": fill_value}
    )
    for r in ranks:
        tax[r] = tax[r].str.replace(r"^[a-z]__+", "", regex=True)

    X = tax.to_numpy()
    n = X.shape[0]
    K = np.zeros((n, n), dtype=float)

    for k in range(len(ranks)):
        col = X[:, k]
        same = (col[:, None] == col[None, :]) & (col[:, None] != fill_value)
        K += weights[k] * same.astype(float)

    np.fill_diagonal(K, 1.0)

    if normalize_rows:
        rs = K.sum(axis=1, keepdims=True)
        rs[rs == 0] = 1.0
        K = K / rs

    K = pd.DataFrame(K, index=common, columns=common)
    return K, tax


def assign_crispr(cri_big: pd.DataFrame, cri_s: pd.DataFrame) -> pd.DataFrame:
    """Fill matching rows/columns of a target matrix from a source CRISPR matrix.

    Args:
        cri_big: Target (full-size) matrix.
        cri_s: Source CRISPR matrix (subset).

    Returns:
        Updated copy of cri_big.
    """
    full_cri = cri_big.copy()
    overlap_rows = full_cri.index.intersection(cri_s.index)
    overlap_cols = full_cri.columns.intersection(cri_s.columns)
    if len(overlap_rows) > 0 and len(overlap_cols) > 0:
        full_cri.loc[overlap_rows, overlap_cols] = cri_s.loc[overlap_rows, overlap_cols]
    return full_cri


def build_binary_crispr_matrix(
    raw_crispr_path: str,
    bacteria_features,
    virus_features,
    transpose_after_load: bool = True,
) -> pd.DataFrame:
    """Load a raw CRISPR network and build a binary bacteria x viruses matrix.

    Args:
        raw_crispr_path: Path to the CSV of the raw CRISPR network.
        bacteria_features: Ordered list of bacteria feature IDs.
        virus_features: Ordered list of virus feature IDs.
        transpose_after_load: If True, transpose the loaded matrix (contigs were rows).

    Returns:
        Binary bacteria x viruses DataFrame.
    """
    df_crispr = pd.read_csv(raw_crispr_path, index_col=0)
    if transpose_after_load:
        df_crispr = df_crispr.T

    crispr = df_crispr.T
    crispr.index = [str(i) for i in crispr.index]

    target = pd.DataFrame(
        np.zeros((len(bacteria_features), len(virus_features))),
        index=pd.Index(bacteria_features),
        columns=pd.Index(virus_features),
    )
    full = assign_crispr(target, crispr)
    full[full != 0] = 1
    return full


def get_hierarchies(df_b1: pd.DataFrame, df_v1: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Prepare bacteria and virus taxonomy hierarchies for CRISPR smoothing.

    Args:
        df_b1: Bacteria taxonomy DataFrame with progenomes_taxid_genus column.
        df_v1: Virus taxonomy DataFrame with lev0 column.

    Returns:
        Tuple (bacteria_tax, virus_tax) with cleaned indexes.
    """
    df_b = df_b1.copy()
    df_v = df_v1.copy()

    df_v = df_v.set_index("lev0")

    silva_str = [str(i) for i in df_b["progenomes_taxid_genus"]]
    silva_str1 = [i.split(".")[0] for i in silva_str]
    df_b["progenomes_taxid_genus"] = silva_str1
    df_b = df_b.set_index("progenomes_taxid_genus")

    keep_index = [i for i in df_b.index if "nan" not in i]
    df_b = df_b.loc[keep_index]
    df_b = df_b.drop_duplicates()
    df_b = df_b.loc[~df_b.index.duplicated(keep="first")]
    df_v = df_v.loc[~df_v.index.duplicated(keep="first")]
    df_b["progenomes_taxid_genus"] = df_b.index
    df_v["lev0"] = df_v.index

    return df_b, df_v


def build_smoothed_crispr_for_study(
    raw_crispr_path: str,
    bacteria_features,
    virus_features,
    tax_bac_path: str,
    tax_vir_path: str,
    bacterial_ranks,
    viral_ranks,
    bacterial_weights,
    viral_weights,
    alpha: float = 0.95,
    transpose_after_load: bool = True,
) -> dict[str, pd.DataFrame]:
    """Build a smoothed CRISPR matrix for a single study.

    Args:
        raw_crispr_path: Path to the raw CRISPR network CSV.
        bacteria_features: Bacteria feature IDs (from processed abundance).
        virus_features: Virus feature IDs (from processed abundance).
        tax_bac_path: Path to bacteria taxonomy CSV.
        tax_vir_path: Path to virus taxonomy CSV.
        bacterial_ranks: Bacterial taxonomy rank columns for kernel.
        viral_ranks: Viral taxonomy rank columns for kernel.
        bacterial_weights: Per-rank weights for bacteria kernel.
        viral_weights: Per-rank weights for virus kernel.
        alpha: CRISPR smoothing weight.
        transpose_after_load: Passed to build_binary_crispr_matrix.

    Returns:
        Dict with crispr_binary, K_bac, K_vir, crispr_smooth, bac_tax_aligned, vir_tax_aligned.
    """
    tax_bac = pd.read_csv(tax_bac_path, index_col=0)
    tax_vir = pd.read_csv(tax_vir_path, index_col=0)
    tax_bac_aligned, tax_vir_aligned = get_hierarchies(tax_bac, tax_vir)

    crispr_binary = build_binary_crispr_matrix(
        raw_crispr_path=raw_crispr_path,
        bacteria_features=bacteria_features,
        virus_features=virus_features,
        transpose_after_load=transpose_after_load,
    )

    K_bac, bac_tax_aligned = build_taxonomy_kernel_from_shared_ranks(
        ids=crispr_binary.index,
        tax_df=tax_bac_aligned,
        ranks=bacterial_ranks,
        weights=bacterial_weights,
        normalize_rows=True,
        strict=False,
    )

    K_vir, vir_tax_aligned = build_taxonomy_kernel_from_shared_ranks(
        ids=crispr_binary.columns,
        tax_df=tax_vir_aligned,
        ranks=viral_ranks,
        weights=viral_weights,
        normalize_rows=True,
        strict=False,
    )

    crispr_smooth = smooth_crispr_bac_vir(
        crispr_df=crispr_binary,
        K_bac=K_bac,
        K_vir=K_vir,
        alpha=alpha,
        preserve_original=True,
    )

    return {
        "crispr_binary": crispr_binary,
        "K_bac": K_bac,
        "K_vir": K_vir,
        "crispr_smooth": crispr_smooth,
        "bac_tax_aligned": bac_tax_aligned,
        "vir_tax_aligned": vir_tax_aligned,
    }


def crispr_matrix_aggregate_viruses(
    df_crispr: pd.DataFrame,
    vir_tax: pd.DataFrame,
    *,
    bac_col: int = 0,
    vir_col: int = 1,
    vir_rank: str = "lev0",
    vir_id_col=None,
    dropna_vir: bool = True,
    dtype=int,
) -> pd.DataFrame:
    """Parse a SpacePHARER output TSV and aggregate viruses by taxonomy rank.

    Args:
        df_crispr: Raw SpacePHARER predictions DataFrame (no header).
        vir_tax: Virus taxonomy DataFrame.
        bac_col: Column index for bacterial spacer IDs.
        vir_col: Column index for viral contig IDs.
        vir_rank: Viral taxonomy rank column to aggregate to.
        vir_id_col: Optional viral ID column in vir_tax; uses index if None.
        dropna_vir: Drop viruses not found in taxonomy.
        dtype: Output matrix dtype.

    Returns:
        Bacteria x viral_groups crosstab matrix.
    """
    def parse_bac_taxid(s):
        if pd.isna(s):
            return None
        s = str(s)
        try:
            i1 = s.split(">")[1]
            return int(i1.split(".")[0])
        except Exception:
            return None

    df1 = df_crispr[[bac_col, vir_col]].copy()
    df1["bac_taxid"] = df1[bac_col].map(parse_bac_taxid)
    df1["contig"] = df1[vir_col].astype(str).str.strip()
    df1 = df1.dropna(subset=["bac_taxid", "contig"])
    df1["bac_taxid"] = df1["bac_taxid"].astype(int)

    M = pd.crosstab(df1["bac_taxid"], df1["contig"]).astype(dtype)

    if vir_rank not in vir_tax.columns:
        raise KeyError(f"vir_tax missing column {vir_rank!r}. Available: {list(vir_tax.columns)[:20]}")

    if vir_id_col is None:
        vt = vir_tax[[vir_rank]].copy()
        vt.index = vt.index.astype(str).str.strip()
        vt = vt.loc[~vt.index.duplicated(keep="first")]
        contig_to_group = vt[vir_rank]
        labels = pd.Series(M.columns.astype(str), index=M.columns).map(contig_to_group)
    else:
        if vir_id_col not in vir_tax.columns:
            raise KeyError(f"vir_tax missing id column {vir_id_col!r}")
        vt = vir_tax[[vir_id_col, vir_rank]].copy()
        vt = vt.dropna(subset=[vir_id_col])
        vt[vir_id_col] = vt[vir_id_col].astype(str).str.strip()
        vt = vt.drop_duplicates(subset=[vir_id_col], keep="first")
        contig_to_group = vt.set_index(vir_id_col)[vir_rank]
        labels = pd.Series(M.columns.astype(str), index=M.columns).map(contig_to_group)

    valid = labels.notna() & (labels.astype(str).str.strip() != "")
    if dropna_vir:
        M = M.loc[:, valid.values]
        labels = labels.loc[valid]
    else:
        labels = labels.fillna("Unassigned")

    out = M.T.groupby(labels.astype(str).str.strip()).sum().T
    out.index = out.index.astype(str)
    return out


# ── Abundance helpers ──────────────────────────────────────────────────────────

def aggregate_otu_columns_by_rank_skip_nan(
    otu: pd.DataFrame,
    tax: pd.DataFrame,
    rank: str,
) -> pd.DataFrame:
    """Aggregate OTU/ASV columns to a taxonomy rank, skipping NaN labels.

    Args:
        otu: Samples x ASVs abundance DataFrame.
        tax: ASVs x ranks taxonomy DataFrame.
        rank: Taxonomy column to aggregate to.

    Returns:
        Samples x rank-groups abundance DataFrame.
    """
    common_asvs = otu.columns.intersection(tax.index)
    if len(common_asvs) == 0:
        raise ValueError("No common ASV/OTU IDs between abundance columns and taxonomy index.")

    otu2 = otu.loc[:, common_asvs]
    labels = tax.loc[common_asvs, rank]

    keep = labels.notna() & (labels.astype(str).str.strip() != "")
    otu2 = otu2.loc[:, keep.values]
    labels = labels.loc[keep].astype(str).str.strip()

    return otu2.T.groupby(labels).sum().T


def prepare_bacteria_genus_abundance(
    otu: pd.DataFrame,
    tax: pd.DataFrame,
    rank: str = "target_taxids",
) -> pd.DataFrame:
    """Aggregate bacteria OTUs to the selected taxonomy rank and sanitize feature names.

    Args:
        otu: Samples x ASVs raw abundance.
        tax: ASV taxonomy DataFrame with a column matching rank.
        rank: Taxonomy column to aggregate to.

    Returns:
        Samples x genus-level abundance DataFrame with sanitized column names.
    """
    genus_abund = aggregate_otu_columns_by_rank_skip_nan(otu, tax, rank)
    genus_abund = rename_clostridium_sensu_stricto(genus_abund)
    genus_abund.columns = sanitize_index(genus_abund.columns)
    genus_abund = apply_custom_renames(genus_abund)
    return genus_abund


def remove_disease_columns_from_virus_abundance(V: pd.DataFrame) -> pd.DataFrame:
    """Remove phenotype/metadata columns accidentally stored in viral abundance tables.

    Args:
        V: Viral abundance DataFrame.

    Returns:
        Cleaned copy without non-feature columns.
    """
    V = V.copy()
    bad_cols = {
        "disease", "Disease", "disease_original", "disease_binary",
        "phenotype", "Phenotype", "label", "Label",
        "class", "Class", "group", "Group",
        "sample_id", "SampleID", "subject_id", "SubjectNo",
        "Reads", "reads",
    }
    cols_to_drop = [c for c in V.columns if str(c) in bad_cols]
    if cols_to_drop:
        print("Dropping non-feature columns from viral abundance:", cols_to_drop)
        V = V.drop(columns=cols_to_drop)
    return V


def prevalence_filter_df(df: pd.DataFrame, prevalence: float = 0.10, verbose: bool = True) -> pd.DataFrame:
    """Keep only features present in at least prevalence * n_samples samples.

    Args:
        df: Samples x features DataFrame.
        prevalence: Minimum fractional prevalence threshold.
        verbose: Print kept/total feature count.

    Returns:
        Filtered DataFrame.
    """
    n_samples = df.shape[0]
    min_samples = max(int(prevalence * n_samples), 1)
    keep = (df > 0).sum(axis=0) >= min_samples
    if verbose:
        print(
            f"prevalence={prevalence} -> min_samples={min_samples}/{n_samples}; "
            f"kept {int(keep.sum())}/{df.shape[1]} features"
        )
    return df.loc[:, keep]


def align_abundance_from_metadata(
    virus_abundance: pd.DataFrame,
    bacteria_abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    keep_col: str = "keep_for_analysis",
    virus_id_col: str = "virus_sample_id",
    bacteria_id_col: str = "bacteria_sample_id",
    final_index_col: str = "virus_sample_id",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Align viral and bacterial abundance using standardized metadata.

    Args:
        virus_abundance: Viral abundance DataFrame.
        bacteria_abundance: Bacterial abundance DataFrame.
        metadata: Metadata DataFrame with sample ID and keep columns.
        keep_col: Column used to filter metadata rows.
        virus_id_col: Metadata column for viral sample IDs.
        bacteria_id_col: Metadata column for bacterial sample IDs.
        final_index_col: Metadata column to use as the aligned index.

    Returns:
        Tuple (V, B, meta_aligned) with aligned, filtered abundance and metadata.
    """
    V_raw = clean_df_ids(virus_abundance)
    B_raw = clean_df_ids(bacteria_abundance)
    meta = metadata.copy()

    for col in [virus_id_col, bacteria_id_col]:
        if col not in meta.columns:
            raise ValueError(f"metadata is missing required column: {col}")

    if keep_col is not None and keep_col in meta.columns:
        keep = parse_bool_series(meta[keep_col])
        print(f"metadata filter {keep_col}: kept {int(keep.sum())} / {len(keep)} rows")
        meta = meta.loc[keep].copy()
    else:
        print(f"metadata filter {keep_col!r} not found; keeping all {meta.shape[0]} rows")

    if "analysis_order" in meta.columns:
        meta = meta.sort_values("analysis_order").copy()

    for col in [virus_id_col, bacteria_id_col, final_index_col, "sample_id", "subject_id"]:
        if col in meta.columns:
            meta[col] = clean_index_ids(meta[col])

    virus_ids = list(meta[virus_id_col].astype(str))
    bacteria_ids = list(meta[bacteria_id_col].astype(str))

    missing_v = sorted(set(virus_ids) - set(V_raw.index.astype(str)))
    missing_b = sorted(set(bacteria_ids) - set(B_raw.index.astype(str)))

    if missing_v:
        raise ValueError(f"{len(missing_v)} viral sample IDs missing in viral abundance. Examples: {missing_v[:10]}")
    if missing_b:
        raise ValueError(f"{len(missing_b)} bacterial sample IDs missing in bacterial abundance. Examples: {missing_b[:10]}")

    V = V_raw.loc[virus_ids].copy()
    B = B_raw.loc[bacteria_ids].copy()

    final_index = list(meta[final_index_col].astype(str)) if final_index_col in meta.columns else virus_ids
    V.index = final_index
    B.index = final_index

    meta = meta.reset_index(drop=True)
    return V, B, meta


def summarize_df(name: str, df: pd.DataFrame) -> None:
    """Print a short summary of a DataFrame shape and index uniqueness.

    Args:
        name: Label for the printout.
        df: DataFrame to summarize.
    """
    print(f"{name}: shape={df.shape}, index_unique={df.index.is_unique}, columns_unique={df.columns.is_unique}")


def study_outdir(study: str, subdir: str, output_root: Path) -> Path:
    """Construct the output directory path for a study.

    Args:
        study: Study identifier string.
        subdir: Sub-directory name within the study folder.
        output_root: Root output directory.

    Returns:
        Path to the study sub-directory.
    """
    base = output_root / study
    return base / subdir if subdir else base


def out_path(study: str, subdir: str, filename: str, output_root: Path) -> Path:
    """Construct the full path to an output file within a study directory.

    Args:
        study: Study identifier string.
        subdir: Sub-directory name.
        filename: File name.
        output_root: Root output directory.

    Returns:
        Full Path to the output file.
    """
    return study_outdir(study, subdir, output_root) / filename


def orient_W_viruses_by_bacteria(
    W_df: pd.DataFrame,
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    verbose: bool = True,
) -> pd.DataFrame:
    """Detect and correct the orientation of a CRISPR matrix to be viruses x bacteria.

    Args:
        W_df: CRISPR matrix (may be bacteria x viruses or viruses x bacteria).
        V_df: Viral abundance (samples x viruses).
        B_df: Bacterial abundance (samples x bacteria).
        verbose: Print overlap counts.

    Returns:
        W_df in viruses x bacteria orientation.
    """
    n_vir_row = len(set(W_df.index.astype(str)) & set(map(str, V_df.columns)))
    n_bac_col = len(set(W_df.columns.astype(str)) & set(map(str, B_df.columns)))
    n_bac_row = len(set(W_df.index.astype(str)) & set(map(str, B_df.columns)))
    n_vir_col = len(set(W_df.columns.astype(str)) & set(map(str, V_df.columns)))

    if verbose:
        print("W rows overlapping viruses:", n_vir_row, "/", V_df.shape[1])
        print("W cols overlapping bacteria:", n_bac_col, "/", B_df.shape[1])
        print("W rows overlapping bacteria:", n_bac_row, "/", B_df.shape[1])
        print("W cols overlapping viruses:", n_vir_col, "/", V_df.shape[1])

    if n_bac_row > n_vir_row and n_vir_col > n_bac_col:
        print("Detected bacteria x viruses orientation; transposing W.")
        W_df = W_df.T

    return W_df
