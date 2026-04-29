"""Numerical transformations: CLR, closure, message-passing, shrinkage."""

from __future__ import annotations

import numpy as np
import pandas as pd


def closure(X: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Row-wise closure: each row sums to 1.

    Args:
        X: Nonnegative matrix of shape (n_samples, n_features).
        eps: Small constant to avoid division by zero.

    Returns:
        Closed matrix with row sums approximately equal to 1.
    """
    X = np.asarray(X, dtype=float)
    row_sums = X.sum(axis=1, keepdims=True)
    return X / (row_sums + eps)


def clr(X: np.ndarray, pseudocount: float = 1e-6, eps: float = 1e-12) -> np.ndarray:
    """Row-wise centered log-ratio transform.

    Args:
        X: Nonnegative matrix.
        pseudocount: Added after clipping to handle zeros safely.
        eps: Numerical stability constant inside the log.

    Returns:
        CLR-transformed matrix.
    """
    X = np.asarray(X, dtype=float)
    X = np.clip(X, 0.0, None) + pseudocount
    log_x = np.log(X + eps)
    return log_x - log_x.mean(axis=1, keepdims=True)


def double_clr_transform(
    V_df: pd.DataFrame,
    B_df: pd.DataFrame,
    pseudocount: float = 1e-6,
    eps: float = 1e-12,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Apply closure and CLR separately to viruses and bacteria.

    Args:
        V_df: Samples x viruses abundance matrix.
        B_df: Samples x bacteria abundance matrix.
        pseudocount: Pseudocount for CLR.
        eps: Numerical stability constant.

    Returns:
        Tuple of (V_clr_df, B_clr_df, X_clr_df).
    """
    V_closed = closure(V_df.to_numpy(), eps=eps)
    B_closed = closure(B_df.to_numpy(), eps=eps)

    V_clr = clr(V_closed, pseudocount=pseudocount, eps=eps)
    B_clr = clr(B_closed, pseudocount=pseudocount, eps=eps)

    V_clr_df = pd.DataFrame(V_clr, index=V_df.index, columns=V_df.columns)
    B_clr_df = pd.DataFrame(B_clr, index=B_df.index, columns=B_df.columns)
    X_clr_df = pd.concat([V_clr_df, B_clr_df], axis=1)

    return V_clr_df, B_clr_df, X_clr_df


def row_normalize(W: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Row-normalize a matrix so each row sums to 1.

    Args:
        W: Matrix to normalize.
        eps: Small constant to avoid division by zero.

    Returns:
        Row-normalized matrix.
    """
    W = np.asarray(W, dtype=float)
    row_sums = W.sum(axis=1, keepdims=True)
    return W / (row_sums + eps)


def normalize_columns(X: pd.DataFrame) -> pd.DataFrame:
    """Normalize each column of a DataFrame to sum to 1.

    Args:
        X: Input DataFrame (features x samples convention from old notebook).

    Returns:
        Column-normalized DataFrame.
    """
    col_sums = X.sum(axis=0)
    return X.divide(col_sums.replace(0, np.nan), axis=1).fillna(0.0)


def geometric_mean_safe(x: pd.Series, eps: float = 1e-12) -> float:
    """Compute the geometric mean of a series, clipping zeros.

    Args:
        x: Input series.
        eps: Floor value for clipping.

    Returns:
        Geometric mean as a float.
    """
    x = np.asarray(x, dtype=float)
    x = np.clip(x, eps, None)
    return float(np.exp(np.log(x).mean()))


def old_style_clr_transform(df: pd.DataFrame) -> pd.DataFrame:
    """CLR transform matching the original notebook's shrinkage preprocessing.

    Input is samples x features. The transform adds 1, transposes to features x samples,
    normalizes samples, divides by sample geometric mean, takes log, and transposes back.

    Args:
        df: Samples x features DataFrame.

    Returns:
        CLR-transformed DataFrame (samples x features).
    """
    X = (df + 1.0).T
    X = normalize_columns(X)
    g = X.apply(geometric_mean_safe, axis=0)
    Z = np.log(X.divide(g.replace(0, np.nan), axis=1)).replace([np.inf, -np.inf], 0).fillna(0)
    return Z.T


def schaefer_strimmer_corr(X: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Schäfer-Strimmer shrinkage correlation estimator.

    Args:
        X: Samples x features DataFrame (at least 3 samples required).

    Returns:
        Tuple of (shrinkage_correlation_DataFrame, diagnostics_dict).
        diagnostics_dict contains keys: n_samples, n_features, lambda_var, lambda_corr.

    Raises:
        ValueError: If fewer than 3 samples are provided.
    """
    values = X.to_numpy(dtype=float)
    n, p = values.shape
    if n < 3:
        raise ValueError(f"Need at least 3 samples for shrinkage correlation, got {n}")

    w = (values - values.mean(axis=0)) ** 2
    w_bar = np.mean(w, axis=0)
    var_unb = (n / (n - 1)) * w_bar
    var_s = (n / (n - 1) ** 3) * np.sum((w - w_bar) ** 2, axis=0)

    med_var = np.median(var_unb)
    denom_var = np.sum((var_unb - med_var) ** 2)
    lambda_var = 1.0 if denom_var <= 1e-15 else min(1.0, float(np.sum(var_s) / denom_var))

    sd = np.std(values, axis=0, ddof=1)
    X_st = np.divide(values, sd, out=np.zeros_like(values), where=sd > 1e-12)
    X_c_st = X_st - X_st.mean(axis=0)

    w_st = X_c_st.T @ X_c_st
    w_st_sq = (X_c_st ** 2).T @ (X_c_st ** 2)
    w_bar_st = w_st / n
    var_s_st = (n / (n - 1) ** 3) * (w_st_sq - 2 * w_bar_st * w_st + n * w_bar_st ** 2)

    corr_unb_st = (n / (n - 1)) * w_bar_st
    denom_corr = np.sum(corr_unb_st ** 2) - np.sum(np.diag(corr_unb_st) ** 2)
    numer_corr = np.sum(var_s_st) - np.sum(np.diag(var_s_st))
    lambda_corr = 1.0 if denom_corr <= 1e-15 else min(1.0, float(numer_corr / denom_corr))

    corr_X = np.nan_to_num(np.corrcoef(values.T), nan=0.0, posinf=0.0, neginf=0.0)
    corr_shrink = (1.0 - lambda_corr) * corr_X
    np.fill_diagonal(corr_shrink, 1.0)

    out = pd.DataFrame(corr_shrink, index=X.columns, columns=X.columns)
    diag = {
        "n_samples": int(n),
        "n_features": int(p),
        "lambda_var": float(lambda_var),
        "lambda_corr": float(lambda_corr),
    }
    return out, diag
