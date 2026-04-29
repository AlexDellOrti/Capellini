"""I/O helpers: file reading, writing, subprocess execution."""

from __future__ import annotations

import gzip
import subprocess
from pathlib import Path
from typing import IO


def open_maybe_gzip(path: str | Path, mode: str = "rt", encoding: str = "utf-8") -> IO:
    """Open a file, transparently decompressing if it ends with .gz.

    Args:
        path: Path to the file.
        mode: File open mode.
        encoding: Text encoding.

    Returns:
        File-like object.
    """
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, mode, encoding=encoding, errors="replace")
    return open(p, mode, encoding=encoding, errors="replace")


def read_table(path: str | Path, index_col: int = 0, **kwargs):
    """Read a CSV, TSV, or Excel file into a DataFrame, including gzip-compressed variants.

    Args:
        path: Path to the tabular file.
        index_col: Column to use as the row index.
        **kwargs: Additional keyword arguments forwarded to pandas.

    Returns:
        pd.DataFrame with the file contents.
    """
    import pandas as pd

    path = Path(path)
    suffixes = "".join(path.suffixes).lower()

    if suffixes.endswith(".xlsx") or suffixes.endswith(".xls"):
        return pd.read_excel(path, index_col=index_col, **kwargs)
    if suffixes.endswith(".tsv") or suffixes.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t", index_col=index_col, **kwargs)
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz") or suffixes.endswith(".gz"):
        return pd.read_csv(path, index_col=index_col, **kwargs)
    raise ValueError(f"Unsupported file type: {path}")


def write_df(df, path: str | Path, *, overwrite: bool = True, verbose: bool = False, **to_csv_kwargs) -> Path:
    """Write a DataFrame to CSV, creating parent directories as needed.

    Args:
        df: DataFrame to write.
        path: Destination file path.
        overwrite: Skip if file exists and overwrite is False.
        verbose: Print a message on skip or save.
        **to_csv_kwargs: Forwarded to DataFrame.to_csv.

    Returns:
        Path to the written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        if verbose:
            print("skip existing:", path)
        return path
    df.to_csv(path, **to_csv_kwargs)
    if verbose:
        print("saved:", path, getattr(df, "shape", ""))
    return path


def sh(cmd: str, desc: str = "") -> subprocess.CompletedProcess:
    """Run a shell command, printing description and command, raising on failure.

    Args:
        cmd: Shell command string.
        desc: Human-readable description printed before execution.

    Returns:
        CompletedProcess result.

    Raises:
        RuntimeError: If the command exits with a non-zero return code,
            including full stdout and stderr in the message.
    """
    if desc:
        print(desc)
    print(f"Executing command: {cmd}")
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        msg = (
            f"Command failed with code {e.returncode}\n"
            f"STDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}"
        )
        raise RuntimeError(msg) from e
    if r.stdout:
        print(r.stdout)
    if r.stderr.strip():
        print("STDERR (tool messages):\n" + r.stderr)
    return r
