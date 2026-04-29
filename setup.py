"""Build-time guard: refuse to install capellini if external CLI tools are missing.

All Python metadata still comes from pyproject.toml. This file exists only to
register a custom ``build`` command that runs before any wheel/editable install
and aborts with a clear error message listing what to install and how.

Bypass with ``CAPELLINI_SKIP_DEP_CHECK=1`` (intended for CI / source archives
where the binaries aren't yet on PATH but will be at runtime).
"""

from __future__ import annotations

import os
import shutil
import sys

from setuptools import setup

# (display name, candidate executables, install hint)
REQUIRED_TOOLS: list[tuple[str, list[str], str]] = [
    ("spacepharer", ["spacepharer"], "conda install -c bioconda spacepharer"),
    ("mmseqs2",     ["mmseqs", "mmseqs2"], "conda install -c bioconda mmseqs2"),
    ("minced",      ["minced"], "conda install -c bioconda minced"),
    ("prodigal",    ["prodigal"], "conda install -c bioconda prodigal"),
    ("Rscript (with the dada2 Bioconductor package)",
        ["Rscript"],
        "brew install r   # then in R:  BiocManager::install('dada2')"),
    ("micro",       ["micro"], "brew install micro"),
]


def _missing_tools() -> list[tuple[str, str]]:
    missing: list[tuple[str, str]] = []
    for name, candidates, hint in REQUIRED_TOOLS:
        if not any(shutil.which(c) for c in candidates):
            missing.append((name, hint))
    return missing


def _enforce_external_deps() -> None:
    if os.environ.get("CAPELLINI_SKIP_DEP_CHECK") == "1":
        return
    missing = _missing_tools()
    if not missing:
        return
    print(
        "\n"
        "============================================================\n"
        "  capellini cannot be installed: missing external tools.\n"
        "  Install the following and re-run `pip install`:\n"
        "============================================================",
        file=sys.stderr,
    )
    for name, hint in missing:
        print(f"\n  • {name}\n      install: {hint}", file=sys.stderr)
    print(
        "\n"
        "  (set CAPELLINI_SKIP_DEP_CHECK=1 to bypass this check —\n"
        "   recommended only for CI / containers where the tools\n"
        "   will be present at runtime.)\n",
        file=sys.stderr,
    )
    sys.exit(1)


_enforce_external_deps()
setup()
