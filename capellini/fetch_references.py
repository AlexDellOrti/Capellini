"""Download the bundled reference FASTAs from the GitHub release.

These files are too large to ship inside the source repo, so we host them as
release assets on GitHub and download them post-install on demand.

Override the source tag with the ``CAPELLINI_REFERENCES_TAG`` env var if you
need to pin a different release.
"""

from __future__ import annotations

import importlib.resources as pkg_resources
import logging
import os
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Iterable
from urllib.error import HTTPError, URLError

logger = logging.getLogger(__name__)

DEFAULT_TAG = "v0.1.0"
RELEASES_BASE = "https://github.com/AlexDellOrti/Capellini/releases/download"

# (asset_filename, destination relative to capellini/data/)
ASSETS: tuple[tuple[str, str], ...] = (
    ("progenome16S.fasta", "references/progenome16S.fasta"),
    ("spacers_CompleteCollection.fasta", "references/spacers/spacers_CompleteCollection.fasta"),
)


def _data_root() -> Path:
    """Return the absolute path to the installed ``capellini/data`` directory."""
    ref = pkg_resources.files("capellini.data").joinpath(".")
    return Path(str(ref))


def _release_tag() -> str:
    return os.environ.get("CAPELLINI_REFERENCES_TAG", DEFAULT_TAG).strip() or DEFAULT_TAG


def _asset_url(tag: str, filename: str) -> str:
    return f"{RELEASES_BASE}/{tag}/{filename}"


def _human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def _download_one(url: str, dest: Path) -> None:
    """Stream-download ``url`` to ``dest`` with a simple progress bar."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")

    def _progress(block_num: int, block_size: int, total: int) -> None:
        if total <= 0:
            return
        downloaded = block_num * block_size
        pct = min(100.0, downloaded * 100.0 / total)
        bar = "█" * int(pct // 2) + "·" * (50 - int(pct // 2))
        sys.stdout.write(
            f"\r  [{bar}] {pct:5.1f}%  {_human_size(downloaded)} / {_human_size(total)}"
        )
        sys.stdout.flush()

    print(f"  → {url}")
    try:
        urllib.request.urlretrieve(url, str(tmp), reporthook=_progress)
    except (HTTPError, URLError) as exc:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        raise RuntimeError(f"Download failed for {url}: {exc}") from exc
    print()  # newline after progress bar
    shutil.move(str(tmp), str(dest))


def fetch_references(
    tag: str | None = None,
    overwrite: bool = False,
    assets: Iterable[tuple[str, str]] = ASSETS,
) -> list[Path]:
    """Download the bundled reference FASTAs from the GitHub release.

    Args:
        tag: Release tag to pull from. Defaults to the env var
            ``CAPELLINI_REFERENCES_TAG`` if set, otherwise ``DEFAULT_TAG``.
        overwrite: If True, re-download even when the file already exists.
        assets: Iterable of ``(filename, dest_relative_to_data_root)`` pairs.

    Returns:
        List of paths to the downloaded (or already-present) files.
    """
    tag = tag or _release_tag()
    data_root = _data_root()

    print(f"CAPELLINI references — release tag: {tag}")
    print(f"Target directory:  {data_root}\n")

    paths: list[Path] = []
    for filename, rel in assets:
        dest = data_root / rel
        if dest.exists() and not overwrite:
            print(f"✓ {rel}  (already present, skipping)")
            paths.append(dest)
            continue
        print(f"↓ {rel}")
        _download_one(_asset_url(tag, filename), dest)
        size = dest.stat().st_size
        print(f"  saved: {dest}  ({_human_size(size)})\n")
        paths.append(dest)

    print("Done.")
    return paths


def main(argv: list[str] | None = None) -> int:
    """Console-script entry point: ``capellini-fetch-references``."""
    import argparse

    parser = argparse.ArgumentParser(
        prog="capellini-fetch-references",
        description=(
            "Download the bundled CAPELLINI reference FASTAs "
            "(progenome16S.fasta, spacers_CompleteCollection.fasta) "
            "from the GitHub release into the installed package."
        ),
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="Release tag to pull from (default: the value of "
             "CAPELLINI_REFERENCES_TAG, else the package default).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-download even if the file already exists.",
    )
    args = parser.parse_args(argv)

    try:
        fetch_references(tag=args.tag, overwrite=args.overwrite)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
