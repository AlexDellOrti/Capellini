"""Interactive terminal UI for CAPELLINI."""

from __future__ import annotations

import importlib.resources as pkg_resources
import os
import shutil
import subprocess
import sys
import threading
from collections import deque
from pathlib import Path

import questionary
from rich.console import Console, Group
from rich.panel import Panel
from rich.live import Live
from rich.table import Table
from rich.text import Text

from capellini.config import CapelliniConfig
from capellini.pipeline import CapelliniPipeline

CONSOLE = Console()
LAST_CONFIG_FILE = Path.home() / ".capellini" / "last_config"

LOGO = r"""
              /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /
             /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   
            /   /    ▗▄▄▖ ▗▄▖ ▗▄▄▖ ▗▄▄▄▖▗▖   ▗▖   ▗▄▄▄▖▗▖  ▗▖▗▄▄▄▖      /   /   /
           /   /    ▐▌   ▐▌ ▐▌▐▌ ▐▌▐▌   ▐▌   ▐▌     █  ▐▛▚▖▐▌  █       /   /   /
          /   /     ▐▌   ▐▛▀▜▌▐▛▀▘ ▐▛▀▀▘▐▌   ▐▌     █  ▐▌ ▝▜▌  █      /   /   /
         /   /      ▝▚▄▄▖▐▌ ▐▌▐▌   ▐▙▄▄▖▐▙▄▄▖▐▙▄▄▖▗▄█▄▖▐▌  ▐▌▗▄█▄▖   /   /   /
        /   /   /                                                   /   /   /
       /   /   /       CRISPR-Abundance Phage-Evidence Linkage     /   /   /
      /   /   /                for Network Inference              /   /   /
     /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /                
    /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   /   / 
"""


# ── Section grouping for "Show config" ────────────────────────────────────────
CONFIG_SECTIONS: list[tuple[str, list[str]]] = [
    ("Global", [
        "base", "download_path",
        "input_fasta_folder", "dada2_folder", "mmseq_folder", "sp_folder",
        "procs_folder", "enhanced_networks_folder",
        "silva_ref_path", "silva_taxmap_path", "full_ncbi_taxonomy_path",
        "bacterial_raw_fasta_folder", "virus_fasta_name", "metadata_path",
        "species_level", "fresh_start", "ref_removal",
        "genes_reference_url", "bacContigs_reference_url", "protein_reference_url",
    ]),
    ("DADA2", ["direction", "bacteria_fasta_name", "fasta_generation"]),
    ("MMSeqs2", [
        "isolate_ref_16S", "mapping_saving", "min_bitscore", "max_matches",
        "add_taxonomy", "extend_taxonomy",
    ]),
    ("SpacePHARER", [
        "min_n_spacers", "min_length", "max_length", "fdr",
        "keep_spacers_collection", "remove_decomp_fasta",
    ]),
    ("ProCs", [
        "proteins_extraction_path", "clustering_path", "matrix_type",
        "save_single_bacgenome_collection", "keep_coords",
        "filter_1bac_1vir", "remove_collections", "batch_size",
    ]),
    ("Network — global", [
        "OUTPUT_ROOT", "OVERWRITE", "VERBOSE",
        "RUN_COMMON_ABUNDANCE", "RUN_SHRINKAGE_CORRELATIONS",
        "RUN_RAW_CRISPR_NETWORKS", "RUN_SMOOTH_CRISPR", "RUN_XSTAR",
    ]),
    ("Network — common abundance", [
        "PREVALENCE", "KEEP_COLUMN", "BACTERIA_TAXONOMY_RANK",
    ]),
    ("Network — CRISPR smoothing", [
        "BACTERIAL_RANKS", "BACTERIAL_WEIGHTS",
        "CRISPR_SMOOTH_ALPHA", "TRANSPOSE_RAW_CRISPR_AFTER_LOAD",
    ]),
    ("Network — X*", ["PSEUDOCOUNT", "LAM", "N_STEPS", "PRESERVE_SCALE"]),
    ("Network — study inputs", [
        "STUDY", "virus_abundance_raw", "bacteria_otu", "bacteria_taxonomy",
        "phage_host_predictions", "tax_bac_for_smoothing", "tax_vir",
        "viral_ranks", "viral_weights", "aggregate_viral_rank",
    ]),
]


# ── Config persistence ────────────────────────────────────────────────────────

def _read_last_config_path() -> Path | None:
    """Return the path of the last config the user loaded, or None."""
    if LAST_CONFIG_FILE.exists():
        try:
            p = LAST_CONFIG_FILE.read_text().strip()
            if p:
                return Path(p)
        except OSError:
            pass
    return None


def _write_last_config_path(path: Path) -> None:
    LAST_CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    LAST_CONFIG_FILE.write_text(str(path))


# ── Editor ────────────────────────────────────────────────────────────────────

def _edit_file(path: Path) -> None:
    if shutil.which("micro") is None:
        CONSOLE.print(
            "[red]The 'micro' editor is required but was not found on PATH.[/red]\n"
            "Install it: https://github.com/zyedidia/micro  (e.g. `brew install micro`)."
        )
        return
    CONSOLE.print(
        "[dim]micro shortcuts:[/dim] "
        "[bold]Ctrl+S[/bold] save  ·  "
        "[bold]Ctrl+Q[/bold] quit  ·  "
        "[bold]Ctrl+G[/bold] full help"
    )
    _pause("Press Enter to open the editor…")
    subprocess.call(["micro", str(path)])


# ── Display helpers ───────────────────────────────────────────────────────────

def _show_logo() -> None:
    CONSOLE.print(f"[cyan]{LOGO}[/cyan]")


def _refresh_screen() -> None:
    """Clear the terminal and re-show the logo before drawing a new menu."""
    CONSOLE.clear()
    _show_logo()


def _show_config(cfg: CapelliniConfig) -> None:
    for section, fields in CONFIG_SECTIONS:
        table = Table(title=section, show_header=True, header_style="bold cyan")
        table.add_column("Field", style="cyan", no_wrap=True)
        table.add_column("Value", style="white")
        for f in fields:
            if hasattr(cfg, f):
                table.add_row(f, str(getattr(cfg, f)))
        CONSOLE.print(table)


def _bundled_data_dir() -> Path:
    """Return the path of the installed capellini.data directory.

    Falls back gracefully when the ``references`` subfolder is not a Python
    package or doesn't exist (which can happen if the user deleted it).
    """
    return Path(str(pkg_resources.files("capellini.data").joinpath(".")))


def _bundled_reference_paths() -> tuple[Path, Path]:
    data = _bundled_data_dir()
    return (
        data / "references" / "progenome16S.fasta",
        data / "references" / "spacers" / "spacers_CompleteCollection.fasta",
    )


def _check_inputs(cfg: CapelliniConfig) -> tuple[list[tuple[str, bool, str]], list[tuple[str, bool, str]]]:
    """Run all input/dependency checks. Returns (paths, deps) row tuples."""
    paths: list[tuple[str, bool, str]] = []

    virus = Path(cfg.input_fasta_folder) / cfg.virus_fasta_name if cfg.virus_fasta_name else None
    paths.append(("Virus FASTA", virus is not None and virus.exists(), str(virus or "")))
    paths.append(("SILVA reference",
                  bool(cfg.silva_ref_path) and Path(cfg.silva_ref_path).exists(),
                  cfg.silva_ref_path or ""))
    paths.append(("SILVA taxmap",
                  bool(cfg.silva_taxmap_path) and Path(cfg.silva_taxmap_path).exists(),
                  cfg.silva_taxmap_path or ""))
    paths.append(("Bacterial raw FASTA folder",
                  bool(cfg.bacterial_raw_fasta_folder) and Path(cfg.bacterial_raw_fasta_folder).is_dir(),
                  cfg.bacterial_raw_fasta_folder or ""))
    paths.append(("Metadata",
                  bool(cfg.metadata_path) and Path(cfg.metadata_path).exists(),
                  cfg.metadata_path or ""))

    try:
        bundled_16s, bundled_spacers = _bundled_reference_paths()
        paths.append(("Bundled progenome16S.fasta", bundled_16s.exists(), str(bundled_16s)))
        paths.append(("Bundled spacers_CompleteCollection.fasta",
                      bundled_spacers.exists(), str(bundled_spacers)))
    except (ModuleNotFoundError, FileNotFoundError):
        paths.append(("Bundled progenome16S.fasta", False, "<not installed>"))
        paths.append(("Bundled spacers_CompleteCollection.fasta", False, "<not installed>"))

    deps: list[tuple[str, bool, str]] = []
    for tool in ("spacepharer", "minced", "Rscript", "prodigal", "micro"):
        deps.append((tool, shutil.which(tool) is not None, ""))
    mmseqs_ok = shutil.which("mmseqs") is not None or shutil.which("mmseqs2") is not None
    deps.append(("mmseqs/mmseqs2", mmseqs_ok, ""))

    return paths, deps


_REFERENCE_LABELS = {
    "Bundled progenome16S.fasta",
    "Bundled spacers_CompleteCollection.fasta",
}


def _preflight_full_pipeline(session: "_Session") -> bool:
    """Validate inputs before launching the full pipeline.

    Returns True if the pipeline can proceed. Returns False if the user
    must go back to the main menu (missing inputs/deps), or if they declined
    to fix missing references.
    """
    if not session.ensure_loaded():
        _pause()
        return False
    cfg = session.config
    paths, deps = _check_inputs(cfg)

    missing_refs = [(n, d) for (n, ok, d) in paths if not ok and n in _REFERENCE_LABELS]
    missing_paths = [(n, d) for (n, ok, d) in paths if not ok and n not in _REFERENCE_LABELS]
    missing_deps = [(n, d) for (n, ok, d) in deps if not ok]

    # Case 1: only the bundled references are missing → offer to download.
    if missing_refs and not missing_paths and not missing_deps:
        _refresh_screen()
        CONSOLE.print("[yellow]References not found:[/yellow]")
        for n, _ in missing_refs:
            CONSOLE.print(f"  • {n}")
        CONSOLE.print()
        if _confirm("Download the missing references now?", default=True):
            _refresh_screen()
            from capellini.fetch_references import fetch_references
            try:
                fetch_references(overwrite=False)
            except RuntimeError as exc:
                CONSOLE.print(f"[red]{exc}[/red]")
                _pause()
                return False
            _pause()
            # Re-check after the download
            return _preflight_full_pipeline(session)
        return False

    # Case 2: anything else is missing → abort with a summary.
    if missing_refs or missing_paths or missing_deps:
        _refresh_screen()
        CONSOLE.print("[red]Cannot start the full pipeline — missing items:[/red]\n")
        if missing_refs:
            CONSOLE.print("[bold]References:[/bold]")
            for n, d in missing_refs:
                CONSOLE.print(f"  • {n}  [dim]{d}[/dim]")
            CONSOLE.print()
        if missing_paths:
            CONSOLE.print("[bold]Input paths:[/bold]")
            for n, d in missing_paths:
                CONSOLE.print(f"  • {n}  [dim]{d or '<empty>'}[/dim]")
            CONSOLE.print()
        if missing_deps:
            CONSOLE.print("[bold]External dependencies:[/bold]")
            for n, _ in missing_deps:
                CONSOLE.print(f"  • {n}  [dim]not on PATH[/dim]")
            CONSOLE.print()
        _select(
            "",
            choices=[questionary.Choice(" » Back to main menu", "back")],
            default="back",
        )
        return False

    return True


def _validate_inputs(cfg: CapelliniConfig) -> None:
    paths, deps = _check_inputs(cfg)

    def make_table(title: str) -> Table:
        t = Table(title=title)
        t.add_column("Check", style="cyan", no_wrap=True)
        t.add_column("Result")
        return t

    def add_row(t: Table, name: str, ok: bool, detail: str = "") -> None:
        marker = "[green]OK[/green]" if ok else "[red]MISSING[/red]"
        t.add_row(name, f"{marker} {detail}".strip())

    paths_table = make_table("Input paths")
    for n, ok, d in paths:
        add_row(paths_table, n, ok, d)

    deps_table = make_table("External dependencies")
    for n, ok, d in deps:
        add_row(deps_table, n, ok, d)

    CONSOLE.print(paths_table)
    CONSOLE.print()
    CONSOLE.print(deps_table)


_TAIL_MAX_LINES = 30


def _live_progress(stages: list[str], runner) -> None:
    statuses = {s: "waiting" for s in stages}
    tail: deque[str] = deque(maxlen=_TAIL_MAX_LINES)
    current_stage: dict[str, str] = {"name": ""}
    lock = threading.Lock()

    labels = CapelliniPipeline.STAGE_LABELS

    def render_table() -> Table:
        t = Table(title="CAPELLINI pipeline")
        t.add_column("Stage", style="cyan")
        t.add_column("Status")
        for s in stages:
            mark = {
                "waiting": "[grey50]waiting[/grey50]",
                "running": "[yellow]running[/yellow]",
                "done": "[green]✓ done[/green]",
                "failed": "[red]✗ failed[/red]",
            }[statuses[s]]
            t.add_row(labels.get(s, s), mark)
        return t

    def render() -> Group:
        with lock:
            body_text = "\n".join(tail) if tail else "[dim](no output yet)[/dim]"
            title = labels.get(current_stage["name"], current_stage["name"]) or "output"
        panel = Panel(
            Text.from_markup(body_text),
            title=f"[bold]{title}[/bold] (last {_TAIL_MAX_LINES} lines)",
            border_style="grey50",
            height=_TAIL_MAX_LINES + 2,
        )
        return Group(render_table(), Text(""), panel)

    class _Refreshable:
        def __rich__(self):
            return render()

    # Save real terminal fds before redirecting
    saved_stdout_fd = os.dup(1)
    saved_stderr_fd = os.dup(2)
    real_term = os.fdopen(saved_stdout_fd, "w", buffering=1)
    live_console = Console(file=real_term, force_terminal=True)

    # Pipe to capture all writes to fd 1 / fd 2 (Python prints + subprocess output)
    r_fd, w_fd = os.pipe()
    saved_py_stdout = sys.stdout
    saved_py_stderr = sys.stderr

    def reader_loop() -> None:
        buf = bytearray()
        while True:
            try:
                chunk = os.read(r_fd, 4096)
            except OSError:
                return
            if not chunk:
                return
            buf.extend(chunk)
            while True:
                nl = buf.find(b"\n")
                if nl < 0:
                    break
                line = bytes(buf[:nl]).decode("utf-8", errors="replace").rstrip("\r")
                del buf[: nl + 1]
                with lock:
                    tail.append(line)

    reader = threading.Thread(target=reader_loop, daemon=True)

    try:
        os.dup2(w_fd, 1)
        os.dup2(w_fd, 2)
        os.close(w_fd)
        sys.stdout = os.fdopen(1, "w", buffering=1, closefd=False)
        sys.stderr = os.fdopen(2, "w", buffering=1, closefd=False)
        reader.start()

        with Live(_Refreshable(), console=live_console, refresh_per_second=8) as live:  # noqa: F841
            for s in stages:
                with lock:
                    tail.clear()
                    current_stage["name"] = s
                statuses[s] = "running"
                try:
                    runner(s)
                except Exception:
                    statuses[s] = "failed"
                    raise
                statuses[s] = "done"
    finally:
        try:
            sys.stdout.flush()
            sys.stderr.flush()
        except Exception:
            pass
        os.dup2(saved_stdout_fd, 1)
        os.dup2(saved_stderr_fd, 2)
        sys.stdout = saved_py_stdout
        sys.stderr = saved_py_stderr
        os.close(saved_stderr_fd)
        # writer end of the pipe is now fully closed (fd 1 & fd 2 restored);
        # the reader sees EOF and exits.
        reader.join(timeout=2)
        try:
            os.close(r_fd)
        except OSError:
            pass
        try:
            real_term.close()  # closes saved_stdout_fd
        except Exception:
            pass


# ── Config session ────────────────────────────────────────────────────────────

class _Session:
    def __init__(self) -> None:
        self.config_path: Path | None = None
        self.config: CapelliniConfig | None = None

    def ensure_loaded(self) -> bool:
        if self.config is not None:
            return True
        path = self._resolve_initial_path()
        if path is None:
            return False
        return self._load(path)

    def _resolve_initial_path(self) -> Path | None:
        path = _read_last_config_path()
        if path is not None and path.exists():
            return path
        CONSOLE.print(
            "[yellow]No previous config remembered.[/yellow] "
            "Use Settings → Load config to point CAPELLINI at one."
        )
        return None

    def _load(self, path: Path) -> bool:
        if not path.exists():
            CONSOLE.print(f"[red]Config not found:[/red] {path}")
            return False
        try:
            self.config = CapelliniConfig.from_yaml(path)
        except Exception as exc:  # noqa: BLE001
            CONSOLE.print(f"[red]Failed to load config:[/red] {exc}")
            return False
        self.config_path = path
        CONSOLE.print(f"[green]Loaded config:[/green] {path}")
        return True

    def reload(self) -> None:
        if self.config_path is not None:
            self._load(self.config_path)

    def switch(self, path: Path) -> bool:
        return self._load(path)


# ── Settings sub-menu ─────────────────────────────────────────────────────────

def _select(message: str, choices: list, default=None):
    """Arrow-key selector that returns None if the user cancels with Ctrl-C."""
    try:
        return questionary.select(message, choices=choices, default=default, qmark="»").ask()
    except KeyboardInterrupt:
        return None


def _ask_path(message: str, default: str) -> Path | None:
    answer = questionary.path(message, default=default).ask()
    if answer is None or answer.strip() == "":
        return None
    return Path(answer).expanduser()


def _confirm(message: str, default: bool = True) -> bool:
    return bool(questionary.confirm(message, default=default).ask())


def _pause(message: str = "Press Enter to return…") -> None:
    try:
        input(f"\n{message}")
    except (KeyboardInterrupt, EOFError):
        pass


def _settings_menu(session: _Session) -> None:
    while True:
        _refresh_screen()
        current = session.config_path
        CONSOLE.print(
            f"[dim]Current config:[/dim] {current if current else '[red]<none>[/red]'}\n"
        )
        choice = _select(
            "Settings",
            choices=[
                questionary.Choice("Load config (set / change current)", "load"),
                questionary.Choice("Edit current config", "edit"),
                questionary.Choice("Show current config", "show"),
                questionary.Choice("Validate inputs", "validate"),
                questionary.Separator(),
                questionary.Choice(" » Back to main menu", "back"),
            ],
            default="back",
        )
        if choice in (None, "back"):
            return

        if choice == "load":
            default_str = str(current) if current else str(Path.cwd())
            target = _ask_path("Config path", default=default_str)
            if target is None:
                continue
            if not target.exists():
                CONSOLE.print(f"[red]Not found:[/red] {target}")
                _pause()
                continue
            if session.switch(target):
                _write_last_config_path(target)
                CONSOLE.print("[green]Config remembered for next time.[/green]")

        elif choice == "edit":
            if current is None:
                CONSOLE.print(
                    "[yellow]No config loaded.[/yellow] "
                    "Use 'Load config' first."
                )
                _pause()
                continue
            _edit_file(current)
            session.reload()

        elif choice in {"show", "validate"}:
            if not session.ensure_loaded():
                _pause()
                continue
            _refresh_screen()
            if choice == "show":
                _show_config(session.config)
            else:
                _validate_inputs(session.config)

        _pause()


# ── Stage sub-menus ───────────────────────────────────────────────────────────

def _selected_stages_menu(session: _Session) -> None:
    if not session.ensure_loaded():
        return
    stages = list(CapelliniPipeline.STAGE_ORDER)
    labels = CapelliniPipeline.STAGE_LABELS
    while True:
        _refresh_screen()
        CONSOLE.print(
            "[dim]Toggle a stage with [space], confirm with [enter].\n"
            "To go back to the main menu, do not select any stage and press enter.[/dim]\n"
        )
        try:
            picked = questionary.checkbox(
                "Select stages to run",
                choices=[questionary.Choice(labels.get(s, s), s) for s in stages],
                qmark="»",
            ).ask()
        except KeyboardInterrupt:
            return
        if picked is None or not picked:
            return
        ordered = [s for s in stages if s in picked]
        _refresh_screen()
        pipeline = CapelliniPipeline(session.config)
        _live_progress(ordered, pipeline.run_stage)
        _pause()


def _single_stage_menu(session: _Session) -> None:
    if not session.ensure_loaded():
        return
    stages = list(CapelliniPipeline.STAGE_ORDER)
    labels = CapelliniPipeline.STAGE_LABELS
    while True:
        _refresh_screen()
        options = [questionary.Choice(labels.get(s, s), s) for s in stages]
        options.extend([questionary.Separator(), questionary.Choice(" » Back to main menu", "back")])

        choice = _select("Run a single stage", choices=options, default="back")
        if choice in (None, "back"):
            return
        _refresh_screen()
        pipeline = CapelliniPipeline(session.config)
        _live_progress([choice], pipeline.run_stage)
        _pause()


def _run_pipeline_menu(session: _Session) -> None:
    while True:
        _refresh_screen()
        choice = _select(
            "Run pipeline",
            choices=[
                questionary.Choice("Run full pipeline", "run_all"),
                questionary.Choice("Run selected stages", "run_selected"),
                questionary.Choice("Run single stage", "run_one"),
                questionary.Separator(),
                questionary.Choice(" » Back to main menu", "back"),
            ],
            default="run_all",
        )
        if choice in (None, "back"):
            return
        if choice == "run_selected":
            _selected_stages_menu(session)
            continue
        if choice == "run_one":
            _single_stage_menu(session)
            continue
        if choice == "run_all":
            if not _preflight_full_pipeline(session):
                continue
            _refresh_screen()
            pipeline = CapelliniPipeline(session.config)
            _live_progress(CapelliniPipeline.STAGE_ORDER, pipeline.run_stage)
            _pause()


def _demo_menu() -> None:
    _refresh_screen()
    CONSOLE.print("[bold]Demo[/bold]\n")
    CONSOLE.print("Nothing here yet! Stay tuned.\n")
    _select(
        "",
        choices=[questionary.Choice(" » Back to main menu", "back")],
        default="back",
    )


# ── Main loop ─────────────────────────────────────────────────────────────────

def main() -> None:
    """Entry point for the `capellini` CLI."""
    # Lightweight subcommand dispatch before launching the interactive UI.
    if len(sys.argv) > 1 and sys.argv[1] in {"fetch-references", "fetch_references"}:
        from capellini.fetch_references import main as _fetch_main
        sys.exit(_fetch_main(sys.argv[2:]))

    _show_logo()
    session = _Session()
    session.ensure_loaded()

    while True:
        _refresh_screen()

        choice = _select(
            "Main menu",
            choices=[
                questionary.Choice("Run pipeline", "run"),
                questionary.Choice("Demo", "demo"),
                questionary.Choice("Settings", "settings"),
                questionary.Choice("Fetch/Update reference FASTAs from GitHub release", "fetch_refs"),
                questionary.Separator(),
                questionary.Choice("Quit", "quit"),
            ],
            default="run",
        )
        if choice in (None, "quit"):
            return
        if choice == "run":
            _run_pipeline_menu(session)
            continue
        if choice == "demo":
            _demo_menu()
            continue
        if choice == "settings":
            _settings_menu(session)
            continue
        if choice == "fetch_refs":
            _refresh_screen()
            from capellini.fetch_references import fetch_references
            overwrite = _confirm(
                "Re-download even if the files already exist?", default=False
            )
            try:
                fetch_references(overwrite=overwrite)
            except RuntimeError as exc:
                CONSOLE.print(f"[red]{exc}[/red]")
            _pause()
            continue


if __name__ == "__main__":
    sys.exit(main() or 0)
