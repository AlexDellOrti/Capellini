"""Top-level pipeline orchestrator."""

from __future__ import annotations

import logging
from typing import Any

from capellini.config import CapelliniConfig
from capellini.stages import (
    dada2 as dada2_stage,
    mmseqs2 as mmseqs2_stage,
    ncbi_mapping as ncbi_stage,
    network as network_stage,
    preflight as preflight_stage,
    procs as procs_stage,
    spacepharer as spacepharer_stage,
)

logger = logging.getLogger(__name__)


def _build_gca_target_set(silva_fixed, species_level: bool) -> set:
    """Replicate the SpacePHARER stage's target-GCA derivation."""
    if species_level:
        return (
            set(silva_fixed["GCA_species"].dropna())
            | set(silva_fixed["GCA_family"].dropna())
        )
    return (
        set(silva_fixed["GCA_genus"].dropna())
        | set(silva_fixed["GCA_family"].dropna())
    )


class CapelliniPipeline:
    """Run the CAPELLINI pipeline stages in order, sharing inter-stage state."""

    STAGE_ORDER = [
        "preflight",
        "dada2",
        "ncbi_mapping",
        "spacepharer",
        "procs",
        "network",
    ]

    STAGE_LABELS = {
        "preflight": "Preflight",
        "dada2": "DADA2",
        "ncbi_mapping": "3-layer NCBI ID Mapping",
        "spacepharer": "SpacePHARER Execution",
        "procs": "Protein Clusters (ProCs) Estimation",
        "network": "Enhanced Networks Estimation",
    }

    def __init__(self, config: CapelliniConfig) -> None:
        self.config = config
        self.state: dict[str, Any] = {}
        if not logging.getLogger().handlers:
            logging.basicConfig(
                level=logging.INFO,
                format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            )

    # ── Public API ───────────────────────────────────────────────────────────

    def run_all(self) -> None:
        """Run every stage in STAGE_ORDER."""
        for stage in self.STAGE_ORDER:
            self.run_stage(stage)

    def run_from(self, name: str) -> None:
        """Run from the given stage to the end."""
        idx = self.STAGE_ORDER.index(name)
        for stage in self.STAGE_ORDER[idx:]:
            self.run_stage(stage)

    def run_stage(self, name: str) -> Any:
        """Dispatch a single stage by name."""
        logger.info("=== Stage: %s ===", name)
        if name == "preflight":
            return preflight_stage.run_preflight(self.config)
        if name == "dada2":
            return dada2_stage.run_dada2(self.config)
        if name == "ncbi_mapping":
            tt = ncbi_stage.run_ncbi_mapping(self.config)
            self.state["taxonomy_table"] = tt
            silva_fixed = mmseqs2_stage.run_mmseqs2(self.config, tt)
            self.state["silva_fixed"] = silva_fixed
            self.state["gca_target_set"] = _build_gca_target_set(
                silva_fixed, self.config.species_level
            )
            return silva_fixed
        if name == "spacepharer":
            silva_fixed = self._require("silva_fixed")
            return spacepharer_stage.run_spacepharer(self.config, silva_fixed)
        if name == "procs":
            gca = self.state.get("gca_target_set")
            if gca is None:
                silva_fixed = self._require("silva_fixed")
                gca = _build_gca_target_set(silva_fixed, self.config.species_level)
                self.state["gca_target_set"] = gca
            pa = procs_stage.run_procs(self.config, gca)
            self.state["pa_matrix"] = pa
            return pa
        if name == "network":
            return network_stage.run_network(self.config)
        raise ValueError(f"Unknown stage: {name!r}")

    # ── Internal helpers ─────────────────────────────────────────────────────

    def _require(self, key: str) -> Any:
        if key not in self.state:
            raise RuntimeError(
                f"Stage requires '{key}' in pipeline state. "
                f"Run upstream stages first or call run_all()."
            )
        return self.state[key]
