<p align="center">
  <img src="docs/_static/logo.png" alt="CAPELLINI logo" width="320">
</p>

# CAPELLINI

**CRISPR-Abundance Phage-Evidence Linkage for Leveraging Interaction Network Inference**

CAPELLINI is an integrative method for inferring phage–host interactions from
paired 16S rRNA amplicon and viral metagenomics data. It combines
**dataset-specific CRISPR evidence** with **abundance-based statistical
network inference** so that each modality informs the other:

1. A **SpacePHARER–proGenomes3 bioinformatics workflow** that links
   16S-derived bacterial taxa to high-quality reference genomes and extracts
   cohort-specific CRISPR spacers, producing a high-confidence virus–host
   adjacency matrix `W`.
2. A **statistical network estimation and infusion framework** that learns
   sparse partial correlations from CLR-transformed bacterial and viral
   abundances (SPIEC-EASI / shrinkage correlation), then propagates abundance
   residuals across the CRISPR-derived network to obtain a host-informed
   abundance representation `Z*`.

Across IBD, colorectal cancer (CRC), and GI-GvHD gut cohorts, propagation
through the CRISPR network increases the agreement between abundance-derived
associations and independent host predictions (iPHoP) by **+13 % to +63 %**,
with enrichment concentrated among the most stable inferred edges.

> Paper: *CAPELLINI: CRISPR-Abundance Phage-Evidence Linkage for Leveraging
> Interaction Network Inference.* Pugno, Dell'Orti, Tran, Müller (ECCB 2026).

---

## Pipeline overview

```
Preflight → DADA2 → 3-layer NCBI ID Mapping → SpacePHARER Execution
         → Protein Clusters (ProCs) Estimation → Enhanced Networks Estimation
```

- **Preflight** — folder layout; optional fresh-start cleanup that preserves
  bundled references and the input virus FASTA.
- **DADA2** — denoise 16S reads to ASVs and assign SILVA taxonomy via the
  bundled `DADA2_Pipe.R` script.
- **3-layer NCBI ID Mapping** — download `names.dmp` and assign real NCBI
  taxids to the SILVA taxonomy table, then run a taxonomy-aware mapping of
  ASVs to proGenomes3 representative genomes via `mmseqs easy-search` with
  three-layer fallback (ASV → genus → family) and derivation of the
  `target_taxids` column.
- **SpacePHARER Execution** — filter the bundled spacer collection to the
  cohort, build SpacePHARER databases, and run `predictmatch` with FDR
  control to obtain the virus–host adjacency `W`.
- **Protein Clusters (ProCs) Estimation** — protein clustering of bacterial
  and viral proteins, building the ProCs presence/count matrix.
- **Enhanced Networks Estimation** — common-abundance preprocessing, CLR
  transformation, Schäfer–Strimmer shrinkage correlations, raw and
  taxonomy-smoothed CRISPR networks (`W̃ = (1−α)W + α K_vir W K_bac`),
  and X* message-passing propagation
  (`Z*_v = Z_v + η (Z_b P_h − Z_v)`, `Z*_b = Z_b + η (Z_v P_v − Z_b)`).

---

## Installation

```bash
cd Package/capellini
pip install -e .
```

External tools must be on `PATH` **before** running `pip install` —
the install will abort otherwise.

| Tool                  | Install                                                          |
|-----------------------|------------------------------------------------------------------|
| `spacepharer`         | `conda install -c bioconda spacepharer`                          |
| `mmseqs` / `mmseqs2`  | `conda install -c bioconda mmseqs2`                              |
| `minced`              | `conda install -c bioconda minced`                               |
| `prodigal`            | `conda install -c bioconda prodigal`                             |
| `Rscript` + `dada2`   | `brew install r` → `R -e "BiocManager::install('dada2')"`        |
| `micro`               | `brew install micro`                                             |

Bypass the build-time check with `CAPELLINI_SKIP_DEP_CHECK=1 pip install -e .`
(intended only for CI / containers where the tools become available later).

### Reference FASTAs

Two large reference files are not shipped with the source repo and must be
fetched separately from the GitHub release:

- `capellini/data/references/progenome16S.fasta` (~76 MB)
- `capellini/data/references/spacers/spacers_CompleteCollection.fasta` (~74 MB)

The package provides a one-shot downloader that places them in the right
location inside the installed package:

```bash
capellini fetch-references
# or, equivalently:
capellini-fetch-references
```

You can also trigger it from the interactive UI:
**Main menu → Fetch reference FASTAs from GitHub release**.

If you skip this step, CAPELLINI will fall back to downloading the full
proGenomes3 source (~45 GB) and rebuilding both files at runtime — slower,
but useful if you intentionally want to regenerate them
(`regenerate_16S_reference: true` / `regenerate_spacers_collection: true`).

To pull from a specific release tag, set `CAPELLINI_REFERENCES_TAG` (default
is the latest tagged release):

```bash
CAPELLINI_REFERENCES_TAG=v0.1.0 capellini fetch-references
```

---

## Quick start

### Terminal UI

```bash
capellini
# or
python -m capellini
```

The arrow-key menu lets you create a default config, run the full pipeline,
run a single stage, run a custom selection of stages, validate inputs, and
inspect the resolved configuration.

### Programmatic API

```python
from capellini import CapelliniConfig, CapelliniPipeline

cfg = CapelliniConfig.from_yaml("capellini_config.yaml")
CapelliniPipeline(cfg).run_all()
```

CAPELLINI does not ship with prebuilt config files: it remembers the path of
the last config you loaded (in `~/.capellini/last_config`) and re-uses it on
the next run. Use **Settings → Load config** to point it at the YAML you
want to use.

---

## Inputs

CAPELLINI expects, per cohort:

- Raw 16S rRNA amplicon FASTQ files (forward, reverse, or paired).
- A viral contig FASTA (e.g. ViroProfiler output).
- A sample metadata table (CSV) used to align bacterial and viral
  abundances and to drop disease-associated columns when needed.
- The SILVA reference (Release 138.1) and SILVA taxmap.

All output folders under `base/` are created automatically.

---

## Key parameters

| Parameter                  | Default            | Meaning                                                                |
|----------------------------|--------------------|------------------------------------------------------------------------|
| `species_level`            | `false`            | Genus-level (`false`) vs species-level (`true`) target resolution      |
| `BACTERIA_TAXONOMY_RANK`   | `target_taxids`    | Rank used to aggregate bacteria for the network stage                  |
| `PREVALENCE`               | `0.10`             | Keep features present in ≥ 10 % of samples                             |
| `CRISPR_SMOOTH_ALPHA`      | `0.95`             | Strength of taxonomy smoothing applied to `W`                          |
| `LAM` (η in the paper)     | `0.5`              | Strength of CRISPR-informed abundance propagation                      |
| `N_STEPS`                  | `1`                | Number of message-passing updates                                      |
| `fdr`                      | `0.05`             | SpacePHARER FDR threshold                                              |
| `min_n_spacers`            | `3`                | MinCED minimum spacers per array                                       |

---

## Outputs

For each study, CAPELLINI writes under `Enhanced Networks/<study>/`:

- `common/` — aligned, prevalence-filtered `V`, `B`, and metadata tables.
- `shrinkage/` — Schäfer–Strimmer shrinkage correlations on the CLR-stacked
  `Z = [B^CLR  V^CLR]`.
- `crispr/` — raw and taxonomy-smoothed CRISPR matrices (`W`, `W̃`).
- `xstar/` — host-informed abundances `Z*` from convex and residual
  message-passing variants.

Intermediate stage outputs live under their respective folders
(`DADA2 output/`, `MMSeqs2 Output/`, `SpacePHARER output/`,
`Procs Estimations/`).

---

## Documentation

Full HTML documentation is built with Sphinx and the Read the Docs theme:

```bash
pip install -e ".[docs]"
cd docs
make html
open build/html/index.html
```

Live-reload while editing: `make livehtml`.

---

## Citation

If you use CAPELLINI, please cite the paper:

```
Pugno D., Dell'Orti A., Tran V., Müller C. L.
CAPELLINI: CRISPR-Abundance Phage-Evidence Linkage for Leveraging
Interaction Network Inference. ECCB 2026.
```

---

## License & contact

Code: <https://github.com/bio-datascience/capellini>
Corresponding author: **Daniele Pugno** — `D.Pugno@lmu.de`
