"""MMSeqs2 stage: 16S reference, easy-search, and 3-layer NCBI/GCA assignment."""

from __future__ import annotations

import bz2
import logging
import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

from capellini.config import CapelliniConfig
from capellini.utils.io import open_maybe_gzip
from capellini.utils.taxonomy import build_rank_to_taxids

logger = logging.getLogger(__name__)

# ── Bundled 16S reference path ─────────────────────────────────────────────────
_BUNDLED_16S = Path(__file__).parent.parent / "data" / "references" / "progenome16S.fasta"

# ── Regex patterns (exact from notebook) ──────────────────────────────────────
RE_TAXID = re.compile(r"(?:\b(?:taxid|tax_id)\b\s*[:=]\s*|\btaxid[|:])(\d+)\b", re.IGNORECASE)
RE_ASSEMBLY = re.compile(r"\b(GC[AF]_\d{9}\.\d+)\b")
RE_BIOSAMPLE = re.compile(r"\b(SAM[END]\d+)\b")
RE_BIOPROJECT = re.compile(r"\b(PRJ[ENDB][A-Z]?\d+)\b")
RE_ACCESSION = re.compile(
    r"\b("
    r"(?:NZ_)?[A-Z]{2}_[0-9A-Z]{4,12}\.\d+"
    r"|"
    r"[A-Z]{1,2}\d{5,9}\.\d+"
    r")\b"
)


def _uniq_preserve(seq: List[str]) -> List[str]:
    seen: set = set()
    out = []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def extract_ids_from_header(header: str) -> Dict[str, List[str]]:
    """Extract all structured IDs from a FASTA header using compiled regex patterns.

    Args:
        header: FASTA record description string.

    Returns:
        Dict with keys taxids, assemblies, biosamples, bioprojects, accessions.
    """
    return {
        "taxids": _uniq_preserve(RE_TAXID.findall(header)),
        "assemblies": _uniq_preserve(RE_ASSEMBLY.findall(header)),
        "biosamples": _uniq_preserve(RE_BIOSAMPLE.findall(header)),
        "bioprojects": _uniq_preserve(RE_BIOPROJECT.findall(header)),
        "accessions": _uniq_preserve(RE_ACCESSION.findall(header)),
    }


def is_16s_gene(description: str) -> bool:
    """Return True if a FASTA description corresponds to a 16S rRNA gene.

    Args:
        description: FASTA record description string.

    Returns:
        True if the record is a 16S rRNA gene.
    """
    if ("16S" in description) or ("16s" in description):
        return True
    if "ribosomal RNA" in description:
        if "product" in description:
            d1 = description.split("product")
            return "16" in d1[1]
        return False
    return False


def get_reference_16s(cfg: CapelliniConfig) -> Path:
    """Return path to the 16S ProGenomes reference, using the bundled file if available.

    Modification 1: checks the bundled progenome16S.fasta first. If not present,
    downloads the full genes FASTA from ProGenomes3 and filters it for 16S records.

    Args:
        cfg: Populated CapelliniConfig instance.

    Returns:
        Path to the ready progenome16S.fasta reference.
    """
    reference_16S_path = Path(cfg.input_fasta_folder) / "progenome16S.fasta"

    # Check bundled file first (unless user explicitly asked to regenerate)
    if _BUNDLED_16S.exists() and not getattr(cfg, "regenerate_16S_reference", False):
        logger.info("Bundled progenome16S.fasta found — using it")
        if not reference_16S_path.exists():
            import shutil as _shutil
            reference_16S_path.parent.mkdir(parents=True, exist_ok=True)
            _shutil.copy2(str(_BUNDLED_16S), str(reference_16S_path))
        return reference_16S_path

    if reference_16S_path.exists():
        logger.info("Isolated 16S records found in input folder — skipping download")
        return reference_16S_path

    # Download and filter
    filename = os.path.basename(cfg.genes_reference_url)
    genes_reference_path = Path(cfg.download_path) / filename

    if not genes_reference_path.exists():
        logger.info("Downloading ProGenomes3 gene reference (~large) ...")
        urllib.request.urlretrieve(cfg.genes_reference_url, str(genes_reference_path))
    else:
        logger.info("ProGenomes3 gene reference found — skipping download")

    logger.info("Filtering for 16S records ...")
    reference_16S_path.parent.mkdir(parents=True, exist_ok=True)
    with bz2.open(str(genes_reference_path), "rt") as fasta_in, open(str(reference_16S_path), "w") as fasta_out:
        for record in tqdm(SeqIO.parse(fasta_in, "fasta"), desc="Filtering 16S"):
            if is_16s_gene(record.description):
                fasta_out.write(record.format("fasta"))
    logger.info("16S reference saved: %s", reference_16S_path)
    return reference_16S_path


def run_mmseqs_easy_search(
    bact_path: Path,
    reference_16s_path: Path,
    mmseq_folder: str,
) -> Path:
    """Run mmseqs easy-search (nucleotide mode) of 16S ASVs against the reference.

    Args:
        bact_path: Query FASTA (16S DADA2 bacteria sequences).
        reference_16s_path: Subject FASTA (progenome16S.fasta).
        mmseq_folder: Output directory for the .m8 file.

    Returns:
        Path to the output.m8 file.
    """
    mmseqs_bin = "mmseqs2" if shutil.which("mmseqs2") else "mmseqs"
    tmp_dir = tempfile.mkdtemp()
    mmseqs_output = Path(mmseq_folder) / "output.m8"

    logger.info("Running mmseqs easy-search")
    subprocess.run(
        [
            mmseqs_bin, "easy-search",
            str(bact_path),
            str(reference_16s_path),
            str(mmseqs_output),
            tmp_dir,
            "--search-type", "3",
        ],
        check=True,
    )
    return mmseqs_output


def parse_mmseqs_output(
    mmseqs_output: Path,
    min_bitscore: int,
    max_matches: int,
) -> pd.DataFrame:
    """Parse an mmseqs .m8 output file and return the top-scoring hits per query.

    Adds NCBI ID, Genome Accession ID, and Gene Index columns parsed from the target header.

    Args:
        mmseqs_output: Path to output.m8.
        min_bitscore: Minimum bitscore threshold.
        max_matches: Maximum number of hits to keep per query (by bitscore).

    Returns:
        Filtered DataFrame with NCBI ID, Genome Accession ID, Gene Index columns.
    """
    cols = ["query", "target", "pident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "tstart", "tend", "evalue", "bitscore"]
    m8_df = pd.read_csv(str(mmseqs_output), sep="\t", names=cols)
    m8_df["bitscore"] = pd.to_numeric(m8_df["bitscore"], errors="coerce")

    filtered_m8_df = m8_df[m8_df["bitscore"] >= min_bitscore].copy()

    topBitScore_df = (
        filtered_m8_df
        .sort_values(["query", "bitscore"], ascending=[True, False])
        .groupby("query", group_keys=False)
        .head(max_matches)
        .copy()
    )

    target_split = topBitScore_df["target"].str.split(".", expand=True)
    topBitScore_df["NCBI ID"] = target_split[0]
    topBitScore_df["Genome Accession ID"] = target_split[2].str.rsplit("_", n=1).str[0]
    topBitScore_df["Gene Index"] = target_split[2].str.rsplit("_", n=1).str[1]

    return topBitScore_df


def extract_ids_from_reference(reference_16s_path: Path) -> tuple[list, pd.DataFrame]:
    """Parse the 16S FASTA reference and extract IDs from every record header.

    Args:
        reference_16s_path: Path to progenome16S.fasta.

    Returns:
        Tuple (per_record_rows list, df DataFrame) where df has columns
        record_id, ncbi_id, assembly, biosample, bioproject, accessions, description.
    """
    per_record_rows = []
    n = 0
    with open_maybe_gzip(str(reference_16s_path)) as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            n += 1
            header = rec.description
            ids = extract_ids_from_header(header)
            per_record_rows.append({
                "record_id": rec.id,
                "taxid": ids["taxids"][0] if ids["taxids"] else None,
                "assembly": ids["assemblies"][0] if ids["assemblies"] else None,
                "biosample": ids["biosamples"][0] if ids["biosamples"] else None,
                "bioproject": ids["bioprojects"][0] if ids["bioprojects"] else None,
                "accessions": ";".join(ids["accessions"]) if ids["accessions"] else None,
                "description": header,
            })
            if n % 200_000 == 0:
                logger.debug("Parsed %s records ...", n)

    logger.debug("Done parsing %s records from 16S reference", n)
    df = pd.DataFrame(per_record_rows)

    tax_ids_pro = []
    for record_id in df["record_id"]:
        tax_ids_pro.append(int(record_id.split(".")[0]))
    df["ncbi_id"] = tax_ids_pro

    return per_record_rows, df


def build_top200_dicts(
    taxonomy_table: pd.DataFrame,
    topBitScore_df: pd.DataFrame,
) -> tuple[dict, dict, dict]:
    """Build the three ranked-hit dictionaries for the 3-layer NCBI ID assignment.

    Layer 1 (species): per-ASV ranked hits for ASVs with a Genus.
    Layer 2 (genus): pooled ranked hits for all ASVs in each genus.
    Layer 3 (family): per-ASV ranked hits for ASVs without a Genus but with a Family.

    Args:
        taxonomy_table: ASV taxonomy DataFrame (index = ASV names, Genus/Family columns).
        topBitScore_df: MMSeqs2 hits with query, NCBI ID, bitscore columns.

    Returns:
        Tuple (top_200_per_asv, top_200_per_genus, top_200_per_family).
    """
    # Layer 1
    top_200_per_asv: dict = {}
    for asv in taxonomy_table.index:
        gen = taxonomy_table.loc[asv, "Genus"] if "Genus" in taxonomy_table.columns else None
        if isinstance(gen, str) and gen.strip():
            mm_sub = topBitScore_df.loc[topBitScore_df["query"] == asv].copy()
            mm_sub = mm_sub.sort_values("bitscore", ascending=False, kind="mergesort")
            top_200_per_asv[asv] = [int(i) for i in mm_sub.head(200)["NCBI ID"].tolist()]

    # Layer 2
    top_200_per_genus: dict = {}
    for genus_name in taxonomy_table["Genus"].unique() if "Genus" in taxonomy_table.columns else []:
        if not isinstance(genus_name, str):
            continue
        sub = taxonomy_table[taxonomy_table["Genus"] == genus_name]
        mm_sub = topBitScore_df.loc[topBitScore_df["query"].isin(sub.index)].copy()
        mm_sub = mm_sub.sort_values(["query", "bitscore"], ascending=[True, False], kind="mergesort")
        top_200_per_genus[genus_name] = [int(i) for i in mm_sub.head(200)["NCBI ID"].tolist()]

    # Layer 3
    top_200_per_family: dict = {}
    for asv in taxonomy_table.index:
        fam = taxonomy_table.loc[asv, "Family"] if "Family" in taxonomy_table.columns else None
        gen = taxonomy_table.loc[asv, "Genus"] if "Genus" in taxonomy_table.columns else None
        if isinstance(fam, str) and not isinstance(gen, str):
            mm_sub = topBitScore_df.loc[topBitScore_df["query"] == asv].copy()
            mm_sub = mm_sub.sort_values("bitscore", ascending=False, kind="mergesort")
            top_200_per_family[asv] = [int(i) for i in mm_sub.head(200)["NCBI ID"].tolist()]

    return top_200_per_asv, top_200_per_genus, top_200_per_family


def pick_bounded(ranked_taxids: list, allowed_space: set, used_set: set):
    """Pick the first ranked taxid that is in allowed_space and not already used.

    Args:
        ranked_taxids: Ordered list of NCBI taxids (best first).
        allowed_space: Set of valid taxids for this rank/layer.
        used_set: Global set of already-assigned taxids.

    Returns:
        Chosen integer taxid, or None if none qualify.
    """
    for t in ranked_taxids:
        t = int(t)
        if t in allowed_space and t not in used_set:
            return t
    return None


def fallback_pick_from_space(allowed_space: set, used_set: set):
    """Deterministic fallback: pick the smallest unused taxid from the allowed space.

    Args:
        allowed_space: Set of valid taxids.
        used_set: Set of already-assigned taxids.

    Returns:
        Smallest unused int taxid, or None if the space is exhausted.
    """
    for t in sorted(allowed_space):
        if t not in used_set:
            return int(t)
    return None


def map_silva_to_progenomes_bounded(
    silva: pd.DataFrame,
    top_200_per_asv: dict,
    top_200_per_genus: dict,
    top_200_per_family: dict,
    genus_to_taxids: dict,
    family_to_taxids: dict,
    ncbi_taxids_by_genus: dict,
    ncbi_taxids_by_family: dict,
    allowed_universe: set,
    topBitScore_df: pd.DataFrame,
    progenomes_ref_df: pd.DataFrame,
    enforce_unique_taxids: bool = True,
    debug: bool = False,
) -> tuple:
    """3-layer NCBI taxid assignment: species, genus, and family resolution.

    Adds six new columns to the input silva DataFrame:
        progenomes_taxid_species, progenomes_taxid_genus, progenomes_taxid_family
        GCA_species, GCA_genus, GCA_family, GCA (consolidated)

    Args:
        silva: ASV taxonomy DataFrame.
        top_200_per_asv: Layer 1 ranked hits (per-ASV).
        top_200_per_genus: Layer 2 ranked hits (per-genus).
        top_200_per_family: Layer 3 ranked hits (per-ASV, family fallback).
        genus_to_taxids: ProGenomes genus name → set of taxids.
        family_to_taxids: ProGenomes family name → set of taxids.
        ncbi_taxids_by_genus: NCBI genus name → set of taxids.
        ncbi_taxids_by_family: NCBI family name → set of taxids.
        allowed_universe: Set of all NCBI taxids in the 16S reference.
        topBitScore_df: MMSeqs2 hits DataFrame (with NCBI ID, Genome Accession ID).
        progenomes_ref_df: DataFrame from extract_ids_from_reference (ncbi_id, record_id).
        enforce_unique_taxids: Enforce uniqueness across genera/families (Layer 2/3).
        debug: Print per-ASV assignment details.

    Returns:
        Tuple (silva_out, species_mapping, genus_mapping, family_mapping, genus_rep).
    """
    silva_out = silva.copy()

    def _is_str(x):
        return isinstance(x, str) and x.strip() != ""

    # GCA lookup tables
    topBitScore_df = topBitScore_df.copy()
    topBitScore_df["NCBI ID int"] = topBitScore_df["NCBI ID"].astype(int)

    _best_per_asv = (
        topBitScore_df
        .sort_values("bitscore", ascending=False)
        .drop_duplicates(subset=["query", "NCBI ID int"], keep="first")
        .set_index(["query", "NCBI ID int"])["Genome Accession ID"]
    )

    _best_per_ncbi = (
        topBitScore_df
        .sort_values("bitscore", ascending=False)
        .drop_duplicates(subset=["NCBI ID int"], keep="first")
        .set_index("NCBI ID int")["Genome Accession ID"]
    )

    _ref_gca = progenomes_ref_df[["ncbi_id", "record_id"]].copy()
    _ref_gca["gca"] = _ref_gca["record_id"].str.split(".").str[2].str.rsplit("_", n=1).str[0]
    _ref_ncbi_to_gca = (
        _ref_gca[_ref_gca["gca"].notna()]
        .assign(ncbi_id_int=lambda x: x["ncbi_id"].astype(int))
        .drop_duplicates(subset=["ncbi_id_int"], keep="first")
        .set_index("ncbi_id_int")["gca"]
    )

    def _resolve_gca(asv, ncbi_id):
        if pd.isna(ncbi_id):
            return pd.NA
        ncbi_int = int(ncbi_id)
        try:
            return _best_per_asv.loc[(asv, ncbi_int)]
        except KeyError:
            pass
        try:
            return _best_per_ncbi.loc[ncbi_int]
        except KeyError:
            pass
        try:
            return _ref_ncbi_to_gca.loc[ncbi_int]
        except KeyError:
            pass
        return pd.NA

    # Pre-build genus allowed spaces
    genus_space: dict = {}
    if "Genus" in silva_out.columns:
        for gen in silva_out["Genus"].dropna().unique():
            if not _is_str(gen):
                continue
            space = (
                set(genus_to_taxids.get(gen, set()))
                .intersection(set(ncbi_taxids_by_genus.get(gen, set())))
                .intersection(set(allowed_universe))
            )
            genus_space[gen] = space

    # ── Layer 1 — progenomes_taxid_species ────────────────────────────────
    species_mapping: dict = {}
    gca_species_mapping: dict = {}

    for asv in silva_out.index:
        gen = silva_out.at[asv, "Genus"] if "Genus" in silva_out.columns else None
        if not _is_str(gen):
            continue
        space = genus_space.get(gen, set())
        if not space:
            continue
        ranked = top_200_per_asv.get(asv, [])
        chosen = pick_bounded(ranked, space, set())
        if chosen is None:
            chosen = fallback_pick_from_space(space, set())
        if chosen is None:
            continue
        species_mapping[asv] = int(chosen)
        gca_species_mapping[asv] = _resolve_gca(asv, chosen)
        if debug:
            print(f"[SPECIES] asv={asv} genus={gen} -> {chosen}")

    # ── Layer 2 — progenomes_taxid_genus ──────────────────────────────────
    used_genus: set = set()
    genus_rep: dict = {}
    genus_rep_gca: dict = {}

    if "Genus" in silva_out.columns:
        for gen in sorted([g for g in silva_out["Genus"].unique() if _is_str(g)]):
            ranked = top_200_per_genus.get(gen, [])
            if not ranked:
                continue
            space = genus_space.get(gen, set())
            if not space:
                continue
            chosen = pick_bounded(ranked, space, used_genus if enforce_unique_taxids else set())
            if chosen is None:
                chosen = fallback_pick_from_space(space, used_genus if enforce_unique_taxids else set())
            if chosen is None:
                continue
            genus_rep[gen] = int(chosen)
            genus_rep_gca[gen] = _resolve_gca(gen, chosen)
            if enforce_unique_taxids:
                used_genus.add(int(chosen))
            if debug:
                print(f"[GENUS]   {gen} -> {chosen}")

    genus_mapping: dict = {}
    gca_genus_mapping: dict = {}
    for asv in silva_out.index:
        gen = silva_out.at[asv, "Genus"] if "Genus" in silva_out.columns else None
        if _is_str(gen) and gen in genus_rep:
            genus_mapping[asv] = genus_rep[gen]
            gca_genus_mapping[asv] = genus_rep_gca[gen]

    # ── Layer 3 — progenomes_taxid_family ─────────────────────────────────
    used_family: set = set(used_genus)
    family_mapping: dict = {}
    gca_family_mapping: dict = {}

    for asv in silva_out.index:
        gen = silva_out.at[asv, "Genus"] if "Genus" in silva_out.columns else None
        fam = silva_out.at[asv, "Family"] if "Family" in silva_out.columns else None

        if asv in genus_mapping:
            continue
        if not _is_str(fam):
            continue

        ranked = top_200_per_family.get(asv, [])
        space = (
            set(family_to_taxids.get(fam, set()))
            .intersection(set(ncbi_taxids_by_family.get(fam, set())))
            .intersection(set(allowed_universe))
        )
        if not space:
            continue

        chosen = pick_bounded(ranked, space, used_family if enforce_unique_taxids else set())
        if chosen is None:
            chosen = fallback_pick_from_space(space, used_family if enforce_unique_taxids else set())
        if chosen is None:
            continue

        family_mapping[asv] = int(chosen)
        gca_family_mapping[asv] = _resolve_gca(asv, chosen)
        if enforce_unique_taxids:
            used_family.add(int(chosen))
        if debug:
            print(f"[FAMILY]  fam={fam} asv={asv} -> {chosen}")

    # ── Write output columns ───────────────────────────────────────────────
    silva_out["progenomes_taxid_species"] = pd.Series(species_mapping).reindex(silva_out.index).astype("Int64")
    silva_out["progenomes_taxid_genus"] = pd.Series(genus_mapping).reindex(silva_out.index).astype("Int64")
    silva_out["progenomes_taxid_family"] = pd.Series(family_mapping).reindex(silva_out.index).astype("Int64")
    silva_out["GCA_species"] = pd.Series(gca_species_mapping).reindex(silva_out.index)
    silva_out["GCA_genus"] = pd.Series(gca_genus_mapping).reindex(silva_out.index)
    silva_out["GCA_family"] = pd.Series(gca_family_mapping).reindex(silva_out.index)
    silva_out["GCA"] = (
        silva_out["GCA_species"]
        .fillna(silva_out["GCA_genus"])
        .fillna(silva_out["GCA_family"])
    )

    for legacy_col in ["progenomes_taxid", "target_taxids", "taxid_matched_rank", "NCBI_taxid"]:
        if legacy_col in silva_out.columns:
            silva_out.drop(columns=[legacy_col], inplace=True)

    # Summary
    n = len(silva_out)
    n_sp = silva_out["progenomes_taxid_species"].notna().sum()
    n_gen = silva_out["progenomes_taxid_genus"].notna().sum()
    n_fam = silva_out["progenomes_taxid_family"].notna().sum()
    n_gca = silva_out["GCA"].notna().sum()

    print(f"\nAssignment summary ({n} ASVs total):")
    print(f"  progenomes_taxid_species  [Layer 1] : {n_sp:>5}  ({100*n_sp/n:.1f}%)")
    print(f"  progenomes_taxid_genus    [Layer 2] : {n_gen:>5}  ({100*n_gen/n:.1f}%)")
    print(f"  progenomes_taxid_family   [Layer 3] : {n_fam:>5}  ({100*n_fam/n:.1f}%)")
    print(f"  GCA assigned (any layer)            : {n_gca:>5}  ({100*n_gca/n:.1f}%)")

    comparison = silva_out[["progenomes_taxid_species", "progenomes_taxid_genus"]].dropna()
    differs = (comparison["progenomes_taxid_species"] != comparison["progenomes_taxid_genus"]).sum()
    print(f"  ASVs where species != genus taxid   : {differs} / {len(comparison)}")

    return silva_out, species_mapping, genus_mapping, family_mapping, genus_rep


def run_mmseqs2(cfg: CapelliniConfig, taxonomy_table: pd.DataFrame) -> pd.DataFrame:
    """Orchestrate the full MMSeqs2 stage and 3-layer NCBI/GCA assignment.

    Steps:
        1. Get/bundle 16S reference (Modification 1).
        2. Run mmseqs easy-search.
        3. Parse .m8 output.
        4. Extract IDs from reference FASTA.
        5. Build per-layer ranked hit dicts.
        6. Run map_silva_to_progenomes_bounded.
        7. Save silva_fixed CSV.
        8. Optionally delete downloaded reference.

    Args:
        cfg: Populated CapelliniConfig instance.
        taxonomy_table: Taxonomy table from the NCBI mapping stage.

    Returns:
        silva_fixed DataFrame with progenomes_taxid and GCA columns.
    """
    logger.info("MMSeqs2: starting")

    # Rename ASV indices to ASV_1, ASV_2 …
    taxonomy_table = taxonomy_table.copy()
    taxonomy_table.index = [f"ASV_{i+1}" for i in range(len(taxonomy_table))]

    reference_16s_path = get_reference_16s(cfg)
    bact_path = Path(cfg.input_fasta_folder) / cfg.bacteria_fasta_name

    mmseqs_output = run_mmseqs_easy_search(bact_path, reference_16s_path, cfg.mmseq_folder)
    topBitScore_df = parse_mmseqs_output(mmseqs_output, cfg.min_bitscore, cfg.max_matches)

    per_record_rows, df = extract_ids_from_reference(reference_16s_path)
    allowed_universe = set(df["ncbi_id"].astype(int))

    # Load full NCBI taxonomy for rank lookups
    logger.info("Loading full NCBI taxonomy for genus/family bounds")
    df_all_ncbis = pd.read_csv(cfg.full_ncbi_taxonomy_path, sep="\t")
    ncbi_taxids_by_genus = build_rank_to_taxids(df_all_ncbis, "genus")
    ncbi_taxids_by_family = build_rank_to_taxids(df_all_ncbis, "family")

    genus_to_taxids = (
        df_all_ncbis.groupby("genus")["taxid"]
        .apply(lambda s: set(map(int, s)))
        .to_dict()
    )
    family_to_taxids = (
        df_all_ncbis.groupby("family")["taxid"]
        .apply(lambda s: set(map(int, s)))
        .to_dict()
    )

    top_200_per_asv, top_200_per_genus, top_200_per_family = build_top200_dicts(
        taxonomy_table, topBitScore_df
    )

    silva_fixed, *_ = map_silva_to_progenomes_bounded(
        silva=taxonomy_table,
        top_200_per_asv=top_200_per_asv,
        top_200_per_genus=top_200_per_genus,
        top_200_per_family=top_200_per_family,
        genus_to_taxids=genus_to_taxids,
        family_to_taxids=family_to_taxids,
        ncbi_taxids_by_genus=ncbi_taxids_by_genus,
        ncbi_taxids_by_family=ncbi_taxids_by_family,
        allowed_universe=set(map(int, allowed_universe)),
        topBitScore_df=topBitScore_df,
        progenomes_ref_df=df,
        enforce_unique_taxids=True,
        debug=False,
    )

    # Build target_taxids column (species- or genus-level, falling back to family)
    if cfg.species_level:
        silva_fixed["target_taxids"] = silva_fixed["progenomes_taxid_species"].fillna(
            silva_fixed["progenomes_taxid_family"]
        ).astype("Int64")
    else:
        silva_fixed["target_taxids"] = silva_fixed["progenomes_taxid_genus"].fillna(
            silva_fixed["progenomes_taxid_family"]
        ).astype("Int64")

    # Save
    suffix = {"forward": "F", "reverse": "R", "paired": "P"}.get(cfg.direction, "F")
    out_path = Path(cfg.dada2_folder) / f"taxonomy_table_{suffix}_progenomeLikeNCBIIDs.csv"
    silva_fixed.to_csv(str(out_path))
    logger.info("Saved silva_fixed: %s", out_path)

    if cfg.ref_removal:
        filename = os.path.basename(cfg.genes_reference_url)
        genes_reference_path = Path(cfg.download_path) / filename
        try:
            os.remove(str(genes_reference_path))
            logger.info("Deleted %s", genes_reference_path.name)
        except FileNotFoundError:
            pass

    logger.info("MMSeqs2 stage complete")
    return silva_fixed
