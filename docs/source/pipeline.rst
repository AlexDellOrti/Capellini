Pipeline overview
=================

.. code-block:: text

   Preflight → DADA2 → 3-layer NCBI ID Mapping → SpacePHARER Execution
            → Protein Clusters (ProCs) Estimation → Enhanced Networks Estimation

Stages
------

**Preflight**
    Folder layout; optional fresh-start cleanup that preserves bundled
    references and the input virus FASTA.

**DADA2**
    Denoise 16S reads to ASVs and assign SILVA taxonomy via the bundled
    ``DADA2_Pipe.R`` script.

**3-layer NCBI ID Mapping**
    Download ``names.dmp`` and assign real NCBI taxids to the SILVA
    taxonomy table, then run a taxonomy-aware mapping of ASVs to
    proGenomes3 representative genomes via ``mmseqs easy-search`` with
    three-layer fallback (ASV → genus → family) and derivation of the
    ``target_taxids`` column.

**SpacePHARER Execution**
    Filter the bundled spacer collection to the cohort, build SpacePHARER
    databases, and run ``predictmatch`` with FDR control to obtain the
    virus–host adjacency :math:`W`.

**Protein Clusters (ProCs) Estimation**
    Protein clustering of bacterial and viral proteins, building the ProCs
    presence/count matrix.

**Enhanced Networks Estimation**
    Common-abundance preprocessing, CLR transformation, Schäfer–Strimmer
    shrinkage correlations, raw and taxonomy-smoothed CRISPR networks

    .. math::

       \tilde{W} = (1 - \alpha) W + \alpha\, K_{\mathrm{vir}}\, W\, K_{\mathrm{bac}},

    and X* message-passing propagation

    .. math::

       Z^*_v = Z_v + \eta (Z_b P_h - Z_v), \quad
       Z^*_b = Z_b + \eta (Z_v P_v - Z_b).
