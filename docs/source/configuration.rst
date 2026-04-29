Configuration
=============

CAPELLINI is driven by a single YAML configuration file that you provide.
The package does not ship any bundled configs: it simply remembers the
path of the last config you loaded — stored at
``~/.capellini/last_config`` — and re-uses it on the next run.

In the terminal UI, point CAPELLINI at your YAML via
**Settings → Load config**. Programmatic users pass the path directly:

.. code-block:: python

   from capellini import CapelliniConfig, CapelliniPipeline

   cfg = CapelliniConfig.from_yaml("/path/to/your_config.yaml")
   CapelliniPipeline(cfg).run_all()

Key parameters
--------------

==========================  ==================  ================================================================
Parameter                   Default             Meaning
==========================  ==================  ================================================================
``species_level``           ``false``           Genus-level (``false``) vs species-level (``true``) target resolution
``BACTERIA_TAXONOMY_RANK``  ``target_taxids``   Rank used to aggregate bacteria for the network stage
``PREVALENCE``              ``0.10``            Keep features present in ≥ 10 % of samples
``CRISPR_SMOOTH_ALPHA``     ``0.95``            Strength of taxonomy smoothing applied to ``W``
``LAM`` (η in the paper)    ``0.5``             Strength of CRISPR-informed abundance propagation
``N_STEPS``                 ``1``               Number of message-passing updates
``fdr``                     ``0.05``            SpacePHARER FDR threshold
``min_n_spacers``           ``3``               MinCED minimum spacers per array
==========================  ==================  ================================================================

Inputs
------

CAPELLINI expects, per cohort:

* Raw 16S rRNA amplicon FASTQ files (forward, reverse, or paired).
* A viral contig FASTA (e.g. ViroProfiler output).
* A sample metadata CSV used to align bacterial and viral abundances.
* The SILVA reference (Release 138.1) and SILVA taxmap.

All output folders under ``base/`` are created automatically.

Outputs
-------

For each study, CAPELLINI writes under ``Enhanced Networks/<study>/``:

* ``common/`` — aligned, prevalence-filtered ``V``, ``B``, and metadata tables.
* ``shrinkage/`` — Schäfer–Strimmer shrinkage correlations on the CLR-stacked
  :math:`Z = [B^{\mathrm{CLR}}\ V^{\mathrm{CLR}}]`.
* ``crispr_raw/`` and ``crispr_smooth/`` — raw and taxonomy-smoothed CRISPR
  matrices (:math:`W`, :math:`\tilde{W}`).
* ``xstar/`` — host-informed abundances :math:`Z^*` from convex and
  residual message-passing variants.
