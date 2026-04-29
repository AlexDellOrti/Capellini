Installation
============

.. code-block:: bash

   cd Package/capellini
   pip install -e .

External tools must be on ``PATH`` **before** running ``pip install`` —
the install will abort otherwise.

External dependencies
---------------------

================== =============================================================
Tool               Install
================== =============================================================
``spacepharer``    ``conda install -c bioconda spacepharer``
``mmseqs2``        ``conda install -c bioconda mmseqs2``
``minced``         ``conda install -c bioconda minced``
``prodigal``       ``conda install -c bioconda prodigal``
``Rscript`` +      ``brew install r`` then ``R -e "BiocManager::install('dada2')"``
``dada2``
``micro``          ``brew install micro``
================== =============================================================

Bypass the build-time check with ``CAPELLINI_SKIP_DEP_CHECK=1 pip install -e .``
(intended only for CI / containers where the tools become available later).

Bundled reference files
-----------------------

CAPELLINI ships with two large reference files inside the package so they
don't need to be re-downloaded or re-derived:

* ``capellini/data/references/progenome16S.fasta``
* ``capellini/data/references/spacers/spacers_CompleteCollection.fasta``

To rebuild them from proGenomes3, set ``regenerate_16S_reference: true`` or
``regenerate_spacers_collection: true`` in the config.
