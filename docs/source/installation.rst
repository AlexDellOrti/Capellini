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

Reference FASTAs
----------------

Two large reference files are not shipped inside the source repo and must
be fetched separately from the GitHub release:

* ``capellini/data/references/progenome16S.fasta`` (~76 MB)
* ``capellini/data/references/spacers/spacers_CompleteCollection.fasta`` (~74 MB)

After installing CAPELLINI, run the bundled downloader to place them in
the right location inside the installed package:

.. code-block:: bash

   capellini fetch-references
   # or equivalently
   capellini-fetch-references

The same action is available from the interactive UI:
**Main menu → Fetch reference FASTAs from GitHub release**.

If you skip this step, CAPELLINI will fall back to downloading the full
proGenomes3 source and rebuilding both files at runtime — slower, but
useful if you intentionally want to regenerate them
(``regenerate_16S_reference: true`` / ``regenerate_spacers_collection: true``).

To pull from a specific release tag, set ``CAPELLINI_REFERENCES_TAG``:

.. code-block:: bash

   CAPELLINI_REFERENCES_TAG=v0.1.0 capellini fetch-references
