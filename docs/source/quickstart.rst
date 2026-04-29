Quick start
===========

Terminal UI
-----------

.. code-block:: bash

   capellini
   # or
   python -m capellini

The arrow-key menu lets you load a config, edit it, run the full pipeline,
run a single stage, run a custom selection of stages, validate inputs, and
inspect the resolved configuration. The path of the last config you loaded
is remembered at ``~/.capellini/last_config`` and re-used on the next
launch.

Programmatic API
----------------

.. code-block:: python

   from capellini import CapelliniConfig, CapelliniPipeline

   cfg = CapelliniConfig.from_yaml("/path/to/your_config.yaml")
   CapelliniPipeline(cfg).run_all()
