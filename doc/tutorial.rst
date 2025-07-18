.. _tutorials-and-examples:

Tutorials and Examples
======================

Overview
--------

**SynTemp** extracts graph-based reaction templates (partial Imaginary Transition State graphs) from reaction databases and clusters them hierarchically to reveal common transformation patterns. This section walks you through running **SynTemp** both as a Python API and from the command line, explains key parameters and outputs, and offers tips for adapting the workflow to your own data :cite:`phan2025syntemp`.

Use in a script
---------------

1. **Data preparation**  
   Prepare your reaction dataset as a list of Python dicts, each containing:  
   - **`id`** (e.g. `"R-id"`): a unique identifier for the reaction.  
   - **`reactions`** (e.g. `"reactions"`): the reaction SMILES string in `reactants>>products` format.

   .. code-block:: python

      data = [
          {
              "R-id": 0,
              "reactions": (
                  "COC(=O)C(CCCCNC(=O)OCc1ccccc1)NC(=O)Nc1cc(OC)cc"
                  "(C(C)(C)C)c1O.O>>COC(=O)C(CCCCN)NC(=O)Nc1cc(OC)cc"
                  "(C(C)(C)C)c1O.O=C(O)OCc1ccccc1"
              ),
          },
          # … more entries …
      ]

2. **Instantiate `AutoTemp`**  
   Configure the extraction pipeline via keyword arguments:

   .. code-block:: python

      from syntemp.auto_template import AutoTemp

      auto = AutoTemp(
          rebalancing=True,             # balance atoms across reactants/products
          mapper_types=[
              "rxn_mapper",
              "graphormer",
              "local_mapper",
          ],
          id="R-id",                    # key for reaction identifier
          rsmi="reactions",             # key for reaction SMILES
          n_jobs=1,                     # number of parallel workers
          verbose=2,                    # logging level (0–3)
          batch_size=1,                 # reactions per worker batch
          job_timeout=None,             # no per-job timeout
          safe_mode=False,              # skip strict mechanistic checks
          save_dir=None,                # directory to save outputs
          fix_hydrogen=True,            # infer missing hydrogens in ITS
      )

   Adjust `n_jobs` and `batch_size` for performance on larger datasets; use `safe_mode=True` to skip expensive validation if needed.

3. **Extract templates**  
   Run the extraction and collect results:

   .. code-block:: python

      reaction_dicts, templates, hier_templates, its_correct, uncertain_hydrogen = \
          auto.temp_extract(data, lib_path=None)

   - **`reaction_dicts`**: input records annotated with atom-mapping data.  
   - **`templates`**: lists of extracted ITS templates (as GML graphs).  
   - **`hier_templates`**: hierarchical clusters of templates.  
   - **`its_correct`**: booleans marking validated ITS graphs.  
   - **`uncertain_hydrogen`**: flags for reactions with ambiguous hydrogen placement.

4. **Inspect a core template**  
   View the GML representation of the first reaction’s primary template:

   .. code-block:: python

      core_tpl = templates[0]
      print(core_tpl[0]["gml"])

5. **Forward prediction with `SynReactor`**
   Convert a GML template back to an ITS graph and apply it to a substrate SMILES to predict the forward reaction. The `SynReactor` class from the *synkit* package (installed as a dependency of **SynTemp**) supports both forward and retrosynthetic modes via the `invert` flag.

   .. code-block:: python

      from synkit.IO import gml_to_its
      from synkit.Synthesis.Reactor.syn_reactor import SynReactor

      # 1. Select the substrate SMILES (reactant) from your data entry
      substrate_smiles = data[0]["reactions"].split(">>")[0]

      # 2. Convert the GML template to an ITS graph object
      its_graph = gml_to_its(core_tpl[0]["gml"])

      # 3. Initialize SynReactor for forward prediction
      reactor = SynReactor(
         substrate=substrate_smiles,
         template=its_graph,
         invert=False,      # False: forward direction; True: retrosynthesis
      )

      # 4. Run the reactor and retrieve results
      reaction_smarts    = reactor.smarts      # SMARTS pattern of the predicted reaction


      print("Reaction SMARTS:   ", reaction_smarts)


Use on the command line
-----------------------

You can run the same pipeline without writing Python:

.. code-block:: bash

   printf "R-id,reaction\n0,COC(=O)[C@H](CCCCNC(=O)OCc1ccccc1)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O>>COC(=O)[C@H](CCCCN)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O\n" \
     > test.csv

   python -m syntemp \
       --data_path test.csv \
       --rebalancing \
       --rerun_aam \
       --fix_hydrogen \
       --id R-id \
       --rsmi reaction \
       --mapper_types rxn_mapper graphormer local_mapper \
       --n_jobs 2 \
       --batch_size 10 \
       --log_file ./log.txt \
       --save_dir ./results

- **`--data_path`**: CSV file with header columns matching `--id` and `--rsmi`.  
- **`--rerun_aam`**: re-compute atom maps even if cached.  
- **`--save_dir`**: output directory; subfolders `meta/`, `templates/`, etc., will be created.

Reproduce full template extraction
----------------------------------

To replicate published results on the USPTO-50K dataset, run from the repository root:

.. code-block:: bash

   python -m syntemp \
       --data_path Data/USPTO_50K_original.csv \
       --log_file Data/Test/log.txt \
       --save_dir Data/Test/ \
       --rebalancing \
       --fix_hydrogen \
       --rerun_aam \
       --n_jobs 3 \
       --batch_size 1000 \
       --id ID \
       --rsmi reactions

This will process 50 000 reactions in parallel, infer ensemble atom mappings, complete ITS graphs with hydrogens, detect and extend reaction centers, and hierarchically cluster templates :cite:`phan2025syntemp`.

Tips and Troubleshooting
------------------------

- **Dependency conflicts**  
  If import errors arise (e.g., RDKit, RXNMapper), ensure your `requirements.txt` matches the repository’s pinned versions:  
  ```bash
  rdkit>=2024.3.5
  networkx>=3.3
  synrbl>=1.0.0
  synkit>=0.0.10
  # and, for ensemble AAMs:
  dgl==2.1.0
  dgllife==0.3.2
  localmapper>=0.1.5
  rxn-chem-utils>=1.6.0
  rxn-utils>=2.0.0
  rxnmapper>=0.4.1
  chython==1.78
  chytorch>=1.65
  chytorch-rxnmap>=1.4
  torch==2.2.0
  torchdata==0.7.1
  transformers==4.51.1 #temporary fix conflict

- **Hydrogen placement warnings**  
  Reactions flagged in `uncertain_hydrogen` may have ambiguous protonation states; inspect these manually or disable `fix_hydrogen` to skip automatic inference.

- **Performance tuning**  
  Increase `n_jobs` to utilize more CPU cores and raise `batch_size` for fewer, larger batches; set `safe_mode=True` to skip extensive validation on very large datasets.

For more details, browse the **SynTemp** repository on GitHub: https://github.com/TieuLongPhan/SynTemp/  
Enjoy exploring and extracting reaction templates with **SynTemp**!  

