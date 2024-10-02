import os
import glob
import pandas as pd
from typing import List, Any, Dict, Optional, Tuple
from syntemp.pipeline import (
    rebalance,
    clean,
    run_aam,
    extract_its,
    rule_extract,
    write_gml,
)
from syntemp.SynUtils.utils import (
    prune_branches,
    reindex_data,
    save_database,
    setup_logging,
)


class AutoTemp:
    def __init__(
        self,
        rebalancing: bool = False,
        mapper_types: List[str] = ["local_mapper", "rxn_mapper", "graphormer"],
        id: str = "R-id",
        rsmi: str = "reactions",
        n_jobs: int = 4,
        verbose: int = 0,
        batch_size: int = 1000,
        job_timeout: Optional[int] = None,
        safe_mode: bool = False,
        save_dir: Optional[str] = None,
        fix_hydrogen: bool = True,
        refinement_its: bool = False,
        node_label_names: List[str] = ["element", "charge"],
        node_label_default: List[Any] = ["*", 0],
        edge_attribute: str = "order",
        max_radius: int = 3,
        reindex: bool = True,
        rerun_aam: bool = True,
        log_file: str = None,
        log_level: str = "INFO",
        clean_data: bool = True,
        get_random_hydrogen: bool = False,
        fast_process: bool = False,
    ):
        """
        Initializes the AutoTemp class with specified settings for processing chemical
        reaction data.

        Parameters:
        - rebalancing (bool): Whether to use synrbl to rebalance reactions.
        Defaults to False.
        - mapper_types (List[str]): List of mapper names to use for processing.
        Defaults to ['local_mapper'].
        - id (str): Identifier for reaction IDs. Defaults to 'R-id'.
        - rsmi (str): Identifier for reaction SMILES strings. Defaults to 'reactions'.
        - n_jobs (int): Number of parallel jobs to run. Defaults to 4.
        - verbose (int): Verbosity level for logging. Defaults to 0.
        - batch_size (int): Number of reactions to process in each batch.
        Defaults to 1000.
        - job_timeout (Optional[int]): Timeout for each batch job in seconds.
        None means no timeout.
        - safe_mode (bool): Flag to run processing in a safe mode, which might handle
        exceptions or special cases. Defaults to False.
        - save_dir (Optional[str]): Directory to save results.
        None means results are not saved.
        - fix_hydrogen (bool): Whether to fix hydrogen atoms in the ITS graphs.
        Defaults to True.
        - refinement_its (bool): Whether to refine incorrect ITS graphs.
        Defaults to False.
        - node_label_names (List[str]): Names of node labels in the graph.
        Defaults to ["element", "charge"].
        - node_label_default (List[Any]): Default values for node labels if unspecified.
        Defaults to ["*", 0].
        - edge_attribute (str): Name of the edge attribute in the graph.
        Defaults to "order".
        - max_radius (int): Maximum radius for node connectivity in the graph.
        Defaults to 3.
        - reindex (bool): Reindex the rule or not .Defaults to True,
        - rerun_aam: Run atom map for the whole data or already have the atom map data.
        Defaults to False,

        """
        self.rebalancing = rebalancing
        self.id = id
        self.rsmi = rsmi
        self.mapper_types = mapper_types
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.batch_size = batch_size
        self.job_timeout = job_timeout
        self.safe_mode = safe_mode
        self.save_dir = save_dir
        self.fix_hydrogen = fix_hydrogen
        self.refinement_its = refinement_its
        self.node_label_names = node_label_names
        self.node_label_default = node_label_default
        self.edge_attribute = edge_attribute
        self.max_radius = max_radius
        self.reindex = reindex
        self.rerun_aam = rerun_aam
        self.clean_data = clean_data
        self.get_random_hydrogen = get_random_hydrogen
        self.fast_process = fast_process

        # log_level = getattr(logging, log_level.upper(), None)
        # if not isinstance(log_level, int):
        #     raise ValueError(f"Invalid log level: {log_level}")
        self.logger = setup_logging("INFO", log_file)

    def temp_extract(
        self, data: List[Dict[str, Any]], lib_path: str = None
    ) -> Tuple[Any, List[Dict[str, Any]], List[Dict[str, Any]]]:
        """
        Processes the given chemical reaction data through cleaning, atom-atom mapping,
        ITS extraction, and rule extraction.

        Parameters:
        - data (List[Dict[str, Any]]): Input dataset containing chemical reactions.

        Returns:
        - Tuple containing:
          - Extracted rules.
          - List of ITS data that was incorrectly processed.
          - List of ITS data with uncertain hydrogen adjustments.
        """

        if isinstance(data, pd.DataFrame):
            data = data.to_dict("records")
            self.logger.info("Data converted to list of dictionaries.")

        # Step 1: rebalance and clean the data
        if self.rerun_aam:
            if self.rebalancing:
                data = rebalance(data, self.rsmi, self.id, self.n_jobs, self.batch_size)
            if self.clean_data:
                clean_data = clean(data, self.id, self.rsmi, self.n_jobs)
            else:
                clean_data = data

            # Step 2: Run atom-atom mapping
            aam_data = run_aam(
                clean_data,
                self.mapper_types,
                "reactions",
                self.job_timeout,
                self.safe_mode,
            )
        else:
            aam_data = data

        if self.save_dir:
            save_database(aam_data, f"{self.save_dir}/aam.json.gz")

        # Step 3: Extract ITS graphs and categorize them
        self.logger.info("Extract ITS graphs and categorize them.")
        its_correct, its_incorrect, uncertain_hydrogen = extract_its(
            aam_data,
            self.mapper_types,
            self.batch_size,
            self.verbose,
            self.n_jobs,
            self.fix_hydrogen,
            self.refinement_its,
            self.save_dir,
            get_random_results=self.get_random_hydrogen,
            fast_process=self.fast_process,
        )

        # Step 4: Extract rules from the correct ITS graphs
        self.logger.info("Extract rules from the correct ITS graphs.")
        reaction_dicts, templates, hier_templates = rule_extract(
            its_correct,
            self.node_label_names,
            self.node_label_default,
            self.edge_attribute,
            self.max_radius,
            self.save_dir,
        )
        if lib_path is None:
            self.logger.info("Write Rules.")
            gml_rules = write_gml(templates, self.save_dir, "Cluster_id", "RC", True)
            return (
                gml_rules,
                reaction_dicts,
                templates,
                hier_templates,
                its_incorrect,
                uncertain_hydrogen,
            )
        else:
            from syntemp.lib_isomorphism import LibIsomorphism

            # Generate GML rules without saving to a path
            gml_rules = write_gml(templates, None, "Cluster_id", "RC", True)

            # Check rules for isomorphism and collect IDs that do not pass the check
            radius_0_lib_path = os.path.join(lib_path, "R0")
            radius_0_id = [
                templates[0][i]["Cluster_id"]
                for i, rule in enumerate(gml_rules[0])
                if not LibIsomorphism.lib_isomorphism(
                    rule=rule, lib_path=radius_0_lib_path
                )
            ]

            new_templates = [
                (
                    [t for t in layer if t["Cluster_id"] in radius_0_id]
                    if idx == 0
                    else layer
                )
                for idx, layer in enumerate(templates)
            ]
            new_hier_templates = [
                (
                    [t for t in layer if t["Cluster_id"] in radius_0_id]
                    if idx == 0
                    else layer
                )
                for idx, layer in enumerate(hier_templates)
            ]
            starting_indices = []
            for radius in range(self.max_radius + 1):
                rule_path = os.path.join(lib_path, f"R{radius}")
                num_rules = len(glob.glob(os.path.join(rule_path, "*.gml")))
                starting_indices.append(num_rules + 1)
            # Prune branches from templates and hierarchical templates
            templates = prune_branches(new_templates)
            templates = reindex_data(templates, starting_indices)
            hier_templates = prune_branches(new_hier_templates)
            hier_templates = reindex_data(hier_templates, starting_indices)
            gml_rules = write_gml(templates, self.save_dir, "Cluster_id", "RC", True)
            # Return all relevant data
            return (
                gml_rules,
                reaction_dicts,
                templates,
                hier_templates,
                its_correct,
                uncertain_hydrogen,
            )
