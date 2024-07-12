import logging
from typing import List, Any, Dict, Optional, Tuple
from SynTemp.pipeline import rebalance, clean, run_aam, extract_its, rule_extract


# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class AutoTemp:
    def __init__(
        self,
        rebalancing: bool = False,
        mapper_types: List[str] = ["local_mapper"],
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
    ):
        """
        Initializes the AutoTemp class with specified settings for processing chemical
        reaction data.

        Parameters:
        - rebalancing (bool): Whether to use synrbl to rebalance reactions. Defaults to False.
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
        # Step 1: rebalance and clean the data
        if self.rebalancing:
            data = rebalance(data, self.rsmi, self.id, self.n_jobs, self.batch_size)
            print(data[0])
        clean_data = clean(data, self.id, self.rsmi, self.n_jobs)

        # Step 2: Run atom-atom mapping
        aam_data = run_aam(
            clean_data, self.mapper_types, "reactions", self.job_timeout, self.safe_mode
        )

        # Step 3: Extract ITS graphs and categorize them
        its_correct, its_incorrect, uncertain_hydrogen = extract_its(
            aam_data,
            self.mapper_types,
            self.batch_size,
            self.verbose,
            self.n_jobs,
            self.fix_hydrogen,
            self.refinement_its,
            self.save_dir,
        )
        if lib_path is None:
            # Step 4: Extract rules from the correct ITS graphs
            rules, reaction_dicts, templates, hier_templates = rule_extract(
                its_correct,
                self.node_label_names,
                self.node_label_default,
                self.edge_attribute,
                self.max_radius,
                self.save_dir,
            )

            return (
                rules,
                reaction_dicts,
                templates,
                hier_templates,
                its_incorrect,
                uncertain_hydrogen,
            )
        else:
            pass
