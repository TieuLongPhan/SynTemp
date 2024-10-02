import os
import shutil
import pandas as pd
from typing import List, Any, Dict, Optional, Union, Tuple
from syntemp.SynChemistry.neutralize import Neutralize
from syntemp.SynChemistry.deionize import Deionize
from syntemp.SynAAM.atom_map_consensus import AAMConsensus
from syntemp.SynITS.its_extraction import ITSExtraction
from syntemp.SynITS.its_hadjuster import ITSHAdjuster
from syntemp.SynITS.its_refinement import ITSRefinement
from syntemp.SynRule.hierarchical_clustering import HierarchicalClustering
from syntemp.SynRule.rule_writing import RuleWriting
from syntemp.SynUtils.utils import save_to_pickle, collect_data, setup_logging
from synrbl import Balancer


logger = setup_logging()


def rebalance(
    data: Union[pd.DataFrame, List[Dict[str, Any]]],
    reaction_col: str = "reactions",
    id_col: str = "id",
    n_jobs: int = 1,
    batch_size: int = 10,
) -> List[Dict[str, Any]]:
    """
    Rebalances a dataset of chemical reactions by modifying them to a balanced state
    using the Balancer class.

    Parameters:
    - data (Union[pd.DataFrame, List[Dict[str, Any]]]): Input chemical reaction data
    as a pandas DataFrame or list of dicts.
    - reaction_col (str): Column name for the reaction strings in the data. Default is
    "reactions".
    - id_col (str): Column name for the reaction identifiers in the data. Default is
    "id".
    - n_jobs (int): Number of jobs to run in parallel. Default is 1.
    - batch_size (int): Size of the batch to process at once. Default is 10.

    Returns:
    - List[Dict[str, Any]]: A list of dictionaries containing the balanced reactions.
    """
    logger.info("Starting the rebalancing process.")

    # Convert DataFrame to list of dictionaries if necessary
    if isinstance(data, pd.DataFrame):
        data = data.to_dict("records")
        logger.info("Data converted to list of dictionaries.")

    # Initialize the Balancer class with specified parameters
    synrbl = Balancer(
        reaction_col=reaction_col,
        id_col=id_col,
        n_jobs=n_jobs,
        batch_size=batch_size,
    )

    # Process rebalancing and receive output as a list
    balanced_reactions = synrbl.rebalance(reactions=data, output_dict=False)

    # Update the original data with balanced reactions
    new_data = [
        dict(entry, **{reaction_col: balanced_reactions[idx]})
        for idx, entry in enumerate(data)
    ]

    logger.info("Rebalancing process completed.")
    return new_data


def clean(
    data: Union[pd.DataFrame, List[Dict[str, Any]]],
    id: str = "R-id",
    rsmi: str = "reactions",
    n_jobs: int = 4,
) -> List[Dict[str, Any]]:
    """
    Cleans and standardizes reaction data by converting DataFrame to list,
    neutralizing charges, and applying deionization processes to
    SMILES representations.

    Parameters:
    - data (Union[pd.DataFrame, List[Dict[str, Any]]]): The reaction data to clean,
    either as a DataFrame or a list of dictionaries.
    - id (str, optional): The column/key name for reaction IDs. Defaults to 'R-id'.
    - rsmi (str, optional): The column/key name for reaction SMILES strings. Defaults
    to 'reactions'.
    - n_jobs (int, optional): Number of parallel jobs for the neutralization process.
    Defaults to 4.

    Returns:
    - List[Dict[str, Any]]: A list of dictionaries with cleaned and standardized
    reaction data.
    """
    logger.info("Starting the cleaning process.")
    if isinstance(data, pd.DataFrame):
        data = data.to_dict("records")
        logger.info("Data converted to list of dictionaries.")
    data = [{"R-id": value[id], "reactions": value[rsmi]} for value in data]
    data = Neutralize.parallel_fix_unbalanced_charge(data, "reactions", n_jobs)
    logger.info("Neutralization process completed.")
    data = Deionize.apply_uncharge_smiles_to_reactions(data, Deionize.uncharge_smiles)
    logger.info("Deionization process completed.")
    data = [
        {"R-id": value["R-id"], "reactions": value["standardized_reactions"]}
        for value in data
    ]
    logger.info("Cleaning process completed.")
    return data


def run_aam(
    data: List[Dict[str, Any]],
    mapper_types: List[Any],
    rsmi_column: str = "reactions",
    job_timeout: Optional[int] = None,
    safe_mode: bool = False,
) -> List[Dict[str, Any]]:
    """
    Runs atom-atom mapping on reaction data using a consensus approach among specified
    mappers.

    Parameters:
    - data (List[Dict[str, Any]]): The reaction data.
    - mapper_types (List[Any]): List of mapper objects/types to be used.
    - rsmi_column (str, optional): The key name for reaction SMILES strings in the
    input data. Defaults to "reactions".
    - job_timeout (Optional[int], optional): Timeout for each batch job in seconds.
    If None, there is no timeout. Defaults to None.
    - safe_mode (bool, optional): Whether to use safe mode during mapping,
    which might handle exceptions or special cases. Defaults to False.

    Returns:
    - List[Dict[str, Any]]: Results of the atom mapping process as a list of
    dictionaries.
    """
    logger.info("Starting atom mapping consensus process.")
    aam = AAMConsensus(data, mappers=mapper_types)
    results = aam.batch_consensus(
        data,
        rsmi_column=rsmi_column,
        batch_size=len(data),
        job_timeout=job_timeout,
        safe_mode=safe_mode,
    )
    logger.info("Atom mapping consensus process completed.")
    return results


def extract_its(
    data: List[str],
    mapper_types: Optional[List[str]] = None,
    batch_size: int = 1000,
    verbose: int = 1,
    n_jobs: int = 4,
    fix_hydrogen: bool = True,
    refinement_its: bool = False,
    save_dir: Optional[str] = None,
    data_name: str = "",
    symbol: str = ">>",
    get_random_results: bool = False,
    fast_process: bool = False,
    job_timeout: int = 1,
) -> List[dict]:
    """
    Executes the extraction of ITS graphs from reaction data in batches,
    adjusts hydrogen if necessary, and saves the data.

    Parameters:
    - data (List[str]): List of reaction data as SMILES strings or equivalent.
    - mapper_types (Optional[List[str]]): Optional list of mapper names to
    use for processing. Defaults to common mappers if not specified.
    - batch_size (int): Number of reactions to process in each batch.
    - verbose (int): Level of verbosity in logging output (0, 1, or 2).
    - n_jobs (int): Number of concurrent jobs for parallel processing.
    - fix_hydrogen (bool): Flag to indicate if hydrogen atoms should be adjusted
    in the ITS graphs.
    - save_dir (Optional[str]): Directory path where the processed data should be
    saved.
    If not provided, defaults to the current directory.
    - data_name (str): Base name for the saved data files, which helps in
    identifying different datasets.

    Returns:
    - List[dict]: A list of dictionaries containing the processed ITS graph data,
    with keys for correct, incorrect, and uncertain mappings.
    """
    logger.info(f"Extracting ITS graph with {n_jobs} CPUs.")

    if not mapper_types:
        mapper_types = ["rxn_mapper", "graphormer", "local_mapper"]
        logger.info("No mapper_types provided. Using default mappers.")

    temp_dir = os.path.join(save_dir, "temp_batches") if save_dir else "temp_batches"
    os.makedirs(temp_dir, exist_ok=True)

    num_batches = len(data) // batch_size + (len(data) % batch_size > 0)
    for i in range(num_batches):
        logger.info(f"Processing batch {i + 1}/{num_batches}")
        start_index = i * batch_size
        end_index = min((i + 1) * batch_size, len(data))
        batch_data = data[start_index:end_index]

        batch_correct, batch_incorrect = ITSExtraction.parallel_process_smiles(
            batch_data,
            mapper_types,
            n_jobs=n_jobs,
            verbose=verbose,
            export_full=False,
            check_method="RC",
            symbol=symbol,
        )

        if fix_hydrogen:
            if i == 1 or (i % 10 == 0 and i >= 10):
                logger.info(f"Fixing hydrogen for batch {i + 1}/{num_batches}.")
            batch_processed = ITSHAdjuster.process_graph_data_parallel(
                batch_correct,
                "ITSGraph",
                n_jobs=n_jobs,
                verbose=verbose,
                get_random_results=get_random_results,
                fast_process=fast_process,
                job_timeout=job_timeout,
            )

            uncertain_hydrogen = [
                value for value in batch_processed if value["ITSGraph"] is None
            ]
            batch_correct = [
                value for value in batch_processed if value["ITSGraph"] is not None
            ]

            # Save uncertain hydrogen data for the batch
            save_to_pickle(
                uncertain_hydrogen,
                os.path.join(temp_dir, f"uncertain_hydrogen_{i}.pkl"),
            )

        # Save batch results to temporary files
        save_to_pickle(batch_correct, os.path.join(temp_dir, f"batch_correct_{i}.pkl"))
        save_to_pickle(
            batch_incorrect, os.path.join(temp_dir, f"batch_incorrect_{i}.pkl")
        )

    # Combine saved batch data
    logger.info("Combine batch data.")

    logger.info("Processing equivalent ITS correct")
    its_correct = collect_data(num_batches, temp_dir, "batch_correct_{}.pkl")
    logger.info("Processing unequivalent ITS correct")
    its_incorrect = collect_data(num_batches, temp_dir, "batch_incorrect_{}.pkl")
    try:
        all_uncertain_hydrogen = []
        if fix_hydrogen:
            logger.info("Processing ambiguous hydrogen-ITS")
            all_uncertain_hydrogen = collect_data(
                num_batches, temp_dir, "uncertain_hydrogen_{}.pkl"
            )
    except Exception as e:
        logger.error(f"{e}")
        all_uncertain_hydrogen = []

    # logger.info(f"Number of correct mappers before refinement: {len(its_correct)}")
    if refinement_its:
        logger.info("Refining unequivalent ITS correct")
        its_refine = ITSRefinement.process_graphs_in_parallel(
            its_incorrect, mapper_types, n_jobs, verbose
        )
        its_refine = [value for value in its_refine if value]
        refined_ids = {value["R-id"] for value in its_refine}
        its_incorrect = [
            value for value in its_incorrect if value["R-id"] not in refined_ids
        ]

        its_correct.extend(its_refine)

    logger.info(f"Number of correct mappers: {len(its_correct)}")
    logger.info(f"Number of incorrect mappers: {len(its_incorrect)}")
    logger.info(
        "Number of uncertain hydrogen:"
        + f"{len(data)-len(its_correct)-len(its_incorrect)}"
    )
    if save_dir:
        logger.info("Combining and saving data")

        save_to_pickle(
            its_correct, os.path.join(save_dir, f"{data_name}_its_correct.pkl.gz")
        )

        save_to_pickle(
            its_incorrect,
            os.path.join(save_dir, f"{data_name}_its_incorrect.pkl.gz"),
        )

        save_to_pickle(
            all_uncertain_hydrogen,
            os.path.join(save_dir, f"{data_name}_uncertain_hydrogen.pkl.gz"),
        )
    # Clean up temporary files
    shutil.rmtree(temp_dir)

    return its_correct, its_incorrect, all_uncertain_hydrogen


def rule_extract(
    data: List[Any],
    node_label_names: List[str] = ["element", "charge"],
    node_label_default: List[Any] = ["*", 0],
    edge_attribute: str = "order",
    max_radius: int = 3,
    save_path: Optional[str] = None,
) -> Tuple[Optional[List[Any]], Optional[List[Any]], Optional[List[Any]]]:
    """
    Automates the generation of templates and rules from hierarchical clustering of
    given data. Handles potential errors gracefully.

    Parameters:
    - data (List[Any]): Input data to process.
    - node_label_names (List[str]): Labels for the nodes in the hierarchical clustering.
    - node_label_default (List[Any]): Default values for node labels if not specified.
    - edge_attribute (str): The attribute used for defining edges.
    - max_radius (int): Maximum radius for the clustering algorithm.
    - save_path (Optional[str]): Path to save output files.

    Returns:
    - Tuple of Optional[List[Any]]: A tuple containing lists of rules, templates,
      and hierarchical templates extracted, or None for each if an error occurs.
    """
    try:
        hier_cluster = HierarchicalClustering(
            node_label_names, node_label_default, edge_attribute, max_radius
        )
        logger.info("Hierarchical clustering initialized successfully.")
        reaction_dicts, templates, hier_templates = hier_cluster.fit(data)
        logger.info("Clustering completed and data extracted.")

        if save_path:
            for obj, name in zip(
                [reaction_dicts, templates, hier_templates],
                ["data_cluster.pkl.gz", "templates.pkl.gz", "hier_templates.pkl.gz"],
            ):
                save_to_pickle(obj, f"{save_path}/{name}")
                logger.info(f"{name} successfully saved.")

        return reaction_dicts, templates, hier_templates
    except Exception as e:
        logger.error("An error occurred during template generation: %s", e)
        return None, None, None


def write_gml(
    template_data: List[Any],
    save_path: Optional[str] = None,
    id_column: str = "Cluster_id",
    rule_column: str = "RC",
    reindex: bool = True,
) -> List:
    """
    Process templates to extract and save rules.

    Parameters:
    - template_data (List[Any]): List of template data to process.
    - save_path (Optional[str]): Base directory where directories will be created and
    files saved.

    Returns:
    - List: A list of results from the rule extraction process.
    """
    rules = []
    for radius, template in enumerate(template_data):
        directory_path = None
        if save_path:
            directory_path = os.path.join(save_path, f"R{radius}")
            try:
                # Ensure directory exists
                os.makedirs(directory_path, exist_ok=True)
                logger.info(f"Ensured directory exists at {directory_path}")
            except Exception as e:
                logger.error(f"Failed to create directory {directory_path}: {e}")
                continue  # Skip this iteration on failure to create directory

        try:
            write = RuleWriting.auto_extraction(
                template,
                id_column=id_column,
                rule_column=rule_column,
                reindex=reindex,
                save_path=directory_path,
            )
            rules.append(write)
            logger.info(f"Rules extracted for template at radius {radius}")
        except Exception as e:
            logger.error(f"Error extracting rules for radius {radius}: {e}")

    return rules
