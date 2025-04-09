import os
import shutil
import random
import pandas as pd
from joblib import Parallel, delayed
from typing import List, Any, Dict, Optional, Union, Tuple

from synkit.Chem.Reaction.deionize import Deionize
from synkit.Chem.Reaction.neutralize import Neutralize
from synkit.Chem.Reaction.standardize import Standardize
from synkit.IO.nx_to_gml import NXToGML
from synkit.IO.debug import setup_logging
from synkit.IO.data_io import save_to_pickle, collect_data, save_list_to_file


from syntemp.SynAAM.atom_map_consensus import AAMConsensus
from syntemp.SynITS.its_extraction import ITSExtraction
from syntemp.SynITS.its_hadjuster import ITSHAdjuster
from syntemp.SynRule.hierarchical_clustering import HierarchicalClustering
from syntemp.SynRule.rule_writing import RuleWriting

from synrbl import Balancer


logger = setup_logging()
std = Standardize()


def normalize_rsmi_dict(
    data: Dict[str, Any], reaction_col: str = "reactions"
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Normalizes the reaction SMILES in a dictionary using a specified standardizer.

    Parameters:
    - data (Dict[str, Any]): The dictionary containing the reaction data.
    - reaction_col (str): The key in the dictionary under which reaction data is stored.

    Returns:
    - (Dict[str, Any], Dict[str, Any]): A tuple of dictionaries where
    the first item is the dictionary with the normalized reaction data,
    and the second item is a dictionary of the original data if normalization failed.
    """
    original_data = data.copy()
    try:
        data[reaction_col] = std.fit(data[reaction_col])
    except Exception as e:
        logger.error(f"Error occurred during standardization: {e}")
        # Return the original data as the issue when an error occurs
        return {}, original_data

    # Remove keys with None values, if any remain
    data = {k: v for k, v in data.items() if v is not None}
    return data, {}


def normalize_rsmi_list(
    data: Union[pd.DataFrame, List[Dict[str, Any]]],
    reaction_col: str = "reactions",
    n_jobs: int = 1,
    verbose: int = 0,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Normalizes a list of dictionaries or a DataFrame containing reaction SMILES.

    Parameters:
    - data (Union[pd.DataFrame, List[Dict[str, Any]]]): The data to normalize,
      can be a list of dictionaries or a DataFrame.
    - reaction_col (str): The column or key where the reaction data is stored.
    - n_jobs (int): The number of jobs to run in parallel (default 1).
    - verbose (int): The verbosity level of the parallel processing (default 0).

    Returns:
    - (List[Dict[str, Any]], List[Dict[str, Any]]): A tuple of lists where the
    first list contains dictionaries with normalized reaction data, and the second list
    contains dictionaries of the original data where normalization failed.
    """
    if isinstance(data, pd.DataFrame):
        data = data.to_dict(orient="records")

    results = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(normalize_rsmi_dict)(d, reaction_col) for d in data
    )

    normalized_data, issues = zip(*results)  # Unpack results into separate lists

    # Filter out empty dictionaries from normalized_data
    normalized_data = [d for d in normalized_data if d]
    issues = [d for d in issues if d]

    return list(normalized_data), list(issues)


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

    df = []
    for value in data:
        try:
            value[reaction_col] = Standardize().fit(value[reaction_col])
            df.append(value)
        except Exception as e:
            logger.error(f"Error occurred during standardization: {e}")
            pass
    data = df
    # Initialize the Balancer class with specified parameters
    synrbl = Balancer(
        reaction_col=reaction_col,
        id_col=id_col,
        n_jobs=n_jobs,
        batch_size=batch_size,
    )

    # Rebalance the reactions
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
    save_dir: Optional[str] = None,
    data_name: str = "",
    job_timeout: int = 60,
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
        )

        if fix_hydrogen:
            if i == 1 or (i % 10 == 0 and i >= 10):
                logger.info(f"Fixing hydrogen for batch {i + 1}/{num_batches}.")
            batch_processed = ITSHAdjuster().process_graph_data_parallel(
                batch_correct,
                "ITSGraph",
                n_jobs=1,
                verbose=verbose,
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

    logger.info(f"Number of correct mappers: {len(its_correct)}")
    logger.info(f"Number of incorrect mappers: {len(its_incorrect)}")
    logger.info(
        "Number of uncertain hydrogen:"
        + f"{len(data)-len(its_correct)-len(its_incorrect)}"
    )
    if save_dir:
        logger.info("Combining and saving data")
        meta_folder = os.path.join(save_dir, "meta")
        os.makedirs(meta_folder, exist_ok=True)

        save_to_pickle(
            its_correct, os.path.join(meta_folder, f"{data_name}_its_correct.pkl.gz")
        )

        save_to_pickle(
            its_incorrect,
            os.path.join(meta_folder, f"{data_name}_its_incorrect.pkl.gz"),
        )

        save_to_pickle(
            all_uncertain_hydrogen,
            os.path.join(meta_folder, f"{data_name}_uncertain_hydrogen.pkl.gz"),
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
        for key, value in enumerate(templates):
            templates[key] = temp_list(value, "RC", "cls_id", "gml")

        reaction_dicts = temp_list(reaction_dicts, "ITSGraph", None, "its_gml")
        reaction_dicts = temp_list(reaction_dicts, "GraphRules", None, "rc_gml")

        if save_path:
            temp_folder = os.path.join(save_path, "template")
            os.makedirs(temp_folder, exist_ok=True)
            for obj, name in zip(
                [reaction_dicts, templates, hier_templates],
                [
                    "data_cluster.list.json",
                    "templates.list.json",
                    "hier_templates.list.json",
                ],
            ):
                save_list_to_file(obj, f"{temp_folder}/{name}")
                logger.info(f"{name} successfully saved.")

        return reaction_dicts, templates, hier_templates
    except Exception as e:
        logger.error("An error occurred during template generation: %s", e)
        return None, None, None


def write_gml(
    template_data: List[Any],
    save_path: Optional[str] = None,
    id_column: str = "cls_id",
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


def temp_dict(
    data: Dict[str, Any], rc_key: str, id_key: str = None, rule_key: str = "gml"
) -> Dict[str, Any]:
    """
    Transforms a dictionary entry representing a reaction by converting
    networkx graph representations into GML format and removing
    the original reaction key.

    Parameters:
    - data (Dict[str, Any]): Dictionary containing reaction components.
    - rc_key (str): The key in the dictionary that holds the reaction component to transform.
    - id_key (str): The key in the dictionary that holds the identifier used in the transformation.
    - rule_key (str): The key under which the transformed data will be stored (default 'gml').

    Returns:
    - Dict[str, Any]: Updated dictionary with the transformed data and
    the original reaction key removed.
    """
    if id_key is None:
        rule_name = random.randint(-1000, 1000)
    else:
        rule_name = data[id_key]
    try:
        transformer = NXToGML()
        data[rule_key] = transformer.transform(data[rc_key], rule_name, reindex=True)
        data.pop(rc_key, None)  # Remove the original reaction component key
    except Exception as e:
        logger.error(f"Failed to transform data for key {rc_key} with ID {id_key}: {e}")
        data.pop(rc_key, None)  # Ensure the original key is removed even on failure
        data[rule_key] = None  # Set the rule key to None to indicate failure
    return data


def temp_list(
    data: List[Dict[str, Any]],
    rc_key: str = "RC",
    id_key: str = "R-id",
    rule_key: str = "gml",
    n_jobs: int = 1,
    verbose: int = 0,
) -> List[Dict[str, Any]]:
    """
    Processes a list of dictionaries by transforming each dictionary's specified reaction component
    into the GML format using parallel processing.

    Parameters:
    - data (List[Dict[str, Any]]): List of dictionaries to process.
    - rc_key (str): The key in the dictionaries for the reaction component to transform.
    - id_key (str): The identifier key used in the transformation process.
    - rule_key (str): The key under which the transformed data will be stored (default 'gml').
    - n_jobs (int): Number of parallel jobs to run (default 1).
    - verbose (int): Level of verbosity for parallel processing (default 0).

    Returns:
        List[Dict[str, Any]]: A list of dictionaries with updated transformation data.
    """
    result = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(temp_dict)(d, rc_key, id_key, rule_key) for d in data
    )
    return result
