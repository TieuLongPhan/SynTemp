import os
from typing import List, Dict, Any, Set, Tuple, Optional, Generator
import json
import pickle
import random
import subprocess
import logging
import pandas as pd
from sklearn.model_selection import train_test_split


def setup_logging(log_level: str = "INFO", log_filename: str = None) -> logging.Logger:
    """
    Configures the logging for an application. It sets up logging to either the console
    or a file based on whether a log filename is provided. The function adjusts the
    logging level dynamically according to the specified log level.

    Parameters
    ----------
    log_level : str, optional
        Specifies the logging level. Accepted values are 'DEBUG', 'INFO', 'WARNING',
        'ERROR', 'CRITICAL'. Default is 'INFO'.
    log_filename : str, optional
        Specifies the filename of the log file. If provided, logs will be written
        to this file. If None, logs will be written to the console. Default is None.

    Returns
    -------
    logging.Logger
        The configured logger object.

    Raises
    ------
    ValueError
        If the specified log_level is not recognized as a valid logging level.
    """
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {log_level}")

    logger = logging.getLogger()
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    if log_filename:
        os.makedirs(os.path.dirname(log_filename), exist_ok=True)
        logging.basicConfig(
            level=numeric_level, format=log_format, filename=log_filename, filemode="w"
        )
    else:
        logging.basicConfig(level=numeric_level, format=log_format)

    return logger


def save_database(database: list[dict], pathname: str = "./Data/database.json") -> None:
    """
    Save a database (a list of dictionaries) to a JSON file.

    Args:
        database: The database to be saved.
        pathname: The path where the database will be saved.
                    Defaults to './Data/database.json'.

    Raises:
        TypeError: If the database is not a list of dictionaries.
        ValueError: If there is an error writing the file.
    """
    if not all(isinstance(item, dict) for item in database):
        raise TypeError("Database should be a list of dictionaries.")

    try:
        with open(pathname, "w") as f:
            json.dump(database, f)
    except IOError as e:
        raise ValueError(f"Error writing to file {pathname}: {e}")


def load_database(pathname: str = "./Data/database.json") -> List[Dict]:
    """
    Load a database (a list of dictionaries) from a JSON file.

    Args:
        pathname: The path from where the database will be loaded.
                    Defaults to './Data/database.json'.

    Returns:
        The loaded database.

    Raises:
        ValueError: If there is an error reading the file.
    """
    try:
        with open(pathname, "r") as f:
            database = json.load(f)  # Load the JSON data from the file
        return database
    except IOError as e:
        raise ValueError(f"Error reading to file {pathname}: {e}")


def save_to_pickle(data: List[Dict[str, Any]], filename: str) -> None:
    """
    Save a list of dictionaries to a pickle file.

    Parameters:
    data (List[Dict[str, Any]]): A list of dictionaries to be saved.
    filename (str): The name of the file where the data will be saved.
    """
    with open(filename, "wb") as file:
        pickle.dump(data, file)


def load_from_pickle(filename: str) -> List[Any]:
    """
    Load data from a pickle file.

    Parameters:
    filename (str): The name of the pickle file to load data from.

    Returns:
    List[Any]: The data loaded from the pickle file.
    """
    with open(filename, "rb") as file:
        return pickle.load(file)


def stratified_random_sample(
    data: List[Dict], property_key: str, samples_per_class: int, seed: int = None
) -> List[Dict]:
    """
    Stratifies and samples data from a list of dictionaries based on a specified property.

    Parameters:
    - data (List[Dict]): The data to sample from, a list of dictionaries.
    - property_key (str): The key in the dictionaries to stratify by.
    - samples_per_class (int): The number of samples to take from each class.
    - seed (int): The seed for the random number generator for reproducibility.

    Returns:
    - List[Dict]: A list of sampled dictionaries.
    """

    if seed is not None:
        random.seed(seed)

    # Group data by the specified property
    stratified_data = {}
    for item in data:
        key = item.get(property_key)
        if key in stratified_data:
            stratified_data[key].append(item)
        else:
            stratified_data[key] = [item]

    # Sample data from each group
    sampled_data = []
    for key, items in stratified_data.items():
        if len(items) >= samples_per_class:
            sampled_data.extend(random.sample(items, samples_per_class))
        else:
            raise ValueError(
                f"Not enough data to sample {samples_per_class} items for class {key}"
            )

    return sampled_data


def ensure_directory_exists(directory_name):
    """
    Checks if a directory with the given name exists, and if not, creates it.

    Parameters:
    - directory_name (str): The name of the directory to check and potentially create.
    """
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
        print(f"Directory '{directory_name}' created.")
    else:
        print(f"Directory '{directory_name}' already exists.")


def run_shell_command(command="mod_post"):
    """
    Executes a shell command and prints its output.

    Parameters:
    - command (str): The shell command to execute.
    """
    try:
        # Execute the command and capture the output
        result = subprocess.run(
            command,
            shell=True,
            text=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print("Command output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error occurred while executing the command.")
        print("Error Code:", e.returncode)
        print("Error Message:", e.stderr)


def create_unique_value_dict(original_dict):
    """
    Creates a dictionary where each unique value from the original dictionary
    is mapped to the first key associated with that value.

    Parameters:
        original_dict (dict): The original dictionary with potentially duplicate values.

    Returns:
        dict: A new dictionary with unique values and their first corresponding keys.
    """
    new_dict = {}
    for key, value in original_dict.items():
        if value not in new_dict:
            new_dict[value] = key

    return new_dict


def train_val_test_split_df(
    df: pd.DataFrame,
    target: str,
    train_size: float = 0.8,
    val_size: float = 0.1,
    test_size: float = 0.1,
    random_state: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Splits a DataFrame into train, validation, and test sets with stratification.

    Parameters:
    - df (pd.DataFrame): The DataFrame to split.
    - target (str): The name of the column to stratify by.
    - train_size (float): Proportion of the dataset to include in the train split.
    - val_size (float): Proportion of the dataset to include in the validation split.
    - test_size (float): Proportion of the dataset to include in the test split.
    - random_state (int): The seed used by the random number generator.

    Returns:
    - Tuple[pd.DataFrame, pd.UtilityModel, pd.DataFrame]: A tuple containing the training
    set, validation set, and test set as DataFrames.

    Raises:
    - AssertionError: If the sum of train_size, val_size, and test_size does not equal 1.
    """
    # Ensure the sizes sum to 1
    assert (
        train_size + val_size + test_size == 1
    ), "The sum of train_size, val_size, and test_size must equal 1"

    # Split the data into training and a temporary set
    train_df, temp_df = train_test_split(
        df, train_size=train_size, stratify=df[target], random_state=random_state
    )

    # Calculate the relative size of the test set from the remaining data
    test_size_relative = test_size / (val_size + test_size)

    # Further split the temporary set into validation and test sets
    val_df, test_df = train_test_split(
        temp_df,
        train_size=1 - test_size_relative,
        stratify=temp_df[target],
        random_state=random_state,
    )

    return train_df, val_df, test_df


def modify_smiles(smiles):
    """
    Modify the input SMILES string by replacing '.[HH]' with '[H][H]'.

    Args:
    smiles (str): The input SMILES string representing a chemical reaction.

    Returns:
    str: The modified SMILES string.
    """
    return smiles.replace("[HH]", "[H][H]")


def prune_branches(
    data: List[List[Dict[str, Optional[int]]]]
) -> List[List[Dict[str, Optional[int]]]]:
    """
    Prunes branches in a hierarchical structure where any parent does not exist in the
    immediately preceding layer. Each layer is represented as a list of nodes, and
    each node is a dictionary with keys including 'Cluster_id' (integer) for the node's
    identifier and 'Parent' (integer or None for root nodes) for the node's parent ID
    in the previous layer.

    Args:
    - data (List[List[Dict[str, Optional[int]]]]): Hierarchical data where
    each sub-list represents a layer.

    Returns:
    - List[List[Dict[str, Optional[int]]]]: Pruned hierarchical data containing
    only nodes with existing parents.

    Each node's existence in a given layer depends on the presence of its 'Parent' ID in
    the set of 'Cluster_id's from the immediately previous layer. Root nodes in the
    first layer should have their 'Parent' set to None.
    """
    # Dictionary to keep track of existing cluster IDs for each layer
    existing_ids: Dict[int, Set[int]] = {0: {item["Cluster_id"] for item in data[0]}}

    # Process each layer after the first
    for layer_index in range(1, len(data)):
        layer = data[layer_index]
        # Store current layer's cluster IDs temporarily
        current_ids: Set[int] = set()

        new_layer: List[Dict[str, Optional[int]]] = []
        for item in layer:
            # Check if the single parent exists in the previous layer's IDs
            if item["Parent"] in existing_ids[layer_index - 1]:
                new_layer.append(item)
                current_ids.add(item["Cluster_id"])

        # Update the data layer with pruned nodes
        data[layer_index] = new_layer
        # Update existing IDs for this layer
        existing_ids[layer_index] = current_ids

    return data


def reindex_data(
    data: List[List[Dict[str, Optional[int]]]], starting_indices: List[int]
) -> List[List[Dict[str, Optional[int]]]]:
    """
    Reindexes 'Cluster_id', 'Parent', and optionally 'Child' fields in hierarchical data.
    First, each layer is normalized to start indexing from 0. Then, a layer-specific
    starting index from 'starting_indices' is added to each 'Cluster_id' in the layer.

    Parameters:
    - data (List[List[Dict[str, Optional[int]]]]): Hierarchical data
    where each sub-list represents a layer.
    - starting_indices (List[int]): The list of starting indices for 'Cluster_id'
    in each layer.

    Returns:
    - List[List[Dict[str, Optional[int]]]]: Hierarchical data with updated indices.
    """
    reindexed_data = []

    # Initialize id_map for reindexing all layers from 0, then apply starting_indices
    id_map = []
    for layer in data:
        current_map = {item["Cluster_id"]: idx for idx, item in enumerate(layer)}
        id_map.append(current_map)

    for layer_index, layer in enumerate(data):
        new_layer = []
        for item in layer:
            new_item = item.copy()
            # Apply the starting index for the current layer
            new_item["Cluster_id"] = (
                id_map[layer_index][item["Cluster_id"]] + starting_indices[layer_index]
            )

            # Update the 'Parent' field to the new index if it's not None
            if item["Parent"] is not None and layer_index > 0:
                new_item["Parent"] = (
                    id_map[layer_index - 1][item["Parent"]]
                    + starting_indices[layer_index - 1]
                )

            if "Child" in item and (layer_index + 1) < len(data):
                new_item["Child"] = [
                    id_map[layer_index + 1][child] + starting_indices[layer_index + 1]
                    for child in item["Child"]
                    if child in id_map[layer_index + 1]
                ]

            new_layer.append(new_item)
        reindexed_data.append(new_layer)

    return reindexed_data


def load_from_pickle_generator(file_path: str) -> Generator[Any, None, None]:
    """
    A generator that yields items from a pickle file where each pickle load returns a list
    of dictionaries.

    Paremeters:
    - file_path (str): The path to the pickle file to load.

    - Yields:
    Any: Yields a single item from the list of dictionaries stored in the pickle file.
    """
    with open(file_path, "rb") as file:
        while True:
            try:
                batch_items = pickle.load(file)
                for item in batch_items:
                    yield item
            except EOFError:
                break


def collect_data(num_batches: int, temp_dir: str, file_template: str) -> List[Any]:
    """
    Collects and aggregates data from multiple pickle files into a single list.

    Paremeters:
    - num_batches (int): The number of batch files to process.
    - temp_dir (str): The directory where the batch files are stored.
    - file_template (str): The template string for batch file names, expecting an integer
    formatter.

    Returns:
    List[Any]: A list of aggregated data items from all batch files.
    """
    collected_data: List[Any] = []
    for i in range(num_batches):
        file_path = os.path.join(temp_dir, file_template.format(i))
        for item in load_from_pickle_generator(file_path):
            collected_data.append(item)
    return collected_data
