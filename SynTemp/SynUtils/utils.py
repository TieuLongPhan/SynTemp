import os
from typing import List, Dict, Set, Any
from typing import Optional, Union, Callable, Tuple
import json
import pickle
import random
import subprocess


def save_database(database: list[dict], pathname: str = "./Data/database.json") -> None:
    """
    Save a database (a list of dictionaries) to a JSON file.

    Args:
        database: The database to be saved.
        pathname: The path where the database will be saved. Defaults to './Data/database.json'.

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
        pathname: The path from where the database will be loaded. Defaults to './Data/database.json'.

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
