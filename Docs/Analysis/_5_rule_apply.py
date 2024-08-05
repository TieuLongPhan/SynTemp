import copy
import glob
from typing import List, Dict, Tuple
import pathlib
import sys
import pathlib
import argparse
import logging


def setup_logging(log_dir, log_level):
    log_path = pathlib.Path(log_dir)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Create a custom logger
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)

    # Create handlers
    fh = logging.FileHandler(log_dir, mode="w")
    fh.setLevel(getattr(logging, log_level))

    # Create formatters and add it to handlers
    f_format = logging.Formatter("%(asctime)s - %(message)s")
    fh.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(fh)

    return logger


root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))


from syntemp.SynUtils.chemutils import (
    categorize_reactions,
    standardize_rsmi,
)
from syntemp.SynUtils.utils import load_from_pickle, load_database, save_database
from syntemp.SynRule.reduce_reactions import ReduceReactions
from syntemp.SynRule.rule_engine import RuleEngine
import logging


class RuleBenchmark:
    """
    The RuleBenchmark class encapsulates functionalities for reaction modeling using the
    MØD toolkit.It provides methods for forward and backward prediction based on templates
    library.
    """

    @staticmethod
    def get_rule_files(entry, rule_file_path):
        """
        Determines the rule files to use based on the provided entry and rule class
        column.

        Parameters:
        - entry (dict): The dictionary containing rule class information.
        - rule_class_col (str or None): The key in the entry dict that holds the rule
        class ID(s).
        - rule_file_path (str): The base path to the rule files.

        Returns:
        - list of str: A list of file paths to the rule files.
        """

        if len(entry["rule_class"]) > 0:
            rule_ids = entry.get("rule_class", [])
            rule_ids = rule_ids if isinstance(rule_ids, list) else [rule_ids]
            rule_files = [f"{rule_file_path}/{rule_id}.gml" for rule_id in rule_ids]
        else:
            rule_files = glob.glob(f"{rule_file_path}/*.gml")

        return rule_files

    @staticmethod
    def reproduce_reactions(
        database: List[Dict],
        rule_file_path: str,
        original_rsmi_col: str = "reactions",
        verbosity: int = 0,
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Simulates chemical reactions for each entry in a molecular database, processing
        them in both forward and backward directions. It categorizes reactions based on
        matching the original reaction SMILES string and updates the database entries with
        the simulation results.

        Parameters:
        - database (List[Dict]): A list of dictionaries representing the entries to be
        processed.
        - rule_class_col (str): The key in the dictionaries identifying the rule file(s)
        for each entry.
        - rule_file_path (str): Path to the directory containing the rule files.
        - original_rsmi_col (str, optional): Key for the original reaction SMILES string
        in the dictionaries. Defaults to 'reactions'.
        - repeat_times (int, optional): The number of times to simulate the reaction for
        each entry. Defaults to 1.
        - use_specific_rules (bool, optional): If True, uses specific rule files
        identified by 'rule_class'. Otherwise, uses all rule files.
        - job_timeout (int): Timeout

        Returns:
        - Tuple[List[Dict], List[Dict]]: Two lists of updated dictionaries for forward and
        backward reactions, respectively.
        """

        updated_database_forward = copy.deepcopy(database)
        updated_database_backward = copy.deepcopy(database)

        for reaction_direction, updated_database in (
            ("forward", updated_database_forward),
            ("backward", updated_database_backward),
        ):
            for entry in updated_database:
                logging.info(f"Process reaction {entry[original_rsmi_col]}")
                entry["positive_reactions"] = []
                entry["unrank"] = []
                entry["unrank_raw"] = []

                reaction_side_index = 0 if reaction_direction == "forward" else 1
                initial_smiles_list = (
                    entry[original_rsmi_col].split(">>")[reaction_side_index].split(".")
                )

                rule_files = RuleBenchmark.get_rule_files(entry, rule_file_path)

                result_temp = {}
                for rule_file in rule_files:
                    rule_id = rule_file.split("/")[-1].split(".")[0]
                    try:
                        reactions = RuleEngine.perform_reaction(
                            rule_file_path=rule_file,
                            initial_smiles=initial_smiles_list,
                            prediction_type=reaction_direction,
                            verbosity=verbosity,
                        )
                    except FileNotFoundError as fnf_error:
                        reactions = []
                        logging.error(
                            f"File not found: {rule_file} - {fnf_error}"
                            + f"for {initial_smiles_list}"
                        )
                    except Exception as e:
                        reactions = []
                        logging.error(
                            f"Error processing file {rule_file}: {e}"
                            + f"for {initial_smiles_list}"
                        )
                    reactions = {standardize_rsmi(rxn) for rxn in reactions if rxn}
                    reactions = [rxn for rxn in reactions if rxn is not None]
                    if len(reactions) > 0:
                        result_temp[rule_id] = reactions

                    matched_reactions, _ = categorize_reactions(
                        reactions, entry[original_rsmi_col]
                    )

                    # Accumulate reactions
                    if matched_reactions:
                        entry["positive_reactions"].extend(matched_reactions)
                    entry["unrank"].extend(reactions)
                    entry["unrank_raw"].append(result_temp)
                    entry["positive_reactions"] = list(set(entry["positive_reactions"]))
                if len(entry["positive_reactions"]) > 0:
                    entry["positive_reactions"] = entry["positive_reactions"][0]
                else:
                    entry["positive_reactions"] = None
                entry["unrank"] = ReduceReactions.process_list_of_rsmi(
                    list(set(entry["unrank"]))
                )
        return updated_database_forward, updated_database_backward


def benchmark(args: argparse.Namespace) -> None:
    """Run benchmark with provided command line arguments."""
    logger = setup_logging(args.log_dir, args.log_level)
    logger.info("Start process...")
    logger.info("Loading database...")

    try:
        database = load_database(args.data_dir)
    except Exception as e:
        logger.error(f"Failed to load database: {e}")
        return

    batch_size = 500
    results_fw, results_bw = [], []

    for i in range(0, len(database), batch_size):
        logger.info(f"Processing from {i} to {i + batch_size}")
        batch_reactions = database[i : i + batch_size]

        try:
            fw, bw = RuleBenchmark.reproduce_reactions(
                database=batch_reactions,
                rule_file_path=args.rule_file_path,
                original_rsmi_col=args.original_rsmi_col,
            )
            results_fw.extend(fw)
            results_bw.extend(bw)
        except Exception as e:
            logger.error(f"Error processing batch from {i} to {i + batch_size}: {e}")

    try:
        save_database(results_fw, f"{args.save_dir}/fw_{args.radius}.json.gz")
        save_database(results_bw, f"{args.save_dir}/bw_{args.radius}.json.gz")
        logger.info("Results successfully saved.")
    except Exception as e:
        logger.error(f"Cannot save results due to: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run the synthesis rules benchmarking with configurable paths and options."
    )
    parser.add_argument(
        "--log_dir",
        type=str,
        default="log.txt",
        help="File path for logging output. Default is 'log.txt'.",
    )
    parser.add_argument(
        "--data_dir",
        type=str,
        default="data/test.json.gz",
        help="Path to the data directory containing the database file. Default is 'data/test.json.gz'.",
    )
    parser.add_argument(
        "--rule_file_path",
        type=str,
        default="rules.json",
        help="Path to the JSON file containing synthesis rules. Default is 'rules.json'.",
    )
    parser.add_argument(
        "--original_rsmi_col",
        type=str,
        default="reactions",
        help="Column name in the database for the original SMILES notation of reactions. Default is 'reactions'.",
    )
    parser.add_argument(
        "--log_level",
        type=str,
        default="INFO",
        help="Logging level to use. Options are DEBUG, INFO, WARNING, ERROR, CRITICAL. Default is 'INFO'.",
    )
    parser.add_argument(
        "--save_dir",
        type=str,
        default="data/",
        help="Path to the directory where data will be saved. Default is 'data/'.",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=0,
        help="Radius of reaction center used. Default is 0.",
    )

    args = parser.parse_args()
    benchmark(args)
