import copy
import glob
from typing import List, Dict, Tuple
from SynTemp.SynUtils.chemutils import (
    standardize_rsmi,
    remove_stereochemistry_from_reaction_smiles,
)
from SynTemp.SynRule.rule_engine import RuleEngine
from SynTemp.SynChemistry.sf_factory import SFFactory
from SynTemp.SynChemistry.sf_similarity import SFSimilarity
from SynTemp.SynChemistry.sf_maxfrag import SFMaxFrag


class RuleBenchmark:
    """
    The MØDModeling class encapsulates functionalities for reaction modeling using the MØD toolkit.
    It provides methods for forward and backward prediction based on templates library.
    """

    @staticmethod
    def categorize_reactions(
        reactions: List[str], target_reaction: str
    ) -> Tuple[List[str], List[str]]:
        """
        Sorts a list of reaction SMILES strings into two groups based on their match with a specified target reaction. The
        categorization process distinguishes between reactions that align with the target reaction and those that do not.

        Parameters:
        - reactions (List[str]): The array of reaction SMILES strings to be categorized.
        - target_reaction (str): The SMILES string of the target reaction used as the benchmark for categorization.

        Returns:
        - Tuple[List[str], List[str]]: A pair of lists, where the first contains reactions matching the target and the second
                                    comprises non-matching reactions.
        """
        match, not_match = [], []
        target_reaction = standardize_rsmi(target_reaction)
        for reaction_smiles in reactions:
            if reaction_smiles == target_reaction:
                match.append(reaction_smiles)
            else:
                not_match.append(reaction_smiles)
        return match, list(set(not_match))

    @staticmethod
    def reproduce_reactions(
        database: List[Dict],
        id_col: str,
        rule_file_path: str,
        original_rsmi_col: str = "reactions",
        repeat_times: int = 1,
        prior: bool = True,
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Simulates reactions from a database of molecular structures in both forward and backward directions, categorizes
        the reactions based on matching with the original reaction SMILES, and updates the database entries with the
        results of these simulations.

        The process involves interpreting the SMILES strings from the database entries, using them to simulate
        reactions based on specified rules, and then updating the entries with the outcomes. This is done separately
        for both forward and backward reaction directions, producing two updated databases.

        Parameters:
        - database (List[Dict]): A list of dictionaries representing the database entries to be processed.
        - id_col (str): The key in the dictionaries used to identify the rule file for each entry.
        - rule_file_path (str): The base path to the directory containing the rule files.
        - original_rsmi_col (str, optional): The key in the dictionaries for the original reaction SMILES string. Defaults to 'reactions'.
        - repeat_times (int, optional): The number of times to repeat the reaction simulation for each entry. Defaults to 1.
        - prior (bool, optional): If True, uses specific rule files identified by 'id_col'. If False, uses all rule files in the directory.

        Returns:
        - Tuple[List[Dict], List[Dict]]: A tuple containing two lists of dictionaries, with the first list representing
                                         the updated database for forward reactions and the second for backward reactions.
        """
        updated_database_forward = copy.deepcopy(database)
        updated_database_backward = copy.deepcopy(database)

        for reaction_direction, updated_database in (
            ("forward", updated_database_forward),
            ("backward", updated_database_backward),
        ):
            for entry in updated_database:
                entry["positive_reactions"] = []
                # entry["negative_reactions"] = []
                entry["unrank"] = []

                if prior:
                    rule_files = [f"{rule_file_path}/{entry[id_col]}.gml"]
                else:
                    rule_files = glob.glob(f"{rule_file_path}/*.gml")

                for rule_file in rule_files:
                    reaction_side_index = 0 if reaction_direction == "forward" else 1
                    initial_smiles_list = (
                        entry[original_rsmi_col]
                        .split(">>")[reaction_side_index]
                        .split(".")
                    )

                    # Process reactions for each rule file
                    reactions = RuleEngine.perform_reaction(
                        rule_file_path=rule_file,
                        initial_smiles=initial_smiles_list,
                        repeat_times=repeat_times,
                        prediction_type=reaction_direction,
                    )

                    reactions = list(
                        set([standardize_rsmi(value) for value in reactions])
                    )
                    matched_reactions, unmatched_reactions = (
                        RuleBenchmark.categorize_reactions(
                            reactions, entry[original_rsmi_col]
                        )
                    )

                    # Accumulate reactions
                    if matched_reactions:
                        entry["positive_reactions"].extend(matched_reactions)
                    # entry["negative_reactions"].extend(unmatched_reactions)
                    entry["unrank"].extend(reactions)
                entry["positive_reactions"] = list(set(entry["positive_reactions"]))
                if len(entry["positive_reactions"]) > 0:
                    entry["positive_reactions"] = entry["positive_reactions"][0]
                else:
                    entry["positive_reactions"] = None
                entry["unrank"] = list(set(entry["unrank"]))

        return updated_database_forward, updated_database_backward

    @staticmethod
    def TopKAccuracy(
        list_of_dicts: List[Dict[str, List[str]]],
        ground_truth_key: str,
        rank_list_key: str,
        k: int,
        ignore_stero: bool = False,
        scoring_function: str = SFSimilarity(['FCFP6']),
    ) -> float:
        """
        Calculates the top-k accuracy from a list of dictionaries based on the specified ground truth and ranking list keys.

        This function evaluates each dictionary in the list to determine if the ground truth value is within the top 'k' values
        of the rank list provided in each dictionary. Each dictionary is then marked with a 'top_{k}_correct' boolean indicating
        whether the prediction was correct within the top k predictions.

        Parameters:
        - list_of_dicts (List[Dict[str, List[str]]]): List of dictionaries, each containing at least the keys for ground truth
        and rank list.
        - ground_truth_key (str): The key in each dictionary which corresponds to the ground truth value. Assumes ground truth
        value is standardized using `standardize_rsmi`.
        - rank_list_key (str): The key in each dictionary which holds the list of ranked predictions.
        - k (int): The number of top elements in the rank list to consider for determining if the ground truth is correctly predicted.

        Returns:
        - float: The top-k accuracy as a percentage of the total dictionaries evaluated, rounded to two decimal places.

        Raises:
        - ValueError: If any dictionary does not contain the required keys.
        """
        correct = 0
        total = len(list_of_dicts)

        # if scoring_function == "Similarity":
        #     list_of_dicts = SFSimilarity.process_list_of_dicts(
        #         list_of_dicts, "unrank", list_fp
        #     )

        factory = SFFactory(scoring_function=scoring_function)
        list_of_dicts = factory.process_list_of_dicts(list_of_dicts, 'unrank')
        # Check if the required keys are present
        for item in list_of_dicts:
            if ground_truth_key not in item or rank_list_key not in item:
                raise ValueError(
                    f"Dictionary is missing required keys: {ground_truth_key} or {rank_list_key}"
                )

            top_k_predictions = item[rank_list_key][
                :k
            ]  # Get the top-k elements from the rank list

            if ignore_stero:
                base_rsmi = remove_stereochemistry_from_reaction_smiles(
                    item[ground_truth_key]
                )
                top_k_predictions = [
                    remove_stereochemistry_from_reaction_smiles(x)
                    for x in top_k_predictions
                ]
            else:
                base_rsmi = standardize_rsmi(item[ground_truth_key])

            if base_rsmi in top_k_predictions:
                item[f"top_{k}_correct"] = True
                correct += 1
            else:
                item[f"top_{k}_correct"] = False

        accuracy = (correct / total) * 100 if total > 0 else 0
        return round(accuracy, 2)
