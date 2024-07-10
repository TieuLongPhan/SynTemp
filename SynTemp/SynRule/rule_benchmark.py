# import copy
# import os
# import glob
# from typing import List, Dict, Tuple
# from SynTemp.SynUtils.chemutils import (
#     categorize_reactions,
#     standardize_rsmi,
#     remove_stereochemistry_from_reaction_smiles,
# )
# from SynTemp.SynUtils.utils import load_from_pickle
# from SynTemp.SynRule.rule_engine import RuleEngine
# from SynTemp.SynRule.hier_engine import HierEngine
# from SynTemp.SynChemistry.sf_factory import SFFactory
# from SynTemp.SynChemistry.sf_similarity import SFSimilarity
# from SynTemp.SynChemistry.reduce_reactions import ReduceReactions
# import logging


# class RuleBenchmark:
#     """
#     The MØDModeling class encapsulates functionalities for reaction modeling using the MØD toolkit.
#     It provides methods for forward and backward prediction based on templates library.
#     """

#     @staticmethod
#     def reproduce_reactions(
#         database: List[Dict],
#         rule_class_col: str,
#         rule_file_path: str,
#         original_rsmi_col: str = "reactions",
#         repeat_times: int = 1,
#         use_specific_rules: bool = False,
#         verbosity: int = 0,
#         hierarchical=False,
#         max_radius: int = 3,
#         max_solutions: int = 100,
#         prune: bool = True,
#         prune_size: int = 1,
#         templates_threshold: int = 0.00,
#     ) -> Tuple[List[Dict], List[Dict]]:
#         """
#         Simulates chemical reactions for each entry in a molecular database, processing them in both forward
#         and backward directions. It categorizes reactions based on matching the original reaction SMILES string
#         and updates the database entries with the simulation results.

#         Parameters:
#             database (List[Dict]): A list of dictionaries representing the entries to be processed.
#             rule_class_col (str): The key in the dictionaries identifying the rule file(s) for each entry.
#             rule_file_path (str): Path to the directory containing the rule files.
#             original_rsmi_col (str, optional): Key for the original reaction SMILES string in the dictionaries. Defaults to 'reactions'.
#             repeat_times (int, optional): The number of times to simulate the reaction for each entry. Defaults to 1.
#             use_specific_rules (bool, optional): If True, uses specific rule files identified by 'rule_class'. Otherwise, uses all rule files.
#             job_timeout (int): Timeout

#         Returns:
#             Tuple[List[Dict], List[Dict]]: Two lists of updated dictionaries for forward and backward reactions, respectively.
#         """
#         updated_database_forward = copy.deepcopy(database)
#         updated_database_backward = copy.deepcopy(database)

#         for reaction_direction, updated_database in (
#             ("forward", updated_database_forward),
#             ("backward", updated_database_backward),
#         ):
#             for entry in updated_database:
#                 logging.info(f"Process reaction {entry[original_rsmi_col]}")
#                 entry["positive_reactions"] = []
#                 entry["unrank"] = []
#                 entry["unrank_raw"] = []

#                 reaction_side_index = 0 if reaction_direction == "forward" else 1
#                 initial_smiles_list = (
#                     entry[original_rsmi_col].split(">>")[reaction_side_index].split(".")
#                 )

#                 # Determine the rule files to use
#                 if use_specific_rules:
#                     rule_files = [
#                         f"{rule_file_path}/{rule}.gml"
#                         for rule in (
#                             entry[rule_class_col]
#                             if isinstance(entry[rule_class_col], list)
#                             else [entry[rule_class_col]]
#                         )
#                     ]
#                 else:
#                     rule_files = glob.glob(f"{rule_file_path}/*.gml")

#                     if hierarchical:
#                         root_rule_file_path = os.path.dirname(rule_file_path)
#                         hier_temp = load_from_pickle(
#                             f"{root_rule_file_path}/hier_rules.pkl.gz"
#                         )

#                         reactions, result_temp = HierEngine.hier_rule_apply(
#                             initial_smiles=initial_smiles_list,
#                             hier_temp=hier_temp,
#                             rule_file_path=root_rule_file_path,
#                             max_radius=max_radius,
#                             prediction_type=reaction_direction,
#                             max_solutions=max_solutions,
#                             prune=prune,
#                             prune_size=prune_size,
#                             templates_threshold=templates_threshold,
#                         )

#                         reactions = {standardize_rsmi(rxn) for rxn in reactions if rxn}
#                         reactions = [rxn for rxn in reactions if rxn is not None]

#                         for key, values in result_temp.items():
#                             standardized_reactions = {
#                                 standardize_rsmi(value) for value in values if value
#                             }
#                             standardized_reactions = [
#                                 rxn for rxn in standardized_reactions if rxn is not None
#                             ]
#                             reduced_reactions = ReduceReactions.process_list_of_rsmi(
#                                 standardized_reactions
#                             )
#                             result_temp[key] = reduced_reactions

#                         matched_reactions, _ = categorize_reactions(
#                             reactions, entry[original_rsmi_col]
#                         )
#                         if matched_reactions:
#                             entry["positive_reactions"] = matched_reactions
#                             # entry["negative_reactions"].extend(unmatched_reactions)
#                         entry["unrank"].extend(reactions)
#                         entry["unrank_raw"].append(result_temp)

#                     else:
#                         for rule_file in rule_files:
#                             try:
#                                 reactions = RuleEngine.perform_reaction(
#                                     rule_file_path=rule_file,
#                                     initial_smiles=initial_smiles_list,
#                                     repeat_times=repeat_times,
#                                     prediction_type=reaction_direction,
#                                     verbosity=verbosity,
#                                 )
#                             except FileNotFoundError as fnf_error:
#                                 reactions = []
#                                 logging.error(
#                                     f"File not found: {rule_file} - {fnf_error} for {initial_smiles_list}"
#                                 )
#                             except Exception as e:
#                                 reactions = []
#                                 logging.error(
#                                     f"Error processing file {rule_file}: {e} for {initial_smiles_list}"
#                                 )
#                             reactions = {
#                                 standardize_rsmi(rxn) for rxn in reactions if rxn
#                             }
#                             reactions = [rxn for rxn in reactions if rxn is not None]

#                             matched_reactions, _ = categorize_reactions(
#                                 reactions, entry[original_rsmi_col]
#                             )

#                             # Accumulate reactions
#                             if matched_reactions:
#                                 entry["positive_reactions"].extend(matched_reactions)
#                             entry["unrank"].extend(reactions)
#                         entry["positive_reactions"] = list(
#                             set(entry["positive_reactions"])
#                         )
#                 if len(entry["positive_reactions"]) > 0:
#                     entry["positive_reactions"] = entry["positive_reactions"][0]
#                 else:
#                     entry["positive_reactions"] = None
#                 entry["unrank"] = ReduceReactions.process_list_of_rsmi(
#                     list(set(entry["unrank"]))
#                 )
#         return updated_database_forward, updated_database_backward

#     @staticmethod
#     def TopKAccuracy(
#         list_of_dicts: List[Dict[str, List[str]]],
#         ground_truth_key: str,
#         rank_list_key: str,
#         k: int,
#         ignore_stero: bool = True,
#         scoring_function: str = SFSimilarity(["FCFP6"]),
#     ) -> float:
#         """
#         Calculates the top-k accuracy from a list of dictionaries based on the specified ground truth and ranking list keys.

#         This function evaluates each dictionary in the list to determine if the ground truth value is within the top 'k' values
#         of the rank list provided in each dictionary. Each dictionary is then marked with a 'top_{k}_correct' boolean indicating
#         whether the prediction was correct within the top k predictions.

#         Parameters:
#         - list_of_dicts (List[Dict[str, List[str]]]): List of dictionaries, each containing at least the keys for ground truth
#         and rank list.
#         - ground_truth_key (str): The key in each dictionary which corresponds to the ground truth value. Assumes ground truth
#         value is standardized using `standardize_rsmi`.
#         - rank_list_key (str): The key in each dictionary which holds the list of ranked predictions.
#         - k (int): The number of top elements in the rank list to consider for determining if the ground truth is correctly predicted.

#         Returns:
#         - float: The top-k accuracy as a percentage of the total dictionaries evaluated, rounded to two decimal places.

#         Raises:
#         - ValueError: If any dictionary does not contain the required keys.
#         """
#         correct = 0
#         total = len(list_of_dicts)
#         factory = SFFactory(scoring_function=scoring_function)
#         list_of_dicts = factory.process_list_of_dicts(list_of_dicts, "unrank")
#         # Check if the required keys are present
#         for item in list_of_dicts:
#             if ground_truth_key not in item or rank_list_key not in item:
#                 raise ValueError(
#                     f"Dictionary is missing required keys: {ground_truth_key} or {rank_list_key}"
#                 )

#             top_k_predictions = item[rank_list_key][
#                 :k
#             ]  # Get the top-k elements from the rank list

#             if ignore_stero:
#                 base_rsmi = remove_stereochemistry_from_reaction_smiles(
#                     item[ground_truth_key]
#                 )
#                 top_k_predictions = [
#                     remove_stereochemistry_from_reaction_smiles(x)
#                     for x in top_k_predictions
#                 ]
#             else:
#                 base_rsmi = standardize_rsmi(item[ground_truth_key])

#             if base_rsmi in top_k_predictions:
#                 item[f"top_{k}_correct"] = True
#                 correct += 1
#             else:
#                 item[f"top_{k}_correct"] = False

#         accuracy = (correct / total) * 100 if total > 0 else 0
#         return round(accuracy, 2)
