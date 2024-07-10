# import pandas as pd
# import glob
# import copy
# from typing import Union, List, Dict
# from SynTemp.SynUtils.chemutils import standardize_rsmi
# from SynTemp.SynUtils.utils import run_shell_command
# from SynTemp.SynRule.rule_engine import RuleEngine


# class RuleExecutor:
#     @staticmethod
#     def reaction_prediction(
#         input_smiles: Union[str, List[str]],
#         rule_file_path: str,
#         prediction_type: str = "forward",
#         repeat_times: int = 1,
#         print_results: bool = False,
#         rule_class: Union[List[str], str] = None,
#         use_specific_rules: bool = False,
#         verbosity: int = 0,
#     ) -> List[str]:
#         """
#         Performs chemical reaction predictions on given SMILES strings using specified rule files. This function supports both forward
#         and backward prediction directions and can operate in single file mode or directory mode depending on the setup of the rule files.

#         Parameters:
#             - input_smiles (Union[str, List[str]]): Input SMILES or a list of SMILES strings.
#                         If a single string with '>>', it can indicate reactants and products.
#             - rule_file_path (str): Path to a GML file or a directory with multiple GML files for reaction rules.
#             - prediction_type (str, optional): Direction of the reaction prediction ('forward' or 'backward'). Defaults to 'forward'.
#             - repeat_times (int, optional): Number of iterations for rule application on the input. Defaults to 1.
#             - print_results (bool): Whether to print the results (for example in LaTeX format). Defaults to False.
#             - rule_class (Union[List[str], str], optional): Specific rule(s) to apply if using specific rules.
#             - use_specific_rules (bool, optional): Flag to use specific rules from rule_class or all rules in the directory.

#         Returns:
#             List[str]: Deduplicated list of predicted reaction SMILES strings, standardized.
#         """
#         # Prepare the input smiles list
#         if isinstance(input_smiles, str):
#             # Split by '>>' to separate reactants from products if needed
#             reaction_parts = input_smiles.split(">>")
#             smiles_list = (
#                 reaction_parts[0].split(".")
#                 if prediction_type == "forward"
#                 else reaction_parts[-1].split(".")
#             )
#         elif isinstance(input_smiles, list):
#             smiles_list = input_smiles
#         else:
#             raise TypeError(
#                 "input_smiles must be either a string or a list of strings."
#             )

#         predictions = []
#         # Define rule files based on parameters
#         if rule_file_path.endswith(".gml"):
#             # Single rule file scenario
#             rule_files = [rule_file_path]
#         elif use_specific_rules and rule_class:
#             # Specific rule files scenario
#             rule_files = [
#                 f"{rule_file_path}/{rule}.gml"
#                 for rule in (
#                     rule_class if isinstance(rule_class, list) else [rule_class]
#                 )
#             ]
#         else:
#             # Use all rule files in the directory
#             rule_files = glob.glob(f"{rule_file_path}/*.gml")

#             for rule_file in rule_files:
#                 reactions = RuleEngine.perform_reaction(
#                     rule_file_path=rule_file,
#                     initial_smiles=smiles_list,
#                     repeat_times=repeat_times,
#                     prediction_type=prediction_type,
#                     print_results=print_results,
#                     verbosity=verbosity,
#                 )
#                 predictions.extend(reactions)

#         # Standardize and deduplicate the predictions
#         predictions = [standardize_rsmi(value) for value in predictions]
#         if print_results:
#             run_shell_command(command="mod_post")
#         return list(set(predictions))

#     @staticmethod
#     def reaction_database_prediction(
#         database: Union[pd.DataFrame, List[Dict]],
#         rule_file_path: str,
#         original_rsmi_col: str = "reactions",
#         prediction_type: str = "forward",
#         repeat_times: int = 1,
#         print_results: bool = False,
#         rule_class_col: str = "class",
#         use_specific_rules: bool = False,
#         verbosity: int = 0,
#     ) -> Union[pd.DataFrame, List[Dict]]:
#         """
#         Applies reaction predictions across a database of entries using specified rule files, accommodating both forward
#         and backward directions. This method extracts SMILES strings from each database entry, applies reaction rules, and
#         aggregates the results back into the database.

#         Parameters:
#         - database (Union[pd.DataFrame, List[Dict]]): The database to process, which can be either a Pandas DataFrame or a list of dictionaries.
#         - rule_file_path (str): The file path to a GML file or a directory containing multiple GML files for reaction prediction.
#         - original_rsmi_col (str, optional): The key in the database entries that identifies where the original reaction SMILES strings are stored. Defaults to 'reactions'.
#         - prediction_type (str, optional): Specifies the direction of the reaction, either 'forward' or 'backward'. Defaults to 'forward'.
#         - repeat_times (int, optional): Specifies the number of iterations for the reaction rule application. Defaults to 1.
#         - print_results (bool): Print results in latex or not. Defaults to False.
#         - rule_class_col (str): Columns contain specific rule class want to apply.
#         - use_specific_rules (bool, optional): Flag to use specific rules from rule_class or all rules in the directory.
#         Returns:
#         - Union[pd.DataFrame, List[Dict]]: The updated database with each entry containing the results of the reaction predictions.
#         """
#         if isinstance(database, pd.DataFrame):
#             database_copy = database.copy().to_dict(orient="records")
#         elif isinstance(database, list):
#             database_copy = copy.deepcopy(database)
#         else:
#             raise ValueError(
#                 "Unsupported database type. Expected a Pandas DataFrame or a list of dictionaries."
#             )

#         for entry in database_copy:
#             input_smiles = entry.get(original_rsmi_col, "")
#             rule_class = entry.get(rule_class_col, "")
#             predictions = RuleExecutor.reaction_prediction(
#                 input_smiles,
#                 rule_file_path,
#                 prediction_type,
#                 repeat_times,
#                 print_results,
#                 rule_class,
#                 use_specific_rules,
#                 verbosity,
#             )
#             entry["predictions"] = predictions
#             entry["number_predictions"] = len(predictions)

#         if isinstance(database, pd.DataFrame):
#             return pd.DataFrame(database_copy)
#         return database_copy
