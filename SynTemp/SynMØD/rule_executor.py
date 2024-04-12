import pandas as pd
import glob
import copy
from typing import Union, List, Dict
from SynTemp.SynUtils.chemutils import standardize_rsmi
from SynTemp.SynUtils.utils import run_shell_command, ensure_directory_exists
from SynTemp.SynMØD.rule_engine import RuleEngine


class RuleExecutor:
    @staticmethod
    def reaction_prediction(
        input_smiles: Union[str, List[str]],
        rule_file_path: str,
        prediction_type: str = "forward",
        repeat_times: int = 1,
        print_results: bool = False,
    ) -> List[str]:
        """
        Performs reaction predictions on a given SMILES input using specified rule files, accommodating both forward
        and backward directions. This function processes the input SMILES to determine the relevant reaction pathway
        and executes the reaction using the provided rule files, either in a single file mode or directory mode.

        Parameters:
        - input_smiles (Union[str, List[str]]): The input SMILES or a list of SMILES strings. If a single SMILES string
          is provided, it is expected to possibly contain '>>' to separate reactants and products.
        - rule_file_path (str): The file path to a GML file or a directory containing multiple GML files for reaction prediction.
        - prediction_type (str, optional): Specifies the direction of the reaction, either 'forward' or 'backward'. Defaults to 'forward'.
        - repeat_times (int, optional): Specifies the number of times the reaction should be attempted, allowing for multiple
          iterations over the rule application. Defaults to 1.
        - print_results (bool): Print results in latex or not. Defaults to False.

        Returns:
        - List[str]: A list of predicted reaction SMILES strings, standardized and deduplicated.
        """
        if print_results:
            ensure_directory_exists("./out")

        if isinstance(input_smiles, str):
            split_smiles = input_smiles.split(">>")
            smiles_list = (
                split_smiles[0].split(".")
                if prediction_type == "forward"
                else split_smiles[-1].split(".")
            )
        elif isinstance(input_smiles, list):
            smiles_list = input_smiles
        else:
            raise ValueError(
                "Unsupported input_smiles type. Expected a string or a list of strings."
            )

        predictions = []
        if rule_file_path.endswith(".gml"):
            # Single rule file specified directly
            reactions = RuleEngine.perform_reaction(
                rule_file_path=rule_file_path,
                initial_smiles=smiles_list,
                repeat_times=repeat_times,
                prediction_type=prediction_type,
                print_results=print_results,
            )
            predictions.extend(reactions)
        else:
            # Assume it's a directory containing multiple .gml files
            for rule_file in glob.glob(f"{rule_file_path}/*.gml"):
                reactions = RuleEngine.perform_reaction(
                    rule_file_path=rule_file,
                    initial_smiles=smiles_list,
                    repeat_times=repeat_times,
                    prediction_type=prediction_type,
                    print_results=print_results,
                )
                predictions.extend(reactions)

        # Standardize and deduplicate the predictions
        predictions = [standardize_rsmi(value) for value in predictions]
        if print_results:
            run_shell_command(command="mod_post")

        return list(set(predictions))

    @staticmethod
    def reaction_database_prediction(
        database: Union[pd.DataFrame, List[Dict]],
        rule_file_path: str,
        original_rsmi_col: str = "reactions",
        prediction_type: str = "forward",
        repeat_times: int = 1,
        print_results: bool = False,
    ) -> Union[pd.DataFrame, List[Dict]]:
        """
        Applies reaction predictions across a database of entries using specified rule files, accommodating both forward
        and backward directions. This method extracts SMILES strings from each database entry, applies reaction rules, and
        aggregates the results back into the database.

        Parameters:
        - database (Union[pd.DataFrame, List[Dict]]): The database to process, which can be either a Pandas DataFrame or a list of dictionaries.
        - rule_file_path (str): The file path to a GML file or a directory containing multiple GML files for reaction prediction.
        - original_rsmi_col (str, optional): The key in the database entries that identifies where the original reaction SMILES strings are stored. Defaults to 'reactions'.
        - prediction_type (str, optional): Specifies the direction of the reaction, either 'forward' or 'backward'. Defaults to 'forward'.
        - repeat_times (int, optional): Specifies the number of iterations for the reaction rule application. Defaults to 1.
        - print_results (bool): Print results in latex or not. Defaults to False.

        Returns:
        - Union[pd.DataFrame, List[Dict]]: The updated database with each entry containing the results of the reaction predictions.
        """
        if isinstance(database, pd.DataFrame):
            database_copy = database.copy().to_dict(orient="records")
        elif isinstance(database, list):
            database_copy = copy.deepcopy(database)
        else:
            raise ValueError(
                "Unsupported database type. Expected a Pandas DataFrame or a list of dictionaries."
            )

        for entry in database_copy:
            input_smiles = entry.get(original_rsmi_col, "")
            predictions = RuleExecutor.reaction_prediction(
                input_smiles,
                rule_file_path,
                prediction_type,
                repeat_times,
                print_results,
            )
            entry["predictions"] = predictions
            entry["number_predictions"] = len(predictions)

        if isinstance(database, pd.DataFrame):
            return pd.DataFrame(database_copy)
        return database_copy
