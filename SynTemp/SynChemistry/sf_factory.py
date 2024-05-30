from typing import List, Dict, Any
import copy


class SFFactory:
    """
    A factory class for processing chemical reaction data with a given scoring function.
    The scoring function is expected to rank reactions based on some criteria, typically
    calculated from molecular fingerprints or other chemical informatics methods.
    """

    def __init__(self, scoring_function: Any):
        """
        Initializes the SFFactory with a specific scoring function used for ranking reactions.

        Parameters:
        - scoring_function (Any): An object that has a method `sort_reactions` which
                                  takes a list of reaction SMILES strings and returns
                                  a sorted list based on some criteria.
        """
        self.scoring_function = scoring_function

    def process_list_of_dicts(
        self, database: List[Dict[str, List[str]]], col_name: str
    ) -> List[Dict[str, List[str]]]:
        """
        Processes a list of dictionaries containing reaction data. Each dictionary is expected
        to contain a list of reaction SMILES strings identified by a specific column name.
        The method ranks these reactions using the scoring function and updates each dictionary
        to include this ranked list.

        Parameters:
        - database (List[Dict[str, List[str]]]): A list of dictionaries, each containing reaction data.
        - col_name (str): The key in each dictionary that indicates where the reaction SMILES are stored.

        Returns:
        - List[Dict[str, List[str]]]: The updated list of dictionaries with reactions replaced by
          a ranked list according to the scoring function.
        """
        list_of_dicts = copy.deepcopy(database)
        for item in list_of_dicts:
            reactions = item.get(col_name, [])
            if reactions:
                sorted_reactions = self.scoring_function.sort_reactions(reactions)
                item["ranked_reactions"] = sorted_reactions
                item.pop(
                    col_name, None
                )  # Optionally remove the original reactions list
            else:
                item["ranked_reactions"] = []

        return list_of_dicts
