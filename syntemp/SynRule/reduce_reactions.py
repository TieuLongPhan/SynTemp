import copy
from typing import List, Dict, Tuple
from syntemp.SynUtils.chemutils import get_combined_molecular_formula


class ReduceReactions:

    @staticmethod
    def check_balance(reactants: str, products: str) -> bool:
        """
        Checks if the reaction is balanced based on the molecular formulas of reactants
        and products.

        Parameters:
        - reactants (str): SMILES string for the reactants.
        - products (str): SMILES string for the products.

        Returns:
        - bool: True if the reaction is balanced, otherwise False.
        """
        reactant_formula = get_combined_molecular_formula(reactants)
        product_formula = get_combined_molecular_formula(products)
        return reactant_formula == product_formula

    @staticmethod
    def process_rsmi(reaction_smiles: str) -> Tuple[str, bool]:
        """
        Processes a reaction SMILES string to check if the reaction is balanced.

        Parameters:
        - reaction_smiles (str): The reaction SMILES string formatted as
        'reactants>>products'.

        Returns:
        - Tuple[str, bool]: A tuple containing the reaction SMILES and a boolean
        indicating if it is balanced.
        """
        reactants, products = reaction_smiles.split(">>")
        balanced = ReduceReactions.check_balance(reactants, products)
        return balanced

    @staticmethod
    def process_list_of_rsmi(reactions: List[str]) -> List[Tuple[str, bool]]:
        """
        Processes a list of reaction SMILES strings and filters out those that are
        balanced.

        Parameters:
        - reactions (List[str]): A list of reaction SMILES strings.

        Returns:
        - List[Tuple[str, bool]]: A list of tuples, each containing a reaction SMILES
        string and a boolean for balance.
        """
        return [rsmi for rsmi in reactions if ReduceReactions.process_rsmi(rsmi)]

    @classmethod
    def process_list_of_dicts(cls, database: List[Dict], col_name: str) -> List[Dict]:
        """
        Processes a list of dictionaries, each containing a list of reaction SMILES
        strings, and updates them with balanced reactions.

        Parameters:
        - database (List[Dict]): The list of dictionaries containing reaction data.
        - col_name (str): The key in each dictionary that holds the reaction SMILES
        strings.

        Returns:
        - List[Dict]: The updated list of dictionaries with balanced reactions.
        """
        list_of_dicts = copy.deepcopy(database)
        for item in list_of_dicts:
            balanced_reactions = cls.process_list_of_rsmi(item[col_name])
            item[col_name] = balanced_reactions
        return list_of_dicts
