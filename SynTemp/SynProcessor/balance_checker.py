from rdkit import Chem
from collections import defaultdict
from joblib import Parallel, delayed
from typing import List, Dict, Union, Tuple
from rdkit.Chem import rdMolDescriptors


class BalanceReactionCheck:
    """
    A class to check the balance of chemical reactions given in SMILES format.
    It supports parallel execution and maintains the input format in the output.
    """

    def __init__(
        self,
        input_data: Union[str, List[Union[str, Dict[str, str]]]],
        rsmi_column: str = "reactions",
        n_jobs: int = -1,
        verbose: int = 0,
    ):
        """
        Initializes the class with given input data, the column name for reactions in the input,
        number of jobs for parallel processing, and verbosity level.

        Parameters:
        - input_data (Union[str, List[Union[str, Dict[str, str]]]]): A single SMILES string,
          a list of SMILES strings, or a list of dictionaries with 'reactions' keys.
        - rsmi_column (str): The key/column name for reaction SMILES strings in the input data.
        - n_jobs (int): The number of parallel jobs to run for balance checking (default: -1, using all processors).
        - verbose (int): The verbosity level of joblib parallel execution (default: 0).
        """
        self.rsmi_column = rsmi_column
        self.reactions = self.parse_input(input_data)
        self.n_jobs = n_jobs
        self.verbose = verbose

    def parse_input(
        self, input_data: Union[str, List[Union[str, Dict[str, str]]]]
    ) -> List[Dict[str, str]]:
        """
        Parses the input data into a standardized list containing dictionaries for each reaction.

        Parameters:
        - input_data (Union[str, List[Union[str, Dict[str, str]]]]): The input data to be processed.

        Returns:
        - List[Dict[str, str]]: A list of dictionaries with reaction SMILES strings.
        """
        standardized_input = []
        if isinstance(input_data, str):
            standardized_input.append({self.rsmi_column: input_data})
        elif isinstance(input_data, list):
            for item in input_data:
                if isinstance(item, str):
                    standardized_input.append({self.rsmi_column: item})
                elif isinstance(item, dict) and self.rsmi_column in item:
                    standardized_input.append(item)
        else:
            raise ValueError("Unsupported input type")
        return standardized_input

    @staticmethod
    def parse_reaction(reaction_smiles: str) -> Tuple[List[str], List[str]]:
        """
        Splits a reaction SMILES string into reactants and products.

        Parameters:
        - reaction_smiles (str): A SMILES string representing a chemical reaction.

        Returns:
        - Tuple[List[str], List[str]]: Lists of SMILES strings for reactants and products.
        """
        reactants_smiles, products_smiles = reaction_smiles.split(">>")
        return reactants_smiles, products_smiles

    @staticmethod
    def mol_to_molecular_formula(mol):
        """
        Converts an RDKit molecule object to its molecular formula.

        Args:
        - mol (rdkit.Chem.Mol): An RDKit molecule object.

        Returns:
        - str: The molecular formula of the molecule.
        """
        return rdMolDescriptors.CalcMolFormula(mol)

    @staticmethod
    def is_balanced(
        reaction_dict: Dict[str, str], rsmi_column: str
    ) -> Dict[str, Union[bool, str]]:
        """
        Checks if a single reaction (in SMILES format) is balanced, maintaining the input format.

        Parameters:
        - reaction_dict (Dict[str, str]): A dictionary containing the reaction SMILES string.

        Returns:
        - Dict[str, Union[bool, str]]: A dictionary indicating if the reaction is balanced,
          along with the original reaction data.
        """
        reaction_smiles = reaction_dict[rsmi_column]
        reactants_smiles, products_smiles = BalanceReactionCheck.parse_reaction(
            reaction_smiles
        )
        reactants_mols = Chem.MolFromSmiles(reactants_smiles)
        products_mols = Chem.MolFromSmiles(products_smiles)
        reactants_forumula = BalanceReactionCheck.mol_to_molecular_formula(
            reactants_mols
        )
        products_forumula = BalanceReactionCheck.mol_to_molecular_formula(products_mols)
        if reactants_forumula != products_forumula:
            return {"balanced": False, **reaction_dict}

        return {"balanced": True, **reaction_dict}

    def check_balances(
        self,
    ) -> Tuple[List[Dict[str, Union[bool, str]]], List[Dict[str, Union[bool, str]]]]:
        """
        Checks the balance of all reactions in the input data.

        Returns:
        - Tuple[List[Dict[str, Union[bool, str]]], List[Dict[str, Union[bool, str]]]]: Two lists containing dictionaries
          of balanced and unbalanced reactions, respectively.
        """
        results = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
            delayed(self.is_balanced)(reaction, self.rsmi_column)
            for reaction in self.reactions
        )

        balanced_reactions = [reaction for reaction in results if reaction["balanced"]]
        unbalanced_reactions = [
            reaction for reaction in results if not reaction["balanced"]
        ]

        return balanced_reactions, unbalanced_reactions
