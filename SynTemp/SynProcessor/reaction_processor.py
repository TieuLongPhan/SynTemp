from typing import List, Dict, Union
from rdkit import Chem
from joblib import Parallel, delayed


class ReactionProcessor:
    """
    A class for processing chemical reactions, including labeling and canonicalization
    of SMILES strings using RDKit.
    """

    @staticmethod
    def label_reactions(reaction_dict: Dict) -> Dict:
        """
        Labels chemical reactions based on their reactants, indicating whether they
        are oxidation or reduction reactions, and canonicalizes the SMILES strings.

        Parameters:
        - reaction_list (List[Dict]): A list of dictionaries, each representing a reaction
          with keys 'R-id' and 'new_reaction'.

        Returns:
        - List[Dict]: A list of dictionaries, each augmented with a 'label', 'reactants',
          and 'products' keys, where 'reactants' and 'products' are canonicalized SMILES.
        """

        label = "unspecified"
        r_id = reaction_dict.get("R-id", "N/A")
        new_reaction = reaction_dict.get("new_reaction", "")

        try:
            reactants, products = new_reaction.split(">>", 1)
        except ValueError:
            reactants, products = "", ""

        labeling_criteria = {
            ".[O]": "Oxidation",
            ".[H]": "Reduction",
        }

        for marker, reaction_label in labeling_criteria.items():
            if marker in reactants:
                label = reaction_label
                break

        reactants_smiles = Chem.CanonSmiles(reactants) if reactants else ""
        products_smiles = Chem.CanonSmiles(products) if products else ""

        new_dict = {
            "R-id": r_id,
            "new_reaction": new_reaction,
            "label": label,
            "reactants": reactants_smiles,
            "products": products_smiles,
        }

        return new_dict

    @classmethod
    def process_reactions_parallel(
        cls, reaction_list: List[Dict], n_jobs: int = -1
    ) -> List[Dict]:
        """
        Processes a list of chemical reactions in parallel, labeling them and
        canonicalizing their SMILES strings.

        Parameters:
        - reaction_list (List[Dict]): A list of dictionaries, each representing a reaction.
        - n_jobs (int): The number of jobs to run in parallel. -1 means using all processors.

        Returns:
        - List[Dict]: A processed list of reactions with labels and canonicalized SMILES.
        """
        results = Parallel(n_jobs=n_jobs)(
            delayed(cls.label_reactions)(reaction) for reaction in reaction_list
        )

        return results

    @staticmethod
    def calculate_charge(smiles: str) -> int:
        """
        Calculates the formal charge of a given molecule represented by a SMILES string.

        Parameters:
        - smiles (str): A SMILES string representing a molecule.

        Returns:
        - int: The formal charge of the molecule.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0  # Return 0 if the molecule cannot be parsed
        return Chem.rdmolops.GetFormalCharge(mol)

    @staticmethod
    def process_reaction(reaction: Dict[str, str]) -> Dict[str, Union[str, int]]:
        """
        Calculates and adds the total charge of products in a single reaction.

        Parameters:
        - reaction (Dict[str, str]): A dictionary representing a reaction with keys
          'R-id' and 'new_reaction'.

        Returns:
        - Dict[str, Union[str, int]]: The same reaction dictionary, with an added key
          'total_charge_in_products' indicating the sum of formal charges in its products.
        """
        products = reaction.get("products", "").split(".")
        total_charge = sum(
            ReactionProcessor.calculate_charge(product) for product in products
        )
        reaction["total_charge_in_products"] = total_charge
        return reaction

    @classmethod
    def sum_of_charge_in_products(
        cls, reaction_list: List[Dict], n_jobs: int = -1
    ) -> List[Dict]:
        """
        Calculates the sum of formal charges in the products of each reaction in the
        given list and adds this information to the reaction dictionaries.

        Parameters:
        - reaction_list (List[Dict]): A list of dictionaries, each representing a reaction.
        - n_jobs (int): The number of jobs to run in parallel. -1 means using all processors.

        Returns:
        - List[Dict]: The same list of reactions, with an added key 'total_charge_in_products'
          for each reaction, indicating the sum of formal charges in its products.
        """

        # Parallel processing of reactions to calculate the sum of charges in products
        processed_list = Parallel(n_jobs=n_jobs)(
            delayed(cls.process_reaction)(reaction) for reaction in reaction_list
        )

        return processed_list
