from typing import Callable, List, Tuple, Union, Dict, Any
import pandas as pd
import copy
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Avalon import pyAvalonTools as fpAvalon
from rdkit import DataStructs
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


class SimilarityRanking:
    def __init__(self, fingerprint_types: List[str]):
        """
        Initialize the SimilarityRanking class with specific fingerprint types used for molecular analysis.

        Parameters:
        - fingerprint_types (List[str]): A list of fingerprint types such as ['ECFP4', 'RDK5'].
        """
        self.fingerprint_types = fingerprint_types

    @classmethod
    def parse_fingerprint_settings(cls, fptype: str) -> Callable:
        """
        Parses the fingerprint type and returns a callable function for generating the specified molecular fingerprint.

        Parameters:
        - fptype (str): The fingerprint type, e.g., 'ECFP4', 'RDK5', 'MACCS', or 'Avalon'.

        Returns:
        - Callable: A function that takes an RDKit molecule and returns its fingerprint.
        """
        if fptype.startswith("ECFP"):
            radius = int(fptype[4]) // 2
            nBits = 2 ** (10 + radius)
            return lambda mol: AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)
        if fptype.startswith("FCFP"):
            radius = int(fptype[4]) // 2
            nBits = 2 ** (10 + radius)
            return lambda mol: AllChem.GetMorganFingerprintAsBitVect(
                mol, radius, nBits, useFeatures=True
            )
        elif fptype.startswith("RDK"):
            maxPath = int(fptype[3])
            fpSize = 2048 if maxPath < 7 else 4096
            return lambda mol: Chem.RDKFingerprint(mol, maxPath=maxPath, fpSize=fpSize)
        elif fptype == "MACCS":
            return lambda mol: MACCSkeys.GenMACCSKeys(mol)
        elif fptype == "Avalon":
            return lambda mol: fpAvalon.GetAvalonFP(mol, 1024)
        else:
            raise ValueError(f"Unsupported fingerprint type: {fptype}")

    @classmethod
    def calculate_fingerprint(
        cls, mol: Chem.Mol, fingerprint_type: str
    ) -> DataStructs.cDataStructs.ExplicitBitVect:
        """
        Generates a fingerprint for a given molecule based on the specified fingerprint type.

        Parameters:
        - mol (Chem.Mol): The molecule for which to generate the fingerprint.
        - fingerprint_type (str): The type of fingerprint to generate.

        Returns:
        - DataStructs.cDataStructs.ExplicitBitVect: The generated fingerprint.
        """
        fingerprint_func = cls.parse_fingerprint_settings(fingerprint_type)
        return fingerprint_func(mol)

    @staticmethod
    def calculate_tanimoto_similarity(
        fp1: DataStructs.cDataStructs.ExplicitBitVect,
        fp2: DataStructs.cDataStructs.ExplicitBitVect,
    ) -> float:
        """
        Calculates the Tanimoto similarity between two molecular fingerprints.

        Parameters:
        - fp1 (DataStructs.cDataStructs.ExplicitBitVect): The first fingerprint.
        - fp2 (DataStructs.cDataStructs.ExplicitBitVect): The second fingerprint.

        Returns:
        - float: The Tanimoto similarity score.
        """
        return DataStructs.TanimotoSimilarity(fp1, fp2)

    @staticmethod
    def smiles_to_mol(smiles: str) -> Union[Chem.Mol, None]:
        """
        Converts a SMILES string to an RDKit molecule object.

        Parameters:
        - smiles (str): The SMILES string to convert.

        Returns:
        - Chem.Mol: The RDKit molecule object, or None if the conversion fails.
        """
        return Chem.MolFromSmiles(smiles)

    @staticmethod
    def get_combined_molecular_formula(smiles: str) -> str:
        """
        Computes the molecular formula for a molecule represented by a SMILES string.

        Parameters:
        - smiles (str): The SMILES string of the molecule.

        Returns:
        - str: The molecular formula, or an empty string if the molecule is invalid.
        """
        mol = SimilarityRanking.smiles_to_mol(smiles)
        if not mol:
            return ""
        return CalcMolFormula(mol)

    @classmethod
    def check_balance(cls, reactants: str, products: str) -> bool:
        """
        Checks if the reaction is balanced based on the molecular formulas of reactants and products.

        Parameters:
        - reactants (str): SMILES string for the reactants.
        - products (str): SMILES string for the products.

        Returns:
        - bool: True if the reaction is balanced, otherwise False.
        """
        reactant_formula = cls.get_combined_molecular_formula(reactants)
        product_formula = cls.get_combined_molecular_formula(products)
        return reactant_formula == product_formula

    @classmethod
    def process_reaction_smiles(
        cls, reaction_smiles: str, fingerprint_types
    ) -> Tuple[float, bool]:
        """
        Processes a reaction SMILES string to calculate the summed Tanimoto similarity and checks if the reaction is balanced.

        Parameters:
        - reaction_smiles (str): The reaction SMILES string formatted as 'reactants>>products'.
        - fingerprint_types (List[str]): The list of fingerprint types to apply.

        Returns:
        - Tuple[float, bool]: The summed Tanimoto similarity and a boolean indicating if the reaction is balanced.
        """
        reactants, products = reaction_smiles.split(">>")
        reactant_mol = cls.smiles_to_mol(reactants)
        product_mol = cls.smiles_to_mol(products)

        total_similarity = 0
        for fptype in fingerprint_types:
            reactant_fp = cls.calculate_fingerprint(reactant_mol, fptype)
            product_fp = cls.calculate_fingerprint(product_mol, fptype)
            total_similarity += cls.calculate_tanimoto_similarity(
                reactant_fp, product_fp
            )

        balanced = cls.check_balance(reactants, products)
        return total_similarity, balanced

    @classmethod
    def fit(
        cls, reactions: List[str], fingerprint_types
    ) -> Tuple[List[str], List[float]]:
        """
        Processes a list of reaction SMILES strings, returning lists of reaction SMILES and their summed similarities if balanced, sorted by similarity.

        Parameters:
        - reactions (List[str]): A list of reaction SMILES strings.
        - fingerprint_types (List[str]): The list of fingerprint types to apply.

        Returns:
        - Tuple[List[str], List[float]]: Lists of reaction SMILES and their summed similarities, sorted by similarity in descending order.
        """
        results = []
        for rsmi in reactions:
            total_similarity, balanced = cls.process_reaction_smiles(
                rsmi, fingerprint_types
            )
            if balanced:
                results.append((rsmi, total_similarity))
        results.sort(key=lambda x: x[1], reverse=True)  # Sort by similarity descending
        rsmi_list, similarity_list = zip(*results) if results else ([], [])
        return list(rsmi_list), list(similarity_list)

    @classmethod
    def process_list_of_dicts(
        cls, database: List[Dict], col_name: str, fingerprint_types
    ) -> List[Dict]:
        """
        Processes a list of dictionaries containing reaction SMILES strings, returning the processed data
        with reactions ranked by their calculated Tanimoto similarity.

        Parameters:
        - database (List[Dict]): List of dictionaries containing the reaction data.
        - col_name (str): The key in each dictionary where reaction SMILES are stored.
        - fingerprint_types (List[str]): List of fingerprint types used for calculating the similarity.

        Returns:
        - List[Dict]: The list of dictionaries updated to include 'rank' and 'similarity' keys.
        """
        list_of_dicts = copy.deepcopy(database)
        for item in list_of_dicts:
            reactions = []
            for rsmi in item[col_name]:
                total_similarity, balanced = cls.process_reaction_smiles(
                    rsmi, fingerprint_types
                )
                if balanced:
                    reactions.append((rsmi, total_similarity))

            # Sort reactions by similarity in descending order
            reactions.sort(key=lambda x: x[1], reverse=True)

            # Prepare the output to include 'rank' and 'similarity' keys
            item["rank"] = [rsmi for rsmi, sim in reactions]
            item["similarity"] = [sim for rsmi, sim in reactions]

            # Optionally remove the original reactions list if no longer needed
            item.pop(col_name, None)

        return list_of_dicts
