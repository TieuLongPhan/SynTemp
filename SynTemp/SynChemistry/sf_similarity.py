from typing import Callable, List, Tuple, Dict
import copy
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Avalon import pyAvalonTools as fpAvalon
from rdkit import DataStructs
from SynTemp.SynUtils.chemutils import mol_from_smiles


class SFSimilarity:

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
        reactant_mol = mol_from_smiles(reactants)
        product_mol = mol_from_smiles(products)

        total_similarity = 0
        for fptype in fingerprint_types:
            reactant_fp = cls.calculate_fingerprint(reactant_mol, fptype)
            product_fp = cls.calculate_fingerprint(product_mol, fptype)
            total_similarity += cls.calculate_tanimoto_similarity(
                reactant_fp, product_fp
            )

        return total_similarity

    def sort_reactions(
        self,
        reactions: List[str],
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
            total_similarity = self.process_reaction_smiles(
                rsmi, self.fingerprint_types
            )
            results.append((rsmi, total_similarity))
        results.sort(key=lambda x: x[1], reverse=True)  # Sort by similarity descending
        rsmi_list, _ = zip(*results) if results else ([], [])
        return rsmi_list
