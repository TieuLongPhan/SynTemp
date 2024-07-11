from rdkit import Chem
from joblib import Parallel, delayed
from typing import Dict, List, Tuple, Any


class AMMPostprocessor:
    """
    A class to validate Consensus Atom-Atom Mapping (AAM) represented in SMILES format.
    """

    def __init__(self):
        pass

    @staticmethod
    def extract_mappings_and_count_atoms(smiles: str) -> Tuple[List[int], int]:
        """
        Takes a SMILES string with atom mapping as input and returns a sorted
        list of atom mapping numbers based on their appearance in the molecule
        along with the total count of atoms.

        Parameters:
        - smiles (str): A SMILES string with atom mapping representing the molecule.

        Returns:
        - tuple: A tuple where the first element is a list of atom mapping numbers
                sorted based on their appearance, and the second element is the
                total count of atoms in the molecule.
        """
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return ("Invalid SMILES", 0)

        atom_map_numbers = []
        for atom in molecule.GetAtoms():
            atom_map_num = atom.GetAtomMapNum()
            if atom_map_num:
                atom_map_numbers.append(atom_map_num)

        sorted_atom_map_numbers = sorted(atom_map_numbers)
        total_atom_count = len(molecule.GetAtoms())

        return (sorted_atom_map_numbers, total_atom_count)

    @staticmethod
    def is_consistent_mapping(reaction_smiles: str) -> bool:
        """
        Determines if the atom mapping numbers and atom counts in the reactants
        and products of a given reaction SMILES are consistent.

        Parameters:
        - reaction_smiles (str): A reaction SMILES string with atom mapping.

        Returns:
        - bool: True if the mapping is consistent and the
                atom counts match, False otherwise.
        """
        reactants_smiles, products_smiles = reaction_smiles.split(">>")
        mapping_reactants, atom_count_reactants = (
            AMMPostprocessor.extract_mappings_and_count_atoms(reactants_smiles)
        )
        mapping_products, atom_count_products = (
            AMMPostprocessor.extract_mappings_and_count_atoms(products_smiles)
        )

        if mapping_reactants != mapping_products:
            return False
        return (
            len(mapping_reactants) == atom_count_reactants
            and len(mapping_products) == atom_count_products
        )

    @staticmethod
    def postprocess(
        mapped_smiles: Dict[str, str], mapper_names: List[str], threshold: int
    ) -> Dict[str, Any]:
        """
        Post-processes the mapped SMILES based on consistency checks and a threshold.

        Parameters:
        - mapped_smiles (Dict[str, str]): Dictionary of mapped SMILES strings.
        - mapper_names (List[str]): List of mapper names to process.
        - threshold (int): Threshold for determining validity.

        Returns:
        - Dict[str, Any]: Updated dictionary with validity status.
        """
        valid_count = sum(
            AMMPostprocessor.is_consistent_mapping(mapped_smiles[mapper_name])
            for mapper_name in mapper_names
        )
        mapped_smiles["Valid"] = valid_count == threshold
        return mapped_smiles

    @staticmethod
    def parallel_postprocess(
        mapped_smiles_list: List[Dict[str, str]],
        mapper_names: List[str],
        threshold: int,
        n_jobs: int = -1,
        verbose: int = 10,
    ) -> List[Dict[str, Any]]:
        """
        Processes a list of mapped SMILES strings in parallel.

        Parameters:
        - mapped_smiles_list (List[Dict[str, str]]): A list of dictionaries of
                                                        mapped SMILES strings.
        - mapper_names (List[str]): List of mapper names to process.
        - threshold (int): Threshold for determining validity.
        - n_jobs (int): Number of jobs to run in parallel.
                        Defaults to -1 (use all processors).
        - verbose (int): Verbosity level.

        Returns:
        - List[Dict[str, Any]]: List of dictionaries with updated validity status.
        """
        return Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(AMMPostprocessor.postprocess)(
                mapped_smiles, mapper_names, threshold
            )
            for mapped_smiles in mapped_smiles_list
        )
