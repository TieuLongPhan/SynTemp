import re
from typing import List, Union
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from collections import Counter
from itertools import combinations_with_replacement
from joblib import Parallel, delayed


class MØDPostprocess:
    """
    A class for post-processing molecule operations, specifically focusing on generating combinations of SMILES strings
    based on carbon count or molecular formula criteria.

    Methods:
    - count_carbons: Counts the total number of carbon atoms across all molecules represented by a list of SMILES strings.
    - get_combined_molecular_formula: Computes the combined molecular formula for a list of molecules represented by SMILES strings.
    - generate_smiles_combinations: Generates all combinations of SMILES strings from the given list that match the specified target criteria.
    """
    
    @staticmethod
    def filter_valid_molecules(smiles_list: List[str]) -> List[Chem.Mol]:
        """
        Filters and returns valid RDKit molecule objects from a list of SMILES strings.
        
        Parameters:
        - smiles_list (List[str]): A list of SMILES strings.
        
        Returns:
        - List[Chem.Mol]: A list of valid RDKit molecule objects.
        """
        valid_molecules = []
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:  # Check if the molecule is valid
                valid_molecules.append(mol)
        return valid_molecules

    @staticmethod
    def standardize_rsmi(rsmi: str) -> str:
        """
        Standardizes a reaction SMILES (rSMI) by ensuring that all reactants and products are valid molecules with atoms
        and that the SMILES strings within the reactants and products are in ascending order.
        
        The function splits the reaction into reactants and products, filters and validates them, sorts them in ascending order,
        and then assembles them back into a standardized reaction SMILES string.
        
        Parameters:
        - rsmi (str): The reaction SMILES string to be standardized.
        
        Returns:
        - str: The standardized reaction SMILES string with valid, non-empty, and sorted reactants and products.
        """
        reactants, products = rsmi.split('>>')
        reactant_molecules = MØDPostprocess.filter_valid_molecules(reactants.split('.'))
        product_molecules = MØDPostprocess.filter_valid_molecules(products.split('.'))

        # Convert molecules back to SMILES, sort them, and assemble the standardized reaction SMILES string
        standardized_reactants = '.'.join(sorted(Chem.MolToSmiles(mol) for mol in reactant_molecules))
        standardized_products = '.'.join(sorted(Chem.MolToSmiles(mol) for mol in product_molecules))

        return f"{standardized_reactants}>>{standardized_products}"

    @staticmethod
    def count_carbons(smiles_list: List[str]) -> int:
        """
        Counts the total number of carbon atoms across all molecules represented by a list of SMILES strings.

        Args:
            smiles_list (List[str]): A list of SMILES strings representing molecules.

        Returns:
            int: The total number of carbon atoms across all valid SMILES strings in the list.
        """
        total_carbon_count = 0
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:  # Valid SMILES string
                total_carbon_count += sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        return total_carbon_count

    @staticmethod
    def get_combined_molecular_formula(smiles_list: List[str]) -> str:
        """
        Computes the combined molecular formula for a list of molecules represented by SMILES strings.

        Args:
            smiles_list (List[str]): A list of SMILES strings representing separate molecules.

        Returns:
            str: The combined molecular formula. Returns an empty string if any SMILES string is invalid.
        """
        total_formula = {}
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:  # Invalid SMILES string
                return ''  
            formula = CalcMolFormula(mol)
            for element in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
                symbol, count = element
                count = int(count) if count else 1
                total_formula[symbol] = total_formula.get(symbol, 0) + count
        return ''.join(f"{symbol}{total_formula[symbol]}" for symbol in sorted(total_formula.keys()))

    @staticmethod
    def generate_smiles_combinations(smiles_list: List[str], target_smiles: str, use_formula: bool = False) -> List[List[str]]:
        """
        Generates all combinations of SMILES strings from the given list that match the target specified by a target SMILES string.

        Args:
            smiles_list (List[str]): A list of SMILES strings.
            target_smiles (str): The SMILES string of the target molecule.
            use_formula (bool): If True, uses the molecular formula for comparison. If False, uses the total number of carbon atoms.

        Returns:
            List[List[str]]: A list of combinations, where each combination is a list of SMILES strings that together match the target
                             specified by either the total number of carbon atoms or the molecular formula.
        """
        if use_formula:
            target_formula = MØDPostprocess.get_combined_molecular_formula([target_smiles])
        else:
            target_carbons = MØDPostprocess.count_carbons([target_smiles])

        valid_combinations = []
        smiles_counts = Counter(smiles_list)

        for r in range(1, sum(smiles_counts.values()) + 1):
            for combo in combinations_with_replacement(smiles_list, r):
                if Counter(combo) <= smiles_counts:
                    if use_formula:
                        combo_formula = MØDPostprocess.get_combined_molecular_formula(combo)
                        if combo_formula == target_formula:
                            valid_combinations.append(list(combo))
                    else:  # When `use_formula` is False
                        total_carbons = MØDPostprocess.count_carbons(combo)
                        if total_carbons == target_carbons:
                            valid_combinations.append(list(combo))

        return valid_combinations
