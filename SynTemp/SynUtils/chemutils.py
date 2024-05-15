import re
from rdkit import Chem
from rdkit.Chem.MolStandardize import normalize, tautomer, charge
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from typing import List, Optional, Tuple
from collections import Counter
from itertools import combinations_with_replacement


def normalize_molecule(mol: Chem.Mol) -> Chem.Mol:
    """
    Normalize a molecule using RDKit's Normalizer.

    Args:
        mol (Chem.Mol): RDKit Mol object to be normalized.

    Returns:
        Chem.Mol: Normalized RDKit Mol object.
    """
    return normalize.Normalizer().normalize(mol)


def canonicalize_tautomer(mol: Chem.Mol) -> Chem.Mol:
    """
    Canonicalize the tautomer of a molecule using RDKit's TautomerCanonicalizer.

    Args:
    - mol (Chem.Mol): RDKit Mol object.

    Returns:
    - Chem.Mol: Mol object with canonicalized tautomer.
    """
    return tautomer.TautomerCanonicalizer().canonicalize(mol)


def salts_remover(mol: Chem.Mol) -> Chem.Mol:
    """
    Remove salt fragments from a molecule using RDKit's SaltRemover.

    Args:
    - mol (Chem.Mol): RDKit Mol object.

    Returns:
    - Chem.Mol: Mol object with salts removed.
    """
    remover = SaltRemover()
    return remover.StripMol(mol)


def reionize_charges(mol: Chem.Mol) -> Chem.Mol:
    """
    Adjust molecule to its most likely ionic state using RDKit's Reionizer.

    Args:
    - mol: RDKit Mol object.

    Returns:
    - Mol object with reionized charges.
    """
    return charge.Reionizer().reionize(mol)


def uncharge_molecule(mol: Chem.Mol) -> Chem.Mol:
    """
    Neutralize a molecule by removing counter-ions using RDKit's Uncharger.

    Args:
        mol: RDKit Mol object.

    Returns:
        Neutralized Mol object.
    """
    uncharger = rdMolStandardize.Uncharger()
    return uncharger.uncharge(mol)


def assign_stereochemistry(
    mol: Chem.Mol, cleanIt: bool = True, force: bool = True
) -> None:
    """
    Assigns stereochemistry to a molecule using RDKit's AssignStereochemistry.

    Args:
        mol: RDKit Mol object.
        cleanIt: Flag indicating whether to clean the molecule.
                Default is True.
        force: Flag indicating
            whether to force stereochemistry assignment.
            Default is True.

    Returns:
        None
    """
    Chem.AssignStereochemistry(mol, cleanIt=cleanIt, force=force)


def fragments_remover(mol: Chem.Mol) -> Chem.Mol:
    """
    Remove small fragments from a molecule, keeping only the largest one.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        Chem.Mol: Mol object with small fragments removed.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, default=None, key=lambda m: m.GetNumAtoms())


def remove_hydrogens_and_sanitize(mol: Chem.Mol) -> Chem.Mol:
    """
    Remove explicit hydrogens and sanitize a molecule.

    Args:
        mol (Chem.Mol): RDKit Mol object.

    Returns:
        Chem.Mol: Mol object with explicit hydrogens removed and sanitized.
    """
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    return mol


def remove_atom_mapping(smiles: str) -> str:
    """
    Removes atom mapping numbers and simplifies atomic notation in a SMILES string.

    This function processes a SMILES string to:
    1. Remove any atom mapping numbers denoted by ':'
        followed by one or more digits.
    2. Simplify the atomic notation by removing square
        brackets around atoms that do not need them.

    Parameters:
    - smiles (str): The SMILES string to be processed.

    Returns:
    - str: The processed SMILES string with atom mappings
            removed and simplified atomic notations.
    """
    # Remove atom mapping numbers
    pattern = re.compile(r":\d+")
    smiles = pattern.sub("", smiles)
    # Simplify atomic notation by removing unnecessary square brackets
    pattern = re.compile(r"\[(?P<atom>(B|C|N|O|P|S|F|Cl|Br|I){1,2})(?:H\d?)?\]")
    smiles = pattern.sub(r"\g<atom>", smiles)
    return smiles


def mol_from_smiles(smiles: str) -> Optional[Chem.Mol]:
    """
    Converts a SMILES string to an RDKit Mol object,
    with error handling for invalid strings.

    Parameters:
    - smiles (str): The SMILES string to be converted.

    Returns:
    - Chem.Mol: An RDKit Mol object created from the given SMILES string.
                None if conversion fails.

    Raises:
    - ValueError: If the SMILES string is invalid and cannot be converted to a Mol object.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return mol


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


def standardize_rsmi(rsmi: str, stereo: bool = False) -> str:
    """
    Standardizes a reaction SMILES (rSMI) by ensuring that all reactants and products
    are valid molecules with atoms and that the SMILES strings within the reactants
    and products are in ascending order. Optionally ignores stereochemistry.

    Parameters:
    - rsmi (str): The reaction SMILES string to be standardized.
    - ignore_stereo (bool): If True, stereochemical information is
                            ignored in the SMILES representation.

    Returns:
    - str: The standardized reaction SMILES string with valid,
            non-empty, and sorted reactants and products.
    """
    reactants, products = rsmi.split(">>")
    reactant_molecules = filter_valid_molecules(reactants.split("."))
    product_molecules = filter_valid_molecules(products.split("."))

    standardized_reactants = ".".join(
        sorted(
            Chem.MolToSmiles(mol, isomericSmiles=stereo) for mol in reactant_molecules
        )
    )
    standardized_products = ".".join(
        sorted(
            Chem.MolToSmiles(mol, isomericSmiles=stereo) for mol in product_molecules
        )
    )

    return f"{standardized_reactants}>>{standardized_products}"


def categorize_reactions(
    reactions: List[str], target_reaction: str
) -> Tuple[List[str], List[str]]:
    """
    Sorts a list of reaction SMILES strings into two groups based on
    their match with a specified target reaction. The categorization process
    distinguishes between reactions that align with the target reaction
    and those that do not.

    Parameters:
    - reactions (List[str]): The array of reaction SMILES strings to be categorized.
    - target_reaction (str): The SMILES string of the target reaction
                                used as the benchmark for categorization.

    Returns:
    - Tuple[List[str], List[str]]: A pair of lists, where the first contains
                                reactions matching the target and the second
                                comprises non-matching reactions.
    """
    match, not_match = [], []
    target_reaction = standardize_rsmi(target_reaction, stereo=False)
    for reaction_smiles in reactions:
        if reaction_smiles == target_reaction:
            match.append(reaction_smiles)
        else:
            not_match.append(reaction_smiles)
    return match, not_match


def count_carbons(smiles_list: List[str]) -> int:
    """
    Counts the total number of carbon atoms across all molecules
    represented by a list of SMILES strings.

    Args:
        smiles_list (List[str]): A list of SMILES strings representing molecules.

    Returns:
        int: The total number of carbon atoms across all valid SMILES strings in the list.
    """
    total_carbon_count = 0
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:  # Valid SMILES string
            total_carbon_count += sum(
                1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C"
            )
    return total_carbon_count


def get_combined_molecular_formula(smiles: str) -> str:
    """
    Computes the molecular formula for a molecule represented by a SMILES string.

    Parameters:
    - smiles (str): The SMILES string of the molecule.

    Returns:
    - str: The molecular formula, or an empty string if the molecule is invalid.
    """
    mol = mol_from_smiles(smiles)
    if not mol:
        return ""
    return CalcMolFormula(mol)


@staticmethod
def generate_smiles_combinations(
    smiles_list: List[str], target_smiles: str, use_formula: bool = False
) -> List[List[str]]:
    """
    Generates all combinations of SMILES strings from the given list
    that match the target specified by a target SMILES string.

    Args:
        smiles_list (List[str]): A list of SMILES strings.
        target_smiles (str): The SMILES string of the target molecule.
        use_formula (bool): If True, uses the molecular formula for comparison.
                            If False, uses the total number of carbon atoms.

    Returns:
        List[List[str]]: A list of combinations, where each combination is a list of
        SMILES strings that together match the target specified by either the
        total number of carbon atoms or the molecular formula.
    """
    if use_formula:
        target_formula = get_combined_molecular_formula([target_smiles])
    else:
        target_carbons = count_carbons([target_smiles])

    valid_combinations = []
    smiles_counts = Counter(smiles_list)

    for r in range(1, sum(smiles_counts.values()) + 1):
        for combo in combinations_with_replacement(smiles_list, r):
            if Counter(combo) <= smiles_counts:
                if use_formula:
                    combo_formula = get_combined_molecular_formula(combo)
                    if combo_formula == target_formula:
                        valid_combinations.append(list(combo))
                else:  # When `use_formula` is False
                    total_carbons = count_carbons(combo)
                    if total_carbons == target_carbons:
                        valid_combinations.append(list(combo))

        return valid_combinations


def remove_stereochemistry_from_reaction_smiles(reaction_smiles: str) -> str:
    """
    Removes stereochemical information from a reaction SMILES string.

    Parameters:
    - reaction_smiles (str): A reaction SMILES string possibly
                    containing stereochemical information.

    Returns:
    - str: A reaction SMILES string with stereochemistry removal.
    """
    # Split the reaction SMILES into reactants, agents (optional), and products
    parts = reaction_smiles.split(">>")
    if len(parts) not in [2, 3]:
        raise ValueError("Invalid reaction SMILES format.")

    # Process each part to remove stereochemistry
    def process_part(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Chem.RemoveStereochemistry(mol)
            return Chem.MolToSmiles(mol)
        return smiles

    # Apply the processing to each part
    processed_parts = [process_part(part) for part in parts]

    # Rejoin the parts into a complete reaction SMILES string
    return ">>".join(processed_parts)
