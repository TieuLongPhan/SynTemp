import re
from rdkit import Chem
from rdkit.Chem.MolStandardize import normalize, tautomer, charge
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import rdmolops
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


def remove_atom_mapping(reaction_smiles):
    # Split the reaction SMILES into reactants and products
    parts = reaction_smiles.split(">>")
    if len(parts) != 2:
        raise ValueError("Invalid reaction SMILES format.")

    # Function to remove atom mappings from a SMILES string
    def clean_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to an RDKit mol object
        if mol is None:
            raise ValueError("Invalid SMILES string.")
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)  # Remove atom mapping
        return Chem.MolToSmiles(mol, True)  # Convert mol back to SMILES

    # Apply the cleaning function to both reactants and products
    reactants_clean = clean_smiles(parts[0])
    products_clean = clean_smiles(parts[1])

    # Combine the cleaned reactants and products back into a reaction SMILES
    return f"{reactants_clean}>>{products_clean}"


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


def remove_hydrogen(mol, atom_index):
    edited_mol = Chem.RWMol(mol)
    atom = edited_mol.GetAtomWithIdx(atom_index)
    # Ensure there is at least one explicit hydrogen to remove
    if atom.GetNumExplicitHs() > 0:
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
    return edited_mol.GetMol()


# Process a single component of the molecule
def process_component(component):
    for atom_index in range(component.GetNumAtoms()):
        try:
            test_mol = remove_hydrogen(component, atom_index)
            rdmolops.SanitizeMol(test_mol)
            print("Fuck")
            return test_mol
        except Exception:
            continue
    return None


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
        # try:
        mol = Chem.MolFromSmiles(smiles)
        # except:
        #     mol = Chem.MolFromSmiles(smiles, sanitize=False)
        #     mol = process_component(mol)
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
    if reactant_molecules and product_molecules:
        standardized_reactants = ".".join(
            sorted(
                Chem.MolToSmiles(mol, isomericSmiles=stereo)
                for mol in reactant_molecules
            )
        )
        standardized_products = ".".join(
            sorted(
                Chem.MolToSmiles(mol, isomericSmiles=stereo)
                for mol in product_molecules
            )
        )

        return f"{standardized_reactants}>>{standardized_products}"
    else:
        return None


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


def enumerate_tautomers(reaction_smiles: str) -> Optional[List[str]]:
    """
    Enumerates possible tautomers for reactants while canonicalizing the products in a
    reaction SMILES string. This function first splits the reaction SMILES string into
    reactants and products. It then generates all possible tautomers for the reactants and
    canonicalizes the product molecule. The function returns a list of reaction SMILES
    strings for each tautomer of the reactants combined with the canonical product.

    Parameters:
    - reaction_smiles (str): A SMILES string of the reaction formatted as
    'reactants>>products'.

    Returns:
    - List[str] | None: A list of SMILES strings for the reaction, with each string
    representing a different
    - tautomer of the reactants combined with the canonicalized products. Returns None if
    an error occurs or if invalid SMILES strings are provided.

    Raises:
    - ValueError: If the provided SMILES strings cannot be converted to molecule objects,
    indicating invalid input.
    """
    try:
        # Split the input reaction SMILES string into reactants and products
        reactants_smiles, products_smiles = reaction_smiles.split(">>")

        # Convert SMILES strings to molecule objects
        reactants_mol = Chem.MolFromSmiles(reactants_smiles)
        products_mol = Chem.MolFromSmiles(products_smiles)

        if reactants_mol is None or products_mol is None:
            raise ValueError(
                "Invalid SMILES string provided for reactants or products."
            )

        # Initialize tautomer enumerator

        enumerator = rdMolStandardize.TautomerEnumerator()

        # Enumerate tautomers for the reactants and canonicalize the products
        try:
            reactants_can = enumerator.Enumerate(reactants_mol)
        except Exception as e:
            print(f"An error occurred: {e}")
            reactants_can = [reactants_mol]
        products_can = products_mol

        # Convert molecule objects back to SMILES strings
        reactants_can_smiles = [Chem.MolToSmiles(i) for i in reactants_can]
        products_can_smiles = Chem.MolToSmiles(products_can)

        # Combine each reactant tautomer with the canonical product in SMILES format
        rsmi_list = [i + ">>" + products_can_smiles for i in reactants_can_smiles]
        if len(rsmi_list) == 0:
            return [reaction_smiles]
        else:
            # rsmi_list.remove(reaction_smiles)
            rsmi_list.insert(0, reaction_smiles)
            return rsmi_list

    except Exception as e:
        print(f"An error occurred: {e}")
        return [reaction_smiles]


def mapping_success_rate(list_mapping_data):
    """
    Calculate the success rate of entries containing atom mappings in a list of data
    strings.

    Parameters:
    - list_mapping_in_data (list of str): List containing strings to be searched for atom
    mappings.

    Returns:
    - float: The success rate of finding atom mappings in the list as a percentage.

    Raises:
    - ValueError: If the input list is empty.
    """
    atom_map_pattern = re.compile(r":\d+")
    if not list_mapping_data:
        raise ValueError("The input list is empty, cannot calculate success rate.")

    success = sum(
        1 for entry in list_mapping_data if re.search(atom_map_pattern, entry)
    )
    rate = 100 * (success / len(list_mapping_data))

    return round(rate, 2)


def generate_reaction_smiles(
    temp_results: List[str], base_smiles: str, is_forward: bool = True
) -> List[str]:
    """
    Constructs reaction SMILES strings from intermediate results using a base SMILES
    string, indicating whether the process is a forward or backward reaction. This
    function iterates over a list of intermediate SMILES strings, combines them with the
    base SMILES, and formats them into complete reaction SMILES strings.

    Parameters:
    - temp_results (List[str]): Intermediate SMILES strings resulting from partial
    reactions or combinations.
    - base_smiles (str): The SMILES string representing the starting point of the
    reaction, either as reactants or products, depending on the reaction direction.
    - is_forward (bool, optional): Flag to determine the direction of the reaction; 'True'
    for forward reactions where 'base_smiles' are reactants, and 'False' for backward
    reactions where 'base_smiles' are products. Defaults to True.

    Returns:
    - List[str]: A list of complete reaction SMILES strings, formatted according to the
    specified reaction direction.
    """
    results = []
    for comb in temp_results:
        if comb:
            joined_smiles = ".".join(comb)
            reaction_smiles = (
                f"{base_smiles}>>{joined_smiles}"
                if is_forward
                else f"{joined_smiles}>>{base_smiles}"
            )
            results.append(reaction_smiles)
    return results
