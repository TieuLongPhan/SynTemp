from rdkit import Chem
from rdkit.Chem.MolStandardize import normalize, tautomer, charge
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import List, Optional
import re

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

def assign_stereochemistry(mol: Chem.Mol, cleanIt: bool = True, force: bool = True) -> None:
    """
    Assigns stereochemistry to a molecule using RDKit's AssignStereochemistry.

    Args:
        mol: RDKit Mol object.
        cleanIt: Flag indicating whether to clean the molecule. Default is True.
        force: Flag indicating whether to force stereochemistry assignment. Default is True.

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
    1. Remove any atom mapping numbers denoted by ':' followed by one or more digits.
    2. Simplify the atomic notation by removing square brackets around atoms that do not need them.
    
    Parameters:
    - smiles (str): The SMILES string to be processed.
    
    Returns:
    - str: The processed SMILES string with atom mappings removed and simplified atomic notations.
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
    Converts a SMILES string to an RDKit Mol object, with error handling for invalid strings.
    
    Parameters:
    - smiles (str): The SMILES string to be converted.
    
    Returns:
    - Chem.Mol: An RDKit Mol object created from the given SMILES string. None if conversion fails.
    
    Raises:
    - ValueError: If the SMILES string is invalid and cannot be converted to a Mol object.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return mol

def filter_valid_molecules(smiles_list: List[str]) -> List[Chem.Mol]:
    """
    Filters a list of SMILES strings, converting them to RDKit Mol objects, while ignoring invalid or empty molecules.
    
    Parameters:
    - smiles_list (List[str]): A list of SMILES strings to be processed.
    
    Returns:
    - List[Chem.Mol]: A list of RDKit Mol objects derived from valid, non-empty SMILES strings in the input list.
    """
    valid_molecules = []
    for smiles in smiles_list:
        try:
            mol = mol_from_smiles(smiles)
            if mol.GetNumAtoms() > 0:
                valid_molecules.append(mol)
        except ValueError:
            continue
    return valid_molecules

def standardize_rsmi(rsmi: str) -> str:
    """
    Standardizes a reaction SMILES (rSMI) by ensuring that all reactants and products are valid molecules with atoms.
    
    The function splits the reaction into reactants and products, filters and validates them, and then
    assembles them back into a standardized reaction SMILES string.
    
    Parameters:
    - rsmi (str): The reaction SMILES string to be standardized.
    
    Returns:
    - str: The standardized reaction SMILES string with valid and non-empty reactants and products.
    """
    reactants, products = rsmi.split('>>')
    reactant_molecules = filter_valid_molecules(reactants.split('.'))
    product_molecules = filter_valid_molecules(products.split('.'))

    # Convert molecules back to SMILES and assemble the standardized reaction SMILES string
    standardized_reactants = '.'.join(Chem.MolToSmiles(mol) for mol in reactant_molecules)
    standardized_products = '.'.join(Chem.MolToSmiles(mol) for mol in product_molecules)

    return f"{standardized_reactants}>>{standardized_products}"
