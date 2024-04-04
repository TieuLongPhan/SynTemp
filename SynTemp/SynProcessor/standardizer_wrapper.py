from rdkit import Chem
from rdkit.Chem.MolStandardize import normalize, tautomer, charge
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize


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
