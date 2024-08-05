import unittest
from rdkit import Chem
from pathlib import Path
from syntemp.SynUtils.chemutils import (
    normalize_molecule,
    canonicalize_tautomer,
    salts_remover,
    reionize_charges,
    uncharge_molecule,
    assign_stereochemistry,
    fragments_remover,
    remove_hydrogens_and_sanitize,
)

root_dir = Path(__file__).parents[2]


class TestMoleculeProcessing(unittest.TestCase):
    def setUp(self):
        # Example molecule: caffeine
        self.caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        self.mol = Chem.MolFromSmiles(self.caffeine_smiles)

    def test_normalize_molecule(self):
        normalized_mol = normalize_molecule(self.mol)
        self.assertIsInstance(normalized_mol, Chem.Mol)

    def test_canonicalize_tautomer(self):
        canonical_tautomer = canonicalize_tautomer(self.mol)
        self.assertIsInstance(canonical_tautomer, Chem.Mol)

    def test_salts_remover(self):
        # Using a molecule with a known salt
        mol_with_salt = Chem.MolFromSmiles("CC(=O)O.[Na+]")
        salt_removed_mol = salts_remover(mol_with_salt)
        self.assertNotEqual(
            Chem.MolToSmiles(salt_removed_mol), Chem.MolToSmiles(mol_with_salt)
        )

    def test_reionize_charges(self):
        reionized_mol = reionize_charges(self.mol)
        self.assertIsInstance(reionized_mol, Chem.Mol)

    def test_uncharge_molecule(self):
        charged_mol = Chem.AddHs(Chem.MolFromSmiles("[NH4+].[Cl-]"))
        uncharged_mol = uncharge_molecule(charged_mol)
        self.assertEqual(Chem.rdmolops.GetFormalCharge(uncharged_mol), 0)

    def test_assign_stereochemistry(self):
        mol_with_stereo = Chem.MolFromSmiles("C[C@H](O)[C@@H](O)C")
        assign_stereochemistry(mol_with_stereo)
        self.assertEqual(
            len(Chem.FindMolChiralCenters(mol_with_stereo, includeUnassigned=True)), 2
        )

    def test_fragments_remover(self):
        mol_with_fragments = Chem.MolFromSmiles("CCO.OCC")
        largest_fragment = fragments_remover(mol_with_fragments)
        self.assertEqual(
            largest_fragment.GetNumAtoms(), 3
        )  # Expecting the ethyl alcohol fragment

    def test_remove_hydrogens_and_sanitize(self):
        mol_with_h = Chem.AddHs(self.mol)
        cleaned_mol = remove_hydrogens_and_sanitize(mol_with_h)
        self.assertEqual(cleaned_mol.GetNumAtoms(), self.mol.GetNumAtoms())


if __name__ == "__main__":
    unittest.main()
