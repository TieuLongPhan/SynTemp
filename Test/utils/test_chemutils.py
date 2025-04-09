import unittest
from rdkit import Chem
from pathlib import Path
from syntemp.utils.chemutils import (
    remove_hydrogens_and_sanitize,
)

root_dir = Path(__file__).parents[2]


class TestMoleculeProcessing(unittest.TestCase):
    def setUp(self):
        # Example molecule: caffeine
        self.caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
        self.mol = Chem.MolFromSmiles(self.caffeine_smiles)

    def test_remove_hydrogens_and_sanitize(self):
        mol_with_h = Chem.AddHs(self.mol)
        cleaned_mol = remove_hydrogens_and_sanitize(mol_with_h)
        self.assertEqual(cleaned_mol.GetNumAtoms(), self.mol.GetNumAtoms())


if __name__ == "__main__":
    unittest.main()
