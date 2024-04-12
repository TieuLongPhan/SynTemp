import unittest
from rdkit import Chem
import unittest
import sys
import pandas as pd
from pathlib import Path

root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynProcessor.smiles_standardizer import SMILESStandardizer


class TestSMILESStandardizer(unittest.TestCase):  ## add Test case
    def setUp(self):
        self.standardizer = SMILESStandardizer()
        self.test_smiles = (
            "C[N+](C)(C)C.[Cl-]"  # Example molecule: Tetramethylammonium chloride
        )

    def test_standardize_mol(self):
        mol = Chem.MolFromSmiles(self.test_smiles)
        standardized_mol = self.standardizer.standardize_mol(
            mol, remove_salts=True, uncharge=True
        )
        standardized_smiles = Chem.MolToSmiles(standardized_mol)
        self.assertNotIn("[Cl-]", standardized_smiles)
        self.assertNotIn("[NH+]", standardized_smiles)

    def test_standardize_smiles(self):
        standardized_smiles, _ = self.standardizer.standardize_smiles(
            self.test_smiles, remove_salts=True, uncharge=True
        )
        self.assertNotIn("[Cl-]", standardized_smiles)
        self.assertNotIn("[NH+]", standardized_smiles)

    def test_standardize_dict_smiles(self):
        data = [{"smiles": self.test_smiles, "id": 1}]
        standardized_data = self.standardizer.standardize_dict_smiles(
            data, keys=["smiles"], remove_salts=True, uncharge=True, parallel=False
        )
        for item in standardized_data:
            self.assertNotIn("[Cl-]", item["standardized_smiles"])
            self.assertNotIn("[NH+]", item["standardized_smiles"])

    def test_parallel_standardize_dict_smiles(self):
        data = pd.DataFrame(
            {"smiles": [self.test_smiles for _ in range(10)], "id": range(10)}
        )
        standardized_data = self.standardizer.standardize_dict_smiles(
            data,
            keys=["smiles"],
            remove_salts=True,
            uncharge=True,
            parallel=True,
            n_jobs=2,
        )
        for standardized_smiles in standardized_data["standardized_smiles"]:
            self.assertNotIn("[Cl-]", standardized_smiles)
            self.assertNotIn("[NH+]", standardized_smiles)


if __name__ == "__main__":
    unittest.main()
