import unittest
import sys
from pathlib import Path

root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynProcessor.balance_checker import BalanceReactionCheck


class TestBalanceReactionCheck(unittest.TestCase):
    def test_single_smiles_balanced(self):
        """Test a single balanced reaction in SMILES format."""
        smiles = "Clc1cnc2nc1Nc1ccc(OCCC3CCNCC3)c(c1)CCc1cncc(c1)N2.O=C=NCc1ccco1>>O=C(NCc1ccco1)N1CCC(CCOc2ccc3cc2CCc2cncc(c2)Nc2ncc(Cl)c(n2)N3)CC1"
        checker = BalanceReactionCheck(smiles)
        balanced, unbalanced = checker.check_balances()
        self.assertEqual(len(balanced), 1)
        self.assertEqual(len(unbalanced), 0)

    def test_single_smiles_unbalanced(self):
        """Test a single unbalanced reaction in SMILES format."""
        smiles = "CC(=O)O.CCO>>CC(=O)OCC"
        checker = BalanceReactionCheck(smiles)
        balanced, unbalanced = checker.check_balances()
        self.assertEqual(len(balanced), 0)
        self.assertEqual(len(unbalanced), 1)

    def test_list_of_smiles_mixed(self):
        """Test a list of mixed balanced and unbalanced reactions in SMILES format."""
        reactions = [
            "Clc1cnc2nc1Nc1ccc(OCCC3CCNCC3)c(c1)CCc1cncc(c1)N2.O=C=NCc1ccco1>>O=C(NCc1ccco1)N1CCC(CCOc2ccc3cc2CCc2cncc(c2)Nc2ncc(Cl)c(n2)N3)CC1",
            "CC(=O)O.CCO>>CC(=O)OCC",
        ]
        checker = BalanceReactionCheck(reactions)
        balanced, unbalanced = checker.check_balances()
        self.assertEqual(len(balanced), 1)
        self.assertEqual(len(unbalanced), 1)

    def test_list_of_dicts_mixed(self):
        """Test a list of dictionaries with mixed balanced and unbalanced reactions."""
        reactions = [
            {
                "R-id": "test_1",
                "reactions": "Clc1cnc2nc1Nc1ccc(OCCC3CCNCC3)c(c1)CCc1cncc(c1)N2.O=C=NCc1ccco1>>O=C(NCc1ccco1)N1CCC(CCOc2ccc3cc2CCc2cncc(c2)Nc2ncc(Cl)c(n2)N3)CC1",
            },
            {"R-id": "test_2", "reactions": "CC(=O)O.CCO>>CC(=O)OCC"},
        ]
        checker = BalanceReactionCheck(reactions, rsmi_column="reactions")
        balanced, unbalanced = checker.check_balances()
        self.assertEqual(len(balanced), 1)
        self.assertEqual(len(unbalanced), 1)
        self.assertEqual(balanced[0]["R-id"], "test_1")
        self.assertEqual(unbalanced[0]["R-id"], "test_2")


if __name__ == "__main__":
    unittest.main()
