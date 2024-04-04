import unittest
import sys
from pathlib import Path
root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynITSG.SynProcessor.uncharge_reaction import UnchargeReaction  # Replace 'your_module' with the actual module name
from rdkit import Chem

class TestUnchargeReaction(unittest.TestCase):
    def test_mol_to_molecular_formula(self):
        mol = Chem.MolFromSmiles("CCO")
        formula = UnchargeReaction.mol_to_molecular_formula(mol)
        self.assertEqual(formula, "C2H6O")

    def test_uncharge_anion(self):
        smiles = "[OH-]"
        uncharged_smiles = UnchargeReaction.uncharge_anion(smiles)
        self.assertEqual(uncharged_smiles, "O")

    def test_uncharge_cation(self):
        smiles = "[NH4+]"
        uncharged_smiles = UnchargeReaction.uncharge_cation(smiles)
        self.assertEqual(uncharged_smiles, "[NH4]") # can not uncharged N+

    def test_ammonia_hydroxide_standardize(self):
        reaction_smiles = "[NH4+].[OH-]>>N.O"
        standardized_smiles = UnchargeReaction.ammonia_hydroxide_standardize(reaction_smiles)
        self.assertTrue("N.O" in standardized_smiles or "O.N" in standardized_smiles)

    def test_uncharge_smiles(self):
        charge_smiles = "[NH4+].[OH-]"
        uncharged_smiles = UnchargeReaction.uncharge_smiles(charge_smiles)
        self.assertTrue("N" in uncharged_smiles and "O" in uncharged_smiles)

    def test_apply_uncharge_smiles_to_reactions(self):
        reactions = [
            {"reactants": "[NH4+].[OH-]", "products": "N.O"}
        ]
        uncharged_reactions = UnchargeReaction.apply_uncharge_smiles_to_reactions(
            reactions, UnchargeReaction.uncharge_smiles
        )
        for reaction in uncharged_reactions:
            self.assertTrue("N" in reaction["new_reactants"] and "O" in reaction["new_reactants"])
            self.assertTrue(reaction["success"])

if __name__ == "__main__":
    unittest.main()
