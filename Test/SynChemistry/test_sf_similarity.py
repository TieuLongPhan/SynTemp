import unittest
from rdkit import Chem
from rdkit.Chem import AllChem
from SynTemp.SynChemistry.sf_similarity import SFSimilarity


class TestSimilarityRanking(unittest.TestCase):
    def setUp(self):
        # Setup can include creating instances of the class, etc.
        self.similarity_ranking = SFSimilarity(["ECFP4", "MACCS"])
        self.data = [
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2ccc(NC(=O)C(c3ccc4c(c3NC(=O)C=C(C)C)CCN(C(=O)OCC)CC4)C(C)C)cc2CC1",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2cc3c(cc2CC1)C(C(C)C)C(=O)N3",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2ccc(NC(=O)CC(C)(C)c3ccc4c(c3NC(=O)C=C(C)C)CCN(C(=O)OCC)CC4)cc2CC1",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2ccc(NC(=O)C(c3cc4c(cc3NC(=O)C=C(C)C)CCN(C(=O)OCC)CC4)C(C)C)cc2CC1",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2ccc3c(c2CC1)NC(=O)C3C(C)C",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2ccc(NC(=O)CC(C)(C)c3cc4c(cc3NC(=O)C=C(C)C)CCN(C(=O)OCC)CC4)cc2CC1",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2cc3c(cc2CC1)C(C)(C)CC(=O)N3",
            "CCOC(=O)N1CCc2ccc(NC(=O)C=C(C)C)cc2CC1>>CCOC(=O)N1CCc2ccc3c(c2CC1)NC(=O)CC3(C)C",
        ]

    def test_parse_fingerprint_settings_ecfp(self):
        # Test the fingerprint parsing for ECFP4
        result = self.similarity_ranking.parse_fingerprint_settings("ECFP4")
        # Test that the result is a callable (function)
        self.assertTrue(callable(result))
        # Further checks can include checking the type of fingerprint etc.

    # def test_calculate_fingerprint(self):
    #     # Test fingerprint calculation
    #     mol = Chem.MolFromSmiles('c1ccccc1')
    #     fp = self.similarity_ranking.calculate_fingerprint(mol, 'ECFP4')
    #     # Check that a fingerprint is returned
    #     self.assertIsInstance(fp, Chem.rdchem.ExplicitBitVect)

    # def test_calculate_tanimoto_similarity(self):
    #     # Create two fingerprints and calculate their similarity
    #     mol1 = Chem.MolFromSmiles('c1ccccc1')
    #     mol2 = Chem.MolFromSmiles('c1ccccn1')
    #     fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
    #     fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
    #     similarity = self.similarity_ranking.calculate_tanimoto_similarity(fp1, fp2)
    #     # Check that similarity is a float and between 0 and 1
    #     self.assertIsInstance(similarity, float)
    #     self.assertTrue(0 <= similarity <= 1)

    # def test_check_balance(self):
    #     # Test the check balance functionality
    #     reactants = 'CCO'
    #     products = 'OCC'
    #     result = self.similarity_ranking.check_balance(reactants, products)
    #     # Check that balance is reported correctly
    #     self.assertTrue(result)

    # def test_process_reaction_smiles(self):
    #     # Test processing of a reaction SMILES string
    #     reaction_smiles = 'CCO>>OCC'
    #     fptypes = ['ECFP4']
    #     similarity, balanced = self.similarity_ranking.process_reaction_smiles(reaction_smiles, fptypes)
    #     # Check types and correctness of return values
    #     self.assertIsInstance(similarity, float)
    #     self.assertIsInstance(balanced, bool)
    #     self.assertTrue(balanced)

    # def test_fit(self):
    #     # Test the fit function with a set of reactions
    #     reactions = ['CCO>>OCC', 'CC>>CC']
    #     fptypes = ['ECFP4']
    #     rsmi_list, similarity_list = self.similarity_ranking.fit(reactions, fptypes)
    #     # Check that outputs are lists and have correct content
    #     self.assertIsInstance(rsmi_list, list)
    #     self.assertIsInstance(similarity_list, list)
    #     self.assertEqual(len(rsmi_list), len(similarity_list))


if __name__ == "__main__":
    unittest.main()
