import unittest
from syntemp.SynAAM.aam_postprocess import AMMPostprocessor


class TestAMMPostprocessor(unittest.TestCase):

    def setUp(self):
        self.smiles_with_mapping = "[Cl:36][c:35]1[cH:34][n:33][c:32]2[n:38][c:37]1[NH:39][c:20]1[cH:19][cH:18][c:17]([O:16][CH2:15][CH2:14][CH:13]3[CH2:12][CH2:11][NH:10][CH2:41][CH2:40]3)[c:22]([cH:21]1)[CH2:23][CH2:24][c:25]1[cH:26][n:27][cH:28][c:29]([cH:30]1)[NH:31]2.[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1"
        self.reaction_smiles_with_mapping = "[Cl:36][c:35]1[cH:34][n:33][c:32]2[n:38][c:37]1[NH:39][c:20]1[cH:19][cH:18][c:17]([O:16][CH2:15][CH2:14][CH:13]3[CH2:12][CH2:11][NH:10][CH2:41][CH2:40]3)[c:22]([cH:21]1)[CH2:23][CH2:24][c:25]1[cH:26][n:27][cH:28][c:29]([cH:30]1)[NH:31]2.[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)[N:10]1[CH2:11][CH2:12][CH:13]([CH2:14][CH2:15][O:16][c:17]2[cH:18][cH:19][c:20]3[cH:21][c:22]2[CH2:23][CH2:24][c:25]2[cH:26][n:27][cH:28][c:29]([cH:30]2)[NH:31][c:32]2[n:33][cH:34][c:35]([Cl:36])[c:37]([n:38]2)[NH:39]3)[CH2:40][CH2:41]1"

    def test_extract_mappings_and_count_atoms(self):
        _, count = AMMPostprocessor.extract_mappings_and_count_atoms(
            self.smiles_with_mapping
        )
        self.assertEqual(count, 41)

    def test_is_consistent_mapping(self):
        # Testing with a consistent reaction mapping
        self.assertTrue(
            AMMPostprocessor.is_consistent_mapping(self.reaction_smiles_with_mapping)
        )

        # Testing with an inconsistent reaction mapping
        inconsistent_reaction_smiles = "[CH3:1][CH2:2][OH:3]>>[CH2:1]=[CH2:2]"
        self.assertFalse(
            AMMPostprocessor.is_consistent_mapping(inconsistent_reaction_smiles)
        )

    def test_postprocess(self):
        mapped_smiles = {
            "mapper1": self.reaction_smiles_with_mapping,
            "mapper2": self.reaction_smiles_with_mapping,
        }
        mapper_names = ["mapper1", "mapper2"]
        result = AMMPostprocessor.postprocess(mapped_smiles, mapper_names, 2)
        self.assertTrue(result["Valid"])

    def test_parallel_postprocess(self):
        mapped_smiles_list = [
            {
                "mapper1": self.reaction_smiles_with_mapping,
                "mapper2": self.reaction_smiles_with_mapping,
            },
            {
                "mapper1": self.reaction_smiles_with_mapping,
                "mapper2": self.reaction_smiles_with_mapping,
            },
            {
                "mapper1": self.reaction_smiles_with_mapping,
                "mapper2": self.reaction_smiles_with_mapping,
            },
            {
                "mapper1": self.reaction_smiles_with_mapping,
                "mapper2": self.reaction_smiles_with_mapping,
            },
            {
                "mapper1": self.reaction_smiles_with_mapping,
                "mapper2": self.reaction_smiles_with_mapping,
            },
        ]
        mapper_names = ["mapper1", "mapper2"]
        results = AMMPostprocessor.parallel_postprocess(
            mapped_smiles_list, mapper_names, 2, n_jobs=4, verbose=1
        )
        self.assertTrue(results[0]["Valid"])


if __name__ == "__main__":
    unittest.main()
