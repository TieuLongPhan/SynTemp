import unittest
from syntemp.SynAAM.rxn_mapper_wrapper import (
    map_with_rxn_mapper,
    map_with_rxn_mapper_batch,
)


class TestRXNMapper(unittest.TestCase):
    def setUp(self):

        self.reaction_smiles = "C=O.N#[C-]~[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1"
        self.mapped_reaction = "[CH2:3]=[O:4].[N:1]#[C-:2]~[K+].[O:6]=[C:5](Cl)[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[N:1]#[C:2][CH:3]([OH:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1"
        self.reaction_list = [
            "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
            "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
            "CCCC1(CC(=O)OCC)OCCc2c1[nH]c1c(C)c(C(=O)OCc3ccccc3)cc(Br)c21.N#[C-]~[Cu+]>>CCCC1(CC(=O)OCC)OCCc2c1[nH]c1c(C)c(C(=O)OCc3ccccc3)cc(C#N)c21",
            "Cc1cc(=O)c2c(N)ccc(C)c2o1.N#[C-]~[Cu+].O=N[O-]>>Cc1cc(=O)c2c(C#N)ccc(C)c2o1",
        ]
        self.expected_list = [
            "[CH2:3]=[O:4].[N:1]#[C-:2].[K+].[O:6]=[C:5](Cl)[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1>>[N:1]#[C:2][CH:3]([OH:4])[C:5](=[O:6])[c:7]1[cH:8][cH:9][cH:10][cH:11][cH:12]1",
            "F[c:4]1[cH:5][cH:6][c:7]([Br:8])[cH:9][c:10]1[C:11]([F:12])([F:13])[F:14].N#[C-:2]~[K+].[OH2:1].[OH2:3]>>[O:1]=[C:2]([OH:3])[c:4]1[cH:5][cH:6][c:7]([Br:8])[cH:9][c:10]1[C:11]([F:12])([F:13])[F:14]",
            "[CH3:1][CH2:2][CH2:3][C:4]1([CH2:5][C:6](=[O:7])[O:8][CH2:9][CH3:10])[O:11][CH2:12][CH2:13][c:14]2[c:15]1[nH:16][c:17]1[c:18]([CH3:19])[c:20]([C:21](=[O:22])[O:23][CH2:24][c:25]3[cH:26][cH:27][cH:28][cH:29][cH:30]3)[cH:31][c:32](Br)[c:35]21.[N:34]#[C-:33]~[Cu+]>>[CH3:1][CH2:2][CH2:3][C:4]1([CH2:5][C:6](=[O:7])[O:8][CH2:9][CH3:10])[O:11][CH2:12][CH2:13][c:14]2[c:15]1[nH:16][c:17]1[c:18]([CH3:19])[c:20]([C:21](=[O:22])[O:23][CH2:24][c:25]3[cH:26][cH:27][cH:28][cH:29][cH:30]3)[cH:31][c:32]([C:33]#[N:34])[c:35]21",
            "[CH3:1][c:2]1[cH:3][c:4](=[O:5])[c:6]2[c:7](N)[cH:10][cH:11][c:12]([CH3:13])[c:14]2[o:15]1.[N:9]#[C-:8]~[Cu+].O=N[O-]>>[CH3:1][c:2]1[cH:3][c:4](=[O:5])[c:6]2[c:7]([C:8]#[N:9])[cH:10][cH:11][c:12]([CH3:13])[c:14]2[o:15]1",
        ]

    def test_map_with_rxn_mapper_single_success(self):
        result = map_with_rxn_mapper(self.reaction_smiles)
        print(result)
        self.assertEqual(result, self.mapped_reaction)

    def test_map_with_rxn_mapper_single_failure(self):
        fail_reaction = "fail>>C"
        result = map_with_rxn_mapper(fail_reaction)  # will return original reaction
        self.assertEqual(result, fail_reaction)

    def test_map_with_rxn_mapper_batch(self):
        results = map_with_rxn_mapper_batch(self.reaction_list, batch_size=4)
        self.assertEqual(results, self.expected_list)


if __name__ == "__main__":
    unittest.main()
