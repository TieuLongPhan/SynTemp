import unittest
from syntemp.SynAAM.graphormer_wrapper import (
    map_with_graphormer,
    map_with_graphormer_batch,
)


class TestGraphormer(unittest.TestCase):
    def setUp(self):

        self.reaction_smiles = "C=O.N#[C-]~[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1"
        self.mapped_reaction = "[K+:5]~[C-:4]#[N:3].[O:2]=[CH2:1].[cH:10]1[cH:11][cH:12][cH:13][cH:14][c:9]1[C:7]([Cl:8])=[O:6]>>[cH:11]1[cH:10][c:9]([cH:14][cH:13][cH:12]1)[C:7](=[O:6])[CH:1]([C:4]#[N:3])[OH:2]"
        self.reaction_list = [
            "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
            "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
            "CCCC1(CC(=O)OCC)OCCc2c1[nH]c1c(C)c(C(=O)OCc3ccccc3)cc(Br)c21.N#[C-]~[Cu+]>>CCCC1(CC(=O)OCC)OCCc2c1[nH]c1c(C)c(C(=O)OCc3ccccc3)cc(C#N)c21",
            "Cc1cc(=O)c2c(N)ccc(C)c2o1.N#[C-]~[Cu+].O=N[O-]>>Cc1cc(=O)c2c(C#N)ccc(C)c2o1",
        ]
        self.expected_list = [
            "[C-:4]#[N:3].[K+:5].[O:2]=[CH2:1].[cH:10]1[cH:11][cH:12][cH:13][cH:14][c:9]1[C:7]([Cl:8])=[O:6]>>[cH:11]1[cH:10][c:9]([cH:14][cH:13][cH:12]1)[C:7](=[O:6])[CH:1]([C:4]#[N:3])[OH:2]",
            "[F:10][C:9]([F:11])([F:12])[c:8]1[cH:7][c:5]([Br:6])[cH:4][cH:3][c:2]1[F:1].[K+:15]~[C-:14]#[N:13].[OH2:16].[OH2:17]>>[F:10][C:9]([F:11])([F:12])[c:8]1[cH:7][c:5]([Br:6])[cH:4][cH:3][c:2]1[C:14](=[O:17])[OH:16]",
            "[Cu+:37]~[C-:36]#[N:35].[cH:27]1[cH:26][c:25]([cH:30][cH:29][cH:28]1)[CH2:24][O:23][C:21]([c:20]1[c:18]([CH3:19])[c:17]2[nH:16][c:15]3[C:4]([O:11][CH2:12][CH2:13][c:14]3[c:34]2[c:32]([Br:33])[cH:31]1)([CH2:5][C:6](=[O:7])[O:8][CH2:9][CH3:10])[CH2:3][CH2:2][CH3:1])=[O:22]>>[cH:27]1[cH:26][c:25]([cH:30][cH:29][cH:28]1)[CH2:24][O:23][C:21]([c:20]1[cH:31][c:32]([C:36]#[N:35])[c:34]2[c:17]([nH:16][c:15]3[c:14]2[CH2:13][CH2:12][O:11][C:4]3([CH2:5][C:6](=[O:7])[O:8][CH2:9][CH3:10])[CH2:3][CH2:2][CH3:1])[c:18]1[CH3:19])=[O:22]",
            "[Cu+:17]~[C-:16]#[N:15].[N:19]([O-:20])=[O:18].[NH2:8][c:7]1[cH:9][cH:10][c:11]([c:13]2[o:14][c:2]([cH:3][c:4]([c:6]12)=[O:5])[CH3:1])[CH3:12]>>[C:16](#[N:15])[c:7]1[c:6]2[c:4]([cH:3][c:2]([CH3:1])[o:14][c:13]2[c:11]([CH3:12])[cH:10][cH:9]1)=[O:5]",
        ]

    def test_map_with_graphormer_single_success(self):
        result = map_with_graphormer(self.reaction_smiles)
        print(result)
        self.assertEqual(result, self.mapped_reaction)

    def test_map_with_graphormer_single_failure(self):
        fail_reaction = "fail>>C"
        result = map_with_graphormer(fail_reaction)  # will return original reaction
        self.assertEqual(result, fail_reaction)

    def test_map_with_graphormer_batch(self):
        results = map_with_graphormer_batch(self.reaction_list, batch_size=4)
        self.assertEqual(results, self.expected_list)


if __name__ == "__main__":
    unittest.main()
