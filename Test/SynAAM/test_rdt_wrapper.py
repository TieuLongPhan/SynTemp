import unittest
import importlib.resources
from syntemp.SynAAM.rdt_wrapper import map_with_rdt, map_with_rdt_batch


class TestRDT(unittest.TestCase):

    def setUp(self):
        self.working_dir = "./"
        try:
            self.rdt_jar_path = importlib.resources.files("syntemp.SynAAM").joinpath(
                "RDT_2.4.1.jar"
            )
        except FileNotFoundError:
            raise FileNotFoundError(
                "The RDT jar file is not found in the specified package resources."
            )

        self.reaction_smiles = "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1"
        self.mapped_reaction = "[O:1]=[C:2]([Cl:13])[c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1.[O:9]=[CH2:10].[C-:11]#[N:12].[K+:14]>>[cH:8]1[cH:7][cH:6][cH:5][cH:4][c:3]1[C:2]([CH:10]([C:11]#[N:12])[OH:9])=[O:1]"
        self.reaction_list = [
            "C=O.N#[C-].[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1",
            "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
            "CCCC1(CC(=O)OCC)OCCc2c1[nH]c1c(C)c(C(=O)OCc3ccccc3)cc(Br)c21.N#[C-].[Cu+]>>CCCC1(CC(=O)OCC)OCCc2c1[nH]c1c(C)c(C(=O)OCc3ccccc3)cc(C#N)c21",
            "Cc1cc(=O)c2c(N)ccc(C)c2o1.N#[C-]~[Cu+].O=N[O-]>>Cc1cc(=O)c2c(C#N)ccc(C)c2o1",
        ]
        self.expected_list = [
            "[O:1]=[C:2]([Cl:13])[c:3]1[cH:4][cH:5][cH:6][cH:7][cH:8]1.[O:9]=[CH2:10].[C-:11]#[N:12].[K+:14]>>[cH:8]1[cH:7][cH:6][cH:5][cH:4][c:3]1[C:2]([CH:10]([C:11]#[N:12])[OH:9])=[O:1]",
            "Fc1ccc(Br)cc1C(F)(F)F.N#[C-]~[K+].O.O>>O=C(O)c1ccc(Br)cc1C(F)(F)F",
            "[O:1]=[C:2]([O:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1)[c:11]2[cH:12][c:13]([Br:36])[c:14]3[c:15]([nH:16][c:17]4[c:18]3[CH2:19][CH2:20][O:21][C:22]4([CH2:23][C:24](=[O:25])[O:26][CH2:27][CH3:28])[CH2:29][CH2:30][CH3:31])[c:32]2[CH3:33].[C-:34]#[N:35].[Cu+:37]>>[CH3:33][c:32]1[c:11]([cH:12][c:13]([C:34]#[N:35])[c:14]2[c:15]1[nH:16][c:17]3[c:18]2[CH2:19][CH2:20][O:21][C:22]3([CH2:29][CH2:30][CH3:31])[CH2:23][C:24](=[O:25])[O:26][CH2:27][CH3:28])[C:2](=[O:1])[O:3][CH2:4][c:5]4[cH:6][cH:7][cH:8][cH:9][cH:10]4",
            "Cc1cc(=O)c2c(N)ccc(C)c2o1.N#[C-]~[Cu+].O=N[O-]>>Cc1cc(=O)c2c(C#N)ccc(C)c2o1",
        ]

    def test_map_with_rdt_single_success(self):
        result = map_with_rdt(self.reaction_smiles, self.rdt_jar_path, self.working_dir)
        self.assertEqual(result, self.mapped_reaction)

    def test_map_with_rdt_single_failure(self):
        fail_reaction = "fail>>C"
        result = map_with_rdt(
            fail_reaction, self.rdt_jar_path, self.working_dir
        )  # will return original reaction
        self.assertEqual(result, fail_reaction)

    def test_map_with_rdt_single_failure_smiles(self):
        fail_reaction = (
            "C=O.N#[C-]~[K+].O=C(Cl)c1ccccc1>>N#CC(O)C(=O)c1ccccc1"  # ~ in smiles
        )
        result = map_with_rdt(
            fail_reaction, self.rdt_jar_path, self.working_dir
        )  # will return original reaction
        self.assertEqual(result, fail_reaction)

    def test_map_with_rdt_batch(self):
        results = map_with_rdt_batch(
            self.reaction_list, self.rdt_jar_path, self.working_dir, 4
        )
        self.assertEqual(results, self.expected_list)


if __name__ == "__main__":
    unittest.main()
