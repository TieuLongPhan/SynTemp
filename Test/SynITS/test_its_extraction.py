import unittest
from syntemp.SynITS.its_extraction import ITSExtraction
from syntemp.SynITS.its_construction import ITSConstruction


class TestITSExtraction(unittest.TestCase):

    def setUp(self):
        # Simple SMILES strings for testing
        self.smiles1 = "CC(C)C"
        self.smiles2 = "CC"
        self.mapped_smiles_list = [
            {
                "balanced": True,
                "R-id": "USPTO_50K_26",
                "reactions": "Clc1cnc2nc1Nc1ccc(OCCC3CCNCC3)c(c1)CCc1cncc(c1)N2.O=C=NCc1ccco1>>O=C(NCc1ccco1)N1CCC(CCOc2ccc3cc2CCc2cncc(c2)Nc2ncc(Cl)c(n2)N3)CC1",
                "local_mapper": "[Cl:36][c:35]1[cH:34][n:33][c:32]2[n:38][c:37]1[NH:39][c:20]1[cH:19][cH:18][c:17]([O:16][CH2:15][CH2:14][CH:13]3[CH2:12][CH2:11][NH:10][CH2:41][CH2:40]3)[c:22]([cH:21]1)[CH2:23][CH2:24][c:25]1[cH:26][n:27][cH:28][c:29]([cH:30]1)[NH:31]2.[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)[N:10]1[CH2:11][CH2:12][CH:13]([CH2:14][CH2:15][O:16][c:17]2[cH:18][cH:19][c:20]3[cH:21][c:22]2[CH2:23][CH2:24][c:25]2[cH:26][n:27][cH:28][c:29]([cH:30]2)[NH:31][c:32]2[n:33][cH:34][c:35]([Cl:36])[c:37]([n:38]2)[NH:39]3)[CH2:40][CH2:41]1",
                "rxn_mapper": "[Cl:36][c:35]1[cH:34][n:33][c:32]2[n:38][c:37]1[NH:39][c:20]1[cH:19][cH:18][c:17]([O:16][CH2:15][CH2:14][CH:13]3[CH2:12][CH2:11][NH:10][CH2:41][CH2:40]3)[c:22]([cH:21]1)[CH2:23][CH2:24][c:25]1[cH:26][n:27][cH:28][c:29]([cH:30]1)[NH:31]2.[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)[N:10]1[CH2:11][CH2:12][CH:13]([CH2:14][CH2:15][O:16][c:17]2[cH:18][cH:19][c:20]3[cH:21][c:22]2[CH2:23][CH2:24][c:25]2[cH:26][n:27][cH:28][c:29]([cH:30]2)[NH:31][c:32]2[n:33][cH:34][c:35]([Cl:36])[c:37]([n:38]2)[NH:39]3)[CH2:40][CH2:41]1",
                "graphormer": "[C:34](=[O:33])=[N:35][CH2:36][c:37]1[cH:38][cH:39][cH:40][o:41]1.[CH2:17]1[CH2:18][NH:19][CH2:20][CH2:21][CH:16]1[CH2:15][CH2:14][O:13][c:12]1[cH:11][cH:10][c:9]2[cH:23][c:22]1[CH2:24][CH2:25][c:26]1[cH:31][c:30]([NH:32][c:5]3[n:6][c:7]([c:2]([cH:3][n:4]3)[Cl:1])[NH:8]2)[cH:29][n:28][cH:27]1>>[CH2:17]1[CH2:18][N:19]([CH2:20][CH2:21][CH:16]1[CH2:15][CH2:14][O:13][c:12]1[cH:11][cH:10][c:9]2[cH:23][c:22]1[CH2:24][CH2:25][c:26]1[cH:31][c:30]([NH:32][c:5]3[n:6][c:7]([c:2]([cH:3][n:4]3)[Cl:1])[NH:8]2)[cH:29][n:28][cH:27]1)[C:34]([NH:35][CH2:36][c:37]1[cH:38][cH:39][cH:40][o:41]1)=[O:33]",
            },
            {
                "balanced": True,
                "R-id": "inequi",
                "reactions": "CC(C)(C)[Si](C)(C)OC(COCC=1C=CC=CC=1)CCC(CC=C)OS(=O)(=O)C.O>>C1=CC(=CC=C1)COCC2OC(CC2)CC=C.CC(C)(C)[Si](C)(C)O.CS(=O)(=O)O",
                "local_mapper": "[CH3:18][C:19]([CH3:20])([CH3:21])[Si:22]([CH3:23])([CH3:24])[O:11][CH:10]([CH2:9][O:8][CH2:7][c:3]1[cH:2][cH:1][cH:6][cH:5][cH:4]1)[CH2:14][CH2:13][CH:12]([CH2:15][CH:16]=[CH2:17])[O:30][S:27](=[O:28])(=[O:29])[CH3:26].[OH2:25]>>[cH:1]1[cH:2][c:3]([CH2:7][O:8][CH2:9][CH:10]2[O:11][CH:12]([CH2:15][CH:16]=[CH2:17])[CH2:13][CH2:14]2)[cH:4][cH:5][cH:6]1.[CH3:18][C:19]([CH3:20])([CH3:21])[Si:22]([CH3:23])([CH3:24])[OH:25].[CH3:26][S:27](=[O:28])(=[O:29])[OH:30]",
                "rxn_mapper": "[CH3:18][C:19]([CH3:20])([CH3:21])[Si:22]([CH3:23])([CH3:24])[O:25][CH:10]([CH2:9][O:8][CH2:7][C:5]1=[CH:6][CH:1]=[CH:2][CH:3]=[CH:4]1)[CH2:14][CH2:13][CH:12]([CH2:15][CH:16]=[CH2:17])[O:30][S:27](=[O:28])(=[O:29])[CH3:26].[OH2:11]>>[CH:1]1=[CH:2][C:3]([CH2:7][O:8][CH2:9][CH:10]2[O:11][CH:12]([CH2:15][CH:16]=[CH2:17])[CH2:13][CH2:14]2)=[CH:4][CH:5]=[CH:6]1.[CH3:18][C:19]([CH3:20])([CH3:21])[Si:22]([CH3:23])([CH3:24])[OH:25].[CH3:26][S:27](=[O:28])(=[O:29])[OH:30]",
                "graphormer": "[CH3:1][C:2]([CH3:3])([CH3:4])[Si:5]([CH3:6])([CH3:7])[O:8][CH:9]([CH2:10][O:11][CH2:12][C:13]=1[CH:14]=[CH:15][CH:16]=[CH:17][CH:18]=1)[CH2:19][CH2:20][CH:21]([CH2:22][CH:23]=[CH2:24])[O:25][S:26](=[O:27])(=[O:28])[CH3:29].[OH2:30]>>[CH3:1][C:2]([CH3:3])([CH3:4])[Si:5]([CH3:6])([CH3:7])[OH:8].[CH:17]1=[CH:18][C:13](=[CH:14][CH:15]=[CH:16]1)[CH2:12][O:11][CH2:10][CH:9]1[O:30][CH:21]([CH2:20][CH2:19]1)[CH2:22][CH:23]=[CH2:24].[O:27]=[S:26](=[O:28])([OH:25])[CH3:29]",
            },
        ]
        self.mapper_names = ["local_mapper", "rxn_mapper", "graphormer"]

    def test_graph_from_smiles(self):
        graph = ITSExtraction.graph_from_smiles(self.smiles1)
        self.assertEqual(len(graph.nodes()), 4)
        self.assertEqual(len(graph.edges()), 3)

    def test_check_equivariant_graph(self):
        react_local_mapper, prod_local_mapper = self.mapped_smiles_list[0][
            "local_mapper"
        ].split(">>")
        G_local = ITSExtraction.graph_from_smiles(react_local_mapper)
        H_local = ITSExtraction.graph_from_smiles(prod_local_mapper)
        ITS_local = ITSConstruction.ITSGraph(G_local, H_local)

        react_rxn_mapper, prod_rxn_mapper = self.mapped_smiles_list[0][
            "rxn_mapper"
        ].split(">>")
        G_rxn = ITSExtraction.graph_from_smiles(react_rxn_mapper)
        H_rxn = ITSExtraction.graph_from_smiles(prod_rxn_mapper)
        ITS_rxn = ITSConstruction.ITSGraph(G_rxn, H_rxn)

        react_graphormer, prod_graphormer = self.mapped_smiles_list[0][
            "graphormer"
        ].split(">>")
        G_graphormer = ITSExtraction.graph_from_smiles(react_graphormer)
        H_graphormer = ITSExtraction.graph_from_smiles(prod_graphormer)
        ITS_graphormer = ITSConstruction.ITSGraph(G_graphormer, H_graphormer)

        classified, equivariant = ITSExtraction.check_equivariant_graph(
            [ITS_local, ITS_rxn, ITS_graphormer]
        )
        self.assertEqual(equivariant, 2)  # Expect 2 isomorphic graphs
        self.assertEqual(classified, [(0, 1), (0, 2)])  # matching pairs

    def test_process_mapped_smiles(self):
        graphs_by_map_correct, _ = ITSExtraction.process_mapped_smiles(
            self.mapped_smiles_list[0],
            self.mapper_names,
        )
        self.assertIsNotNone(graphs_by_map_correct["ITSGraph"])
        self.assertIsNotNone(graphs_by_map_correct["GraphRules"])

    def test_parallel_process_smiles(self):
        # Measure execution time with a single job

        results, results_wrong = ITSExtraction.parallel_process_smiles(
            self.mapped_smiles_list, self.mapper_names, n_jobs=2, verbose=0
        )
        print(results_wrong)
        # Equivalent AAM
        self.assertIsNotNone(results[0]["ITSGraph"])
        self.assertIsNotNone(results[0]["GraphRules"])

        # Inequivalent AAM
        self.assertEqual(results_wrong[0]["equivariant"], 0)


if __name__ == "__main__":
    unittest.main()
