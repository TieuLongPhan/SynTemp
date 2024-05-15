import unittest
from pathlib import Path


from SynTemp.SynAAM.atom_mappers import (
    map_with_rxn_mapper,
    map_with_graphormer,
    map_with_local_mapper,
    map_with_rdt,
    remove_atom_mapping,
    normalize_smiles,
)
from rxnmapper import RXNMapper
from localmapper import localmapper
from rdkit import Chem

root_dir = Path(__file__).parents[2]


class TestMappingFunctions(unittest.TestCase):

    def setUp(self):
        # This SMILES string will be used in each test method
        self.smiles = "O=C=NCc1ccco1>>O=C(NCc1ccco1)N1CCC(CCOc2ccc3cc2CCc2cncc(c2)Nc2ncc(Cl)c(n2)N3)CC1"
        self.rxn_mapper = RXNMapper()

    def test_map_with_rxn_mapper(self):
        #
        result = map_with_rxn_mapper(self.smiles, self.rxn_mapper)
        mapped_result = "[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)[N:10]1[CH2:11][CH2:12][CH:13]([CH2:14][CH2:15][O:16][c:17]2[cH:18][cH:19][c:20]3[cH:21][c:22]2[CH2:23][CH2:24][c:25]2[cH:26][n:27][cH:28][c:29]([cH:30]2)[NH:31][c:32]2[n:33][cH:34][c:35]([Cl:36])[c:37]([n:38]2)[NH:39]3)[CH2:40][CH2:41]1"
        self.assertEqual(result, mapped_result)

    def test_map_with_graphormer(self):

        result = map_with_graphormer(self.smiles)

        mapped_result = "[C:2](=[O:1])=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[CH2:12]1[CH2:11][N:10]([CH2:41][CH2:40][CH:13]1[CH2:14][CH2:15][O:16][c:17]1[cH:18][cH:19][c:20]2[cH:21][c:22]1[CH2:23][CH2:24][c:25]1[cH:30][c:29]([NH:31][c:32]3[n:38][c:37]([c:35]([cH:34][n:33]3)[Cl:36])[NH:39]2)[cH:28][n:27][cH:26]1)[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)=[O:1]"
        self.assertEqual(result, mapped_result)

    def test_map_with_local_mapper(self):

        result = map_with_local_mapper(self.smiles, mapper=localmapper())
        print(result)
        mapped_result = "[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)N1CCC(CCOc2ccc3cc2CCc2cncc(c2)Nc2ncc(Cl)c(n2)N3)CC1"
        self.assertEqual(result, mapped_result)

    def test_map_with_rdt(self):

        # Ensure the RDT jar path and working directory are correctly set
        rdt_jar_path = f"{root_dir}/Data/RDT_2.4.1.jar"
        result = map_with_rdt(self.smiles, rdt_jar_path, working_dir=root_dir)
        mapped_result = "[O:1]=[C:2]=[N:3][CH2:4][c:5]1[o:6][cH:7][cH:8][cH:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[o:6][cH:7][cH:8][cH:9]1)[N:10]2[CH2:11][CH2:12][CH:13]([CH2:14][CH2:15][O:16][c:17]3[cH:18][cH:19][c:20]4[cH:21][c:22]3[CH2:23][CH2:24][c:25]5[cH:26][n:27][cH:28][c:29]([cH:30]5)[NH:31][c:32]6[n:33][cH:34][c:35]([Cl:36])[c:37]([n:38]6)[NH:39]4)[CH2:40][CH2:41]2"
        self.assertEqual(result, mapped_result)

    def test_remove_atom_mapping(self):

        unmapped_smiles = remove_atom_mapping(
            "[O:1]=[C:2]=[N:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1>>[O:1]=[C:2]([NH:3][CH2:4][c:5]1[cH:6][cH:7][cH:8][o:9]1)[N:10]1[CH2:11][CH2:12][CH:13]([CH2:14][CH2:15][O:16][c:17]2[cH:18][cH:19][c:20]3[cH:21][c:22]2[CH2:23][CH2:24][c:25]2[cH:26][n:27][cH:28][c:29]([cH:30]2)[NH:31][c:32]2[n:33][cH:34][c:35]([Cl:36])[c:37]([n:38]2)[NH:39]3)[CH2:40][CH2:41]1"
        )
        self.assertEqual(
            normalize_smiles(unmapped_smiles), normalize_smiles(self.smiles)
        )

    def test_normalize_smiles(self):
        # Test the normalization of the specific reaction SMILES string
        normalized_smiles = normalize_smiles(self.smiles)
        self.assertIsInstance(
            Chem.rdChemReactions.ReactionFromSmarts(normalized_smiles),
            Chem.rdChemReactions.ChemicalReaction,
        )


if __name__ == "__main__":
    unittest.main()
