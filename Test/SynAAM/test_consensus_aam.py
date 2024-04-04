import unittest
import sys
from pathlib import Path
root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynAAM.consensus_aam import ConsensusAAM
from rxnmapper import RXNMapper
import os

#TODO, add real test cases
class TestConsensusAAM(unittest.TestCase):

    def setUp(self):
        self.reactions_list = [
            {'reactions': 'C=CCC(N)(c1ccccc1)c1ccccc1.O=CC(C(=O)N1CCOCC1)c1ccccc1>>C=CCC(N=C(c1ccccc1)c1ccccc1)C(C(=O)N1CCOCC1)c1ccccc1'},
            {'reactions': 'C=CCSCC(C)=O.C=CC(O)C(C)N>>CC(=O)CSCC=CC(O)C(C)N'}  
        ]
        self.consensus_aam = ConsensusAAM(self.reactions_list, rsmi_column='reactions', save_dir='/tmp', 
                                          mapper_types=['rxn_mapper', 'graphormer', 'local_mapper', 'rdt'])
        self.rxnmapper = RXNMapper()
        self.rdt_jar_path = f'{root_dir}/Data/RDT_2.4.1.jar'
       

    def test_extract_smiles(self):
        smiles_list = self.consensus_aam.extract_smiles()
        self.assertEqual(smiles_list, ['C=CCC(N)(c1ccccc1)c1ccccc1.O=CC(C(=O)N1CCOCC1)c1ccccc1>>C=CCC(N=C(c1ccccc1)c1ccccc1)C(C(=O)N1CCOCC1)c1ccccc1', 
                                       'C=CCSCC(C)=O.C=CC(O)C(C)N>>CC(=O)CSCC=CC(O)C(C)N'])

    def test_process_batch(self):
        # Note: This test will run the actual mapping, affecting the reactions_list.
        batch = self.consensus_aam.extract_smiles()
        self.consensus_aam.process_batch(batch, self.rxnmapper, self.rdt_jar_path, './')

        # Check if the mapping is applied. The exact result depends on the mapper's implementation.
        for reaction_dict in self.consensus_aam.reactions_list:
            self.assertIn('rxn_mapper', reaction_dict)
            self.assertIn('graphormer', reaction_dict)
            self.assertIn('local_mapper', reaction_dict)
            self.assertIn('rdt', reaction_dict)

        # Ensure the file is saved. This doesn't check the file's content, only its existence.
        self.assertTrue(os.path.exists(self.consensus_aam.filename))

    def test_run(self):
        # Running this will process the entire list. It's more of an integration test.
        updated_list = self.consensus_aam.fit(batch_size=1, rxn_mapper=self.rxnmapper, rdt_jar_path=self.rdt_jar_path, working_dir='./')
        self.assertEqual(len(updated_list), len(self.reactions_list))

        # Each reaction should now include the 'rxn_mapper' key.
        for reaction_dict in updated_list:
            self.assertIn('rxn_mapper', reaction_dict)
            self.assertIn('graphormer', reaction_dict)
            self.assertIn('local_mapper', reaction_dict)
            self.assertIn('rdt', reaction_dict)

        # Clean up the generated file
        if os.path.exists(self.consensus_aam.filename):
            os.remove(self.consensus_aam.filename)

if __name__ == '__main__':
    unittest.main()
