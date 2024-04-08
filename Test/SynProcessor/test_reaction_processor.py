import unittest
import sys
from pathlib import Path
root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynProcessor.reaction_processor import ReactionProcessor

class TestReactionProcessor(unittest.TestCase):
    def test_label_reactions(self):
        reaction_dict = {
            'R-id': 'R1',
            'new_reaction': 'CCO.[O]>>CCO[O]'
        }
        expected_label = 'Oxidation'
        processed_reaction = ReactionProcessor.label_reactions(reaction_dict)
        self.assertEqual(processed_reaction['label'], expected_label)
        self.assertTrue('reactants' in processed_reaction and 'products' in processed_reaction)

    def test_calculate_charge(self):
        smiles = "[NH4+]"
        charge = ReactionProcessor.calculate_charge(smiles)
        self.assertEqual(charge, 1)

    def test_process_reaction(self):
        reaction = {
            'R-id': 'R1',
            'new_reaction': 'CCO>>CCO[O]',
            'products': 'CCO[O]'
        }
        processed_reaction = ReactionProcessor.process_reaction(reaction)
        self.assertEqual(processed_reaction['total_charge_in_products'], 0)  # Assuming CCO[O] is neutral for this example

    def test_process_reactions_parallel(self): # add new test case
        reaction_list = [
            {'R-id': 'R1', 'new_reaction': 'CCO.[O]>>CCO[O]'},
            {'R-id': 'R2', 'new_reaction': 'C=C.[H]>>CC'}
        ]
        processed_reactions = ReactionProcessor.process_reactions_parallel(reaction_list, n_jobs=2)
        self.assertEqual(len(processed_reactions), len(reaction_list))
        self.assertTrue(all('label' in reaction for reaction in processed_reactions))

    def test_sum_of_charge_in_products(self):
        reaction_list = [
            {'R-id': 'R1', 'new_reaction': 'CCO>>CCO[O]', 'products': 'CCO[O]'},
            {'R-id': 'R2', 'new_reaction': 'C=C.[H]>>CC', 'products': 'CC'}
        ]
        processed_reactions = ReactionProcessor.sum_of_charge_in_products(reaction_list, n_jobs=2)
        self.assertTrue(all('total_charge_in_products' in reaction for reaction in processed_reactions))
        # Assuming both products are neutral for this example
        self.assertTrue(all(reaction['total_charge_in_products'] == 0 for reaction in processed_reactions))

if __name__ == "__main__":
    unittest.main()
