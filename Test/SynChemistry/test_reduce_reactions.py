import unittest
from SynTemp.SynChemistry.reduce_reactions import ReduceReactions


class TestReduceReactions(unittest.TestCase):

    def test_check_balance(self):
        # Test cases where reactions are balanced and not balanced
        self.assertTrue(ReduceReactions.check_balance("CC", "CC"))  # Balanced
        self.assertFalse(ReduceReactions.check_balance("CC", "CCC"))  # Not balanced

    def test_process_rsmi(self):
        # Check correct processing and balancing check
        self.assertEqual(ReduceReactions.process_rsmi("CC>>CC"), True)  # Balanced
        self.assertEqual(ReduceReactions.process_rsmi("CC>>CCC"), False)  # Not balanced

    def test_process_list_of_rsmi(self):
        # Check processing of multiple reactions
        reactions = ["CC>>CC", "CC>>CCC", "C>>C"]
        expected_results = ["CC>>CC", "C>>C"]  # Filtered and balanced
        self.assertEqual(
            ReduceReactions.process_list_of_rsmi(reactions), expected_results
        )

    def test_process_list_of_dicts(self):
        # Check processing on list of dictionaries
        database = [
            {"reactions": ["CC>>CC", "CC>>CCC", "C>>C"]},
            {"reactions": ["CCC>>CCC"]},
        ]
        expected = [{"reactions": ["CC>>CC", "C>>C"]}, {"reactions": ["CCC>>CCC"]}]
        self.assertEqual(
            ReduceReactions.process_list_of_dicts(database, "reactions"), expected
        )


if __name__ == "__main__":
    unittest.main()
