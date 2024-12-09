import unittest
import networkx as nx
from synutility.SynIO.data_type import load_from_pickle
from syntemp.SynITS.hydrogen_utils import (
    check_explicit_hydrogen,
    check_hcount_change,
    get_cycle_member_rings,
    get_priority,
)


class TestGraphFunctions(unittest.TestCase):

    def setUp(self):
        # Create a test graph for the tests
        self.data = load_from_pickle("./Data/Testcase/hydrogen_test.pkl.gz")

    def test_check_explicit_hydrogen(self):
        # Test the check_explicit_hydrogen function
        # Note, usually only appear in reactants (+H2 reactions)
        count_r, hydrogen_nodes_r = check_explicit_hydrogen(
            self.data[20]["ITSGraph"][0]
        )
        self.assertEqual(count_r, 2)
        self.assertEqual(hydrogen_nodes_r, [45, 46])

    def test_check_hcount_change(self):
        # Test the check_hcount_change function
        max_change = check_hcount_change(
            self.data[20]["ITSGraph"][0], self.data[20]["ITSGraph"][0]
        )
        self.assertEqual(max_change, 2)

    def test_get_cycle_member_rings_minimal(self):
        # Test get_cycle_member_rings with 'minimal' cycles
        member_rings = get_cycle_member_rings(self.data[1]["GraphRules"][2], "minimal")
        self.assertEqual(member_rings, [4])  # Cycles of size 4 and 3

    def test_get_priority(self):
        # Create a test graph for the tests
        self.graph = nx.Graph()
        self.graph.add_nodes_from(
            [
                (1, {"element": "H", "hcount": 2}),
                (2, {"element": "C", "hcount": 1}),
                (3, {"element": "H", "hcount": 1}),
            ]
        )
        self.graph.add_edges_from([(1, 2), (2, 3)])

        # Create another graph for `check_hcount_change` tests
        self.prod_graph = nx.Graph()
        self.prod_graph.add_nodes_from(
            [
                (1, {"element": "H", "hcount": 1}),
                (2, {"element": "C", "hcount": 1}),
                (3, {"element": "H", "hcount": 2}),
            ]
        )
        self.prod_graph.add_edges_from([(1, 2), (2, 3)])

        # Create a more complex graph for cycle tests
        self.complex_graph = nx.Graph()
        self.complex_graph.add_edges_from(
            [
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 1),  # A simple square cycle
                (3, 5),
                (5, 6),
                (6, 3),  # Another cycle
            ]
        )
        reaction_centers = [self.graph, self.prod_graph, self.complex_graph]

        # Get priority indices
        priority_indices = get_priority(reaction_centers)

        self.assertEqual(priority_indices, [0, 1])


if __name__ == "__main__":
    unittest.main()
