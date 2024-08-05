import unittest
from pathlib import Path
import networkx as nx
from syntemp.SynUtils.graph_utils import (
    is_acyclic_graph,
    is_single_cyclic_graph,
    is_complex_cyclic_graph,
    check_graph_type,
    get_cycle_member_rings,
)

root_dir = Path(__file__).parents[2]


class TestUncertainRefinement(unittest.TestCase):

    def test_is_acyclic_graph(self):
        # Create an acyclic graph (tree)
        tree_graph = nx.balanced_tree(2, 3)

        # Create a cyclic graph
        cyclic_graph = nx.cycle_graph(4)

        # Test is_acyclic_graph method
        self.assertTrue(is_acyclic_graph(tree_graph))
        self.assertFalse(is_acyclic_graph(cyclic_graph))

    def test_is_single_cyclic_graph(self):
        # Create a single cyclic graph (cycle)
        single_cycle_graph = nx.cycle_graph(5)

        # Create a graph with two cycles
        two_cycle_graph = nx.Graph()
        two_cycle_graph.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)])

        # Test is_single_cyclic_graph method
        self.assertTrue(is_single_cyclic_graph(single_cycle_graph))
        self.assertFalse(is_single_cyclic_graph(two_cycle_graph))

    def test_is_complex_cyclic_graph(self):
        # Create a graph with two cycles
        two_cycle_graph = nx.Graph()
        two_cycle_graph.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)])

        # Create a graph with three cycles
        three_cycle_graph = nx.Graph()
        three_cycle_graph.add_edges_from(
            [(1, 2), (2, 3), (3, 1), (1, 4), (4, 5), (5, 1)]
        )

        # Test is_complex_cyclic_graph method
        self.assertTrue(is_complex_cyclic_graph(two_cycle_graph))
        self.assertTrue(is_complex_cyclic_graph(three_cycle_graph))

    def test_check_graph_type(self):
        # Create an acyclic graph (tree)
        tree_graph = nx.balanced_tree(2, 3)

        # Create a single cyclic graph (cycle)
        single_cycle_graph = nx.cycle_graph(5)

        # Create a graph with two cycles
        two_cycle_graph = nx.Graph()
        two_cycle_graph.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)])

        # Test check_graph_type method
        self.assertEqual(check_graph_type(tree_graph), "Acyclic")
        self.assertEqual(check_graph_type(single_cycle_graph), "Single Cyclic")
        self.assertEqual(check_graph_type(two_cycle_graph), "Combinatorial Cyclic")

    def test_get_cycle_member_rings(self):
        # Create a test graph with cycles of different sizes
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 1)])  # Cycle of size 4
        G.add_edges_from([(9, 10), (10, 11), (11, 9)])  # Cycle of size 3

        # Call the function to get the cycle member rings
        member_rings = get_cycle_member_rings(G)

        # Define the expected result
        expected_result = [3, 4]

        # Assert that the result matches the expected result
        self.assertEqual(member_rings, expected_result)

        G2 = nx.Graph()
        G2.add_edges_from([(5, 6), (6, 7), (7, 8), (8, 5), (5, 7)])
        member_rings = get_cycle_member_rings(G2)
        self.assertEqual(member_rings, [3, 3])


if __name__ == "__main__":
    unittest.main()
