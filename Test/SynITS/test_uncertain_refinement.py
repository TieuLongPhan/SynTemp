import unittest
import sys
from pathlib import Path

root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
import time
import networkx as nx
from SynTemp.SynITS.uncertain_refinement import UncertainRefinement


class TestUncertainRefinement(unittest.TestCase):

    def test_is_acyclic_graph(self):
        # Create an acyclic graph (tree)
        tree_graph = nx.balanced_tree(2, 3)

        # Create a cyclic graph
        cyclic_graph = nx.cycle_graph(4)

        # Test is_acyclic_graph method
        self.assertTrue(UncertainRefinement.is_acyclic_graph(tree_graph))
        self.assertFalse(UncertainRefinement.is_acyclic_graph(cyclic_graph))

    def test_is_single_cyclic_graph(self):
        # Create a single cyclic graph (cycle)
        single_cycle_graph = nx.cycle_graph(5)

        # Create a graph with two cycles
        two_cycle_graph = nx.Graph()
        two_cycle_graph.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)])

        # Test is_single_cyclic_graph method
        self.assertTrue(UncertainRefinement.is_single_cyclic_graph(single_cycle_graph))
        self.assertFalse(UncertainRefinement.is_single_cyclic_graph(two_cycle_graph))

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
        self.assertTrue(UncertainRefinement.is_complex_cyclic_graph(two_cycle_graph))
        self.assertTrue(UncertainRefinement.is_complex_cyclic_graph(three_cycle_graph))

    def test_check_graph_type(self):
        # Create an acyclic graph (tree)
        tree_graph = nx.balanced_tree(2, 3)

        # Create a single cyclic graph (cycle)
        single_cycle_graph = nx.cycle_graph(5)

        # Create a graph with two cycles
        two_cycle_graph = nx.Graph()
        two_cycle_graph.add_edges_from([(1, 2), (2, 3), (3, 1), (3, 4), (4, 5), (5, 3)])

        # Test check_graph_type method
        self.assertEqual(UncertainRefinement.check_graph_type(tree_graph), "Acyclic")
        self.assertEqual(
            UncertainRefinement.check_graph_type(single_cycle_graph), "Single Cyclic"
        )
        self.assertEqual(
            UncertainRefinement.check_graph_type(two_cycle_graph), "Complex Cyclic"
        )


if __name__ == "__main__":
    unittest.main()
