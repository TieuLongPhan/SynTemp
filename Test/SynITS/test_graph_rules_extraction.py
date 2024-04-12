import unittest
import networkx as nx
<<<<<<< HEAD
import unittest
import sys
from pathlib import Path

root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
import time
=======
>>>>>>> bd12d73c8bed4cee5fafcc379f9f05958c1c21e0
from SynTemp.SynITS.graph_rules_extraction import GraphRuleExtraction


class TestGraphRuleExtraction(unittest.TestCase):

    def test_find_unequal_order_edges(self):
        # Create a test graph
        G = nx.Graph()
        G.add_edge(1, 2, order=(1, 2))  # change
        G.add_edge(2, 3, order=(1, 1))  # not change
        G.add_edge(3, 4, order=(2, 3))  # change

        # Expected nodes with unequal order edges
        expected_nodes = [1, 2, 3, 4]

        # Test find_unequal_order_edges method
        result_nodes = GraphRuleExtraction.find_unequal_order_edges(G)
        self.assertCountEqual(result_nodes, expected_nodes)

    def test_find_nearest_neighbors(self):
        # Create a test graph
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5)])

        # Test find_nearest_neighbors method
        center_nodes = [1]
        expected_neighbors = {1, 2, 3}  # Assuming n_knn = 2
        result_neighbors = GraphRuleExtraction.find_nearest_neighbors(
            G, center_nodes, n_knn=2
        )
        self.assertEqual(result_neighbors, expected_neighbors)

    def test_extract_subgraph(self):
        # Create a test graph
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3), (3, 4)])

        # Node indices for subgraph extraction
        node_indices = [2, 3]

        # Extract subgraph and test
        subgraph = GraphRuleExtraction.extract_subgraph(G, node_indices)
        self.assertTrue(nx.is_isomorphic(subgraph, nx.Graph([(2, 3)])))

    # def test_rules_extraction(self):
    #     # This test would depend on the structure of your 'reaction_dict'
    #     # and would require setting up a test case with a known outcome.
    #     pass

    # def test_process_rules_extraction(self):
    #     # Similar to test_rules_extraction, this would require setting up
    #     # multiple reaction_dicts and verifying the parallel processing outcomes.
    #     pass


if __name__ == "__main__":
    unittest.main()
