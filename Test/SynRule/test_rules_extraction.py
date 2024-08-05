import unittest
import networkx as nx
from syntemp.SynRule.rules_extraction import RuleExtraction


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
        result_nodes = RuleExtraction.find_unequal_order_edges(G)
        self.assertCountEqual(result_nodes, expected_nodes)

    def test_find_nearest_neighbors(self):
        # Create a test graph
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5)])

        # Test find_nearest_neighbors method
        center_nodes = [1]
        expected_neighbors = {1, 2, 3}  # Assuming n_knn = 2
        result_neighbors = RuleExtraction.find_nearest_neighbors(
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
        subgraph = RuleExtraction.extract_subgraph(G, node_indices)
        self.assertTrue(nx.is_isomorphic(subgraph, nx.Graph([(2, 3)])))


if __name__ == "__main__":
    unittest.main()
