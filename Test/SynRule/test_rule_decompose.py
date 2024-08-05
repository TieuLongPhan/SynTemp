import unittest
import networkx as nx
from syntemp.SynRule.rule_decompose import RuleDecompose


class TestRuleDecompose(unittest.TestCase):

    def setUp(self):
        # Setup basic graph structures for testing
        self.G = nx.Graph()
        self.G.add_nodes_from([1, 2, 3, 4], element="C", charge=0, aromatic=False)
        self.G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 1)], standard_order=1)

        self.H = nx.Graph()
        self.H.add_nodes_from([1, 2, 3], element="C", charge=0, aromatic=False)
        self.H.add_edges_from([(1, 2), (2, 3), (3, 1)], standard_order=1)

    def test_get_key_by_value(self):
        test_dict = {"a": 1, "b": 2, "c": 3}
        self.assertEqual(RuleDecompose.get_key_by_value(test_dict, 2), "b")
        self.assertIsNone(RuleDecompose.get_key_by_value(test_dict, 4))

    def test_is_connected(self):
        self.assertTrue(RuleDecompose.is_connected(self.G))
        self.G.add_node(5)
        self.assertFalse(RuleDecompose.is_connected(self.G))

    def test_remove_disconnected_part(self):
        self.G.add_node(5)  # Disconnected node
        connected_G = RuleDecompose.remove_disconnected_part(self.G)
        self.assertFalse(5 in connected_G.nodes())

    def test_node_match(self):
        node1_attrs = {"element": "C", "charge": 0, "aromatic": False}
        node2_attrs = {"element": "C", "charge": 0, "aromatic": False}
        node3_attrs = {"element": "O", "charge": 0, "aromatic": False}
        self.assertTrue(RuleDecompose.node_match(node1_attrs, node2_attrs))
        self.assertFalse(RuleDecompose.node_match(node1_attrs, node3_attrs))

    def test_edge_match(self):
        edge1_attrs = {"standard_order": 1}
        edge2_attrs = {"standard_order": -1}
        edge3_attrs = {"standard_order": 1}
        self.assertFalse(RuleDecompose.edge_match(edge1_attrs, edge2_attrs))
        self.assertTrue(RuleDecompose.edge_match(edge1_attrs, edge3_attrs))

    def test_find_maximum_common_subgraph(self):
        # Parent graph setup
        parent_graph = nx.Graph()
        parent_graph.add_nodes_from([1, 2, 3, 4], element="C", charge=0, aromatic=False)
        parent_graph.add_edges_from(
            [(1, 2), (2, 3), (3, 4), (4, 1), (3, 1)], standard_order=1
        )

        # Child graph that is a subgraph of the parent graph
        child_graph_sub = nx.Graph()
        child_graph_sub.add_nodes_from([1, 2, 3], element="C", charge=0, aromatic=False)
        child_graph_sub.add_edges_from([(1, 2), (2, 3), (3, 1)], standard_order=1)

        # Child graph that is not a subgraph but shares some nodes and edges
        child_graph_nonsub = nx.Graph()
        child_graph_nonsub.add_nodes_from(
            [2, 3, 5], element="C", charge=0, aromatic=True
        )
        child_graph_nonsub.add_edges_from([(2, 3), (3, 5)], standard_order=1)

        # Find isomorphisms
        isomorph_sub = RuleDecompose.find_maximum_common_subgraph(
            parent_graph,
            child_graph_sub,
            node_match=RuleDecompose.node_match,
            edge_match=RuleDecompose.edge_match,
        )
        isomorph_nonsub = RuleDecompose.find_maximum_common_subgraph(
            parent_graph,
            child_graph_nonsub,
            node_match=RuleDecompose.node_match,
            edge_match=RuleDecompose.edge_match,
        )

        # Assertions
        self.assertIsNotNone(isomorph_sub)
        self.assertEqual(len(isomorph_sub), 3)  # Expecting 3 nodes to match
        self.assertIsNone(isomorph_nonsub)  # No subgraph isomorphism should be found

    def test_remove_maximum_common_subgraph_edges(self):
        # Setup parent and child graph from previous test
        parent_graph = nx.Graph()
        parent_graph.add_nodes_from([1, 2, 3, 4], element="C", charge=0, aromatic=False)
        parent_graph.add_edges_from(
            [(1, 2), (2, 3), (3, 4), (4, 1), (3, 1)], standard_order=1
        )

        child_graph_sub = nx.Graph()
        child_graph_sub.add_nodes_from([1, 2, 3], element="C", charge=0, aromatic=False)
        child_graph_sub.add_edges_from([(1, 2), (2, 3), (3, 1)], standard_order=1)

        isomorph_sub = RuleDecompose.find_maximum_common_subgraph(
            parent_graph,
            child_graph_sub,
            node_match=RuleDecompose.node_match,
            edge_match=RuleDecompose.edge_match,
        )

        # Remove edges based on the isomorphism
        modified_graph = RuleDecompose.remove_maximum_common_subgraph_edges(
            parent_graph, child_graph_sub, isomorph_sub, "standard_order"
        )

        # Assertions
        self.assertFalse(modified_graph.has_edge(1, 2))  # Edge (1, 2) should be removed
        self.assertTrue(
            modified_graph.has_edge(3, 4)
        )  # Edge (3, 4) should remain as it's not part of the subgraph
        self.assertEqual(
            modified_graph.number_of_edges(), 2
        )  # Only two edges should remain


if __name__ == "__main__":
    unittest.main()
