import unittest
from unittest.mock import patch
import networkx as nx
from SynTemp.SynRule.rule_cluster import RuleCluster


class TestRuleCluster(unittest.TestCase):

    def test_cluster_graphs_single_graph(self):
        # Test clustering with a single graph
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3])
        graph.add_edges_from([(1, 2), (2, 3)])

        clusterer = RuleCluster()
        clusters = clusterer.cluster_graphs([graph])

        self.assertEqual(len(clusters), 1)
        self.assertEqual(list(clusters[0]), [0])

    def test_cluster_graphs_multiple_graphs(self):
        # Test clustering with multiple identical graphs
        graph1 = nx.Graph()
        graph1.add_nodes_from([1, 2, 3])
        graph1.add_edges_from([(1, 2), (2, 3)])

        graph2 = nx.Graph()
        graph2.add_nodes_from([1, 2, 3])
        graph2.add_edges_from([(1, 2), (2, 3)])

        clusterer = RuleCluster()
        clusters = clusterer.cluster_graphs([graph1, graph2])

        self.assertEqual(len(clusters), 1)
        self.assertEqual(list(clusters[0]), [0, 1])

    def test_get_cluster_indices(self):
        # Test getting cluster indices
        graph = nx.Graph()
        graph.add_nodes_from([1, 2, 3])
        graph.add_edges_from([(1, 2), (2, 3)])

        clusterer = RuleCluster()
        indices = clusterer.get_cluster_indices([graph])

        self.assertEqual(indices, [0])

    def test_process_rules_clustering(self):
        # Test processing rules clustering
        # Setting up a graph that would be deterministic for the functions 'check_graph_type' and 'get_cycle_member_rings'
        graph = nx.Graph([(1, 2), (2, 3)])
        reaction_dicts = [{"rules": ({}, {}, graph), "GraphRules": ({}, {}, graph)}]

        clusterer = RuleCluster()
        updated_dicts = clusterer.process_rules_clustering(reaction_dicts, "rules")

        self.assertEqual(len(updated_dicts), 1)
        self.assertIn("naive_cluster", updated_dicts[0])
        self.assertTrue(isinstance(updated_dicts[0]["naive_cluster"], int))


if __name__ == "__main__":
    unittest.main()
