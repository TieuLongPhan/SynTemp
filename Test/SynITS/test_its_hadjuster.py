# import unittest
# import networkx as nx
# from SynTemp.SynITS.its_hadjuster import ITSHAdjuster


# class TestITSHAdjuster(unittest.TestCase):

#     def create_mock_graph(self, hcounts: dict) -> nx.Graph:
#         """Utility function to create a mock graph with specified
#         hydrogen counts for nodes."""
#         graph = nx.Graph()
#         for node_id, hcount in hcounts.items():
#             graph.add_node(node_id, hcount=hcount)
#         return graph

#     def test_check_hcount_change(self):
#         # Mock reactant and product graphs with specified hydrogen counts
#         react_graph = self.create_mock_graph({1: 1, 2: 2})
#         prod_graph = self.create_mock_graph({1: 0, 2: 3})

#         # Expected: one hydrogen formation (node 1) and one hydrogen break (node 2)
#         max_hydrogen_change = ITSHAdjuster.check_hcount_change(react_graph, prod_graph)
#         self.assertEqual(max_hydrogen_change, 1)

#     def test_add_hydrogen_nodes(self):
#         # Mock reactant and product graphs with specified hydrogen counts
#         react_graph = self.create_mock_graph({1: 1})
#         prod_graph = self.create_mock_graph({1: 0})

#         # Add hydrogen nodes to reactant and product graphs
#         updated_react_graph, _ = ITSHAdjuster.add_hydrogen_nodes(
#             react_graph, prod_graph
#         )

#         # Verify that hydrogen nodes have been added correctly
#         self.assertIn(
#             max(updated_react_graph.nodes), updated_react_graph.nodes
#         )  # Hydrogen node added to reactant graph
#         self.assertEqual(
#             updated_react_graph.nodes[max(updated_react_graph.nodes)]["element"], "H"
# )  # Check element of added node

# def test_add_hydrogen_nodes_multiple(self):
#     # Mock reactant and product graphs with specified hydrogen counts
#     react_graph = self.create_mock_graph({1: 2, 2: 1})
#     prod_graph = self.create_mock_graph({1: 0, 2: 2})

#     # Generate updated graph pairs with multiple hydrogen nodes added
#     updated_graph_pairs = ITSHAdjuster.add_hydrogen_nodes_multiple(
#         react_graph, prod_graph
#     )

#     # Verify that multiple updated graph pairs are generated
#     self.assertTrue(len(updated_graph_pairs) > 1)  # Multiple permutations generated
#     for react_graph, prod_graph in updated_graph_pairs:
#         self.assertIn(
#             max(react_graph.nodes), react_graph.nodes
#         )  # Hydrogen node added to reactant graph
#         self.assertIn(
#             max(prod_graph.nodes), prod_graph.nodes
#         )  # Hydrogen node added to product graph


# if __name__ == "__main__":
#     unittest.main()
