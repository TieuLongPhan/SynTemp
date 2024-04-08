import unittest
import networkx as nx
import unittest
import sys
from pathlib import Path
root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
import time
from SynTemp.SynMØD.MØD_rules import MØDRules
\
class TestMODRules(unittest.TestCase):

    def test_convert_graph_to_gml(self):
        # Create a simple graph to test
        graph = nx.Graph()
        graph.add_node(1, element='H')
        graph.add_node(2, element='O')
        graph.add_edge(1, 2, order=1)

        # Generate GML string
        gml_str = MØDRules.convert_graph_to_gml(graph, "left")
        
        # Check if the GML string is formatted correctly
        self.assertIn("node [ id 1 label \"H\" ]", gml_str)
        self.assertIn("node [ id 2 label \"O\" ]", gml_str)
        self.assertIn("edge [ source 1 target 2 label \"-\" ]", gml_str)

    # def test_RulesGrammar(self):
    #     # Create simple graphs for L, R, and K
    #     L = nx.Graph()
    #     R = nx.Graph()
    #     K = nx.Graph()
    #     L.add_node(1, element='H')
    #     R.add_node(1, element='H')
    #     K.add_node(1, element='H')

    #     # Generate GML string for the rule
    #     gml_str = MØDRules.RulesGrammar(L, R, K, "test_rule")

    #     # Check if the rule is formatted correctly
    #     self.assertIn("ruleID \"test_rule\"", gml_str)
    #     self.assertIn("left [", gml_str)
    #     self.assertIn("right [", gml_str)
    #     self.assertIn("context [", gml_str)

    # def test_process_graph_rules(self):
    #     # Define a rule in dictionary format
    #     graph_rule = {
    #         'R-id': 'rule1',
    #         'graph_rules': (nx.Graph(), nx.Graph(), nx.Graph())
    #     }
    #     graph_rule['graph_rules'][0].add_node(1, element='H')  # Adding a node to L graph

    #     # Process the graph rule
    #     result = MØDRules.process_graph_rules(graph_rule)

    #     # Check if the result is a dictionary with the correct rule ID
    #     self.assertIsInstance(result, str)
    #     self.assertIn("ruleID \"rule1\"", result)

    # def test_auto_extraction(self):
    #     # Define a list of graph rules
    #     data_dicts = [{
    #         'R-id': 'rule1',
    #         'graph_rules': (nx.Graph(), nx.Graph(), nx.Graph())
    #     }]
    #     data_dicts[0]['graph_rules'][0].add_node(1, element='H')  # Adding a node to L graph of the first rule

    #     # Test auto_extraction
    #     results = MØDRules.auto_extraction(data_dicts, save_path=None)  # Not testing file saving here

    #     # Check if the results are correct
    #     self.assertIsInstance(results, list)
    #     self.assertEqual(len(results), 1)
    #     self.assertIn("ruleID \"rule1\"", results[0])

if __name__ == '__main__':
    unittest.main()
