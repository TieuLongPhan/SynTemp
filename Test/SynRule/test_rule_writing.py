import unittest
import networkx as nx
from syntemp.SynUtils.utils import load_from_pickle
from syntemp.SynRule.rule_writing import RuleWriting


class TestRuleWriting(unittest.TestCase):

    def setUp(self) -> None:
        self.data = load_from_pickle("Data/Testcase/templates.pkl.gz")[0]

    def test_charge_to_string(self):
        self.assertEqual(RuleWriting.charge_to_string(3), "3+")
        self.assertEqual(RuleWriting.charge_to_string(-2), "2-")
        self.assertEqual(RuleWriting.charge_to_string(0), "")

    def test_convert_graph_to_gml_context(self):
        G = nx.Graph()
        G.add_node(1, element="C")
        G.add_node(2, element="H")
        changed_node_ids = [2]
        gml_str = RuleWriting.convert_graph_to_gml(G, "context", changed_node_ids)
        expected_str = '   context [\n      node [ id 1 label "C" ]\n   ]\n'
        self.assertEqual(gml_str, expected_str)

    def test_convert_graph_to_gml_left_right(self):
        G = nx.Graph()
        G.add_node(1, element="C", charge=1)
        G.add_node(2, element="H", charge=0)
        G.add_edge(1, 2, order=2)
        changed_node_ids = [1]
        gml_str = RuleWriting.convert_graph_to_gml(G, "left", changed_node_ids)
        expected_str = (
            '   left [\n      edge [ source 1 target 2 label "=" ]'
            + '\n      node [ id 1 label "C+" ]\n   ]\n'
        )
        self.assertEqual(gml_str, expected_str)

    def test_rules_grammar(self):
        L = self.data[0]["RC"][0]
        R = self.data[0]["RC"][1]
        K = self.data[0]["RC"][2]
        changed_node_ids = RuleWriting.find_changed_nodes(L, R, ["charge"])
        rule_name = "test_rule"
        gml_str = RuleWriting.RulesGrammar(L, R, K, rule_name, changed_node_ids)
        expected_str = (
            "rule [\n"
            '   ruleID "test_rule"\n'
            "   left [\n"
            '      edge [ source 21 target 4 label "-" ]\n'
            '      edge [ source 3 target 20 label "-" ]\n'
            "   ]\n"
            "   context [\n"
            '      node [ id 21 label "H" ]\n'
            '      node [ id 3 label "C" ]\n'
            '      node [ id 4 label "N" ]\n'
            '      node [ id 20 label "Br" ]\n'
            "   ]\n"
            "   right [\n"
            '      edge [ source 21 target 20 label "-" ]\n'
            '      edge [ source 3 target 4 label "-" ]\n'
            "   ]\n"
            "]"
        )
        self.assertEqual(gml_str, expected_str)

    def test_find_changed_nodes(self):
        G1 = nx.Graph()
        G1.add_node(1, element="C", charge=0)
        G2 = nx.Graph()
        G2.add_node(1, element="C", charge=1)
        changed_nodes = RuleWriting.find_changed_nodes(G1, G2, ["charge"])
        self.assertEqual(changed_nodes, [1])

    def test_process_graph_rules(self):
        L = self.data[0]["RC"][0]
        R = self.data[0]["RC"][1]
        K = self.data[0]["RC"][2]
        graph_rules = {"R-id": "test_rule", "GraphRules": (L, R, K)}
        gml_str = RuleWriting.process_graph_rules(graph_rules, reindex=True)
        print(gml_str)
        expected_str = (
            "rule [\n"
            '   ruleID "test_rule"\n'
            "   left [\n"
            '      edge [ source 1 target 3 label "-" ]\n'
            '      edge [ source 2 target 4 label "-" ]\n'
            "   ]\n"
            "   context [\n"
            '      node [ id 1 label "H" ]\n'
            '      node [ id 2 label "C" ]\n'
            '      node [ id 3 label "N" ]\n'
            '      node [ id 4 label "Br" ]\n'
            "   ]\n"
            "   right [\n"
            '      edge [ source 1 target 4 label "-" ]\n'
            '      edge [ source 2 target 3 label "-" ]\n'
            "   ]\n"
            "]"
        )
        self.assertEqual(gml_str, expected_str)

    def test_auto_extraction(self):
        write = RuleWriting.auto_extraction(
            self.data, id_column="Cluster_id", rule_column="RC", reindex=True
        )
        expected_str = (
            "rule [\n"
            '   ruleID "0"\n'
            "   left [\n"
            '      edge [ source 1 target 3 label "-" ]\n'
            '      edge [ source 2 target 4 label "-" ]\n'
            "   ]\n"
            "   context [\n"
            '      node [ id 1 label "H" ]\n'
            '      node [ id 2 label "C" ]\n'
            '      node [ id 3 label "N" ]\n'
            '      node [ id 4 label "Br" ]\n'
            "   ]\n"
            "   right [\n"
            '      edge [ source 1 target 4 label "-" ]\n'
            '      edge [ source 2 target 3 label "-" ]\n'
            "   ]\n"
            "]"
        )
        self.assertEqual(write[0], expected_str)


if __name__ == "__main__":
    unittest.main()
