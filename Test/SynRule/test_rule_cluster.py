import unittest
from syntemp.SynRule.rule_cluster import RuleCluster
from syntemp.SynUtils.utils import load_from_pickle


class TestRuleCluster(unittest.TestCase):

    def setUp(self) -> None:

        self.cluster = RuleCluster()

        self.data = load_from_pickle("Data/Testcase/its_graph.pkl.gz")
        self.data_rc = [value["GraphRules"] for value in self.data]

    def test_fit(self):
        cluster_indices, templates = self.cluster.fit(self.data_rc)
        assert any(isinstance(x, int) for x in cluster_indices)
        assert "Cluster_id" in templates[0].keys()
        assert "RC" in templates[0].keys()
        assert "Parent" in templates[0].keys()
        assert "Percentage" in templates[0].keys()


if __name__ == "__main__":
    unittest.main()
