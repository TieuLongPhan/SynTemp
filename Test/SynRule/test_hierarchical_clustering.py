import unittest
from syntemp.SynRule.hierarchical_clustering import HierarchicalClustering
from syntemp.SynUtils.utils import load_from_pickle


class TestRuleCluster(unittest.TestCase):

    def setUp(self) -> None:

        self.cluster = HierarchicalClustering(max_radius=3)

        self.data = load_from_pickle("Data/Testcase/its_graph.pkl.gz")

    def test_fit(self):
        reaction_dicts, templates, hier_templates = self.cluster.fit(self.data)

        # assert data
        assert "Cluster_R0" in reaction_dicts[0].keys()
        assert "Cluster_R1" in reaction_dicts[0].keys()
        assert "Cluster_R2" in reaction_dicts[0].keys()
        assert "Cluster_R3" in reaction_dicts[0].keys()
        self.assertEqual(reaction_dicts[0]["Reaction Type"], "Elementary")
        self.assertEqual(reaction_dicts[0]["Topo Type"], "Single Cyclic")
        self.assertEqual(reaction_dicts[0]["Rings"][0], 4)  # 0 1 2 3
        self.assertEqual(reaction_dicts[0]["Reaction Step"], 1)  # 0 1 2 3

        # assert templates
        self.assertEqual(len(templates), 4)  # 0 1 2 3
        assert "RC" in templates[0][0].keys()
        assert "Parent" in templates[0][0].keys()
        assert "Percentage" in templates[0][0].keys()

        # assert hier templates
        assert isinstance(hier_templates[0][0]["Cluster_id"], int)
        assert isinstance(hier_templates[0][0]["Percentage"], float)
        assert "Parent" in hier_templates[0][0].keys()
        assert "Child" in hier_templates[0][0].keys()


if __name__ == "__main__":
    unittest.main()
