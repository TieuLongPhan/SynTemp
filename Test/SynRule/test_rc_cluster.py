import time
import unittest
from syntemp.SynRule.rc_cluster import RCCluster
from synutility.SynIO.data_type import load_from_pickle
from synutility.SynGraph.Descriptor.graph_descriptors import GraphDescriptor


class TestRCCluster(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Load data once for all tests
        cls.graphs = load_from_pickle("Data/Testcase/its_graph.pkl.gz")
        for value in cls.graphs:
            value["RC"] = value["GraphRules"][2]
            value["ITS"] = value["ITSGraph"][2]
            value = GraphDescriptor.get_descriptors(value)
        cls.clusterer = RCCluster()

    def test_initialization(self):
        """Test the initialization and configuration of the RCCluster."""
        self.assertIsInstance(self.clusterer.nodeLabelNames, list)
        self.assertEqual(self.clusterer.edgeAttribute, "order")
        self.assertEqual(
            len(self.clusterer.nodeLabelNames), len(self.clusterer.nodeLabelDefault)
        )

    def test_auto_cluster(self):
        """Test the auto_cluster method functionality."""
        rc = [value["RC"] for value in self.graphs]
        cycles = [value["cycle"] for value in self.graphs]
        signature = [value["signature_rc"] for value in self.graphs]
        atom_count = [value["atom_count"] for value in self.graphs]
        for att in [None, cycles, signature, atom_count]:
            clusters, graph_to_cluster = self.clusterer.auto_cluster(
                rc,
                att,
                nodeMatch=self.clusterer.nodeMatch,
                edgeMatch=self.clusterer.edgeMatch,
            )
            self.assertIsInstance(clusters, list)
            self.assertIsInstance(graph_to_cluster, dict)
            self.assertEqual(len(clusters), 37)

    def test_auto_cluster_wrong_isomorphism(self):
        rc = [value["RC"] for value in self.graphs]
        cycles = [value["cycle"] for value in self.graphs]
        signature = [value["signature_rc"] for value in self.graphs]
        atom_count = [value["atom_count"] for value in self.graphs]

        # cluster all
        clusters, _ = self.clusterer.auto_cluster(
            rc, None, nodeMatch=None, edgeMatch=None
        )
        self.assertEqual(len(clusters), 10)  # wrong value

        # cluster with cycle
        clusters, _ = self.clusterer.auto_cluster(
            rc, cycles, nodeMatch=None, edgeMatch=None
        )
        self.assertEqual(len(clusters), 10)  # wrong value

        # cluster with atom_count
        clusters, _ = self.clusterer.auto_cluster(
            rc, atom_count, nodeMatch=None, edgeMatch=None
        )
        self.assertEqual(len(clusters), 32)  # wrong value but almost correct

        # cluster with signature
        clusters, _ = self.clusterer.auto_cluster(
            rc, signature, nodeMatch=None, edgeMatch=None
        )
        self.assertEqual(len(clusters), 37)  # correct by some magic. No proof for this

    def test_fit(self):
        """Test the fit method to ensure it correctly updates data entries with cluster indices."""

        clustered_data = self.clusterer.fit(self.graphs, attribute="atom_count")
        max_class = 0
        for item in clustered_data:
            print(item["class"])
            max_class = item["class"] if item["class"] >= max_class else max_class
            # print(max_class)
            self.assertIn("class", item)
        self.assertEqual(max_class, 36)  # 37 classes start from 0 so max is 36

    def test_fit_time_compare(self):
        attributes = {
            "None": None,
            "Cycles": "cycle",
            "Signature": "signature_rc",
            "Atom_count": "atom_count",
        }

        results = {}

        for name, attr in attributes.items():
            start_time = time.time()
            clustered_data = self.clusterer.fit(self.graphs, attribute=attr)
            elapsed_time = time.time() - start_time

            # Optionally print out class information or verify correctness
            max_class = max(item["class"] for item in clustered_data if "class" in item)

            results[name] = elapsed_time

            # Basic verification that 'class' is assigned and max class is as expected
            self.assertTrue(all("class" in item for item in clustered_data))
            self.assertEqual(
                max_class, 36
            )  # Ensure the maximum class index is as expected

        # Compare results to check which attribute took the least/most time
        min_time_attr = min(results, key=results.get)
        max_time_attr = max(results, key=results.get)
        self.assertIn(min_time_attr, ["atom_count", "Signature"])
        self.assertIn(max_time_attr, ["None", "Cycles"])


if __name__ == "__main__":
    unittest.main()
