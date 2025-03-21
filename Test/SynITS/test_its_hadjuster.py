import unittest
import networkx as nx
from copy import deepcopy
from synkit.IO.data_io import load_from_pickle
from syntemp.SynITS.its_hadjuster import ITSHAdjuster


class TestITSHAdjuster(unittest.TestCase):

    def setUp(self):
        """Setup before each test."""
        # Create sample graphs
        self.data = load_from_pickle("./Data/Testcase/hydrogen_test.pkl.gz")

    def test_process_single_graph_data_success(self):
        """Test the process_single_graph_data method."""
        processed_data = ITSHAdjuster.process_single_graph_data(
            self.data[0], "ITSGraph"
        )
        for value in processed_data["ITSGraph"]:
            self.assertTrue(isinstance(value, nx.Graph))
        for value in processed_data["GraphRules"]:
            self.assertTrue(isinstance(value, nx.Graph))

    def test_process_single_graph_data_fail(self):
        """Test the process_single_graph_data method."""
        processed_data = ITSHAdjuster.process_single_graph_data(
            self.data[16], "ITSGraph"
        )
        self.assertIsNone(processed_data["ITSGraph"])
        self.assertIsNone(processed_data["GraphRules"])

    def test_process_single_graph_data_empty_graph(self):
        """Test that an empty graph results in empty ITSGraph and GraphRules."""
        empty_graph_data = {
            "ITSGraph": [None, None, None],
            "GraphRules": [None, None, None],
        }

        processed_data = ITSHAdjuster.process_single_graph_data(
            empty_graph_data, "ITSGraph"
        )

        # Ensure the result is None or empty as expected for an empty graph
        self.assertIsNone(processed_data["ITSGraph"])
        self.assertIsNone(processed_data["GraphRules"])

    def test_process_single_graph_data_safe(self):
        """Test the process_single_graph_data method."""
        processed_data = ITSHAdjuster.process_single_graph_data_safe(
            self.data[0], "ITSGraph", job_timeout=0.0001
        )
        self.assertIsNone(processed_data["ITSGraph"])
        self.assertIsNone(processed_data["GraphRules"])

    def test_process_graph_data_parallel(self):
        """Test the process_graph_data_parallel method."""
        result = ITSHAdjuster().process_graph_data_parallel(
            self.data, "ITSGraph", n_jobs=1, verbose=0, get_priority_graph=True
        )
        result = [value for value in result if value["ITSGraph"]]
        # Check if the result matches the input data structure
        self.assertEqual(len(result), 48)

    def test_process_graph_data_parallel_safe(self):
        """Test the process_graph_data_parallel method."""
        result = ITSHAdjuster().process_graph_data_parallel(
            self.data,
            "ITSGraph",
            n_jobs=1,
            verbose=0,
            get_priority_graph=True,
            safe=True,
            job_timeout=0.0001,  # lower timeout will fail all process
        )
        result = [value for value in result if value["ITSGraph"]]
        # Check if the result matches the input data structure
        self.assertEqual(len(result), 0)

    def test_process_multiple_hydrogens(self):
        """Test the process_multiple_hydrogens method."""
        graphs = deepcopy(self.data[0])
        react_graph, prod_graph, _ = graphs["ITSGraph"]

        result = ITSHAdjuster.process_multiple_hydrogens(
            graphs, react_graph, prod_graph, ignore_aromaticity=False, balance_its=True
        )

        for value in result["ITSGraph"]:
            self.assertTrue(isinstance(value, nx.Graph))
        for value in result["GraphRules"]:
            self.assertTrue(isinstance(value, nx.Graph))


if __name__ == "__main__":
    unittest.main()
