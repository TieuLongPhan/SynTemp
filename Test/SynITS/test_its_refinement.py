# import unittest
# import importlib.resources
# from SynTemp.SynUtils.utils import load_from_pickle
# from SynTemp.SynITS.its_refinement import ITSRefinement


# class TestUncertainRefinement(unittest.TestCase):
#     def setUp(self) -> None:
#         data_path = importlib.resources.files("Data").joinpath(
#             "DPO/USPTO_50K/Hydrogen/USPTO_50K_its_incorrect.pkl.gz"
#         )

#         self.data = load_from_pickle(data_path)[:5]
#         print(self.data)

#     def test_process_and_check_graph(self):
#         _, type_of_graph = ITSRefinement.process_and_check_graph(
#             self.data[0], "local_mapper"
#         )
#         self.assertEqual(type_of_graph, "Single Cyclic")

#         _, type_of_graph = ITSRefinement.process_and_check_graph(
#             self.data[0], "rxn_mapper"
#         )
#         self.assertEqual(type_of_graph, "Complex Cyclic")

#         _, type_of_graph = ITSRefinement.process_and_check_graph(
#             self.data[1], "graphormer"
#         )
#         self.assertEqual(type_of_graph, "None")

#         _, type_of_graph = ITSRefinement.process_and_check_graph(
#             self.data[82], "rxn_mapper"
#         )  # Reaction center has no nodes
#         self.assertEqual(type_of_graph, "None")

#     def test_process_dict(self):
#         data = self.data[0]
#         processed_data = ITSRefinement.process_dict(data, "ITSGraph")
#         self.assertTrue("ITSGraph" and "GraphRules" in processed_data)

#         _, type_of_graph = ITSRefinement.process_and_check_graph(
#             processed_data, "ITSGraph"
#         )
#         self.assertEqual(type_of_graph, "Single Cyclic")

#     def test_process_graphs_in_parallel(self):
#         processed_data = ITSRefinement.process_graphs_in_parallel(self.data[0:100])
#         self.assertTrue("ITSGraph" and "GraphRules" in processed_data[0])


# if __name__ == "__main__":
#     unittest.main()
