# import unittest
# from syntemp.SynRule.auto_cat import AutoCat
# from syntemp.SynRule.rc_cluster import RCCluster
# from synkit.Graph.Feature.graph_descriptors import GraphDescriptor
# from synkit.IO.data_io import load_from_pickle


# class TestRCCluster(unittest.TestCase):

#     @classmethod
#     def setUpClass(cls):
#         # Load data once for all tests
#         cls.graphs = load_from_pickle("Data/Testcase/its_graph.pkl.gz")
#         for value in cls.graphs:
#             value["RC"] = value["GraphRules"][2]
#             value["ITS"] = value["ITSGraph"][2]
#             value = GraphDescriptor.get_descriptors(value)
#         cls.batch_1 = cls.graphs[:20]
#         cls.batch_2 = cls.graphs[20:40]
#         cls.batch_3 = cls.graphs[40:60]
#         cls.batch_4 = cls.graphs[60:80]
#         cls.batch_5 = cls.graphs[80:100]
#         cls.clusterer = RCCluster()

#     def test_fit(self):
#         """Test the fit method"""

#         cluster = RCCluster()
#         templates = None
#         result = []
#         for batch in [
#             self.batch_1,
#             self.batch_2,
#             self.batch_3,
#             self.batch_4,
#             self.batch_5,
#         ]:
#             if templates is None:
#                 templates = cluster.fit(
#                     self.batch_1, RC_col="RC", attribute="signature_rc"
#                 )
#                 result.extend(templates)
#             else:
#                 at = AutoCat(templates)
#                 batch, templates = at.fit(batch, attribute="signature_rc")
#                 result.extend(batch)
#         max_class = max([value["class"] for value in result])
#         self.assertEqual(max_class, 36)
#         self.assertEqual(len(result), len(self.graphs))


# if __name__ == "__main__":
#     unittest.main()
