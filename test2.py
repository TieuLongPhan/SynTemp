from SynTemp.SynUtils.utils import load_from_pickle, save_to_pickle
from SynTemp.SynRule.hierarchical_clustering import HierarchicalClustering
data = load_from_pickle('./Data/test..pkl.gz')

node_label_names = ["element", "charge"]
hcl = HierarchicalClustering(node_label_names=node_label_names,
                            node_label_default=["*", 0],
                            edge_attribute="order",
                            max_radius=3)

reaction_dicts, templates = hcl.fit(data, 'ITSGraph', templates=None, update_template=True, root_sample=100)
save_to_pickle(reaction_dicts, './cluster.pkl.gz')
save_to_pickle(templates, './temp.pkl.gz')