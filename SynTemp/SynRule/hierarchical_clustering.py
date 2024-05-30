from typing import List, Any, Dict
import pandas as pd
import copy
from joblib import Parallel, delayed
from SynTemp.SynRule.rules_extraction import RuleExtraction
from SynTemp.SynRule.rule_cluster import RuleCluster
from SynTemp.SynUtils.graph_utils import check_graph_type, get_cycle_member_rings


class HierarchicalClustering(RuleCluster):
    def __init__(
        self,
        node_label_names: List[str] = [
            "element",
            "aromatic",
            "hcount",
            "charge",
            "typesGH",
        ],
        node_label_default: List[Any] = ["*", False, 0, 0, ()],
        edge_attribute: str = "order",
        max_radius: int = 3,
    ):
        """
        Initializes the HierarchicalClustering with customization options for node and edge matching functions.

        Parameters:
            node_label_names (List[str]): Node attribute names for matching.
            node_label_default (List[Any]): Default values for node attributes.
            edge_attribute (str): Edge attribute name for matching.
            max_radius (int): Maximum number of hierarchical levels.
        """
        super().__init__()
        self.radius = list(range(max_radius + 1))
        self.nodeLabelDefault = node_label_default
        self.nodeLabelNames = node_label_names
        self.edgeAttribute = edge_attribute

    def fit(
        self,
        original_reaction_dicts: List[Dict[str, Any]],
        its_column: str = "ITSGraph",
        n_jobs: int = 1,
        verbose: int = 0
    ) -> List[Dict[str, Any]]:
        """
        Fit the hierarchical clustering model to the data.

        Parameters:
            original_reaction_dicts (List[Dict[str, Any]]): List of reaction dictionaries.
            its_column (str): Column name for the ITS graph.

        Returns:
            List[Dict[str, Any]]: Updated reaction dictionaries with clustering information.
        """
        n_jobs = max(n_jobs, len(self.radius))
        try:
            reaction_dicts = copy.deepcopy(original_reaction_dicts)
            its_graphs = [value[its_column] for value in reaction_dicts]

            def process_level(k):
                rc_graphs = [
                    RuleExtraction.extract_reaction_rules(*value, extend=True, n_knn=k)[
                        2
                    ]
                    for value in its_graphs
                ]
                return RuleCluster(
                    node_label_names=self.nodeLabelNames,
                    node_label_default=self.nodeLabelDefault,
                    edge_attribute=self.edgeAttribute,
                ).get_cluster_indices(rc_graphs)

            cluster_indices = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(process_level)(k) for k in self.radius
            )
            cluster_df = pd.DataFrame(
                {f"Cluster_R{k}": idx for k, idx in zip(self.radius, cluster_indices)}
            ).to_dict("records")

            for key, value in enumerate(reaction_dicts):
                value.update(cluster_df[key])
                value["Reaction Type"] = check_graph_type(value["GraphRules"][2])
                value["Rings"] = get_cycle_member_rings(value["GraphRules"][2])
            return reaction_dicts
        except Exception as e:
            print(f"An error occurred: {e}")
            raise
