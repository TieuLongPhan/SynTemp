from typing import List, Any, Dict, Tuple
import pandas as pd
import copy
from joblib import Parallel, delayed
from SynTemp.SynRule.rules_extraction import RuleExtraction
from SynTemp.SynRule.rule_cluster import RuleCluster
from SynTemp.SynUtils.graph_utils import check_graph_type, get_cycle_member_rings
import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


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

    @staticmethod
    def process_level(
        its_graphs: List[Any],
        k: int,
        nodeLabelNames: List[str],
        nodeLabelDefault: Any,
        edgeAttribute: str,
        templates: Dict,
        update_template: bool,
    ) -> Tuple[Dict, Dict]:
        """
        Processes a level in a graph by extracting rules from the input graphs and clustering them.

        Args:
            its_graphs (List[Any]): A list of input graphs to process.
            k (int): The number of nearest neighbors to consider in the k-NN algorithm for rule extraction.
            nodeLabelNames (List[str]): A list of labels for the nodes in the graph.
            nodeCountDefault (Any): The default value for node labels if no label is specified.
            edgeAttribute (str): The attribute name of the edges used in rule extraction.
            templates (Dict): A dictionary of templates used for clustering rules.
            update_template (bool): A flag to determine whether to update the templates after clustering.

        Returns:
            Tuple[Dict, Dict]: A tuple containing the mapping of graphs to clusters and the potentially updated templates.
        """
        logging.info(f"Processing templates with {k}:")
        rc_graphs = [
            RuleExtraction.extract_reaction_rules(*value, extend=True, n_knn=k)[2]
            for value in its_graphs
        ]

        # Fit the rule clusters with the extracted graphs and templates
        cluster_indices, templates = RuleCluster(
            node_label_names=nodeLabelNames,
            node_label_default=nodeLabelDefault,
            edge_attribute=edgeAttribute,
        ).fit(rc_graphs, templates, update_template)

        return cluster_indices, templates

    def fit(
        self,
        original_reaction_dicts: List[Dict[str, Any]],
        its_column: str = "ITSGraph",
        templates: List[Dict] = None,
        update_template: bool = True,
        root_sample: int = 100,
    ) -> List[Dict[str, Any]]:
        """
        Fit the hierarchical clustering model to the data.

        Parameters:
            original_reaction_dicts (List[Dict[str, Any]]): List of reaction dictionaries.
            its_column (str): Column name for the ITS graph.

        Returns:
            List[Dict[str, Any]]: Updated reaction dictionaries with clustering information.
        """
        try:
            reaction_dicts = copy.deepcopy(original_reaction_dicts)
            its_graphs = [value[its_column] for value in reaction_dicts]
            if templates:
                logging.info("Processing with templates")
                results = [
                    self.process_level(
                        its_graphs,
                        k,
                        self.nodeLabelNames,
                        self.nodeLabelDefault,
                        self.edgeAttribute,
                        templates,
                        update_template,
                    )
                    for k in self.radius
                ]

                cluster_indices = [value[0] for value in results]
                templates = [value[1] for value in results]
            else:
                logging.info("Processing without templates")
                root_length = min(root_sample, len(its_graphs))

                its_root = its_graphs[:root_length]
                its_left = its_graphs[root_length:]

                logging.info(f"Processing {root_length} data to get templates")
                results = [
                    self.process_level(
                        its_root,
                        k,
                        self.nodeLabelNames,
                        self.nodeLabelDefault,
                        self.edgeAttribute,
                        None,
                        update_template,
                    )
                    for k in self.radius
                ]

                cluster_indices_root = [value[0] for value in results]
                templates_root = [value[1] for value in results]

                logging.info("Processing other data with new templates")
                results_left = [
                    self.process_level(
                        its_left,
                        k,
                        self.nodeLabelNames,
                        self.nodeLabelDefault,
                        self.edgeAttribute,
                        templates_root[k],
                        update_template,
                    )
                    for k in self.radius
                ]
                cluster_indices_left = [value[0] for value in results_left]
                templates = [value[1] for value in results_left]

                cluster_indices = [
                    cluster_indices_root[key] + cluster_indices_left[key]
                    for key, _ in enumerate(cluster_indices_root)
                ]

            cluster_df = pd.DataFrame(
                {f"Cluster_R{k}": idx for k, idx in zip(self.radius, cluster_indices)}
            ).to_dict("records")

            for key, value in enumerate(reaction_dicts):
                value.update(cluster_df[key])
                value["Reaction Type"] = check_graph_type(value["GraphRules"][2])
                value["Rings"] = get_cycle_member_rings(value["GraphRules"][2])

            return reaction_dicts, templates

        except Exception as e:
            print(f"An error occurred: {e}")
