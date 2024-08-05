from typing import List, Any, Dict, Tuple
import pandas as pd
import copy
from syntemp.SynRule.rules_extraction import RuleExtraction
from syntemp.SynRule.rule_cluster import RuleCluster
from syntemp.SynUtils.graph_utils import (
    add_child_ids,
    get_descriptors,
)
import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class HierarchicalClustering(RuleCluster):
    def __init__(
        self,
        node_label_names: List[str] = [
            "element",
            "charge",
        ],
        node_label_default: List[Any] = ["*", 0],
        edge_attribute: str = "order",
        max_radius: int = 3,
    ):
        """
        Initializes the HierarchicalClustering with customization options for node and
        edge matching functions.

        Parameters:
        - node_label_names (List[str]): Node attribute names for matching.
        - node_label_default (List[Any]): Default values for node attributes.
        - edge_attribute (str): Edge attribute name for matching.
        - max_radius (int): Maximum number of hierarchical levels.
        """
        super().__init__()
        self.radius = list(range(max_radius + 1))
        self.nodeLabelDefault = node_label_default
        self.nodeLabelNames = node_label_names
        self.edgeAttribute = edge_attribute

    @staticmethod
    def split_graphs_by_class_and_indices(
        graphs: List[Any], class_labels: List[int]
    ) -> Tuple[Dict[int, List[Any]], Dict[int, List[int]]]:
        """
        Splits a list of graphs and their indices into separate lists based on the
        provided class labels.

        Parameters:
        - graphs (List[Any]): The list of graphs to be split.
        - class_labels (List[int]): The list containing class labels corresponding to each
        graph.

        Returns:
        - Tuple[Dict[int, List[Any]], Dict[int, List[int]]]: A tuple containing two
        dictionaries:
            1. A dictionary where keys are class labels and values are lists of graphs
            corresponding to each class.
            2. A dictionary where keys are class labels and values are lists of indices
            corresponding to each class.
        """
        if len(graphs) != len(class_labels):
            raise ValueError(
                "The length of 'graphs' and 'class_labels' must be the same."
            )

        class_dict = {}
        index_dict = {}
        for index, (label, graph) in enumerate(zip(class_labels, graphs)):
            if label not in class_dict:
                class_dict[label] = []
                index_dict[label] = []
            class_dict[label].append(graph)
            index_dict[label].append(index)

        return class_dict, index_dict

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
        Processes a level in a graph by extracting rules from the input graphs and
        clustering them.

        Parameters:
        - its_graphs (List[Any]): A list of input graphs to process.
        - k (int): The number of nearest neighbors to consider in the k-NN algorithm for
        rule extraction.
        - nodeLabelNames (List[str]): A list of labels for the nodes in the graph.
        - nodeCountDefault (Any): The default value for node labels if no label is
        specified.
        - edgeAttribute (str): The attribute name of the edges used in rule extraction.
        - templates (Dict): A dictionary of templates used for clustering rules.
        - update_template (bool): A flag to determine whether to update the templates
        after clustering.

        Returns:
        - Tuple[Dict, Dict]: A tuple containing the mapping of graphs to clusters and the
        potentially updated templates.
        """
        if k > 0:
            rc_graphs = [
                RuleExtraction.extract_reaction_rules(*value, extend=True, n_knn=k)
                for value in its_graphs
            ]
        else:
            rc_graphs = [
                RuleExtraction.extract_reaction_rules(*value, extend=False)
                for value in its_graphs
            ]

        # Fit the rule clusters with the extracted graphs and templates
        cluster_indices, templates = RuleCluster(
            node_label_names=nodeLabelNames,
            node_label_default=nodeLabelDefault,
            edge_attribute=edgeAttribute,
        ).fit(rc_graphs, templates, update_template)

        return cluster_indices, templates

    def process_child_level(
        self,
        parent_graphs,
        parent_cluster_indices,
        node_label_names,
        radius=1,
        nodeLabelDefault=["*", 0],
        edgeAttribute="order",
        templates=None,
        update_template=False,
    ):
        """
        Process graphs by clusters, updating templates and indices based on the specified
        node label names.

        Parameters:
        - graphs (list): A list of graph structures to be processed.
        - indices (list): A list of indices representing cluster identifications.
        - node_label_names (list): A list of node label names used in the graph
        processing.

        Returns:
        - tuple:
            - templates (list): A list of template dictionaries generated during
            processing.
            - cluster_indices_all (list): Updated list of all cluster indices after
            processing.
        """
        graph_dict, index_dict = self.split_graphs_by_class_and_indices(
            parent_graphs, parent_cluster_indices
        )
        templates = []
        cluster_indices_all = ["a"] * len(parent_graphs)
        max_index_template = 0

        for key, value in graph_dict.items():
            cluster_indices_batch, new_templates = self.process_level(
                value,
                radius,
                nodeLabelNames=node_label_names,
                nodeLabelDefault=nodeLabelDefault,
                edgeAttribute=edgeAttribute,
                templates=None,
                update_template=update_template,
            )

            cluster_indices_batch = [
                i + max_index_template for i in cluster_indices_batch
            ]
            new_templates = [
                {
                    "Cluster_id": (
                        template["Cluster_id"] + max_index_template
                        if key is not None
                        else None
                    ),
                    "RC": template["RC"],
                    "Parent": key,
                    "Percentage": template["Percentage"],
                }
                for template in new_templates
            ]
            max_index_template += len(new_templates)
            templates.extend(new_templates)

            for i, j in enumerate(index_dict[key]):
                if key is not None:
                    cluster_indices_all[j] = cluster_indices_batch[i]
                else:
                    cluster_indices_all[j] = None

        return cluster_indices_all, templates

    def fit(
        self,
        original_reaction_dicts: List[Dict[str, Any]],
        its_column: str = "ITSGraph",
        templates: List[Dict] = None,
        update_template: bool = True,
    ) -> List[Dict[str, Any]]:
        """
        Fit the hierarchical clustering model to the data.

        Parameters:
        - original_reaction_dicts (List[Dict[str, Any]]): List of reaction dictionaries.
        - its_column (str): Column name for the ITS graph.

        Returns:
        - List[Dict[str, Any]]: Updated reaction dictionaries with clustering information.
        """
        try:
            reaction_dicts = copy.deepcopy(original_reaction_dicts)
            its_graphs = [value[its_column] for value in reaction_dicts]

            logging.info("Processing with templates")
            logging.info("Parent level")
            cluster_indices_0, templates_0 = self.process_level(
                its_graphs,
                0,
                self.nodeLabelNames,
                self.nodeLabelDefault,
                self.edgeAttribute,
                templates,
                update_template,
            )
            cluster_indices = [cluster_indices_0]
            templates = [get_descriptors(templates_0)]

            parent_cluster_indices = cluster_indices_0
            for k in self.radius:
                if k > 0:
                    logging.info(f"Child level with radius {k}")
                    cluster_indices_k, templates_k = self.process_child_level(
                        its_graphs,
                        parent_cluster_indices,
                        self.nodeLabelNames,
                        k,
                        self.nodeLabelDefault,
                        self.edgeAttribute,
                        None,
                        update_template,
                    )
                    cluster_indices.append(cluster_indices_k)
                    templates.append(templates_k)
                    parent_cluster_indices = cluster_indices_k

            cluster_df = pd.DataFrame(
                {f"Cluster_R{k}": idx for k, idx in zip(self.radius, cluster_indices)}
            ).to_dict("records")
            for key, value in enumerate(reaction_dicts):
                value.update(cluster_df[key])
            reaction_dicts = get_descriptors(
                reaction_dicts, reaction_centers="GraphRules"
            )

            hier_templates = add_child_ids(templates)

            return reaction_dicts, templates, hier_templates

        except Exception as e:
            print(f"An error occurred: {e}")
