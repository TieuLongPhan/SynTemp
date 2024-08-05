import networkx as nx
from typing import List, Set, Dict, Any, Callable, Optional, Tuple
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match
from operator import eq
from syntemp.SynUtils.utils import create_unique_value_dict


class RuleCluster:
    def __init__(
        self,
        node_label_names: List[str] = [
            "element",
            "aromatic",
            "charge",
        ],
        node_label_default: List[Any] = ["*", False, 0],
        edge_attribute: str = "order",
    ):
        """
        Initializes the NaiveClusterer with customization options for node and edge
        matching functions. This class is designed to facilitate clustering of graph nodes
        and edges based on specified attributes and their matching criteria.

        Parameters:
        - node_label_names (List[str]): A list of node attribute names to be considered
        for matching. Each attribute name corresponds to a property of the nodes in the
        graph. Default is ["element", "aromatic", "charge"].
        - node_label_default (List[Any]): Default values for each of the node attributes
        specified in `node_label_names`. These are used where node attributes are missing.
        The length and order of this list should match `node_label_names`.
        Default is ["*", False, 0].
        - edge_attribute (str): The name of the edge attribute to consider for matching
        edges. This attribute is used to assess edge similarity. Default is "order".

        Raises:
        - ValueError: If the lengths of `node_label_names` and `node_label_default` do not
        match.
        """

        self.nodeLabelNames: List[str] = node_label_names
        self.edgeAttribute: str = edge_attribute
        self.nodeLabelDefault: List[Any] = node_label_default
        self.nodeLabelOperator: List[Callable[[Any, Any], bool]] = [
            eq for _ in node_label_names
        ]
        self.nodeMatch: Callable = generic_node_match(
            self.nodeLabelNames, self.nodeLabelDefault, self.nodeLabelOperator
        )
        self.edgeMatch: Callable = generic_edge_match(self.edgeAttribute, 1, eq)

    @staticmethod
    def calculate_cluster_percentages(input_dict: dict) -> dict:
        """
        Calculates the percentage of times each cluster ID appears in the input dictionary
        and returns a new dictionary with cluster IDs as keys and their percentages as
        values, sorted by these percentages in descending order.

        Parameters:
        - input_dict (dict): A dictionary where keys are indices or identifiers and values
        are cluster IDs (integers).

        Returns:
        - dict: A dictionary mapping each cluster ID to the percentage of its appearance,
        sorted by percentage in descending order.
        """
        total_count = len(input_dict)
        cluster_count = {}

        # Count occurrences of each cluster ID
        for cluster_id in input_dict.values():
            cluster_count[cluster_id] = cluster_count.get(cluster_id, 0) + 1

        # Calculate percentage for each cluster ID
        cluster_percentages = {
            cluster_id: round((count / total_count) * 100, 2)
            for cluster_id, count in cluster_count.items()
        }

        # Sort the percentages dictionary by percentage in descending order
        sorted_percentages = dict(
            sorted(cluster_percentages.items(), key=lambda item: item[1], reverse=True)
        )

        return sorted_percentages

    @staticmethod
    def auto_cluster(
        graphs: List[nx.Graph], nodeMatch=None, edgeMatch=None
    ) -> List[Set[int]]:
        """
        Clusters the graphs based on isomorphism, using the predefined node and edge match
        functions.

        Parameters:
        - graphs (List[nx.Graph]): A list of NetworkX graph objects to be clustered.

        Returns:
        - List[Set[int]]: A list of sets, where each set contains the indices of graphs in
        the same cluster.
        """
        visited: Set[int] = set()
        clusters = []
        graph_to_cluster = {}

        for i, graph_i in enumerate(graphs):
            if i in visited:
                continue

            cluster: Set[int] = {i}
            visited.add(i)
            graph_to_cluster[i] = len(clusters)

            for j, graph_j in enumerate(graphs[i + 1 :], start=i + 1):
                if j not in visited and nx.is_isomorphic(
                    graph_i,
                    graph_j,
                    node_match=nodeMatch,
                    edge_match=edgeMatch,
                ):
                    cluster.add(j)
                    visited.add(j)
                    graph_to_cluster[j] = len(clusters)

            clusters.append(cluster)

        return clusters, graph_to_cluster

    @staticmethod
    def library_check(
        graphs: List[nx.Graph],
        templates: List[Dict[str, Any]],
        nodeMatch: Optional[callable] = None,
        edgeMatch: Optional[callable] = None,
    ) -> List[Any]:
        """
        Matches each graph in the provided list against a set of template graphs to
        identify corresponding clusters, and avoids further processing for a graph once a
        match is found.

        Parameters:
        - graphs (List[nx.Graph]): A list of NetworkX graph objects to be checked.
        - templates (List[Dict[str, Any]]): A list of dictionaries, each containing a
        template graph ('RC') and associated cluster ID ('ClusterId').
        - nodeMatch (callable, optional): A function to match nodes between graphs.
        Defaults to None.
        - edgeMatch (callable, optional): A function to match edges between
        - graphs. Defaults to Other.

        Returns:
        - List[Any]: A list containing the 'ClusterId' for each graph in 'graphs' if a
        match is found; otherwise, 'Undefined' for graphs without a match.
        """
        results = []

        for graph in graphs:
            found_match = False

            for template in templates:
                rc = template["RC"]
                if nx.is_isomorphic(
                    graph, rc, node_match=nodeMatch, edge_match=edgeMatch
                ):
                    results.append(template["Cluster_id"])
                    found_match = True
                    break

            if not found_match:
                results.append("Undefined")

        return results

    @staticmethod
    def get_templates(
        graphs: List[List[nx.Graph]],
        single_graphs: List[nx.Graph],
        graph_to_cluster_dict: Dict[int, int],
        sorted_templates: Dict[int, float],
        max_index_template: int = 0,
    ) -> List[Dict[str, Any]]:
        """
        Generates a list of templates from graphs based on cluster mappings, offsetting
        cluster indices by a maximum template index and incorporating percentage values.

        Parameters:
        - graphs (List[List[nx.Graph]]): A list of graph object pairs from which
        templates are derived.
        - single_graphs (List[nx.Graph]): A list of individual graph objects for
        additional data association.
        - graph_to_cluster_dict (Dict[int, int]): A dictionary mapping graph indices to
        cluster IDs.
        - sorted_templates (Dict[int, float]): A dictionary mapping cluster IDs to
        percentage values.
        - max_index_template (int): The maximum index used to offset cluster IDs.

        Returns:
        - List[Dict[str, Any]]: A list of dictionaries each representing a template with a
        'Cluster_id', the associated graph ('RC'), and a 'percentage'.
        """
        # Adjust cluster IDs based on the maximum template index
        temp_graph_to_cluster = {
            key: item + max_index_template
            for key, item in graph_to_cluster_dict.items()
        }
        updated_graphs = []
        for index, graph in enumerate(graphs):
            updated_graph = (graph[0], graph[1], single_graphs[index])
            updated_graphs.append(updated_graph)

        # Create a dictionary with unique values and their first corresponding keys
        unique_temp = create_unique_value_dict(temp_graph_to_cluster)

        template = [
            {
                "Cluster_id": key,
                "RC": updated_graphs[value],
                "Parent": [],
                "Percentage": sorted_templates.get(
                    key, 0
                ),  # Default to 0 if key not found
            }
            for key, value in unique_temp.items()
        ]

        return template

    def fit(
        self,
        graphs: List[nx.Graph],
        templates: Optional[List[Dict[str, Any]]] = None,
        update_template=True,
    ) -> Tuple[List[int], List[Dict[str, Any]]]:
        """
        Automatically clusters the graphs and determines their cluster indices,
        potentially using provided templates for clustering, or generating new templates.

        Parameters:
        - graphs (List[nx.Graph, nx.Graph, nx.Graph]): A list of NetworkX graph objects to
        determine cluster indices for.
        - templates (Optional[List[Dict[str, Any]]]): Optional list of templates used for
        clustering.
        - update_template (bool) : Update new template or not

        Returns:
        - tuple:
            - List[int]: The list of cluster indices for each graph, aligned with the
            order of the input list.
            - List[Dict[str, Any]]: Updated or newly created list of templates.
        """
        if templates is None:
            single_graphs = [value[2] for value in graphs]

            # Perform clustering without predefined templates
            _, graph_to_cluster_dict = self.auto_cluster(
                single_graphs, self.nodeMatch, self.edgeMatch
            )
            sorted_templates = RuleCluster.calculate_cluster_percentages(
                graph_to_cluster_dict
            )
            templates = self.get_templates(
                graphs, single_graphs, graph_to_cluster_dict, sorted_templates, 0
            )

            cluster_indices = [
                graph_to_cluster_dict.get(i, None) for i in range(len(graphs))
            ]

        else:
            # Use existing templates to check graph clusters
            cluster_indices = RuleCluster.library_check(graphs, templates)
            undefined_keys = [
                i for i, cid in enumerate(cluster_indices) if cid == "Undefined"
            ]
            if update_template:
                undefined_keys = [
                    i for i, cid in enumerate(cluster_indices) if cid == "Undefined"
                ]

                if undefined_keys:
                    # Handle undefined clusters by re-clustering them
                    undefined_graphs = [graphs[i] for i in undefined_keys]
                    _, new_graph_to_cluster_dict = self.auto_cluster(
                        undefined_graphs, self.nodeMatch, self.edgeMatch
                    )
                    max_index_template = (
                        max(t["Cluster_id"] for t in templates) if templates else 0
                    )

                    new_templates = self.get_templates(
                        undefined_graphs, new_graph_to_cluster_dict, max_index_template
                    )
                    templates += new_templates

                    # Update cluster indices for re-clustered graphs
                    new_cluster_indices = [
                        new_graph_to_cluster_dict[i] + max_index_template
                        for i in range(len(undefined_graphs))
                    ]
                    for i, key in enumerate(undefined_keys):
                        cluster_indices[key] = new_cluster_indices[i]

        return cluster_indices, templates
