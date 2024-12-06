import networkx as nx
from operator import eq
from collections import OrderedDict
from typing import List, Set, Dict, Any, Callable, Tuple
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match


class RCCluster:
    def __init__(
        self,
        node_label_names: List[str] = [
            "element",
            "charge",
        ],
        node_label_default: List[Any] = ["*", 0],
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
    def auto_cluster(
        graphs: List[nx.Graph], attribute: List[Any], nodeMatch=None, edgeMatch=None
    ) -> List[Set[int]]:
        """
        Clusters the graphs based on isomorphism and cycles similarity,
        using specified node and edge matching functions.

        Parameters:
        - graphs (List[nx.Graph]): List of NetworkX graph objects to be clustered.
        - attribute (List[Any]): Attributes associated with each graph for
        preliminary comparison.
        - nodeMatch (Optional[Callable]): Custom node matching function;
        if None, uses the class's default.
        - edgeMatch (Optional[Callable]): Custom edge matching function;
        if None, uses the class's default.

        Returns:
        - Tuple[List[Set[int]], Dict[int, int]]: A tuple containing a list of sets (
        clusters), where each set contains indices of graphs in the same cluster,
        and a dictionary mapping each graph index to its cluster index.
        """
        if attribute is None:
            att_sorted = [1] * len(graphs)
        else:
            if isinstance(attribute[0], str):
                att_sorted = attribute
            elif isinstance(attribute, List):
                att_sorted = [sorted(value) for value in attribute]
            elif isinstance(attribute, OrderedDict):
                att_sorted = [OrderedDict(sorted(value.items())) for value in attribute]

        visited: Set[int] = set()
        clusters = []
        graph_to_cluster = {}

        for i, graph_i in enumerate(graphs):
            if i in visited:
                continue

            cluster: Set[int] = {i}
            visited.add(i)
            graph_to_cluster[i] = len(clusters)

            # fmt: off
            for j, graph_j in enumerate(graphs[i + 1:], start=i + 1):
                # fmt: on
                if att_sorted[i] == att_sorted[j]:
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

    def fit(
        self, data: List[Dict], RC_col: str = "RC", attribute: str = "cycle"
    ) -> Tuple[List[int], List[Dict[str, Any]]]:
        """
        Automatically clusters the graphs and determines their cluster indices,
        potentially using provided templates for clustering, or generating new templates.

        Parameters:
        - graphs (List[nx.Graph, nx.Graph, nx.Graph]): A list of NetworkX graph objects to
        determine cluster indices for.
        - attribute (List[Any]): Attributes associated with each graph for
        preliminary comparison.


        Returns:
        - tuple:
            - List[int]: The list of cluster indices for each graph, aligned with the
            order of the input list.
            - List[Dict[str, Any]]: Updated or newly created list of templates.
        """
        new_data = data
        data_rc = [value[RC_col] for value in new_data]
        if attribute is None:
            att = None
        else:
            att = [value[attribute] for value in new_data]

        # Perform clustering without predefined templates
        _, graph_to_cluster_dict = self.auto_cluster(
            data_rc, att, self.nodeMatch, self.edgeMatch
        )

        cluster_indices = [graph_to_cluster_dict.get(i, None) for i in range(len(data))]

        for key, value in enumerate(new_data):
            value["class"] = cluster_indices[key]

        return new_data

    def fit_graphs(
        self, graphs_data: List[nx.Graph], attribute: List[str]
    ) -> Tuple[List[int], List[Dict[str, Any]]]:
        """
        Automatically clusters the input graphs based on provided attributes and
        determines the cluster indices for each graph. The method may use existing
        templates for clustering or generate new ones.

        Parameters:
        - graphs_data (List[nx.Graph]): A list of NetworkX graph objects to cluster.
        - attribute (List[str]): A list of attributes associated with each graph,
        used for clustering comparisons.

        Returns:
        - Tuple[List[int], List[Dict[str, Any]]]:
            - A list of cluster indices (List[int]) corresponding to each graph in
            the input list, indicating which cluster each graph belongs to.
            - A list of updated or newly generated templates (List[Dict[str, Any]]),
            representing the cluster information or other relevant data for each
            group of graphs.
        """

        # Perform clustering without predefined templates
        _, graph_to_cluster_dict = self.auto_cluster(
            graphs_data, attribute, self.nodeMatch, self.edgeMatch
        )

        # Generate the cluster indices based on the clustering result
        cluster_indices = [
            graph_to_cluster_dict.get(i, None) for i in range(len(graphs_data))
        ]

        return cluster_indices
