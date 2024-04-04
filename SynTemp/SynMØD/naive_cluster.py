import networkx as nx
from typing import List, Set, Dict, Any, Callable
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match
from operator import eq

class NaiveCluster:
    def __init__(self, node_label_names: List[str] = ["element", "aromatic", "hcount", "charge", "typesGH"],
                 node_label_default: List[Any] = ["*", False, 0, 0, ()], edge_attribute: str = "order"):
        """
        Initializes the NaiveClusterer with customization options for node and edge matching functions.
        
        Parameters:
            node_label_names (List[str]): A list of node attribute names to be considered for matching.
            node_label_default (List[Any]): Default values for node attributes, aligned with `node_label_names`.
            edge_attribute (str): The name of the edge attribute to be considered for matching.
        """
    
        self.nodeLabelNames: List[str] = node_label_names
        self.edgeAttribute: str = edge_attribute
        self.nodeLabelDefault: List[Any] = node_label_default
        self.nodeLabelOperator: List[Callable[[Any, Any], bool]] = [eq for _ in node_label_names]
        self.nodeMatch: Callable = generic_node_match(self.nodeLabelNames, self.nodeLabelDefault, self.nodeLabelOperator)
        self.edgeMatch: Callable = generic_edge_match(self.edgeAttribute, 1, eq)
        self.clusters: List[Set[int]] = []
        self.graph_to_cluster: Dict[int, int] = {}

    def cluster_graphs(self, graphs: List[nx.Graph]) -> List[Set[int]]:
        """
        Clusters the graphs based on isomorphism, using the predefined node and edge match functions.
        
        Parameters:
            graphs (List[nx.Graph]): A list of NetworkX graph objects to be clustered.
            
        Returns:
            List[Set[int]]: A list of sets, where each set contains the indices of graphs in the same cluster.
        """
        visited: Set[int] = set()
        self.clusters = []  # Reset clusters for each call
        self.graph_to_cluster = {}  # Reset graph_to_cluster mapping

        for i, graph_i in enumerate(graphs):
            if i in visited:
                continue

            cluster: Set[int] = {i}
            visited.add(i)
            self.graph_to_cluster[i] = len(self.clusters)

            for j, graph_j in enumerate(graphs[i+1:], start=i+1):
                if j not in visited and nx.is_isomorphic(graph_i, graph_j, node_match=self.nodeMatch, edge_match=self.edgeMatch):
                    cluster.add(j)
                    visited.add(j)
                    self.graph_to_cluster[j] = len(self.clusters)
            
            self.clusters.append(cluster)
        
        return self.clusters

    def get_cluster_indices(self, graphs: List[nx.Graph]) -> List[int]:
        """
        Returns a list where each element is the cluster index for the corresponding graph in the input list.
        This method automatically clusters the graphs before determining their indices.
        
        Parameters:
            graphs (List[nx.Graph]): A list of NetworkX graph objects to determine cluster indices for.
            
        Returns:
            List[int]: The list of cluster indices for each graph, aligned with the order of the input list.
        """
        self.cluster_graphs(graphs)  # Ensure clusters are updated based on the current graphs list
        return [self.graph_to_cluster[i] for i in range(len(graphs))]

    def process_rules_clustering(self, reaction_dicts: List[Dict[str, Any]], rule_column: str = 'rules') -> List[Dict[str, Any]]:
        """
        Processes clustering based on rules extracted from reaction dictionaries, specifically clustering ITS graphs.
        
        Parameters:
            reaction_dicts (List[Dict[str, Any]]): A list of dictionaries, each representing a reaction.
            rule_column (str): The key in the dictionaries where the ITS graph is stored.
            
        Returns:
            List[Dict[str, Any]]: The updated list of reaction dictionaries, each augmented with a 'cluster' key indicating its cluster index.
        """
        # Extract ITS graphs from the reaction dictionaries
        rules_graphs = [reaction[rule_column][2] for reaction in reaction_dicts]

        # Cluster the ITS graphs and get cluster indices for each ITS graph
        self.cluster_graphs(rules_graphs)
        cluster_indices = self.get_cluster_indices(rules_graphs)

        # Update the reaction dictionaries with cluster information
        for i, reaction_dict in enumerate(reaction_dicts):
            reaction_dict['naive_cluster'] = cluster_indices[i]

        return reaction_dicts
