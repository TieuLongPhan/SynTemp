import networkx as nx
from typing import List, Set, Dict, Any
from joblib import Parallel, delayed
from itertools import chain


class RuleExtraction:

    @staticmethod
    def find_unequal_order_edges(G: nx.Graph) -> List[int]:
        """
        Identifies reaction center nodes in a graph based on unequal order edges.

        Parameters:
        - G (nx.Graph): The graph in which to identify reaction center nodes.

        Returns:
        - List[int]: Nodes that are part of edges with unequal orders, indicating reaction
        centers.
        """
        reaction_center_nodes: Set[int] = set()
        for u, v, data in G.edges(data=True):
            order = data.get("order", (1, 1))
            if isinstance(order, tuple) and order[0] != order[1]:
                if data.get("standard_order") != 0:
                    reaction_center_nodes.add(u)
                    reaction_center_nodes.add(v)
        return list(reaction_center_nodes)

    @staticmethod
    def find_nearest_neighbors(
        G: nx.Graph, center_nodes: List[int], n_knn: int = 1
    ) -> Set[int]:
        """
        Finds up to n_knn levels of nearest neighbors for specified center nodes within a
        graph.

        Parameters:
        - G (nx.Graph): The graph to search within.
        - center_nodes (List[int]): Nodes for which to find nearest neighbors.
        - n_knn (int): Levels of nearest neighbors to include.

        Returns:
        - Set[int]: Nodes including the center nodes and their up to n_knn-nearest
        neighbors.
        """
        extended_nodes = set(center_nodes)
        for _ in range(n_knn):
            neighbors = set(
                chain.from_iterable(G.neighbors(node) for node in extended_nodes)
            )
            extended_nodes.update(neighbors)
        return extended_nodes

    @staticmethod
    def extract_subgraph(G: nx.Graph, node_indices: List[int]) -> nx.Graph:
        """
        Extracts a subgraph based on specified node indices from a given graph.

        Parameters:
        - G (nx.Graph): The original graph.
        - node_indices (List[int]): Node indices to include in the subgraph.

        Returns:
        - nx.Graph: The extracted subgraph.
        """

        return G.subgraph(node_indices).copy()

    @staticmethod
    def extract_reaction_rules(
        reactants_graph, products_graph, its_graph, extend: bool = False, n_knn: int = 1
    ) -> Dict[str, Any]:
        """
        Extracts transformation rules for a chemical reaction represented by graphs of
        reactants, products, and intermediate transition states (ITS).

        This method identifies the reaction center nodes, optionally extends them based on
        nearest neighbors, and constructs a tuple of subgraphs representing the rules.

        Parameters:
        - reactants_graph (Graph): Graph representing the reactants.
        - products_graph (Graph): Graph representing the products.
        - its_graph (Graph): Graph representing the intermediate transition states.
        - extend (bool, optional): If True, extends the reaction center nodes by including
        nearest neighbors. Defaults to True.
        - n_knn (int, optional): Specifies the number of nearest neighbors to consider
        when extending the reaction center nodes. Defaults to 1.

        Returns:
        - Tuple: A tuple contains the subgraphs for the reactants, products, and ITS,
        representing the extracted rules.
        """
        reaction_center_nodes = RuleExtraction.find_unequal_order_edges(its_graph)
        if extend:
            reaction_center_nodes = RuleExtraction.find_nearest_neighbors(
                its_graph, reaction_center_nodes, n_knn
            )

        rules_network = tuple(
            RuleExtraction.extract_subgraph(graph, reaction_center_nodes)
            for graph in (reactants_graph, products_graph, its_graph)
        )
        rules_graph = rules_network[2]
        if n_knn == 0:
            rules_graph = RuleExtraction.remove_normal_egdes(
                rules_graph, "standard_order"
            )
        rules_network = (rules_network[0], rules_network[1], rules_graph)
        return rules_network

    @staticmethod
    def rules_extraction(
        reaction_dict: Dict[str, Any],
        mapper_type: str = "ITSGraph",
        extend: bool = True,
        n_knn: int = 1,
    ) -> Dict[str, Any]:
        """
        Extracts rules for a single reaction using a specified mapper, optionally
        extending reaction centers.

        Parameters:
        - reaction_dict (Dict[str, Any]): Single reaction representation.
        - mapper_type (str): Mapper type for rule extraction.
        - extend (bool): Flag to extend reaction center nodes.
        - n_knn (int): Nearest neighbor levels to include if extending.

        Returns:
        - Dict[str, Any]: Reaction dictionary updated with extracted rules.
        'GraphRules': Tuple of reactants graph (L graph), products graph
        (R graph), ITS graph (K graph)
        """
        reactants_graph, products_graph, its_graph = reaction_dict[mapper_type]
        reaction_center_nodes = RuleExtraction.find_unequal_order_edges(its_graph)

        if extend:
            reaction_center_nodes = RuleExtraction.find_nearest_neighbors(
                its_graph, reaction_center_nodes, n_knn
            )

        rules_network = tuple(
            RuleExtraction.extract_subgraph(graph, list(reaction_center_nodes))
            for graph in (reactants_graph, products_graph, its_graph)
        )

        reaction_dict["GraphRules"] = rules_network
        return reaction_dict

    @classmethod
    def process_rules_extraction(
        cls,
        reaction_dicts: List[Dict[str, Any]],
        mapper_type: str = "ITSGraph",
        n_jobs: int = 4,
        verbose: int = 0,
        extend: bool = True,
        n_knn: int = 1,
    ) -> List[Dict[str, Any]]:
        """
        Parallelizes rules extraction for multiple reactions, with options for extension
        and parallelism control.

        Parameters:
        - reaction_dicts (List[Dict[str, Any]]): Reactions to process.
        - mapper_type (str): Mapper type for rule extraction.
        - n_jobs (int): Parallel jobs count (-1 for all CPUs).
        - verbose (int): Verbosity level.
        - extend (bool): Flag to extend reaction center nodes.
        - n_knn (int): Nearest neighbor levels to include if extending.

        Returns:
        - List[Dict[str, Any]]: Updated reactions with extracted rules.
        """
        results = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(cls.rules_extraction)(reaction_dict, mapper_type, extend, n_knn)
            for reaction_dict in reaction_dicts
        )
        return results

    @staticmethod
    def remove_normal_egdes(graph: nx.Graph, property_key: str) -> nx.Graph:
        """
        Create a copy of the input graph and remove edges based on a specified property.

        This function iterates over the edges of the graph and removes any edge for which
        the specified property's values are the same at both the start and end nodes.

        Parameters:
        graph (nx.Graph): The input graph from which edges will be filtered.
        property_key (str): The key for the edge attribute that should be examined.

        Returns:
        nx.Graph: A new graph with the specified edges removed. The original graph is not
        modified.
        """
        # Create a copy of the graph to avoid inplace modification
        filtered_graph = graph.copy()

        # Iterate over the edges and check the specified property
        for edge in list(filtered_graph.edges(data=True)):
            # Unpack edge data
            start_node, end_node, attributes = edge

            if property_key in attributes and attributes[property_key] == 0:
                # Remove the edge if the property values are equal
                filtered_graph.remove_edge(start_node, end_node)

        return filtered_graph
