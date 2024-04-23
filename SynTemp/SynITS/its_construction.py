import networkx as nx
from typing import Tuple, Dict, Any
from copy import deepcopy


class ITSConstruction:
    @staticmethod
    def ITSGraph(
        G: nx.Graph,
        H: nx.Graph,
        ignore_aromaticity: bool = False,
        attributes_defaults: Dict[str, Any] = None,
    ) -> nx.Graph:
        """
        Creates a Combined Graph Representation (CGR) from two input graphs G and H.

        This function merges the nodes of G and H, preserving their attributes. Edges are added
        based on their presence in G and/or H, with special labeling for edges unique to one graph.

        Parameters:
        - G (nx.Graph): The first input graph.
        - H (nx.Graph): The second input graph.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in the graphs. Defaults to False.
        - attributes_defaults (Dict[str, Any]): A dictionary of default attributes to use for nodes that are not present in either G or H.

        Returns:
        - nx.Graph: The Combined Graph Representation as a new graph instance.
        """
        # Create a null graph from a copy of G to preserve attributes
        if len(G.nodes()) > len(H.nodes()):
            ITS = deepcopy(H)
        else:
            ITS = deepcopy(G)
        ITS.remove_edges_from(list(ITS.edges()))

        # Initialize a dictionary to hold node types
        typesDict = dict()

        # Add typeG and typeH attributes, or default attributes for "*" unknown elements
        for v in list(ITS.nodes()):
            # Check if v is in both G and H
            if v not in G.nodes() or v not in H.nodes():
                continue
            else:
                typesG = ITSConstruction.get_node_attributes_with_defaults(
                    G, v, attributes_defaults
                )  # node attribute in reactant graph
                typesH = ITSConstruction.get_node_attributes_with_defaults(
                    H, v, attributes_defaults
                )  # node attribute in product graph
                typesDict[v] = (typesG, typesH)

        nx.set_node_attributes(ITS, typesDict, "typesGH")

        # Add edges from G and H
        ITS = ITSConstruction.add_edges_to_ITS(ITS, G, H, ignore_aromaticity)

        return ITS

    # @staticmethod
    # def ITSGraph(G: nx.Graph, H: nx.Graph, ignore_aromaticity: bool = False, attributes_defaults: Dict[str, Any] = None) -> nx.Graph:
    #     """
    #     Creates a Combined Graph Representation (CGR) from two input graphs G and H by merging their nodes and edges.
    #     Nodes preserve their attributes, and edges are labeled based on their presence in G and/or H.

    #     Parameters:
    #     - G (nx.Graph): The first input graph.
    #     - H (nx.Graph): The second input graph.
    #     - ignore_aromaticity (bool, optional): If True, aromaticity in the graphs is ignored. Defaults to False.
    #     - attributes_defaults (Dict[str, Any], optional): Default attributes for nodes not present in either G or H.

    #     Returns:
    #     - nx.Graph: The Combined Graph Representation as a new graph instance.
    #     """
    #     # Use the smaller graph as the base for ITS
    #     ITS = deepcopy(H if len(G.nodes()) > len(H.nodes()) else G)
    #     ITS.clear_edges()  # Remove all edges while retaining nodes and their attributes

    #     # Iterate over nodes in ITS and get attributes from G and H, or set defaults
    #     for node in ITS:
    #         attributes_in_G = ITSConstruction.get_node_attributes_with_defaults(G, node, attributes_defaults) if node in G else None
    #         attributes_in_H = ITSConstruction.get_node_attributes_with_defaults(H, node, attributes_defaults) if node in H else None

    #         # Combine attributes from G and H or use defaults
    #         ITS.nodes[node]['typesGH'] = (attributes_in_G, attributes_in_H)

    #     # Add edges from G and H with special labeling for unique edges
    #     ITS = ITSConstruction.add_edges_to_ITS(ITS, G, H, ignore_aromaticity)

    #     return ITS

    @staticmethod
    def get_node_attribute(graph: nx.Graph, node: int, attribute: str, default):
        """
        Retrieves a specific attribute for a node in a graph, returning a default value if the attribute is missing.

        Parameters:
        - graph (nx.Graph): The graph from which to retrieve the node attribute.
        - node (int): The node identifier.
        - attribute (str): The attribute to retrieve.
        - default: The default value to return if the attribute is missing.

        Returns:
        - The value of the node attribute, or the default value if the attribute is missing.
        """
        try:
            return graph.nodes[node][attribute]
        except KeyError:
            return default

    @staticmethod
    def get_node_attributes_with_defaults(
        graph: nx.Graph, node: int, attributes_defaults: Dict[str, Any] = None
    ) -> Tuple:
        """
        Retrieves node attributes from a graph, assigning default values if they are missing. Allows
        for an optional dictionary of attribute-default value pairs to specify custom attributes and defaults.

        Parameters:
        - graph (nx.Graph): The graph from which to retrieve node attributes.
        - node (int): The node identifier.
        - attributes_defaults (Dict[str, Any], optional): A dictionary specifying attributes and their default values.

        Returns:
        - Tuple: A tuple containing the node attributes in the order specified by attributes_defaults.
        """
        if attributes_defaults is None:
            attributes_defaults = {
                "element": "*",
                "aromatic": False,
                "hcount": 0,
                "charge": 0,
            }

        return tuple(
            ITSConstruction.get_node_attribute(graph, node, attr, default)
            for attr, default in attributes_defaults.items()
        )

    @staticmethod
    def add_edges_to_ITS(
        ITS: nx.Graph, G: nx.Graph, H: nx.Graph, ignore_aromaticity: bool = False
    ) -> nx.Graph:
        """
        Adds edges to the Combined Graph Representation (ITS) based on the edges of G and H,
        and returns a new graph without modifying the original ITS.

        Parameters:
        - ITS (nx.Graph): The initial combined graph representation.
        - G (nx.Graph): The first input graph.
        - H (nx.Graph): The second input graph.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in the graphs. Defaults to False.

        Returns:
        - nx.Graph: The updated graph with added edges.
        """
        new_ITS = deepcopy(ITS)

        # Add edges from G and H
        for graph_from, graph_to, reverse in [(G, H, False), (H, G, True)]:
            for u, v in graph_from.edges():
                if not new_ITS.has_edge(u, v):
                    if graph_to.has_edge(u, v) or graph_to.has_edge(v, u):
                        edge_label = (
                            (graph_from[u][v]["order"], graph_to[u][v]["order"])
                            if graph_to.has_edge(u, v)
                            else (
                                (graph_from[v][u]["order"], graph_to[v][u]["order"])
                                if reverse
                                else (
                                    graph_from[u][v]["order"],
                                    graph_to[v][u]["order"],
                                )
                            )
                        )
                        new_ITS.add_edge(u, v, order=edge_label)
                    else:
                        edge_label = (
                            (graph_from[u][v]["order"], 0)
                            if not reverse
                            else (0, graph_from[u][v]["order"])
                        )
                        new_ITS.add_edge(u, v, order=edge_label)
        nodes_to_remove = [node for node in new_ITS.nodes() if not new_ITS.nodes[node]]
        new_ITS.remove_nodes_from(nodes_to_remove)
        new_ITS = ITSConstruction.add_standard_order_attribute(
            new_ITS, ignore_aromaticity
        )
        return new_ITS

    @staticmethod
    def add_standard_order_attribute(
        graph: nx.Graph, ignore_aromaticity: bool = False
    ) -> nx.Graph:
        """
        Adds a 'standard_order' attribute to each edge in the given NetworkX graph based on the 'order' attribute.
        The 'standard_order' is computed by subtracting the second element of the 'order' tuple from the first.
        If an element of the tuple is not an integer (e.g., '*'), it is treated as 0 for the computation.

        :param graph: A NetworkX graph whose edges have 'order' attributes as tuples.
        :return: A new NetworkX graph with 'standard_order' attributes added to each edge.
        """
        new_graph = graph.copy()

        for u, v, data in new_graph.edges(data=True):
            if "order" in data and isinstance(data["order"], tuple):
                # Extract order values, replacing non-ints with 0
                first_order = data["order"][0]
                second_order = data["order"][1]
                # Compute standard order
                standard_order = first_order - second_order
                if ignore_aromaticity:
                    if abs(standard_order) < 1:  # to ignore aromaticity
                        standard_order = 0
                # Update the edge data with a new attribute 'standard_order'
                new_graph[u][v]["standard_order"] = standard_order
            else:
                # If 'order' attribute is missing or not a tuple, set 'standard_order' to 0
                new_graph[u][v]["standard_order"] = 0

        return new_graph
