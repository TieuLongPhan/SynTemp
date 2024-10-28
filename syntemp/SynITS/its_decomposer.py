import networkx as nx


def its_decompose(its_graph: nx.Graph, nodes_share="typesGH", edges_share="order"):
    """
    Decompose an ITS graph into two separate graphs G and H based on shared
    node and edge attributes.

    Parameters:
    - its_graph (nx.Graph): The integrated transition state (ITS) graph.
    - nodes_share (str): Node attribute key that stores tuples with node attributes
    or G and H.
    - edges_share (str): Edge attribute key that stores tuples with edge attributes
    for G and H.

    Returns:
    - Tuple[nx.Graph, nx.Graph]: A tuple containing the two graphs G and H.
    """
    G = nx.Graph()
    H = nx.Graph()

    # Decompose nodes
    for node, data in its_graph.nodes(data=True):
        if nodes_share in data:
            node_attr_g, node_attr_h = data[nodes_share]
            # Unpack node attributes for G
            G.add_node(
                node,
                element=node_attr_g[0],
                aromatic=node_attr_g[1],
                hcount=node_attr_g[2],
                charge=node_attr_g[3],
                neighbors=node_attr_g[4],
            )
            # Unpack node attributes for H
            H.add_node(
                node,
                element=node_attr_h[0],
                aromatic=node_attr_h[1],
                hcount=node_attr_h[2],
                charge=node_attr_h[3],
                neighbors=node_attr_h[4],
            )

    # Decompose edges
    for u, v, data in its_graph.edges(data=True):
        if edges_share in data:
            order_g, order_h = data[edges_share]
            if order_g > 0:  # Assuming 0 means no edge in G
                G.add_edge(u, v, order=order_g)
            if order_h > 0:  # Assuming 0 means no edge in H
                H.add_edge(u, v, order=order_h)

    return G, H


def compare_graphs(
    graph1: nx.Graph,
    graph2: nx.Graph,
    node_attrs: list = ["element", "aromatic", "hcount", "charge", "neighbors"],
    edge_attrs: list = ["order"],
) -> bool:
    """
    Compare two graphs based on specified node and edge attributes.

    Parameters:
    - graph1 (nx.Graph): The first graph to compare.
    - graph2 (nx.Graph): The second graph to compare.
    - node_attrs (list): A list of node attribute names to include in the comparison.
    - edge_attrs (list): A list of edge attribute names to include in the comparison.

    Returns:
    - bool: True if both graphs are identical with respect to the specified attributes,
    otherwise False.
    """
    # Compare node sets
    if set(graph1.nodes()) != set(graph2.nodes()):
        return False

    # Compare nodes based on attributes
    for node in graph1.nodes():
        if node not in graph2:
            return False
        node_data1 = {attr: graph1.nodes[node].get(attr, None) for attr in node_attrs}
        node_data2 = {attr: graph2.nodes[node].get(attr, None) for attr in node_attrs}
        if node_data1 != node_data2:
            return False

    # Compare edge sets with sorted tuples
    if set(tuple(sorted(edge)) for edge in graph1.edges()) != set(
        tuple(sorted(edge)) for edge in graph2.edges()
    ):
        return False

    # Compare edges based on attributes
    for edge in graph1.edges():
        # Sort the edge for consistent comparison
        sorted_edge = tuple(sorted(edge))
        if sorted_edge not in graph2.edges():
            return False
        edge_data1 = {attr: graph1.edges[edge].get(attr, None) for attr in edge_attrs}
        edge_data2 = {
            attr: graph2.edges[sorted_edge].get(attr, None) for attr in edge_attrs
        }
        if edge_data1 != edge_data2:
            return False

    return True
