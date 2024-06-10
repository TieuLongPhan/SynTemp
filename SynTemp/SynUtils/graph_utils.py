import networkx as nx
import copy
from typing import List, Dict, Any


def is_acyclic_graph(G: nx.Graph) -> bool:
    """
    Determines if the given graph is acyclic.

    Parameters:
    - G (nx.Graph): The graph to be checked.

    Returns:
    - bool: True if the graph is acyclic, False otherwise.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    return nx.is_tree(G)


def is_single_cyclic_graph(G: nx.Graph) -> bool:
    """
    Determines if the given graph is a single cyclic graph,
    which means the graph has exactly one cycle
    and all nodes in the graph are part of that cycle.

    Parameters:
    - G (nx.Graph): The graph to be checked.

    Returns:
    - bool: True if the graph is single cyclic, False otherwise.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    if not nx.is_connected(G):
        return False

    cycles = nx.cycle_basis(G)

    if cycles:
        nodes_in_cycles = set(node for cycle in cycles for node in cycle)
        if (
            nodes_in_cycles == set(G.nodes())
            and G.number_of_edges() == G.number_of_nodes()
        ):
            return True

    return False


def is_complex_cyclic_graph(G: nx.Graph) -> bool:
    """
    Determines if the given graph is a complex cyclic graph,
    which means all nodes are part of cycles,
    there are multiple cycles, and there are no acyclic parts.

    Parameters:
    - G (nx.Graph): The graph to be checked.

    Returns:
    - bool: True if the graph is complex cyclic, False otherwise.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    # Check if the graph is connected and has at least one cycle
    if not nx.is_connected(G) or not any(nx.cycle_basis(G)):
        return False

    # Get a list of cycles that form a cycle basis for G
    cycles = nx.cycle_basis(G)

    # If there's only one cycle in the basis, it might not be a complex cyclic graph
    if len(cycles) <= 1:
        return False

    # Decompose cycles into a list of nodes, allowing for node overlap between cycles
    nodes_in_cycles = set(node for cycle in cycles for node in cycle)

    # Check if all nodes in G are covered by the nodes in cycles
    return nodes_in_cycles == set(G.nodes())


def check_graph_type(G: nx.Graph) -> str:
    """
    Determines if the given graph is acyclic, single cyclic, or complex cyclic.

    Parameters:
    - G (nx.Graph): The graph to be checked.

    Returns:
    - str: A string indicating if the graph is "Acyclic",
            "Single Cyclic", or "Complex Cyclic".

    Raises:
    - TypeError: If the input G is not a networkx Graph.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    if is_acyclic_graph(G):
        return "Acyclic"
    elif is_single_cyclic_graph(G):
        return "Single Cyclic"
    elif is_complex_cyclic_graph(G):
        return "Complex Cyclic"
    else:
        return "None"


def get_cycle_member_rings(G: nx.Graph) -> List[int]:
    """
    Identifies all cycles in the given graph using cycle bases to ensure no overlap
    and returns a list of the sizes of these cycles (member rings),
    sorted in ascending order.

    Parameters:
    - G (nx.Graph): The NetworkX graph to be analyzed.

    Returns:
    - List[int]: A sorted list of cycle sizes (member rings) found in the graph.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    # Find cycle basis for the graph which gives non-overlapping cycles
    cycles = nx.cycle_basis(G)

    # Determine the size of each cycle (member ring)
    member_rings = [len(cycle) for cycle in cycles]

    # Sort the sizes in ascending order
    member_rings.sort()

    return member_rings


def load_gml_as_text(gml_file_path):
    """
    Load the contents of a GML file as a text string.

    Parameters:
    - gml_file_path (str): The file path to the GML file.

    Returns:
    - str: The text content of the GML file.
    """
    try:
        with open(gml_file_path, "r") as file:
            return file.read()
    except FileNotFoundError:
        print(f"File not found: {gml_file_path}")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def add_child_ids(df: List[List[Dict[str, Any]]]) -> List[List[Dict[str, Any]]]:
    """
    Processes hierarchical data to assign child IDs based on parent-cluster relationships.
    
    Each node in the hierarchy should have a Cluster_id and optionally a Parent. 
    This function will add a Child field to each node, which is a list of Cluster_ids
    from the child nodes in the subsequent layer.

    Args:
    data (List[List[Dict[str, Any]]]): A list of layers, where each layer is a list of dictionaries 
                                       containing at least the keys 'Cluster_id' and 'Parent'.

    Returns:
    List[List[Dict[str, Any]]]: The modified data with each node dictionary containing only 'Cluster_id',
                                'Parent', and 'Child' keys.
    """
    data = copy.deepcopy(df)
    # Initialize a dictionary to hold each node with its unique key as the combination of layer index and Cluster_id
    node_dict = {}

    # Process each layer to assign children based on the parent references
    for layer_index, layer in enumerate(data):
        for node in layer:
            # Generate a unique identifier combining layer index with Cluster_id for internal mapping
            unique_id = f"{layer_index}-{node['Cluster_id']}"
            
            # Initialize the list to hold the cluster_ids of child nodes
            node['Child'] = []
            
            # Store the node in the dictionary with its unique identifier
            node_dict[unique_id] = node

    # Link children to their parents
    for layer_index, layer in enumerate(data):
        if layer_index == 0:
            continue  # Skip the first layer as it has no parents

        for node in layer:
            # Normalize the Parent field to always be a list
            if isinstance(node['Parent'], int):
                node['Parent'] = [node['Parent']]
            
            # Handle the parent list
            for parent_cluster_id in node['Parent']:
                # Compute the parent's unique identifier from the previous layer
                parent_unique_id = f"{layer_index-1}-{parent_cluster_id}"

                # Link the child node's Cluster_id to the parent node's Child list
                if parent_unique_id in node_dict:
                    parent_node = node_dict[parent_unique_id]
                    parent_node['Child'].append(node['Cluster_id'])
                                                     

    # Prepare the final data structure with only the required keys
    for layer in data:
        for node in layer:
            # Retain only Cluster_id, Parent, and Child in each node
            keys_to_keep = {'Cluster_id', 'Parent', 'Child'}
            for key in list(node.keys()):
                if key not in keys_to_keep:
                    del node[key]

    return data