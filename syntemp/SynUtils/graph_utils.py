import networkx as nx
import copy
from typing import List, Dict, Any, Tuple


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
        return "Combinatorial Cyclic"
    else:
        return "Complex Cyclic"


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


def get_descriptors(data: List[Dict], reaction_centers: str = "RC") -> List[Dict]:
    """
    Enhance data with topology type and reaction type descriptors.

    Parameters:
    - data (List[Dict]): List of dictionaries containing reaction data.
    - reaction_centers (str): Key for accessing the reaction centers in the data
    dictionaries.

    Returns:
    - List[Dict]: Enhanced list of dictionaries with added descriptors.
    """
    for value in data:
        graph = value[reaction_centers][2]

        # Determine topology type and cycle member rings
        value["Topo Type"] = check_graph_type(graph)
        value["Rings"] = get_cycle_member_rings(graph)

        # Determine reaction type based on topology type
        if value["Topo Type"] in ["Single Cyclic", "Acyclic"]:
            value["Reaction Type"] = "Elementary"
        else:
            value["Reaction Type"] = "Complicated"

        # Set "Rings" based on topology type
        if value["Topo Type"] == "Acyclic":
            value["Rings"] = [0]
        elif value["Topo Type"] == "Complex Cyclic":
            value["Rings"] = [0] + value["Rings"]

        value["Reaction Step"] = len(value["Rings"])

    return data


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

    Parameters:
    - data (List[List[Dict[str, Any]]]): A list of layers, where each layer is a list of
    dictionaries containing at least the keys 'Cluster_id' and 'Parent'.

    Returns:
    - List[List[Dict[str, Any]]]: The modified data with each node dictionary containing
    'Cluster_id', 'Percentage', 'Parent', and 'Child' keys.
    """
    data = copy.deepcopy(df)
    node_dict = {}

    for layer_index, layer in enumerate(data):
        for node in layer:
            unique_id = f"{layer_index}-{node['Cluster_id']}"
            node["Child"] = []  # Initialize the Child list
            node_dict[unique_id] = node

    # Link children to their parents
    for layer_index, layer in enumerate(data):
        if layer_index == 0:
            continue  # Skip the first layer as it has no parents

        for node in layer:
            parents = node.get("Parent", [])
            parents = (
                [parents] if isinstance(parents, int) else parents
            )  # Ensure parents is always a list

            for parent_cluster_id in parents:
                parent_unique_id = f"{layer_index - 1}-{parent_cluster_id}"
                if parent_unique_id in node_dict:
                    node_dict[parent_unique_id]["Child"].append(node["Cluster_id"])

    for layer in data:
        for node in layer:
            keys_to_keep = {"Cluster_id", "Parent", "Child", "Percentage"}
            for key in list(node.keys()):
                if key not in keys_to_keep:
                    del node[key]

    return data


def check_explicit_hydrogen(graph: nx.Graph) -> tuple:
    """
    Counts the explicit hydrogen nodes in the given graph and collects their IDs.

    Args:
    graph (nx.Graph): The graph to inspect.

    Returns:
    tuple: A tuple containing the number of hydrogen nodes and a list of their node IDs.
    """
    hydrogen_nodes = [
        node_id
        for node_id, attr in graph.nodes(data=True)
        if attr.get("element") == "H"
    ]
    return len(hydrogen_nodes), hydrogen_nodes


def check_hcount_change(react_graph: nx.Graph, prod_graph: nx.Graph) -> int:
    """
    Computes the maximum change in hydrogen count ('hcount') between corresponding nodes
    in the reactant and product graphs. It considers both hydrogen formation and breakage.

    Args:
    react_graph (nx.Graph): The graph representing reactants.
    prod_graph (nx.Graph): The graph representing products.

    Returns:
    int: The maximum hydrogen change observed across all nodes.
    """
    # max_hydrogen_change = 0
    hcount_break, _ = check_explicit_hydrogen(react_graph)
    hcount_form, _ = check_explicit_hydrogen(prod_graph)

    for node_id, attrs in react_graph.nodes(data=True):
        react_hcount = attrs.get("hcount", 0)
        if node_id in prod_graph:
            prod_hcount = prod_graph.nodes[node_id].get("hcount", 0)
        else:
            prod_hcount = 0

        if react_hcount >= prod_hcount:
            hcount_break += react_hcount - prod_hcount
        else:
            hcount_form += prod_hcount - react_hcount

        max_hydrogen_change = max(hcount_break, hcount_form)

    return max_hydrogen_change


def check_graph_connectivity(graph):
    """
    Check the connectivity of a graph.

    Parameters:
    - graph: a NetworkX graph object

    Returns:
    - A string indicating whether the graph is connected.
    """

    if nx.is_connected(graph):
        return "Connected"
    else:
        return "Disconnected."


def get_priority(
    its_list: List[Any],
    reaction_centers: List[Any],
    priority_ring: List[int] = [4, 5, 6],  # Standard priority rings
    priority_pair: List[int] = [3, 5],  # Special priority requiring both rings
    not_priority_ring: List[int] = [3],  # Non-priority ring that disqualifies alone
) -> Tuple[List[Any], List[Any]]:
    """
    Filters reaction centers based on their connectivity and specific ring sizes,
    including those with both rings in the priority pair, and excluding those with non-
    priority ring sizes,
    unless a specific pair condition is met (e.g., 3 must appear with 5).

    Parameters:
    - its_list (List[Any]): List of identifiers for the reaction centers.
    - reaction_centers (List[Any]): List of reaction centers to evaluate.
    - priority_ring (List[int], optional): List of ring sizes given priority.
    Defaults to [4, 6].
    - priority_pair (List[int], optional): List of two ring sizes that must both appear
    together to qualify. Defaults to [3, 5].
    - not_priority_ring (List[int], optional): List of ring sizes that disqualify a center
    unless paired appropriately. Defaults to [3].

    Returns:
    - Tuple[List[Any], List[Any]]: Tuple containing two lists:
        - The first list contains the identifiers from its_list that meet all criteria.
        - The second list contains the corresponding reaction centers that meet the
        criteria.
    """
    priority_set = set(priority_ring)
    not_priority_set = set(not_priority_ring)
    priority_pair_set = set(priority_pair)

    # Filter to include only connected reaction centers
    connected_centers = []
    connected_its_list = []
    for index, center in enumerate(reaction_centers):
        if check_graph_connectivity(center) == "Connected":
            connected_centers.append(center)
            connected_its_list.append(its_list[index])

    cyclic = [get_cycle_member_rings(center) for center in connected_centers]
    # Filter indices based on priority and non-priority ring sizes
    final_indices = []
    for i, rings in enumerate(cyclic):
        ring_set = set(rings)
        # Check for priority conditions and special conditions for non-priority rings
        if (
            ring_set.intersection(priority_set) or (priority_pair_set <= ring_set)
        ) and not (
            ring_set.intersection(not_priority_set)
            and not (5 in ring_set and 3 in ring_set)
        ):
            final_indices.append(i)

    # Retrieve final lists based on filtered indices
    final_its_list = [connected_its_list[i] for i in final_indices]
    final_centers = [connected_centers[i] for i in final_indices]

    return final_its_list, final_centers