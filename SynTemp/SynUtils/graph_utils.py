import networkx as nx
from typing import List


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
