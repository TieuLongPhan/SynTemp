import networkx as nx
from networkx.algorithms.cycles import simple_cycles
from networkx.algorithms.planarity import check_planarity
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
    Determines if the given graph is a single cyclic graph, which means the graph has exactly one cycle.
    
    Parameters:
    - G (nx.Graph): The graph to be checked.
    
    Returns:
    - bool: True if the graph has exactly one cycle, False otherwise.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    # The graph must be connected and have exactly one cycle to be single cyclic
    if nx.is_connected(G):
        cycles = list(simple_cycles(G))
        return len(cycles) == 1 and all(len(cycle) > 2 for cycle in cycles)
    return False

def is_complex_cyclic_graph(G: nx.Graph) -> bool:
    """
    Determines if the given graph is a complex cyclic graph, which means the graph has more than one cycle.

    Parameters:
    - G (nx.Graph): The graph to be checked.

    Returns:
    - bool: True if the graph has more than one cycle, False otherwise.
    """
    if not isinstance(G, nx.Graph):
        raise TypeError("Input must be a networkx Graph object.")

    # Check for more than one cycle
    if sum(1 for _ in simple_cycles(G)) > 1:
        return True
    else:
        return False

def check_graph_type(G: nx.Graph) -> str:
    """
    Determines if the given graph is acyclic, single cyclic, or complex cyclic.

    Parameters:
    - G (nx.Graph): The graph to be checked.

    Returns:
    - str: A string indicating if the graph is "Acyclic", "Single Cyclic", or "Complex Cyclic".

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
    Identifies all cycles in the given graph using cycle bases to ensure no overlap and returns a list of the sizes of these cycles (member rings),
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