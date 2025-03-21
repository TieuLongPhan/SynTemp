import networkx as nx
from operator import eq
from typing import List, Any, Tuple
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match

from synkit.Graph.Feature.graph_descriptors import GraphDescriptor


def check_equivariant_graph(
    its_graphs: List[nx.Graph],
) -> Tuple[List[Tuple[int, int]], int]:
    """
    Checks for isomorphism among a list of ITS graphs.

    Parameters:
    - its_graphs (List[nx.Graph]): A list of ITS graphs.

    Returns:
    - List[Tuple[int, int]]: A list of tuples representing pairs of indices of
    isomorphic graphs.
    """
    nodeLabelNames = ["typesGH"]
    nodeLabelDefault = [()]
    nodeLabelOperator = [eq]
    nodeMatch = generic_node_match(nodeLabelNames, nodeLabelDefault, nodeLabelOperator)
    edgeMatch = generic_edge_match("order", 1, eq)

    classified = []

    for i in range(1, len(its_graphs)):
        # Compare the first graph with each subsequent graph
        if nx.is_isomorphic(
            its_graphs[0], its_graphs[i], node_match=nodeMatch, edge_match=edgeMatch
        ):
            classified.append((0, i))
    return classified, len(classified)


def check_explicit_hydrogen(graph: nx.Graph) -> tuple:
    """
    Counts the explicit hydrogen nodes in the given graph and collects their IDs.

    Parameters:
    - graph (nx.Graph): The graph to inspect.

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

    Parameters:
    - react_graph (nx.Graph): The graph representing reactants.
    - prod_graph (nx.Graph): The graph representing products.

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


def get_cycle_member_rings(G: nx.Graph, type="minimal") -> List[int]:
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

    if type == "minimal":
        cycles = nx.minimum_cycle_basis(G)
    else:
        cycles = nx.cycle_basis(G)
    member_rings = [len(cycle) for cycle in cycles]

    member_rings.sort()

    return member_rings


def get_priority(reaction_centers: List[Any]) -> List[int]:
    """
    Evaluate reaction centers for specific graph characteristics, selecting indices based
    on the shortest reaction paths and maximum ring sizes, and adjusting for certain
    graph types by modifying the ring information.

    Parameters:
    - reaction_centers: List[Any], a list of reaction centers where each center should be
    capable of being analyzed for graph type and ring sizes.

    Returns:
    - List[int]: A list of indices from the original list of reaction centers that meet
    the criteria of having the shortest reaction steps and/or the largest ring sizes.
    Returns indices with minimum reaction steps if no indices meet both criteria.
    """
    # Extract topology types and ring sizes from reaction centers
    topo_type = [
        GraphDescriptor.check_graph_type(center) for center in reaction_centers
    ]
    cyclic = [
        get_cycle_member_rings(center, "fundamental") for center in reaction_centers
    ]

    # Adjust ring information based on the graph type
    for index, graph_type in enumerate(topo_type):
        if graph_type in ["Acyclic", "Complex Cyclic"]:
            cyclic[index] = [0] + cyclic[index]

    # Determine minimum reaction steps
    reaction_steps = [len(rings) for rings in cyclic]
    min_reaction_step = min(reaction_steps)

    # Filter indices with the minimum reaction steps
    indices_shortest = [
        i for i, steps in enumerate(reaction_steps) if steps == min_reaction_step
    ]

    # Filter indices with the maximum ring size
    max_size = max(
        max(rings) for rings in cyclic if rings
    )  # Safeguard against empty sublists
    prior_indices = [i for i, rings in enumerate(cyclic) if max(rings) == max_size]

    # Combine criteria for final indices
    final_indices = [index for index in prior_indices if index in indices_shortest]

    # Fallback to shortest indices if no indices meet both criteria
    if not final_indices:
        return indices_shortest

    return final_indices
