import networkx as nx
from typing import Dict, List, Tuple, Iterable
import itertools
from joblib import Parallel, delayed
from copy import deepcopy
from SynTemp.SynRule.rules_extraction import RuleExtraction
from SynTemp.SynITS.its_construction import ITSConstruction
from SynTemp.SynUtils.graph_utils import (
    check_hcount_change,
    check_explicit_hydrogen,
    get_priority,
)
from SynTemp.SynITS.its_extraction import ITSExtraction

import logging


class ITSHAdjuster:

    @staticmethod
    def process_single_graph_data(
        graph_data: Dict,
        column: str,
        return_all: bool = False,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
    ) -> Dict:
        """
        Processes a single dictionary containing graph information by applying modifications based on hcount changes.
        Optionally handles aromaticity and provides different return behaviors based on the `return_all` flag.

        Args:
        graph_data (Dict): Dictionary containing graph information.
        column (str): Key where the graph tuple is stored.
        return_all (bool): Flag to return all data or just the most relevant. Default is False.
        ignore_aromaticity (bool): Flag to ignore aromaticity in the ITS construction. Default is False.

        Returns:
        Dict: Updated dictionary with new ITS and GraphRules if applicable.
        """
        graphs = deepcopy(graph_data)
        react_graph, prod_graph, its = graphs[column]

        hcount_change = check_hcount_change(react_graph, prod_graph)
        if hcount_change == 0:
            graph_data = ITSHAdjuster.update_graph_data(
                graphs, react_graph, prod_graph, its
            )
        elif hcount_change < 5:
            graph_data = ITSHAdjuster.process_multiple_hydrogens(
                graphs,
                react_graph,
                prod_graph,
                its,
                ignore_aromaticity,
                return_all,
                balance_its,
            )
        else:
            graph_data = ITSHAdjuster.process_high_hcount_change(
                graphs,
                react_graph,
                prod_graph,
                its,
                ignore_aromaticity,
                return_all,
                balance_its,
            )

        return graph_data

    @staticmethod
    def update_graph_data(graph_data, react_graph, prod_graph, its):
        """
        Update graph data dictionary with new ITS and GraphRules based on the graphs provided.

        Args:
        graph_data (Dict): Existing graph data.
        react_graph (nx.Graph), prod_graph (nx.Graph), its: Graphs and ITS to use.

        Returns:
        Dict: Updated graph data dictionary.
        """
        graph_data["ITSGraph"] = (react_graph, prod_graph, its)
        graph_data["GraphRules"] = RuleExtraction.extract_reaction_rules(
            react_graph, prod_graph, its, extend=False, n_knn=1
        )
        return graph_data

    @staticmethod
    def process_multiple_hydrogens(
        graph_data,
        react_graph,
        prod_graph,
        its,
        ignore_aromaticity,
        return_all,
        balance_its,
    ):
        """
        Handles cases with hydrogen count changes between 2 and 4, inclusive.
        Manages the creation of multiple hydrogen node scenarios and evaluates their equivalence.

        Args:
        graph_data, react_graph, prod_graph, its, ignore_aromaticity, return_all as described.

        Returns:
        Dict: Updated graph data.
        """
        combinations_solution = ITSHAdjuster.add_hydrogen_nodes_multiple(
            react_graph, prod_graph
        )
        its_list = [
            ITSConstruction.ITSGraph(
                i[0], i[1], ignore_aromaticity, balance_its=balance_its
            )
            for i in combinations_solution
        ]
        _, equivariant = ITSExtraction.check_equivariant_graph(its_list)
        pairwise_combinations = len(its_list) - 1
        if equivariant == pairwise_combinations:
            graph_data = ITSHAdjuster.update_graph_data(
                graph_data, *combinations_solution[0], its_list[0]
            )
        else:
            graph_data = ITSHAdjuster.process_high_hcount_change(
                graph_data,
                react_graph,
                prod_graph,
                its,
                ignore_aromaticity,
                return_all,
                balance_its,
            )
        return graph_data

    @staticmethod
    def process_high_hcount_change(
        graph_data,
        react_graph,
        prod_graph,
        its,
        ignore_aromaticity,
        return_all,
        balance_its: bool = True,
    ):
        """
        Handles cases with hydrogen count changes of 5 or more.
        Similar to `process_multiple_hydrogens` but tailored for higher counts.

        Args:
        Same as `process_multiple_hydrogens`.

        Returns:
        Dict: Updated graph data.
        """
        combinations_solution = ITSHAdjuster.add_hydrogen_nodes_multiple(
            react_graph, prod_graph
        )
        its_list = [
            ITSConstruction.ITSGraph(
                i[0], i[1], ignore_aromaticity, balance_its=balance_its
            )
            for i in combinations_solution
        ]
        reaction_centers = [
            RuleExtraction.extract_reaction_rules(react_graph, prod_graph, i)[2]
            for i in its_list
        ]

        its_list, rc_list = get_priority(its_list, reaction_centers)
        _, equivariant = ITSExtraction.check_equivariant_graph(rc_list)
        pairwise_combinations = len(its_list) - 1
        if equivariant == pairwise_combinations:
            graph_data = ITSHAdjuster.update_graph_data(
                graph_data, *combinations_solution[0], its_list[0]
            )
        else:
            if return_all:
                graph_data = ITSHAdjuster.update_graph_data(
                    graph_data, react_graph, prod_graph, its
                )
            else:
                graph_data["ITSGraph"], graph_data["GraphRules"] = None, None
        return graph_data

    @staticmethod
    def process_graph_data_parallel(
        graph_data_list: List[Dict],
        column: str,
        n_jobs: int,
        verbose: int,
        return_all: bool = False,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
    ) -> List[Dict]:
        """
        Processes a list of dictionaries containing graph information in parallel.

        Parameters:
        - graph_data_list (List[Dict]): A list of dictionaries containing graph information.
        - column (str): The key in the dictionary where the graph tuple is stored.
        - n_jobs (int): The number of concurrent jobs.
        - verbose (int): The verbosity level.

        Returns:
        - List[Dict]: A list of dictionaries with the updated graph data.
        """
        processed_data = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(ITSHAdjuster.process_single_graph_data)(
                graph_data, column, return_all, ignore_aromaticity, balance_its
            )
            for graph_data in graph_data_list
        )

        return processed_data

    @staticmethod
    def add_hydrogen_nodes_multiple_utils(
        graph: nx.Graph,
        node_id_pairs: Iterable[Tuple[int, int]],
        atom_map_update: bool = True,
    ) -> nx.Graph:
        """
        Creates and returns a new graph with added hydrogen nodes based on the input graph and node ID pairs.

        Parameters:
        - graph (nx.Graph): The base graph to which the nodes will be added.
        - node_id_pairs (Iterable[Tuple[int, int]]): Pairs of node IDs (original node, new hydrogen node) to link with hydrogen.
        - atom_map_update (bool): If True, update the 'atom_map' attribute with the new hydrogen node ID; otherwise, retain the original node's 'atom_map'.

        Returns:
        - nx.Graph: A new graph instance with the added hydrogen nodes.
        """
        new_graph = deepcopy(graph)
        for node_id, new_hydrogen_node_id in node_id_pairs:
            atom_map_val = (
                new_hydrogen_node_id
                if atom_map_update
                else new_graph.nodes[node_id].get("atom_map", 0)
            )
            new_graph.add_node(
                new_hydrogen_node_id,
                charge=0,
                hcount=0,
                aromatic=False,
                element="H",
                atom_map=atom_map_val,
                isomer="N",
                partial_charge=0,
                hybridization=0,
                in_ring=False,
                explicit_valence=0,
                implicit_hcount=0,
            )
            new_graph.add_edge(
                node_id,
                new_hydrogen_node_id,
                order=1.0,
                ez_isomer="N",
                bond_type="SINGLE",
                conjugated=False,
                in_ring=False,
            )
        return new_graph

    @staticmethod
    def add_hydrogen_nodes_multiple(
        react_graph: nx.Graph,
        prod_graph: nx.Graph,
    ) -> List[Tuple[nx.Graph, nx.Graph]]:
        """
        Adds hydrogen nodes to the copies of the reactant and product graphs based on the difference in hcount between them.
        Hydrogen nodes are added or removed to represent the breaking and forming of hydrogen bonds.
        The function generates multiple graph pairs, each with a different permutation of the added hydrogen nodes in the product graph.

        Parameters:
        - react_graph (nx.Graph): The reactant graph.
        - prod_graph (nx.Graph): The product graph.

        Returns:
        - List[Tuple[nx.Graph, nx.Graph]]: A list of tuples, each containing a pair of updated reactant and product graphs.
        """
        react_graph_copy = deepcopy(react_graph)
        prod_graph_copy = deepcopy(prod_graph)
        react_explicit_h, _ = check_explicit_hydrogen(react_graph_copy)
        prod_explicit_h, _ = check_explicit_hydrogen(prod_graph_copy)
        hydrogen_nodes_form, hydrogen_nodes_break = [], []

        primary_graph = (
            react_graph_copy if react_explicit_h <= prod_explicit_h else prod_graph_copy
        )
        for node_id in primary_graph.nodes:
            try:
                # Calculate the difference in hydrogen counts
                hcount_diff = react_graph_copy.nodes[node_id].get(
                    "hcount", 0
                ) - prod_graph_copy.nodes[node_id].get("hcount", 0)
            except KeyError:
                # Handle cases where node_id does not exist in opposite_graph
                continue

            # Decide action based on hcount_diff
            if hcount_diff > 0:
                hydrogen_nodes_break.extend([node_id] * hcount_diff)
            elif hcount_diff < 0:
                hydrogen_nodes_form.extend([node_id] * -hcount_diff)

        max_index = max(
            max(react_graph_copy.nodes, default=0),
            max(prod_graph_copy.nodes, default=0),
        )
        permutations = list(
            itertools.permutations(
                range(
                    max_index + 1 - react_explicit_h,
                    max_index + 1 + len(hydrogen_nodes_form) - react_explicit_h,
                )
            )
        )
        permutations_seed = list(
            itertools.permutations(
                range(
                    max_index + 1 - prod_explicit_h,
                    max_index + 1 + len(hydrogen_nodes_break) - prod_explicit_h,
                )
            )
        )[0]

        updated_graphs = []
        for permutation in permutations:
            current_react_graph, current_prod_graph = deepcopy(
                react_graph_copy
            ), deepcopy(prod_graph_copy)

            # Correctly form the list for new hydrogen node IDs by adding react_explicit_h to each element in permutations_seed
            new_hydrogen_node_ids = [i for i in permutations_seed]

            # Use `zip` to pair `hydrogen_nodes_break` with the new IDs
            node_id_pairs = zip(hydrogen_nodes_break, new_hydrogen_node_ids)
            # Call the method with the formed pairs and specify atom_map_update as False
            current_react_graph = ITSHAdjuster.add_hydrogen_nodes_multiple_utils(
                current_react_graph, node_id_pairs, atom_map_update=False
            )
            # Varied hydrogen nodes in the product graph based on permutation
            current_prod_graph = ITSHAdjuster.add_hydrogen_nodes_multiple_utils(
                current_prod_graph, zip(hydrogen_nodes_form, permutation)
            )
            updated_graphs.append((current_react_graph, current_prod_graph))
        return updated_graphs
