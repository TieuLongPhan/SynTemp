import itertools
import networkx as nx
from operator import eq
from copy import deepcopy
from multiprocessing import Pool
from joblib import Parallel, delayed
from typing import Dict, List, Tuple, Iterable
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match

from synutility.SynIO.debug import setup_logging
from synutility.SynAAM.its_construction import ITSConstruction
from synutility.SynGraph.Descriptor.graph_signature import GraphSignature

from syntemp.SynRule.rules_extraction import RuleExtraction
from syntemp.SynITS.hydrogen_utils import (
    check_hcount_change,
    check_explicit_hydrogen,
    get_priority,
)


logger = setup_logging()


class ITSHAdjuster:
    """
    A class for infering hydrogen to complete reaction center or ITS graph.
    """

    @staticmethod
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
        nodeMatch = generic_node_match(
            nodeLabelNames, nodeLabelDefault, nodeLabelOperator
        )
        edgeMatch = generic_edge_match("order", 1, eq)

        classified = []

        for i in range(1, len(its_graphs)):
            # Compare the first graph with each subsequent graph
            if nx.is_isomorphic(
                its_graphs[0], its_graphs[i], node_match=nodeMatch, edge_match=edgeMatch
            ):
                classified.append((0, i))
        return classified, len(classified)

    @staticmethod
    def update_graph_data(
        graph_data: Dict,
        react_graph: nx.Graph,
        prod_graph: nx.Graph,
        its: nx.Graph,
    ) -> Dict:
        """
        Updates the graph data dictionary with new ITS and GraphRules based on the provided graphs.

        Parameters:
        - graph_data (Dict): Existing graph data dictionary.
        - react_graph (nx.Graph): Reactant graph.
        - prod_graph (nx.Graph): Product graph.
        - its (nx.Graph): Imaginary Transition State graph

        Returns:
        - Dict: Updated graph data dictionary with new ITS and GraphRules.
        """
        graph_data["ITSGraph"] = (react_graph, prod_graph, its)
        graph_data["GraphRules"] = RuleExtraction.extract_reaction_rules(
            react_graph, prod_graph, its, extend=False, n_knn=1
        )
        return graph_data

    @staticmethod
    def process_single_graph_data(
        graph_data: Dict,
        column: str,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
        get_priority_graph: bool = False,
    ) -> Dict:
        """
        Processes a single graph data dictionary, applying modifications based on
        hydrogen count changes. Optionally handles aromaticity and
        returns modified graph data.

        Parameters:
        - graph_data (Dict): Dictionary containing graph.
        - column (str): The key in the dictionary where the graph tuple is stored.
        - ignore_aromaticity (bool): Flag to indicate if aromaticity should be ignored.
        Default is False.
        - balance_its (bool): Flag to balance the ITS. Default is True.
        - get_priority_graph (bool): Flag to determine if priority graphs should be
        considered. Default is False.

        Returns:
        - Dict: Updated graph data dictionary, reflecting changes based on
        hydrogen counts and aromaticity.
        """
        graphs = deepcopy(graph_data)
        react_graph, prod_graph, its = graphs[column]
        is_empty_graph_present = any(
            (not isinstance(graph, nx.Graph) or graph.number_of_nodes() == 0)
            for graph in graphs[column]
        )

        if is_empty_graph_present:
            graphs["ITSGraph"], graphs["GraphRules"] = None, None
            return graphs

        hcount_change = check_hcount_change(react_graph, prod_graph)
        if hcount_change == 0:
            graph_data = ITSHAdjuster.update_graph_data(
                graphs, react_graph, prod_graph, its
            )
        else:
            graph_data = ITSHAdjuster.process_multiple_hydrogens(
                graphs,
                react_graph,
                prod_graph,
                ignore_aromaticity,
                balance_its,
                get_priority_graph,
            )
        if graph_data["GraphRules"] is not None:
            is_empty_rc_present = any(
                (not isinstance(graph, nx.Graph) or graph.number_of_nodes() == 0)
                for graph in graph_data["GraphRules"]
            )
            if is_empty_rc_present:
                graph_data["ITSGraph"] = None
                graph_data["GraphRules"] = None
        return graph_data

    @staticmethod
    def process_single_graph_data_safe(
        graph_data: Dict,
        column: str,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
        job_timeout: int = 1,
        get_priority_graph: bool = False,
    ) -> Dict:
        """
        Processes a single graph data dictionary asynchronously, handling potential
        timeouts during processing.

        Parameters:
        - graph_data (Dict): Dictionary containing graph data.
        - column (str): Key to access the graph data.
        - ignore_aromaticity (bool): Flag to ignore aromaticity during processing.
        Default is False.
        - balance_its (bool): Flag to balance the ITS. Default is True.
        - job_timeout (int): Timeout in seconds for the asynchronous task.
        Default is 1 second.
        - get_priority_graph (bool): Flag to include priority graph processing.
        Default is False.

        Returns:
        - Dict: Processed graph data dictionary.
        """
        pool = Pool(processes=1)
        try:
            async_result = pool.apply_async(
                ITSHAdjuster.process_single_graph_data,
                (
                    graph_data,
                    column,
                    ignore_aromaticity,
                    balance_its,
                    get_priority_graph,
                ),
            )
            graph_data = async_result.get(job_timeout)
            pool.terminate()
            pool.join()
        except Exception as e:
            logger.error(
                f"Issue processing graph data: {e} with time-out {job_timeout}s",
                exc_info=True,
            )
            graph_data["ITSGraph"], graph_data["GraphRules"] = None, None
            pool.terminate()  # Terminate the problematic pool.
            pool.join()
        finally:
            pool.close()
            pool.join()
        return graph_data

    def process_graph_data_parallel(
        self,
        graph_data_list: List[Dict],
        column: str,
        n_jobs: int,
        verbose: int,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
        job_timeout: int = 5,
        safe: bool = False,
        get_priority_graph: bool = False,
    ) -> List[Dict]:
        """
        Processes a list of graph data dictionaries in parallel, utilizing multiple jobs
        for faster processing.

        Parameters:
        - graph_data_list (List[Dict]): List of dictionaries containing graph data.
        - column (str): Key where the graph data is stored.
        - n_jobs (int): Number of parallel jobs to run.
        - verbose (int): Verbosity level for the parallel process.
        - ignore_aromaticity (bool): Flag to ignore aromaticity during processing.
        Default is False.
        - balance_its (bool): Flag to balance ITS.
        Default is True.
        - job_timeout (int): Timeout for job processing in seconds.
        Default is 5 seconds.
        - safe (bool): Flag to use safe parallel processing (timeout).
        Default is False.
        - get_priority_graph (bool): Flag to prioritize graphs.
        Default is False.

        Returns:
        - List[Dict]: A list of processed graph data dictionaries.
        """
        if safe:
            processed_data = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(self.process_single_graph_data_safe)(
                    graph_data,
                    column,
                    ignore_aromaticity,
                    balance_its,
                    job_timeout,
                    get_priority_graph,
                )
                for graph_data in graph_data_list
            )
        else:
            processed_data = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(self.process_single_graph_data)(
                    graph_data,
                    column,
                    ignore_aromaticity,
                    balance_its,
                    get_priority_graph,
                )
                for graph_data in graph_data_list
            )

        return processed_data

    @staticmethod
    def process_multiple_hydrogens(
        graph_data: Dict,
        react_graph: nx.Graph,
        prod_graph: nx.Graph,
        ignore_aromaticity: bool,
        balance_its: bool,
        get_priority_graph: bool = False,
    ) -> Dict:
        """
        Handles cases where hydrogen counts change significantly between the reactant
        and product graphs. Adjusts hydrogen nodes accordingly and evaluates equivalence.

        Parameters:
        - graph_data (Dict): Dictionary of graph data.
        - react_graph (nx.Graph): Reactant graph.
        - prod_graph (nx.Graph): Product graph.
        - ignore_aromaticity (bool): Flag to ignore aromaticity. Default is False.
        - balance_its (bool): Flag to balance ITS. Default is True.
        - get_priority_graph (bool): Flag to prioritize graph processing. Default is False.

        Returns:
        - Dict: Updated graph data after handling hydrogen node adjustments.
        """
        combinations_solution = ITSHAdjuster.add_hydrogen_nodes_multiple(
            react_graph, prod_graph
        )
        # print("length original:", len(combinations_solution))
        its_list = [
            ITSConstruction.ITSGraph(
                i[0], i[1], ignore_aromaticity, balance_its=balance_its
            )
            for i in combinations_solution
        ]
        # rc_list = [get_rc(i) for i in its_list]
        rc_list = [
            RuleExtraction.extract_reaction_rules(react_graph, prod_graph, i, False, 1)[
                2
            ]
            for i in its_list
        ]
        rc_list = [
            rc
            for rc in rc_list
            if rc is not None and isinstance(rc, nx.Graph) and rc.number_of_nodes() > 0
        ]
        rc_sig = [GraphSignature(i).create_graph_signature() for i in rc_list]
        combinations_solution = [
            comb
            for rc, comb in zip(rc_list, combinations_solution)
            if rc is not None and isinstance(rc, nx.Graph) and rc.number_of_nodes() > 0
        ]
        its_list = [
            ITSConstruction.ITSGraph(
                i[0], i[1], ignore_aromaticity, balance_its=balance_its
            )
            for i in combinations_solution
        ]
        if len(rc_sig) != len(set(rc_sig)):
            _, equivariant = ITSHAdjuster.check_equivariant_graph(rc_list)
        else:
            equivariant = 0

        pairwise_combinations = len(rc_list) - 1
        if equivariant == pairwise_combinations:
            graph_data = ITSHAdjuster.update_graph_data(
                graph_data, *combinations_solution[0], its_list[0]
            )
        else:
            graph_data["ITSGraph"], graph_data["GraphRules"] = None, None
            if get_priority_graph:
                priority_indices = get_priority(rc_list)
                rc_list = [rc_list[i] for i in priority_indices]
                rc_sig = [rc_sig[i] for i in priority_indices]
                its_list = [its_list[i] for i in priority_indices]
                combinations_solution = [
                    combinations_solution[i] for i in priority_indices
                ]
                if len(rc_sig) != len(set(rc_sig)):
                    _, equivariant = ITSHAdjuster.check_equivariant_graph(rc_list)
                pairwise_combinations = len(rc_list) - 1
                if equivariant == pairwise_combinations:
                    graph_data = ITSHAdjuster.update_graph_data(
                        graph_data, *combinations_solution[0], its_list[0]
                    )
        return graph_data

    @staticmethod
    def add_hydrogen_nodes_multiple(
        react_graph: nx.Graph,
        prod_graph: nx.Graph,
    ) -> List[Tuple[nx.Graph, nx.Graph]]:
        """
        Adds hydrogen nodes to the copies of the reactant and product graphs based on the
        difference in hcount between them. Hydrogen nodes are added or removed to
        represent the breaking and forming of hydrogen bonds. The function generates
        multiple graph pairs, each with a different permutation of the added hydrogen
        nodes in the product graph.

        Parameters:
        - react_graph (nx.Graph): The reactant graph.
        - prod_graph (nx.Graph): The product graph.

        Returns:
        - List[Tuple[nx.Graph, nx.Graph]]: A list of tuples, each containing a pair of
        updated reactant and product graphs.
        """
        react_graph_copy = deepcopy(react_graph)
        prod_graph_copy = deepcopy(prod_graph)
        react_explicit_h, hydrogen_nodes = check_explicit_hydrogen(react_graph_copy)
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
        range_implicit_h = range(
            max_index + 1,
            max_index + 1 + len(hydrogen_nodes_form) - react_explicit_h,
        )
        combined_indices = list(range_implicit_h) + hydrogen_nodes
        permutations = list(itertools.permutations(combined_indices))
        permutations_seed = permutations[0]

        updated_graphs = []
        for permutation in permutations:
            current_react_graph, current_prod_graph = deepcopy(
                react_graph_copy
            ), deepcopy(prod_graph_copy)

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

    @staticmethod
    def add_hydrogen_nodes_multiple_utils(
        graph: nx.Graph,
        node_id_pairs: Iterable[Tuple[int, int]],
        atom_map_update: bool = True,
    ) -> nx.Graph:
        """
        Creates and returns a new graph with added hydrogen nodes based on the input graph
        and node ID pairs.

        Parameters:
        - graph (nx.Graph): The base graph to which the nodes will be added.
        - node_id_pairs (Iterable[Tuple[int, int]]): Pairs of node IDs (original node, new
        hydrogen node) to link with hydrogen.
        - atom_map_update (bool): If True, update the 'atom_map' attribute with the new
        hydrogen node ID; otherwise, retain the original node's 'atom_map'.

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
