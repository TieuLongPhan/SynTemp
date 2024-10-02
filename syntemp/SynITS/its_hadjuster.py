import networkx as nx
from typing import Dict, List, Tuple, Iterable
import itertools
from joblib import Parallel, delayed
from copy import deepcopy
from syntemp.SynRule.rules_extraction import RuleExtraction
from syntemp.SynITS.its_construction import ITSConstruction
from syntemp.SynUtils.graph_utils import (
    check_hcount_change,
    check_explicit_hydrogen,
    get_priority,
)
from syntemp.SynITS.its_extraction import ITSExtraction
from syntemp.SynUtils.utils import setup_logging
from multiprocessing import Pool

logger = setup_logging()


class ITSHAdjuster:

    @staticmethod
    def process_single_graph_data(
        graph_data: Dict,
        column: str,
        return_all: bool = False,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
        get_random_results=False,
        fast_process: bool = False,
    ) -> Dict:
        """
        Processes a single dictionary containing graph information by applying
        modifications based on hcount changes.
        Optionally handles aromaticity and provides different return behaviors based on
        the `return_all` flag.

        Parameters:
        - graph_data (Dict): A dictionary containing essential graph information. This
        includes nodes, edges, and other graph-specific data.
        - column (str): The key in the dictionary where the graph tuple is stored,
        typically pointing to the specific data structure to be modified.
        - return_all (bool): A flag that determines the nature of the output. If True, the
        function returns all modified data, otherwise it returns only the most relevant
        changes. The default value is False.
        - ignore_aromaticity (bool): A flag to indicate whether aromaticity should be
        ignored during the graph processing. Ignoring aromaticity may affect the ITS
        construction. The default value is False.

        Returns:
        - Dict: An updated dictionary that includes the new internal topology structure
        (ITS) and any applicable GraphRules, reflecting the modifications made based on
        hydrogen counts and aromaticity considerations.
        """
        graphs = deepcopy(graph_data)
        logger.info(f"{graphs}")
        react_graph, prod_graph, its = graphs[column]
        is_empty_graph_present = any(
            (not isinstance(graph, nx.Graph) or graph.number_of_nodes() == 0)
            for graph in graphs[column]
        )

        if is_empty_graph_present:
            # Update graph data if any graph is empty
            graphs["ITSGraph"], graphs["GraphRules"] = None, None
            return graphs

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
                get_random_results,
            )
        else:
            if fast_process:
                graphs["ITSGraph"], graphs["GraphRules"] = None, None
                return graphs
            else:
                graph_data = ITSHAdjuster.process_high_hcount_change(
                    graphs,
                    react_graph,
                    prod_graph,
                    its,
                    ignore_aromaticity,
                    return_all,
                    balance_its,
                    get_random_results,
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
    def update_graph_data(graph_data, react_graph, prod_graph, its, ignore=False):
        """
        Update graph data dictionary with new ITS and GraphRules based on the graphs
        provided.

        Parameters:
        - graph_data (Dict): Existing graph data.
        - react_graph (nx.Graph), prod_graph (nx.Graph), its: Graphs and ITS to use.

        Returns:
        Dict: Updated graph data dictionary.
        """
        graph_data["ITSGraph"] = (react_graph, prod_graph, its)
        graph_data["GraphRules"] = RuleExtraction.extract_reaction_rules(
            react_graph, prod_graph, its, extend=False, n_knn=1
        )
        if ignore:
            graph_data["ITSGraph"], graph_data["GraphRules"] = None, None
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
        get_random_results=False,
    ):
        """
        Handles cases with hydrogen count changes between 2 and 4, inclusive.
        Manages the creation of multiple hydrogen node scenarios and evaluates their
        equivalence.

        Parameters:
        - graph_data, react_graph, prod_graph, its, ignore_aromaticity,
        return_all as described.

        Returns:
        - Dict: Updated graph data.
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
                get_random_results,
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
        get_random_results=False,
    ):
        """
        Handles cases with hydrogen count changes of 5 or more.
        Similar to `process_multiple_hydrogens` but tailored for higher counts.

        Parameters:
        - Same as `process_multiple_hydrogens`.

        Returns:
        - Dict: Updated graph data.
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

        filtered_reaction_centers = [
            rc
            for rc in reaction_centers
            if rc is not None and isinstance(rc, nx.Graph) and rc.number_of_nodes() > 0
        ]

        filtered_combinations_solution = [
            comb
            for rc, comb in zip(reaction_centers, combinations_solution)
            if rc is not None and isinstance(rc, nx.Graph) and rc.number_of_nodes() > 0
        ]

        # Update the original lists with the filtered results
        reaction_centers, combinations_solution = (
            filtered_reaction_centers,
            filtered_combinations_solution,
        )

        priority_indices = get_priority(reaction_centers)
        rc_list = [reaction_centers[i] for i in priority_indices]
        its_list = [its_list[i] for i in priority_indices]
        combinations_solution = [combinations_solution[i] for i in priority_indices]
        _, equivariant = ITSExtraction.check_equivariant_graph(rc_list)
        pairwise_combinations = len(its_list) - 1

        if equivariant == pairwise_combinations:

            graph_data = ITSHAdjuster.update_graph_data(
                graph_data, *combinations_solution[0], its_list[0]
            )

        else:
            if get_random_results is True:
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
    def process_single_graph_data_safe(
        graph_data: Dict,
        column: str,
        return_all: bool = False,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
        get_random_results=False,
        fast_process: bool = False,
        job_timeout: int = 1,
    ) -> Dict:
        # pool = multiprocessing.pool.ThreadPool(1)
        pool = Pool(processes=1)
        try:
            async_result = pool.apply_async(
                ITSHAdjuster.process_single_graph_data,
                (
                    graph_data,
                    column,
                    return_all,
                    ignore_aromaticity,
                    balance_its,
                    get_random_results,
                    fast_process,
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

    @staticmethod
    def process_graph_data_parallel(
        graph_data_list: List[Dict],
        column: str,
        n_jobs: int,
        verbose: int,
        return_all: bool = False,
        ignore_aromaticity: bool = False,
        balance_its: bool = True,
        get_random_results: bool = False,
        fast_process: bool = False,
        job_timeout: int = 1,
    ) -> List[Dict]:
        """
        Processes a list of dictionaries containing graph information in parallel.

        Parameters:
        - graph_data_list (List[Dict]): A list of dictionaries containing graph
        information.
        - column (str): The key in the dictionary where the graph tuple is stored.
        - n_jobs (int): The number of concurrent jobs.
        - verbose (int): The verbosity level.

        Returns:
        - List[Dict]: A list of dictionaries with the updated graph data.
        """
        processed_data = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(ITSHAdjuster.process_single_graph_data_safe)(
                graph_data,
                column,
                return_all,
                ignore_aromaticity,
                balance_its,
                get_random_results,
                fast_process,
                job_timeout,
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
