import networkx as nx
from typing import Dict, List, Tuple, Iterable
from joblib import Parallel, delayed
from copy import deepcopy
from SynTemp.SynITS.graph_rules_extraction import GraphRuleExtraction
from SynTemp.SynITS.its_construction import ITSConstruction
from SynTemp.SynITS.its_extraction import ITSExtraction
import itertools
import math

class ITSHAdjuster:
    
    @staticmethod
    def check_hcount_change(react_graph: nx.Graph, prod_graph: nx.Graph) -> int:
        """
        Checks the 'hcount' attribute of nodes across reactant and product graphs to determine the number
        of hydrogen changes. Hydrogen formation is considered when the hcount decreases by 1 from reactant
        to product, and a hydrogen break is considered when the hcount increases by 1.

        Parameters:
        - react_graph (nx.Graph): The reactant graph.
        - prod_graph (nx.Graph): The product graph.

        Returns:
        - int: Returns the number of hydrogen changes.
        """
        hydrogen_formations = 0
        hydrogen_breaks = 0
        
        for node in react_graph.nodes(data=True):
            node_id = node[0]
            react_hcount = node[1].get('hcount', 0)
            
            if node_id in prod_graph:
                prod_hcount = prod_graph.nodes[node_id].get('hcount', 0)
                hcount_diff = react_hcount - prod_hcount
                    
                if hcount_diff == -1:
                    hydrogen_formations += 1
                    
                if hcount_diff == 1:
                    hydrogen_breaks += 1

        max_hydrogen_change = max(hydrogen_formations, hydrogen_breaks)
        return max_hydrogen_change


    @staticmethod
    def add_hydrogen_nodes(react_graph: nx.Graph, prod_graph: nx.Graph) -> (nx.Graph, nx.Graph):
        """
        Adds a single hydrogen node to the copies of the reactant and product graph based on the difference in hcount between them.
        Only one hydrogen node is added to each graph to represent the breaking and forming of a hydrogen bond.
        The hydrogen nodes added to both graphs will have the same node index.
        The original graphs are not modified.

        Parameters:
        - react_graph (nx.Graph): The reactant graph.
        - prod_graph (nx.Graph): The product graph.

        Returns:
        - Tuple[nx.Graph, nx.Graph]: The updated copies of the reactant and product graphs.
        """
        react_graph_copy = deepcopy(react_graph)
        prod_graph_copy = deepcopy(prod_graph)

        max_index = max(max(react_graph_copy.nodes, default=0), max(prod_graph_copy.nodes, default=0)) + 1
        react_node_to_add_hydrogen = None
        prod_node_to_add_hydrogen = None

        for node_id, attr in react_graph_copy.nodes(data=True):
            if react_node_to_add_hydrogen and prod_node_to_add_hydrogen:
                break

            react_hcount = attr.get('hcount', 0)
            
            if node_id in prod_graph_copy.nodes:
                prod_hcount = prod_graph_copy.nodes[node_id].get('hcount', 0)
                hcount_diff = react_hcount - prod_hcount

                if hcount_diff == 1 and not react_node_to_add_hydrogen:
                    react_node_to_add_hydrogen = node_id

                elif hcount_diff == -1 and not prod_node_to_add_hydrogen:
                    prod_node_to_add_hydrogen = node_id

        if react_node_to_add_hydrogen:
            react_graph_copy.add_node(max_index, charge=0, hcount=0, aromatic=False, element='H', atom_map=max_index,
                                      isomer='N', partial_charge=0, hybridization=0, in_ring=False,
                                      explicit_valence=0, implicit_hcount=0)
            react_graph_copy.add_edge(react_node_to_add_hydrogen, max_index, order=1.0, ez_isomer='N', bond_type='SINGLE',
                                      conjugated=False, in_ring=False)

        if prod_node_to_add_hydrogen:
            prod_graph_copy.add_node(max_index, charge=0, hcount=0, aromatic=False, element='H', atom_map=max_index,
                                     isomer='N', partial_charge=0, hybridization=0, in_ring=False,
                                     explicit_valence=0, implicit_hcount=0)
            prod_graph_copy.add_edge(prod_node_to_add_hydrogen, max_index, order=1.0, ez_isomer='N', bond_type='SINGLE',
                                     conjugated=False, in_ring=False)

        return react_graph_copy, prod_graph_copy

    @staticmethod
    def process_single_graph_data(graph_data: Dict, column: str, return_all: bool=False, ignore_aromaticity: bool = False) -> Dict:
        """
        Processes a single dictionary containing graph information, applying modifications based on hcount changes.

        Parameters:
        - graph_data (Dict): A dictionary containing graph information.
        - column (str): The key in the dictionary where the graph tuple is stored.

        Returns:
        - Dict: The updated dictionary with new ITS and GraphRules if applicable.
        """
        graphs = graph_data[column]
        react_graph, prod_graph, its = graphs

        hcount_change = ITSHAdjuster.check_hcount_change(react_graph, prod_graph)

        if hcount_change == 0:
            graph_data['ITSGraph']= (react_graph, prod_graph, its)
            graph_data['GraphRules'] = GraphRuleExtraction.extract_reaction_rules(react_graph, prod_graph, its, extend=False, n_knn=1)

        elif hcount_change == 1:
            react_graph, prod_graph = ITSHAdjuster.add_hydrogen_nodes(react_graph.copy(), prod_graph.copy())
            its = ITSConstruction.ITSGraph(react_graph, prod_graph, ignore_aromaticity)
            graph_data['ITSGraph'] = (react_graph, prod_graph, its)
            graph_data['GraphRules'] = GraphRuleExtraction.extract_reaction_rules(react_graph, prod_graph, its, extend=False, n_knn=1)
        elif 1 < hcount_change < 5:
            combinations_solution = ITSHAdjuster.add_hydrogen_nodes_multiple(react_graph, prod_graph)
            its_list = [ITSConstruction.ITSGraph(i[0], i[1], ignore_aromaticity) for i in combinations_solution]
            _, equivariant = ITSExtraction.check_equivariant_graph(its_list)
            #pairwise_combinations = math.comb(len(its_list),2)
            pairwise_combinations = len(its_list)-1
            if equivariant==pairwise_combinations:
                graph_data['ITSGraph']= (combinations_solution[0][0], combinations_solution[0][1], its_list[0])
                graph_data['GraphRules'] = GraphRuleExtraction.extract_reaction_rules(combinations_solution[0][0], combinations_solution[0][1], its_list[0], extend=False)
            else:
                if return_all:
                    graph_data['ITSGraph']= (react_graph, prod_graph, its)
                    graph_data['GraphRules'] = GraphRuleExtraction.extract_reaction_rules(react_graph, prod_graph, its, extend=False, n_knn=1)
                    return graph_data
                else:
                    graph_data['ITSGraph'] = None
                    graph_data['GraphRules'] = None
        else:
            if return_all:
                graph_data['ITSGraph']= (react_graph, prod_graph, its)
                graph_data['GraphRules'] = GraphRuleExtraction.extract_reaction_rules(react_graph, prod_graph, its, extend=False, n_knn=1)
                return graph_data
            else:
                graph_data['ITSGraph'] = None
                graph_data['GraphRules'] = None

        return graph_data

    @staticmethod
    def process_graph_data_parallel(graph_data_list: List[Dict], column: str, n_jobs: int, verbose: int,
                                    return_all: bool=False, ignore_aromaticity: bool = False) -> List[Dict]:
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
            delayed(ITSHAdjuster.process_single_graph_data)(graph_data, column, return_all, ignore_aromaticity) for graph_data in graph_data_list
        )

        return processed_data
    
    @staticmethod
    def add_hydrogen_nodes_multiple_utils(graph: nx.Graph, node_id_pairs: Iterable[Tuple[int, int]], atom_map_update: bool = True) -> nx.Graph:
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
            atom_map_val = new_hydrogen_node_id if atom_map_update else new_graph.nodes[node_id].get('atom_map', 0)
            new_graph.add_node(new_hydrogen_node_id, charge=0, hcount=0, aromatic=False, element='H', atom_map=atom_map_val,
                            isomer='N', partial_charge=0, hybridization=0, in_ring=False,
                            explicit_valence=0, implicit_hcount=0)
            new_graph.add_edge(node_id, new_hydrogen_node_id, order=1.0, ez_isomer='N', bond_type='SINGLE',
                            conjugated=False, in_ring=False)
        return new_graph

    @staticmethod
    def add_hydrogen_nodes_multiple(react_graph: nx.Graph, prod_graph: nx.Graph) -> List[Tuple[nx.Graph, nx.Graph]]:
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

        hydrogen_nodes_form, hydrogen_nodes_break = [], []

        for node_id in react_graph_copy.nodes:
            hcount_diff = react_graph_copy.nodes[node_id].get('hcount', 0) - prod_graph_copy.nodes[node_id].get('hcount', 0)
            if hcount_diff > 0:
                hydrogen_nodes_break.extend([node_id] * hcount_diff)
            elif hcount_diff < 0:
                hydrogen_nodes_form.extend([node_id] * -hcount_diff)

        max_index = max(max(react_graph_copy.nodes, default=0), max(prod_graph_copy.nodes, default=0))
        permutations = list(itertools.permutations(range(max_index + 1, max_index + 1 + len(hydrogen_nodes_form))))
        permutations_seed = permutations[0]  # Used for consistent hydrogen nodes in the reactant graph

        updated_graphs = []

        for permutation in permutations:
            current_react_graph, current_prod_graph = deepcopy(react_graph_copy), deepcopy(prod_graph_copy)

            # Consistent hydrogen nodes in the reactant graph
            current_react_graph = ITSHAdjuster.add_hydrogen_nodes_multiple_utils(current_react_graph, zip(hydrogen_nodes_break, permutations_seed), atom_map_update=False)

            # Varied hydrogen nodes in the product graph based on permutation
            current_prod_graph = ITSHAdjuster.add_hydrogen_nodes_multiple_utils(current_prod_graph, zip(hydrogen_nodes_form, permutation))

            updated_graphs.append((current_react_graph, current_prod_graph))

        return updated_graphs
