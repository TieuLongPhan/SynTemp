import networkx as nx
import networkx as nx
from networkx.algorithms.cycles import simple_cycles
from networkx.algorithms.planarity import check_planarity

from joblib import Parallel, delayed
from typing import Dict, Tuple, List, Optional
from SynTemp.SynITS.its_hadjuster import ITSHAdjuster

class UncertainRefinement:
    @staticmethod
    def process_and_check_graph(graph_data: Dict, graph_type: str) -> Tuple[Dict, str]:
        """
        Process the graph data and check its type.

        Parameters:
        - graph_data (Dict): The graph data to be processed.
        - graph_type (str): The type of the graph.

        Returns:
        - Tuple[Dict, str]: Processed graph data and its type.
        """
        processed_data = ITSHAdjuster.process_single_graph_data(graph_data, graph_type, return_all=True)
        graph_type = UncertainRefinement.check_graph_type(processed_data['GraphRules'][2])
        return processed_data, graph_type

    @staticmethod
    def process_dict(input_graph: Dict, graph_type: str = 'ITSGraph', mapper_types: List[str] = ['rxn_mapper', 'graphormer', 'local_mapper']) -> Optional[Dict]:
        """
        Process a graph dictionary to find a 'Single Cyclic' graph.

        Parameters:
        - input_graph (Dict): The graph dictionary to be processed.

        Returns:
        - Dict: The processed graph data with 'Single Cyclic' type or None if not found.
        """
        uncertain_graph = []
        results = []

        for mapper in mapper_types:
            graph_data = {graph_type: input_graph.get(mapper)}
            if graph_data[graph_type]:
                try:
                    processed_data, graph_type = UncertainRefinement.process_and_check_graph(graph_data, graph_type)
                    uncertain_graph.append(processed_data)
                    results.append(graph_type)
                    if graph_type == 'Single Cyclic':
                        processed_data['R-id'] = input_graph.get('R-id', '')
                        return processed_data
                except:
                    return None
        return None

    @staticmethod
    def process_graphs_in_parallel(graph_dicts: List[Dict], n_jobs: int = 4, verbose: int = 1) -> List[Optional[Dict]]:
        """
        Process multiple graph dictionaries in parallel to find 'Single Cyclic' graphs.

        Parameters:
        - graph_dicts (List[Dict]): A list of graph dictionaries to be processed.
        - n_jobs (int): The number of concurrent jobs.
        - verbose (int): The verbosity level.

        Returns:
        - List[Optional[Dict]]: A list of processed graph data.
        """
        results = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(UncertainRefinement.process_dict)(graph_dict) for graph_dict in graph_dicts)
        return results

    @staticmethod
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
  
    @staticmethod
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

    @staticmethod
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

    @staticmethod
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

        if UncertainRefinement.is_acyclic_graph(G):
            return "Acyclic"
        elif UncertainRefinement.is_single_cyclic_graph(G):
            return "Single Cyclic"
        elif UncertainRefinement.is_complex_cyclic_graph(G):
            return "Complex Cyclic"
        else:
            return "None"
