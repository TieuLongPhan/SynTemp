from joblib import Parallel, delayed
from typing import Dict, Tuple, List, Optional
from SynTemp.SynITS.its_hadjuster import ITSHAdjuster
from SynTemp.SynUtils.graph_utils import check_graph_type

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
        graph_type = check_graph_type(processed_data['GraphRules'][2])
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