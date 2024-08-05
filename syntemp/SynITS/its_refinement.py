from joblib import Parallel, delayed
from typing import Dict, Tuple, List, Optional
from syntemp.SynITS.its_hadjuster import ITSHAdjuster
from syntemp.SynRule.rules_extraction import RuleExtraction
from syntemp.SynUtils.graph_utils import check_graph_type


class ITSRefinement:
    @staticmethod
    def process_and_check_graph(graph_data: Dict, column: str) -> Tuple[Dict, str]:
        """
        Process the graph data and check its type.

        Parameters:
        - graph_data (Dict): The graph data to be processed.
        - column (str): The column name of the graph data.

        Returns:
        - Tuple[Dict, str]: Processed graph data and its type.
        """
        processed_data = ITSHAdjuster.process_single_graph_data(
            graph_data, column, return_all=True
        )
        if len(processed_data["GraphRules"][2].nodes) > 0:
            type_of_graph = check_graph_type(processed_data["GraphRules"][2])
        else:
            type_of_graph = "None"
        return processed_data, type_of_graph

    @staticmethod
    def process_dict(
        input_graph: Dict,
        column: str = "ITSGraph",
        mapper_types: List[str] = ["rxn_mapper", "graphormer", "local_mapper"],
    ) -> Optional[Dict]:
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
            graph_data = {column: input_graph.get(mapper)}

            if graph_data[column]:
                try:
                    _, _, rc = RuleExtraction.extract_reaction_rules(
                        graph_data[column][0],
                        graph_data[column][1],
                        graph_data[column][2],
                        extend=False,
                        n_knn=0,
                    )
                    if len(rc.edges) > 20:
                        return None
                    processed_data, type_of_graph = (
                        ITSRefinement.process_and_check_graph(graph_data, column)
                    )
                    uncertain_graph.append(processed_data)
                    results.append(type_of_graph)
                    if type_of_graph == "Single Cyclic":
                        processed_data["R-id"] = input_graph.get("R-id", "")
                        return processed_data
                except ValueError:
                    return None
        return None

    @staticmethod
    def process_graphs_in_parallel(
        graph_dicts: List[Dict],
        mapper_types: List[str] = ["rxn_mapper", "graphormer", "local_mapper"],
        n_jobs: int = 4,
        verbose: int = 1,
    ) -> List[Optional[Dict]]:
        """
        Process multiple graph dictionaries in parallel to find 'Single Cyclic' graphs.

        Parameters:
        - graph_dicts (List[Dict]): A list of graph dictionaries to be processed.
        - n_jobs (int): The number of concurrent jobs.
        - verbose (int): The verbosity level.

        Returns:
        - List[Optional[Dict]]: A list of processed graph data.
        """
        results = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(ITSRefinement.process_dict)(graph_dict, "ITSGraph", mapper_types)
            for graph_dict in graph_dicts
        )
        return results
