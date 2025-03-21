from copy import copy
import networkx as nx
from operator import eq
from joblib import Parallel, delayed
from typing import Dict, List, Tuple
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match

from synkit.IO.debug import setup_logging
from synkit.Graph.ITS.its_construction import ITSConstruction
from synkit.IO.chem_converter import rsmi_to_graph
from synkit.Chem.Reaction.standardize import Standardize
from synkit.Graph.Feature.graph_signature import GraphSignature


from syntemp.SynRule.rules_extraction import RuleExtraction

logger = setup_logging()


class ITSExtraction:
    def __init__(self):
        pass

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
        # If there's only one graph, no comparison is possible
        if len(its_graphs) == 1:
            return [], 0

        # Define node and edge matchers for graph isomorphism check
        nodeLabelNames = ["typesGH"]
        nodeLabelDefault = [()]
        nodeLabelOperator = [eq]
        nodeMatch = generic_node_match(
            nodeLabelNames, nodeLabelDefault, nodeLabelOperator
        )
        edgeMatch = generic_edge_match("order", 1, eq)

        # List to store classified isomorphic pairs
        classified = []

        # Compare each graph to the first one
        for i in range(1, len(its_graphs)):
            if not nx.is_isomorphic(
                its_graphs[0], its_graphs[i], node_match=nodeMatch, edge_match=edgeMatch
            ):
                return [], -1  # Early exit if no isomorphism is found

            classified.append((0, i))

        return classified, len(classified)

    @staticmethod
    def process_mapped_smiles(
        mapped_smiles: Dict[str, str],
        mapper_names: List[str],
        check_method="RC",  # or ITS
        id_column: str = "R-id",
        ignore_aromaticity: bool = False,
        confident_mapper: str = "graphormer",
        sanitize: bool = True,
    ) -> Dict[str, any]:
        """
        Processes mapped SMILES strings representing chemical reactions by constructing
        graph representations for reactants and products, generating ITS (Imaginary
        Transition State) graphs, and checking for isomorphism among these graphs. The
        function also evaluates whether the number of equivariant graphs equals a
        specified threshold to determine a specific graph output.

        Parameters:
        - mapped_smiles (Dict[str, str]): A dictionary where keys are mapper names and
        values are SMILES strings of reactions.
        - mapper_names (List[str]): A list of mapper names to be processed.
        - check_method (str): A method to check for isomorphism among the ITS graphs.
        Either 'RC' or 'ITS'. Defaults to 'RC'.
        - id_column (str): The name of the column in the dataframe that contains the
        reaction ID. Defaults to 'R-id'.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in the reaction
        graphs. Defaults to False.
        - confident_mapper (str): The name of the mapper that was used to generate the
        reaction graphs. Defaults to 'graphormer'.
        - symbol (str): The symbol used to separate reactants and products in the
        reaction SMILES string. Defaults to '>>'.
        - sanitize (bool): Whether to sanitize the molecule(s).
        - id_column (str): The name of the column in the dataframe that contains the
        reaction ID. Defaults to 'R-id'.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in the reaction
        graphs. Defaults to False.
        - confident_mapper (str): The name of the mapper that was used to generate the
        reaction graphs. Defaults to 'graphormer'.
        - symbol (str): The symbol used to separate reactants and products in the
        reaction SMILES string. Defaults to '>>'.
        - sanitize (bool): Whether to sanitize the molecule(s).

        Returns:
        - Dict[str, any]: A dictionary containing graph representations for each reaction
        (as tuples of reactants graph, products graph,
        and ITS graph), ITS graphs, and isomorphism results. Additionally, it includes
        either the first ITS graph or None under the
        'ISTGraph' key, depending on the equivalence of the number of equivariant graphs
        to the threshold.
        """
        threshold = len(mapper_names) - 1
        graphs_by_map = {id_column: mapped_smiles.get(id_column, "N/A")}
        for mapper in mapper_names:
            graphs_by_map[mapper] = (None, None, None)
        rules_by_map = {id_column: mapped_smiles.get(id_column, "N/A")}
        its_graphs = []
        rules_graphs = []
        sig = []
        for mapper in mapper_names:
            try:
                # reactants_side, products_side = mapped_smiles[mapper].split(symbol)
                G, H = rsmi_to_graph(
                    mapped_smiles[mapper],
                    drop_non_aam=True,
                    light_weight=True,
                    sanitize=sanitize,
                )

                # Construct the ITS graph
                ITS = ITSConstruction.ITSGraph(G, H, ignore_aromaticity)
                its_graphs.append(ITS)
                graph_rules = RuleExtraction.extract_reaction_rules(
                    G, H, ITS, extend=False
                )
                rc = graph_rules[2]
                sig_rc = GraphSignature(rc).create_graph_signature()
                if len(sig) > 0:
                    if sig[-1] != sig_rc:
                        rules_graphs = []
                        break

                # Store graphs and ITS
                graphs_by_map[mapper] = (G, H, ITS)

                # Extract reaction rules
                rules_by_map[mapper] = graph_rules
                _, _, rules = rules_by_map[mapper]
                rules_graphs.append(rules)

            except Exception as e:
                logger.info(f"Error processing {mapper}: {e}")
                rules_graphs = []
                break

        if len(rules_graphs) == len(mapper_names):
            if check_method == "RC":
                _, equivariant = ITSExtraction.check_equivariant_graph(rules_graphs)
            elif check_method == "ITS":
                _, equivariant = ITSExtraction.check_equivariant_graph(its_graphs)
        else:
            equivariant = -1
        graphs_by_map["equivariant"] = equivariant
        graphs_by_map_correct = copy(graphs_by_map)
        graphs_by_map_incorrect = copy(graphs_by_map)
        is_empty_graph_present = any(
            (not isinstance(graph, nx.Graph) or graph.number_of_nodes() == 0)
            for graph in rules_graphs
        )
        # Determine the target dictionary based on conditions
        target_dict = (
            graphs_by_map_incorrect
            if is_empty_graph_present or equivariant < threshold
            else graphs_by_map_correct
        )

        # Check if mapper_names is not empty to avoid IndexError
        if mapper_names:
            if "[O]" in Standardize().remove_atom_mapping(
                mapped_smiles[mapper_names[0]]
            ):
                target_dict["ITSGraph"] = graphs_by_map.get(mapper_names[0], None)
                target_dict["GraphRules"] = rules_by_map.get(mapper_names[0], None)
            else:
                if confident_mapper in mapper_names:
                    target_dict["ITSGraph"] = graphs_by_map.get(confident_mapper, None)
                    target_dict["GraphRules"] = rules_by_map.get(confident_mapper, None)
                else:
                    target_dict["ITSGraph"] = graphs_by_map.get(mapper_names[0], None)
                    target_dict["GraphRules"] = rules_by_map.get(mapper_names[0], None)

        return graphs_by_map_correct, graphs_by_map_incorrect

    @staticmethod
    def parallel_process_smiles(
        mapped_smiles_list: List[Dict[str, str]],
        mapper_names: List[str],
        n_jobs: int = 1,
        verbose: int = 0,
        id_column: str = "R-id",
        check_method="RC",
        export_full=False,
        ignore_aromaticity: bool = False,
        confident_mapper: str = "graphormer",
        sanitize: bool = True,
    ) -> List[Dict[str, any],]:
        """
        Processes a list of mapped SMILES strings in parallel.

        Parameters:
        - mapped_smiles_list (List[Dict[str, str]]): A list where each element is a
        dictionary of mapped SMILES strings.
        - mapper_names (List[str]): A list of mapper names to process.
        - n_jobs (int): The number of jobs to run in parallel. Defaults to -1, which means
        using all processors.
        - verbose (int): The verbosity level of the parallel processing.
        - check_method (str): A method to check for isomorphism among the ITS graphs.
        Either 'RC' or 'ITS'. Defaults to 'RC'.
        - export_full (bool): Whether to export the full results. Defaults to False.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in the graph.
        Defaults to False.
        - confident_mapper (str): The mapper name to use if the check_method is 'RC'.
        Defaults to 'graphormer'.
        - symbol (str): The symbol to use if the check_method is 'RC'. Defaults to '>>'.
        - sanitize (bool): Whether to sanitize the molecule(s). Defaults to True.
        - export_full (bool): Whether to export the full results. Defaults to False.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in the graph.
        Defaults to False.
        - confident_mapper (str): The mapper name to use if the check_method is 'RC'.
        Defaults to 'graphormer'.
        - symbol (str): The symbol to use if the check_method is 'RC'. Defaults to '>>'.
        - sanitize (bool): Whether to sanitize the molecule(s). Defaults to True.

        Returns:
        - List[Dict[str, any]]: A list of dictionaries containing graph representations
        for each reaction, ITS graphs, and isomorphism results for each mapped SMILES
        string.
        """
        results = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(ITSExtraction.process_mapped_smiles)(
                mapped_smiles,
                mapper_names,
                check_method,
                id_column,
                ignore_aromaticity,
                confident_mapper,
                sanitize,
            )
            for mapped_smiles in mapped_smiles_list
        )
        results_correct = [result[0] for result in results]
        results_incorrect = [result[1] for result in results]
        if export_full:
            new_results_correct = [d for d in results_correct if "ITSGraph" in d.keys()]
            new_results_incorrect = [
                d for d in results_incorrect if "ITSGraph" in d.keys()
            ]

        else:
            new_results_correct = [
                {
                    id_column: d[id_column],
                    "ITSGraph": d["ITSGraph"],
                    "GraphRules": d["GraphRules"],
                }
                for d in results_correct
                if "ITSGraph" in d.keys()
            ]
            new_results_incorrect = [
                d for d in results_incorrect if "ITSGraph" in d.keys()
            ]
        return new_results_correct, new_results_incorrect
