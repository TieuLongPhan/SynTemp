import networkx as nx
from typing import Dict, List, Tuple
from copy import deepcopy
from rdkit import Chem
from joblib import Parallel, delayed
from operator import eq
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match
from syntemp.SynITS.its_construction import ITSConstruction
from syntemp.SynChemistry.mol_to_graph import MolToGraph
from syntemp.SynRule.rules_extraction import RuleExtraction
from syntemp.SynUtils.chemutils import remove_atom_mapping
from syntemp.SynUtils.utils import setup_logging

logger = setup_logging()


class ITSExtraction:
    def __init__(self):
        pass

    @staticmethod
    def graph_from_smiles(smiles: str) -> nx.Graph:
        """
        Constructs a graph representation from a SMILES string.

        Parameters:
        - smiles (str): A SMILES string representing a molecule or a set of molecules.

        Returns:
        - nx.Graph: A graph representation of the molecule(s).
        """

        mol = Chem.MolFromSmiles(smiles)
        graph = MolToGraph().mol_to_graph(mol, drop_non_aam=True)
        return graph

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
    def process_mapped_smiles(
        mapped_smiles: Dict[str, str],
        mapper_names: List[str],
        check_method="RC",  # or ITS
        id_column: str = "R-id",
        ignore_aromaticity: bool = False,
        confident_mapper: str = "graphormer",
        symbol: str = ">>",
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
        rules_by_map = {id_column: mapped_smiles.get(id_column, "N/A")}
        its_graphs = []
        rules_graphs = []

        for mapper in mapper_names:
            try:
                reactants_side, products_side = mapped_smiles[mapper].split(symbol)

                # Get reactants graph G
                G = ITSExtraction.graph_from_smiles(reactants_side)

                # Get products graph H
                H = ITSExtraction.graph_from_smiles(products_side)

                # Construct the ITS graph
                ITS = ITSConstruction.ITSGraph(G, H, ignore_aromaticity)
                its_graphs.append(ITS)

                # Store graphs and ITS
                graphs_by_map[mapper] = (G, H, ITS)

                # Extract reaction rules
                rules_by_map[mapper] = RuleExtraction.extract_reaction_rules(
                    G, H, ITS, extend=False
                )
                _, _, rules = rules_by_map[mapper]
                rules_graphs.append(rules)

            except Exception as e:
                logger.info(f"Error processing {mapper}: {e}")

                # Fallback: Create a one-node graph for ITS and Rules
                one_node_graph = nx.Graph()
                one_node_graph.add_node(0)  # Create a graph with a single node

                # Use the one-node graph for ITS and Rules
                its_graphs.append(one_node_graph)
                graphs_by_map[mapper] = (one_node_graph, one_node_graph, one_node_graph)
                rules_by_map[mapper] = (one_node_graph, one_node_graph, one_node_graph)
                rules_graphs.append(one_node_graph)
        if len(rules_graphs) > 1:
            if check_method == "RC":
                _, equivariant = ITSExtraction.check_equivariant_graph(rules_graphs)
            elif check_method == "ITS":
                _, equivariant = ITSExtraction.check_equivariant_graph(its_graphs)
        else:
            equivariant = 0
        # graphs_by_map['check_equivariant'] = classified
        graphs_by_map["equivariant"] = equivariant

        graphs_by_map_correct = deepcopy(graphs_by_map)
        graphs_by_map_incorrect = deepcopy(graphs_by_map)
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
            if "[O]" in remove_atom_mapping(mapped_smiles[mapper_names[0]]):
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
        n_jobs: int = -1,
        verbose: int = 10,
        id_column: str = "R-id",
        check_method="RC",
        export_full=False,
        ignore_aromaticity: bool = False,
        confident_mapper: str = "graphormer",
        symbol: str = ">>",
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
                symbol,
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
