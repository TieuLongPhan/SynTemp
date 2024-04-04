import networkx as nx
import pandas as pd
from typing import Dict, List, Tuple, Union, Any
from rdkit import Chem
from operator import eq
from joblib import Parallel, delayed
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match
from SynTemp.SynITS.its_construction import ITSConstruction
from SynTemp.SynProcessor.mol_to_graph import MolToGraph
from SynTemp.SynITS.graph_rules_extraction import GraphRuleExtraction
from itertools import combinations
class AMMValidator:
    def __init__(self):
        """Initializes the AMMValidator class."""
        pass

    @staticmethod
    def graph_from_smiles(smiles: str) -> nx.Graph:
        """
        Constructs a graph representation from a SMILES string.

        Parameters:
            smiles (str): A SMILES string representing a molecule or a set of molecules.

        Returns:
            nx.Graph: A graph representation of the molecule(s).
        """
        mol = Chem.MolFromSmiles(smiles)
        graph = MolToGraph().mol_to_graph(mol, drop_non_aam=True)
        return graph

    @staticmethod
    def check_equivariant_graph(its_graphs: List[nx.Graph]) -> Tuple[List[Tuple[int, int]], int]:
        """
        Checks for isomorphism among a list of ITS graphs and identifies all pairs of isomorphic graphs.

        Parameters:
        - its_graphs (List[nx.Graph]): A list of ITS graphs.

        Returns:
        - List[Tuple[int, int]]: A list of tuples representing pairs of indices of isomorphic graphs.
        - int: The count of unique isomorphic graph pairs found.
        """
        nodeLabelNames = ["typesGH"]
        nodeLabelDefault = ["*", False, 0, 0, ()]
        nodeLabelOperator = [eq, eq, eq, eq, eq]
        nodeMatch = generic_node_match(nodeLabelNames, nodeLabelDefault, nodeLabelOperator)
        edgeMatch = generic_edge_match("order", 1, eq)

        classified = []

        # Use combinations to check each unique pair of graphs without repetition
        for i, j in combinations(range(len(its_graphs)), 2):
            if nx.is_isomorphic(its_graphs[i], its_graphs[j], node_match=nodeMatch, edge_match=edgeMatch):
                classified.append((i, j))

        return classified, len(classified)

    @staticmethod
    def smiles_check(mapped_smile: str, 
                     ground_truth: str, 
                     check_method: str = 'RC', # or 'ITS'
                     ignore_aromaticity: bool = False) -> bool:
        """
        Checks the equivalence of mapped SMILES against ground truth using reaction center (RC) or ITS graph method.

        Parameters:
            mapped_smile (str): The mapped SMILES string.
            ground_truth (str): The ground truth SMILES string.
            check_method (str): The method used for validation ('RC' or 'ITS').
            ignore_aromaticity (bool): Flag to ignore aromaticity in ITS graph construction.

        Returns:
            bool: True if the mapped SMILES is equivalent to the ground truth, False otherwise.
        """
        its_graphs = []
        rules_graphs = []

        for rsmi in [mapped_smile, ground_truth]:
            reactants_side, products_side = rsmi.split(">>")
            G = AMMValidator.graph_from_smiles(reactants_side)  # Reactants graph
            H = AMMValidator.graph_from_smiles(products_side)  # Products graph

            ITS = ITSConstruction.ITSGraph(G, H, ignore_aromaticity)
            its_graphs.append(ITS)
              
            rules = GraphRuleExtraction.extract_reaction_rules(G, H, ITS, extend=False)
            rules_graphs.append(rules[2])
 
        _, equivariant = AMMValidator.check_equivariant_graph(rules_graphs if check_method == 'RC' else its_graphs)
        
        return equivariant == 1

    @staticmethod
    def validate_smiles(data: Union[pd.DataFrame, List[Dict[str, str]]],
                        ground_truth_col: str = 'ground_truth',
                        mapped_cols: List[str] = ['rxn_mapper', 'graphormer', 'local_mapper'],
                        check_method: str = 'RC',
                        ignore_aromaticity: bool = False,
                        n_jobs: int = 1,
                        verbose: int = 0) -> List[Dict[str, Union[str, float, List[bool]]]]:
        """
        Validates collections of mapped SMILES against their ground truths for multiple mappers and calculates the accuracy.

        Parameters:
            data (Union[pd.DataFrame, List[Dict[str, str]]]): The input data containing mapped and ground truth SMILES.
            ground_truth_col (str): The name of the column or key containing the ground truth SMILES.
            mapped_cols (List[str]): The list of columns or keys containing the mapped SMILES for different mappers.
            check_method (str): The method used for validation ('RC' or 'ITS').
            ignore_aromaticity (bool): Whether to ignore aromaticity in ITS graph construction.
            n_jobs (int): The number of parallel jobs to run. -1 means using all processors.
            verbose (int): The verbosity level for joblib's parallel execution.

        Returns:
            List[Dict[str, Union[str, float, List[bool]]]]: A list of dictionaries, each containing the mapper name,
            accuracy, and individual results for each SMILES pair.
        """
        validation_results = []

        for mapped_col in mapped_cols:
            def check_pair(mapping: Dict[str, str]) -> bool:
                return AMMValidator.smiles_check(
                    mapping[mapped_col], mapping[ground_truth_col], check_method, ignore_aromaticity
                )

            if isinstance(data, pd.DataFrame):
                mappings = data.to_dict('records')
            elif isinstance(data, list):
                mappings = data
            else:
                raise ValueError("Data must be either a pandas DataFrame or a list of dictionaries.")

            # Use joblib to parallelize the validation checks
            results = Parallel(n_jobs=n_jobs, verbose=verbose)(delayed(check_pair)(mapping) for mapping in mappings)
            accuracy = sum(results) / len(mappings) if mappings else 0

            # Store the results for each mapper in the list
            validation_results.append({
                'mapper': mapped_col,
                'accuracy': accuracy,
                'results': results
            })

        return validation_results

