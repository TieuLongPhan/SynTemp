import networkx as nx
from syntemp.SynITS.its_extraction import ITSExtraction
from syntemp.SynITS.its_hadjuster import ITSHAdjuster
from syntemp.SynRule.rules_extraction import RuleExtraction
from synutility.SynAAM.its_construction import ITSConstruction
from synutility.SynIO.Format.smi_to_graph import rsmi_to_graph
from synutility.SynIO.debug import setup_logging
from typing import Dict, List, Tuple
from syntemp.SynRule.rc_cluster import RCCluster
from synutility.SynGraph.Descriptor.graph_signature import GraphSignature
from syntemp.SynITS.hydrogen_utils import (
    check_hcount_change,
)
from joblib import Parallel, delayed

logger = setup_logging()


class ITSArbitrary:
    def __init__(self):
        pass

    @staticmethod
    def process_equivalent_map(
        react_graph: nx.Graph,
        prod_graph: nx.Graph,
        ignore_aromaticity: bool,
        balance_its: bool,
    ) -> Tuple[List[nx.Graph], List[nx.Graph]]:
        """
        Process equivalent maps by adding hydrogen nodes and constructing ITS graphs.

        Parameters:
        - react_graph (nx.Graph): The reactant graph.
        - prod_graph (nx.Graph): The product graph.
        - ignore_aromaticity (bool): Whether to ignore aromaticity in graph construction.
        - balance_its (bool): Whether to balance the ITS graph.

        Returns:
        - Tuple of (List[nx.Graph], List[nx.Graph]): Lists of reaction graphs and
        ITS graphs.
        """
        hcount_change = check_hcount_change(react_graph, prod_graph)
        if hcount_change == 0:
            its_list = [ITSConstruction().ITSGraph(react_graph, prod_graph)]
            rc_list = [
                RuleExtraction.extract_reaction_rules(
                    react_graph, prod_graph, i, False, 1
                )[2]
                for i in its_list
            ]
            return list(rc_list), list(its_list)

        combinations_solution = ITSHAdjuster.add_hydrogen_nodes_multiple(
            react_graph, prod_graph
        )

        # Create ITS graphs for each combination solution
        its_list = [
            ITSConstruction.ITSGraph(
                i[0], i[1], ignore_aromaticity, balance_its=balance_its
            )
            for i in combinations_solution
        ]

        # Extract reaction rules for each ITS graph
        rc_list = [
            RuleExtraction.extract_reaction_rules(react_graph, prod_graph, i, False, 1)[
                2
            ]
            for i in its_list
        ]

        # Filter valid reaction graphs and ITS graphs
        valid_rc_its = [
            (rc, its)
            for rc, its in zip(rc_list, its_list)
            if rc is not None and isinstance(rc, nx.Graph) and rc.number_of_nodes() > 0
        ]

        # Unzip valid results
        rc_list, its_list = zip(*valid_rc_its) if valid_rc_its else ([], [])
        return list(rc_list), list(its_list)

    @staticmethod
    def process_non_equivalent_map(
        data: Dict[str, str],
        mapped_key: List[str],
        sanitize: bool,
        ignore_aromaticity: bool,
        balance_its: bool,
    ) -> Tuple[List[nx.Graph], List[nx.Graph]]:
        """
        Process non-equivalent maps and construct their corresponding reaction and
        ITS graphs.

        Parameters:
        - data (Dict[str, str]): Dictionary of mapped SMILES strings.
        - mapped_key (List[str]): List of mapper names to process.
        - sanitize (bool): Whether to sanitize the molecule(s).
        - ignore_aromaticity (bool): Whether to ignore aromaticity in graph construction.
        - balance_its (bool): Whether to balance the ITS graph.

        Returns:
        - Tuple of (List[nx.Graph], List[nx.Graph]): Lists of reaction graphs and
        ITS graphs.
        """
        rc_list, its_list = [], []
        for mapper in mapped_key:
            try:
                # Convert SMILES to graphs
                G, H = rsmi_to_graph(
                    data[mapper],
                    drop_non_aam=True,
                    light_weight=True,
                    sanitize=sanitize,
                )
                # Process equivalent maps
                rc, its = ITSArbitrary.process_equivalent_map(
                    G, H, ignore_aromaticity, balance_its
                )
                # print(rc)
                rc_list.extend(rc)
                its_list.extend(its)
            except Exception as e:
                logger.warning(f"Error processing {mapper}: {e}")
        return rc_list, its_list

    @staticmethod
    def get_unique_graphs_for_clusters(
        graphs: List[nx.Graph], cluster_indices: List[int]
    ) -> List[nx.Graph]:
        """
        Get a unique graph for each cluster from a list of graphs.

        Parameters:
        - graphs (List[nx.Graph]): List of networkx graphs.
        - cluster_indices (List[int]): List of indices that represent cluster assignments
        for each graph.

        Returns:
        - List[nx.Graph]: List of unique graphs, one per cluster.
        """
        # Create a dictionary to store graphs by cluster index
        cluster_graphs = {}

        for idx, cluster_id in enumerate(cluster_indices):
            # Add graph to the appropriate cluster
            if cluster_id not in cluster_graphs:
                cluster_graphs[cluster_id] = []
            cluster_graphs[cluster_id].append(graphs[idx])

        # Now, select one unique graph per cluster (e.g., the first graph)
        unique_graphs = []
        for cluster_id, graphs_in_cluster in cluster_graphs.items():
            unique_graphs.append(graphs_in_cluster[0])

        return unique_graphs

    @staticmethod
    def its_expand(
        data: Dict[str, str],
        mapped_key: List[str],
        check_method: str = "RC",
        id_column: str = "R-id",
        ignore_aromaticity: bool = False,
        confident_mapper: str = "graphormer",
        sanitize: bool = True,
        balance_its: bool = True,
    ) -> Tuple[List[nx.Graph], List[nx.Graph]]:
        """
        Expand ITS graphs by checking equivalence and processing accordingly.

        Parameters:
        - data (Dict[str, str]): Dictionary of mapped SMILES strings.
        - mapped_key (List[str]): List of mapper names to process.
        - check_method (str): Method to check for isomorphism, "RC" or "ITS".
        - id_column (str): Column name for reaction ID.
        - ignore_aromaticity (bool): Whether to ignore aromaticity.
        - confident_mapper (str): Mapper to use when confident.
        - symbol (str): Reaction symbol separator.
        - sanitize (bool): Whether to sanitize molecules.
        - balance_its (bool): Whether to balance ITS graphs.

        Returns:
        - Tuple of (List[nx.Graph], List[nx.Graph]): Lists of reaction graphs
        and ITS graphs.
        """
        try:
            # Process the mapped SMILES strings and check equivalence
            good, _ = ITSExtraction.process_mapped_smiles(
                data,
                mapped_key,
                check_method,
                id_column,
                ignore_aromaticity,
                confident_mapper,
                sanitize,
            )

            # Ensure equivalence check is valid
            if "equivariant" not in good:
                raise ValueError(
                    "Equivalence check result 'equivariant' not found in the response."
                )
            # print(good)
            # Process based on equivalence check
            if good["equivariant"] != (len(mapped_key) - 1):
                # print(1)
                rc_list, its_list = ITSArbitrary.process_non_equivalent_map(
                    data, mapped_key, sanitize, ignore_aromaticity, balance_its
                )
            else:
                r, p = good[confident_mapper][0], good[confident_mapper][1]
                rc_list, its_list = ITSArbitrary.process_equivalent_map(
                    r, p, ignore_aromaticity, balance_its
                )

            sig = [GraphSignature(i).create_graph_signature() for i in rc_list]
            cluster_indices = RCCluster().fit_graphs(rc_list, sig)

            new_rc = ITSArbitrary.get_unique_graphs_for_clusters(
                rc_list, cluster_indices
            )
            new_its = ITSArbitrary.get_unique_graphs_for_clusters(
                its_list, cluster_indices
            )
            return new_rc, new_its

        except Exception as e:
            # Log error and re-raise exception
            logger.error(f"Error in ITSArbitrary.its_expand: {str(e)}")
            return [], []

    def parallel_its_expand(
        self,
        data: List[Dict],
        mapped_key: List[str],
        check_method: str = "RC",
        id_column: str = "R-id",
        ignore_aromaticity: bool = False,
        confident_mapper: str = "graphormer",
        sanitize: bool = True,
        balance_its: bool = True,
        n_jobs: int = 1,
        verbose: int = 0,
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Expands ITS graphs in parallel for a list of reaction data.

        Parameters:
        - data (List[Dict]): List of dictionaries containing mapped reaction data.
        - mapped_key (List[str]): List of keys for processing each reaction data.
        - check_method (str): Method for checking graph equivalence ("RC" or "ITS").
        Default is "RC".
        - id_column (str): The column in the dictionary representing reaction ID.
        - ignore_aromaticity (bool): Whether to ignore aromaticity when
        constructing chemical graphs.
        - confident_mapper (str): Mapper to use when confident.
        Default: "graphormer".
        - sanitize (bool): Whether to sanitize the molecules during graph construction.
        - balance_its (bool): Whether to balance ITS graphs.
        - n_jobs (int): Number of parallel jobs (default: 1).
        - verbose (int): Verbosity level for parallel processing.

        Returns:
        - Tuple of (List[Dict], List[Dict]): Updated list of dictionaries with expanded
        RC and ITS graphs for each reaction.
        """

        logger.info(f"Starting parallel ITS graph expansion with {n_jobs} jobs.")

        try:
            results = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(self.its_expand)(
                    graph_data,
                    mapped_key,
                    check_method,
                    id_column,
                    ignore_aromaticity,
                    confident_mapper,
                    sanitize,
                    balance_its,
                )
                for graph_data in data
            )
        except Exception as e:
            logger.error(f"Error occurred during parallel processing: {e}")
            return [], []  # In case of failure, return empty lists for RC and ITS

        # Process and store the results into the data dictionary
        for key, result in enumerate(results):
            try:
                rc, its = (
                    result  # Deconstruct the tuple (RC graph list, ITS graph list)
                )
                data[key]["RC"] = rc
                data[key]["ITS"] = its
                logger.info(
                    f"Processed reaction {data[key].get(id_column)} successfully."
                )
            except Exception as e:
                logger.warning(f"Error processing reaction at index {key}: {e}")
                data[key]["RC"] = []
                data[key]["ITS"] = []

        # Return the updated data with RC and ITS graphs
        return data
