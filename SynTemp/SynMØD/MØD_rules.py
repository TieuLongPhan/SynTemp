import networkx as nx
from typing import Tuple, Dict, List, Optional
from joblib import Parallel, delayed
import os
class MØDRules:
    @staticmethod
    def convert_graph_to_gml(graph: nx.Graph, section: str) -> str:
        """
        Convert a NetworkX graph to a GML string representation, focusing on nodes for the 'context' section
        and on nodes and edges for the 'left' or 'right' sections.

        Args:
        graph (nx.Graph): The NetworkX graph to be converted.
        section (str): The section name in the GML output, typically "left", "right", or "context".

        Returns:
        str: The GML string representation of the graph for the specified section.
        """
        order_to_label = {1: '-', 1.5: ':', 2: '=', 3: '#'}
        gml_str = f"   {section} [\n"

        if section == 'context':
            for node in graph.nodes(data=True):
                element = node[1].get('element', 'X') 
                gml_str += f'      node [ id {node[0]} label "{element}" ]\n'

        if section != 'context':
            for edge in graph.edges(data=True):
                label = order_to_label.get(edge[2].get('order', 1), '-')
                gml_str += f'      edge [ source {edge[0]} target {edge[1]} label "{label}" ]\n'

        gml_str += "   ]\n"
        return gml_str

    @staticmethod
    def RulesGrammar(L: nx.Graph, R: nx.Graph, K: nx.Graph, rule_name: str) -> str:
        """
        Generate a GML string representation for a chemical rule, including its left, context, and right graphs.

        Args:
        L (nx.Graph): The left graph.
        R (nx.Graph): The right graph.
        K (nx.Graph): The context graph.
        rule_name (str): The name of the rule.

        Returns:
        str: The GML string representation of the rule.
        """
        gml_str = "rule [\n"
        gml_str += f'   ruleID "{rule_name}"\n'
        gml_str += MØDRules.convert_graph_to_gml(L, "left")
        gml_str += MØDRules.convert_graph_to_gml(K, "context")
        gml_str += MØDRules.convert_graph_to_gml(R, "right")
        gml_str += "]"
        return gml_str

    @staticmethod
    def process_graph_rules(graph_rules: Dict[str, Tuple[nx.Graph, nx.Graph, nx.Graph]], id_column : str = 'R-id',
                            rule_column :str = 'GraphRules', reindex: bool = False) -> Dict[str, str]:
        """
        Process a dictionary of graph rules to generate GML strings for each rule, with an option to reindex nodes and edges.

        Args:
        graph_rules (Dict[str, Tuple[nx.Graph, nx.Graph, nx.Graph]]): A dictionary mapping rule names to tuples of (L, R, K) graphs.
        reindex (bool): If true, reindex node IDs based on the L graph sequence.

        Returns:
        Dict[str, str]: A dictionary mapping rule names to their GML string representations.
        """
        
        rule_name = graph_rules[id_column]
        L,K,R = graph_rules[rule_column]
        if reindex:
            # Create an index mapping from L graph
            index_mapping = {old_id: new_id for new_id, old_id in enumerate(L.nodes(), 1)}
            
            # Apply the mapping to L, R, and K graphs
            L = nx.relabel_nodes(L, index_mapping)
            R = nx.relabel_nodes(R, index_mapping)
            K = nx.relabel_nodes(K, index_mapping)
        
        rule_grammar = MØDRules.RulesGrammar(L, R, K, rule_name)
        return rule_grammar
    
    @classmethod
    def auto_extraction(cls, data_dicts: List[Dict], id_column: str = 'R-id',
                        n_jobs: int = 4, verbose: int = 1, rule_column: str = 'GraphRules',
                        reindex: bool = False, save_path: Optional[str] = None) -> List[str]:
        """
        Automatically extract and process a list of graph rules, optionally saving the results to files.

        Args:
            data_dicts (List[Dict]): A list of dictionaries, each representing a rule with its associated graphs.
            id_column (str): The dictionary key where the rule ID is stored.
            n_jobs (int): The number of concurrent jobs to run.
            verbose (int): The verbosity level of the job execution.
            rule_column (str): The key in the dictionary that maps to the tuple of (L, R, K) graphs.
            reindex (bool): If True, reindex node IDs based on the L graph sequence.
            save_path (Optional[str]): The directory path where the GML files will be saved. If None, files are not saved.

        Returns:
            List[str]: A list of GML string representations for each rule.
        """
        results = Parallel(n_jobs=n_jobs, verbose=verbose)(
            delayed(cls.process_graph_rules)(data_dict, id_column=id_column,
                                             rule_column=rule_column, reindex=reindex)
            for data_dict in data_dicts
        )

        if save_path and not os.path.exists(save_path):
            os.makedirs(save_path)

        if save_path:
            for i, gml in enumerate(results):
                rule_id = data_dicts[i].get(id_column, f"rule_{i}")
                file_path = os.path.join(save_path, f"{rule_id}.gml")
                with open(file_path, 'w') as file:
                    file.write(gml)

        return results
    
    @staticmethod
    def check_property_changes_in_nodes(left_graph: nx.Graph, right_graph: nx.Graph, properties: list[str]) -> list:
        """
        Compares specified properties of nodes between two graphs and identifies nodes with property changes.

        Parameters:
        - left_graph (nx.Graph): The first graph to compare, serving as the reference.
        - right_graph (nx.Graph): The second graph to compare, checked against the reference.
        - properties (list[str]): A list of property names (strings) to check for changes on each node.

        Returns:
        - List[int]: A list of node indices where any of the specified properties have changed between the two graphs.
        """
        
        if not isinstance(left_graph, nx.Graph) or not isinstance(right_graph, nx.Graph):
            raise ValueError("Both left_graph and right_graph must be NetworkX graph instances.")

        changed_nodes = []
        for node in left_graph.nodes:
            if node in right_graph.nodes:
                for prop in properties:
                    if left_graph.nodes[node].get(prop) != right_graph.nodes[node].get(prop):
                        changed_nodes.append(node)
                        break  

        return changed_nodes

        