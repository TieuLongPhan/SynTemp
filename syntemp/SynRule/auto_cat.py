from operator import eq
import networkx as nx
from networkx.algorithms.isomorphism import generic_node_match, generic_edge_match
from typing import List, Dict, Any, Callable, Optional, Tuple

from synkit.Graph.Feature.graph_descriptors import GraphDescriptor


class AutoCat:
    def __init__(
        self,
        templates: List[Dict],
        node_label_names: List[str] = ["element", "charge"],
        node_label_default: List[Any] = ["*", 0],
        edge_attribute: str = "order",
    ):
        """
        Initializes an AutoCat instance which uses graph templates for categorizing
        new graphs.

        Parameters:
        - templates (List[Dict]): A list of template dictionaries containing
        graph information and classification.
        - node_label_names (List[str]): Names of the node attributes to use in
        isomorphism checks.
        - node_label_default (List[Any]): Default values for node attributes if they are
        missing in the graph data.
        - edge_attribute (str): The edge attribute to consider when checking isomorphism
        between graphs.

        Raises:
        - ValueError: If the lengths of `node_label_names` and `node_label_default`
        do not match.
        """
        if len(node_label_names) != len(node_label_default):
            raise ValueError(
                "The lengths of `node_label_names` and `node_label_default` must match."
            )

        self.templates = templates
        self.nodeLabelNames = node_label_names
        self.nodeLabelDefault = node_label_default
        self.edgeAttribute = edge_attribute
        self.nodeMatch = generic_node_match(
            self.nodeLabelNames, self.nodeLabelDefault, [eq] * len(node_label_names)
        )
        self.edgeMatch = generic_edge_match(edge_attribute, 1, eq)

    def lib_check(
        self,
        data: Dict,
        RC_col: str = "RC",
        nodeMatch: Optional[Callable] = None,
        edgeMatch: Optional[Callable] = None,
        attribute: str = "cycle",
    ) -> Dict:
        """
        Checks and classifies a graph based on existing templates. If no match is found,
        a new class is assigned.

        Parameters:
        - data (Dict): A dictionary representing a graph with its attributes and
        classification.
        - nodeMatch (Optional[Callable]): A function to match nodes.
        Defaults to a predefined generic_node_match.
        - edgeMatch (Optional[Callable]): A function to match edges.
        Defaults to a predefined generic_edge_match.

        Returns:
        - Dict: The updated dictionary with its classification.
        """
        cycle = data[attribute]
        sub_temp = [temp for temp in self.templates if temp[attribute] == cycle]

        for template in sub_temp:
            if nx.is_isomorphic(
                template[RC_col],
                data[RC_col],
                node_match=nodeMatch or self.nodeMatch,
                edge_match=edgeMatch or self.edgeMatch,
            ):
                data["class"] = template["class"]
                break
        else:
            new_class = max((temp["class"] for temp in self.templates), default=-1) + 1
            data["class"] = new_class
            self.templates.append(
                data.copy()
            )  # Make sure to append a copy to avoid reference issues

        return data

    def fit(
        self, data: List[Dict], RC_col: str = "RC", attribute: str = "cycle"
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Processes a list of graph data entries, classifying each based on
        existing templates.

        Parameters:
        - data (List[Dict]): A list of dictionaries, each representing a graph
        to be classified.

        Returns:
        - Tuple[List[Dict], List[Dict]]: A tuple containing the list of classified data
        and the updated templates.
        """

        data = GraphDescriptor.process_entries_in_parallel(data, RC_col)
        for entry in data:
            self.lib_check(entry, RC_col, self.nodeMatch, self.edgeMatch, attribute)
        return data, self.templates
