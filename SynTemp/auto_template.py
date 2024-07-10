import os
import logging
from typing import List, Any, Optional
from SynTemp.SynRule.hierarchical_clustering import HierarchicalClustering
from SynTemp.SynRule.rule_writing import RuleWriting
from SynTemp.SynUtils.utils import save_to_pickle

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def AutoTemp(data, 
             node_label_names: List[str] = ["element", "charge"],
             node_label_default: List[Any] = ["*", 0], 
             edge_attribute: str = "order", 
             max_radius: int = 3,
             save_path: Optional[str] = None) -> List:
    """
    Automates the generation of templates and rules from hierarchical clustering of given 
    data.

    Parameters:
    - data: Input data to process.
    - node_label_names (List[str]): Labels for the nodes in the hierarchical clustering.
    - node_label_default (List[Any]): Default values for node labels if not specified.
    - edge_attribute (str): The attribute used for defining edges.
    - max_radius (int): Maximum radius for the clustering algorithm.
    - save_path (Optional[str]): Path to save output files.

    Returns:
    - List: A list of rules extracted from the generated templates.
    """
    try:
        hier_cluster = HierarchicalClustering(node_label_names, node_label_default, 
                                              edge_attribute, max_radius)
        logging.info("Hierarchical clustering initialized successfully.")
        reaction_dicts, templates, hier_templates = hier_cluster.fit(data)
        logging.info("Clustering completed and data extracted.")
        
        rules = []
        for radius, template in enumerate(templates):
            directory_path = None
            if save_path:
                directory_path = os.path.join(save_path, f'R{radius}')
                if not os.path.exists(directory_path):
                    os.makedirs(directory_path)
                    logging.info(f"Directory created at {directory_path}")
                else:
                    logging.info(f"Directory {directory_path} already exists")

            write = RuleWriting.auto_extraction(template, id_column='Cluster_id', 
                                                rule_column='RC', reindex=True, 
                                                save_path=directory_path)
            rules.append(write)
            logging.info(f"Rules extracted for template at radius {radius}")

        if save_path:
            save_to_pickle(reaction_dicts, f"{save_path}/data_cluster.pkl.gz")
            save_to_pickle(templates, f"{save_path}/templates.pkl.gz")
            save_to_pickle(hier_templates, f"{save_path}/hier_templates.pkl.gz")
            logging.info("All data successfully saved.")

        return rules
    except Exception as e:
        logging.error("An error occurred during template generation: %s", e)
        raise

