import pathlib
import time
import logging
import os
import sys
import shutil
from joblib import Parallel, delayed
import tempfile

root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))
import pandas as pd
from SynTemp.SynUtils.utils import (
    # load_database,
    # save_database,
    save_to_pickle,
    load_from_pickle,
)
from SynTemp.SynRule.rule_cluster import RuleCluster


def configure_logging(save_dir: str, verbose: int, data_name: str) -> logging.Logger:
    """
    Configures logging to a file and console with configurable verbosity.

    Args:
        save_dir (str): Directory to save log file. If None, only log to console.
        verbose (int): Verbosity level. Higher value means more logging.
        data_name (str): Name of the data, used for naming the log file.

    Returns:
        logging.Logger: Configured logger.
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)

    # Clear existing handlers
    logger.handlers = []

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler
    if save_dir:
        log_file = os.path.join(save_dir, f"{data_name}_alignment.log")
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def rule_align(data_name, save_dir, verbose=0, n_jobs=4):
    logger = configure_logging(save_dir, verbose, data_name)
    logger.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
    logger.info(f"Processing data {data_name} with {n_jobs} cpus")
    process_graph_data = load_from_pickle(
        f"{root_dir}/Data/{data_name}/{data_name}_its_correct.pkl.gz"
    )
    logger.info(f"Data all for alignment: {len(process_graph_data)}")
    logger.info("Clustering ITS graphs.")
    node_label_names = ["element", "charge"]
    naive_cluster = RuleCluster(
        node_label_names=node_label_names,
        node_label_default=["*", 0],
        edge_attribute="order",
    )
    its_graph_rules_cluster = naive_cluster.process_rules_clustering(
        process_graph_data, rule_column="GraphRules"
    )
    if save_dir:
        save_to_pickle(
            its_graph_rules_cluster,
            f"{save_dir}/{data_name}_its_graph_rules_cluster.pkl.gz",
        )
    logger.info(
        f"Number of clusters: {pd.DataFrame(its_graph_rules_cluster)['naive_cluster'].value_counts()}"
    )
    return its_graph_rules_cluster


if __name__ == "__main__":
    folder_name = "USPTO_50K"
    save_dir = f"{root_dir}/Data/{folder_name}"
    rule_align(folder_name, save_dir)
