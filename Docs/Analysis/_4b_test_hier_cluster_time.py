import pathlib
import logging
import sys

# Setup the root directory based on the script's location
root_dir = pathlib.Path(__file__).resolve().parents[2]
sys.path.append(str(root_dir))
# Configure logging
logging.basicConfig(
    filename=f"{root_dir}/Docs/Analysis/hier_cluster_time.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# Importing necessary functions and classes
from syntemp.SynUtils.utils import load_from_pickle

from syntemp.SynRule.hierarchical_clustering import HierarchicalClustering
from syntemp.SynRule.rules_extraction import RuleExtraction

# Load data
data = load_from_pickle(f"{root_dir}/Data/Temp/_its_correct.pkl.gz")
its_graphs = [value["ITSGraph"] for value in data]
cluster = HierarchicalClustering()
logging.info(f"Processing templates")
cluster.fit(data)
