import pathlib
import logging
import time
import sys

# Setup the root directory based on the script's location
root_dir = pathlib.Path(__file__).resolve().parents[2]
sys.path.append(str(root_dir))
# Configure logging
logging.basicConfig(
    filename=f"{root_dir}/Docs/Analysis/cluster_time.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

# Importing necessary functions and classes
from syntemp.SynUtils.utils import load_from_pickle
from syntemp.SynRule.rule_cluster import RuleCluster
from syntemp.SynRule.rules_extraction import RuleExtraction

# Load data
data = load_from_pickle(f"{root_dir}/Data/Temp/_its_correct.pkl.gz")
its_graphs = [value["ITSGraph"] for value in data]
cluster = RuleCluster()

# Process the data for different values of k and log the processing time
for k in range(4):
    start_time = time.time()  # Start time measurement

    logging.info(f"Processing templates with k={k}")

    if k > 0:
        # Extract reaction rules with extension and k-nearest neighbors if k > 0
        rc_graphs = [
            RuleExtraction.extract_reaction_rules(*value, extend=True, n_knn=k)
            for value in its_graphs
        ]
    else:
        # Extract reaction rules without extension if k = 0
        rc_graphs = [
            RuleExtraction.extract_reaction_rules(*value, extend=False)
            for value in its_graphs
        ]

    # Fit the rule clusters with the extracted graphs
    cluster_indices, templates = cluster.fit(
        rc_graphs, templates=None, update_template=True
    )

    end_time = time.time()  # End time measurement
    processing_time = end_time - start_time  # Calculate processing time

    # Log the processing time
    logging.info(f"Finished processing for k={k} in {processing_time:.2f} seconds")
