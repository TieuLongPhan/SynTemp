import pathlib
import sys
import time
import logging
import pandas as pd


root_dir = pathlib.Path(__file__).parents[2]
# Setup logging
logging.basicConfig(
    filename=f"{root_dir}/Data/AAM/unbalance/aam_processing.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

logger = logging.getLogger("AAMProcessing")

sys.path.append(str(root_dir))
from syntemp.SynAAM.atom_map_consensus import AAMConsensus
from syntemp.SynUtils.utils import load_database, save_database

list_data = {
    "ecoli": 273,
    "recon3d": 382,
    "uspto_3k": 3000,
    "golden": 1758,
    "natcomm": 491,
}
mapper_types = ["rxn_mapper", "graphormer", "local_mapper", "rdt"]

rdt_jar_path = f"{root_dir}/Data/RDT_2.4.1.jar"
working_dir = f"{root_dir}"

# Create a dictionary to store time taken for each mapper type for each data
time_dict = {mapper: {} for mapper in mapper_types}


for folder_name, num_columns in list_data.items():
    save_dir = f"{root_dir}/Data/AAM/unbalance/{folder_name}"
    save_path = f"{root_dir}/Data/AAM/unbalance/{folder_name}/{folder_name}_aam_reactions.json.gz"
    data = load_database(
        f"{root_dir}/Data/AAM/unbalance/{folder_name}/{folder_name}_reactions.json.gz"
    )
    logger.info(f"Loaded {len(data)} reactions from {folder_name}")
    for mapper in mapper_types:
        aam = AAMConsensus(data, mappers=[mapper])
        start_time = time.time()
        logger.info(f"Starting consensus mapping for {folder_name} using {mapper}")

        results = aam.batch_consensus(
            data,
            rsmi_column="reactions",
            batch_size=len(data),
            job_timeout=None,
            safe_mode=False,
            rdt_jar_path=rdt_jar_path,
            working_dir=working_dir,
        )
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(
            f"Completed consensus mapping for {folder_name} using {mapper} in {elapsed_time:.2f} seconds"
        )

        # Store the time taken in the dictionary
        time_dict[mapper][folder_name] = elapsed_time
        for key, value in enumerate(data):
            data[key][mapper] = results[key][mapper]
    save_database(data, save_path)
    logger.info(f"Saved AAM {folder_name} results to {save_path}")


# Create a DataFrame from the time dictionary
df = pd.DataFrame(time_dict)

# Log the DataFrame to file
logger.info("\n" + df.to_string())

# Optionally, save the DataFrame to a file for further analysis
df.to_csv(f"{root_dir}/Data/AAM/results_benchmark/mapping_times.csv")
