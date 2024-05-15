import pathlib
import sys

root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))
import time
import logging
import pandas as pd
from SynTemp.SynAAM.run_consensus import run_consensus_aam
from SynTemp.SynUtils.utils import load_database


# Configure logging
logging.basicConfig(
    filename="./Data/AAM/aam_processing.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

list_data = {
    "ecoli": 273,
    "recon3d": 382,
    "uspto_3k": 3000,
    "golden": 1758,
    "natcomm": 491,
}
mapper_types = ["rxn_mapper", "graphormer", "local_mapper", "rdt"]

rdt_jar_path = f"{root_dir}/Data/RDT_2.4.1.jar"

# Create a dictionary to store time taken for each mapper type for each data
time_dict = {mapper_type: {} for mapper_type in mapper_types}

# Record number of data columns
num_data_columns = pd.Series(list_data)

# Logging number of data columns
logging.info("Number of data columns:\n" + num_data_columns.to_string())

for mapper_type in mapper_types:
    for folder_name, num_columns in list_data.items():
        start_time = time.time()  
        save_dir = f"{root_dir}/Data/{folder_name}"
        data = load_database(
            f"{root_dir}/Data/AAM/{folder_name}/{folder_name}_reactions.json.gz"
        )
        mapped_reactions = run_consensus_aam(
            data=data,
            rsmi_column="reactions",
            save_dir=None,
            mapper_types=[mapper_type],
            data_name=folder_name,
            batch_size=50,
            check_balance=False,
            verbose=0,
            rdt_jar_path=rdt_jar_path,
            working_dir="./",
        )
        end_time = time.time()  # End time for individual data processing
        elapsed_time = round(
            end_time - start_time, 2
        )  # Round elapsed time to two decimal places
        time_dict[mapper_type][folder_name] = elapsed_time  # Update time dictionary
        logging.info(
            f"{folder_name} ({num_columns}) - {mapper_type}: {elapsed_time} seconds"
        )

# Create DataFrame from time dictionary
df = pd.DataFrame(time_dict)

# Log the DataFrame to file
logging.info("\n" + df.to_string())
