import argparse
import pandas as pd
import logging
import traceback
from synrbl import Balancer
from synrbl.rsmi_utils import load_database, save_database

def setup_logging(log_file=None):
    """Set up logging configuration."""
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"
    if log_file:
        logging.basicConfig(
            filename=log_file,
            filemode='w',
            level=logging.INFO,
            format=log_format,
            datefmt=date_format,
        )
    else:
        logging.basicConfig(level=logging.INFO, format=log_format, datefmt=date_format)

def main(input_file, output_file, log_file, n_jobs=4, batch_size=5000):
    setup_logging(log_file)

    synrbl = Balancer(reaction_col="reactions", id_col="id", n_jobs=n_jobs, batch_size=batch_size)

    try:
        logging.info(f"Processing file {input_file}")
        try:
            reactions = load_database(input_file)
        except Exception as e:
            logging.error(f"Failed to load database, attempting to read CSV. Error: {e}")
            reactions = pd.read_csv(input_file).to_dict('records')
        
        balanced_reactions = synrbl.rebalance(
            reactions=reactions, output_dict=False
        )

        for i, new_reaction in enumerate(balanced_reactions):
            reactions[i]["reactions"] = new_reaction

        save_database(
            database=reactions, pathname=output_file
        )
        logging.info("Reactions balanced and saved successfully.")

    except Exception as e:
        logging.error(
            f"An error occurred while processing {input_file}: {e}\n"
            + f"{traceback.format_exc()}"
        )

    logging.info("Processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script processes chemical reaction files, "
        + "rebalances them using SynRBL, and logs the processing progress. "
        + "It requires specifying input and output files and "
        + "optionally a log file.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Path to the reaction files.\n"
        + "Example: --input_file ./path/to/input/reactions.csv",
    )
    parser.add_argument(
        "--output_file",
        type=str,
        required=True,
        help="Path to the output file where balanced reaction files "
        + "will be saved.\nExample: --output_file ./path/to/output/reactions_balance.csv",
    )
    parser.add_argument(
        "--log",
        type=str,
        help="Optional log file name to record the processing progress.\n"
        + "Example: --log process_log.txt",
    )

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.log)
