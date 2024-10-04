import argparse
import logging
import os
import pandas as pd
from syntemp.SynUtils.utils import load_database
from syntemp.auto_template import AutoTemp


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run the AutoTemp extraction process.")
    parser.add_argument(
        "--data_path", type=str, default=None, help="Path to the data file"
    )
    parser.add_argument("--rebalancing", action="store_true", help="Enable rebalancing")
    parser.add_argument(
        "--mapper_types",
        nargs="+",
        default=["local_mapper", "rxn_mapper", "graphormer"],
        help="Types of atom map techniques used",
    )
    parser.add_argument("--id", type=str, default="R-id", help="ID column")
    parser.add_argument(
        "--rsmi", type=str, default="reactions", help="Reaction SMILES column"
    )
    parser.add_argument(
        "--n_jobs", type=int, default=1, help="Number of jobs to run in parallel"
    )
    parser.add_argument("--verbose", type=int, default=2, help="Verbosity level")
    parser.add_argument(
        "--batch_size", type=int, default=1000, help="Batch size for processing"
    )
    parser.add_argument("--safe_mode", action="store_true", help="Enable safe mode")
    parser.add_argument(
        "--job_timeout", type=int, default=None, help="Timeout for jobs"
    )
    parser.add_argument(
        "--fix_hydrogen", action="store_true", help="Enable fixing hydrogen"
    )
    parser.add_argument(
        "--refinement_its", action="store_true", help="Refine non-equivalent ITS"
    )
    parser.add_argument(
        "--fast_process",
        action="store_true",
        help="Ignore hydrogen adjustment" + " with num_h >= 5",
    )
    parser.add_argument(
        "--rerun_aam", action="store_true", help="Run AAM based on mapper types input"
    )
    parser.add_argument("--lib_path", type=str, default=None, help="Library path")
    parser.add_argument("--save_dir", type=str, default=None, help="Save directory")
    parser.add_argument(
        "--log_file", type=str, default=None, help="File to log the process"
    )
    parser.add_argument(
        "--log_level", type=str, default="INFO", help="File to log the process"
    )
    parser.add_argument(
        "--get_random_hydrogen",
        action="store_true",
        help="Get random full ITS hydrogen",
    )
    return parser.parse_args()


def read_data(filepath):
    """Load data from a JSON or CSV file."""
    file_ext = os.path.splitext(filepath)[1].lower()
    try:
        if file_ext == ".gz":
            return load_database(filepath)
        elif file_ext == ".csv":
            return pd.read_csv(filepath)
        else:
            raise ValueError(f"Unsupported file type: {file_ext}")
    except Exception as e:
        logging.error(f"Error loading data from {filepath}: {e}")
        raise  # Re-raise the exception to handle it further up in the call stack


def main():
    args = parse_arguments()

    data = read_data(args.data_path)

    try:
        auto = AutoTemp(
            rebalancing=args.rebalancing,
            mapper_types=args.mapper_types,
            id=args.id,
            rsmi=args.rsmi,
            n_jobs=args.n_jobs,
            verbose=args.verbose,
            batch_size=args.batch_size,
            job_timeout=args.job_timeout,
            safe_mode=args.safe_mode,
            save_dir=args.save_dir,
            fix_hydrogen=args.fix_hydrogen,
            refinement_its=args.refinement_its,
            rerun_aam=args.rerun_aam,
            log_file=args.log_file,
            log_level=args.log_level,
            get_random_hydrogen=args.get_random_hydrogen,
            fast_process=args.fast_process,
        )
        auto.temp_extract(data, lib_path=args.lib_path)
        logging.info("Extraction successful.")
    except Exception as e:
        logging.error(f"An error occurred during processing: {e}")


if __name__ == "__main__":
    main()
