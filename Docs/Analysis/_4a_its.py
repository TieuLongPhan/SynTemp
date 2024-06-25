import pathlib
import time
import logging
import os
import sys
import shutil
from joblib import Parallel, delayed
import tempfile
import argparse



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
        log_file = os.path.join(save_dir, f"{data_name}_pipeline.log")
        fh = logging.FileHandler(log_file, mode="w")
        fh.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def run_synitsg_pipeline(
    data,
    mapper_name=None,
    batch_size=1000,
    verbose=1,
    n_jobs=4,
    fix_hydrogen=True,
    save_dir=None,
    data_name="",
):
    """
    Runs the Synthetic ITS Graph pipeline.

    Parameters:
        data (list): List of reaction data.
        mapper_name (list): List of mapper names.
        batch_size (int): Batch size for processing.
        verbose (int): Verbosity level. 0: no output, 1: summary output, 2: detailed output.
        n_jobs (int): Number of parallel jobs to run.
        check_valid (bool): Whether to check the validity of the data.
        fix_hydrogen (bool): Whether to fix hydrogen in the ITS graphs.
        curate_uncertain_mapping (bool): Whether to curate uncertain mappings.
        alignment (bool): Whether to align ITS graphs.
        save_dir (str): Directory to save results.
        data_name (str): Name of the data, used for naming the log file.

    Returns:
        list: List of processed ITS graph data.
    """
    from SynTemp.SynUtils.utils import (
        load_database,
        save_database,
        save_to_pickle,
        load_from_pickle,
    )
    from SynTemp.SynITS.its_extraction import ITSExtraction
    from SynTemp.SynITS.its_hadjuster import ITSHAdjuster


    logger = configure_logging(save_dir, verbose, data_name)
    logger.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
    logger.info(f"Processing data {data_name} with {n_jobs} cpus")
    start_time = time.time()

    if not mapper_name:
        mapper_name = ["rxn_mapper", "graphormer", "local_mapper"]
        logger.info("No mapper_name provided. Using default mappers.")

    temp_dir = os.path.join(save_dir, "temp_batches") if save_dir else "temp_batches"
    os.makedirs(temp_dir, exist_ok=True)

    num_batches = len(data) // batch_size + (len(data) % batch_size > 0)
    for i in range(num_batches):
        logger.debug(f"Processing batch {i+1}/{num_batches}")
        start_index = i * batch_size
        end_index = min((i + 1) * batch_size, len(data))
        batch_data = data[start_index:end_index]

        batch_correct, batch_incorrect = ITSExtraction.parallel_process_smiles(
            batch_data,
            mapper_name,
            n_jobs=n_jobs,
            verbose=verbose,
            export_full=False,
            check_method="RC",
        )

        if fix_hydrogen:
            logger.info(f"Fixing hydrogen for batch {i+1}/{num_batches}.")
            batch_processed = ITSHAdjuster.process_graph_data_parallel(
                batch_correct, "ITSGraph", n_jobs=n_jobs, verbose=verbose
            )
            uncertain_hydrogen = [
                value for value in batch_processed if value["ITSGraph"] is None
            ]
            batch_correct = [
                value for value in batch_processed if value["ITSGraph"] is not None
            ]
            # Save uncertain hydrogen data for the batch
            save_to_pickle(
                uncertain_hydrogen,
                os.path.join(temp_dir, f"uncertain_hydrogen_{i}.pkl"),
            )

        # Save batch results to temporary files
        save_to_pickle(batch_correct, os.path.join(temp_dir, f"batch_correct_{i}.pkl"))
        save_to_pickle(
            batch_incorrect, os.path.join(temp_dir, f"batch_incorrect_{i}.pkl")
        )

    # Combine saved batch data
    its_correct, its_incorrect, all_uncertain_hydrogen = [], [], []
    for i in range(num_batches):

        its_correct.extend(
            load_from_pickle(os.path.join(temp_dir, f"batch_correct_{i}.pkl"))
        )
        its_incorrect.extend(
            load_from_pickle(os.path.join(temp_dir, f"batch_incorrect_{i}.pkl"))
        )
        if fix_hydrogen:
            all_uncertain_hydrogen.extend(
                load_from_pickle(os.path.join(temp_dir, f"uncertain_hydrogen_{i}.pkl"))
            )

    # Clean up temporary files
    shutil.rmtree(temp_dir)

   
    # Save final combined results
    logger.info(f"Combining and save data")
    if save_dir:
        logger.info(f"Number of correct mappers:{len(its_correct)}")
        logger.info(f"Number of incorrect mappers:{len(its_incorrect)}")
        logger.info(f"Number of uncertain hydrogen:{len(all_uncertain_hydrogen)}")
        save_to_pickle(
            its_correct, os.path.join(save_dir, f"{data_name}_its_correct.pkl.gz")
        )
        save_to_pickle(
            its_incorrect, os.path.join(save_dir, f"{data_name}_its_incorrect.pkl.gz")
        )
        save_to_pickle(
            all_uncertain_hydrogen,
            os.path.join(save_dir, f"{data_name}_uncertain_hydrogen.pkl.gz"),
        )

    elapsed_time = time.time() - start_time
    logger.info(f"Execution time: {elapsed_time:.2f} seconds")
    

def main():
    root_dir = pathlib.Path(__file__).parents[2]
    sys.path.append(str(root_dir))
    from SynTemp.SynUtils.utils import (
        load_database,
    )
    
    parser = argparse.ArgumentParser(description="Run the Synthetic ITS Graph pipeline.")
    parser.add_argument("--mapper_name", nargs='+', default=["rxn_mapper", "graphormer", "local_mapper"], help="List of mapper names")
    parser.add_argument("--batch_size", type=int, default=500, help="Batch size for processing")
    parser.add_argument("--verbose", type=int, default=1, help="Verbosity level")
    parser.add_argument("--n_jobs", type=int, default=4, help="Number of parallel jobs")
    parser.add_argument("--fix_hydrogen", type=bool, default=False, help="Whether to fix hydrogen")
    parser.add_argument("--data_name", type=str, default="", help="Name of the data")
    parser.add_argument("--rule_folder", type=str, default="", help="Name of folder to store rules")

    args = parser.parse_args()

    data = load_database(os.path.join(root_dir, 'Data', 'DPO', args.data_name, 'train.json.gz'))

    args.save_dir = os.path.join(root_dir, 'Data', 'DPO', args.data_name, args.rule_folder)
    run_synitsg_pipeline(
        data=data,
        mapper_name=args.mapper_name,
        batch_size=args.batch_size,
        verbose=args.verbose,
        n_jobs=args.n_jobs,
        fix_hydrogen=args.fix_hydrogen,
        save_dir=args.save_dir,
        data_name=args.data_name)


if __name__ == "__main__":

    main()
    
    #python Docs/Analysis/_4a_its.py --batch_size 1000 --data_name USPTO_unbalance --rule_folder Raw
    #python Docs/Analysis/_4a_its.py --batch_size 1000 -fix_hydrogen True --data_name USPTO_unbalance --rule_folder Complete
    #python Docs/Analysis/_4a_its.py --batch_size 1000 --fix_hydrogen True --data_name USPTO_balance --rule_folder Complete