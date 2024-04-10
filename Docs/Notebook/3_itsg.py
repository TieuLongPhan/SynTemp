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
from SynTemp.SynUtils.utils import load_database, save_database, save_to_pickle, load_from_pickle
from SynTemp.SynAAM.aam_postprocess import AMMPostprocessor
from SynTemp.SynITS.its_extraction import ITSExtraction
from SynTemp.SynITS.its_hadjuster import ITSHAdjuster
from SynTemp.SynMØD.naive_cluster import NaiveCluster
from SynTemp.SynITS.uncertain_refinement import UncertainRefinement

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
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # File handler
    if save_dir:
        log_file = os.path.join(save_dir, f'{data_name}_pipeline_log.log')
        fh = logging.FileHandler(log_file, mode = 'w')
        fh.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger

def process_batch(batch, batch_index, temp_dir, n_jobs, verbose):
    result = UncertainRefinement.process_graphs_in_parallel(batch, n_jobs=n_jobs, verbose=verbose)
    result = [value for value in result if value is not None]
    # Save the result to a file in the temporary directory
    temp_file_path = os.path.join(temp_dir, f"batch_{batch_index}.pkl.gz")
    save_to_pickle(result, temp_file_path)


def run_synitsg_pipeline(data, mapper_name=None, batch_size=1000, verbose=1, n_jobs=4, check_valid=False, fix_hydrogen=True, 
                         curate_uncertain_mapping=True, alignment=True, save_dir=None, data_name = ''):
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
    logger = configure_logging(save_dir, verbose, data_name)
    logger.setLevel(logging.DEBUG if verbose > 1 else logging.INFO)
    logger.info(f"Processing data {data_name} with {n_jobs} cpus")
    start_time = time.time()
    
    if not mapper_name:
        mapper_name = ['rxn_mapper', 'graphormer', 'local_mapper']
        logger.info("No mapper_name provided. Using default mappers.")
    threshold = len(mapper_name)-1
    if check_valid:
        logger.info("Checking and filtering valid data.")
        check_valid_data = AMMPostprocessor.parallel_postprocess(data, mapper_name, threshold=3, n_jobs=n_jobs, verbose=verbose)
        valid_data = [reaction for reaction in check_valid_data if reaction.get('Valid')]
        unvalid_data = [reaction for reaction in check_valid_data if not reaction.get('Valid')]
        logger.info(f"Valid data count: {len(valid_data)}, Invalid data count: {len(unvalid_data)}")
        data = valid_data 

        
    temp_dir = os.path.join(save_dir, "temp_batches") if save_dir else "temp_batches"
    os.makedirs(temp_dir, exist_ok=True)

    num_batches = len(data) // batch_size + (len(data) % batch_size > 0)
    for i in range(num_batches):
        logger.debug(f"Processing batch {i+1}/{num_batches}")
        start_index = i * batch_size
        end_index = min((i + 1) * batch_size, len(data))
        batch_data = data[start_index:end_index]

        batch_correct, batch_incorrect = ITSExtraction.parallel_process_smiles(batch_data, mapper_name, threshold=threshold, n_jobs=n_jobs, verbose=verbose, export_full=False, check_method='RC')
        
        if fix_hydrogen:
            logger.info(f"Fixing hydrogen for batch {i+1}/{num_batches}.")
            batch_processed = ITSHAdjuster.process_graph_data_parallel(batch_correct, 'ITSGraph', n_jobs=n_jobs, verbose=verbose)
            uncertain_hydrogen = [value for value in batch_processed if value['ITSGraph'] is None]
            batch_correct = [value for value in batch_processed if value['ITSGraph'] is not None]
            # Save uncertain hydrogen data for the batch
            save_to_pickle(uncertain_hydrogen, os.path.join(temp_dir, f'uncertain_hydrogen_{i}.pkl'))

        # Save batch results to temporary files
        save_to_pickle(batch_correct, os.path.join(temp_dir, f'batch_correct_{i}.pkl'))
        save_to_pickle(batch_incorrect, os.path.join(temp_dir, f'batch_incorrect_{i}.pkl'))

    # Combine saved batch data
    its_correct, its_incorrect, all_uncertain_hydrogen = [], [], []
    for i in range(num_batches):
        
        its_correct.extend(load_from_pickle(os.path.join(temp_dir, f'batch_correct_{i}.pkl')))
        its_incorrect.extend(load_from_pickle(os.path.join(temp_dir, f'batch_incorrect_{i}.pkl')))
        all_uncertain_hydrogen.extend(load_from_pickle(os.path.join(temp_dir, f'uncertain_hydrogen_{i}.pkl')))

    # Clean up temporary files
    shutil.rmtree(temp_dir)

    process_graph_data = its_correct
    # Save final combined results
    logger.info(f"Combining and save data")
    if save_dir:
        logger.info(f"Number of correct mappers:{len(its_correct)}")
        logger.info(f"Number of incorrect mappers:{len(its_incorrect)}")
        logger.info(f"Number of uncertain hydrogen:{len(all_uncertain_hydrogen)}")
        save_to_pickle(its_correct, os.path.join(save_dir, f'{data_name}_its_correct.pkl.gz'))
        save_to_pickle(its_incorrect, os.path.join(save_dir, f'{data_name}_its_incorrect.pkl.gz'))
        save_to_pickle(all_uncertain_hydrogen, os.path.join(save_dir, f'{data_name}_uncertain_hydrogen.pkl.gz'))

    if curate_uncertain_mapping:
        logger.info("Curating incorrect mappings.")

        temp_dir = os.path.join(save_dir, "temp_batches") if save_dir else "temp_batches"
        os.makedirs(temp_dir, exist_ok=True)

        # Determine batch size, for example, 100 items per batch
        batch_size = min(batch_size, 500)
        batches = [its_incorrect[i:i + batch_size] for i in range(0, len(its_incorrect), batch_size)]

        # Process each batch
        for i, batch in enumerate(batches):
            logger.info(f"Processing batch {i+1}/{len(batches)}")
            process_batch(batch, i, temp_dir, n_jobs, verbose)

        # Combine results from all batches
        uncertain_fix = []
        for i in range(len(batches)):
            uncertain_fix.extend(load_from_pickle(os.path.join(temp_dir, f'batch_{i}.pkl.gz')))
      
        process_graph_data.extend(uncertain_fix)
        logger.info(f"incorrect mappings fixed: {len(uncertain_fix)}")

        shutil.rmtree(temp_dir)

    if alignment:
        logger.info(f"Data all for alignment: {len(process_graph_data)}")
        logger.info("Clustering ITS graphs.")
        node_label_names = ["element", "charge"]
        naive_cluster = NaiveCluster(node_label_names=node_label_names, node_label_default=["*", 0], edge_attribute="order")
        its_graph_rules_cluster = naive_cluster.process_rules_clustering(process_graph_data, rule_column='GraphRules')
        if save_dir:
            save_to_pickle(its_graph_rules_cluster, f'{save_dir}/{data_name}_its_graph_rules_cluster.pkl.gz')
        logger.info(f"Number of clusters: {pd.DataFrame(its_graph_rules_cluster)['naive_cluster'].value_counts()}")
    else:
        its_graph_rules_cluster = process_graph_data

    elapsed_time = time.time() - start_time
    logger.info(f"Execution time: {elapsed_time:.2f} seconds")
    return its_graph_rules_cluster


if __name__ == "__main__":
    mapper_name = ['rxn_mapper', 'graphormer', 'local_mapper']
    #folder_names = ['uspto', 'jaworski', 'golden', 'ecoli']
    folder_name = 'USPTO_50K'
    save_dir = f'{root_dir}/Data/{folder_name}'

    data = load_database(f'{root_dir}/Data/{folder_name}/{folder_name}_aam_reactions.json.gz')[:]
    run_synitsg_pipeline(data, mapper_name, batch_size=500, verbose=1, n_jobs=4, check_valid=False, 
                         curate_uncertain_mapping=True, fix_hydrogen=True, 
                         alignment=True, save_dir=save_dir, data_name=folder_name)
  


