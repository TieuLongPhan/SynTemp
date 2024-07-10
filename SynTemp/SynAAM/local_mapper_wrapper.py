from multiprocessing import Pool, TimeoutError
from localmapper import localmapper
import logging

logger = logging.getLogger("LocalMapper")
logging.basicConfig(level=logging.INFO)


def map_with_local_mapper(batch):
    """
    Process a batch of reactions using the local mapper,
    handling each reaction individually. If a single reaction is passed,
    it returns a single mapped reaction string. If a list is passed,
    it returns a list of mapped reactions.

    Parameters:
    - batch (str or list): A single reaction string or a batch of reactions
    to be processed.

    Returns:
    - str or list: The mapping results corresponding to the input format,
    with original reaction retained if an error occurs.
    """
    mapper = localmapper()

    single_input = False
    if not isinstance(batch, list):
        batch = [batch]
        single_input = True

    results = []
    for reaction in batch:
        try:
            result = mapper.get_atom_map(reaction)
            results.append(result)
        except Exception as e:
            logger.error(
                f"Failed to process reaction '{reaction}' due to error: {e}",
                exc_info=True,
            )
            results.append(reaction)

    return results[0] if single_input else results


def map_with_local_mapper_batch(reaction_list, batch_size=200, job_timeout=None):
    """
    Maps reactions in batches using a thread pool.

    Parameters:
    - reaction_list (list): List of SMILES strings representing reactions.
    - batch_size (int): Number of reactions to process in a batch. Defaults to 200.
    - job_timeout (int): Timeout in seconds for processing a batch. Defaults to 100.

    Returns:
    - list: List of mapped SMILES strings.
    """
    max_size = len(reaction_list)
    results_all = []
    if not job_timeout:
        job_timeout = int(batch_size / 2) + 5
    pool = Pool(processes=1)  # Initialize the pool outside the loop
    try:
        for i in range(0, max_size, batch_size):
            current_batch = reaction_list[i : i + batch_size]
            try:
                async_result = pool.apply_async(map_with_local_mapper, (current_batch,))
                results = async_result.get(job_timeout)
                results_all.extend(results)
                logger.info(f"Successfully processed Internal Batch {i}")
            except TimeoutError:
                logger.error(
                    f"Timeout: Batch {i} processing exceeded {job_timeout} seconds."
                )
                pool.terminate()  # Terminate the problematic pool.
                pool.join()  # Ensure all workers are cleaned up properly
                for reaction in current_batch:
                    with Pool(processes=1) as fallback_pool:
                        try:
                            result = fallback_pool.apply_async(
                                map_with_local_mapper, (reaction,)
                            ).get(5)
                            results_all.append(result)
                        except TimeoutError:
                            fallback_pool.terminate()
                            fallback_pool.join()
                            logger.error(
                                "Timeout: Individual reaction processing"
                                + f"exceeded limits for {reaction}."
                            )
                            results_all.append(reaction)
                logger.info(
                    f"Successfully processed Internal Batch from {i} to {i+batch_size}"
                )
                pool = Pool(processes=1)
    finally:
        pool.close()
        pool.join()

    return results_all
