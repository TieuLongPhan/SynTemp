from SynTemp.SynAAM.atom_map_consensus import AAMConsensus
from SynTemp.SynITS.its_extraction import ITSExtraction


def run_aam(data, mapper_types):
    """
    Executes the atom mapping consensus algorithm on provided data using specified mapper types.

    Args:
    - data (list): List of data records to process.
    - mapper_types (list): List of mapper types to use.

    Returns:
    - results (list): Results from the atom mapping consensus process.
    """
    aam = AAMConsensus(data, mappers=mapper_types)
    results = aam.batch_consensus(
        data,
        rsmi_column="record",
        batch_size=len(data),
        job_timeout=None,
        safe_mode=False,
    )
    return results


def low_cost_its(
    data, mapper_name, batch_size=1000, verbose=1, n_jobs=4, check_method="ITS"
):
    """
    Processes input data to extract ITS using specified mapper, handling data in batches for efficiency.

    Args:
    - data (list): Data to process.
    - mapper_name (str): Name of the mapper to use.
    - batch_size (int): Size of each data batch.
    - verbose (int): Verbosity level of the process.
    - n_jobs (int): Number of jobs to run in parallel.
    - check_method (str): Method used to check ITS correctness.

    Returns:
    - Tuple containing updated data with 'equivalent' flags, correct ITS, and incorrect ITS.
    """
    num_batches = (len(data) + batch_size - 1) // batch_size
    its_correct, its_incorrect = [], []

    for batch_start in range(0, len(data), batch_size):
        batch_end = batch_start + batch_size
        batch_data = data[batch_start:batch_end]
        batch_correct, batch_incorrect = ITSExtraction.parallel_process_smiles(
            batch_data,
            mapper_name,
            n_jobs=n_jobs,
            verbose=verbose,
            export_full=False,
            check_method=check_method,
        )
        its_correct.extend(batch_correct)
        its_incorrect.extend(batch_incorrect)

    # Update data with equivalence information
    equivalent_ids = {item["R-id"] for item in its_correct}
    for item in data:
        item["equivalent"] = item["R-id"] in equivalent_ids

    return data, its_correct, its_incorrect
