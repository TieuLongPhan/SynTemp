import sys
import pathlib

root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import save_to_pickle, load_from_pickle
from SynTemp.SynITS.its_refinement import ITSRefinement

its_wrong = load_from_pickle(
    f"{root_dir}/Data/DPO/USPTO_balance/Complete/USPTO_balance_its_incorrect.pkl.gz"
)
its_correct = load_from_pickle(
    f"{root_dir}/Data/DPO/USPTO_balance/Complete/USPTO_balance_its_correct.pkl.gz"
)


def batch_process_items(data_items, batch_size):
    """
    Processes items in parallel in batches.

    :param data_items: List of items to be processed.
    :param batch_size: The number of items in each batch.
    """
    # Process each batch
    total_results = []
    for i in range(0, len(data_items), batch_size):
        batch = data_items[i : i + batch_size]
        print(
            f"Processing batch {i // batch_size + 1}/{(len(data_items) + batch_size - 1) // batch_size}"
        )
        result = ITSRefinement.process_graphs_in_parallel(batch, n_jobs=1)
        total_results.extend(result)
        # Optionally, handle the result here, such as saving it or further processing.
    return total_results


batch_size = 500  # Define the size of each batch
total_results = batch_process_items(its_wrong[:], batch_size)

new_data = [value for value in total_results if value]
its_correct.extend(new_data)
save_to_pickle(
    its_correct,
    f"{root_dir}/Data/DPO/USPTO_balance/Expand/USPTO_balance_its_correct.pkl.gz",
)
