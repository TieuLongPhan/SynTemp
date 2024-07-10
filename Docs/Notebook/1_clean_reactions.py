import pathlib
import sys
import time
import logging

root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import load_database, save_database
import warnings
from SynTemp.SynStandardizer.neutralize import Neutralize
from SynTemp.SynStandardizer.deionize import Deionize


warnings.filterwarnings("ignore")


def main():
    data = load_database(
        f"{root_dir}/Data/DPO/USPTO_balance/USPTO_50K_bad_reactions.json.gz"
    )[:]
    df = [{"R-id": value["id"], "reactions": value["reactions"]} for value in data]
    data = Neutralize.parallel_fix_unbalanced_charge(df, "reactions", 4)
    data = Deionize.apply_uncharge_smiles_to_reactions(data, Deionize.uncharge_smiles)
    data = [
        {"R-id": value["R-id"], "reactions": value["standardized_reactions"]}
        for value in data
    ]

    save_database(
        data, f"{root_dir}/Data/DPO/USPTO_balance/USPTO_50K_reactions.json.gz"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    start_time = time.time()
    main()
    end_time = time.time()

    elapsed_time = end_time - start_time
    logging.info(f"Execution time: {elapsed_time:.2f} seconds")
