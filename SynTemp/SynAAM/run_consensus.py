import pathlib
import sys
from typing import List, Optional, Dict
import pandas as pd
import warnings
from joblib import Parallel, delayed
import transformers
from rxnmapper import RXNMapper
from SynTemp.SynChemistry.balance_checker import BalanceReactionCheck
from SynTemp.SynAAM.consensus_aam import ConsensusAAM
from SynTemp.SynUtils.utils import save_database

root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))

# Configure warnings and logging
warnings.filterwarnings("ignore")
transformers.logging.set_verbosity_error()


def map_batch(
    batch_data: List[Dict],
    rxn_mapper: RXNMapper,
    rsmi_column: str = "reactions",
    mapper_types: List[str] = ["rxn_mapper", "graphormer", "local_mapper"],
    rdt_jar_path: Optional[str] = None,
    working_dir: Optional[str] = None,
) -> pd.DataFrame:
    """
    Map a batch of reactions using the given mappers.

    Parameters:
    - batch_data: List[Dict], a list of dictionaries each containing a reaction.
    - rxn_mapper: RXNMapper, the RXNMapper object used for reaction mapping.
    - mapper_types: List[str], list of mapper types to use for consensus AAM.
    - rdt_jar_path: Optional[str], the file path to the
                                RDT JAR file if RDT mapper is used.
    - working_dir: Optional[str], the working directory to execute
                    external tools like RDT.

    Returns:
    - pd.DataFrame: A DataFrame with the batch of mapped reactions.
    """
    consensus_aam = ConsensusAAM(
        batch_data,
        rsmi_column=rsmi_column,
        save_dir=f"{root_dir}/Data",
        mapper_types=mapper_types,
    )
    return consensus_aam.fit(len(batch_data), rxn_mapper, rdt_jar_path, working_dir)


def run_consensus_aam(
    data: List[Dict],
    rsmi_column: str = "reactions",
    save_dir: Optional[str] = None,
    mapper_types: List[str] = ["rxn_mapper", "graphormer", "local_mapper"],
    data_name: str = "",
    batch_size: int = 1000,
    check_balance: bool = True,
    n_jobs: int = 1,
    verbose: int = 0,
    rdt_jar_path: Optional[str] = None,
    working_dir: Optional[str] = None,
) -> List[Dict]:
    """
    Automatically maps a dataset of chemical reactions using
    Atom Atom-Mapping (AAM) and returns the mapped reactions.

    Parameters:
    - data: List[Dict], the input list of dictionaries
    - rsmi_column: str, the name of reaction smiles key in the dictionaries/
    - save_dir: Optional[str], the directory to save the mapped
                                reactions to. If None, do not save.
    - mapper_types: List[str], the types of mappers to use for consensus AAM.
    - data_name: str, the name to give the saved file.
    - batch_size: int, the batch size to use when running the mapper.
    - check_balance: bool, if True, checks the reaction balance before running the mapper.
    - n_jobs: int, the number of parallel jobs to run.
    - verbose: int, the verbosity level.
    - rdt_jar_path: Optional[str], the file path to the
                    RDT JAR file if RDT mapper is used.
    - working_dir: Optional[str], the working directory
                to execute external tools like RDT.

    Returns:
    - List[Dict]: A list of dictionaries containing the mapped reactions.
    """
    rxn_mapper = RXNMapper()

    if check_balance:
        df_data = pd.DataFrame(data)
        checker = BalanceReactionCheck(
            df_data, rsmi_column=rsmi_column, n_jobs=5, verbose=verbose
        )
        balanced_reactions_df, _ = checker.check_balances()
        balanced_reactions = balanced_reactions_df.to_dict("records")
    else:
        balanced_reactions = data

    # Split the data into batches
    batches = [
        balanced_reactions[i: i + batch_size]
        for i in range(0, len(balanced_reactions), batch_size)
    ]

    # Parallel map the batches
    mapped_reactions = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(map_batch)(
            batch, rxn_mapper, rsmi_column, mapper_types, rdt_jar_path, working_dir
        )
        for batch in batches
    )

    mapped_reactions = [item for sublist in mapped_reactions for item in sublist]

    if save_dir:
        save_database(
            mapped_reactions,
            f"{save_dir}/{data_name}_aam_reactions.json.gz",
        )

    return mapped_reactions
