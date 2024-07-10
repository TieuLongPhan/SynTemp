import os
import shutil
from uuid import uuid4
from contextlib import contextmanager
from typing import List
import logging

logger = logging.getLogger("RDT")


@contextmanager
def temporary_change_dir(target_directory):
    """
    Context manager to temporarily change the working directory.
    """
    original_directory = os.getcwd()
    try:
        os.chdir(target_directory)
        yield
    finally:
        os.chdir(original_directory)


def map_with_rdt(reaction_smiles: str, rdt_jar_path: str, working_dir: str) -> str:
    """
    Maps a reaction using the RDT (Retro-Data) tool.

    Parameters:
    - reaction_smiles (str): The SMILES string of the reaction to be mapped.
    - rdt_jar_path (str): The file path to the RDT JAR file.
    - working_dir (str): The working directory to execute RDT.

    Returns:
    - str: The mapped reaction SMILES string.
    """
    unique_id = uuid4()
    unique_dir = os.path.join(working_dir, f"RDT_{unique_id}")
    output_file = "ECBLAST_smiles_AAM.txt"

    try:
        os.makedirs(unique_dir)

        with temporary_change_dir(unique_dir):
            command = (
                f'java -jar "{rdt_jar_path}" -Q SMI -q "{reaction_smiles}"'
                + f" -c -j AAM -f TEXT > {output_file}"
            )
            os.system(command)

            # Read the output from the file
            with open(output_file, "r") as inputFile:
                RDTaam = inputFile.read().splitlines()[3]

        mapped_smiles = RDTaam.split(" ")[0] if " " in RDTaam else RDTaam
        if ">>" in mapped_smiles:
            return mapped_smiles
        else:
            logger.error("Error in original SMILES, could have ~ in reaction smiles")
            return reaction_smiles
    except Exception as e:
        logger.error(f"RDT mapping failed: {e}")
        return reaction_smiles
    finally:
        shutil.rmtree(unique_dir, ignore_errors=True)


def map_with_rdt_batch(
    reaction_list: List[str], rdt_jar_path: str, working_dir: str, batch_size: int = 200
) -> List[str]:
    """
    Maps a batch of reactions using the RDT (Retro-Data) tool with
    a specified batch size and level of concurrency.

    Parameters:
    - reaction_list (List[str]): A list of SMILES strings of reactions to be mapped.
    - rdt_jar_path (str): The path to the RDT JAR file.
    - working_dir (str): The working directory for RDT executions.
    - batch_size (int): Number of reactions to process in a sub-batch.
    - processes (int): Number of concurrent processes to use for mapping.

    Returns:
    - List[str]: A list of mapped reaction SMILES strings, retains original on failure.
    """

    mapped_rxns = []
    total_reactions = len(reaction_list)
    # Split the reaction list into manageable sub-batches
    for i in range(0, total_reactions, batch_size):
        current_batch = reaction_list[i : i + batch_size]
        for reaction_smiles in current_batch:
            try:
                mapped_rxn = map_with_rdt(reaction_smiles, rdt_jar_path, working_dir)
                mapped_rxns.append(mapped_rxn)
            except Exception as e:
                logger.error(
                    f"RDT mapping failed for reaction '{reaction_smiles}': {e}",
                    exc_info=True,
                )
                mapped_rxns.append(reaction_smiles)
        logger.info(f"Successfully processed Internal Batch from {i} to {i+batch_size}")

    print("Completed mapping all reactions.")
    return mapped_rxns
