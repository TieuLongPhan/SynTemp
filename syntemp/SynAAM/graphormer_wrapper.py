import logging
from chython import smiles

logger = logging.getLogger("graphormer")
logging.basicConfig(level=logging.INFO)


def map_with_graphormer(reaction_smiles: str) -> str:
    """
    Maps a reaction using the Graphormer.

    Parameters:
    - reaction_smiles (str): The SMILES string of the reaction to be mapped.

    Returns:
        str: The mapped reaction SMILES string or the original string if mapping fails.
    """
    try:
        # Initialize Graphormer with reaction SMILES
        myrxn = smiles(reaction_smiles)
        myrxn.reset_mapping()
        # Get mapped reaction
        mapped_rxn = format(myrxn, "m")
        return mapped_rxn.split(" ")[0] if " " in mapped_rxn else mapped_rxn
    except Exception as e:
        logger.error(
            f"Graphormer mapping failed for reaction '{reaction_smiles}': {e}",
            exc_info=True,
        )
        return reaction_smiles


def map_with_graphormer_batch(reaction_list: list, batch_size: int = 200) -> list:
    """
    Maps a batch of reactions using the Graphormer with a specified batch size.

    Parameters:
    - reaction_list (list): A list of SMILES strings of reactions to be mapped.
    - batch_size (int): Number of reactions to process in a sub-batch.

    Returns:
    - list: A list of mapped reaction SMILES strings, retains original on failure.
    """
    mapped_rxns = []
    total_reactions = len(reaction_list)
    logger.info(
        f"Starting mapping of {total_reactions} reactions in batches of {batch_size}."
    )

    for i in range(0, total_reactions, batch_size):
        current_batch = reaction_list[i : i + batch_size]
        logger.info(
            f"Processing batch from index {i} to {min(i + batch_size, total_reactions)}."
        )

        for reaction_smiles in current_batch:
            try:
                mapped_rxn = map_with_graphormer(reaction_smiles)
                mapped_rxns.append(mapped_rxn)
            except Exception as e:
                logger.error(
                    f"Graphormer mapping failed for reaction '{reaction_smiles}': {e}",
                    exc_info=True,
                )
                mapped_rxns.append(reaction_smiles)
        logger.info(f"Successfully processed Internal Batch from {i} to {i+batch_size}")

    logger.info("Completed mapping all reactions.")
    return mapped_rxns
