from rxnmapper import RXNMapper
import logging

logger = logging.getLogger("rxnmapper")
logging.basicConfig(level=logging.INFO)

rxn_mapper = RXNMapper()


def map_with_rxn_mapper(
    reaction_smiles: str, rxn_mapper: RXNMapper = rxn_mapper
) -> str:
    """
    Maps a reaction using the RXNMapper.

    Parameters:
    - reaction_smiles (str): The SMILES string of the reaction to be mapped.
    - rxn_mapper (RXNMapper): The instance of RXNMapper used for mapping.

    Returns:
    - str: The mapped reaction SMILES string, or the original string if mapping fails.
    """
    try:
        mapped_rxn = rxn_mapper.get_attention_guided_atom_maps(
            [reaction_smiles], canonicalize_rxns=False
        )[0]["mapped_rxn"]
        return mapped_rxn.split(" ")[0] if " " in mapped_rxn else mapped_rxn
    except Exception as e:
        logger.error(
            f"RXNMapper mapping failed for reaction '{reaction_smiles}': {e}",
            exc_info=True,
        )
        return reaction_smiles


def map_with_rxn_mapper_batch(
    reaction_list: list, rxn_mapper: RXNMapper = rxn_mapper, batch_size: int = 200
) -> list:
    """
    Maps a batch of reactions using the RXNMapper with a specified batch size.

    Parameters:
    - reaction_list (list): A list of SMILES strings of reactions to be mapped.
    - rxn_mapper (RXNMapper): The instance of RXNMapper used for mapping.
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
                mapped_rxn = rxn_mapper.get_attention_guided_atom_maps(
                    [reaction_smiles], canonicalize_rxns=False
                )[0]["mapped_rxn"]
                mapped_rxns.append(
                    mapped_rxn.split(" ")[0] if " " in mapped_rxn else mapped_rxn
                )
            except Exception as e:
                logger.error(
                    f"RXNMapper mapping failed for reaction '{reaction_smiles}': {e}",
                    exc_info=True,
                )
                mapped_rxns.append(
                    reaction_smiles
                )  # Append the original SMILES on failure
        logger.info(f"Successfully processed Internal Batch from {i} to {i+batch_size}")
    logger.info("Completed mapping all reactions.")
    return mapped_rxns
