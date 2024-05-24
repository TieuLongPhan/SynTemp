import logging
from typing import List, Dict, Tuple, Callable
from importlib import import_module

logger = logging.getLogger("Consensus")


class AAMConsensus:
    """
    A class to handle consensus mapping of atom mappings using various mapping tools.

    Attributes:
        data (List[Dict]): The dataset containing reaction information.
        mappers (List[str]): The list of names of the mappers to use.

    Methods:
        import_mapper(mapper_name: str): Dynamically imports the mapper functions.
        single_consensus(reaction_dict: Dict, rsmi_column: str = 'reactions'): Processes a single reaction dictionary using specified mappers.
        batch_consensus(reaction_dicts: List[Dict], rsmi_column: str = 'reactions', batch_size: int = 200, job_timeout: int = 100): Processes a batch of reactions using specified mappers.
    """

    def __init__(self, data: List[Dict], mappers: List[str] = None):
        """
        Initializes the AAMConsensus object with the provided data and mappers.

        Parameters:
            data (List[Dict]): The list of reaction dictionaries to process.
            mappers (List[str], optional): A list of strings specifying which mappers to use. Defaults to ['local_mapper', 'rxn_mapper', 'graphormer'] if None is provided.
        """
        self.data = data
        # Default mappers if none provided
        self.mappers = (
            mappers
            if mappers is not None
            else ["local_mapper", "rxn_mapper", "graphormer"]
        )

    def import_mapper(self, mapper_name: str) -> Tuple[Callable, Callable]:
        """
        Dynamically imports the necessary mapper functions based on the mapper name.

        Parameters:
            mapper_name (str): The name of the mapper to import.

        Returns:
            Tuple[Callable, Callable]: A tuple containing the single and batch mapping functions.

        Raises:
            ValueError: If an unsupported mapper name is provided.
        """
        if mapper_name == "local_mapper":
            mapper_module = import_module("SynTemp.SynAAM.local_mapper_wrapper")
            return (
                mapper_module.map_with_local_mapper,
                mapper_module.map_with_local_mapper_batch,
            )
        elif mapper_name == "rxn_mapper":
            mapper_module = import_module("SynTemp.SynAAM.rxn_mapper_wrapper")
            return (
                mapper_module.map_with_rxn_mapper,
                mapper_module.map_with_rxn_mapper_batch,
            )
        elif mapper_name == "graphormer":
            mapper_module = import_module("SynTemp.SynAAM.graphormer_wrapper")
            return (
                mapper_module.map_with_graphormer,
                mapper_module.map_with_graphormer_batch,
            )
        elif mapper_name == "rdt":
            mapper_module = import_module("SynTemp.SynAAM.rdt_wrapper")
            return mapper_module.map_with_rdt, mapper_module.map_with_rdt_batch
        else:
            raise ValueError("Unsupported mapper name")

    def single_consensus(
        self,
        reaction_dict: Dict,
        rsmi_column: str = "reactions",
        rdt_jar_path: str = None,
        working_dir: str = None,
    ) -> Dict:
        """
        Processes a single reaction dictionary using the specified mappers and returns the results.

        Parameters:
            reaction_dict (Dict): The reaction dictionary to process.
            rsmi_column (str): The key in the dictionary where the reaction SMILES string is stored.

        Returns:
            Dict: The reaction dictionary updated with the mapping results from each mapper.
        """
        results_dict = reaction_dict.copy()
        rsmi = reaction_dict[rsmi_column]
        for mapper_name in self.mappers:
            mapper_func, _ = self.import_mapper(mapper_name)
            try:
                if mapper_name == "rdt":
                    results_dict[mapper_name] = mapper_func(
                        reaction_smiles=rsmi,
                        rdt_jar_path=rdt_jar_path,
                        working_dir=working_dir,
                    )
                else:
                    results_dict[mapper_name] = mapper_func(rsmi)
            except Exception as e:
                logger.error(f"Error mapping with {mapper_name}: {e}")
                results_dict[mapper_name] = (
                    rsmi  # Fallback to original SMILES on failure
                )
        return results_dict

    def batch_consensus(
        self,
        reaction_dicts: List[Dict],
        rsmi_column: str = "reactions",
        batch_size: int = 200,
        job_timeout: int = 100,
        safe_mode=True,
        rdt_jar_path: str = None,
        working_dir: str = None,
    ) -> List[Dict]:
        """
        Processes a batch of reactions using the specified mappers and returns the results.

        Parameters:
            reaction_dicts (List[Dict]): A list of reaction dictionaries to process.
            rsmi_column (str): The key where the reaction SMILES string is stored in each dictionary.
            batch_size (int): The number of reactions to process in each batch.
            job_timeout (int): The timeout in seconds for the local_mapper batch processing.

        Returns:
            List[Dict]: The list of reaction dictionaries updated with the mapping results from each mapper.
        """
        rsmis = [d[rsmi_column] for d in reaction_dicts]
        results = {}
        for mapper_name in self.mappers:
            logger.info(f"{mapper_name} process...")
            mapper_func, mapper_batch_func = self.import_mapper(mapper_name)
            try:
                if mapper_name == "local_mapper":
                    if safe_mode:
                        results[mapper_name] = mapper_batch_func(
                            rsmis, batch_size=batch_size, job_timeout=job_timeout
                        )
                    else:
                        results[mapper_name] = mapper_func(rsmis)
                elif mapper_name == "rdt":
                    results[mapper_name] = mapper_batch_func(
                        rsmis,
                        rdt_jar_path=rdt_jar_path,
                        working_dir=working_dir,
                        batch_size=batch_size,
                    )
                else:
                    results[mapper_name] = mapper_batch_func(
                        rsmis, batch_size=batch_size
                    )

            except Exception as e:
                logger.error(f"Error mapping batch with {mapper_name}: {e}")
                results[mapper_name] = [
                    d[rsmi_column] for d in reaction_dicts
                ]  # Fallback to original SMILES on failure

        # Apply mapping results to each reaction dictionary
        for rd in reaction_dicts:
            for mapper_name in self.mappers:
                rd[mapper_name] = results[mapper_name][reaction_dicts.index(rd)]
        return reaction_dicts
