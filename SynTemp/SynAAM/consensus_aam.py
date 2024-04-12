import logging
import os
import sys
import pathlib
import pandas as pd
from typing import List, Dict, Optional
from SynTemp.SynAAM.atom_mappers import (
    map_with_rxn_mapper,
    map_with_graphormer,
    map_with_local_mapper,
    map_with_rdt,
)
from SynTemp.SynUtils.utils import load_database, save_database
from rxnmapper import RXNMapper
from localmapper import localmapper
from SynTemp.SynProcessor.balance_checker import BalanceReactionCheck

root_dir = pathlib.Path(__file__).parents[2]
sys.path.append(str(root_dir))


class ConsensusAAM:
    """
    A class to apply atom-atom mapping (AAM) consensus from multiple mappers to a list of chemical reactions.

    Attributes:
        reactions_list (List[Dict]): A list of dictionaries containing reaction information.
        return_confidence (bool): Flag to indicate whether to return confidence scores.
        root_dir (str): Root directory for the local mapper.
        save_dir (str): Directory where processed batches will be saved.
        mapper_types (List[str]): List of mapper types to be used.
    """

    def __init__(
        self,
        reactions_list: List[Dict],
        rsmi_column: str = "reactions",
        save_dir: str = "./",
        mapper_types: List[str] = ["rxn_mapper", "graphormer", "local_mapper"],
    ):
        """
        Constructs all the necessary attributes for the ConsensusAAM object.

        Parameters:
            reactions_list (List[Dict]): List of dictionaries containing reaction information.
            return_confidence (bool): Flag to return confidence scores (default is False).
            root_dir (str): Root directory for the local mapper (default is '../SynITSG/LocalMapper').
            save_dir (str): Directory to save processed batches (default is current directory).
            mapper_types (List[str]): List of mapper types to be used (default includes rxn_mapper, graphormer, and local_mapper).
        """
        self.reactions_list = reactions_list
        self.rsmi_column = rsmi_column
        self.save_dir = save_dir
        self.mapper_types = mapper_types
        self.smiles_list = self.extract_smiles()
        self.filename = os.path.join(self.save_dir, f"aam_reactions.json.gz")

    def extract_smiles(self) -> List[str]:
        """Extracts SMILES strings from the reactions list."""
        return [
            reaction_dict[self.rsmi_column] for reaction_dict in self.reactions_list
        ]

    def process_batch(
        self, batch: List[str], rxn_mapper, rdt_jar_path: str, working_dir: str
    ) -> None:
        """
        Applies mapping functions to a batch of SMILES and updates the original dictionary with the results.
        Also saves the updated reactions_list after processing the batch.

        Parameters:
            batch (List[str]): Batch of SMILES strings to process.
            rxn_mapper: The reaction mapper object.
        """
        results = {}
        if "local_mapper" in self.mapper_types:
            results["local_mapper"] = [map_with_local_mapper(rxn) for rxn in batch]
        if "rxn_mapper" in self.mapper_types:
            results["rxn_mapper"] = [
                map_with_rxn_mapper(rxn, rxn_mapper) for rxn in batch
            ]
        if "graphormer" in self.mapper_types:
            results["graphormer"] = [map_with_graphormer(rxn) for rxn in batch]
        if "rdt" in self.mapper_types:
            results["rdt"] = [
                map_with_rdt(rxn, rdt_jar_path, working_dir) for rxn in batch
            ]

        batch_dict = {smiles: {} for smiles in batch}
        for mapper_type, mapper_results in results.items():
            for rxn, result in zip(batch, mapper_results):
                batch_dict[rxn][mapper_type] = result

        for reaction_dict in self.reactions_list:
            smiles = reaction_dict[self.rsmi_column]
            if smiles in batch_dict:
                reaction_dict.update(batch_dict[smiles])

        # Save the updated reactions_list
        save_database(self.reactions_list, self.filename)

    def fit(
        self,
        batch_size: int,
        rxn_mapper,
        rdt_jar_path: str = None,
        working_dir: str = None,
    ) -> List[Dict]:
        """
        Processes the SMILES list in batches and returns the updated list of reaction dictionaries.

        Parameters:
            batch_size (int): Number of reactions to process in each batch.
            rxn_mapper: The reaction mapper object.

        Returns:
            List[Dict]: Updated list of reaction dictionaries with mapping results.
        """
        for i in range(0, len(self.smiles_list), batch_size):
            batch = self.smiles_list[i : i + batch_size]
            self.process_batch(batch, rxn_mapper, rdt_jar_path, working_dir)
            batch_number = i // batch_size + 1
            logging.info(
                f"Processed batch {batch_number}/{len(self.smiles_list) // batch_size}"
            )
        return self.reactions_list
