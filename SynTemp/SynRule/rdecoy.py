from typing import Optional, List, Tuple, Dict
from SynTemp.SynUtils.chemutils import standardize_rsmi
import copy
from SynTemp.SynRule.rule_engine import RuleEngine


class RDecoy:
    @staticmethod
    def categorize_reactions(
        reactions: List[str], target_reaction: str
    ) -> Tuple[List[str], List[str]]:
        """
        Sorts a list of reaction SMILES strings into two groups based on their match with a specified target reaction. The
        categorization process distinguishes between reactions that align with the target reaction and those that do not.

        Parameters:
        - reactions (List[str]): The array of reaction SMILES strings to be categorized.
        - target_reaction (str): The SMILES string of the target reaction used as the benchmark for categorization.

        Returns:
        - Tuple[List[str], List[str]]: A pair of lists, where the first contains reactions matching the target and the second
                                    comprises non-matching reactions.
        """
        match, not_match = [], []
        target_reaction = standardize_rsmi(target_reaction, stereo=False)
        for reaction_smiles in reactions:
            if reaction_smiles == target_reaction:
                match.append(reaction_smiles)
            else:
                not_match.append(reaction_smiles)
        return match, not_match

    @staticmethod
    def generate_reactions(
        database: List[Dict],
        rule_file_path: str,
        original_rsmi_col: str = "reactions",
        repeat_times: int = 1,
        cluster_col: Optional[str] = "Cluster",
        reaction_side_index: int = 0,  # Assuming you want the index to be passed or set a default value
    ) -> Tuple[List[Dict], List[Dict]]:

        updated_database_forward = copy.deepcopy(database)
        for entry in updated_database_forward:
            entry[original_rsmi_col] = standardize_rsmi(
                entry[original_rsmi_col], stereo=False
            )
            entry["positive_reactions"] = []
            entry["negative_reactions"] = []

            rule_files = [f"{rule_file_path}/{entry[cluster_col]}.gml"]

            for rule_file in rule_files:
                initial_smiles_list = (
                    entry[original_rsmi_col].split(">>")[reaction_side_index].split(".")
                )

                reactions = RuleEngine.perform_reaction(
                    rule_file_path=rule_file,
                    initial_smiles=initial_smiles_list,
                    repeat_times=repeat_times,
                    prediction_type="forward",
                )

                reactions = list(
                    set([standardize_rsmi(value, stereo=True) for value in reactions])
                )
                matched_reactions, unmatched_reactions = RDecoy.categorize_reactions(
                    reactions, entry[original_rsmi_col]
                )

                if matched_reactions:
                    entry["positive_reactions"].extend(matched_reactions)
                # Optional: Uncomment and handle unmatched_reactions
                if unmatched_reactions:
                    entry["negative_reactions"].extend(unmatched_reactions)
                #
            entry["positive_reactions"] = list(set(entry["positive_reactions"]))
            entry["positive_reactions"] = (
                entry["positive_reactions"][0] if entry["positive_reactions"] else None
            )

            entry["negative_reactions"] = list(set(entry["negative_reactions"]))
            entry["negative_reactions"] = (
                entry["negative_reactions"] if entry["negative_reactions"] else None
            )

        return updated_database_forward
