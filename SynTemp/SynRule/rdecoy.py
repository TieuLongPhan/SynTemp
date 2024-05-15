import copy
from typing import Optional, List, Tuple, Dict
from SynTemp.SynUtils.chemutils import standardize_rsmi, categorize_reactions
from SynTemp.SynRule.rule_engine import RuleEngine


class RDecoy:

    @staticmethod
    def generate_reactions(
        database: List[Dict],
        rule_file_path: str,
        original_rsmi_col: str = "reactions",
        repeat_times: int = 1,
        cluster_col: Optional[str] = "Cluster",
        reaction_side_index: int = 0,
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
                matched_reactions, unmatched_reactions = categorize_reactions(
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
