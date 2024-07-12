import os
import glob
import logging
from typing import List
from SynTemp.SynComp.valence_constrain import ValenceConstrain
from SynTemp.SynUtils.graph_utils import load_gml_as_text
from mod import RCMatch, ruleGMLString

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class RuleCompose:
    def __init__(self) -> None:
        pass

    @staticmethod    
    def _compose(rule_1, rule_2):
        """
        Compose two rules and filter the results based on chemical valence constraints.

        Parameters:
        - rule_1: First rule object to compose.
        - rule_2: Second rule object to compose.

        Returns:
        - list: List of 'good' modifications where the resulting rules pass the
        valence checks.
        """
        try:
            # Attempt to match and compose the rules
            m = RCMatch(rule_1, rule_2)
            modRes = m.composeAll()
            valence_check = ValenceConstrain()
            goodMod, _ = valence_check.split(modRes)
            return goodMod
        except Exception as e:
            print(e)
            return []  # Return an empty list in case of failure

    @staticmethod
    def _process_compose(rule_1_id, rule_2_id, rule_path, rule_path_compose):
        """
        Process and compose two rules based on their GML files.

        Parameters:
        - rule_1_id (str): Identifier for the first rule.
        - rule_2_id (str): Identifier for the second rule.
        - rule_path (str): Directory path where the original GML files are stored.
        - rule_path_compose (str): Directory path where the composed GML files
        will be saved.

        Returns:
        - list: Composed rules from the two provided rules.
        """
        rule_1 = load_gml_as_text(f"{rule_path}/{rule_1_id}.gml")
        rule_1 = ruleGMLString(rule_1)
        rule_2 = ruleGMLString(load_gml_as_text(f"{rule_path}/{rule_2_id}.gml"))
        rules_compose = RuleCompose._compose(rule_1, rule_2)
        if rule_path_compose:
            for key, value in enumerate(rules_compose):
                filepath = f"{rule_path_compose}/p_{rule_1_id}_{rule_2_id}_r{key}.gml"
                RuleCompose.save_gml_from_text(
                    value.getGMLString(), filepath, key, [rule_1_id, rule_2_id]
                )
        return rules_compose

    @staticmethod
    def _auto_compose(rule_path, rule_path_compose):
        """
        Automatically find all GML files in the given directory and compose them pairwise.

        Parameters:
        - rule_path (str): Directory path where the GML files are stored.
        - rule_path_compose (str): Directory path where the composed GML files will
        be saved.

        Returns:
        - None: Composed rules are saved directly to the filesystem.
        """
        # Get all gml file names in the directory
        gml_files = [os.path.basename(f) for f in glob.glob(f"{rule_path}/*.gml")]
        gml_ids = [
            os.path.splitext(f)[0] for f in gml_files
        ]  # Strip the .gml extension to get IDs

        # Compose each pair of rules once (i.e., (rule1, rule2) but not (rule2, rule1))
        # Calculate the total number of compositions for progress logging
        num_files = len(gml_ids)
        total_compositions = num_files * (num_files - 1) // 2
        current_composition = 0
        for i in range(len(gml_ids)):
            for j in range(i + 1, len(gml_ids)):
                RuleCompose._process_compose(
                    gml_ids[i], gml_ids[j], rule_path, rule_path_compose
                )
                current_composition += 1
                if current_composition % 100 == 0:
                    logging.info(
                        f"Progress: {current_composition}/{total_compositions}"
                        + "compositions completed."
                    )

    @staticmethod
    def save_gml_from_text(
        gml_content: str, gml_file_path: str, rule_id: str, parent_ids: List[str]
    ) -> bool:
        """
        Save a text string to a GML file by modifying the 'ruleID' line to include parent
        rule names. This function parses the given GML content, identifies any lines
        starting with 'ruleID', and replaces these lines with a new ruleID that
        incorporates identifiers from parent rules.

        Parameters:
        - gml_content (str): The content to be saved to the GML file. This should be the
        entire textual content of a GML file.
        - gml_file_path (str): The file path where the GML file should be saved. If the
        path does not exist or is inaccessible, the function will return False and print
        an error message.
        - rule_id (str): The original rule ID from the content. This is the identifier
        that will be modified to include parent IDs in the new ruleID.
        - parent_ids (List[str]): List of parent rule IDs to prepend to the original rule
        ID. These are combined into a new identifier to reflect the hierarchical
        relationship in rule IDs.

        Returns:
        - bool: True if the file was successfully saved, False otherwise. The function
        attempts to write the modified GML content to the specified file path.
        """
        try:
            parent_ids = [str(i) for i in parent_ids]
            rule_id = str(rule_id)
            # Create the new ruleID by concatenating parent IDs with the original rule ID
            new_rule_id = (
                "p_" + "_".join(parent_ids) + "_r_" + rule_id
                if parent_ids
                else "r_" + rule_id
            )

            # Initialize a list to hold the modified lines
            modified_lines = []

            # Iterate through each line and replace the 'ruleID' line as needed
            for line in gml_content.splitlines():
                if line.strip().startswith("ruleID"):
                    # Replace the whole line with the new ruleID
                    modified_lines.append(f'\truleID "{new_rule_id}"')
                else:
                    modified_lines.append(line)

            # Join all lines back into a single string
            modified_content = "\n".join(modified_lines)

            # Write the modified content to the file
            with open(gml_file_path, "w") as file:
                file.write(modified_content)
            return True
        except FileNotFoundError:
            print(f"Unable to access the file path: {gml_file_path}")
            return False
        except Exception as e:
            print(f"An error occurred while writing to the file: {e}")
            return False
