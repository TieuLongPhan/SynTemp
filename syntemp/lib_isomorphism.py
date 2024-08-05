import os
from syntemp.SynUtils.graph_utils import load_gml_as_text
from mod import ruleGMLString


class LibIsomorphism:

    @staticmethod
    def rule_isomorphism(rule_1: str, rule_2: str) -> bool:
        """
        Determines if two GML rules are isomorphic. This function converts two rule
        representations from string format to a GML rule object and then checks for
        isomorphism between them.

        Parameters:
        rule_1 (str): The GML string representation of the first rule.
        rule_2 (str): The GML string representation of the second rule.

        Returns:
        bool: True if the rules are isomorphic, False otherwise.
        """
        rule_1 = ruleGMLString(rule_1)
        rule_2 = ruleGMLString(rule_2)
        isomorphic = rule_1.isomorphism(rule_2) == 1
        return isomorphic

    @staticmethod
    def lib_isomorphism(rule: str, lib_path: str) -> bool:
        """
        Checks if a given GML rule string is isomorphic to any rule in a specified library
        directory. This function searches for GML files in the given directory and checks
        each for isomorphism with the provided rule.

        Paremeters:
        - rule (str): The GML string representation of the rule to check.
        - lib_path (str): The directory path containing GML files to check against.

        Returns:
        - bool: True if any file in the library is isomorphic to the given rule, False
        otherwise.
        """
        # Construct full file path list for GML files in the specified library path
        gml_files = [
            os.path.join(lib_path, f)
            for f in os.listdir(lib_path)
            if f.endswith(".gml")
        ]

        # Check each file for isomorphism
        for file_path in gml_files:
            file_content = load_gml_as_text(file_path)
            if LibIsomorphism.rule_isomorphism(rule, file_content):
                return True
        return False
