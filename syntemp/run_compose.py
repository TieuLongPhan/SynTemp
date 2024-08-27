import glob
from typing import List, Tuple
import argparse
from mod import ruleGMLString
from syntemp.SynComp.rule_compose import RuleCompose
from syntemp.SynUtils.graph_utils import load_gml_as_text


def load_compose_rules(compose_rule_path: str) -> List[Tuple[str, object]]:
    """
    Load and convert all GML files from a given directory to rule objects.

    Parameters:
    - compose_rule_path (str): Path to the directory containing GML files.

    Returns:
    - List[Tuple[str, object]]: A list of tuples, each containing the path
    and the corresponding rule object loaded from a GML file.
    """
    compose_files = glob.glob(f"{compose_rule_path}/*.gml")
    return [(path, ruleGMLString(load_gml_as_text(path))) for path in compose_files]


def detect_isomorphic_rules(
    double_rule_path: str, compose_rules: List[Tuple[str, object]]
) -> List[str]:
    """
    Detects and collects file paths that contain rules isomorphic to any rules in the p
    rovided set.

    Parameters:
    - double_rule_path (str): Path to the directory containing GML files to be checked.
    - compose_rules (List[Tuple[str, object]]): A list of tuples containing paths
    and rule objects against which to check for isomorphism.

    Returns:
        List[str]: A list of paths to GML files that contain isomorphic rules.
    """
    detected_files = []
    for double_file in glob.glob(f"{double_rule_path}/*.gml"):
        rule_double = ruleGMLString(load_gml_as_text(double_file))
        for _, rule in compose_rules:
            if rule_double.isomorphism(rule) == 1:
                detected_files.append(double_file)
                print(double_file)
                print("Composed_rule", rule)
    return detected_files


def main(args):
    """
    Main function to handle rule composition and detection of isomorphic rules.
    """
    if args.compose:
        rule_compose = RuleCompose()
        rule_compose._auto_compose(args.single_rule_path, args.compose_rule_path)

    compose_rules = load_compose_rules(args.compose_rule_path)
    detected_files = detect_isomorphic_rules(args.double_rule_path, compose_rules)
    print(f"Detected {len(detected_files)} isomorphic files.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and compare graph rules.")
    parser.add_argument(
        "-s",
        "--single_rule_path",
        type=str,
        required=True,
        help="Path to single expansion rules",
    )
    parser.add_argument(
        "-c",
        "--compose_rule_path",
        type=str,
        required=True,
        help="Path to composite rules",
    )
    parser.add_argument(
        "-d", "--double_rule_path", type=str, required=True, help="Path to double rules"
    )
    parser.add_argument(
        "-C",
        "--compose",
        action="store_true",
        help="Flag to trigger auto-composition of rules",
    )

    args = parser.parse_args()
    main(args)


# python run_compose.py -s Data/Temp/RuleComp/Single/R0
# -c Data/Temp/RuleComp/Compose -d Data/Temp/RuleComp/Double/R0/
