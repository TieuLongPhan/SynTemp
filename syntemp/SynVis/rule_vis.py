import os
import glob
from syntemp.SynUtils.graph_utils import load_gml_as_text
from mod import ruleGMLString, inputRules


def rule_vis(rule_path, save_path=None, verbose=False):
    """
    Visualizes a reaction rule from a GML file.

    Parameters:
    - rule_path (str): Path to the GML file containing the reaction rule.
    - save_path (str, optional): Directory to save the visualized rule outputs.
                                 Defaults to 'out' if None provided.
    - verbose (bool, optional): Flag to enable detailed logging of the process.

    Returns:
    - None: Outputs are saved directly to files in the specified directory.
    """
    try:
        # Load the rule from the GML file
        rule = load_gml_as_text(rule_path)
        if verbose:
            print(f"Loaded rule from {rule_path}")

        # Ensure the save directory exists
        output_dir = save_path if save_path else "out"
        os.makedirs(output_dir, exist_ok=True)
        if verbose:
            print(f"Output directory set to {output_dir}")

        # Process the rule using custom logic from 'mod'
        ruleGMLString(rule)
        if verbose:
            print("Processed the rule using ruleGMLString.")

        # Iterate over rules, possibly meant for further processing or output
        for rule in inputRules:
            # This may involve saving or further manipulating rules
            rule.print(second=True)

    except Exception as e:
        print(f"An error occurred: {e}")
        if verbose:
            raise


def rules_vis(rule_paths, save_path=None, verbose=False):
    # Ensure the save directory exists
    output_dir = save_path if save_path else "out"
    os.makedirs(output_dir, exist_ok=True)
    if verbose:
        print(f"Output directory set to {output_dir}")
    for rule_path in rule_paths:
        rule = load_gml_as_text(rule_path)
        ruleGMLString(rule)
    for rule in inputRules:
        # This may involve saving or further manipulating rules
        rule.print(second=True)


def auto_rules_vis(path_to_gml, save_path=None, verbose=False):
    # Ensure the save directory exists
    rule_paths = [os.path.basename(f) for f in glob.glob(f"{path_to_gml}/*.gml")]
    output_dir = save_path if save_path else "out"
    os.makedirs(output_dir, exist_ok=True)
    if verbose:
        print(f"Output directory set to {output_dir}")
    for rule_path in rule_paths:
        rule = load_gml_as_text(f"{path_to_gml}/{rule_path}")
        ruleGMLString(rule)
    for rule in inputRules:
        # This may involve saving or further manipulating rules
        rule.print(second=True)
