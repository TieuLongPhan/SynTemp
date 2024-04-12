from typing import List
from SynTemp.SynUtils.graph_utils import load_gml_as_text
from SynTemp.SynUtils.utils import ensure_directory_exists, run_shell_command
from mod import *


class RuleEngine:
    """
    The MØDModeling class encapsulates functionalities for reaction modeling using the MØD toolkit.
    It provides methods for forward and backward prediction based on templates library.
    """

    @staticmethod
    def generate_reaction_smiles(
        temp_results: List[str], base_smiles: str, is_forward: bool = True
    ) -> List[str]:
        """
        Constructs reaction SMILES strings from intermediate results using a base SMILES string, indicating whether the
        process is a forward or backward reaction. This function iterates over a list of intermediate SMILES strings,
        combines them with the base SMILES, and formats them into complete reaction SMILES strings.

        Parameters:
        - temp_results (List[str]): Intermediate SMILES strings resulting from partial reactions or combinations.
        - base_smiles (str): The SMILES string representing the starting point of the reaction, either as reactants
                            or products, depending on the reaction direction.
        - is_forward (bool, optional): Flag to determine the direction of the reaction; 'True' for forward reactions
                                    where 'base_smiles' are reactants, and 'False' for backward reactions where
                                    'base_smiles' are products. Defaults to True.

        Returns:
        - List[str]: A list of complete reaction SMILES strings, formatted according to the specified reaction direction.
        """
        results = []
        for comb in temp_results:
            joined_smiles = ".".join(comb)
            reaction_smiles = (
                f"{base_smiles}>>{joined_smiles}"
                if is_forward
                else f"{joined_smiles}>>{base_smiles}"
            )
            results.append(reaction_smiles)
        return results

    @staticmethod
    def perform_reaction(
        rule_file_path: str,
        initial_smiles: List[str],
        repeat_times: int = 1,
        prediction_type: str = "forward",
        print_results: bool = False,
    ) -> List[str]:
        """
        Applies a specified reaction rule, loaded from a GML file, to a set of initial molecules represented by SMILES strings.
        The reaction can be simulated in forward or backward direction and repeated multiple times.

        Parameters:
        - rule_file_path (str): Path to the GML file containing the reaction rule.
        - initial_smiles (List[str]): Initial molecules represented as SMILES strings.
        - repeat_times (int, optional): Number of iterations for the reaction. Defaults to 1.
        - type (str, optional): Direction of the reaction ('forward' for forward, 'backward' for backward). Defaults to 'forward'.
        - print_results (bool): Print results in latex or not. Defaults to False.

        Returns:
        - List[str]: SMILES strings of the resulting molecules or reactions.
        """

        # Determine the rule inversion based on reaction type
        invert_rule = prediction_type == "backward"

        # Convert SMILES strings to molecule objects, avoiding duplicate conversions
        initial_molecules = [smiles(smile) for smile in set(initial_smiles)]
        max_vertices = sum(molecule.numVertices for molecule in initial_molecules)

        # Load the reaction rule from the GML file
        gml_content = load_gml_as_text(rule_file_path)
        reaction_rule = ruleGMLString(gml_content, invert=invert_rule)

        # Define the reaction strategy
        universe = addUniverse(initial_molecules[0])
        subsets = (
            addSubset(initial_molecules[1:])
            if len(initial_molecules) > 1
            else addSubset(initial_molecules[0])
        )
        strategy = universe >> subsets >> repeat[repeat_times]([reaction_rule])

        # Initialize the derivation graph and execute the strategy
        dg = DG(graphDatabase=initial_molecules)
        dg.build().execute(strategy)
        if print_results:
            dg.print()

        temp_results = []
        for e in dg.edges:
            productSmiles = [v.graph.smiles for v in e.targets]
            temp_results.append(productSmiles)

        # Define a dictionary to map reaction types to their respective processing functions
        reaction_processing_map = {
            "forward": lambda smiles: RuleEngine.generate_reaction_smiles(
                temp_results, ".".join(smiles), is_forward=True
            ),
            "backward": lambda smiles: RuleEngine.generate_reaction_smiles(
                temp_results, ".".join(smiles), is_forward=False
            ),
        }

        # Use the reaction type to select the appropriate processing function and apply it
        if prediction_type in reaction_processing_map:
            return reaction_processing_map[prediction_type](initial_smiles)
        else:
            return ""
