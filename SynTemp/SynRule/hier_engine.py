from typing import List, Dict, Any
import logging
from SynTemp.SynUtils.graph_utils import load_gml_as_text
from SynTemp.SynUtils.chemutils import generate_reaction_smiles
import glob
from mod import *

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


class HierEngine:
    @staticmethod
    def rule_apply(initial_molecules, reaction_rule):
        """
        Apply a chemical reaction rule to a set of initial molecules and return the resulting product SMILES strings.

        Args:
            initial_molecules (List[Any]): A list of molecule objects to which the reaction rule is applied.
            reaction_rule (Any): The reaction rule object, typically parsed from a GML file.

        Returns:
            List[List[str]]: A nested list where each sublist contains the SMILES strings of the reaction products
            for one molecule. Returns an empty list if no products are formed or an error occurs.
        """
        try:
            # Initialize the DG with the given initial molecules
            dg = DG(graphDatabase=initial_molecules)
            config.dg.doRuleIsomorphismDuringBinding = False
            dg.build().apply(initial_molecules, reaction_rule)
            products = []
            for e in dg.edges:
                productSmiles = [v.graph.smiles for v in e.targets]
                products.append(productSmiles)
            return products if products else []
        except Exception as e:
            print(f"An error occurred: {e}")
            return []

    @staticmethod
    def hier_child_level(
        initial_molecules: List[Any],
        rule_file_path: str,
        hier_temp: Dict[int, List[Dict[str, Any]]],
        radius: int,
        rule_id: int = 0,
        max_radius: int = 3,
        invert_rule: bool = False,
        prune: bool = True,
        prune_size: int = 1,
    ) -> List[Any]:
        """
        Recursively apply hierarchical chemical reaction rules starting from a given reaction radius level.
        Args:
            initial_molecules: List of initial molecule structures.
            rule_file_path: Path to the folder containing hierarchical reaction rules.
            hier_temp: Hierarchical template dictating the application of rules.
            radius: Current radius level for rule application.
            rule_id: Identifier for the specific rule to apply (optional).
            max_radius: Maximum radius to apply rules up to.
            invert_rule: Whether to invert the rule during application, depending on the reaction direction.
            prune: Option to stop recursion and return early if fewer than `prune_size` results are produced.
            prune_id: Threshold number of results below which pruning happens.
        Returns:
            A flattened list of unique products from the applied hierarchical rules.
        """
        try:
            rule_path = f"{rule_file_path}/R{radius}/{rule_id}.gml"
            gml_content = load_gml_as_text(rule_path)
            reaction_rule = ruleGMLString(gml_content, invert=invert_rule, add=False)

            temp_results = HierEngine.rule_apply(initial_molecules, reaction_rule)

            # if radius > max_radius - 1:
            #     return temp_results

            # if prune and len(temp_results) < prune_size:
            #     return temp_results

            if radius > max_radius - 1 or (prune and len(temp_results) < prune_size):
                return temp_results

            new_rule_ids = [
                entry["Child"]
                for entry in hier_temp[radius]
                if entry["Cluster_id"] == rule_id
            ][0]
            results = []
            for entry in new_rule_ids:
                result = HierEngine.hier_child_level(
                    initial_molecules,
                    rule_file_path,
                    hier_temp,
                    radius + 1,
                    entry,
                    max_radius,
                    invert_rule,
                )
                if result:
                    results.extend(result)

            return results

        except FileNotFoundError:
            logging.error(f"File not found: {rule_path}")
            return []
        except Exception as e:
            logging.error(f"Error in processing at hier_child_level: {e}")
            return []

    @staticmethod
    def hier_rule_apply(
        initial_smiles: List[Any],
        hier_temp: Dict[int, List[Dict[str, Any]]],
        rule_file_path: str,
        prediction_type: str = "forward",
        max_radius: int = 3,
        max_solutions: int = 1000,
        prune: bool = True,
        prune_size: int = 1,
        templates_threshold: int = 0.00
    ) -> List[List[Any]]:
        """
        Apply hierarchical chemical reaction rules to a dataset of molecules based on their SMILES strings.
        Args:
            initial_smiles: List of SMILES strings representing the molecules to be processed.
            hier_temp: Hierarchical template for applying rules.
            rule_file_path: Path to the directory containing rule files.
            prediction_type: Type of prediction, "forward" for forward reactions or "backward" for retrosynthesis.
            max_radius: Maximum radial depth for applying the rules.
            max_solutions: Maximum number of solutions to consider from one rule application.
            prune: Option to stop recursion and return early based on the number of results.
            prune_size: Threshold number of results below which pruning happens.
        Returns:
            A list containing the results from applying the hierarchical rules to all entries, processed
            according to the reaction type.
        """

        initial_molecules = sorted(
            (smiles(smile, add=False) for smile in initial_smiles),
            key=lambda molecule: molecule.numVertices,
        )
        invert_rule = prediction_type == "backward"

        hier_0 = [value for value in hier_temp[0] if value['Percentage'] >= templates_threshold]
        hier_0_id = [value['Cluster_id'] for value in hier_0]

        # Create a list of file patterns to search for, using each Cluster_id
        file_patterns = [f"{rule_file_path}/R0/{a}.gml" for a in hier_0_id]

        # Use glob.glob to find all files that match the patterns in the list
        rule_files = []
        for pattern in file_patterns:
            rule_files.extend(glob.glob(pattern))

        #rule_files = glob.glob(f"{rule_file_path}/R0/*.gml")
        temp_results = []
        for rule_file in rule_files:
            result = HierEngine.hier_child_level(
                initial_molecules,
                rule_file_path,
                hier_temp,
                radius=0,
                rule_id=int(rule_file.split("/")[-1].replace(".gml", "")),
                max_radius=max_radius,
                invert_rule=invert_rule,
                prune=prune,
                prune_size=prune_size,
            )
            if len(result) < max_solutions:
                temp_results.extend(result)
            else:
                temp_results.extend(result[:max_solutions])

        reaction_processing_map = {
            "forward": lambda smiles: generate_reaction_smiles(
                temp_results, ".".join(smiles), is_forward=True
            ),
            "backward": lambda smiles: generate_reaction_smiles(
                temp_results, ".".join(smiles), is_forward=False
            ),
        }

        return reaction_processing_map.get(prediction_type, lambda x: [])(
            initial_smiles
        )
