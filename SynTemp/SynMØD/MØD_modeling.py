import glob
import copy
from typing import List, Any, Dict, Tuple
from SynTemp.SynMØD.MØD_postprocess import MØDPostprocess
from mod import DG, ruleGMLString, smiles, graphGMLString, addUniverse, addSubset, repeat, rightPredicate
from SynTemp.SynUtils.graph_utils import load_gml_as_text

class MØDModeling:
    """
    The MØDModeling class encapsulates functionalities for reaction modeling using the MØD toolkit.
    It provides methods for forward and backward prediction based on templates library.
    """
    @staticmethod
    def smilesFromProduct(product: Any) -> str:
        """
        Converts a product object from a reaction into a SMILES string representation.

        The product object, typically a graph representing a molecular structure, is traversed to construct a GML graph
        string. This GML string is then converted into a graph object, from which the SMILES representation is derived.

        Parameters:
        - product (Any): The product object to be converted into SMILES. The object structure is assumed to have vertices
                         and edges compatible with molecular graph representations.

        Returns:
        - str: The SMILES string representation of the product.
        """
        graphString = "graph [\n"
        for i, v in enumerate(product.vertices):
            graphString += f'  node [ id {i} label "{v.stringLabel}" ]\n'
            for e in v.incidentEdges:
                if v.id < e.target.id:
                    graphString += f'  edge [ source {v.id} target {e.target.id} label "{e.stringLabel}" ]\n'
        graphString += "]"
        graph = graphGMLString(graphString, name="Nan")  
        return graph.smiles
    
    @staticmethod
    def categorize_reactions(reactions: List[str], target_reaction: str) -> Tuple[List[str], List[str]]:
        """
        Categorizes the reactions into matches and not matches based on the target reaction.

        Args:
            reactions (List[str]): A list of reaction SMILES strings to categorize.
            target_reaction (str): The target reaction SMILES string to compare against.

        Returns:
            Tuple[List[str], List[str]]: A tuple containing two lists: matched reactions and non-matched reactions.
        """
        match, not_match = [], []
        target_reaction = MØDPostprocess.standardize_rsmi(target_reaction)
        for reaction_smiles in reactions:
            reaction_smiles = MØDPostprocess.standardize_rsmi(reaction_smiles)
            if reaction_smiles == target_reaction:
                match.append(reaction_smiles)
            else:
                not_match.append(reaction_smiles)
        return match, list(set(not_match))
    
    @staticmethod
    def reproduce_reactions(database: List[Dict], id_col: str, rule_file_path: str,
                            original_rsmi_col: str = 'reactions', 
                            repeat_times: int = 1, max_solutions: int = 10) -> List[Dict]:
        """
        Processes the database by performing reactions, categorizing them, and updating the database entries.

        Args:
            database (List[Dict]): The database to process, represented as a list of dictionaries.

        Returns:
            List[Dict]: The updated database after processing.
        """
        database_fw = copy.deepcopy(database)
        for entry in database_fw:
            rule_name = entry[id_col]
            rule_file = f'{rule_file_path}/{rule_name}.gml'
            initial_smiles_list = entry[original_rsmi_col].split('>>')[0].split('.')

            # Process reactions
            reactions = MØDModeling.perform_reaction(rule_file_path=rule_file, invert_rule=False, 
                                                    initial_smiles=initial_smiles_list, repeat_times=repeat_times, 
                                                    type='fw', max_solutions = max_solutions)

            # Categorize reactions
            match, not_match = MØDModeling.categorize_reactions(reactions, entry[original_rsmi_col])

            # Update database entry
            entry['positive_reactions'] = match[0] if match else None
            entry['negative_reactions'] = not_match

        database_bw = copy.deepcopy(database)
        for entry in database_bw:
            rule_name = entry[id_col]
            rule_file = f'{rule_file_path}/{rule_name}.gml'
            initial_smiles_list = entry[original_rsmi_col].split('>>')[1].split('.')

            # Process reactions
            reactions = MØDModeling.perform_reaction(rule_file_path=rule_file, invert_rule=True, 
                                                    initial_smiles=initial_smiles_list, repeat_times=repeat_times, 
                                                    type='bw', max_solutions = max_solutions)

            # Categorize reactions
            match, not_match = MØDModeling.categorize_reactions(reactions, entry[original_rsmi_col])

            # Update database entry
            entry['positive_reactions'] = match[0] if match else None
            entry['negative_reactions'] = not_match

        return database_fw, database_bw
    
    @staticmethod
    def forward_prediction(database: List[Dict], rule_file_path: str,
                        original_rsmi_col: str = 'reactions', 
                        repeat_times: int = 1, max_solutions: int = 10) -> List[Dict]:
        """
        Processes the database by performing forward reaction predictions for each entry, using the specified rule files.
        The function creates a copy of the input database and adds the forward predictions to each entry in the copy, 
        thereby not modifying the input database in-place.

        Args:
            database (List[Dict]): The database to process, represented as a list of dictionaries.
            rule_file_path (str): The file path where the rule files are stored.
            original_rsmi_col (str, optional): The column name in the database that contains the original reaction SMILES. Defaults to 'reactions'.
            repeat_times (int, optional): The number of times to repeat the reaction process. Defaults to 1.

        Returns:
            List[Dict]: A new database list, with each entry updated to include forward reaction predictions.
        """
        # Create a deep copy of the database to avoid modifying the input in-place
        database_copy = copy.deepcopy(database)

        for entry in database_copy:
            if '>>' in entry.get(original_rsmi_col, ''):
                reactants_list = entry.get(original_rsmi_col, '').split('>>')[0].split('.')
            else:
                reactants_list = entry.get(original_rsmi_col, '').split('.')
            predictions = []
            for rule_file in glob.glob(f'{rule_file_path}/*.gml'):
                # Process reactions using each rule file
                reactions = MØDModeling.perform_reaction(rule_file_path=rule_file, invert_rule=False, 
                                                        initial_smiles=reactants_list, repeat_times=repeat_times,
                                                        max_solutions = max_solutions)
                predictions.extend(reactions)
            
            # Update the entry with aggregated predictions from all rule files
            predictions = [MØDPostprocess.standardize_rsmi(value) for value in predictions]
            predictions = list(set(predictions))
            entry['forward_predictions'] = predictions
            entry['number_predictions'] = len(predictions)

        return database_copy

    @staticmethod
    def backward_prediction(database: List[Dict], rule_file_path: str,
                        original_rsmi_col: str = 'reactions', 
                        repeat_times: int = 1, max_solutions: int = 10) -> List[Dict]:
        """
        Processes the database by performing forward reaction predictions for each entry, using the specified rule files.
        The function creates a copy of the input database and adds the forward predictions to each entry in the copy, 
        thereby not modifying the input database in-place.

        Args:
            database (List[Dict]): The database to process, represented as a list of dictionaries.
            rule_file_path (str): The file path where the rule files are stored.
            original_rsmi_col (str, optional): The column name in the database that contains the original reaction SMILES. Defaults to 'reactions'.
            repeat_times (int, optional): The number of times to repeat the reaction process. Defaults to 1.

        Returns:
            List[Dict]: A new database list, with each entry updated to include forward reaction predictions.
        """
        # Create a deep copy of the database to avoid modifying the input in-place
        database_copy = copy.deepcopy(database)

        for entry in database_copy:
            if '>>' in entry.get(original_rsmi_col, ''):
                products_list = entry.get(original_rsmi_col, '').split('>>')[1].split('.')
            else:
                products_list = entry.get(original_rsmi_col, '').split('.')

            predictions = []
            for rule_file in glob.glob(f'{rule_file_path}/*.gml'):
                # Process reactions using each rule file
                reactions = MØDModeling.perform_reaction(rule_file_path=rule_file, invert_rule=True, 
                                                        initial_smiles=products_list, repeat_times=repeat_times, 
                                                        type ='bw', max_solutions=max_solutions)
                predictions.extend(reactions)
            predictions = [MØDPostprocess.standardize_rsmi(value) for value in predictions]
            predictions = list(set(predictions))
            # Update the entry with aggregated predictions from all rule files
            entry['backward_predictions'] = predictions
            entry['number_predictions'] = len(predictions)


        return database_copy
    
    @staticmethod
    def generate_reaction_smiles(temp_results: List[str], base_smiles: str, is_forward: bool = True) -> List[str]:
        """
        Generate reaction SMILES strings based on the temporary results, given the base SMILES string.

        Parameters:
            temp_results (List[str]): List of temporary result SMILES strings.
            base_smiles (str): Base SMILES string representing the reactants or products.
            is_forward (bool, optional): Indicates whether the reaction is forward (True) or backward (False). Defaults to True.

        Returns:
            List[str]: List of reaction SMILES strings.
        """
        results = []
        for comb in MØDPostprocess.generate_smiles_combinations(temp_results, base_smiles, True):
            joined_smiles = '.'.join(comb)
            reaction_smiles = f"{base_smiles}>>{joined_smiles}" if is_forward else f"{joined_smiles}>>{base_smiles}"
            results.append(reaction_smiles)
        return results


    @staticmethod
    def perform_reaction(rule_file_path: str, invert_rule: bool, initial_smiles: List[str], 
                     repeat_times: int = 1, type: str = 'fw', max_solutions: int = 10) -> List[str]:
        """
        Loads a reaction rule from a GML file, applies it to specified SMILES strings, and generates the resulting products or reaction SMILES.

        Parameters:
        - rule_file_path (str): Path to the GML file containing the reaction rule.
        - invert_rule (bool): Whether to invert the reaction rule. Useful for backward reactions.
        - initial_smiles (List[str]): List of initial molecules represented as SMILES strings.
        - repeat_times (int): Number of times to repeat the reaction rule. Defaults to 1.
        - type (str): Types of prediction: forward (fw) or backward (bw). Defaults to 'fw'.
        - max_solutions (int): maximum number of solutions 

        Returns:
        - List[str]: SMILES strings of the resulting molecules or reactions after applying the rule.
        """

        # Convert SMILES strings to molecule objects
        initial_molecules = []
        max_vertice = 0
        for smile in initial_smiles:
            initial_molecules.append(smiles(smile))
            max_vertice += smiles(smile).numVertices
            #print(smiles(smile).numVertices)
        # Load the rule from the GML file
        gml_content = load_gml_as_text(rule_file_path)
        reaction_rule = ruleGMLString(gml_content, invert=invert_rule)

        if len (initial_molecules) > 1:
        #Define the strategy
            # strategy = (addUniverse(initial_molecules[0]) 
            #             >> addSubset(initial_molecules[1:]) 
            #             >> rightPredicate[
            #                 lambda derivation: all(g.numVertices <= max_vertice for g in derivation.right)
            #             ]
            #             >> repeat[repeat_times]([reaction_rule]))
            strategy = (addUniverse(initial_molecules[0]) >> addSubset(initial_molecules[1:]) >> repeat[repeat_times]([reaction_rule]))
        else:
            strategy = (addUniverse(initial_molecules[0]) >> addSubset(initial_molecules[0]) >> repeat[repeat_times]([reaction_rule]))

        # Initialize the derivation graph with the initial molecules
        dg = DG(graphDatabase=initial_molecules)
        dg.build().execute(strategy)


        # Collect resulting products or reactions
        temp_results = [MØDModeling.smilesFromProduct(graph) for graph in dg.products]
        if type == 'fw':
            reactant_smiles = '.'.join(initial_smiles)
            if len(temp_results) <= max_solutions:
                return MØDModeling.generate_reaction_smiles(temp_results, reactant_smiles)
            else:
                # Handle cases with more than 10 temporary results separately
                return [f"{reactant_smiles}>>{smiles}" for smiles in temp_results
                        if MØDPostprocess.get_combined_molecular_formula([smiles]) ==
                        MØDPostprocess.get_combined_molecular_formula([reactant_smiles])]
        elif type == 'bw':
            product_smiles = '.'.join(initial_smiles)
            if len(temp_results) <= max_solutions:
                return MØDModeling.generate_reaction_smiles(temp_results, product_smiles, is_forward=False)
            else:
                # Handle cases with more than 10 temporary results separately
                return [f"{smiles}>>{product_smiles}" for smiles in temp_results
                        if MØDPostprocess.get_combined_molecular_formula([smiles]) ==
                        MØDPostprocess.get_combined_molecular_formula([product_smiles])]

        return []  # Return an empty list for unsupported reaction types
