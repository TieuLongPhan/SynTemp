from typing import List, Any
from mod import DG, ruleGMLString, smiles, graphGMLString, addUniverse, addSubset, repeat
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
    def perform_reaction(rule_file_path: str, invert_rule: bool, initial_smiles: List[str], 
                     return_type: str = 'product', repeat_times: int = 1) -> List[str]:
        """
        Loads a reaction rule from a GML file, applies it to specified SMILES strings, and generates the resulting products or reaction SMILES.

        Parameters:
        - rule_file_path (str): Path to the GML file containing the reaction rule.
        - invert_rule (bool): Whether to invert the reaction rule. Useful for backward reactions.
        - initial_smiles (List[str]): List of initial molecules represented as SMILES strings.
        - return_type (str): Specifies the type of SMILES to return ('product' for product SMILES, 'reaction' for reaction SMILES).
        - repeat_times (int): Number of times to repeat the reaction rule. Defaults to 1.

        Returns:
        - List[str]: SMILES strings of the resulting molecules or reactions after applying the rule.
        """

        # Convert SMILES strings to molecule objects
        initial_molecules = [smiles(smile) for smile in initial_smiles]

        # Load the rule from the GML file
        gml_content = load_gml_as_text(rule_file_path)
        reaction_rule = ruleGMLString(gml_content, invert=invert_rule)

        if len (initial_molecules) > 1:
        # Define the strategy
            strategy = (addUniverse(initial_molecules[0]) >> addSubset(initial_molecules[1:]) >> repeat[repeat_times]([reaction_rule]))
        else:
            strategy = (addUniverse(initial_molecules[0]) >> addSubset(initial_molecules[0]) >> repeat[repeat_times]([reaction_rule]))

        # Initialize the derivation graph with the initial molecules
        dg = DG(graphDatabase=initial_molecules)
        dg.build().execute(strategy)

        # Collect the resulting products or reactions
        results = []
        if return_type == 'product':
            for graph in dg.products:
                product_smiles = MØDModeling.smilesFromProduct(graph)
                results.append(product_smiles)
                print(f"Product {graph.id}: {product_smiles}")
        elif return_type == 'reaction':
            reactant_smiles = '.'.join(initial_smiles)
            for graph in dg.products:
                product_smiles = MØDModeling.smilesFromProduct(graph)
                reaction_smiles = f"{reactant_smiles}>>{product_smiles}"
                results.append(reaction_smiles)
        return results