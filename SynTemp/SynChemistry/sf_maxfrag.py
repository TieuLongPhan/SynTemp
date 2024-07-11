# from rdkit import Chem
# from typing import List, Tuple


# class SFMaxFrag:
#     """
#     A class to process a list of chemical reaction SMILES and sort them based on the size
#     of the largest non-Hydrogen atom fragment in the reactants and products.

#     Attributes:
#         reactions (List[str]): A list of reaction SMILES strings.
#     """

#     def __init__(self):
#         """
#         Initializes the SFMaxFrag object with a list of reaction SMILES strings.

#         Parameters:
#             reactions (List[str]): A list of chemical reaction SMILES strings.
#         """
#         pass

#     @staticmethod
#     def get_largest_fragment(smiles: str) -> str:
#         """
#         Extracts and returns the largest fragment by non-Hydrogen atom count from a SMILES string.

#         Parameters:
#             smiles (str): A SMILES string representing a molecule.

#         Returns:
#             str: The SMILES string of the largest fragment.
#         """
#         mol = Chem.MolFromSmiles(smiles)
#         frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
#         max_count = 0
#         largest_frag = None
#         for frag in frags:
#             count = sum(1 for atom in frag.GetAtoms() if atom.GetSymbol() != "H")
#             if count > max_count:
#                 max_count = count
#                 largest_frag = frag
#         return Chem.MolToSmiles(largest_frag)

#     @staticmethod
#     def count_non_hydrogen_atoms(smiles: str) -> int:
#         """
#         Counts the number of non-Hydrogen atoms in a molecule from its SMILES string.

#         Parameters:
#             smiles (str): The SMILES string of the molecule.

#         Returns:
#             int: The count of non-Hydrogen atoms.
#         """
#         mol = Chem.MolFromSmiles(smiles)
#         return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() != "H")

#     @staticmethod
#     def process_reaction(reaction: str) -> Tuple[str, int, int]:
#         """
#         Processes a single reaction SMILES string and returns the reaction along with the count of
#         non-Hydrogen atoms in the largest fragments of the reactants and products.

#         Parameters:
#             reaction (str): A reaction SMILES string formatted as 'reactants>>products'.

#         Returns:
#             Tuple[str, int, int]: A tuple containing the reaction SMILES string, the count of
#                                   non-Hydrogen atoms in the largest reactant fragment, and the
#                                   count in the largest product fragment.
#         """
#         reactants, products = reaction.split(">>")
#         largest_reactant = SFMaxFrag.get_largest_fragment(reactants)
#         largest_product = SFMaxFrag.get_largest_fragment(products)
#         reactant_count = SFMaxFrag.count_non_hydrogen_atoms(largest_reactant)
#         product_count = SFMaxFrag.count_non_hydrogen_atoms(largest_product)
#         return (reaction, reactant_count, product_count)

#     def sort_reactions(self, reactions) -> List[str]:
#         """
#         Sorts the stored reactions by the non-Hydrogen atom counts in the largest fragments of
#         the reactants and products, in descending order.

#         Returns:
#             List[str]: A list of sorted reaction SMILES strings.
#         """
#         processed_reactions = [SFMaxFrag.process_reaction(r) for r in reactions]
#         processed_reactions.sort(key=lambda x: (x[2], x[1]), reverse=True)
#         sorted_reactions = [item[0] for item in processed_reactions]
#         return sorted_reactions
