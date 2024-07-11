from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
from typing import Dict
import random


class MolToGraph:
    """
    A class for converting molecules from SMILES strings to graph representations using
    RDKit and NetworkX.
    """

    def __init__(self):
        """Initialize the MolToGraphConverter class."""
        pass

    @staticmethod
    def add_partial_charges(mol: Chem.Mol) -> None:
        """
        Computes and assigns Gasteiger partial charges to each atom in the given molecule.

        Parameters:
        - mol (Chem.Mol): An RDKit molecule object.
        """
        AllChem.ComputeGasteigerCharges(mol)

    @staticmethod
    def get_stereochemistry(atom: Chem.Atom) -> str:
        """
        Determines the stereochemistry (R/S configuration) of a given atom.

        Parameters:
        - atom (Chem.Atom): An RDKit atom object.

        Returns:
        - str: The stereochemistry ('R', 'S', or 'N' for non-chiral).
        """
        chiral_tag = atom.GetChiralTag()
        if chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return "S"
        elif chiral_tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return "R"
        return "N"

    @staticmethod
    def get_bond_stereochemistry(bond: Chem.Bond) -> str:
        """
        Determines the stereochemistry (E/Z configuration) of a given bond.

        Parameters:
        - bond (Chem.Bond): An RDKit bond object.

        Returns:
        - str: The stereochemistry ('E', 'Z', or 'N' for non-stereospecific
        or non-double bonds).
        """
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            return "N"
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            return "E"
        elif stereo == Chem.BondStereo.STEREOZ:
            return "Z"
        return "N"

    @staticmethod
    def has_atom_mapping(mol: Chem.Mol) -> bool:
        """
        Check if the given molecule has any atom mapping numbers.

        Atom mapping numbers are used in chemical reactions to track the correspondence
        between atoms in reactants and products.

        Parameters:
        - mol (Chem.Mol): An RDKit molecule object.

        Returns:
        - bool: True if any atom in the molecule has a mapping number, False otherwise.
        """
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                return True
        return False

    @staticmethod
    def random_atom_mapping(mol: Chem.Mol) -> Chem.Mol:
        """
        Assigns a random atom mapping number to each atom in the given molecule.

        This method iterates over all atoms in the molecule and assigns a random
        mapping number between 1 and the total number of atoms to each atom.

        Parameters:
        - mol (Chem.Mol): An RDKit molecule object.

        Returns:
        - Chem.Mol: The RDKit molecule object with random atom mapping numbers assigned.
        """
        atom_indices = list(range(1, mol.GetNumAtoms() + 1))
        random.shuffle(atom_indices)
        for atom, idx in zip(mol.GetAtoms(), atom_indices):
            atom.SetProp("molAtomMapNumber", str(idx))
        return mol

    @classmethod
    def mol_to_graph(cls, mol: Chem.Mol, drop_non_aam: bool = False) -> nx.Graph:
        """
        Converts an RDKit molecule object to a NetworkX graph with specified atom and bond
        attributes.
        Optionally excludes atoms without atom mapping numbers if drop_non_aam is True.

        Parameters:
        - mol (Chem.Mol): An RDKit molecule object.
        - drop_non_aam (bool, optional): If True, nodes without atom mapping numbers will
        be dropped.

        Returns:
        - nx.Graph: A NetworkX graph representing the molecule.
        """
        cls.add_partial_charges(mol)
        graph = nx.Graph()
        index_to_class: Dict[int, int] = {}
        if cls.has_atom_mapping(mol) is False:
            mol = cls.random_atom_mapping(mol)

        for atom in mol.GetAtoms():
            atom_map = atom.GetAtomMapNum()

            if drop_non_aam and atom_map == 0:
                continue
            gasteiger_charge = round(float(atom.GetProp("_GasteigerCharge")), 3)
            index_to_class[atom.GetIdx()] = atom_map
            graph.add_node(
                atom_map,
                charge=atom.GetFormalCharge(),
                hcount=atom.GetTotalNumHs(),
                aromatic=atom.GetIsAromatic(),
                element=atom.GetSymbol(),
                atom_map=atom_map,
                isomer=cls.get_stereochemistry(atom),
                partial_charge=gasteiger_charge,
                hybridization=str(atom.GetHybridization()),
                in_ring=atom.IsInRing(),
                explicit_valence=atom.GetExplicitValence(),
                implicit_hcount=atom.GetNumImplicitHs(),
                neighbors=sorted(
                    [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
                ),
            )

        for bond in mol.GetBonds():
            begin_atom_class = index_to_class.get(bond.GetBeginAtomIdx())
            end_atom_class = index_to_class.get(bond.GetEndAtomIdx())
            if begin_atom_class is None or end_atom_class is None:
                continue
            graph.add_edge(
                begin_atom_class,
                end_atom_class,
                order=bond.GetBondTypeAsDouble(),
                ez_isomer=cls.get_bond_stereochemistry(bond),
                bond_type=str(bond.GetBondType()),
                conjugated=bond.GetIsConjugated(),
                in_ring=bond.IsInRing(),
            )

        return graph
