from rdkit import Chem
import networkx as nx
from typing import Dict


class GraphToMol:
    """
    Converts a NetworkX graph representation of a molecule back into an RDKit molecule
    object, taking into account specific node and edge attributes for the construction of
    the molecule. An option to ignore bond orders can be specified,
    treating all bonds as single bonds.
    """

    def __init__(
        self, node_attributes: Dict[str, str], edge_attributes: Dict[str, str]
    ):
        """
        Initializes the GraphToMol converter with mappings for node and edge attributes
        to their corresponding RDKit atom and bond properties.
        """
        self.node_attributes = node_attributes
        self.edge_attributes = edge_attributes

    def graph_to_mol(
        self, graph: nx.Graph, ignore_bond_order: bool = False
    ) -> Chem.Mol:
        """
        Converts a NetworkX graph into an RDKit molecule object by interpreting node and
        edge attributes according to the provided mappings. Optionally ignores bond
        orders, treating all bonds as single bonds.

        Parameters:
        - graph (nx.Graph): The NetworkX graph representation of the molecule to be
        converted.
        - ignore_bond_order (bool): If True, all bonds are treated as single bonds
        regardless of their specified order in the graph.

        Returns:
        - Chem.Mol: An RDKit molecule object constructed based on the graph
        representation.
        """
        mol = Chem.RWMol()

        # Map for tracking RDKit atom indices corresponding to NetworkX nodes
        node_to_idx: Dict[int, int] = {}

        # Add atoms to the molecule based on node attributes
        for node, data in graph.nodes(data=True):
            element = data.get(
                self.node_attributes["element"], "C"
            )  # Defaults to Carbon
            charge = data.get(self.node_attributes["charge"], 0)
            atom = Chem.Atom(element)
            atom.SetFormalCharge(charge)
            if "atom_class" in data:  # Set atom map number if available
                atom.SetAtomMapNum(data["atom_class"])
            idx = mol.AddAtom(atom)
            node_to_idx[node] = idx

        # Add bonds between atoms based on edge attributes
        for u, v, data in graph.edges(data=True):
            if ignore_bond_order:
                bond_order = 1  # Treat all bonds as single bonds
            else:
                bond_order = abs(
                    data.get(self.edge_attributes["order"], 1)
                )  # Use absolute value of bond order
            bond_type = self.get_bond_type_from_order(bond_order)
            mol.AddBond(node_to_idx[u], node_to_idx[v], bond_type)

        # Attempt to sanitize the molecule to ensure its chemical validity
        Chem.SanitizeMol(mol)

        return mol

    @staticmethod
    def get_bond_type_from_order(order: float) -> Chem.BondType:
        """
        Converts a numerical bond order into the corresponding RDKit BondType enum.
        """
        if order == 1:
            return Chem.BondType.SINGLE
        elif order == 2:
            return Chem.BondType.DOUBLE
        elif order == 3:
            return Chem.BondType.TRIPLE
        return Chem.BondType.AROMATIC
