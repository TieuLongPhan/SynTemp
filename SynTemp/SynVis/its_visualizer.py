from collections import namedtuple
from typing import List, Tuple, Dict
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions


class ITSVisualizer:
    """
    Class for visualizing and analyzing chemical reactions using RDKit.
    """

    AtomInfo = namedtuple(
        "AtomInfo", ("mapnum", "reactant", "reactantAtom", "product", "productAtom")
    )
    BondInfo = namedtuple(
        "BondInfo", ("product", "productAtoms", "productBond", "status")
    )

    def __init__(self, reaction_smiles: str) -> None:
        """
        Initialize the ITSVisualizer object with a reaction SMILES string.
        """
        self.reaction_smiles: str = reaction_smiles
        self.rxn: rdChemReactions.ChemicalReaction = rdChemReactions.ReactionFromSmarts(
            self.reaction_smiles, useSmiles=True
        )
        self.rxn.Initialize()
        self.atms: List[ITSVisualizer.AtomInfo]
        self.bnds: List[ITSVisualizer.BondInfo]
        self.atms, self.bnds = self.find_modifications_in_products(self.rxn)

    @staticmethod
    def map_reacting_atoms_to_products(
        rxn: rdChemReactions.ChemicalReaction, reactingAtoms: List[List[int]]
    ) -> List[AtomInfo]:
        """
        Maps reacting atoms to their corresponding products.
        """
        res: List[ITSVisualizer.AtomInfo] = []
        for ridx, reacting in enumerate(reactingAtoms):
            reactant = rxn.GetReactantTemplate(ridx)
            for raidx in reacting:
                mapnum = reactant.GetAtomWithIdx(raidx).GetAtomMapNum()
                for pidx, product in enumerate(rxn.GetProducts()):
                    for paidx, patom in enumerate(product.GetAtoms()):
                        if patom.GetAtomMapNum() == mapnum:
                            res.append(
                                ITSVisualizer.AtomInfo(mapnum, ridx, raidx, pidx, paidx)
                            )
                            break
        return res

    @staticmethod
    def get_mapped_neighbors(atom: Chem.Atom) -> Dict[Tuple[int, int], Tuple[int, int]]:
        """
        Finds all mapped neighbors of a mapped atom.
        """
        res: Dict[Tuple[int, int], Tuple[int, int]] = {}
        amap = atom.GetAtomMapNum()
        for nbr in atom.GetNeighbors():
            nmap = nbr.GetAtomMapNum()
            if nmap:
                key = tuple(sorted((amap, nmap)))
                res[key] = (atom.GetIdx(), nbr.GetIdx())
        return res

    @staticmethod
    def find_modifications_in_products(
        rxn: rdChemReactions.ChemicalReaction,
    ) -> Tuple[List[AtomInfo], List[BondInfo]]:
        """
        Identifies modifications in products compared to reactants.
        """
        reactingAtoms = rxn.GetReactingAtoms()
        amap = ITSVisualizer.map_reacting_atoms_to_products(rxn, reactingAtoms)
        res: List[ITSVisualizer.AtomInfo] = []
        seen: set = set()
        for _, ridx, raidx, pidx, paidx in amap:
            reactant = rxn.GetReactantTemplate(ridx)
            ratom = reactant.GetAtomWithIdx(raidx)
            rnbrs = ITSVisualizer.get_mapped_neighbors(ratom)
            product = rxn.GetProductTemplate(pidx)
            patom = product.GetAtomWithIdx(paidx)
            pnbrs = ITSVisualizer.get_mapped_neighbors(patom)
            for tpl in pnbrs:
                if tpl not in rnbrs:
                    pbond = product.GetBondBetweenAtoms(*pnbrs[tpl])
                    if (pidx, pbond.GetIdx()) not in seen:
                        seen.add((pidx, pbond.GetIdx()))
                        res.append(
                            ITSVisualizer.BondInfo(
                                pidx, pnbrs[tpl], pbond.GetIdx(), "New"
                            )
                        )
                else:
                    pbond = product.GetBondBetweenAtoms(*pnbrs[tpl])
                    # present in both reactants and products, check to see if it changed
                    rbond = reactant.GetBondBetweenAtoms(*rnbrs[tpl])
                    if rbond.GetBondType() != pbond.GetBondType():
                        res.append(
                            ITSVisualizer.BondInfo(
                                pidx, pnbrs[tpl], pbond.GetIdx(), "Changed"
                            )
                        )
        return amap, res

    def draw_product_with_modified_bonds(
        self, productIdx: int = None, showAtomMaps: bool = False, img_size=(1200, 500)
    ) -> str:
        """
        Draws the product molecule highlighting the modified bonds and atoms.
        """
        rxn, atms, bnds = self.rxn, self.atms, self.bnds
        if productIdx is None:
            pcnts = [x.GetNumAtoms() for x in rxn.GetProducts()]
            largestProduct = list(sorted(zip(pcnts, range(len(pcnts))), reverse=True))[
                0
            ][1]
            productIdx = largestProduct
        d2d = Draw.rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])

        pmol = Chem.Mol(rxn.GetProductTemplate(productIdx))
        Chem.SanitizeMol(pmol)
        if not showAtomMaps:
            for atom in pmol.GetAtoms():
                atom.SetAtomMapNum(0)
        bonds_to_highlight = []
        highlight_bond_colors = {}
        atoms_seen = set()
        for binfo in bnds:
            if binfo.product == productIdx and binfo.status == "New":
                bonds_to_highlight.append(binfo.productBond)
                atoms_seen.update(binfo.productAtoms)
                highlight_bond_colors[binfo.productBond] = (0.4, 0.4, 1)
            if binfo.product == productIdx and binfo.status == "Changed":
                bonds_to_highlight.append(binfo.productBond)
                atoms_seen.update(binfo.productAtoms)
                highlight_bond_colors[binfo.productBond] = (1, 0.4, 0.4)
        atoms_to_highlight = set()
        for ainfo in atms:
            if ainfo.product != productIdx or ainfo.productAtom in atoms_seen:
                continue
            atoms_to_highlight.add(ainfo.productAtom)

        d2d.drawOptions().useBWAtomPalette()
        d2d.drawOptions().continuousHighlight = False
        d2d.drawOptions().highlightBondWidthMultiplier = 24
        d2d.drawOptions().setHighlightColour((0.9, 0.9, 0))
        d2d.drawOptions().fillHighlights = False
        d2d.drawOptions().bondLineWidth = 3
        d2d.drawOptions().scale = 1
        d2d.drawOptions().atomLabelFontSize = 24
        d2d.drawOptions().padding = 0.1

        atoms_to_highlight.update(atoms_seen)
        d2d.DrawMolecule(
            pmol,
            highlightAtoms=atoms_to_highlight,
            highlightBonds=bonds_to_highlight,
            highlightBondColors=highlight_bond_colors,
        )
        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def draw_highlighted_subgraph(
        self, productIdx: int = None, showAtomMaps: bool = False, img_size=(1200, 500)
    ) -> str:
        """
        Draws only the subgraph of the product molecule containing modified
        atoms and bonds, with specified colors for highlighted bonds.
        """
        rxn, _, bnds = self.rxn, self.atms, self.bnds
        if productIdx is None:
            pcnts = [x.GetNumAtoms() for x in rxn.GetProducts()]
            largestProduct = max(enumerate(pcnts), key=lambda x: x[1])[0]
            productIdx = largestProduct

        d2d = Draw.rdMolDraw2D.MolDraw2DCairo(
            img_size[0], img_size[1]
        )  # Improved resolution and anti-aliasing

        pmol = Chem.Mol(rxn.GetProductTemplate(productIdx))
        Chem.SanitizeMol(pmol)
        if not showAtomMaps:
            for atom in pmol.GetAtoms():
                atom.SetAtomMapNum(0)

        # Identifying the atoms and bonds to highlight
        atoms_to_highlight = set()
        bonds_to_highlight = {}
        for binfo in bnds:
            if binfo.product == productIdx:
                for atom in binfo.productAtoms:
                    atoms_to_highlight.add(atom)
                if binfo.status == "New":
                    bonds_to_highlight[binfo.productBond] = (0.4, 0.4, 1)
                elif binfo.status == "Changed":
                    bonds_to_highlight[binfo.productBond] = (1, 0.4, 0.4)

        # Creating a subgraph with the highlighted atoms and bonds
        submol = Chem.PathToSubmol(pmol, list(bonds_to_highlight.keys()))
        atom_map = {
            atom_idx: submol.GetAtomWithIdx(i).GetIdx()
            for i, atom_idx in enumerate(atoms_to_highlight)
        }

        highlight_atom_ids = [
            atom_map[atom_idx]
            for atom_idx in atoms_to_highlight
            if atom_idx in atom_map
        ]
        highlight_bond_colors = {}
        for bond_idx, color in bonds_to_highlight.items():
            bond = pmol.GetBondWithIdx(bond_idx)
            if bond:
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if begin_idx in atom_map and end_idx in atom_map:
                    sub_bond = submol.GetBondBetweenAtoms(
                        atom_map[begin_idx], atom_map[end_idx]
                    )
                    if sub_bond:
                        highlight_bond_colors[sub_bond.GetIdx()] = color

        # Drawing options for the subgraph
        d2d.drawOptions().useBWAtomPalette()
        d2d.drawOptions().continuousHighlight = False
        d2d.drawOptions().highlightBondWidthMultiplier = 24
        d2d.drawOptions().setHighlightColour((0.9, 0.9, 0))
        d2d.drawOptions().fillHighlights = False
        d2d.drawOptions().bondLineWidth = 3
        d2d.drawOptions().scale = 1
        d2d.drawOptions().atomLabelFontSize = 24
        d2d.drawOptions().padding = 0.1

        d2d.DrawMolecule(
            submol,
            highlightAtoms=highlight_atom_ids,
            highlightBonds=list(highlight_bond_colors.keys()),
            highlightBondColors=highlight_bond_colors,
        )
        d2d.FinishDrawing()
        return d2d.GetDrawingText()
