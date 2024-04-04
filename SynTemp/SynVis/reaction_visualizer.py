import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image
from typing import List, Dict, Tuple, Optional


class ReactionVisualizer:
    """
    Plot and optionally compare two chemical reactions for visualization using data from a given source.

    This method allows the visualization of chemical reactions represented as strings. It can plot a single reaction or compare two reactions side-by-side, depending on the 'compare' parameter. The reactions are visualized using RDKit and plotted using Matplotlib.

    Parameters
    ----------
    data : dict
        A dictionary containing the reaction data.
    old_reaction_col : str
        The key in 'data' for the old reaction string.
    new_reaction_col : str
        The key in 'data' for the new reaction string.
    compare : bool, optional
        If True, both the old and new reactions are plotted side by side for comparison. Default is False.
    orientation : str, optional
        The layout orientation of the plots, either 'vertical' or 'horizontal'. Default is 'vertical'.
    savefig : bool, optional
        If True, the figure is saved to the specified pathname. Default is False.
    pathname : str, optional
        The pathname where the figure will be saved, if 'savefig' is True.
    dpi : int, optional
        The resolution of the figure in dots per inch. Default is 300.

    Examples
    --------
    # Example usage of the plot_reactions method
    visualizer = ReactionVisualizer()
    reaction_data = {
        'old_reaction': 'C1=CC=CC=C1.CCO>>C1=CC=CC=C1OCCO',
        'new_reaction': 'C1=CC=CC=C1>>CCO'
    }
    visualizer.plot_reactions(reaction_data, 'old_reaction', 'new_reaction', compare=True, orientation='horizontal')
    
    # To save the plot as an image
    visualizer.plot_reactions(reaction_data, 'old_reaction', 'new_reaction', compare=True, savefig=True, pathname='reaction_comparison.png')
    """

    def __init__(
        self,
        compare: bool = True,
        orientation: str = "vertical",
        figsize: Tuple[int, int] = (10, 5),
        label_position: str = "below",
        dpi: int = 300,
        bond_line_width: float = 6,
        atom_label_font_size: int = 50,
        padding: float = 0.001,
    ) -> None:
        """
        Initialize the ReactionVisualizer.

        Parameters:
        - compare (bool): Whether to compare old and new reactions side by side. Default is True.
        - orientation (str): Orientation of comparison ('vertical' or 'horizontal'). Default is 'vertical'.
        - figsize (Tuple[int, int]): Figure size (width, height) in inches. Default is (10, 5).
        - label_position (str): Position of labels ('above' or 'below') the reaction images. Default is 'below'.
        - dpi (int): Dots per inch for image resolution. Default is 300.
        - bond_line_width (float): Width of bond lines in the reaction image. Default is 6.
        - atom_label_font_size (int): Font size for atom labels in the reaction image. Default is 50.
        - padding (float): Padding around the drawing in the reaction image. Default is 0.001.
        """
        self.compare = compare
        self.orientation = orientation
        self.figsize = figsize
        self.label_position = label_position
        self.dpi = dpi
        self.bond_line_width = bond_line_width
        self.atom_label_font_size = atom_label_font_size
        self.padding = padding

    
    @staticmethod
    def draw_molecule_with_atom_numbers(mol):
        """
        Draw a molecule with atom numbers annotated.
        
        Parameters
        ----------
        mol : RDKit Molecule object
            The molecule to be drawn with atom numbers.
        
        Returns
        -------
        mol : RDKit Molecule object
            The molecule with atom numbers.
        """
        mol_with_atom_numbers = Chem.Mol(mol)
        for atom in mol_with_atom_numbers.GetAtoms():
            atom.SetProp('atomLabel', str(atom.GetIdx()))
        return mol_with_atom_numbers


    def visualize_reaction(self, reaction_str: str, show_atom_numbers: bool = False) -> Image:
        """
        Visualize a single chemical reaction and return the image.

        Parameters
        ----------
        reaction_str : str
            A string representation of the reaction (e.g., 'C1=CC=CC=C1>>CCO').

        Returns
        -------
        PIL.Image.Image
            An image of the chemical reaction.
        """
        # Parse reactants and products from the reaction string
        reactants_str, products_str = reaction_str.split('>>')
        reactants = [Chem.MolFromSmiles(smiles) for smiles in reactants_str.split('.')]
        products = [Chem.MolFromSmiles(smiles) for smiles in products_str.split('.')]

        if show_atom_numbers:
            reactants = [self.draw_molecule_with_atom_numbers(mol) for mol in reactants]
            products = [self.draw_molecule_with_atom_numbers(mol) for mol in products]
            
        
        rxn = AllChem.ChemicalReaction()
        for reactant in reactants:
            rxn.AddReactantTemplate(reactant)
        for product in products:
            rxn.AddProductTemplate(product)

        # Set up RDKit drawer with customizable parameters
        drawer = rdMolDraw2D.MolDraw2DCairo(2000, 600)  # Adjust canvas size as needed
        opts = drawer.drawOptions()
        opts.bondLineWidth = self.bond_line_width  # Increase bond line width
        opts.atomLabelFontSize = self.atom_label_font_size  # Increase font size for atom labels
        opts.padding = self.padding  # Adjust padding around the drawing

        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        return Image.open(io.BytesIO(drawer.GetDrawingText()))


    def plot_reactions(
        self, 
        data: Dict[str, str], 
        old_reaction_col: str, 
        new_reaction_col: str, 
        compare: bool = False, 
        orientation: str = 'vertical', 
        savefig: bool = False, 
        pathname: Optional[str] = None, 
        dpi: int = 300,
        show_atom_numbers: bool = False
        ) -> None:
        """
        Plot one or two chemical reactions for visualization.

        Parameters
        ----------
        data : Dict[str, str]
            A dictionary containing the data with keys as column names and values as reaction strings.
        old_reaction_col : str
            The column name containing the string representation of the old reaction.
        new_reaction_col : str
            The column name containing the string representation of the new reaction.
        compare : bool, optional
            If True, both the old and new reactions are plotted. Default is False.
        orientation : str, optional
            The layout orientation of the plots ('vertical' or 'horizontal'). Default is 'vertical'.
        savefig : bool, optional
            If True, saves the figure to the specified pathname. Default is False.
        pathname : str, optional
            Pathname to save the figure, if savefig is True.
        dpi : int, optional
            Resolution of the figure in dots per inch. Default is 300.
        """
        # Get the old and new reaction strings from the data
        old_reaction_str = data[old_reaction_col]
        new_reaction_str = data[new_reaction_col]
        
        # Create reaction images
        new_reaction_image = self.visualize_reaction(new_reaction_str, show_atom_numbers)
        old_reaction_image = self.visualize_reaction(old_reaction_str, show_atom_numbers) if compare else None

        # Set up plot layout
        if compare:
            nrows, ncols = (2, 1) if orientation == 'vertical' else (1, 2)
            fig, axs = plt.subplots(nrows, ncols, figsize=self.figsize, dpi=dpi)
            ax_new = axs[1] if orientation == 'vertical' else axs[0]
            ax_old = axs[0] if orientation == 'vertical' else axs[1]
        else:
            fig, ax = plt.subplots(figsize=self.figsize, dpi=dpi)

        # Plotting logic
        if compare:
            ax_old.imshow(old_reaction_image)
            ax_old.axis('off')
            ax_new.imshow(new_reaction_image)
            ax_new.axis('off')
        else:
            ax.imshow(new_reaction_image)
            ax.axis('off')

        # Setting titles
        label_y_position = -0.1 if self.label_position == 'below' else 1.1
        if compare:
            ax_old.set_title('Old Reaction', position=(0.5, label_y_position), weight='bold')
            ax_new.set_title('New Reaction', position=(0.5, label_y_position), weight='bold')
        else:
            ax.set_title('New Reaction', position=(0.5, label_y_position), weight='bold')

        # Saving the figure
        if savefig and pathname:
            fig.savefig(pathname, dpi=dpi)

        #plt.tight_layout()
        return fig