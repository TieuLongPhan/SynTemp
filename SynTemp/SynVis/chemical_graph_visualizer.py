import random
import networkx as nx
import matplotlib.pyplot as plt
from typing import Optional, Dict, Any

color_scheme = {
    "H": "#FFFFFF",  # White
    "C": "#909090",  # Gray
    "N": "#3050F8",  # Blue
    "O": "#FF0D0D",  # Red
    "F": "#90E050",  # Green
    "Cl": "#1FF01F",  # Green
    "Br": "#A62929",  # Dark Red/Brown
    "I": "#940094",  # Purple
    "P": "#FF8000",  # Orange
    "S": "#FFFF30",  # Yellow
    "Na": "#E0E0E0",  # Light Gray
    "K": "#8F40D4",  # Light Purple
    "Ca": "#3DFF00",  # Light Green
    "Mg": "#8AFF00",  # Light Green
    "Fe": "#B7410E",  # Rust Red
    "Zn": "#7D80B0",  # Light Blue
    "Cu": "#C88033",  # Copper Red/Orange
    "Ag": "#C0C0C0",  # Light Gray
    "Au": "#FFD123",  # Gold Yellow
    "Hg": "#B8B8D0",  # Silver Gray
    "Pb": "#575961",  # Dark Gray
    "Al": "#BFA6A6",  # Light Gray
    "Si": "#F0C8A0",  # Light Brown
    "B": "#FFA1A1",  # Pink
    "As": "#BD80E3",  # Light Gray
    "Sb": "#9E63B5",  # Dark Gray
    "Se": "#FFA100",  # Light Pink
    "Te": "#D47A00",  # Gray
    "Cd": "#FFD98F",  # Light Blue/Gray
    "Ti": "#BFC2C7",  # Light Gray
    "V": "#A6A6AB",  # Light Gray/Blue
    "Cr": "#8A99C7",  # Steel Gray
    "Mn": "#9C7AC7",  # Gray
    "Co": "#FF7A00",  # Light Pink
    "Ni": "#4DFF4D",  # Light Green
}


class ChemicalGraphVisualizer:
    def __init__(
        self,
        seed: Optional[int] = None,
        element_colors: Optional[Dict[str, str]] = color_scheme,
    ):
        """
        Initialize the visualizer with optional seed and color scheme.

        Parameters:
        seed (int, optional): Seed for random number generator for reproducibility.
        element_colors (dict, optional): Dictionary mapping elements to their color codes.
        """
        # Define a popular color scheme in chemistry if not provided
        if element_colors is None:
            self.element_colors = {
                "H": "#FFFFFF",  # White
                "C": "#909090",  # Gray
                "N": "#3050F8",  # Blue
                "O": "#FF0D0D",  # Red
                "F": "#90E050",  # Green
                "Cl": "#1FF01F",  # Green
                # Additional elements can be added here
            }
        else:
            self.element_colors = element_colors
        self.seed = seed

    def graph_vis(
        self,
        G: nx.Graph,
        node_size: int = 100,
        visualize_edge_weight: bool = False,
        edge_font_size: int = 10,
        show_node_labels: bool = False,
        node_label_font_size: int = 12,
        ax: Optional[plt.Axes] = None,
    ) -> None:
        """
        Visualize a NetworkX graph with standard representation.

        Parameters:
        G (nx.Graph): The graph to visualize.
        node_size (int): The size of the nodes.
        visualize_edge_weight (bool): Whether to display edge weights.
        edge_font_size (int): Font size for edge labels.
        show_node_labels (bool): Whether to show labels on the nodes.
        node_label_font_size (int): Font size for node labels.
        """
        # Set random seed for reproducibility
        if self.seed is not None:
            random.seed(self.seed)

        # Get colors for each node
        node_colors = [
            self.element_colors.get(G.nodes[node]["element"], "#000000")
            for node in G.nodes()
        ]

        # Draw the graph
        pos = nx.spring_layout(G, seed=self.seed)  # Use spring layout

        if ax is None:
            ax = plt.gca()  # Get current axes if not provided

        if show_node_labels:
            node_labels = {node: G.nodes[node]["element"] for node in G.nodes()}
            nx.draw(
                G,
                pos,
                ax=ax,
                with_labels=True,
                labels=node_labels,
                node_color=node_colors,
                node_size=node_size,
                font_size=node_label_font_size,
                font_weight="bold",
            )
        else:
            nx.draw(
                G,
                pos,
                ax=ax,
                with_labels=False,
                node_color=node_colors,
                node_size=node_size,
                font_weight="bold",
            )

        # Get edge labels if needed
        if visualize_edge_weight:
            edge_labels = {(u, v): G.edges[u, v]["order"] for u, v in G.edges()}
            nx.draw_networkx_edge_labels(
                G, pos, ax=ax, edge_labels=edge_labels, font_size=edge_font_size
            )

    def its_vis(
        self,
        G: nx.Graph,
        node_size: int = 100,
        show_node_labels: bool = False,
        node_label_font_size: int = 12,
        ax: Optional[plt.Axes] = None,
    ) -> None:
        """
        Visualize a NetworkX graph with edge colors indicating bond changes.

        Parameters:
        G (nx.Graph): The graph to visualize.
        node_size (int): The size of the nodes.
        show_node_labels (bool): Whether to show labels on the nodes.
        node_label_font_size (int): Font size for node labels.
        """
        # Set random seed for reproducibility
        if self.seed is not None:
            random.seed(self.seed)

        # Draw the graph
        pos = nx.spring_layout(G, seed=self.seed)  # Use spring layout
        # Get colors for each node
        node_colors = [
            self.element_colors.get(G.nodes[node]["element"], "#000000")
            for node in G.nodes()
        ]

        # Determine edge colors based on 'order'
        edge_colors = []
        for u, v, data in G.edges(data=True):
            order = data.get("standard_order", 0)
            # order = tuple(0 if isinstance(x, str) else float(x) for x in order)
            if order == 0:
                edge_colors.append("black")  # Normal bond
            elif order < 0:
                edge_colors.append("blue")  # Increasing bond
            else:
                edge_colors.append("red")  # Breaking bond

        if ax is None:
            ax = plt.gca()  # Get current axes if not provided

        if show_node_labels:
            node_labels = {node: G.nodes[node]["element"] for node in G.nodes()}
            nx.draw(
                G,
                pos,
                ax=ax,
                with_labels=True,
                labels=node_labels,
                node_color=node_colors,
                node_size=node_size,
                font_size=node_label_font_size,
                font_weight="bold",
                edge_color=edge_colors,
            )
        else:
            nx.draw(
                G,
                pos,
                ax=ax,
                with_labels=False,
                node_color=node_colors,
                node_size=node_size,
                font_weight="bold",
                edge_color=edge_colors,
            )

    def visualize_all(
        self,
        graph_tuple,
        figsize=(15, 5),
        left_graph_title="Reactants",
        k_graph_title="ITS Graph",
        right_graph_title="Products",
        save_path=None,
    ):
        """
        Visualize reactants, ITS graph, and products in one figure.

        Parameters:
        graph_tuple (tuple): Tuple of NetworkX graphs (reactants, products, ITS graph).
        """
        # Unpack the tuple
        reactants_graph, products_graph, its_graph = graph_tuple

        # Create a figure with subplots
        fig, axs = plt.subplots(1, 3, figsize=figsize)  # Adjust figsize as needed

        # Initialize visualizer

        # Visualize each graph on its respective subplot
        self.graph_vis(reactants_graph, ax=axs[0], show_node_labels=True)
        self.its_vis(its_graph, ax=axs[1], show_node_labels=True)
        self.graph_vis(products_graph, ax=axs[2], show_node_labels=True)

        # Set titles for subplots
        axs[0].set_title(left_graph_title, fontsize=24, weight="bold")
        axs[1].set_title(k_graph_title, fontsize=24, weight="bold")
        axs[2].set_title(right_graph_title, fontsize=24, weight="bold")

        plt.tight_layout()
        if save_path is not None:
            plt.savefig(save_path, dpi=600)
        plt.show()

    def multigraph_vis(
        self,
        G: nx.MultiGraph,
        node_size: int = 100,
        show_node_labels: bool = False,
        node_label_font_size: int = 12,
        ax: Optional[Any] = None,
    ) -> None:
        if ax is None:
            _, ax = plt.subplots()

        pos = nx.spring_layout(G, seed=self.seed)
        # Get colors for each node
        node_colors = [
            self.element_colors.get(G.nodes[node]["element"], "#000000")
            for node in G.nodes()
        ]

        # Draw nodes
        nx.draw_networkx_nodes(
            G,
            pos,
            ax=ax,
            node_size=node_size,
            node_color=[
                self.element_colors.get(G.nodes[node]["element"], "#000000")
                for node in G.nodes()
            ],
        )

        # Draw edges based on bond order
        for u, v, data in G.edges(data=True):
            order = data.get("order", 1)
            if order == 1:
                style = "solid"
            elif order == 1.5:
                style = "dotted"
            elif order == 2:
                style = "double"
            else:  # Assuming order 3 for triple bond
                style = "triple"

            if style == "double":
                nx.draw_networkx_edges(
                    G, pos, ax=ax, edgelist=[(u, v)], style=":", width=2
                )
                nx.draw_networkx_edges(
                    G,
                    pos,
                    ax=ax,
                    edgelist=[(u, v)],
                    style="solid",
                    width=2,
                    edge_color="white",
                    alpha=1,
                )
            elif style == "triple":
                nx.draw_networkx_edges(
                    G, pos, ax=ax, edgelist=[(u, v)], style="solid", width=4
                )
            elif style == "dotted":
                nx.draw_networkx_edges(
                    G, pos, ax=ax, edgelist=[(u, v)], style="dotted", width=2
                )
            else:  # single or dotted
                nx.draw_networkx_edges(G, pos, ax=ax, edgelist=[(u, v)], style=style)

        # Draw node labels if needed
        if show_node_labels:
            node_labels = {node: G.nodes[node]["element"] for node in G.nodes()}
            nx.draw(
                G,
                pos,
                ax=ax,
                with_labels=True,
                labels=node_labels,
                node_color=node_colors,
                node_size=node_size,
                font_size=node_label_font_size,
                font_weight="bold",
            )
        else:
            nx.draw(
                G,
                pos,
                ax=ax,
                with_labels=False,
                node_color=node_colors,
                node_size=node_size,
                font_weight="bold",
            )

        ax.axis("off")
