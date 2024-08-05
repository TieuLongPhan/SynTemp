import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from typing import List, Optional


def plot_reaction_times(df: pd.DataFrame, ax: Axes) -> None:
    """
    Plot average reaction times for different datasets with error bars.

    This function visualizes the average reaction times for 'All', 'Biochemical',
    and 'Chemical' datasets, using error bars to indicate variability. It expects
    a DataFrame with reaction times across these categories and an Axes object
    to draw the plot.

    Parameters:
    - df (pd.DataFrame): DataFrame containing the reaction time data.
      It should have columns for each mapper (e.g., 'Mapper1', 'Mapper2', etc.)
      and a 'Number of Reactions' column for normalization.
    - ax (matplotlib.axes.Axes): Matplotlib Axes object where the plot will be drawn.

    Returns:
    - None: The function directly modifies the provided Axes object with the plot.

    Notes:
    - The function uses LaTeX rendering for text, which requires a LaTeX system installed.
    - It sets the y-axis limit based on assumed data range, which may need adjustment.
    """
    # Enable LaTeX rendering in matplotlib
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
    sns.set(style="darkgrid")  # Correct style for background grid

    # Define colors for better differentiation in the plot
    colors = ["#4c92c3", "#ffdd57", "#ff6f61"]

    # Define bar width and positions
    bar_width = 0.25
    x = np.arange(len(df.columns[:-1]))  # Position indexes for mappers

    # Plot each category for mappers with error bars
    for i, mapper in enumerate(df.columns[:-1]):
        all_avg, all_std = _calculate_stats(df, mapper)
        bio_avg, bio_std = _calculate_stats(df.iloc[:2], mapper)
        chem_avg, chem_std = _calculate_stats(df.iloc[2:], mapper)

        # Plotting bars for each category on the provided ax object
        ax.bar(
            x[i] - bar_width,
            all_avg,
            width=bar_width,
            color=colors[0],
            label=r"$\textit{All dataset}$" if i == 0 else "",
            yerr=all_std,
            capsize=5,
        )
        ax.bar(
            x[i],
            bio_avg,
            width=bar_width,
            color=colors[1],
            label=r"$\textit{Biochemical dataset}$" if i == 0 else "",
            yerr=bio_std,
            capsize=5,
        )
        ax.bar(
            x[i] + bar_width,
            chem_avg,
            width=bar_width,
            color=colors[2],
            label=r"$\textit{Chemical dataset}$" if i == 0 else "",
            yerr=chem_std,
            capsize=5,
        )

        # Add text labels for average values above each bar
        ax.text(
            x[i] - bar_width,
            all_avg + 0.05,
            f"{all_avg:.2f}",
            ha="center",
            va="bottom",
            color="black",
            fontsize=18,
        )
        ax.text(
            x[i],
            bio_avg + 0.05,
            f"{bio_avg:.2f}",
            ha="center",
            va="bottom",
            color="black",
            fontsize=18,
        )
        ax.text(
            x[i] + bar_width,
            chem_avg + 0.05,
            f"{chem_avg:.2f}",
            ha="center",
            va="bottom",
            color="black",
            fontsize=18,
        )

    # Set labels, title, and ticks for the plot
    ax.set_ylabel(
        "Average Time per Reaction (seconds)", fontsize=24, weight="bold", color="black"
    )
    ax.set_xticks(x)
    ax.set_xticklabels(
        df.columns[:-1], rotation=45, fontsize=18, weight="bold", color="black"
    )
    ax.set_title(r"A. Processing time benchmarking", fontsize=28, weight="bold")
    ax.tick_params(axis="y", labelsize=18, labelcolor="black")
    ax.legend(fontsize=20, loc="upper left", frameon=True, edgecolor="black")
    ax.set_ylim(0, 7)  # Adjust this limit based on your data range


def _calculate_stats(sub_df: pd.DataFrame, column: str) -> tuple:
    """
    Helper function to calculate the average and standard deviation of reaction times.

    Parameters:
    - sub_df (pd.DataFrame): Subset of the main DataFrame.
    - column (str): Name of the column for which to calculate statistics.

    Returns:
    - tuple: A tuple containing the average and standard deviation.
    """
    avg = sub_df[column].sum() / sub_df["Number of Reactions"].sum()
    std = np.sqrt(
        np.sum((sub_df[column] / sub_df["Number of Reactions"] - avg) ** 2)
        / sub_df["Number of Reactions"].sum()
    )
    return avg, std


def plot_heatmap(data: pd.DataFrame, save_path: Optional[str] = None) -> Figure:
    """
    Plot a heatmap to visualize accuracies across datasets and mappers with optional saving.

    This function generates a heatmap of accuracies, displaying data from different mappers
    and datasets. It assumes that the data is already prepared with 'mapper' as index.
    The plot can optionally be saved to a specified path.

    Parameters:
    - data (pd.DataFrame): A DataFrame containing the heatmap data with mappers as index
      and datasets as columns.
    - save_path (str, optional): Path to save the generated heatmap as a PDF file.
      If None, the plot is not saved.

    Returns:
    - fig (matplotlib.figure.Figure): A Matplotlib Figure object containing the heatmap.

    Raises:
    - ValueError: If the data does not contain the expected index or columns.
    """
    # Set up LaTeX rendering and seaborn theme
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath}")
    sns.set_theme(style="darkgrid")

    # Validate the DataFrame structure
    if "mapper" not in data.index.names:
        data = data.set_index("mapper")

    # Creating a figure object for more control
    fig, ax = plt.subplots(figsize=(16, 8))

    # Create a heatmap
    sns.heatmap(
        data,
        ax=ax,
        annot=True,
        cmap="coolwarm",
        fmt=".1f",
        linewidths=0.3,
        linecolor="white",
        cbar=True,
        cbar_kws={"label": r"Accuracy (\%)", "orientation": "vertical"},
        annot_kws={"size": 18},
    )

    # Customize the color bar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label(r"Accuracy (\%)", size=18)

    # Customizing tick marks
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=18,
        fontweight="bold",
        color="darkred",
    )
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=18,
        fontweight="bold",
        color="darkgreen",
    )

    plt.tight_layout()  # Adjust layout to make sure everything fits without overlap

    # Check if save_path is provided
    if save_path:
        fig.savefig(save_path, dpi=600, bbox_inches="tight", pad_inches=0)

    return fig


def plot_accuracy_success_rate_subplot(
    df: pd.DataFrame,
    accuracy_cols: List[str],
    success_cols: List[str],
    titles: List[str],
    axes: np.ndarray,
    save_path: Optional[str] = None,
) -> None:
    """
    Plots 1x2 subplots for accuracy and success rates using Seaborn from specified columns
    in a DataFrame.

    Parameters:
    ----------
    df : pd.DataFrame
        The DataFrame containing the data. Must include 'mapper' and specified columns for
        accuracy and success rates.
    accuracy_cols : List[str]
        A list of column names in the DataFrame that contain accuracy data.
    success_cols : List[str]
        A list of column names in the DataFrame that contain success rate data.
    titles : List[str]
        Titles for each subplot.
    axes : np.ndarray
        An array of Matplotlib axes objects where the plots will be drawn.
    save_path : str, optional
        File path to save the plot as a PDF. If None, the plot is not saved.
        (default is None)

    Raises:
    -------
    ValueError
        If any specified columns are missing from the DataFrame.
    """

    # Check if all specified columns exist in the DataFrame
    for col in accuracy_cols + success_cols:
        if col not in df.columns:
            raise ValueError(f"Column {col} does not exist in the DataFrame.")

    for idx, ax in enumerate(axes):
        # Ensure 'mapper' is treated as categorical
        df["mapper"] = df["mapper"].astype(str)

        # Create a temporary DataFrame to facilitate plotting with Seaborn
        temp_df = df[["mapper", accuracy_cols[idx], success_cols[idx]]].melt(
            id_vars=["mapper"], var_name="Metric", value_name="Percentage"
        )

        # Mapping column names to more user-friendly names for the legend
        temp_df["Metric"] = temp_df["Metric"].map(
            {accuracy_cols[idx]: "Accuracy", success_cols[idx]: "Success Rate"}
        )

        # Plotting using Seaborn on the specified axes
        sns.barplot(
            x="mapper",
            y="Percentage",
            hue="Metric",
            data=temp_df,
            palette="coolwarm",
            ax=ax,
        )

        # Enhance annotations above bars for clarity
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.annotate(
                    format(height, ".1f"),
                    (p.get_x() + p.get_width() / 2.0, height),
                    ha="center",
                    va="center",
                    xytext=(0, 9),
                    textcoords="offset points",
                    fontsize=18,
                )

        ax.set_title(titles[idx], fontsize=28, fontweight="bold")
        ax.set_xlabel(None)
        ax.set_ylabel(r"Percentage (\%)", fontsize=24, fontweight="semibold")

        ax.set_xticklabels(
            ax.get_xticklabels(), rotation=45, fontsize=18, fontweight="bold"
        )
        ax.set_yticklabels(
            [f"{int(x)}%" for x in ax.get_yticks()], fontsize=18, fontweight="bold"
        )
        ax.grid(
            True, which="both", linestyle="--", linewidth=0.5, color="gray", alpha=0.5
        )
        ax.legend([], [], frameon=False)

    # Handle legend and save plot if path provided
    handles, labels = axes[-1].get_legend_handles_labels()
    fig = plt.gcf()
    fig.legend(
        handles, labels, loc="center left", bbox_to_anchor=(0.45, 0.02), fontsize=20
    )
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches="tight")

    plt.show()
