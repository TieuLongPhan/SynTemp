import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Optional
from matplotlib.axes import Axes

def create_pie_chart(data, column, ax=None, title=None, color_pallet="pastel"):
    """
    Generates a pie chart for the specified column from a list of dictionaries.
    Displays percentage labels inside the slices only and category names in an external
    legend without percentages. Allows customization of the plot title, supporting LaTeX
    formatted strings.

    Parameters:
        data (list of dict): Data to plot.
        column (str): Column name to plot percentages for.
        ax (matplotlib.axes.Axes, optional): Matplotlib axis object to plot on.
        title (str, optional): Title for the pie chart, supports LaTeX formatted strings.

    Returns:
        matplotlib.axes.Axes: The axis with the pie chart.
    """
    # Enable LaTeX formatting for better quality text rendering
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")

    # Convert list of dictionaries to DataFrame
    df = pd.DataFrame(data)

    # Calculate percentage
    percentage = df[column].value_counts(normalize=True) * 100

    # Define a color palette using Seaborn
    colors = sns.color_palette(color_pallet, len(percentage))

    # Create pie plot
    if ax is None:
        fig, ax = plt.subplots()

    wedges, texts, autotexts = ax.pie(
        percentage,
        startangle=90,
        colors=colors,
        autopct="%1.1f%%",
        pctdistance=0.85,
        explode=[0.05] * len(percentage),
    )

    # Draw a circle at the center of pie to make it look like a donut
    centre_circle = plt.Circle((0, 0), 0.70, fc="white")
    ax.add_artist(centre_circle)

    # Equal aspect ratio ensures that pie is drawn as a circle.
    ax.axis("equal")

    # Add legend with category names only
    ax.legend(
        wedges,
        [f"{label}" for label in percentage.index],
        title=column,
        loc="upper right",
        bbox_to_anchor=(0.6, 0.1, 0.6, 1),
        prop={"size": 16},
        title_fontsize=16,
    )  # Set label font size

    # Set title using LaTeX if provided, else default to a generic title
    if title:
        ax.set_title(title, fontsize=24)
    else:
        ax.set_title(f"Pie Chart of {column}", fontsize=32)

    # Enhance the font size and color of the autotexts
    for autotext in autotexts:
        autotext.set_color("black")
        autotext.set_fontsize(18)

    return ax

def count_column_values(data, column):
    """
    Count the occurrences of each unique value in the specified column from a list of 
    dictionaries. Treats all data types, including lists, as single entities by converting 
    lists to strings.

    Parameters:
        data (list of dict): The data to process.
        column (str): The column to count values from.

    Returns:
        dict: A dictionary with keys as unique values (strings if lists) and values as the 
        count of occurrences.
    """
    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data)

    # Handle if the column contains lists
    if df[column].dtype == object and df[column].apply(
        lambda x: isinstance(x, list)).all():
        
        df[column] = df[column].apply(lambda x: str(x))
    
    # Count occurrences of each unique value
    return df[column].value_counts().to_dict()



def plot_rules_distribution(
    rules: Dict[str, int], 
    rule_type: str = 'single', 
    ax: Optional[Axes] = None, 
    title: Optional[str] = None, 
    refinement: bool = False, 
    threshold: float = 1, 
    remove: bool = True, 
    color_pallet: str = 'pastel'
) -> None:
    """
    Plots the distribution of rules in a bar chart, optionally combining all entries under the threshold into a 
    single category 'Under 1%' if `refinement` is True.

    Parameters:
        rules (Dict[str, int]): Dictionary with rule counts keyed by rule name, where the values are counts.
        rule_type (str, optional): Specifies the type of rules to plot ('single' or 'complex'). Default is 'single'.
        ax (matplotlib.axes.Axes, optional): Matplotlib axis object to plot on. If None, a new figure is created.
        title (str, optional): Optional title for the chart. If None, a default title based on `rule_type` is used.
        refinement (bool, optional): If True, combines all percentages under the threshold into one category 'Under 1%'.
                                     Default is False.
        threshold (float, optional): The percentage threshold under which all categories are combined into 'Under 1%' 
                                     if `refinement` is True. Default is 1.
        remove (bool, optional): If True, removes the last category from the plot. Default is True.
        color_pallet (str, optional): Color palette to use for the plot. Default is 'pastel'.

    Returns:
        None: The function directly modifies the `ax` object or creates a new plot.
    """
    # Calculate total counts for the rules
    total_rules = sum(rules.values())

    # Convert counts to percentages and optionally combine small values
    if refinement:
        refined_rules = {}
        small_value_aggregate = 0
        for key, value in rules.items():
            percentage = value / total_rules * 100
            if percentage < threshold:
                small_value_aggregate += percentage
            else:
                refined_rules[key] = percentage
        if small_value_aggregate > 0:
            refined_rules['Under 1%'] = small_value_aggregate
        percentages = list(refined_rules.values())
        types_of_rules = list(refined_rules.keys())
        if remove:
            percentages = percentages[:-1]
            types_of_rules = types_of_rules[:-1]
    else:
        percentages = [value / total_rules * 100 for value in rules.values()]
        types_of_rules = list(rules.keys())

    # Set style
    sns.set(style="whitegrid")

    # Enable LaTeX rendering in matplotlib
    plt.rc('text', usetex=True)
    plt.rc('text.latex', preamble=r'\usepackage{amsmath}')  # Ensure amsmath is loaded

    # Create figure and axis if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

    # Plot the data
    sns.barplot(ax=ax, x=types_of_rules, y=percentages, palette=color_pallet)
    if title:
        ax.set_title(title, fontsize=24)
    else:
        ax.set_title(f'Distribution of {rule_type.capitalize()} Rules', fontsize=16)
    ax.set_xlabel('Type of Rings', fontsize=18)
    ax.set_ylabel(r'Percentage (\%)', fontsize=18)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")


    # Set font size for x-tick and y-tick labels
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    # Add text labels above the bars
    for index, value in enumerate(percentages):
        ax.text(index, value + 0.5, f'{value:.1f}%', ha='center', va='bottom', fontsize=18)

    # Only show plot if ax is not provided (i.e., we created the figure here)
    if ax is None:
        plt.show()



def plot_heatmap(
    data,
    title="Heatmap of Test Counts by Topo Type and Reaction Step",
    color_palette="coolwarm",
    title_fontsize=24,
    label_fontsize=20,
    annot_fontsize=18,
    cbar_label_fontsize=18,
    legend_fontsize=24,
    xtick_fontsize=18,
    ytick_fontsize=18,
    ax=None,
):
    """
    Plots a heatmap based on the provided dataset with options for customization, specific 
    aggregation, and an enhanced legend.

    Parameters:
        data (list of dict): Data to be visualized.
        title (str, optional): Title for the heatmap. Defaults to a generic title if none 
        provided.
        color_palette (str, optional): Color palette for the heatmap. 
        Defaults to 'coolwarm'.
        title_fontsize (int, optional): Font size for the title. Defaults to 16.
        label_fontsize (int, optional): Font size for the axis labels. Defaults to 14.
        annot_fontsize (int, optional): Font size for the annotations. Defaults to 12.
        cbar_label_fontsize (int, optional): Font size for the color bar label. 
        Defaults to 12.
        legend_fontsize (int, optional): Font size for the legend. Defaults to 10.
        ax (matplotlib.axes.Axes, optional): Matplotlib axis object to plot on. 
        If none, a new figure is created.
    """
    # Convert input data to DataFrame
    df = pd.DataFrame(data)
    df["Test"] = 1

    # Define a custom aggregation function to calculate percentages
    def custom_agg(series):
        total = series.sum()
        return (
            total / len(data)
        ) * 100  # Dividing by the total number of data points and multiplying by 100

    # Create pivot table for heatmap using the custom aggregation function
    pivot_table = df.pivot_table(
        index="Topo Type", columns="Reaction Step", values="Test", aggfunc=custom_agg
    )

    # Check if an axis is provided; if not, create a new figure and axis
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))

    # Plot heatmap on the provided or created axis
    heatmap = sns.heatmap(
        pivot_table,
        annot=True,
        cmap=color_palette,
        fmt=".1f",
        ax=ax,
        cbar_kws={"label": r"Percentage (\%)"},
    )

    # Customize the title and axis labels font size
    ax.set_title(title, fontsize=title_fontsize)
    ax.set_ylabel("Topo Type", fontsize=label_fontsize)
    ax.set_xlabel("Reaction Step", fontsize=label_fontsize)
    # ax.set_xticks

    # Customize the font size of the annotations
    for text in heatmap.texts:
        text.set_fontsize(annot_fontsize)

    # Set font size for x-tick and y-tick labels
    ax.tick_params(axis="x", labelsize=xtick_fontsize)
    ax.tick_params(axis="y", labelsize=ytick_fontsize)
    # Customize the font size of the color bar label
    heatmap.figure.axes[-1].yaxis.label.set_size(cbar_label_fontsize)

    # Create a legend with specified font size
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            title="Legend",
            loc="upper right",
            bbox_to_anchor=(1.05, 1),
            fontsize=legend_fontsize,
        )

    if not ax:
        plt.show()
