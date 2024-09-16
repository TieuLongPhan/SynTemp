import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Patch
from matplotlib.axes import Axes
import seaborn as sns
from collections import Counter
from typing import List, Dict, Any, Optional, Union, Sequence


def plot_top_rules_with_seaborn(
    data: List[Dict[str, float]], top_n: int = 10, ax: Optional[Axes] = None
) -> Optional[Axes]:
    """
    Generates a bar plot for the top N rules based on their 'Percentage'
    from a list of dictionaries using Seaborn.

    Parameters:
    - data (List[Dict[str, float]]): A list of dictionaries containing rule data,
    where each dictionary has 'Cluster_id' as a str and 'Percentage' as a float.
    - top_n (int, optional): Number of top entries to plot. Defaults to 10.
    - ax (matplotlib.axes.Axes, optional): Matplotlib axis object on which to plot.
    If None, a new axis will be created. Defaults to None.

    Returns:
    - Optional[matplotlib.axes.Axes]: The axis with the bar plot or None if
    an error occurs.

    Raises:
    - ValueError: If `data` is empty, `top_n` is non-positive, or the required columns
    are missing in the data.
    """
    if not data or top_n <= 0:
        raise ValueError("Invalid data or number of rules to plot.")

    df = pd.DataFrame(data)

    required_columns = {"Cluster_id", "Percentage"}
    if not required_columns.issubset(df.columns):
        raise ValueError("Data must include 'Cluster_id' and 'Percentage' columns.")

    df_sorted = df.sort_values(by="Percentage", ascending=False).head(top_n)

    sns.set_theme(style="whitegrid")

    # Prepare the figure and axes for plotting if ax is not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 6))

    # Use LaTeX for text rendering in plots
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath}")

    # Create the bar plot using seaborn on the specified axes 'ax'
    barplot = sns.barplot(
        x="Cluster_id",
        y="Percentage",
        data=df_sorted,
        palette="coolwarm",
        ax=ax,
        order=df_sorted["Cluster_id"].tolist(),
    )

    # Add labels on top of each bar
    for p in barplot.patches:
        ax.annotate(
            rf"{p.get_height():.1f}%",
            (p.get_x() + p.get_width() / 2.0, p.get_height()),
            ha="center",
            va="center",
            xytext=(0, 9),
            textcoords="offset points",
            fontsize=20,
        )

    # Setting plot labels and titles
    ax.set_xlabel(r"Rule ID", fontsize=24)
    ax.set_ylabel(r"Percentage (\%)", fontsize=24)
    ax.set_title(rf"Top {top_n} Popular Rules", fontsize=32, weight="medium")
    ax.tick_params(axis="both", labelsize=24)
    plt.xticks(rotation=45)  # Rotate x-axis labels for better readability

    # Show the plot if ax was not provided
    if ax is None:
        plt.show()

    return ax


def load_and_title_png(
    file_path: str,
    title: str,
    ax: Optional[Axes] = None,
    save_path: Optional[str] = None,
) -> Axes:
    """
    Load a PNG file, add a title, and optionally save the image to a specified path.

    Parameters:
    - file_path (str): Path to the PNG file.
    - title (str): Title to add to the image.
    - ax (Optional[matplotlib.axes.Axes]): Axis object on which the image is displayed.
    If None, a new figure and axis will be created.
    - save_path (Optional[str]): If provided, the path where the titled image
    will be saved.

    Returns:
    - matplotlib.axes.Axes: The axis with the displayed image and title.

    Raises:
    - FileNotFoundError: If the file at `file_path` cannot be found.
    IOError: If saving the image fails.
    """
    # Load the image
    try:
        img = mpimg.imread(file_path)
    except FileNotFoundError:
        raise FileNotFoundError("The specified file cannot be found.")

    # Use LaTeX for text rendering in plots
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath}")

    # If ax is not provided, create a new figure and axis
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(13, 6))
        created_fig = True

    # Display the image on the specified axis
    ax.imshow(img)
    ax.set_title(
        rf"{title}", fontsize=32, weight="medium"
    )  # Set the title with a specified fontsize
    ax.axis("off")  # Hide the axes

    # Save the image if a save path is provided
    if save_path:
        try:
            fig = ax.get_figure()
            fig.savefig(save_path, format="png", bbox_inches="tight", pad_inches=0.1)
        except IOError as e:
            raise IOError(f"Failed to save the image: {e}")

    # Show the plot if a new figure was created
    if created_fig:
        plt.show()

    return ax


def calculate_value_percentage(
    data: List[Dict[str, Any]], key: str
) -> Dict[Union[str, tuple], float]:
    """
    Calculates the percentage distribution of each unique value for a specified key
    in a list of dictionaries. If the value under the key is a list,
    it is treated as a unique item (tuple) and not flattened.

    Parameters:
    - data (List[Dict[str, Any]]): A list of dictionaries containing data.
    - key (str): The key in the dictionaries for which to calculate
    the percentage distribution.

    Returns:
    - Dict[Union[str, tuple], float]: A dictionary where the keys are the
    unique values or tuples (for list values) for the specified key, and
    the values are their respective percentages, rounded to two decimal places.

    Raises:
    - ValueError: If the specified key is not found in any dictionary within the list.
    """
    # Extract values associated with the key from each dictionary in the list
    try:
        values = [
            tuple(entry[key]) if isinstance(entry[key], list) else entry[key]
            for entry in data
            if key in entry
        ]
    except KeyError:
        raise ValueError(f"Key '{key}' not found in one or more dictionary entries.")

    if not values:
        raise ValueError(f"No data found for key '{key}'.")

    # Count occurrences of each value
    value_counts = Counter(values)
    total = sum(value_counts.values())

    # Calculate percentage for each value
    percentages = {k: round((v / total * 100), 2) for k, v in value_counts.items()}

    return percentages


def create_pie_chart(data, column, ax=None, title=None, color_pallet="pastel"):
    """
    Generates a pie chart for the specified column from a list of dictionaries.
    Displays percentage labels inside the slices only and category names in an external
    legend without percentages. Allows customization of the plot title, supporting LaTeX
    formatted strings.

    Parameters:
    - data (list of dict): Data to plot.
    - column (str): Column name to plot percentages for.
    - ax (matplotlib.axes.Axes, optional): Matplotlib axis object to plot on.
    - title (str, optional): Title for the pie chart, supports LaTeX formatted strings.

    Returns:
    - matplotlib.axes.Axes: The axis with the pie chart.
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


def count_column_values(data: List[Dict[str, Any]], column: str) -> Dict[str, int]:
    """
    Count the occurrences of each unique value in the specified column from a list of
    dictionaries. Treats all data types, including lists, as single entities by converting
    lists to strings.

    Parameters:
    - data (List[Dict[str, Any]]): The data to process, where each item is a dictionary.
    - column (str): The key representing the column in each dictionary from which
    to count values.

    Returns:
    - Dict[str, int]: A dictionary where keys are unique values (converted to strings if
    lists) and values are the count of occurrences for each unique value.

    Raises:
    - KeyError: If the specified column is not found in any of the dictionaries.
    - ValueError: If the data list is empty or the column name is not provided.
    """
    if not data or not column:
        raise ValueError("Data list is empty or column name is not provided.")

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data)

    # Verify if the specified column exists
    if column not in df.columns:
        raise KeyError(f"Column '{column}' not found in data.")

    # Convert lists to strings if the column contains lists
    if df[column].apply(lambda x: isinstance(x, list)).any():
        df[column] = df[column].apply(lambda x: str(x) if isinstance(x, list) else x)

    # Count occurrences of each unique value
    return df[column].value_counts().to_dict()


def plot_rules_distribution(
    rules: Dict[str, int],
    rule_type: str = "single",
    ax: Optional[Axes] = None,
    title: Optional[str] = None,
    refinement: bool = False,
    threshold: float = 1,
    remove: bool = True,
    color_pallet: str = "pastel",
) -> None:
    """
    Plots the distribution of rules in a bar chart, optionally combining all entries under
    the threshold into a single category 'Under 1%' if `refinement` is True.

    Parameters:
    - rules (Dict[str, int]): Dictionary with rule counts keyed by rule name,
    where the values are counts.
    - rule_type (str, optional): Specifies the type of rules to plot
    ('single' or 'complex'). Default is 'single'.
    - ax (matplotlib.axes.Axes, optional): Matplotlib axis object to plot on.
    If None, a new figure is created.
    - title (str, optional): Optional title for the chart. If None,
    a default title based on `rule_type` is used.
    - refinement (bool, optional): If True, combines all percentages under
    the threshold into one category 'Under 1%'. Default is False.
    - threshold (float, optional): The percentage threshold under which all
    categories are combined into 'Under 1%' if `refinement` is True. Default is 1.
    - remove (bool, optional): If True, removes the last category from the plot.
    Default is True.
    - color_pallet (str, optional): Color palette to use for the plot.
    Default is 'pastel'.

    Returns:
    - None: The function directly modifies the `ax` object or creates a new plot.
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
            refined_rules["Under 1%"] = small_value_aggregate
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
    plt.rc("text", usetex=True)
    plt.rc("text.latex", preamble=r"\usepackage{amsmath}")  # Ensure amsmath is loaded

    # Create figure and axis if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6), dpi=120)

    # Plot the data
    sns.barplot(ax=ax, x=types_of_rules, y=percentages, palette=color_pallet)
    if title:
        ax.set_title(title, fontsize=24)
    else:
        ax.set_title(f"Distribution of {rule_type.capitalize()} Rules", fontsize=16)
    ax.set_xlabel("Cycle length", fontsize=18)
    ax.set_ylabel(r"Percentage (\%)", fontsize=18)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    # Set font size for x-tick and y-tick labels
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    # Add text labels above the bars
    for index, value in enumerate(percentages):
        ax.text(
            index, value + 0.5, f"{value:.1f}%", ha="center", va="bottom", fontsize=18
        )

    # Only show plot if ax is not provided (i.e., we created the figure here)
    if ax is None:
        plt.show()


def plot_bar_compare(
    data: List[Dict[str, Any]],
    ax: Axes,
    title: str = "Processing Time by Method and Radius",
    show_values: bool = True,
) -> None:
    """
    Plots a bar chart comparing processing times by method and radius on the provided Axes
    object.

    Paremeters:
    - data (List[Dict[str, Any]]): Data to plot, expected to have keys
    corresponding to radii ('R0', 'R1', etc.) and values as times.
    - ax (Axes): Matplotlib Axes object on which to plot the bar chart.
    - title (str, optional): Title for the plot.
    Defaults to "Processing Time by Method and Radius".
    - show_values (bool, optional): If True, values are displayed above the bars.
    Defaults to True.

    Returns:
    - None: This function modifies the `ax` object in-place and does not return any value.
    """
    plt.rc("text", usetex=True)  # Enable LaTeX rendering
    plt.rc("font", family="serif")  # Optional: set the default font family to serif

    types = [item["Type"] for item in data]
    R_labels = ["R0", "R1", "R2", "R3"]  # Simple labels to match dictionary keys
    formatted_R_labels = [
        r"$R_{0}$",
        r"$R_{1}$",
        r"$R_{2}$",
        r"$R_{3}$",
    ]  # LaTeX formatted labels for display
    colors = ["#3A8EBA", "#F4A582"]  # Reduced color palette
    n_types = len(types)
    index = np.arange(len(R_labels))  # One group per radius
    bar_width = 0.35  # Width of each bar

    # Prepare custom legend handles
    legend_handles = []

    # Loop through each type
    for i, typ in enumerate(types):
        means = [float(data[i][r]) for r in R_labels]
        bars = ax.bar(
            index + i * bar_width,
            means,
            bar_width,
            color=colors[i % len(colors)],
            label=typ,
        )

        # Add custom legend handle
        legend_handles.append(
            Patch(facecolor=colors[i % len(colors)], label=r"{" + typ + "}")
        )  # LaTeX in legend

        # Optionally add values on top of the bars
        if show_values:
            for bar in bars:
                height = bar.get_height()
                ax.annotate(
                    f"{height:.1f}",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=16,
                )

    # Configuring the axis
    ax.set_ylabel("Time (minutes)", fontsize=18)
    ax.set_title(r"{A. " + title + "}", fontsize=20, weight="medium")
    ax.set_xticks(index + bar_width / 2 * (n_types - 1))
    ax.set_xticklabels(
        formatted_R_labels, fontsize=16
    )  # Use LaTeX formatted labels for display
    ax.tick_params(axis="y", labelsize=16)
    ax.legend(
        handles=legend_handles, fontsize=16
    )  # Use custom legend handles with LaTeX formatting


def plot_cumulative_line(
    ax: Axes,
    radius: Sequence[Union[int, float]],
    hier: Sequence[Union[int, float]],
    emp: Sequence[Union[int, float]],
    colors: List[str] = ["#3A8EBA", "#F4A582"],
) -> None:
    """
    Plot cumulative line plots for hierarchical and empirical data sets on the provided
    Axes object.

    Parameters:
    - ax (Axes): The axes on which to plot.
    - radius (Sequence[Union[int, float]]):
    Radius values (e.g., time points or spatial measurements).
    - hier (Sequence[Union[int, float]]): Hierarchical data values
    over the specified radii.
    - emp (Sequence[Union[int, float]]): Empirical data values over the specified radii.
    - colors (List[str], optional): List of colors to use for each dataset.
    Defaults to ['#3A8EBA', '#F4A582'].

    Returns:
    - None: This function modifies the `ax` object in-place and does not return any value.
    """
    plt.rc("text", usetex=True)  # Enable LaTeX rendering
    plt.rc("font", family="serif")  # Optional: use serif font

    # Plot cumulative line plots with specified colors
    sns.lineplot(
        ax=ax,
        x=radius,
        y=np.cumsum(hier),
        color=colors[0],
        label="Hierarchical",
        linewidth=2,
    )
    sns.lineplot(
        ax=ax,
        x=radius,
        y=np.cumsum(emp),
        color=colors[1],
        label="Empirical",
        linewidth=2,
    )

    # Set axis labels and title with LaTeX formatting
    # ax.set_xlabel(r'\textbf{Radius}', fontsize=12)  # Bold label
    ax.set_ylabel(r"{Cumulative Time (minutes)}", fontsize=18)
    ax.set_title(r"{B. Cumulative Processing Time by Radius}", fontsize=20)

    # Set x-ticks and labels with LaTeX formatting
    ax.set_xticks(radius)
    ax.set_xticklabels([r"$R_{0}$", r"$R_{1}$", r"$R_{2}$", r"$R_{3}$"], fontsize=16)
    ax.tick_params(axis="y", labelsize=16)

    # Add legend with font settings
    ax.legend(fontsize=16)
