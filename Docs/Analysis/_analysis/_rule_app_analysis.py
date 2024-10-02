import json
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, Optional, Dict, List, Any


plt.rc("text", usetex=True)
plt.rc("text.latex", preamble=r"\usepackage{amsmath}")


def save_results_to_json(results: Dict, file_path: str) -> None:
    """
    Save the given results dictionary to a JSON file specified by the file path.

    Parameters:
        results (Dict): The results dictionary to be saved.
        file_path (str): The path to the JSON file where the results will be saved.

    Raises:
        IOError: If there is an issue writing to the file.
    """
    try:
        with open(file_path, "w", encoding="utf-8") as file:
            json.dump(results, file, ensure_ascii=False, indent=4)
        print(f"Results successfully saved to {file_path}")
    except IOError as e:
        print(f"Failed to write to file: {e}")


def load_results_from_json(file_path: str) -> Dict:
    """
    Load and return the results dictionary from a specified JSON file.

    Parameters:
        file_path (str): The path to the JSON file from which results will be loaded.

    Returns:
        Dict: The dictionary containing the loaded results.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file is not a valid JSON.
    """
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            data = json.load(file)
        print(f"Results successfully loaded from {file_path}")
        return data
    except FileNotFoundError:
        print(f"No file found at the specified path: {file_path}")
        raise
    except json.JSONDecodeError:
        print(f"Error decoding JSON from the file at {file_path}")
        raise


def load_database(pathname: str = "./Data/database.json") -> List[Dict]:
    """
    Load a database (a list of dictionaries) from a JSON file.

    Args:
        pathname: The path from where the database will be loaded.
                    Defaults to './Data/database.json'.

    Returns:
        The loaded database.

    Raises:
        ValueError: If there is an error reading the file.
    """
    try:
        with open(pathname, "r") as f:
            database = json.load(f)  # Load the JSON data from the file
        return database
    except IOError as e:
        raise ValueError(f"Error reading to file {pathname}: {e}")


def coverage_rate(
    data: List[Dict],
    solutions_column: str = "unrank",
    positive_columns: str = "positive_reactions",
) -> Tuple[float, float, float]:
    """
    Calculate the average solutions count, coverage rate, and false positive rate for
    solutions in the data.

    Parameters:
    - data: List of dictionaries containing the data.
    - solutions_column: The key in the dictionaries representing the solutions.
    Default is 'unrank'.
    - positive_columns: The key in the dictionaries representing the positive reactions.
    Default is 'positive_reactions'.

    Returns:
    - A tuple containing:
      - The average number of solutions per entry.
      - The coverage rate (percentage): the number of entries with positive reactions
      divided by the total number of entries.
      - The average false positive rate: the average ratio of false solutions per solution
      when at least one solution is present.
    """
    if not data:
        return 0.0, 0.0, 0.0

    total_solutions = []
    positive_counts = 0
    false_positive_ratios = []

    for entry in data:
        try:
            solutions = entry.get(solutions_column, [])
            positive = entry.get(positive_columns, False)

            if not isinstance(solutions, list):
                continue  # Skip entries with incorrect format

            num_solutions = len(solutions)
            total_solutions.append(num_solutions)

            if positive:
                positive_counts += 1
                if num_solutions > 0:
                    false_positive_ratios.append((num_solutions - 1) / num_solutions)
            else:
                if num_solutions > 0:
                    false_positive_ratios.append(1)

        except KeyError as e:
            print(f"Key error: {e} in entry: {entry}")

    total_entries = len(data)

    average_solutions = np.mean(total_solutions) if total_solutions else 0.0
    coverage_rate = (
        (positive_counts / total_entries) * 100 if total_entries > 0 else 0.0
    )
    average_fpr = np.mean(false_positive_ratios) * 100 if false_positive_ratios else 0.0

    return round(average_solutions, 2), round(coverage_rate, 2), round(average_fpr, 2)


def automatic_results(
    test_types: List[str],
    temp_types: List[str],
    predict_types: List[str],
    radii: List[int],
    base_path="../../Data/Temp/Benchmark",
) -> Dict[str, Dict[str, Tuple[float, float, float]]]:
    """
    Automatically computes coverage rates for combinations of test type, template type,
    predict type, and radii. Iterates over the provided parameter lists, loads data,
    and computes statistics.

    Parameters:
    - test_types (List[str]): List of test types.
    - temp_types (List[str]): List of template types.
    - predict_types (List[str]): List of prediction types.
    - radii (List[int]): List of radii values.
    - base_path (str): path to data

    Returns:
    - Dict[str, Dict[str, Tuple[float, float, float]]]: A dictionary where the key
    is the test type and the value is another dictionary. The inner dictionary's keys are
    combinations of parameters as strings, and its values are tuples with the results from
    `coverage_rate` (average solutions, coverage rate, false positive rate).
    """
    all_results = {}

    for test in test_types:
        test_results = {}
        for predict in predict_types:
            predict_results = {}
            for temp in temp_types:
                for rad in radii:
                    path = f"{base_path}/{temp}/Output/{test}/{predict}_{rad}.json.gz"
                    name = f"{temp}_{rad}"
                    data = load_database(path)
                    if data:
                        predict_results[name] = coverage_rate(data)
                    else:
                        predict_results[name] = (0.0, 0.0, 0.0)
            test_results[predict] = predict_results
        all_results[test] = test_results

    return all_results


def plot_percentage(
    df: pd.DataFrame,
    ax: plt.Axes,
    column: str,
    title: str = "A",
    color_map: Optional[List[str]] = None,
    fontsettings: Optional[Dict[str, int]] = None,
) -> None:
    """
    Plot a percentage bar chart for different categories and subcategories within the data.

    Parameters:
    df (pd.DataFrame): DataFrame containing the data to plot. Index of the DataFrame
                       should be string labels in the format 'category_subcategory'.
    ax (plt.Axes): Matplotlib Axes object where the chart will be drawn.
    column (str): Column name in df that contains the percentage values to plot.
    title (str, optional): Title of the plot. Default is 'A'.
    color_map (List[str], optional): List of hex color strings for the bars. If None,
                                     a default set of colors will be used.
    fontsettings (Dict[str, int], optional): Dictionary containing font size settings
                                             for various elements of the plot. If None,
                                             default settings are applied.

    Returns:
    None: This function does not return any value but modifies the ax object by drawing a bar chart.

    Example:
    >>> fig, ax = plt.subplots()
    >>> data = pd.DataFrame({'Value': [20, 30, 40, 50]}, index=['Type1_10', 'Type1_20', 'Type2_10', 'Type2_20'])
    >>> plot_percentage(data, ax, 'Value')
    >>> plt.show()
    """
    if fontsettings is None:
        fontsettings = {
            "title_size": 18,
            "label_size": 16,
            "ticks_size": 16,
            "annotation_size": 12,
        }

    # Split the index into template type and radii
    df["Type"] = [i.split("_")[0] for i in df.index]
    df["Radii"] = [int(i.split("_")[1]) for i in df.index]

    # Sort data to group by type and then by radii
    df = df.sort_values(by=["Radii"])

    # Prepare color map for radii using coolwarm
    if color_map is None:
        color_map = ["#3A8EBA", "#92C5DE", "#F4A582", "#D6604D"]

    # Plotting logic with annotations
    total_width = 3  # Total width for group
    width = total_width / len(
        df["Radii"].unique()
    )  # Width for each bar within each type group
    type_positions = np.arange(len(df["Type"].unique())) * (
        len(df["Radii"].unique()) + 1
    )

    for i, t in enumerate(df["Type"].unique()):
        for j, r in enumerate(df["Radii"].unique()):
            # print(t)
            bar_positions = type_positions[i] + j * width
            heights = df[(df["Type"] == t) & (df["Radii"] == r)][column]
            ax.bar(
                bar_positions,
                heights,
                width=width,
                label=f"$R_{{{r}}}$" if i == 0 else "",
                color=color_map[j % len(color_map)],
            )
    # Adding annotations
    for rect in ax.patches:
        height = rect.get_height()
        ax.annotate(
            f"{height:.1f}%",
            xy=(rect.get_x() + rect.get_width() / 2, height),
            xytext=(0, 3),  # 3 points vertical offset
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsettings["annotation_size"],
        )

    # Enhancements like axes labeling, ticks setting, and adding grid
    ax.set_ylabel(rf"$\mathcal{{{column}}} (\%)$", fontsize=fontsettings["label_size"])
    ax.set_title(title, fontsize=fontsettings["title_size"], weight="medium")
    ax.set_xticks(type_positions + total_width / 2 - width / 2)
    ax.set_xticklabels(
        [f"$Q_{{\\text{{{t}}}}}$" for t in df["Type"].unique()],
        fontsize=fontsettings["ticks_size"],
    )
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_yticklabels(
        [f"{i}%" for i in range(0, 101, 20)], fontsize=fontsettings["ticks_size"]
    )
    ax.grid(True, which="major", linestyle="--", linewidth="0.5", color="grey")
    ax.set_axisbelow(True)


def plot_roc_curves(
    df: pd.DataFrame,
    ax: plt.Axes,
    selected_types: Optional[List[str]] = None,
    fontsettings: Optional[Dict[str, int]] = None,
    title: str = "A",
    add_legend: bool = False,
) -> List[Any]:
    """
    Plot ROC curves for specified types from a DataFrame on a given matplotlib Axes.

    Parameters:
    df (pd.DataFrame): DataFrame containing the data for plotting. Must include columns 'Type', 'C' for TPR,
                       and 'FPR' for FPR, where 'Type' differentiates data series.
    ax (plt.Axes): The matplotlib Axes object where the ROC curves will be drawn.
    selected_types (Optional[List[str]]): List of strings representing the types to be included in the plot.
                                          If None, all types in the DataFrame will be plotted.
    fontsettings (Optional[Dict[str, int]]): Dictionary containing font settings for titles, labels,
                                             ticks, and annotations. If None, defaults will be applied.
    title (str): Title of the plot.
    add_legend (bool): If True, add a legend to the plot.

    Returns:
    List[Any]: List containing matplotlib line handles for the legend, useful if further customization
               or reference is needed.

    Raises:
    ValueError: If selected_types is provided and contains non-string elements.

    Example:
    >>> fig, ax = plt.subplots()
    >>> data = pd.DataFrame({
    ...     'Type': ['Type1', 'Type1', 'Type2', 'Type2'],
    ...     'C': [90, 85, 88, 80],
    ...     'FPR': [5, 10, 5, 10]
    ... })
    >>> plot_roc_curves(data, ax, ['Type1', 'Type2'])
    >>> plt.show()
    """
    if selected_types is not None:
        if not all(isinstance(t, str) for t in selected_types):
            raise ValueError("selected_types must be a list of strings.")
        original_types = [t for t in selected_types if t in df["Type"].unique()]
    else:
        original_types = df["Type"].unique()

    types = [f"$Q_{{\\text{{{t}}}}}$" for t in original_types]

    if fontsettings is None:
        fontsettings = {
            "title_size": 28,
            "label_size": 24,
            "ticks_size": 24,
            "annotation_size": 18,
        }

    markers = ["o", "^", "s", "p"]
    markers.reverse()
    marker_labels = [r"$R_{0}$", r"$R_{1}$", r"$R_{2}$", r"$R_{3}$"]
    marker_labels.reverse()
    marker_color = "gray"

    colors = plt.cm.coolwarm(np.linspace(0, 1, len(types)))
    colors = ["#3A8EBA", "#D6604D"]
    legend_handles = []

    for index, type_ in enumerate(original_types):
        type_data = df[df["Type"] == type_]
        tpr = type_data["C"].tolist()
        fpr = type_data["NR"].tolist()
        tpr = [x / 100 for x in tpr]
        fpr = [x / 100 for x in fpr]
        tpr.reverse()
        fpr.reverse()

        (line,) = ax.plot(
            fpr, tpr, linestyle="-", color=colors[index], label=f"{types[index]}"
        )
        legend_handles.append(line)

        for i, (f, t) in enumerate(zip(fpr, tpr)):
            marker = ax.plot(
                f, t, marker=markers[i % len(markers)], color=marker_color
            )[0]
            if index == 1:
                marker_handle = plt.Line2D(
                    [0],
                    [0],
                    marker=markers[i % len(markers)],
                    color="none",
                    markerfacecolor=marker_color,
                    markersize=10,
                    label=marker_labels[i],
                )
                legend_handles.append(marker_handle)

    ax.set_xlabel(r"$\mathcal{NR}\ (\%)$", fontsize=fontsettings["label_size"])
    ax.set_ylabel(r"$\mathcal{C}\ (\%)$", fontsize=fontsettings["label_size"])
    ax.set_title(rf"{title}", fontsize=fontsettings["title_size"], weight="medium")
    ax.tick_params(axis="both", which="major", labelsize=fontsettings["ticks_size"])
    ax.grid(True)

    if add_legend:
        ax.legend(
            handles=legend_handles,
            loc="lower right",
            fancybox=True,
            title_fontsize=fontsettings["label_size"],
            fontsize=fontsettings["annotation_size"],
            ncol=3,
        )

    ax.grid(True)
    return legend_handles


def plot_processing_times(
    times: Dict[str, List[float]], ax: Optional[plt.Axes] = None, title: str = "A"
) -> None:
    """
    Plot processing times for various methods across different stages.

    This function takes a dictionary of processing times, converts them into hours,
    and plots them using a bar chart.

    Parameters:
    times (Dict[str, List[float]]): A dictionary where keys are method names and values
                                    are lists of processing times in seconds for each stage.
    ax (Optional[plt.Axes]): Matplotlib Axes object where the plot will be drawn. If None,
                             the current active Axes will be used.
    title (str): The title of the plot.

    Returns:
    None: The function creates a plot but does not return any value.

    Example:
    >>> times = {
    ...     "Method1": [3600, 7200, 1800, 5400],
    ...     "Method2": [1800, 3600, 900, 2700],
    ... }
    >>> fig, ax = plt.subplots()
    >>> plot_processing_times(times, ax=ax, title="Processing Times Analysis")
    >>> plt.show()
    """
    # Convert to hours
    for key in times:
        times[key] = np.array(times[key]) / 3600

    # Stages
    stages = [r"$R_{0}$", r"$R_{1}$", r"$R_{2}$", r"$R_{3}$"]

    # Create a DataFrame
    df = (
        pd.DataFrame(times, index=stages)
        .reset_index()
        .melt(id_vars="index", var_name="Method", value_name="Time (hours)")
    )
    df.rename(columns={"index": "Stage"}, inplace=True)

    # Create the plot on the provided ax
    if ax is None:
        ax = plt.gca()  # Get current axis if not provided

    custom_colors = ["#5e4fa2", "#3A8EBA", "#D6604D"]
    palette = sns.color_palette(custom_colors[: len(times.keys())])
    bar_plot = sns.barplot(
        x="Stage", y="Time (hours)", hue="Method", data=df, palette=palette, ax=ax
    )

    ax.set_title(title, fontsize=24, weight="bold")
    ax.set_xlabel(None)
    ax.set_ylabel("Time (Hours)", fontsize=20)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=20)
    ax.set_yticklabels([f"{y:.0f}" for y in ax.get_yticks()], fontsize=20)
    ax.legend(
        title="Method",
        title_fontsize="20",
        fontsize="18",
        loc="upper left",
        bbox_to_anchor=(0.01, 1),
    )

    # Add text annotations on the bars
    for p in bar_plot.patches:
        bar_height = p.get_height()
        if bar_height > 0.01:  # Adjust this threshold as needed
            annotation = format(
                p.get_height(), ".1f" if p.get_height() < 100 else ".0f"
            )
            ax.annotate(
                annotation,
                (p.get_x() + p.get_width() / 2, p.get_height()),
                ha="center",
                va="center",
                xytext=(0, 9),
                textcoords="offset points",
                fontsize=16,
            )

    ax.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
