from typing import List, Dict, Optional, Any
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

try:
    from umap import UMAP  # UMAP might not be available in all environments
except ImportError:
    UMAP = None
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from joblib import Parallel, delayed


class ClusterVisualization:
    def __init__(
        self,
        data: List[Dict],
        fps_col: str = "fps",
        cluster_col: str = "Cluster",
        n_jobs: int = 4,
        verbose: int = 1,
    ):
        """
        Initialize the ClusterVisualization class with the dataset and relevant column names.

        Parameters:
        - data (List[Dict]): List of dictionaries containing molecule data, including fingerprints and cluster assignments.
        - fps_col (str): The key in the dictionary for accessing the fingerprint data.
        - cluster_col (str): The key in the dictionary for accessing the cluster assignment data.
        - n_jobs (int): The number of jobs to run in parallel for fingerprint conversion.
        - verbose (int): The verbosity level of job execution.
        """
        self.data = data
        self.fps_col = fps_col
        self.cluster_col = cluster_col
        self.n_jobs = n_jobs
        self.verbose = verbose

    @staticmethod
    def fps_to_array(fp: ExplicitBitVect) -> np.ndarray:
        """
        Converts a fingerprint to a numpy array.

        Parameters:
        - fp (ExplicitBitVect): The fingerprint to convert.

        Returns:
        - np.ndarray: The fingerprint represented as a numpy array.
        """
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr

    def dimension_reduction(
        self, algorithm: str = "pca", combine: bool = False, **kwargs: Any
    ) -> np.ndarray:
        """
        Performs dimensionality reduction on the fingerprint data.
        If combine is True, PCA is used to reduce to 50 components before applying t-SNE or UMAP.

        Parameters:
        - algorithm (str): The algorithm to use for dimensionality reduction ('pca', 'tsne', 'umap').
        - combine (bool): Whether to combine PCA and t-SNE/UMAP for dimensionality reduction.
        - **kwargs: Additional keyword arguments for the dimensionality reduction algorithms.

        Returns:
        - np.ndarray: The dimensionality-reduced components.
        """
        fps = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
            delayed(self.fps_to_array)(d[self.fps_col]) for d in self.data
        )
        arr = np.array(fps)

        if combine:
            pca = PCA(n_components=50)
            pca_components = pca.fit_transform(arr)
            if algorithm == "tsne":
                reducer = TSNE(n_components=2, **kwargs)
            elif algorithm == "umap" and UMAP is not None:
                reducer = UMAP(n_components=2, **kwargs)
            else:
                raise ValueError(
                    "Combine option is only supported with 'tsne' or 'umap'."
                )
            components = reducer.fit_transform(pca_components)
        else:
            if algorithm == "pca":
                reducer = PCA(n_components=2, **kwargs)
            elif algorithm == "tsne":
                reducer = TSNE(n_components=2, **kwargs)
            elif algorithm == "umap" and UMAP is not None:
                reducer = UMAP(n_components=2, **kwargs)
            else:
                raise ValueError(f"Unsupported algorithm: {algorithm}")
            components = reducer.fit_transform(arr)
        return components

    def visualize(
        self,
        components: Optional[np.ndarray] = None,
        algorithm: str = "pca",
        combine: bool = False,
        raw_data: Optional[np.ndarray] = None,
        figsize: tuple = (10, 8),
        title_fontsize: int = 16,
        label_fontsize: int = 14,
        legend_fontsize: int = 12,
        color_palette: str = "husl",
        custom_title: Optional[str] = None,
        **kwargs: Any,
    ) -> None:
        """
        Visualizes the clustered data using dimensionality reduction for layout with enhanced plotting features.

        This method supports visualization enhancements such as placing the legend outside the plot and
        setting the plot title to semi-bold. It integrates Seaborn for aesthetic improvements and allows
        dynamic dimensionality reduction if raw data is provided.

        Parameters:
        - components (Optional[np.ndarray]): Pre-computed components for visualization. Required if raw_data is None.
        - algorithm (str): Dimensionality reduction algorithm ('pca', 'tsne', 'umap'). 'umap' requires UMAP to be installed.
        - combine (bool): If True, combines PCA (50 components) with 'tsne' or 'umap' for dimensionality reduction.
        - raw_data (Optional[np.ndarray]): Raw data for on-the-fly dimensionality reduction. Required if components is None.
        - figsize (tuple): Figure size for the plot.
        - title_fontsize (int): Font size for the plot title.
        - label_fontsize (int): Font size for the axis labels.
        - legend_fontsize (int): Font size for the legend.
        - color_palette (str): Seaborn color palette for the plot.
        - custom_title (Optional[str]): Custom title for the plot. Auto-generated based on the algorithm if None.
        - **kwargs: Additional keyword arguments passed to the dimensionality reduction algorithm or Seaborn scatterplot.
        """
        sns.set_style("whitegrid")
        if components is None and raw_data is None:
            components = self.dimension_reduction(
                algorithm=algorithm, combine=combine, **kwargs
            )

        self.df_cls = pd.DataFrame(components, columns=["x", "y"])
        self.df_cls[self.cluster_col] = [d[self.cluster_col] for d in self.data]

        title = custom_title if custom_title else f"{algorithm.upper()} Visualization"
        if combine:
            title = f"PCA -> {algorithm.upper()} Combined Visualization"

        plt.figure(figsize=figsize)
        plot = sns.scatterplot(
            x="x",
            y="y",
            hue=self.cluster_col,
            palette=sns.color_palette(color_palette, as_cmap=True),
            data=self.df_cls,
            legend="full",
            **kwargs,
        )
        plt.title(title, fontsize=title_fontsize, weight="bold")
        plot.set(xlabel="Component 1", ylabel="Component 2")

        plot.title.set_fontsize(title_fontsize)
        plot.xaxis.label.set_fontsize(label_fontsize)
        plot.yaxis.label.set_fontsize(label_fontsize)
        plt.legend(
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            fontsize=legend_fontsize,
            title="Cluster",
            borderaxespad=0.0,
        )

        plt.tight_layout()
        plt.show()
