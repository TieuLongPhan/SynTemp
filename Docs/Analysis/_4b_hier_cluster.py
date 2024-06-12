import sys
import pathlib
import argparse
def main():
    root_dir = pathlib.Path(__file__).parents[2]
    sys.path.append(str(root_dir))
    from SynTemp.SynUtils.utils import load_from_pickle, save_to_pickle
    from SynTemp.SynRule.hierarchical_clustering import HierarchicalClustering

    parser = argparse.ArgumentParser(description="Run the Hierarchical Clustering pipeline.")
    parser.add_argument("--node_label_names", nargs='+', default=["element", "charge"], help="List of node labels in isomorphism check")
    parser.add_argument("--node_label_defaults", nargs='+', default=["*", 0], help="List of node labels defaults in isomorphism check")
    parser.add_argument("--edge_attribute", type=str, default="order", help="edge attribute using in isomorphism check")
    parser.add_argument("--max_radius", type=int, default=3, help="Maximun number radius extend for reaction center")
    parser.add_argument("--root_sample", type=int, default=100, help="Number of data in first process step")
    parser.add_argument("--folder_name", type=str, default="", help="Name of folder to store rules")
    parser.add_argument("--data_name", type=str, default="", help="Name of input data")

    args = parser.parse_args()
    data_path = f"{root_dir}/Data/DPO/{args.folder_name}/{args.data_name}"

    data = load_from_pickle(data_path)

    hcl = HierarchicalClustering(node_label_names=args.node_label_names,
                                node_label_default=args.node_label_defaults,
                                edge_attribute=args.edge_attribute,
                                max_radius=args.max_radius)

    reaction_dicts, templates = hcl.fit(data, 'ITSGraph', templates=None, update_template=True, root_sample=100)
    
    save_to_pickle(reaction_dicts, f"{root_dir}/Data/DPO/{args.folder_name}/cluster.pkl.gz")

    save_to_pickle(templates, f"{root_dir}/Data/DPO/{args.folder_name}/templates.pkl.gz")

if __name__ == "__main__":

    main()