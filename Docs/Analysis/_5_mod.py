import sys
from pathlib import Path

root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import load_database
from SynTemp.SynRule.rule_benchmark import RuleBenchmark
from SynTemp.SynChemistry.sf_similarity import SFSimilarity
from SynTemp.SynChemistry.sf_maxfrag import SFMaxFrag
if __name__ == "__main__":
    import logging
    import pandas as pd
    folder_name = 'Hydrogen'

    # Set up logging
    logging.basicConfig(filename=f'{root_dir}/Data/DPO/USPTO_50K/Hydrogen/topk_accuracy_good.log', level=logging.INFO, format='%(asctime)s - %(message)s')

    # Load the database
    database = load_database(f"{root_dir}/Data/DPO/USPTO_50K/test.json.gz")[:10]

    # Set the parameters for the experiment
    top_k_values = [1, 3, 5, 10]
    scoring_functions = {
        'MaxFrag': SFMaxFrag(),
        'ECFP6': SFSimilarity(["ECFP6"]),
        'MACCS': SFSimilarity(["MACCS"]),
        'RDK7': SFSimilarity(["RDK7"])
    }

    # Prepare DataFrame to store results
    results_list = []

    # Run benchmark for each scoring function and Top K
    fw, bw = RuleBenchmark.reproduce_reactions(
        database=database,
        rule_class_col="R-id",
        rule_file_path=f"{root_dir}/Data/DPO/USPTO_50K/Hydrogen/Rules",
        original_rsmi_col="reactions",
        repeat_times=1,
        use_specific_rules=False,
    )

    for name, func in scoring_functions.items():
        for k in top_k_values:
            accuracy = RuleBenchmark.TopKAccuracy(
                fw,
                "reactions",
                "ranked_reactions",
                k,
                ignore_stero=True,
                scoring_function=func,
            )
            
            log_message = f"Top {k} accuracy for {name}: {accuracy}"
            logging.info(log_message)
            
            # Append results to the list
            results_list.append({'Scoring Function': name, 'Top K': f'Top {k}', 'Accuracy': accuracy})

    # Convert list to DataFrame
    results_df = pd.DataFrame(results_list)

    # Pivot the DataFrame to get the desired layout
    pivot_df = results_df.pivot(index='Scoring Function', columns='Top K', values='Accuracy')

    # # Save results to CSV
    # pivot_df.to_csv('topk_accuracy_matrix.csv')
    # logging.info("Results matrix saved to topk_accuracy_matrix.csv")

    # Log the pivot table
    pivot_log = pivot_df.to_string()
    logging.info("Results Matrix:\n" + pivot_log)
