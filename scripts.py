from SynTemp.SynUtils.utils import load_database, save_database
from SynTemp.SynChemistry.sf_similarity import SFSimilarity
from SynTemp.SynChemistry.sf_maxfrag import SFMaxFrag
from SynTemp.SynRule.rule_benchmark import RuleBenchmark
import sys
import os
import argparse
import logging
import pandas as pd
from pathlib import Path
def setup_logging(log_dir, log_level):
    logging.basicConfig(filename=log_dir,
                        level=getattr(logging, log_level),
                        format='%(asctime)s - %(message)s')

def main(args):
    #root_dir = Path(__file__).parents[2]
    logging.info("Start process....")
    # Setup logging
    setup_logging(args.log_dir, args.log_level)
    logging.info("Loading database....")
    # Load the database
    database = load_database(args.data_dir)[:]
    folder_path = os.path.dirname(args.data_dir)
    # Scoring functions dictionary
    scoring_functions = {
        'MaxFrag': SFMaxFrag(),
        'ECFP6': SFSimilarity(["ECFP6"]),
        'MACCS': SFSimilarity(["MACCS"]),
        'RDK7': SFSimilarity(["RDK7"])
    }

    # Benchmarking setup
    top_k_values = [1, 3, 5, 10]
    
    fw, bw = RuleBenchmark.reproduce_reactions(
        database=database,
        rule_class_col=args.rule_class_col,
        rule_file_path=args.rule_file_path,
        original_rsmi_col=args.original_rsmi_col,
        repeat_times=1,
        use_specific_rules=args.use_specific_rules,
        verbosity=0
    )
    try:
        save_database(fw, f'{folder_path}/fw_good.json.gz')
        save_database(bw, f'{folder_path}/bw_good.json.gz')
    except:
        logging.error('Cannot save')

    results_list_fw = []
    logging.info("Forward Prediction Validation")
    results_list_fw = []
    for name, func in scoring_functions.items():
        for k in top_k_values:
            accuracy = RuleBenchmark.TopKAccuracy(fw, args.original_rsmi_col, "ranked_reactions", k, ignore_stero=True, scoring_function=func)
            log_message = f"Top {k} accuracy for {name}: {accuracy}"
            logging.info(log_message)
            results_list_fw.append({'Scoring Function': name, 'Top K': f'Top {k}', 'Accuracy': accuracy})
    results_df_fw = pd.DataFrame(results_list_fw)

    pivot_df_fw = results_df_fw.pivot(index='Scoring Function', columns='Top K', values='Accuracy')

    logging.info("Backward Prediction Validation")
    results_list_bw = []
    for name, func in scoring_functions.items():
        for k in top_k_values:
            accuracy = RuleBenchmark.TopKAccuracy(bw, args.original_rsmi_col, "ranked_reactions", k, ignore_stero=True, scoring_function=func)
            log_message = f"Top {k} accuracy for {name}: {accuracy}"
            logging.info(log_message)
            results_list_bw.append({'Scoring Function': name, 'Top K': f'Top {k}', 'Accuracy': accuracy})
    
    results_df_bw = pd.DataFrame(results_list_bw)

    pivot_df_bw = results_df_bw.pivot(index='Scoring Function', columns='Top K', values='Accuracy')
    logging.info("Forward Prediction Results Matrix:\n" + pivot_df_fw.to_string())
    logging.info("Backward Prediction Results Matrix:\n" + pivot_df_bw.to_string())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the synthesis rules benchmarking with configurable paths and options.')
    parser.add_argument('--log_dir', type=str, default='log.txt', help='File path for logging output.')
    parser.add_argument('--data_dir', type=str, default='data/test.json.gz', help='Path to the data directory containing the database file.')
    parser.add_argument('--rule_file_path', type=str, default='rules.json', help='Path to the JSON file containing synthesis rules.')
    parser.add_argument('--rule_class_col', type=str, default='R-id', help='Column name in the database that specifies the rule class identifier.')
    parser.add_argument('--original_rsmi_col', type=str, default='reactions', help='Column name in the database for the original SMILES notation of reactions.')
    parser.add_argument('--use_specific_rules', type=bool, default=False, help='Boolean to decide whether to use specific rules or not. Default is False.')
    parser.add_argument('--log_level', type=str, default='INFO', help='Logging level to use, options include DEBUG, INFO, WARNING, ERROR, CRITICAL.')
    args = parser.parse_args()
    main(args)


# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Hydrogen/rules_good_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Hydrogen/Rules_good