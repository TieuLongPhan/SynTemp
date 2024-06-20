from pathlib import Path
import argparse
import logging
from SynTemp.SynRule.rule_benchmark import RuleBenchmark
from SynTemp.SynUtils.utils import load_database, save_database
from SynTemp.SynChemistry.sf_similarity import SFSimilarity
from SynTemp.SynChemistry.sf_random import SFRandom
import pandas as pd


def setup_logging(log_dir, log_level):
    log_path = Path(log_dir)
    log_path.parent.mkdir(parents=True, exist_ok=True)  # Ensure directory exists
    
    # Create a custom logger
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    
    # Create handlers
    fh = logging.FileHandler(log_dir, mode='w')  # 'w' to overwrite existing logs
    fh.setLevel(getattr(logging, log_level))
    
    # Create formatters and add it to handlers
    f_format = logging.Formatter('%(asctime)s - %(message)s')
    fh.setFormatter(f_format)
    
    # Add handlers to the logger
    logger.addHandler(fh)
    
    return logger  # Return the logger instance if you want to use it later

def main(args):
    #root_dir = Path(__file__).parents[2]

    # Setup logging
    logger = setup_logging(args.log_dir, args.log_level)
    logger.info("Start process....")
    logger.info("Loading database....")
    # Load the database
    database = load_database(args.data_dir)[:]
    # Scoring functions dictionary
    scoring_functions = {
        'Random': SFRandom(),
        'ECFP6': SFSimilarity(["ECFP6"]),
        'MACCS': SFSimilarity(["MACCS"]),
        'RDK7': SFSimilarity(["RDK7"])
    }

    # Benchmarking setup
    top_k_values = [1, 3, 5, 10]

    batch_size = 500
    results_fw = []
    results_bw = []

    for i in range(0, len(database), batch_size):
        logger.info(f"Processing from {i} to {i+batch_size}")
        batch_reactions = database[i:i+batch_size]
        fw, bw = RuleBenchmark.reproduce_reactions(
        database=batch_reactions,
        rule_class_col=args.rule_class_col,
        rule_file_path=args.rule_file_path,
        original_rsmi_col=args.original_rsmi_col,
        repeat_times=1,
        use_specific_rules=args.use_specific_rules,
        verbosity=0,
        hierarchical=args.hierarchical,
        max_radius=args.max_radius,
        max_solutions = args.max_solutions)

        results_fw.extend(fw)
        results_bw.extend(bw)
    try:
        if args.hierarchical:
            save_database(results_fw, f'{args.save_dir}/fw_hier_{args.max_radius}.json.gz')
            save_database(results_bw, f'{args.save_dir}/bw_hier_{args.max_radius}.json.gz')
        else:
            save_database(results_fw, f'{args.save_dir}/fw_{args.radius}.json.gz')
            save_database(results_bw, f'{args.save_dir}/bw_{args.radius}.json.gz')
    except:
        logger.error('Cannot save')

    results_list_fw = []
    logger.info("Forward Prediction Validation")
    results_list_fw = []
    for name, func in scoring_functions.items():
        for k in top_k_values:
            accuracy = RuleBenchmark.TopKAccuracy(results_fw, args.original_rsmi_col, "ranked_reactions", k, ignore_stero=True, scoring_function=func)
            log_message = f"Top {k} accuracy for {name}: {accuracy}"
            logger.info(log_message)
            results_list_fw.append({'Scoring Function': name, 'Top K': f'Top {k}', 'Accuracy': accuracy})
    results_df_fw = pd.DataFrame(results_list_fw)

    pivot_df_fw = results_df_fw.pivot(index='Scoring Function', columns='Top K', values='Accuracy')

    logger.info("Backward Prediction Validation")
    results_list_bw = []
    for name, func in scoring_functions.items():
        for k in top_k_values:
            accuracy = RuleBenchmark.TopKAccuracy(results_bw, args.original_rsmi_col, "ranked_reactions", k, ignore_stero=True, scoring_function=func)
            log_message = f"Top {k} accuracy for {name}: {accuracy}"
            logger.info(log_message)
            results_list_bw.append({'Scoring Function': name, 'Top K': f'Top {k}', 'Accuracy': accuracy})
    
    results_df_bw = pd.DataFrame(results_list_bw)

    pivot_df_bw = results_df_bw.pivot(index='Scoring Function', columns='Top K', values='Accuracy')
    logger.info("Forward Prediction Results Matrix:\n" + pivot_df_fw.to_string())
    logger.info("Backward Prediction Results Matrix:\n" + pivot_df_bw.to_string())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run the synthesis rules benchmarking with configurable paths and options.')
    parser.add_argument('--log_dir', type=str, default='log.txt', help='File path for logging output.')
    parser.add_argument('--data_dir', type=str, default='data/test.json.gz', help='Path to the data directory containing the database file.')
    parser.add_argument('--rule_file_path', type=str, default='rules.json', help='Path to the JSON file containing synthesis rules.')
    parser.add_argument('--rule_class_col', type=str, default='R-id', help='Column name in the database that specifies the rule class identifier.')
    parser.add_argument('--original_rsmi_col', type=str, default='reactions', help='Column name in the database for the original SMILES notation of reactions.')
    parser.add_argument('--use_specific_rules', type=bool, default=False, help='Boolean to decide whether to use specific rules or not. Default is False.')
    parser.add_argument('--log_level', type=str, default='INFO', help='Logging level to use, options include DEBUG, INFO, WARNING, ERROR, CRITICAL.')
    parser.add_argument('--save_dir', type=str, default='data/', help='Path to the save data dir.')
    parser.add_argument('--radius', type=int, default=0, help='Radius of reaction center used' )
    parser.add_argument('--hierarchical', type=bool, default=False, help='Apply hierachical rules or not')
    parser.add_argument('--max_radius', type=int, default=3, help='Max radius of reaction center used' )
    parser.add_argument('--max_solutions', type=int, default=100, help='Maximun solution for one rule' )

    args = parser.parse_args()
    main(args)

# cd Documents/Project/TACsy/SynTemp/
# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Log/rules_r0_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Non_hydrogen/R0 --save_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Output --radius 0

# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Log/rules_r1_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Non_hydrogen/R1 --save_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Output --radius 1
    
# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Log/rules_r2_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Non_hydrogen/R2 --save_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Output --radius 2
    
# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Log/rules_r3_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Non_hydrogen/R3 --save_dir ./Data/DPO/USPTO_50K/Non_hydrogen/Output --radius 3

# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Hydrogen/Log/rules_r0_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Hydrogen/R0 --save_dir ./Data/DPO/USPTO_50K/Hydrogen/Output --radius 0

# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Hydrogen/Log/rules_r1_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Hydrogen/R1 --save_dir ./Data/DPO/USPTO_50K/Hydrogen/Output --radius 1
    
# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Log/rules_r0_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Good_hydrogen/R0 --save_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Output --radius 0

# python scripts.py --log_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Log/rules_r1_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Good_hydrogen/R1 --save_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Output --radius 1
    
# python scripts.py --log_dir Data/DPO/USPTO_50K/Good_hydrogen/Log/hier_rule_0_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Good_hydrogen/R0 --save_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Output --radius 0 --hierarchical True --max_radius 0 --max_solutions 1000
    
# python scripts.py --log_dir Data/DPO/USPTO_50K/Good_hydrogen/Log/hier_rule_1_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Good_hydrogen/R0 --save_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Output --radius 0 --hierarchical True --max_radius 1 --max_solutions 1000
    
# python scripts.py --log_dir Data/DPO/USPTO_50K/Good_hydrogen/Log/hier_rule_2_log.txt --data_dir ./Data/DPO/USPTO_50K/test.json.gz --rule_file_path ./Data/DPO/USPTO_50K/Good_hydrogen/R0 --save_dir ./Data/DPO/USPTO_50K/Good_hydrogen/Output --radius 0 --hierarchical True --max_radius 2 --max_solutions 1000
    
    