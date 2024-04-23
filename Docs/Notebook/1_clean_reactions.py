import pathlib
import sys
import logging
import time
root_dir = pathlib.Path(__file__).parents[1]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import load_database, save_database
from SynTemp.SynChemistry.reaction_processor import ReactionProcessor
from SynTemp.SynChemistry.unbalanced_charge import UnbalancedCharge
from SynTemp.SynChemistry.uncharge_reaction import UnchargeReaction
import warnings
warnings.filterwarnings("ignore")

def main():
    rule_balanced = load_database(f'{root_dir}/Data/rule_based_reactions.json.gz')
    rule_balanced = [{'R-id': d['R-id'], 'new_reaction': d['new_reaction']} for d in rule_balanced if 'R-id' in d and 'new_reaction' in d]

    labeled_list = ReactionProcessor.process_reactions_parallel(rule_balanced, n_jobs=4)
    label_charge_list = ReactionProcessor.sum_of_charge_in_products(labeled_list, n_jobs=4)

    
    fix_charge_balance = UnbalancedCharge.parallel_fix_unbalanced_charge(label_charge_list, charges_column='total_charge_in_products', n_jobs=4)

    uncharge_fix = UnchargeReaction()
    fix_charge_balance_post = uncharge_fix.apply_uncharge_smiles_to_reactions(fix_charge_balance, uncharge_fix.uncharge_smiles, n_jobs=4)

    save_database(fix_charge_balance_post, f'{root_dir}/Data/clean_reactions.json.gz')

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    start_time = time.time()  
    main()
    end_time = time.time() 

    elapsed_time = end_time - start_time
    logging.info(f"Execution time: {elapsed_time:.2f} seconds")