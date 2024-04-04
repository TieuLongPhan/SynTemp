import pathlib
import sys
import logging
import torch
import time
root_dir = pathlib.Path(__file__).parents[1]
sys.path.append(str(root_dir))
from SynTemp.SynProcessor.balance_checker import BalanceReactionCheck
from SynTemp.SynAAM.consensus_aam import ConsensusAAM
from SynTemp.utils import load_database, save_database
import warnings
warnings.filterwarnings("ignore")
import transformers
from rxnmapper import RXNMapper
from localmapper import localmapper
transformers.logging.set_verbosity_error()

def main(data, save_dir=None, data_name = '', batch_size=1000, check_balance=True):
    rxn_mapper = RXNMapper()
    #local_mapper = AtomMapper(device=torch.device('cpu') , dataset='USPTO_FULL', root_dir=f'{root_dir}/SynITSG/LocalMapper')
    
    if check_balance:
        checker = BalanceReactionCheck(data, rsmi_column='reactions', n_jobs=5, verbose=2)
        balanced_reactions, _ = checker.check_balances()
    else:
        balanced_reactions = data
    
    consensus_aam = ConsensusAAM(balanced_reactions, rsmi_column='reactions', save_dir=f'{root_dir}/Data', 
                             mapper_types=['rxn_mapper', 'graphormer', 'local_mapper'])
       
    mapped_reactions = consensus_aam.run(batch_size, rxn_mapper)
    if save_dir:
        save_database(mapped_reactions, f'{save_dir}/{data_name}_aam_reactions.json.gz')
  
    return mapped_reactions

if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO)
    #folder_names = ['uspto', 'jaworski', 'golden', 'ecoli']
    folder_name = 'jaworski'
    start_time = time.time()  
    save_dir = f'{root_dir}/Data/{folder_name}'
    data = load_database(f'{save_dir}/{folder_name}_reactions.json.gz')
    print(len(data))
    mapped_reactions = main(data, save_dir=save_dir, data_name = folder_name, batch_size=50, check_balance=False)
    end_time = time.time() 

    elapsed_time = end_time - start_time
    logging.info(f"Execution time: {elapsed_time:.2f} seconds")