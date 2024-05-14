import sys
from pathlib import Path
root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import load_database, save_database
from SynTemp.SynRule.rdecoy import RDecoy


if __name__ == "__main__": 
    
    data = load_database(f'{root_dir}/Data/USPTO_50K/clean_cluster_uspto.json.gz')
    data = [value for value in data if value['R-id'] !=2027] #issue data
    fw_decoy = RDecoy.generate_reactions(database=data, rule_file_path=f'{root_dir}/Data/USPTO_50K/Rule',
                                        original_rsmi_col='reactions', repeat_times=1, cluster_col='Cluster')
    save_database(fw_decoy, f'{root_dir}/Data/USPTO_50K/rdecoy.json.gz')