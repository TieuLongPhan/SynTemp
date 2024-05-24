import pathlib
import sys
root_dir = pathlib.Path(__file__).parents[3]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import load_database, save_database
from synrbl import Balancer
list_data = ["ecoli", "recon3d", "uspto_3k", "golden", "natcomm"]

for data_name in list_data:
    save_path = f'{root_dir}/Data/AAM/balance/{data_name}/{data_name}_reactions.json.gz'
    data = load_database(f'{root_dir}/Data/AAM/balance/{data_name}_reactions.json.gz')
    synrbl = Balancer(reaction_col="reactions", id_col="R-id", n_jobs=2, confidence_threshold=0.5)
    balanced_reactions = synrbl.rebalance(
                reactions=data, output_dict=False
            )
    for i, new_reaction in enumerate(balanced_reactions):
        data[i]["old_reactions"] = data[i]["reactions"]
        data[i]["reactions"] = new_reaction
    save_database(data, save_path)
    