import sys
from pathlib import Path
root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynUtils.utils import load_database
from SynTemp.SynRule.rule_benchmark import RuleBenchmark
from SynTemp.SynChemistry.sf_similarity import SFSimilarity
from SynTemp.SynChemistry.sf_maxfrag import SFMaxFrag

if __name__ == "__main__": 
    database = load_database(f'{root_dir}/Data/uspto/demo_database.json.gz')
    fw, bw = RuleBenchmark.reproduce_reactions(database=database,  id_col='R-id', rule_file_path=f'{root_dir}/Data/uspto/Rule',
                                            original_rsmi_col='reactions', repeat_times=1, prior=True)



    print("Top 5 accuracy-MaxFrag:", RuleBenchmark.TopKAccuracy(fw, 'reactions','ranked_reactions', 5, ignore_stero=True, scoring_function=SFMaxFrag()))
    print("Top 5 accuracy-ECFP6:", RuleBenchmark.TopKAccuracy(fw, 'reactions','ranked_reactions', 5, ignore_stero=True, scoring_function=SFSimilarity(['ECFP6'])))                                                    