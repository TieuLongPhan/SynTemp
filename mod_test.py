from mod import DG
from mod import ruleGMLString, smiles, graphGMLString, addUniverse, addSubset, repeat
from mod import *
from SynTemp.SynUtils.utils import load_database, save_database

from SynTemp.SynMØD.MØD_modeling import MØDModeling

database = load_database('./test_database.json.gz')

for key, value in enumerate(database):
    rule_name = value['R-id']
    rule_file = f'./Data/uspto/Rule/{rule_name}.gml'
    initial_smiles = value['reactions'].split('>>')[0].split('.')
    reactions = MØDModeling.perform_reaction(rule_file_path=rule_file, invert_rule=False, initial_smiles=initial_smiles, return_type='reaction',
                                 repeat_times=1)
    database[key]['new_reactions'] = reactions

save_database(database, './test_database_out.json.gz')
