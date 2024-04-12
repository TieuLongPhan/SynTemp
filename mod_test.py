# from mod import DG
# from mod import ruleGMLString, smiles, graphGMLString, addUniverse, addSubset, repeat
# from mod import *
# from SynTemp.SynMØD.rule_engine import RuleEngine
#from SynTemp.SynUtils.utils import load_database, save_database
from SynTemp.SynUtils.utils import run_shell_command, ensure_directory_exists

from SynTemp.SynUtils.chemutils import standardize_rsmi
from SynTemp.SynRule.rule_executor import RuleExecutor
ensure_directory_exists('./out')
test = RuleExecutor.reaction_prediction(input_smiles=['C=C1C(=C)C2OC1C1=C2CC(C(C)=O)CC1'],
                                        rule_file_path='./Data/uspto/Rule/USPTO_50K_31.gml',
                                        prediction_type='backward', repeat_times=1, print_results=True)

