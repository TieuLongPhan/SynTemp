from SynTemp.auto_template import AutoTemp
from SynTemp.SynUtils.utils import load_from_pickle

data = load_from_pickle('Data/Testcase/its_graph.pkl.gz')
results = AutoTemp(data, save_path='Data/Testcase')
print(results)