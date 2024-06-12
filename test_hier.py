from SynTemp.SynRule.rule_benchmark import RuleBenchmark
from SynTemp.SynUtils.utils import load_database, save_to_pickle
data = load_database('./Data/DPO/USPTO_50K/test.json.gz')


# start_index = 3
# end_index = 4

fw_hier, bw_hier = RuleBenchmark.reproduce_reactions(
        database=data[0:100],
        rule_class_col='R-id',
        rule_file_path='./Data/DPO/USPTO_50K/Good_hydrogen/R0',
        original_rsmi_col='reactions',
        repeat_times=1,
        use_specific_rules=False,
        verbosity=0,
        job_timeout=5,
        hierarchical=True,
        max_radius=3,
        max_solutions = 10
    )
save_to_pickle(fw_hier, './fw_hier.pkl.gz')
save_to_pickle(bw_hier, './bw_hier.pkl.gz')

# fw, bw = RuleBenchmark.reproduce_reactions(
#         database=data[start_index:end_index],
#         rule_class_col='R-id',
#         rule_file_path='./Data/DPO/USPTO_50K/Good_hydrogen/R0',
#         original_rsmi_col='reactions',
#         repeat_times=1,
#         use_specific_rules=False,
#         verbosity=0,
#         job_timeout=5
    # )
# save_to_pickle(fw, './fw.pkl.gz')
# save_to_pickle(bw, './bw.pkl.gz')
print('Hierachical....')

print(fw_hier[0]['positive_reactions'])
print(len(fw_hier[0]['unrank']))

print(bw_hier[0]['positive_reactions'])
print(len(bw_hier[0]['unrank']))


# print('Normal....')
# print(fw[0]['positive_reactions'])
# print(len(fw[0]['unrank']))

# print(bw[0]['positive_reactions'])
# print(len(bw[0]['unrank']))

