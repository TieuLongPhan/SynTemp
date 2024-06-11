from SynTemp.SynRule.rule_benchmark import RuleBenchmark
from SynTemp.SynUtils.utils import load_database
data = load_database('./Data/DPO/USPTO_50K/test.json.gz')
fw_hier, bw_hier = RuleBenchmark.reproduce_reactions(
        database=data[0:1],
        rule_class_col='R-id',
        rule_file_path='./Data/DPO/USPTO_50K/Hydrogen/R0',
        original_rsmi_col='reactions',
        repeat_times=1,
        use_specific_rules=False,
        verbosity=0,
        job_timeout=5,
        hierarchical=True,
        max_radius=2
    )
print(fw_hier)

print(len(fw_hier[0]['unrank']))

print(bw_hier)
print(len(bw_hier[0]['unrank']))

# fw, bw = RuleBenchmark.reproduce_reactions(
#         database=data[0:1],
#         rule_class_col='R-id',
#         rule_file_path='./Data/DPO/USPTO_50K/Good_hydrogen/R0',
#         original_rsmi_col='reactions',
#         repeat_times=1,
#         use_specific_rules=False,
#         verbosity=0,
#         job_timeout=5
#     )


# print(len(fw[0]['unrank']))