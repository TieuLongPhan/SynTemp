# from SynTemp.auto_template import AutoTemp
# from SynTemp.SynUtils.utils import load_from_pickle, load_database
# import logging

# # data = load_from_pickle('Data/Testcase/its_graph.pkl.gz')
# data = load_database("Data/Testcase/demo.json.gz")[:100]
# auto = AutoTemp(
#     rebalancing=True,
#     mapper_types=["rxn_mapper", "graphormer"],
#     id="R-id",
#     rsmi="reactions",
#     n_jobs=1,
#     verbose=2,
#     batch_size=50,
#     job_timeout=None,
#     safe_mode=False,
#     save_dir="Data/Testcase/Test",
#     fix_hydrogen=True,
#     refinement_its=True,
# )
# (
#     rules,
#     reaction_dicts,
#     templates,
#     hier_templates,
#     its_incorrect,
#     uncertain_hydrogen,
# ) = auto.fit(data)

# print(len(rules))

# print(rules[0][0])
