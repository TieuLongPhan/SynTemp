from SynTemp.SynRule.rule_benchmark import RuleBenchmark
from SynTemp.SynUtils.utils import load_database, save_to_pickle
data = load_database('./Data/DPO/USPTO_50K/test.json.gz')


# start_index = 3
# end_index = 4
# ban_list = ['CC(C)(C)OC(=O)N(CCCCC1(c2ccc(N)cc2)SCCS1)C[C@H](O[Si](C)(C)C(C)(C)C)c1ccc(O)c2[nH]c(=O)ccc12.COc1cccc(Nc2c(C(N)=O)cnc3c(C)cc(S(=O)(=O)c4cccc(C(=O)O)c4)cc23)c1>>COc1cccc(Nc2c(C(N)=O)cnc3c(C)cc(S(=O)(=O)c4cccc(C(=O)Nc5ccc(C6(CCCCN(C[C@H](O[Si](C)(C)C(C)(C)C)c7ccc(O)c8[nH]c(=O)ccc78)C(=O)OC(C)(C)C)SCCS6)cc5)c4)cc23)c1.O',
# 'CC(C)[Si](O[C@@H]1CC[C@@H](N2CC[C@@H](Cc3c(Cl)cc(-c4ccc(C(=O)N5CCC(C(F)(F)F)CC5)cc4)cc3Cl)C2=O)CC1)(C(C)C)C(C)C.O>>O=C(c1ccc(-c2cc(Cl)c(C[C@@H]3CCN([C@@H]4CC[C@@H](O)CC4)C3=O)c(Cl)c2)cc1)N1CCC(C(F)(F)F)CC1.CC(C)[Si](O)(C(C)C)C(C)C',
# 'CCOC(=O)CC(O)[C@H]1[C@H](N(Cc2ccc(OC)cc2)C(=O)OC(C)(C)C)[C@@H](OCc2ccccc2)[C@@H](OCc2ccccc2)[C@@H]1OCc1ccccc1.O>>CCOC(=O)CC(O)[C@H]1[C@H](NC(=O)OC(C)(C)C)[C@@H](OCc2ccccc2)[C@@H](OCc2ccccc2)[C@@H]1OCc1ccccc1.COc1ccc(CO)cc1',
# 'COC(=O)c1ccc(CN(C(=O)[C@@H]2Cc3ccccc3CN2C(=O)[C@@H](NC(=O)[C@H](C)N(C)C(=O)OC(C)(C)C)C(C)(C)C)[C@H](C)c2ccc(F)cc2)cc1.O>>C[C@@H](C(=O)N[C@H](C(=O)N1Cc2ccccc2C[C@H]1C(=O)N(Cc1ccc(C(=O)O)cc1)[C@H](C)c1ccc(F)cc1)C(C)(C)C)N(C)C(=O)OC(C)(C)C.CO',
# 'COC(=O)[C@H](CCCCNC(=O)[C@H](CCCCNC(=O)[C@H](CCCCNC(=NC(=O)OC(C)(C)C)NC(=O)OC(C)(C)C)NC(=O)OC(C)(C)C)NC(=O)OC(C)(C)C)NC(=O)OC(C)(C)C.O>>CC(C)(C)OC(=O)N=C(NCCCC[C@H](NC(=O)OC(C)(C)C)C(=O)NCCCC[C@H](NC(=O)OC(C)(C)C)C(=O)NCCCC[C@H](NC(=O)OC(C)(C)C)C(=O)O)NC(=O)OC(C)(C)C.CO',
# 'CC(C)(C)OC(=O)N1CCC(COC2CC=C(B3OC(C)(C)C(C)(C)O3)CC2)CC1.Cc1nc(S(C)(=O)=O)ccc1OS(=O)(=O)C(F)(F)F>>Cc1nc(S(C)(=O)=O)ccc1C1=CCC(OCC2CCN(C(=O)OC(C)(C)C)CC2)CC1.CC1(C)OB(OS(=O)(=O)C(F)(F)F)OC1(C)C',
# 'NCC(=O)OCc1ccccc1.CC(=O)N1[C@H](CC23CC4CC(CC(C4)C2)C3)C(=O)N(Cc2ccccc2)c2ccccc2C(=O)C[C@H]1C(=O)O[Na]>>CC(=O)N1[C@H](CC23CC4CC(CC(C4)C2)C3)C(=O)N(Cc2ccccc2)c2ccccc2C(=O)C[C@H]1C(=O)NCC(=O)OCc1ccccc1.O[Na]',
# 'CCCC[Sn](C#Cc1ccccc1)(CCCC)CCCC.COC(=O)c1cc(NC(C)=O)cc(NC(C)=O)c1Br>>COC(=O)c1cc(NC(C)=O)cc(NC(C)=O)c1C#Cc1ccccc1.CCCC[Sn](Br)(CCCC)CCCC',
# 'CC1(C)CC(=O)N(COCC[Si](C)(C)C)c2c(CCl)cccc21.COc1ccccc1COCCCOc1ccc(C2CCN(C(=O)OC(C)(C)C)CC2O)cc1>>COc1ccccc1COCCCOc1ccc(C2CCN(C(=O)OC(C)(C)C)CC2OCc2cccc3c2N(COCC[Si](C)(C)C)C(=O)CC3(C)C)cc1.Cl',
# 'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C.CCOC(=O)CNCCCn1c(-c2ccc(F)cc2)csc1=Nc1ccc(Cl)cc1OC>>CCOC(=O)CN(CCCn1c(-c2ccc(F)cc2)csc1=Nc1ccc(Cl)cc1OC)C(=O)OC(C)(C)C.CC(C)(C)OC(=O)O',
# 'COCCCNc1nc(C(C)(C)C)ncc1C(=O)N(CC(C)C)[C@H]1C[C@@H](C(=O)N2CCN(C)CC2)CN(C(=O)OC(C)(C)C)C1.O>>COCCCNc1nc(C(C)(C)C)ncc1C(=O)N(CC(C)C)[C@@H]1CNC[C@H](C(=O)N2CCN(C)CC2)C1.CC(C)(C)OC(=O)O',
# 'COCOC[C@@H]1C(C(=O)N(Cc2cccc(Cl)c2Cl)C2CC2)=C(OS(=O)(=O)C(F)(F)F)CCN1C(=O)OC(C)(C)C.OB(O)c1ccc(O)cc1>>COCOC[C@@H]1C(C(=O)N(Cc2cccc(Cl)c2Cl)C2CC2)=C(c2ccc(O)cc2)CCN1C(=O)OC(C)(C)C.O=S(=O)(OB(O)O)C(F)(F)F',
# 'Cc1cc(C[C@@H](OC(=O)N2CCC(N3CCc4ccccc4NC3=O)CC2)C(=O)N2CCC(C3CCN(CC(=O)O)CC3)CC2)cc(C)c1O.O=C1CCCN1CCO>>Cc1cc(C[C@@H](OC(=O)N2CCC(N3CCc4ccccc4NC3=O)CC2)C(=O)N2CCC(C3CCN(CC(=O)OCCN4CCCC4=O)CC3)CC2)cc(C)c1O.O', 'CN(C(=O)OC(C)(C)C)C(C(=O)O)c1cccc(Br)c1.COc1ccccc1-c1nn(COCC[Si](C)(C)C)c2ncc(B3OC(C)(C)C(C)(C)O3)cc12>>COc1ccccc1-c1nn(COCC[Si](C)(C)C)c2ncc(-c3cccc(C(C(=O)O)N(C)C(=O)OC(C)(C)C)c3)cc12.CC1(C)OB(Br)OC1(C)C',
# 'COc1ccc2c(O[C@@H]3C[C@H]4C(=O)N[C@]5(C(=O)O)C[C@H]5C=CCCCCC[C@H](NC(=O)OC(C)(C)C)C(=O)N4C3)cc(-c3ccccc3)nc2c1.[H].[H]>>COc1ccc2c(O[C@@H]3C[C@H]4C(=O)N[C@]5(C(=O)O)C[C@H]5CCCCCCC[C@H](NC(=O)OC(C)(C)C)C(=O)N4C3)cc(-c3ccccc3)nc2c1',
# 'CC(C)(C)OC(=O)N[C@@H](Cc1ccc(B2OC(C)(C)C(C)(C)O2)cc1F)C(=O)N1CCC[C@H]1C#N.CC(C)c1ccc(CN2CC(Oc3ccc(Br)cn3)C2)cc1>>CC(C)c1ccc(CN2CC(Oc3ccc(-c4ccc(C[C@H](NC(=O)OC(C)(C)C)C(=O)N5CCC[C@H]5C#N)c(F)c4)cn3)C2)cc1.CC1(C)OB(Br)OC1(C)C',
# 'CCOC(=O)/C=C/c1ccc2c(=O)n(CC(C)C)c(CNC(=O)OC(C)(C)C)c(-c3ccccc3)c2c1.[H].[H]>>CCOC(=O)CCc1ccc2c(=O)n(CC(C)C)c(CNC(=O)OC(C)(C)C)c(-c3ccccc3)c2c1',
# 'COc1cc(C[C@@]2(CCO[Si](C)(C)C(C)(C)C)OC(C)(C)O[C@H]2C)c(OC)c2c(OC)c3ccccc3c(OC)c12.O>>COc1cc(C[C@@]2(CCO)OC(C)(C)O[C@H]2C)c(OC)c2c(OC)c3ccccc3c(OC)c12.CC(C)(C)[Si](C)(C)O'
            # ]
# data = [value for value in data if value['reactions'] not in ban_list]
# fw_hier, bw_hier = RuleBenchmark.reproduce_reactions(
#         database=data[0:],
#         rule_class_col='R-id',
#         rule_file_path='./Data/DPO/USPTO_50K/Good_hydrogen/R0',
#         original_rsmi_col='reactions',
#         repeat_times=1,
#         use_specific_rules=False,
#         verbosity=0,
#         job_timeout=5,
#         hierarchical=True,
#         max_radius=3,
#         max_solutions = 10
#     )
# save_to_pickle(fw_hier, './fw_hier.pkl.gz')
# save_to_pickle(bw_hier, './bw_hier.pkl.gz')

fw, bw = RuleBenchmark.reproduce_reactions(
        database=data[0:1],
        rule_class_col='R-id',
        rule_file_path='./Data/DPO/USPTO_50K/Non_hydrogen/R1',
        original_rsmi_col='reactions',
        repeat_times=1,
        use_specific_rules=False,
        verbosity=0,
    )

print(fw)
# save_to_pickle(fw, './fw.pkl.gz')
# save_to_pickle(bw, './bw.pkl.gz')
# print('Hierachical....')

# print(fw_hier[0]['positive_reactions'])
# print(len(fw_hier[0]['unrank']))

# print(bw_hier[0]['positive_reactions'])
# print(len(bw_hier[0]['unrank']))


print('Normal....')
print(fw[0]['positive_reactions'])
print(len(fw[0]['unrank']))

print(bw[0]['positive_reactions'])
print(len(bw[0]['unrank']))

