from pathlib import Path
import sys
root_dir = Path(__file__).parents[1]
sys.path.append(root_dir)
print(root_dir)

def main():
    
    data = load_database('./Data/DPO/USPTO_balance/test.json.gz')

    # fw_hier, bw_hier = RuleBenchmark.reproduce_reactions(
    #         database=data[0:200],
    #         rule_class_col='R-id',
    #         rule_file_path='./Data/DPO/USPTO_balance/Expand/R0',
    #         original_rsmi_col='reactions',
    #         repeat_times=1,
    #         use_specific_rules=False,
    #         verbosity=3,
    #         hierarchical=True,
    #         max_radius=3,
    #         max_solutions = 1000,
    #         templates_threshold = 1
    #     )
    # save_to_pickle(fw_hier, './fw_hier.pkl.gz')
    # save_to_pickle(bw_hier, './bw_hier.pkl.gz')

    fw, bw = RuleBenchmark.reproduce_reactions(
            database=data[0:1],
            rule_class_col='R-id',
            rule_file_path='./Data/DPO/USPTO_balance/Complete/R3',
            original_rsmi_col='reactions',
            repeat_times=1,
            use_specific_rules=False,
            verbosity=0,
        )

    print(fw[0])
    # save_to_pickle(fw, './fw.pkl.gz')
    # save_to_pickle(bw, './bw.pkl.gz')
    # print('Hierachical....')
    # print(fw_hier[0]['unrank_raw'])
    # print(fw_hier[0]['positive_reactions'])
    # print(len(fw_hier[0]['unrank']))

    # print(bw_hier[0]['positive_reactions'])
    # print(len(bw_hier[0]['unrank']))


    print('Normal....')
    print(fw[0]['positive_reactions'])
    print(len(fw[0]['unrank']))

    print(bw[0]['positive_reactions'])
    print(len(bw[0]['unrank']))


if __name__ == "__main__":
    sys.path.append(root_dir)
    from SynTemp.SynRule.rule_benchmark import RuleBenchmark
    from SynTemp.SynUtils.utils import load_database, save_to_pickle
    main()