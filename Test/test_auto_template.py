import unittest
from syntemp.auto_template import AutoTemp
from syntemp.SynUtils.utils import load_database
from pathlib import Path

root_dir = Path(__file__).parents[1]


class TestAutoTemp(unittest.TestCase):

    def setUp(self) -> None:
        self.data = load_database(f"{root_dir}/Data/Testcase/demo.json.gz")[:20]
        self.auto = AutoTemp(
            rebalancing=True,
            mapper_types=["rxn_mapper", "graphormer", "local_mapper"],
            id="R-id",
            rsmi="reactions",
            n_jobs=1,
            verbose=2,
            batch_size=50,
            job_timeout=None,
            safe_mode=False,
            save_dir=f"{root_dir}/Data/Testcase/Test",
            fix_hydrogen=True,
            refinement_its=True,
        )

    def test_temp_extract(self):
        (rules, _, _, _, _, _) = self.auto.temp_extract(self.data, lib_path=None)
        self.assertIn("ruleID", rules[0][0])
        self.assertEqual(len(rules[0]), 11)

    def test_temp_extract_lib(self):
        print(f"{root_dir}/Data/Testcase/Compose/SingleRule")
        (rules, _, _, _, _, _) = self.auto.temp_extract(
            self.data, lib_path=f"{root_dir}/Data/Testcase/Compose/SingleRule"
        )  # 1 rules exist
        self.assertIn("ruleID", rules[0][0])
        self.assertEqual(len(rules[0]), 9)


if __name__ == "__main__":
    unittest.main()
