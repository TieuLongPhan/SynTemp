import unittest
from SynTemp.auto_template import AutoTemp
from SynTemp.SynUtils.utils import load_database
from pathlib import Path

root_dir = Path(__file__).parents[1]


class TestRuleCluster(unittest.TestCase):

    def setUp(self) -> None:
        self.data = load_database(f"{root_dir}/Data/Testcase/demo.json.gz")[:10]
        self.auto = AutoTemp(
            rebalancing=True,
            mapper_types=["rxn_mapper", "graphormer"],
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

    def test_fit(self):
        (rules, _, _, _, _, _) = self.auto.fit(self.data)
        self.assertIn("ruleID", rules[0][0])


if __name__ == "__main__":
    unittest.main()
