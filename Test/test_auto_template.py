import unittest
from pathlib import Path
from syntemp.auto_template import AutoTemp
from synkit.IO.data_io import load_database


root_dir = Path(__file__).parents[1]


class TestAutoTemp(unittest.TestCase):

    def setUp(self) -> None:
        self.data = load_database(f"{root_dir}/Data/Testcase/demo.json.gz")[:20]
        self.radii = 3
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
            max_radius=self.radii,
        )

    def test_temp_extract(self):
        (_, templates, _, _, _) = self.auto.temp_extract(self.data, lib_path=None)
        self.assertIn("gml", templates[0][0])
        self.assertEqual(len(templates), self.radii + 1)


if __name__ == "__main__":
    unittest.main()
