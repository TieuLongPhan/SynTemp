import unittest
import pandas as pd
from pathlib import Path
from synkit.IO.data_io import load_database
from syntemp.pipeline import (
    normalize_rsmi_dict,
    normalize_rsmi_list,
    rebalance,
    clean,
    run_aam,
    extract_its,
    rule_extract,
)

root_dir = Path(__file__).parents[1]


class TestPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Suppose this function loads test data from a specified source
        cls.data = load_database(f"{root_dir}/Data/Testcase/demo.json.gz")[:20]
        cls.data_df = pd.DataFrame(cls.data)

    def test_normalize_rsmi_dict(self):
        # Example test case for normalize_rsmi_dict
        result, issue = normalize_rsmi_dict(self.data[0], "reactions")
        self.assertIn("reactions", result)
        self.assertIsInstance(issue, dict)

    def test_normalize_rsmi_list_no_issue(self):
        # Testing with a DataFrame
        normalized_data, issues = normalize_rsmi_list(self.data_df, "reactions")
        self.assertEqual(len(normalized_data), len(self.data))
        self.assertEqual(len(issues), 0)
        self.assertTrue(all(issue for issue in issues))

    def test_normalize_rsmi_list_issue(self):
        # Testing with a DataFrame
        d = [
            {"R-id": 1, "reactions": "C=C>>CC"},
            {"R-id": 2, "reactions": "C=C>>CCO(C)C"},
        ]
        normalized_data, issues = normalize_rsmi_list(d, "reactions")
        self.assertEqual(len(normalized_data), 1)
        self.assertEqual(len(issues), 1)
        self.assertTrue(all(issue for issue in issues))

    def test_rebalance(self):
        # Testing the rebalance function
        balanced_data = rebalance(self.data_df, "reactions", "id")
        self.assertEqual(len(balanced_data), len(self.data))

    def test_clean(self):
        # Test the clean function
        cleaned_data = clean(self.data_df)
        self.assertEqual(len(cleaned_data), len(self.data))

    def test_run_aam(self):
        # Assuming run_aam needs a list of mappers
        mapper_types = ["rxn_mapper"]
        # mapper_types = ["rxn_mapper", "graphormer", "local_mapper"]
        aam_results = run_aam(self.data, mapper_types)
        self.assertEqual(len(aam_results), len(self.data))
        self.assertTrue(any(keys in aam_results[0] for keys in mapper_types))

    def test_extract_its(self):
        mapper_types = ["rxn_mapper"]
        aam_results = run_aam(self.data, mapper_types)
        its_correct, its_incorrect, all_uncertain_hydrogen = extract_its(
            aam_results, mapper_types=mapper_types
        )
        self.assertEqual(len(its_correct), 17)
        self.assertEqual(len(its_incorrect), 0)
        self.assertEqual(len(all_uncertain_hydrogen), 3)

    def test_rule_extract(self):
        mapper_types = ["rxn_mapper"]
        aam_results = run_aam(self.data, mapper_types)
        its_correct, _, _ = extract_its(aam_results, mapper_types=mapper_types)
        reaction_dicts, templates, hier_templates = rule_extract(
            its_correct, ["element", "charge"]
        )
        self.assertEqual(len(reaction_dicts), 17)
        self.assertEqual(len(templates), 4)
        self.assertEqual(len(hier_templates), 4)


if __name__ == "__main__":
    unittest.main()
