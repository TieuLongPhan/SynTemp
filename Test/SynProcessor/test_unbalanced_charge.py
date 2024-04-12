import unittest
import sys
from pathlib import Path

root_dir = Path(__file__).parents[2]
sys.path.append(str(root_dir))
from SynTemp.SynProcessor.unbalanced_charge import UnbalancedCharge


class TestUnbalancedCharge(unittest.TestCase):
    def test_fix_negative_charge(self):
        reaction_dict = {
            "reactants": "H2O",
            "products": "OH-",
            "total_charge_in_products": -1,
            "R-id": "R1",
            "label": "test_label",
        }
        expected = {
            "R-id": "R1",
            "new_reaction": "H2O.[Na+]>>OH-.[Na+]",
            "label": "test_label",
            "reactants": "H2O.[Na+]",
            "products": "OH-.[Na+]",
            "total_charge_in_products": 0,
        }
        self.assertEqual(UnbalancedCharge.fix_negative_charge(reaction_dict), expected)

    def test_fix_positive_charge(self):
        reaction_dict = {
            "reactants": "HCl",
            "products": "H+",
            "total_charge_in_products": 1,
            "R-id": "R2",
            "label": "test_label",
        }
        expected = {
            "R-id": "R2",
            "new_reaction": "HCl.[Cl-]>>H+.[Cl-]",
            "label": "test_label",
            "reactants": "HCl.[Cl-]",
            "products": "H+.[Cl-]",
            "total_charge_in_products": 0,
        }
        self.assertEqual(UnbalancedCharge.fix_positive_charge(reaction_dict), expected)

    def test_fix_unbalanced_charged_no_change(self):
        reaction_dict = {
            "reactants": "NaCl",
            "products": "Na+ + Cl-",
            "total_charge_in_products": 0,
            "R-id": "R3",
            "label": "test_label",
        }
        # Since the charge is already balanced, the original dictionary should be returned unchanged.
        self.assertEqual(
            UnbalancedCharge.fix_unbalanced_charged(reaction_dict), reaction_dict
        )

    def test_parallel_fix_unbalanced_charge(self):
        reaction_dicts = [
            {
                "reactants": "H2O",
                "products": "OH-",
                "total_charge_in_products": -1,
                "R-id": "R1",
                "label": "test_label",
            },
            {
                "reactants": "HCl",
                "products": "H+",
                "total_charge_in_products": 1,
                "R-id": "R2",
                "label": "test_label",
            },
        ]
        fixed_reactions = UnbalancedCharge.parallel_fix_unbalanced_charge(
            reaction_dicts
        )
        # Check the length of the returned list to ensure all reactions were processed
        self.assertEqual(len(fixed_reactions), len(reaction_dicts))
        # Further checks can be added to verify the contents of fixed_reactions


if __name__ == "__main__":
    unittest.main()
