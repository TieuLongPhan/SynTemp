import unittest
from syntemp.SynAAM.atom_map_consensus import AAMConsensus


class TestRDT(unittest.TestCase):

    def setUp(self):
        self.mapper_types = ["local_mapper", "rxn_mapper", "graphormer", "rdt"]
        self.data = [
            {
                "R-id": 0,
                "record": "CC(C)(C)[Si](C)(C)OC(COCC=1C=CC=CC=1)CCC(CC=C)OS(=O)"
                + "(=O)C.O>>C1=CC(=CC=C1)COCC2OC(CC2)CC=C.CC(C)(C)[Si](C)(C)O.CS(=O)(=O)O",
            },
            {
                "R-id": 1,
                "record": "C[Si](C)(C)CCCCCCCCCCCCCCNC=1C=CC(=CC=1)C=CC(OCC)=O.O>>C[Si](C)"
                + "(C)CCCCCCCCCCCCCCNC1=CC=C(C=CC(O)=O)C=C1.CCO",
            },
            {
                "R-id": 2,
                "record": "CCOC(=O)Cc1csc(n1)-c1ccc(cc1C)C(CC)(CC)c1ccc(C=CC(O)(C(F)(F)F)C(F)"
                + "(F)F)c(C)c1.O>>CCC(CC)(c1ccc(C=CC(O)(C(F)(F)F)C(F)(F)F)c(C)c1)c1ccc"
                + "(-c2nc(CC(O)=O)cs2)c(C)c1.CCO",
            },
        ]

    def test_single_consensus(self):
        """Test that single_consensus includes results from all specified mappers."""
        aam = AAMConsensus(self.data, mappers=self.mapper_types)
        result = aam.single_consensus(self.data[0], rsmi_column="record")

        # Check if each mapper type is in the result dictionary
        for mapper in self.mapper_types:
            with self.subTest(mapper=mapper):
                self.assertIn(
                    mapper,
                    result,
                    f"{mapper} results should be present" + "in the consensus output",
                )

    def test_batch_consensus(self):
        aam = AAMConsensus(self.data, mappers=self.mapper_types)
        results = aam.batch_consensus(
            self.data,
            rsmi_column="record",
            batch_size=len(self.data),
            job_timeout=None,
            safe_mode=False,
        )
        # Check if each mapper type is in the result dictionary for each entry
        for i, entry in enumerate(results):
            for mapper in self.mapper_types:
                with self.subTest(mapper=mapper, entry=i):
                    self.assertIn(
                        mapper,
                        entry,
                        f"{mapper} results should be present"
                        + "in the consensus output for entry {i}",
                    )


if __name__ == "__main__":
    unittest.main()
