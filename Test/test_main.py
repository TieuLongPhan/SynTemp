import unittest
import os
import subprocess
import tempfile
import shutil


class TestCMD(unittest.TestCase):

    def setUp(self):
        """Set up a temporary directory and create a test CSV file."""
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

        # Define the path for the test CSV file
        test_csv_path = os.path.join(self.test_dir, "test.csv")

        # Create a test CSV file
        with open(test_csv_path, "w") as file:
            file.write(
                "R-id,reaction\n0,"
                + "COC(=O)[C@H](CCCCNC(=O)OCc1ccccc1)NC(=O)Nc1cc(OC)cc(C(C)"
                + "(C)C)c1O>>COC(=O)[C@H](CCCCN)NC(=O)Nc1cc(OC)cc(C(C)(C)C)c1O"
            )

        # Define paths for output log and save directory
        self.log_path = os.path.join(self.test_dir, "log.txt")
        self.save_dir = os.path.join(self.test_dir, "output")

    def test_cmd(self):
        """Test the command-line interface of the SynTemp module."""
        # Construct the command
        cmd = [
            "python",
            "-m",
            "SynTemp",
            "--data_path",
            os.path.join(self.test_dir, "test.csv"),
            "--rebalancing",
            "--id",
            "R-id",
            "--rsmi",
            "reaction",
            "--rerun_aam",
            "--fix_hydrogen",
            "--log_file",
            self.log_path,
            "--save_dir",
            self.save_dir,
        ]

        # Execute the command
        subprocess.run(cmd, check=True)

        # Check if output exists in the save directory
        self.assertTrue(
            os.path.exists(self.save_dir), "Output directory does not exist."
        )
        self.assertTrue(
            os.path.exists(f"{self.save_dir}/R0"), "Output directory does not exist."
        )

    def tearDown(self):
        """Clean up after tests."""
        shutil.rmtree(self.test_dir)


if __name__ == "__main__":
    unittest.main()
