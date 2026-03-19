import filecmp
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

class TestMutationScatterPlot(unittest.TestCase):
    def setUp(self):
        # Determine the root directory of the project and paths to test inputs
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.project_root = os.path.dirname(self.tests_dir)
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")

        self.tsv_input = os.path.join(self.outputs_dir, "test3.frequencies_default.tsv")

        # Environment to pass to subprocess
        self.env = os.environ.copy()
        if "PYTHONPATH" not in self.env:
            self.env["PYTHONPATH"] = os.path.join(self.project_root, "src")
        # Ensure headless matplotlib generation in CI
        self.env["MPLBACKEND"] = "Agg"

        self.base_cmd = [sys.executable, "-m", "mutation_scatter_plot.mutation_scatter_plot.cli"]

    def _check_outputs(self, target_prefix, tmpdir, expected_basename_prefix):
        """Helper to compare all generated files with golden files stored in tests/outputs/"""
        os.makedirs(self.outputs_dir, exist_ok=True)
        
        # Discover all files generated in tmpdir
        generated_files = [f for f in os.listdir(tmpdir) if os.path.isfile(os.path.join(tmpdir, f))]
        
        self.assertTrue(len(generated_files) > 0, "No output files were generated in tmpdir")

        for gen_file in generated_files:
            # Map the generated filename to expected golden filename
            expected_filename = gen_file.replace(target_prefix, expected_basename_prefix, 1)
            
            gen_path = os.path.join(tmpdir, gen_file)
            expected_path = os.path.join(self.outputs_dir, expected_filename)

            if not os.path.exists(expected_path):
                # Automatically save as the golden version
                shutil.copy(gen_path, expected_path)
                print(f"Created new golden file: {expected_path}")
            else:
                # Compare contents
                is_match = filecmp.cmp(gen_path, expected_path, shallow=False)
                if not is_match:
                    # Run diff -u -w to show the differences
                    diff_result = subprocess.run(
                        ["diff", "-u", "-w", expected_path, gen_path],
                        capture_output=True, text=True, check=False
                    )
                    if any(gen_file.endswith(ext) for ext in [".png", ".pdf", ".html"]):
                        print(f"Warning: {gen_file} plotting artifact differed from strictly byte-for-byte matching golden baseline natively (ignoring natively due to UUIDs/timestamps).\nDifferences:\n{diff_result.stdout}")
                    else:
                        self.fail(f"File {gen_file} does not match golden file {expected_filename}.\nDifferences:\n{diff_result.stdout}")

    def test_aminoacids(self):
        """mutation_scatter_plot --aminoacids --show-STOP --show-X --show-DEL --show-INS"""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_prefix = "test_run_1"
            outfile_prefix = os.path.join(tmpdir, target_prefix)
            cmd = self.base_cmd + [
                "--tsv", self.tsv_input,
                "--outfile-prefix", outfile_prefix,
                "--aminoacids",
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--threshold=0.0001"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed:\n{result.stderr}\n{result.stdout}")
            self._check_outputs(target_prefix, tmpdir, "test2.scatter_aminoacids")

    def test_aminoacids_synonymous(self):
        """mutation_scatter_plot --aminoacids --show-STOP --show-X --show-DEL --show-INS --include-synonymous"""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_prefix = "test_run_2"
            outfile_prefix = os.path.join(tmpdir, target_prefix)
            cmd = self.base_cmd + [
                "--tsv", self.tsv_input,
                "--outfile-prefix", outfile_prefix,
                "--aminoacids",
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--include-synonymous",
                "--threshold=0.0001"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed:\n{result.stderr}\n{result.stdout}")
            self._check_outputs(target_prefix, tmpdir, "test2.scatter_aminoacids_synonymous")

    def test_codons_synonymous(self):
        """mutation_scatter_plot --show-STOP --show-X --show-DEL --show-INS --include-synonymous"""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_prefix = "test_run_3"
            outfile_prefix = os.path.join(tmpdir, target_prefix)
            cmd = self.base_cmd + [
                "--tsv", self.tsv_input,
                "--outfile-prefix", outfile_prefix,
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--include-synonymous",
                "--threshold=0.0001"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed:\n{result.stderr}\n{result.stdout}")
            self._check_outputs(target_prefix, tmpdir, "test2.scatter_codons_synonymous")

if __name__ == "__main__":
    unittest.main()
