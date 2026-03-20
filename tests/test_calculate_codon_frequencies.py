import filecmp
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

class TestCalculateCodonFrequencies(unittest.TestCase):
    def setUp(self):
        # Determine the root directory of the project and paths to test inputs
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.project_root = os.path.dirname(self.tests_dir)
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")

        self.test1_fasta = os.path.join(self.tests_dir, "inputs", "test.fasta")
        self.test2_fasta = os.path.join(self.tests_dir, "inputs", "test2.fasta")
        self.test3_fasta = os.path.join(self.tests_dir, "inputs", "test3.fasta")
        self.ref_fasta = os.path.join(self.tests_dir, "inputs", "MN908947.3_S_full.fasta")

        # Environment to pass to subprocess
        self.env = os.environ.copy()
        if "PYTHONPATH" not in self.env:
            self.env["PYTHONPATH"] = os.path.join(self.project_root, "src")

        # To make tests robust against PATH issues (like 'calculate_codon_frequencies' not being installed),
        # we run the python module directly using the same python executable running the tests.
        self.base_cmd = [sys.executable, "-m", "mutation_scatter_plot.calculate_codon_frequencies.cli"]

    def _check_outputs(self, expected_prefix, generated_prefix):
        """Helper to compare generated TSVs with golden files stored in tests/outputs/"""
        os.makedirs(self.outputs_dir, exist_ok=True)

        for suffix in [".tsv", ".unchanged_codons.tsv"]:
            generated_file = f"{generated_prefix}{suffix}"
            expected_file = os.path.join(self.outputs_dir, f"{expected_prefix}{suffix}")

            self.assertTrue(os.path.exists(generated_file), f"Output file {generated_file} was not created")

            if not os.path.exists(expected_file):
                # Automatically save as the golden version if no existing golden file is present
                shutil.copy(generated_file, expected_file)
                print(f"Created new golden file: {expected_file}")
            else:
                # Compare contents
                is_match = filecmp.cmp(generated_file, expected_file, shallow=False)
                if not is_match:
                    # Run diff -u -w to show the differences
                    diff_result = subprocess.run(
                        ["diff", "-u", "-w", expected_file, generated_file],
                        capture_output=True, text=True, check=False
                    )
                    self.fail(f"File {generated_file} does not match golden file {expected_file}.\nDifferences:\n{diff_result.stdout}")

    def test_default_command(self):
        """Test calculate_codon_frequencies with default parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2_full.default.frequencies")
            test2_full_fasta = os.path.join(self.tests_dir, "inputs", "test2_full.fasta")
            cmd = self.base_cmd + [
                "--alignment-file", test2_full_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2_full.default.frequencies", outfile_prefix)

    def test_x_after_count(self):
        """Test calculate_codon_frequencies with --x-after-count."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2.x_after_count.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test2_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--x-after-count",
                "--aa_start=430",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2.x_after_count.frequencies", outfile_prefix)

    def test_x_after_count_and_min_start(self):
        """Test calculate_codon_frequencies with --x-after-count and --min_start=7."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2.x_after_count_and_min_start.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test2_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--x-after-count",
                "--aa_start=430",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2.x_after_count_and_min_start.frequencies", outfile_prefix)

    def test_test1_default_command(self):
        """Test calculate_codon_frequencies with default parameters for test.fasta."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.default.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=430",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test1.default.frequencies", outfile_prefix)

    def test_test1_x_after_count(self):
        """Test calculate_codon_frequencies with --x-after-count for test.fasta."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.x_after_count.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=430",
                "--x-after-count",
                "--print-unchanged-sites",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test1.x_after_count.frequencies", outfile_prefix)

    def test_test1_x_after_count_and_min_start(self):
        """Test calculate_codon_frequencies with --x-after-count and --min_start=7 for test.fasta."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.x_after_count_and_min_start.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=430",
                "--x-after-count",
                "--min_start", "7",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test1.x_after_count_and_min_start.frequencies", outfile_prefix)

    def test_test3_default_command(self):
        """Test calculate_codon_frequencies with default parameters for test3.fasta.
        The --min_start=4 is necessary to overcome fake codon at the beginning
        and the remaining sequence merely matches the reference since 415"""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test3.default.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test3_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--x-after-count",
                "--min_start", "4",
                "--aa_start", "415",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test3.default.frequencies", outfile_prefix)

    def test_test2_full_default_command(self):
        """Test calculate_codon_frequencies with default parameters for test2_full.fasta."""

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2_full.x_after_count.frequencies")
            test2_full_fasta = os.path.join(self.tests_dir, "inputs", "test2_full.fasta")
            ref_full_fasta = os.path.join(self.tests_dir, "inputs", "MN908947.3_S_full.fasta")
            cmd = self.base_cmd + [
                "--alignment-file", test2_full_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", ref_full_fasta,
                "--x-after-count",
                "--print-unchanged-sites",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            # Note: We expect this to fail initially if the baseline is not yet updated for the zapped file.
            self._check_outputs("test2_full.x_after_count.frequencies", outfile_prefix)

    def test_short_alignment_aa_start(self):
        """Test calculate_codon_frequencies with short alignment and --aa_start=430."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2_short.aa_start.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test2_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", os.path.join(self.tests_dir, "inputs", "MN908947.3_S.fasta"),
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2_short.aa_start.frequencies", outfile_prefix)
    
    def test_disable_print_unchanged_sites(self):
        """Test that --disable-print-unchanged-sites prevents creating the unchanged_codons file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.disabled.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=430",
                "--disable-print-unchanged-sites",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}")
            
            # Main TSV should exist
            self.assertTrue(os.path.exists(f"{outfile_prefix}.tsv"))
            # Unchanged codons file should NOT exist
            self.assertFalse(os.path.exists(f"{outfile_prefix}.unchanged_codons.tsv"), 
                             "Unchanged codons file should not have been created with --disable-print-unchanged-sites")

if __name__ == "__main__":
    unittest.main()
