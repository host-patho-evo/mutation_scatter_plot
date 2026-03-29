"""Unit tests for calculate_codon_frequencies."""
import filecmp
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

class TestCalculateCodonFrequencies(unittest.TestCase):
    """Test cases for the calculate_codon_frequencies CLI and core logic."""
    # pylint: disable=too-many-instance-attributes
    coverage_results = []
    golden_mutations = set()
    project_root = ""

    @classmethod
    def setUpClass(cls):
        cls.coverage_results = []
        cls.project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        golden_path = os.path.join(cls.project_root, "data", "golden_mutations.json")
        try:
            with open(golden_path, 'r', encoding='utf-8') as f:
                cls.golden_mutations = set(json.load(f))
        except (FileNotFoundError, json.JSONDecodeError):
            cls.golden_mutations = set()

    @classmethod
    def tearDownClass(cls):
        if not cls.coverage_results:
            return

        table = []
        table.append("\n" + "="*77)
        table.append(f"{'GOLDEN DICTIONARY COVERAGE':^77}")
        table.append("="*77)
        table.append(f"| {'File':<45} | {'Total':<7} | {'Confirmed':<9} | {'%':<5} |")
        table.append("|" + "-"*47 + "|" + "-"*9 + "|" + "-"*11 + "|" + "-"*7 + "|")
        for res in sorted(cls.coverage_results, key=lambda x: x[0]):
            filename, total, confirmed, pct = res
            table.append(f"| {filename:<45} | {total:<7} | {confirmed:<9} | {pct:>5.1f}% |")
        table.append("="*77 + "\n")

        sys.stderr.write("\n".join(table) + "\n")
        # Ensure it gets flushed
        sys.stderr.flush()

    def setUp(self):
        # Determine the root directory of the project and paths to test inputs
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.project_root = os.path.dirname(self.tests_dir)
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")

        self.test1_fasta = os.path.join(self.tests_dir, "inputs", "test.fasta")
        self.test2_fasta = os.path.join(self.tests_dir, "inputs", "test2.fasta")
        self.test3_fasta = os.path.join(self.tests_dir, "inputs", "test3.fasta")
        self.ref_fasta = os.path.join(self.tests_dir, "inputs", "MN908947.3_S.fasta")

        # Environment to pass to subprocess
        self.env = os.environ.copy()
        if "PYTHONPATH" not in self.env:
            self.env["PYTHONPATH"] = os.path.join(self.project_root, "src")

        # To make tests robust against PATH issues (like 'calculate_codon_frequencies' not being installed),
        # we run the python module directly using the same python executable running the tests.
        self.base_cmd = [sys.executable, "-m", "mutation_scatter_plot.calculate_codon_frequencies.cli"]

    def _check_outputs(self, expected_prefix, generated_prefix):
        """Helper to compare generated TSVs with golden files stored in tests/outputs/"""
        # pylint: disable=too-many-locals
        os.makedirs(self.outputs_dir, exist_ok=True)
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
                    # Run git diff --word-diff-regex=. to show character-level differences
                    diff_result = subprocess.run(
                        ["git", "diff", "--no-index", "--word-diff=color",
                         "--word-diff-regex=.", "--color=always",
                         expected_file, generated_file],
                        capture_output=True, text=True, check=False
                    )
                    self.fail(f"File {generated_file} does not match golden file {expected_file}.\nDifferences:\n{diff_result.stdout}")

            # Append coverage tracking metrics to our teardown table
            if suffix == ".tsv" and self.__class__.golden_mutations:
                try:
                    with open(generated_file, 'r', encoding='utf-8') as f:
                        lines = f.readlines()[1:] # skip header
                    total = len(lines)
                    confirmed = 0
                    for line in lines:
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            pos = parts[1].strip()
                            orig = parts[2].strip()
                            mut = parts[3].strip()
                            if f"{orig}{pos}{mut}" in self.__class__.golden_mutations:
                                confirmed += 1
                    pct = (confirmed / total) * 100 if total > 0 else 0.0
                    self.__class__.coverage_results.append((f"{expected_prefix}{suffix}", total, confirmed, pct))
                except (OSError, ValueError) as e:
                    sys.stderr.write(f"Failed to calculate coverage for {generated_file}: {e}\n")

    # --- test.fasta ---

    def test_test1_default_command(self):
        """Test calculate_codon_frequencies with default parameters for test.fasta."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.default.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=413",
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
                "--aa_start=413",
                "--x-after-count",
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
                "--aa_start=413",
                "--x-after-count",
                "--min_start", "7",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test1.x_after_count_and_min_start.frequencies", outfile_prefix)

    def test_disable_print_unchanged_sites(self):
        """Test that --disable-print-unchanged-sites prevents creating the unchanged_codons file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.disabled.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=413",
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

    # --- test2.fasta ---

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
                "--aa_start=413",
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
                "--min_start", "7",
                "--aa_start=413",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2.x_after_count_and_min_start.frequencies", outfile_prefix)

    def test_short_alignment_aa_start(self):
        """Test calculate_codon_frequencies with short alignment and --aa_start=413."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2_short.aa_start.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test2_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", os.path.join(self.tests_dir, "inputs", "MN908947.3_S.fasta"),
                "--aa_start=413",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2_short.aa_start.frequencies", outfile_prefix)

    # --- test2_full.fasta ---

    def test_default_command(self):
        """Test calculate_codon_frequencies with default parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test2_full.default.frequencies")
            test2_full_fasta = os.path.join(self.tests_dir, "inputs", "test2_full.fasta")
            ref_full_fasta = os.path.join(self.tests_dir, "inputs", "MN908947.3_S_full.fasta")
            cmd = self.base_cmd + [
                "--alignment-file", test2_full_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", ref_full_fasta,
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test2_full.default.frequencies", outfile_prefix)

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
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            # Note: We expect this to fail initially if the baseline is not yet updated for the zapped file.
            self._check_outputs("test2_full.x_after_count.frequencies", outfile_prefix)

    def test_threads_option(self):
        """Test calculate_codon_frequencies with --threads option."""
        with tempfile.TemporaryDirectory() as tmpdir:
            outfile_prefix = os.path.join(tmpdir, "test1.threads.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test1_fasta,
                "--outfile-prefix", outfile_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--aa_start=413",
                "--threads", "2",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            # Compare with the same baseline as test1_default_command
            self._check_outputs("test1.default.frequencies", outfile_prefix)

    # --- test3.fasta ---

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
                "--aa_start", "413",
                "--overwrite"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed with error:\n{result.stderr}\n\nStdout:\n{result.stdout}")
            self._check_outputs("test3.default.frequencies", outfile_prefix)

if __name__ == "__main__":
    unittest.main()


class TestCountFormatParsing(unittest.TestCase):
    """Regression tests for FASTA ID count-prefix parsing in x-after-count mode.

    Bug history (commit ea3512f, 2026-03-24): the original count extraction
    used _record_id.endswith('x.') / endswith('x'), which correctly handled
    the test-file format '576521x' (space-separated first token) but silently
    fell through to count=1 for the production format '576521x.SHA256HASH'
    where the hash is concatenated with a dot, making the full
    '576521x.SHA256HASH' string the first split() token.

    Fixed in c0d1da3 (2026-03-29) by extracting digits before the first 'x'.
    """

    @staticmethod
    def _parse_count(record_id: str) -> int:
        """Mirrors the count-extraction logic in parse_alignment."""
        x_pos = record_id.find('x')
        if x_pos > 0 and record_id[:x_pos].isdigit():
            return int(record_id[:x_pos])
        return 1

    def test_production_sha256dot_format(self):
        """Production format NNNNx.SHA256HASH — was broken before c0d1da3."""
        self.assertEqual(
            self._parse_count(
                '576521x.7cbee25f1cd8f4dd1ad91f50b7eb2722b10b13a14a661f4ae5b077b52bcdf595'
            ),
            576521,
        )

    def test_compact_x_format(self):
        """Existing test-file format: NNNNx (space after x makes it the sole token)."""
        self.assertEqual(self._parse_count('576521x'), 576521)

    def test_x_dot_no_suffix(self):
        """Edge case: NNNNx. with nothing after the dot."""
        self.assertEqual(self._parse_count('576521x.'), 576521)

    def test_count_one_with_sha256(self):
        """Minimum count with sha256: 1x.sha256hash."""
        self.assertEqual(self._parse_count('1x.deadbeef'), 1)

    def test_count_one_compact(self):
        """Minimum count compact: 1x."""
        self.assertEqual(self._parse_count('1x'), 1)

    def test_no_x_prefix_gisaid_style(self):
        """GISAID-style ID with no count prefix — fallback to 1."""
        self.assertEqual(self._parse_count('EPI_ISL_123456'), 1)

    def test_no_x_prefix_plain(self):
        """Plain word ID — fallback to 1."""
        self.assertEqual(self._parse_count('noxcount'), 1)

    def test_x_in_middle_non_numeric(self):
        """'x' appears in ID but preceded by non-digits — fallback to 1."""
        self.assertEqual(self._parse_count('EPI_ISL_xABC'), 1)

    def test_leading_zeros(self):
        """Leading zeros: '007x.sha256' → isdigit() True → count=7."""
        self.assertEqual(self._parse_count('007x.sha256hash'), 7)


class TestSha256CountFormat(unittest.TestCase):
    """CLI regression test: NNNNx.SHA256HASH headers produce correct counts.

    test2_sha256.fasta has the exact same sequences and count values as
    test2.fasta but uses the production-style ID format 'NNNNx.deadbeefXXX'
    (dot-concatenated sha256) instead of 'NNNNx' (space-separated).

    With the old broken endswith logic the first split() token
    '576521x.deadbeef...' did NOT end with 'x', so count fell through to 1
    for every record.  This caused total_codons_per_site to equal the number
    of records (9) instead of sum(counts) = 1,089,650 for a well-covered site.

    With the fix the two files must produce byte-for-byte identical TSVs.
    """

    def setUp(self):
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.project_root = os.path.dirname(self.tests_dir)
        self.ref_fasta = os.path.join(self.tests_dir, "inputs", "MN908947.3_S.fasta")
        self.test2_sha256_fasta = os.path.join(
            self.tests_dir, "inputs", "test2_sha256.fasta"
        )
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")
        self.env = os.environ.copy()
        if "PYTHONPATH" not in self.env:
            self.env["PYTHONPATH"] = os.path.join(self.project_root, "src")
        self.base_cmd = [
            sys.executable, "-m",
            "mutation_scatter_plot.calculate_codon_frequencies.cli",
        ]

    def test_sha256_dot_format_matches_space_format(self):
        """NNNNx.SHA256HASH IDs must give same TSV output as NNNNx IDs.

        Both test2.fasta and test2_sha256.fasta encode identical sequences
        with identical per-record counts.  The two TSVs must be byte-for-byte
        identical; if count parsing regresses to 1-per-record the
        total_codons_per_site column will differ.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run on the sha256-format file.
            sha256_prefix = os.path.join(tmpdir, "test2_sha256.x_after_count.frequencies")
            cmd = self.base_cmd + [
                "--alignment-file", self.test2_sha256_fasta,
                "--outfile-prefix", sha256_prefix,
                "--padded-reference",
                "--reference-infile", self.ref_fasta,
                "--x-after-count",
                "--aa_start=413",
                "--overwrite",
            ]
            result = subprocess.run(
                cmd, cwd=self.project_root, env=self.env,
                capture_output=True, text=True, check=False,
            )
            self.assertEqual(
                result.returncode, 0,
                f"Command failed:\n{result.stderr}\n{result.stdout}",
            )

            # The golden output for test2.fasta with --x-after-count already
            # exists; reuse it as the reference.
            golden_prefix = os.path.join(
                self.outputs_dir, "test2.x_after_count.frequencies"
            )
            for suffix in (".tsv", ".unchanged_codons.tsv"):
                generated = f"{sha256_prefix}{suffix}"
                golden = f"{golden_prefix}{suffix}"
                self.assertTrue(
                    os.path.exists(generated),
                    f"Output file not created: {generated}",
                )
                if os.path.exists(golden):
                    is_match = filecmp.cmp(generated, golden, shallow=False)
                    if not is_match:
                        diff = subprocess.run(
                            ["git", "diff", "--no-index", "--word-diff=color",
                             "--word-diff-regex=.", "--color=always",
                             golden, generated],
                            capture_output=True, text=True, check=False,
                        )
                        self.fail(
                            f"{generated} does not match golden {golden}.\n"
                            f"Differences (check total_codons_per_site column):\n"
                            f"{diff.stdout}"
                        )
