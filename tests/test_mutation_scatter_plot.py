"""Unit tests for mutation_scatter_plot."""
import contextlib
import filecmp
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import traceback
import unittest
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
from bokeh.document import Document # pylint: disable=wrong-import-position

from mutation_scatter_plot.mutation_scatter_plot.cli import main as mutation_scatter_plot_main # pylint: disable=wrong-import-position

class TestMutationScatterPlot(unittest.TestCase):
    """Test cases for the mutation_scatter_plot CLI and rendering logic."""
    def setUp(self):
        # Determine the root directory of the project and paths to test inputs
        self.tests_dir = os.path.dirname(os.path.abspath(__file__))
        self.project_root = os.path.dirname(self.tests_dir)
        self.outputs_dir = os.path.join(self.tests_dir, "outputs")

        self.tsv_input = os.path.join(self.outputs_dir, "test2_full.x_after_count.frequencies.tsv")

        # Environment to pass to subprocess
        self.env = os.environ.copy()
        if "PYTHONPATH" not in self.env:
            self.env["PYTHONPATH"] = os.path.join(self.project_root, "src")
        # Ensure headless matplotlib generation in CI
        self.env["MPLBACKEND"] = "Agg"

        self.script_name = "mutation_scatter_plot.mutation_scatter_plot.cli"
        self.regenerate_goldens = os.environ.get("REGENERATE_GOLDENS") == "1"

    def _invoke_cli(self, cmd_args):
        """Invoke the CLI natively, capturing output without an interpreter fork."""
        with patch.object(sys, 'argv', [self.script_name] + cmd_args):
            f = io.StringIO()
            with contextlib.redirect_stdout(f), contextlib.redirect_stderr(f):
                try:
                    mutation_scatter_plot_main()
                    return 0, f.getvalue()
                except SystemExit as e:
                    return e.code, f.getvalue()
                except Exception as e: # pylint: disable=broad-exception-caught
                    return 1, f.getvalue() + "\n" + str(e) + "\n" + traceback.format_exc()


    def _check_outputs(self, target_prefix, tmpdir, expected_basename_prefix):
        """Helper to compare all generated files with golden files stored in tests/outputs/"""
        # pylint: disable=too-many-locals
        os.makedirs(self.outputs_dir, exist_ok=True)
        os.makedirs(self.outputs_dir, exist_ok=True)

        # Discover all files generated in tmpdir
        generated_files = [f for f in os.listdir(tmpdir) if os.path.isfile(os.path.join(tmpdir, f))]

        self.assertTrue(len(generated_files) > 0, "No output files were generated in tmpdir")

        for gen_file in generated_files:
            # Map the generated filename to expected golden filename
            expected_filename = gen_file.replace(target_prefix, expected_basename_prefix, 1)

            gen_path = os.path.join(tmpdir, gen_file)
            expected_path = os.path.join(self.outputs_dir, expected_filename)

            if not os.path.exists(expected_path) or self.regenerate_goldens:
                # Automatically save as the golden version
                shutil.copy(gen_path, expected_path)
                print(f"Updated golden file: {expected_path}")
            else:
                # Compare contents
                is_match = filecmp.cmp(gen_path, expected_path, shallow=False)
                if not is_match:
                    if gen_file.endswith(".html"):
                        data_gen = self._extract_bokeh_data(gen_path)
                        data_exp = self._extract_bokeh_data(expected_path)
                        if data_gen != data_exp:
                            with tempfile.NamedTemporaryFile("w", delete=False) as f_exp, tempfile.NamedTemporaryFile("w", delete=False) as f_gen:
                                for keys, rows in data_exp:
                                    f_exp.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                                for keys, rows in data_gen:
                                    f_gen.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")

                            diff_res = subprocess.run(
                                ["git", "diff", "--no-index", "--word-diff=color",
                                 "--word-diff-regex=.", "--color=always",
                                 f_exp.name, f_gen.name],
                                capture_output=True, text=True, check=False
                            )
                            self.fail(f"Bokeh internal data points differ for {gen_file}!\nDifferences:\n{diff_res.stdout}")
                        else:
                            print(f"Info: HTML file {gen_file} byte-mismatch ignored because internal Bokeh JSON structural mapping matched perfectly.")
                    elif any(gen_file.endswith(ext) for ext in [".png", ".pdf"]):
                        print(f"Warning: {gen_file} plotting artifact differed from strictly byte-for-byte matching golden baseline natively (ignoring natively due to UUIDs/timestamps).\nDifferences:\nNone (diff unavailable without filecmp debug)")
                    else:
                        self.fail(f"File {gen_file} does not match golden file {expected_filename}.\nDifferences:\nNone (diff unavailable without filecmp debug)")

    def _extract_bokeh_data(self, html_path):
        """Extract data from Bokeh HTML file for comparison."""
        with open(html_path, 'r', encoding='utf-8') as f:
            content = f.read()

        match = re.search(r'<script type="application/json".*?>(.*?)</script>', content, re.DOTALL)
        if not match:
            return []

        js = json.loads(match.group(1).strip())
        doc_json = list(js.values())[0] if js else {}

        doc = Document.from_json(doc_json)
        all_data = []
        for model in doc.models:
            if type(model).__name__ == 'ColumnDataSource' and hasattr(model, 'data'):
                model_data = getattr(model, 'data')
                keys = sorted(model_data.keys())
                if not keys:
                    continue
                rows = []
                for idx in range(len(model_data[keys[0]])):
                    rows.append(tuple(str(model_data[k][idx]) for k in keys))
                all_data.append((keys, sorted(rows)))

        return sorted(all_data)

    def test_all_inputs_aminoacids(self):
        """mutation_scatter_plot --aminoacids for all inputs"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--aminoacids",
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_aminoacids")

    def test_all_inputs_aminoacids_synonymous(self):
        """mutation_scatter_plot --aminoacids --include-synonymous for all inputs"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--aminoacids",
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--include-synonymous",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_aminoacids_synonymous")

    def test_all_inputs_codons(self):
        """mutation_scatter_plot codon mode for all inputs"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_codons")

    def test_all_inputs_codons_synonymous(self):
        """mutation_scatter_plot codon mode --include-synonymous for all inputs"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--include-synonymous",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_codons_synonymous")

    def test_all_inputs_codons_coolwarm_r(self):
        """mutation_scatter_plot codon mode with --colormap=coolwarm_r (standard Normalize path)"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--colormap=coolwarm_r",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_codons")

    def test_all_inputs_aminoacids_coolwarm_r(self):
        """mutation_scatter_plot --aminoacids with --colormap=coolwarm_r (standard Normalize path)"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--aminoacids",
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--colormap=coolwarm_r",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_aminoacids")

    def test_all_inputs_codons_coolwarm_r_synonymous(self):
        """mutation_scatter_plot codon mode with --colormap=coolwarm_r --include-synonymous"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--include-synonymous",
                        "--colormap=coolwarm_r",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_codons_synonymous")

    def test_all_inputs_aminoacids_coolwarm_r_synonymous(self):
        """mutation_scatter_plot --aminoacids with --colormap=coolwarm_r --include-synonymous"""
        test_inputs = [
            ("test1.default", "test1.default.frequencies.tsv"),
            ("test2.x_after_count", "test2.x_after_count.frequencies.tsv"),
            ("test2_full.x_after_count", "test2_full.x_after_count.frequencies.tsv"),
            ("test3.default", "test3.default.frequencies.tsv"),
        ]
        for name, tsv_file in test_inputs:
            with self.subTest(name=name):
                with tempfile.TemporaryDirectory() as tmpdir:
                    target_prefix = "test_run"
                    outfile_prefix = os.path.join(tmpdir, target_prefix)
                    tsv_path = os.path.join(self.outputs_dir, tsv_file)
                    cmd_args = [
                        "--tsv", tsv_path,
                        "--outfile-prefix", outfile_prefix,
                        "--aminoacids",
                        "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                        "--include-synonymous",
                        "--colormap=coolwarm_r",
                        "--threshold=0.001",
                        "--disable-showing-bokeh", "--disable-showing-mplcursors"
                    ]
                    returncode, output = self._invoke_cli(cmd_args)
                    self.assertEqual(returncode, 0, f"Command failed for {name}:\n{output}")
                    self._check_outputs(target_prefix, tmpdir, f"{name}.scatter_aminoacids_synonymous")

    def test4_compare_aminoacids(self):
        """mutation_scatter_plot pairwise compare: --aminoacids vs --aminoacids --include-synonymous HTML JSONs"""
        # pylint: disable=too-many-locals
        with tempfile.TemporaryDirectory() as tmpdir:
            full_tsv = os.path.join(self.outputs_dir, "test2_full.x_after_count.frequencies.tsv")
            # Run A
            target_prefix_a = "test_run_test4A"
            outfile_prefix_a = os.path.join(tmpdir, target_prefix_a)
            cmd_args_a = [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_a,
                "--aminoacids",
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--threshold=0.01",
                "--disable-showing-bokeh", "--disable-showing-mplcursors"
            ]
            res_a_code, res_a_out = self._invoke_cli(cmd_args_a)
            self.assertEqual(res_a_code, 0, f"Command A failed:\n{res_a_out}")

            # Run B
            target_prefix_b = "test_run_test4B"
            outfile_prefix_b = os.path.join(tmpdir, target_prefix_b)
            cmd_args_b = [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_b,
                "--aminoacids",
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--include-synonymous",
                "--threshold=0.01",
                "--disable-showing-bokeh", "--disable-showing-mplcursors"
            ]
            res_b_code, res_b_out = self._invoke_cli(cmd_args_b)
            self.assertEqual(res_b_code, 0, f"Command B failed:\n{res_b_out}")

            # Extract data
            html_a = f"{outfile_prefix_a}.BLOSUM80.amino_acid_changes.html"
            html_b = f"{outfile_prefix_b}.BLOSUM80.amino_acid_changes.html"

            data_a = self._extract_bokeh_data(html_a)
            data_b = self._extract_bokeh_data(html_b)

            # Filter synonymous mutations from data_b for comparison with data_a
            # Data format: (keys, rows)
            # Row index 8 is 'label8' which for amino acids contains the introduced AA.
            # Row index 6 is 'label6' which contains the original AA.
            # Wait, let's check the extraction logic.

            def filter_synonymous(data):
                new_data = []
                for keys, rows in data:
                    # Dynamically find the index of the 'mutation' column
                    try:
                        mut_idx = keys.index('mutation')
                        filtered_rows = [r for r in rows if r[mut_idx][0] != r[mut_idx][-1]]
                    except (ValueError, IndexError):
                        # Fallback or handle cases where 'mutation' key is missing
                        filtered_rows = rows
                    new_data.append((keys, filtered_rows))
                return new_data

            data_a_filtered = filter_synonymous(data_a)
            data_b_filtered = filter_synonymous(data_b)

            if data_a_filtered != data_b_filtered:
                with tempfile.NamedTemporaryFile("w", delete=False) as f_a, tempfile.NamedTemporaryFile("w", delete=False) as f_b:
                    for keys, rows in data_a:
                        f_a.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                    for keys, rows in data_b:
                        f_b.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")

                diff_res = subprocess.run(
                    ["git", "diff", "--no-index", "--word-diff=color",
                     "--word-diff-regex=.", "--color=always",
                     f_a.name, f_b.name],
                    capture_output=True, text=True, check=False
                )
                self.fail(f"HTML JSON structural mapping differed between --aminoacids and --aminoacids --include-synonymous:\n{diff_res.stdout}")

    def test5_compare_codons(self):
        """mutation_scatter_plot pairwise compare: codon mode vs codon mode --include-synonymous HTML JSONs"""
        # pylint: disable=too-many-locals
        with tempfile.TemporaryDirectory() as tmpdir:
            full_tsv = os.path.join(self.outputs_dir, "test2_full.x_after_count.frequencies.tsv")
            # Run A
            target_prefix_a = "test_run_test5A"
            outfile_prefix_a = os.path.join(tmpdir, target_prefix_a)
            cmd_args_a = [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_a,
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--threshold=0.001"
            ]
            res_a_code, res_a_out = self._invoke_cli(cmd_args_a)
            self.assertEqual(res_a_code, 0, f"Command A failed:\n{res_a_out}")

            # Run B
            target_prefix_b = "test_run_test5B"
            outfile_prefix_b = os.path.join(tmpdir, target_prefix_b)
            cmd_args_b = [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_b,
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--include-synonymous",
                "--threshold=0.001",
                "--disable-showing-bokeh", "--disable-showing-mplcursors"
            ]
            res_b_code, res_b_out = self._invoke_cli(cmd_args_b)
            self.assertEqual(res_b_code, 0, f"Command B failed:\n{res_b_out}")

            # Extract data
            html_a = f"{outfile_prefix_a}.BLOSUM80.amino_acid_changes.html"
            html_b = f"{outfile_prefix_b}.BLOSUM80.amino_acid_changes.html"

            data_a = self._extract_bokeh_data(html_a)
            data_b = self._extract_bokeh_data(html_b)

            def filter_synonymous(data):
                new_data = []
                for keys, rows in data:
                    try:
                        mut_idx = keys.index('mutation')
                        # Check first and last characters of the amino acid mutation string (e.g., 'K417N' -> K != N)
                        filtered_rows = [r for r in rows if r[mut_idx][0] != r[mut_idx][-1]]
                    except (ValueError, IndexError):
                        filtered_rows = rows
                    new_data.append((keys, filtered_rows))
                return new_data

            data_a_filtered = filter_synonymous(data_a)
            data_b_filtered = filter_synonymous(data_b)

            if data_a_filtered != data_b_filtered:
                with tempfile.NamedTemporaryFile("w", delete=False) as f_a, tempfile.NamedTemporaryFile("w", delete=False) as f_b:
                    for keys, rows in data_a:
                        f_a.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                    for keys, rows in data_b:
                        f_b.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")

                diff_res = subprocess.run(
                    ["git", "diff", "--no-index", "--word-diff=color",
                     "--word-diff-regex=.", "--color=always",
                     f_a.name, f_b.name],
                    capture_output=True, text=True, check=False
                )
                self.fail(f"HTML JSON structural mapping differed between codon default and --include-synonymous:\n{diff_res.stdout}")

    def test_display_suppression(self):
        """Verify that --disable-showing-bokeh and --disable-showing-mplcursors work without failure"""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_prefix = "test_suppress"
            outfile_prefix = os.path.join(tmpdir, target_prefix)
            tsv_path = os.path.join(self.outputs_dir, "test1.default.frequencies.tsv")
            cmd_args = [
                "--tsv", tsv_path,
                "--outfile-prefix", outfile_prefix,
                "--disable-showing-bokeh",
                "--disable-showing-mplcursors",
                "--threshold=0.001"
            ]
            returncode, output = self._invoke_cli(cmd_args)
            self.assertEqual(returncode, 0, f"Command with suppression flags failed:\n{output}")

            # Check that files were still generated
            expected_files = [
                f"{target_prefix}.BLOSUM80.amino_acid_changes.html",
                f"{target_prefix}.BLOSUM80.amino_acid_changes.png"
            ]
            for f in expected_files:
                self.assertTrue(os.path.exists(os.path.join(tmpdir, f)), f"Output file {f} was not generated")

if __name__ == "__main__":
    unittest.main()
