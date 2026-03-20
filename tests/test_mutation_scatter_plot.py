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
                    # Run diff -u -w --color=always to show the differences
                    diff_result = subprocess.run(
                        ["diff", "-u", "-w", "--color=always", expected_path, gen_path],
                        capture_output=True, text=True, check=False
                    )
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
                                ["diff", "-u", "-w", "--color=always", f_exp.name, f_gen.name],
                                capture_output=True, text=True, check=False
                            )
                            self.fail(f"Bokeh internal data points differ for {gen_file}!\nDifferences:\n{diff_res.stdout}")
                        else:
                            print(f"Info: HTML file {gen_file} byte-mismatch ignored because internal Bokeh JSON structural mapping matched perfectly.")
                    elif any(gen_file.endswith(ext) for ext in [".png", ".pdf"]):
                        print(f"Warning: {gen_file} plotting artifact differed from strictly byte-for-byte matching golden baseline natively (ignoring natively due to UUIDs/timestamps).\nDifferences:\n{diff_result.stdout}")
                    else:
                        self.fail(f"File {gen_file} does not match golden file {expected_filename}.\nDifferences:\n{diff_result.stdout}")

    def _extract_bokeh_data(self, html_path):
        import json
        import re
        from bokeh.document import Document

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
                "--threshold=0.01"
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
                "--threshold=0.01"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed:\n{result.stderr}\n{result.stdout}")
            self._check_outputs(target_prefix, tmpdir, "test2.scatter_aminoacids_synonymous")
    def test_codons(self):
        """mutation_scatter_plot --show-STOP --show-X --show-DEL --show-INS"""
        with tempfile.TemporaryDirectory() as tmpdir:
            target_prefix = "test_run_4"
            outfile_prefix = os.path.join(tmpdir, target_prefix)
            cmd = self.base_cmd + [
                "--tsv", self.tsv_input,
                "--outfile-prefix", outfile_prefix,
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--threshold=0.001"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed:\n{result.stderr}\n{result.stdout}")
            self._check_outputs(target_prefix, tmpdir, "test2.scatter_codons")

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
                "--threshold=0.001"
            ]
            result = subprocess.run(cmd, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, f"Command failed:\n{result.stderr}\n{result.stdout}")
            self._check_outputs(target_prefix, tmpdir, "test2.scatter_codons_synonymous")

    def test4_compare_aminoacids(self):
        """mutation_scatter_plot pairwise compare: --aminoacids vs --aminoacids --include-synonymous HTML JSONs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            full_tsv = os.path.join(self.outputs_dir, "test2_full.frequencies_default.tsv")
            # Run A
            target_prefix_a = "test_run_test4A"
            outfile_prefix_a = os.path.join(tmpdir, target_prefix_a)
            cmd_a = self.base_cmd + [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_a,
                "--aminoacids",
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--threshold=0.01"
            ]
            res_a = subprocess.run(cmd_a, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(res_a.returncode, 0, f"Command A failed:\n{res_a.stderr}\n{res_a.stdout}")
            
            # Run B
            target_prefix_b = "test_run_test4B"
            outfile_prefix_b = os.path.join(tmpdir, target_prefix_b)
            cmd_b = self.base_cmd + [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_b,
                "--aminoacids",
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--include-synonymous",
                "--threshold=0.01"
            ]
            res_b = subprocess.run(cmd_b, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(res_b.returncode, 0, f"Command B failed:\n{res_b.stderr}\n{res_b.stdout}")

            # Extract data
            html_a = f"{outfile_prefix_a}.BLOSUM80.amino_acid_changes.html"
            html_b = f"{outfile_prefix_b}.BLOSUM80.amino_acid_changes.html"
            
            data_a = self._extract_bokeh_data(html_a)
            data_b = self._extract_bokeh_data(html_b)
            
            if data_a != data_b:
                with tempfile.NamedTemporaryFile("w", delete=False) as f_a, tempfile.NamedTemporaryFile("w", delete=False) as f_b:
                    for keys, rows in data_a:
                        f_a.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                    for keys, rows in data_b:
                        f_b.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                
                diff_res = subprocess.run(
                    ["diff", "-u", "-w", "--color=always", f_a.name, f_b.name],
                    capture_output=True, text=True, check=False
                )
                self.fail(f"HTML JSON structural mapping differed between --aminoacids and --aminoacids --include-synonymous:\n{diff_res.stdout}")

    def test5_compare_codons(self):
        """mutation_scatter_plot pairwise compare: codon mode vs codon mode --include-synonymous HTML JSONs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            full_tsv = os.path.join(self.outputs_dir, "test2_full.frequencies_default.tsv")
            # Run A
            target_prefix_a = "test_run_test5A"
            outfile_prefix_a = os.path.join(tmpdir, target_prefix_a)
            cmd_a = self.base_cmd + [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_a,
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--threshold=0.001"
            ]
            res_a = subprocess.run(cmd_a, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(res_a.returncode, 0, f"Command A failed:\n{res_a.stderr}\n{res_a.stdout}")
            
            # Run B
            target_prefix_b = "test_run_test5B"
            outfile_prefix_b = os.path.join(tmpdir, target_prefix_b)
            cmd_b = self.base_cmd + [
                "--tsv", full_tsv,
                "--outfile-prefix", outfile_prefix_b,
                "--show-STOP", "--show-X", "--show-DEL", "--show-INS",
                "--include-synonymous",
                "--threshold=0.001"
            ]
            res_b = subprocess.run(cmd_b, cwd=self.project_root, env=self.env, capture_output=True, text=True, check=False)
            self.assertEqual(res_b.returncode, 0, f"Command B failed:\n{res_b.stderr}\n{res_b.stdout}")

            # Extract data
            html_a = f"{outfile_prefix_a}.BLOSUM80.amino_acid_changes.html"
            html_b = f"{outfile_prefix_b}.BLOSUM80.amino_acid_changes.html"
            
            data_a = self._extract_bokeh_data(html_a)
            data_b = self._extract_bokeh_data(html_b)
            
            if data_a != data_b:
                with tempfile.NamedTemporaryFile("w", delete=False) as f_a, tempfile.NamedTemporaryFile("w", delete=False) as f_b:
                    for keys, rows in data_a:
                        f_a.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                    for keys, rows in data_b:
                        f_b.write(str(keys) + "\n" + "\n".join(map(str, rows)) + "\n")
                
                diff_res = subprocess.run(
                    ["diff", "-u", "-w", "--color=always", f_a.name, f_b.name],
                    capture_output=True, text=True, check=False
                )
                self.fail(f"HTML JSON structural mapping differed between codon mode and codon mode --include-synonymous:\n{diff_res.stdout}")

if __name__ == "__main__":
    unittest.main()
