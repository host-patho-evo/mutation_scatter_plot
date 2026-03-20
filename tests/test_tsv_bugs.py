
import unittest
import os
import subprocess
import pandas as pd

class TestTsvBugs(unittest.TestCase):
    def tearDown(self):
        # The outfile_prefix is modified by the matrix and colormap in the real app
        expected_rendered = "old.BLOSUM80.amino_acid_changes.actually_rendered.tsv"
        files_to_remove = [
            "old.frequencies.tsv", "old.frequencies.unchanged_codons.tsv",
            expected_rendered,
            "old.BLOSUM80.amino_acid_changes.aa.frequencies.colors.tsv",
            "old.BLOSUM80.amino_acid_changes.codon.frequencies.colors.tsv",
            "old.BLOSUM80.amino_acid_changes.png",
            "old.BLOSUM80.amino_acid_changes.pdf",
            "old.BLOSUM80.amino_acid_changes.html"
        ]
        for f in files_to_remove:
            if os.path.exists(f):
                os.remove(f)

    def test_legacy_tsv_loading(self):
        # 6-column headerless TSV
        tsv_content = "1\tM\tV\t0.5\tATG\tGTG\n"
        with open("old.frequencies.tsv", "w") as f:
            f.write(tsv_content)

        with open("old.frequencies.unchanged_codons.tsv", "w") as f:
            f.write("1\tM\tM\t0.5\tATG\tATG\n")

        env = os.environ.copy()
        env["MPLBACKEND"] = "agg"

        cmd = [
            "mutation_scatter_plot",
            "--tsv=old.frequencies.tsv",
            "--outfile-prefix=old",
            "--xmin=1",
            "--xmax=10"
        ]

        subprocess.run(cmd, capture_output=True, text=True, env=env)

        expected_rendered = "old.BLOSUM80.amino_acid_changes.actually_rendered.tsv"
        self.assertTrue(os.path.exists(expected_rendered))

        rendered_df = pd.read_csv(expected_rendered, sep='\t', header=None)
        # Bug 1: header=0 caused loss of the only row
        self.assertEqual(len(rendered_df), 1)
        # Bug 2: 'padded_position' KeyError
        self.assertEqual(rendered_df.iloc[0, 0], 1)
        self.assertEqual(rendered_df.iloc[0, 1], 'M')

if __name__ == "__main__":
    unittest.main()
