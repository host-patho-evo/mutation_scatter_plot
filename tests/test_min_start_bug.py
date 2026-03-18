
import unittest
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestMinStartBug(unittest.TestCase):
    def setUp(self):
        # Reference: ATGGCCGCT (Met Ala Ala)
        self.ref_seq = "ATGGCCGCT" # M A A
        with open("ref.fasta", "w") as f:
            SeqIO.write(SeqRecord(Seq(self.ref_seq), id="ref"), f, "fasta")

        # Alignment: All codons changed to ensure they appear in frequencies.tsv
        # ATG -> GTG (Val), GCC -> GTC (Val), GCT -> CGT (Arg)
        self.aln_seq = "GTGGTCCGT" # V V R
        with open("aln.fasta", "w") as f:
            SeqIO.write(SeqRecord(Seq(self.aln_seq), id="10x"), f, "fasta")

    def tearDown(self):
        for f in ["ref.fasta", "aln.fasta", "test.frequencies.tsv", "test.frequencies.unchanged_codons.tsv", "test.frequencies.count"]:
            if os.path.exists(f):
                os.remove(f)

    def run_calc(self, min_start):
        if os.path.exists("test.frequencies.tsv"):
            os.remove("test.frequencies.tsv")

        cmd = [
            "calculate_codon_frequencies",
            "--reference-infile=ref.fasta",
            "--alignment-file=aln.fasta",
            "--outfile-prefix=test.frequencies",
            "--x-after-count",
            f"--min_start={min_start}",
            "--minimum-alignments-length=5"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return None
        if os.path.exists("test.frequencies.tsv"):
            try:
                return pd.read_csv("test.frequencies.tsv", sep='\t', header=None)
            except pd.errors.EmptyDataError:
                return None
        return None

    def test_min_start_sync(self):
        # min_start=4 is the second codon (index 3, GCC)
        df = self.run_calc(4)
        self.assertIsNotNone(df)

        # Row 0 should be position 2, original_aa 'A'
        # Bug was: reported 'M' (index 0) instead of 'A' (index 1)
        self.assertEqual(df.iloc[0, 2], 'A', f"Expected original_aa 'A' at position 2, got {df.iloc[0, 2]}")
        self.assertEqual(df.iloc[0, 5], 'GCC', f"Expected original_codon 'GCC' at position 2, got {df.iloc[0, 5]}")

if __name__ == "__main__":
    unittest.main()
