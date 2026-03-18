
import unittest
import os
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class TestMinStartGaps(unittest.TestCase):
    def setUp(self):
        # Reference: ---ATGGCC (Gap, Met, Ala)
        # Padded position 1: --- (Gap)
        # Padded position 2: ATG (Met)
        # Padded position 3: GCC (Ala)
        # Depadded position 1: ATG (Met)
        # Depadded position 2: GCC (Ala)
        self.ref_seq = "---ATGGCC"
        with open("ref_gaps.fasta", "w") as f:
            SeqIO.write(SeqRecord(Seq(self.ref_seq), id="ref"), f, "fasta")

        # Alignment: Changed Met->Val, Ala->Val
        self.aln_seq = "---GTGGTC"
        with open("aln_gaps.fasta", "w") as f:
            SeqIO.write(SeqRecord(Seq(self.aln_seq), id="10x"), f, "fasta")

    def tearDown(self):
        for f in ["ref_gaps.fasta", "aln_gaps.fasta", "test_gaps.frequencies.tsv", "test_gaps.frequencies.unchanged_codons.tsv", "test_gaps.frequencies.count"]:
            if os.path.exists(f):
                os.remove(f)

    def run_calc(self, min_start):
        # Ensure it is installed
        cmd = [
            "calculate_codon_frequencies",
            "--reference-infile=ref_gaps.fasta",
            "--alignment-file=aln_gaps.fasta",
            "--outfile-prefix=test_gaps.frequencies",
            "--x-after-count",
            f"--min_start={min_start}",
            "--minimum-alignments-length=1"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if os.path.exists("test_gaps.frequencies.tsv") and os.path.getsize("test_gaps.frequencies.tsv") > 0:
            return pd.read_csv("test_gaps.frequencies.tsv", sep='\t', header=None)
        if os.path.exists("test_gaps.frequencies.unchanged_codons.tsv") and os.path.getsize("test_gaps.frequencies.unchanged_codons.tsv") > 0:
            return pd.read_csv("test_gaps.frequencies.unchanged_codons.tsv", sep='\t', header=None)
        return None

    def test_min_start_gap_sync(self):
        # Start at GCC (nucleotide index 7, so --min_start=7)
        df = self.run_calc(7)
        self.assertIsNotNone(df, "Should have found data in tsv or unchanged_codons.tsv")

        print("Rows with min_start=7:")
        print(df.iloc[:, [0, 1, 2, 5]].values.tolist())

        self.assertEqual(df.iloc[0, 0], 3, "Padded position should be 3")
        self.assertEqual(df.iloc[0, 1], 2, f"Expected depadded position 2, got {df.iloc[0, 1]}")
        self.assertEqual(df.iloc[0, 2], 'A', "Original AA should be 'A'")

if __name__ == "__main__":
    unittest.main()
