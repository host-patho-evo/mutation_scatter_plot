import sys
import unittest
import io

# Mock sys.argv before importing because the script invokes OptionParser.parse_args() at the module level
sys.argv = ['reversecomplement_reads_on_minus.py']

from scripts.reversecomplement_reads_on_minus import (
    string_reverse_complement, 
    shorten_sequence, 
    parse_input,
    myoptions
)

class TestReverseComplementReads(unittest.TestCase):

    def test_string_reverse_complement_basic(self):
        """Test basic native C-level reversing of DNA strings."""
        # Simple
        self.assertEqual(string_reverse_complement("ATGC"), "GCAT")
        self.assertEqual(string_reverse_complement("AAAA"), "TTTT")
        
    def test_string_reverse_complement_ambiguous(self):
        """Test ambiguous IUPAC nucleotide calls natively reversing."""
        # N -> N, Purines(R) <-> Pyrimidines(Y), Strong(S) <-> Strong(S), Weak(W) <-> Weak(W)
        self.assertEqual(string_reverse_complement("N"), "N")
        self.assertEqual(string_reverse_complement("R"), "Y")
        self.assertEqual(string_reverse_complement("Y"), "R")
        self.assertEqual(string_reverse_complement("S"), "S")
        self.assertEqual(string_reverse_complement("W"), "W")
        self.assertEqual(string_reverse_complement("M"), "K")
        self.assertEqual(string_reverse_complement("K"), "M")
        
        # Test full IUPAC geometry string
        orig = "ATCGNRY"
        # Translated: TAGCNYR
        # RC        : RYNCGAT
        expected_rc = "RYNCGAT" 
        self.assertEqual(string_reverse_complement(orig), expected_rc)

    def test_shorten_sequence_padding_no_slice(self):
        """Test sequence padding without strict length boundaries."""
        seq, start, stop = shorten_sequence(
            sequence="ATGC",
            min_start=0,
            number_of_leading_dashes=2,
            number_of_trailing_dashes=3,
            max_stop=0,
            aln_stop=10,
            qseq="ATGC",
            aln_start_qseq=1,
            aln_stop_qseq=4
        )
        self.assertEqual(seq, "--ATGC---")
        
    def test_shorten_sequence_padding_with_slice(self):
        """Test sequence padding perfectly slicing to alignment widths."""
        seq, start, stop = shorten_sequence(
            sequence="ATGC",
            min_start=1,  # Slice off the first character of the resulting dash array
            number_of_leading_dashes=2,
            number_of_trailing_dashes=2,
            max_stop=7,   # Slice ending width
            aln_stop=10,
            qseq="ATGC",
            aln_start_qseq=1,
            aln_stop_qseq=4
        )
        # Native: "--ATGC--" (len 8)
        # Sliced [1:7]: "-ATGC-" (len 6)
        self.assertEqual(seq, "-ATGC-")

    def test_parse_input_minus_strand_flip(self):
        """Test the integration pipeline switching minus orientation and rewriting IDs."""
        # Input simulating a minus strand 
        # >ID ALN_START ALN_STOP MINUS ... QSEQ_START QSEQ_STOP ... SSEQ
        fasta_mock = io.StringIO(
            ">Test1 2000 1000 minus 0 0 0 0 0 0 1 100 0 0 0 0\n"
            "ATGC\n"
        )
        
        # Mock reference sequence
        ref_seq = "A" * 3000
        
        # Generator
        myoptions.x_after_count = False
        myoptions.min_count = 0
        
        sequences = list(parse_input(fasta_mock, ref_seq, "fasta"))
        self.assertEqual(len(sequences), 1)
        
        header, seq = sequences[0]
        # Verify orientation word was successfully overridden from minus natively to plus!
        self.assertIn("plus", header.split())
        self.assertNotIn("minus", header.split())
        
        # Verify sequence was actively reverse complemented
        self.assertIn("GCAT", seq)

    def test_parse_input_plus_strand_passthrough(self):
        """Test that plus strand objects bypass sequence reversion."""
        fasta_mock = io.StringIO(
            ">Test1 1000 2000 plus 0 0 0 0 0 0 1 100 0 0 0 0\n"
            "ATGC\n"
        )
        
        ref_seq = "A" * 3000
        
        sequences = list(parse_input(fasta_mock, ref_seq, "fasta"))
        self.assertEqual(len(sequences), 1)
        
        header, seq = sequences[0]
        # Should stay plus
        self.assertIn("plus", header.split())
        
        # Forward sequence geometry maintained
        self.assertIn("ATGC", seq)

    def test_docstring_development_sequence_minus_1(self):
        """Test the real minus strand read documented in script comments: A00808:1538:HWL7GDRX3:2:2101:7491:1125"""
        seq = "AGAAAGTACTACTAATCTGTATGGTAGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAGAGTAACAATTAGGACCTGCAACACAAATACAAGGTTTGATACCGGCCAGATAGATATCAGTTGAAATATATCTCTCAAAAGGTTTCAGTTTAGACTTCCTAAACAATATATACCAGTAATCATAAATACCACTATGCTTAGAATCAAGCATGTTAGAATTCCAAGCTATATCGCAGCCTTTTAAATCATCTGGTAATTTTTAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATGAGTAAAGCAGTGATCATTTATTTTAGTAGGAGACACTCC"
        rc_seq = string_reverse_complement("".join(seq.split()))
        # Ends with GGAGTGTCTCCTACTAA...
        self.assertTrue(rc_seq.startswith("GGAGTGTCTCCTACTAAAAT"))
        self.assertEqual(len(rc_seq), 402)
    
    def test_docstring_development_sequence_minus_2(self):
        """Test the second real minus strand read documented: A00808:1538:HWL7GDRX3:2:2101:24831:1219"""
        seq = "AGAAAGTACTACTACTCTGTATGGTTGGTGACCAACACCATAAGTGGGTCGGAAACCATATGATTGTAAAGGAAAGTAACAATTAGGACCTTTA---CCTTTACAAGGTTTGTTACCGGCCTGATAGATTTCAGTTGAAATATCTCTCTCAAAAGGTTTGAGTTTAGACTTCCTAAACAATCTATACCAGTAATCATAATTACCACTATGCTTAGAATCAAGCTTGTTAGAATTCCAAGCTATAACGCAGCCTGTAAAATCATCTGGTAATTTATAATTATAATCAGCAATATTTCCAGTTTGCCCTGGAGCGATTTGGCTGACTTCATTACCTTTAATTACAAATGAATCTGCATAGACATTAGTAAAGCAGAGATCATTTAATTTAGTAGGAGACACTCC"
        rc_seq = string_reverse_complement("".join(seq.split()))
        self.assertTrue(rc_seq.startswith("GGAGTGTCTCCTACTAAATT"))

    def test_docstring_development_sequence_plus_1(self):
        """Test the plus strand read: A00808:1538:HWL7GDRX3:2:2101:21296:1204"""
        seq = "GGAGTGTCTCCTACTAAATTAAATGATCTCTGCTTTACTAATGTCTATGCAGATTCATTTGTAATTAAAGGTAATGAAGTCAGCCAAATCGCTCCAGGGCAAACTGGAAATATTGCTGATTATAATTATAAATTACCAGATGATTTTACAGGCTGCGTTATAGCTTGGAATTCTAACAAGCTTGATTCTAAGCATAGTGGTAATTATGATTACTGGTATAGATTGTTTAGGAAGTCTAAACTCAAACCTTTTGAGAGAGATATTTCAACTGAAATCTATCAGGCCGGTAACAAACCTTGTAAAGG---TAAAGGTCCTAATTGTTACTTTCCTTTACAATCATATGGTTTCCGACCCACTTATGGTGTTGGTCACCAACCATACAGAGTAGTAGTACTTTCT"
        # It's plus strand so reverse_complement is NOT called in normal logic, but if manually flipped:
        rc_seq = string_reverse_complement("".join(seq.split()))
        self.assertTrue(rc_seq.endswith("GTAGGAGACACTCC"))

if __name__ == '__main__':
    unittest.main()
