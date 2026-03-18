
import unittest
from mutation_scatter_plot.calculate_codon_frequencies import get_codons

class TestGetCodons(unittest.TestCase):
    def test_padded_multiple_of_3(self):
        # Current logic just slices every 3 chars if len % 3 == 0
        seq = "ATGGCC---"
        expected = ["ATG", "GCC", "---"]
        self.assertEqual(get_codons(seq), expected)

    def test_elif_branch_fix(self):
        # seq length 7 (not multiple of 3), depadded length 6 (multiple of 3)
        # ATG-GCC
        seq = "ATG-GCC"
        # Corrected logic should return codons from depadded sequence
        expected = ["ATG", "GCC"]
        self.assertEqual(get_codons(seq), expected)

    def test_non_multiple_of_3(self):
        seq = "ATGC"
        with self.assertRaises(ValueError):
            get_codons(seq)

if __name__ == "__main__":
    unittest.main()
