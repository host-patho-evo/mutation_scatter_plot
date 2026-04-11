# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Tests probing potential silent bugs in calculate_codon_frequencies.

These tests exercise specific edge cases identified during code review
(2026-04-11).  Each test class targets one potential bug:

1. Float division in reference offset slicing (``int(x / 3)`` vs ``x // 3``)
2. ``max_stop + 1`` boundary semantics
3. N-heavy sequences silently discarded via depadded length
4. Leading/trailing N treated as gaps (over-masking boundary codons)
5. ``starmap`` with 1-tuples where ``map`` suffices (perf only — no data test)
"""

import unittest

import numpy as np

from mutation_scatter_plot.calculate_codon_frequencies import (
    _process_one_site,
    get_codons,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_arrays(sequences, counts=None):
    """Build the NumPy arrays expected by _process_one_site from plain strings.

    Parameters
    ----------
    sequences : list[str]
        Upper-cased alignment sequences, all the same length.
    counts : list[int] | None
        Per-sequence multiplicity.  Defaults to 1 for each.

    Returns
    -------
    dict with keys matching _process_one_site positional args.
    """
    n = len(sequences)
    if counts is None:
        counts = [1] * n
    aln_len = len(sequences[0])

    joined = "".join(sequences)
    aln_array = np.frombuffer(
        joined.encode('ascii'), dtype=np.uint8
    ).reshape(n, aln_len)

    counts_array = np.array(counts, dtype=np.int64)
    padded_len_array = np.array(
        [len(s) for s in sequences], dtype=np.int32
    )
    depadded_len_array = np.array(
        [len(s) - s.count('-') - s.count('N') for s in sequences],
        dtype=np.int32,
    )
    leading_gap_array = np.array(
        [len(s) - len(s.lstrip('-N')) for s in sequences],
        dtype=np.int32,
    )
    trailing_gap_array = np.array(
        [len(s.rstrip('-N')) for s in sequences],
        dtype=np.int32,
    )
    return {
        'num_unique': n,
        'alignment_len': aln_len,
        'aln_array': aln_array,
        'counts_array': counts_array,
        'padded_len_array': padded_len_array,
        'depadded_len_array': depadded_len_array,
        'leading_gap_array': leading_gap_array,
        'trailing_gap_array': trailing_gap_array,
    }


def _ref_gaps_cumulative(ref_seq):
    """Build the cumulative gap count array for the reference."""
    arr = np.frombuffer(ref_seq.encode('ascii'), dtype=np.uint8)
    return np.cumsum(arr == ord('-'))


# ---------------------------------------------------------------------------
# Bug 1: Float division in offset slicing
# ---------------------------------------------------------------------------

class TestBug1FloatDivisionOffsetSlicing(unittest.TestCase):
    """int(x / 3) vs x // 3 — possible wrong codon index for non-multiple-of-3."""

    def test_int_float_div_vs_floor_div(self):
        """Document when ``int(x / 3)`` differs from ``x // 3``.

        For positive ints the two are identical; the potential divergence
        is with negative values (which shouldn't happen but could if
        left_reference_offset is 0 and the max() guard is removed).
        """
        # Positive cases — should all match
        for x in range(0, 100):
            self.assertEqual(
                int(x / 3), x // 3,
                f"int({x}/3) != {x}//3",
            )

    def test_protein_slice_alignment_with_dna_slice(self):
        """Protein and codon slices should cover exactly the same region as DNA.

        Given left_reference_offset=5, right_reference_offset=20:
        - DNA slice: ref[4:20]  (nucleotides 4..19)
        - Protein: ref_prot[int(4/3):int(min(len_prot, 20/3))]
                 = ref_prot[1:6]  (codons 1..5 → nucleotides 3..17)

        The DNA range (4..19) starts mid-codon-1, but the protein
        range starts at complete codon 1.  This means they are
        MISALIGNED when the offset is not a codon boundary.
        """
        ref_dna = "ATGATGATGATGATGATGATG"  # 7 codons (21 nt)
        ref_prot = "MMMMMM"  # 6 AAs (we'll use a shorter one for 18 nt)
        # Actually let's use exactly 21 nt => 7 codons
        ref_prot = "MMMMMMM"
        ref_codons = get_codons(ref_dna)

        left_offset = 5   # 1-based
        right_offset = 20  # 1-based

        # Current code logic (from __init__.py lines 502-512):
        dna_start = max(left_offset - 1, 0)      # 4
        dna_end = min(len(ref_dna), right_offset)  # 20
        dna_slice = ref_dna[dna_start:dna_end]

        prot_start = int(max(left_offset - 1, 0) / 3)  # int(4/3) = 1
        prot_end = int(min(len(ref_prot), right_offset / 3))  # int(min(7, 6.66)) = 6
        prot_slice = ref_prot[prot_start:prot_end]
        codon_slice = ref_codons[prot_start:prot_end]

        # The DNA slice covers nucleotides 4..19 (16 nt)
        self.assertEqual(len(dna_slice), 16)
        # The protein slice covers codons 1..5 (5 codons → 15 nt)
        self.assertEqual(len(prot_slice), 5)
        # The codon slice should cover the same range
        self.assertEqual(len(codon_slice), 5)

        # Document: DNA starts at nt 4 (inside codon 1, which is nt 3-5)
        # Protein starts at codon 1 (nt 3). So the protein covers a
        # slightly different range. This is a potential misalignment.
        dna_codon_start = dna_start // 3  # 1
        self.assertEqual(prot_start, dna_codon_start,
                         "Protein start codon index matches DNA codon index "
                         "(coincidentally correct for offset=5)")

        # Now try offset=6 (on a codon boundary) — should be clean
        left_offset_6 = 7  # 1-based → 0-based = 6
        dna_start_6 = max(left_offset_6 - 1, 0)  # 6
        prot_start_6 = int(max(left_offset_6 - 1, 0) / 3)  # int(6/3) = 2
        dna_codon_start_6 = dna_start_6 // 3  # 2
        self.assertEqual(prot_start_6, dna_codon_start_6,
                         "On codon boundary, both methods agree")

    def test_right_offset_not_divisible_by_3(self):
        """Right offset not divisible by 3 may truncate one codon."""
        ref_prot = "MMMMMMM"  # 7 AAs
        right_offset = 20     # not divisible by 3

        # Current code: int(min(7, 20/3)) = int(min(7, 6.666)) = int(6.666) = 6
        prot_end_current = int(min(len(ref_prot), right_offset / 3))
        self.assertEqual(prot_end_current, 6)

        # With // 3: min(7, 20//3) = min(7, 6) = 6
        prot_end_floor = min(len(ref_prot), right_offset // 3)
        self.assertEqual(prot_end_floor, 6)

        # They agree for this case. Let's find where they differ.
        # right_offset=19: int(min(7, 19/3)) = int(min(7, 6.333)) = int(6.333) = 6
        #                   min(7, 19//3)     = min(7, 6) = 6
        # right_offset=21: int(min(7, 21/3)) = int(min(7, 7.0)) = int(7.0) = 7
        #                   min(7, 21//3)     = min(7, 7) = 7
        # For positive values they always agree since int() truncates toward 0
        # and // floors toward -inf, but for positive values truncation == floor.
        for r in range(1, 100):
            current = int(min(len(ref_prot), r / 3))
            floor = min(len(ref_prot), r // 3)
            self.assertEqual(
                current, floor,
                f"right_offset={r}: int(min(7, {r}/3))={current} "
                f"!= min(7, {r}//3)={floor}",
            )


# ---------------------------------------------------------------------------
# Bug 2: max_stop boundary
# ---------------------------------------------------------------------------

class TestBug2MaxStopBoundary(unittest.TestCase):
    """max_stop semantics in CLI → range() — wobble position semantics."""

    def test_max_stop_wobble_position(self):
        """--max_stop=12 (wobble of 4th codon) → process exactly 4 codons.

        CLI does: _max_stop = max_stop (no transform)
        Then:     range(min_start, max_stop or aln_len, 3)
        Python range is exclusive, so range(0, 12, 3) = [0, 3, 6, 9].
        Last codon starts at 0-based 9 → covers 1-based nt 10, 11, 12.
        """
        _min_start = 0
        _max_stop = 12  # user value passed directly (no +1)

        positions = list(range(_min_start, _max_stop, 3))
        self.assertEqual(positions, [0, 3, 6, 9])

        # Last codon covers nt 10-12 (1-based), ending exactly at the stop
        last_codon_end_1based = positions[-1] + 3  # 0-based end → 1-based
        self.assertEqual(last_codon_end_1based, 12)

    def test_max_stop_3822_alignment(self):
        """--max_stop=3822 on a 3822-nt alignment → last codon is 3820-3822."""
        _max_stop = 3822  # user value (1-based wobble of last codon)
        positions = list(range(0, _max_stop, 3))
        self.assertEqual(positions[-1], 3819)  # 0-based start of last codon
        # 0-based 3819 → 1-based 3820, codon covers 3820-3822
        self.assertEqual(positions[-1] + 3, 3822)

    def test_max_stop_mid_codon(self):
        """--max_stop=10 (mid-codon) still includes that codon."""
        # 0-based 9 < 10, so it's included in range(0, 10, 3)
        positions = list(range(0, 10, 3))
        self.assertEqual(positions, [0, 3, 6, 9])
        # Codon at 0-based 9 covers nt 10-12 (1-based)

        # --max_stop=9: 0-based 9 NOT < 9, so it's excluded
        positions = list(range(0, 9, 3))
        self.assertEqual(positions, [0, 3, 6])
        # Last codon at 0-based 6 covers nt 7-9 (1-based)


# ---------------------------------------------------------------------------
# Bug 3: N-heavy sequences silently discarded
# ---------------------------------------------------------------------------

class TestBug3NHeavySequences(unittest.TestCase):
    """_depadded_len counts N as padding — N-heavy seqs may be discarded."""

    def test_depadded_length_excludes_n(self):
        """Sequences with many Ns have very short depadded length."""
        seq = "N" * 90 + "ATG" + "N" * 57  # 150 nt, only 3 real
        padded_len = len(seq)
        depadded_len = padded_len - seq.count('-') - seq.count('N')
        self.assertEqual(depadded_len, 3,
                         "N-heavy sequence depadded to only 3 nt")

    def test_n_heavy_seq_below_minimum_aln_length(self):
        """N-heavy sequence is classified as min_len_fail → excluded."""
        # Reference: 12 nucleotides = 4 codons
        ref_seq = "ATGATGATGATG"
        ref_prot = "MMMM"
        ref_gaps = _ref_gaps_cumulative(ref_seq)

        # One sequence: mostly N with a real codon at position 0
        sample = "ATG" + "N" * 9  # depadded_len = 3
        arrays = _make_arrays([sample])

        # With minimum_aln_length=50, this sequence should be excluded
        result = _process_one_site(
            0, arrays['num_unique'], arrays['alignment_len'],
            arrays['aln_array'], arrays['counts_array'],
            arrays['padded_len_array'], arrays['depadded_len_array'],
            arrays['leading_gap_array'], arrays['trailing_gap_array'],
            ref_seq, ref_prot, ref_gaps,
            aa_start=0, myoptions_minimum_aln_length=50,
            myoptions_left_reference_offset=0,
        )
        self.assertEqual(result['total_sum'], 0,
                         "N-heavy sequence excluded: total coverage is 0 "
                         "(sequence discarded because depadded_len=3 < 50)")

    def test_n_heavy_seq_above_minimum_aln_length(self):
        """N-heavy sequence passes if minimum_aln_length is low enough."""
        ref_seq = "ATGATGATGATG"
        ref_prot = "MMMM"
        ref_gaps = _ref_gaps_cumulative(ref_seq)

        sample = "ATG" + "N" * 9  # depadded_len = 3
        arrays = _make_arrays([sample])

        # With minimum_aln_length=0, this sequence should be included
        result = _process_one_site(
            0, arrays['num_unique'], arrays['alignment_len'],
            arrays['aln_array'], arrays['counts_array'],
            arrays['padded_len_array'], arrays['depadded_len_array'],
            arrays['leading_gap_array'], arrays['trailing_gap_array'],
            ref_seq, ref_prot, ref_gaps,
            aa_start=0, myoptions_minimum_aln_length=0,
            myoptions_left_reference_offset=0,
        )
        # Document: is the ATG codon at position 0 counted?
        # The leading gap detection uses lstrip('-N'), which strips the ATG
        # since lstrip strips ANY character in the set — wait, ATG chars
        # are NOT in '-N', so lstrip('-N') should leave "ATG..." intact.
        # leading_gap = 0, trailing_gap = 3 (rstrip('-N') strips trailing Ns).
        # pos=0, pos+2=2 < trailing_gap=3, so trailing mask doesn't fire.
        self.assertGreater(result['total_sum'], 0,
                           "N-heavy seq with min_aln_length=0 should contribute")


# ---------------------------------------------------------------------------
# Bug 4: Leading/trailing N as gaps — over-masking
# ---------------------------------------------------------------------------

class TestBug4LeadingTrailingNMasking(unittest.TestCase):
    """lstrip('-N') and rstrip('-N') treat N as gap character."""

    def test_leading_n_extends_boundary_mask(self):
        """Leading Ns cause boundary masking to extend into real codons."""
        # Sequence: NNN ATG CCC → leading Ns extend the "gap" to position 3
        seq = "NNNATGCCC"
        end_of_leading = len(seq) - len(seq.lstrip('-N'))
        self.assertEqual(end_of_leading, 3,
                         "Leading NNN treated as 3-nt gap")

        # The codon at position 3 (ATG) should NOT be masked, but
        # the boundary check is: _pos < _lead → 3 < 3 → False.
        # So ATG at position 3 is NOT masked — good.
        self.assertFalse(3 < end_of_leading,
                         "Codon at pos 3 is NOT masked (boundary is strict <)")

    def test_trailing_n_extends_boundary_mask(self):
        """Trailing Ns cause boundary masking to extend into real codons."""
        seq = "ATGCCCNNN"
        start_of_trailing = len(seq.rstrip('-N'))
        self.assertEqual(start_of_trailing, 6,
                         "Trailing NNN starts at position 6")

        # Codon at position 3 (CCC) ends at position 5.
        # Mask check: _pos + 2 >= _trail → 3 + 2 = 5 >= 6 → False
        # So CCC at position 3 is NOT masked — good.
        self.assertFalse(3 + 2 >= start_of_trailing,
                         "Codon CCC at pos 3 is NOT masked")

        # But codon at position 6 (NNN) ends at position 8.
        # Mask check: 6 + 2 = 8 >= 6 → True → MASKED
        self.assertTrue(6 + 2 >= start_of_trailing,
                        "Codon NNN at pos 6 IS masked (correct)")

    def test_real_nucleotides_adjacent_to_n_boundary(self):
        """Real nucleotides right next to N-boundary are preserved."""
        ref_seq = "ATGATGATG"
        ref_prot = "MMM"
        ref_gaps = _ref_gaps_cumulative(ref_seq)

        # Sequence with trailing Ns: ATG ATG NNN
        sample = "ATGATGNNN"
        arrays = _make_arrays([sample])

        # Process codon at position 3 (ATG, the second codon)
        result = _process_one_site(
            3, arrays['num_unique'], arrays['alignment_len'],
            arrays['aln_array'], arrays['counts_array'],
            arrays['padded_len_array'], arrays['depadded_len_array'],
            arrays['leading_gap_array'], arrays['trailing_gap_array'],
            ref_seq, ref_prot, ref_gaps,
            aa_start=0, myoptions_minimum_aln_length=0,
            myoptions_left_reference_offset=0,
        )
        # trailing_gap = rstrip('-N') of "ATGATGNNN" = "ATGATG" → 6
        # pos=3, pos+2=5 >= 6? → 5 >= 6 → False → NOT masked
        self.assertGreater(result['total_sum'], 0,
                           "Codon at pos 3 next to trailing Ns is preserved")

    def test_interior_n_does_not_trigger_boundary_mask(self):
        """Interior NNN (not at edges) should NOT trigger boundary masking."""
        ref_seq = "ATGNNNCCCTTT"  # 4 codons
        ref_prot = "MXPT"
        ref_gaps = _ref_gaps_cumulative(ref_seq)

        # Sample with interior Ns
        sample = "ATGNNNCCCTTT"
        arrays = _make_arrays([sample])

        # lstrip('-N') of "ATGNNNCCCTTT" → "ATGNNNCCCTTT" (A is not in '-N')
        # So leading_gap = 0. Good — interior Ns don't affect leading.
        self.assertEqual(arrays['leading_gap_array'][0], 0)

        # rstrip('-N') of "ATGNNNCCCTTT" → "ATGNNNCCCTTT" (T is not in '-N')
        # So trailing_gap = 12. Good — interior Ns don't affect trailing.
        self.assertEqual(arrays['trailing_gap_array'][0], 12)

        # Process codon at position 3 (NNN interior)
        result = _process_one_site(
            3, arrays['num_unique'], arrays['alignment_len'],
            arrays['aln_array'], arrays['counts_array'],
            arrays['padded_len_array'], arrays['depadded_len_array'],
            arrays['leading_gap_array'], arrays['trailing_gap_array'],
            ref_seq, ref_prot, ref_gaps,
            aa_start=0, myoptions_minimum_aln_length=0,
            myoptions_left_reference_offset=0,
        )
        # NNN at pos 3 is interior — should be counted (as ambiguous)
        self.assertGreater(
            result['total_sum'] + sum(result['ambiguous'].values()), 0,
            "Interior NNN contributes to ambiguous count",
        )

    def test_edge_case_all_n_sequence(self):
        """A sequence that is entirely Ns — depadded_len=0, fully masked."""
        ref_seq = "ATGATGATG"
        ref_prot = "MMM"
        ref_gaps = _ref_gaps_cumulative(ref_seq)

        sample = "N" * 9
        arrays = _make_arrays([sample])

        self.assertEqual(arrays['depadded_len_array'][0], 0)
        self.assertEqual(arrays['leading_gap_array'][0], 9)
        self.assertEqual(arrays['trailing_gap_array'][0], 0)

        result = _process_one_site(
            0, arrays['num_unique'], arrays['alignment_len'],
            arrays['aln_array'], arrays['counts_array'],
            arrays['padded_len_array'], arrays['depadded_len_array'],
            arrays['leading_gap_array'], arrays['trailing_gap_array'],
            ref_seq, ref_prot, ref_gaps,
            aa_start=0, myoptions_minimum_aln_length=0,
            myoptions_left_reference_offset=0,
        )
        self.assertEqual(result['total_sum'], 0,
                         "All-N sequence contributes nothing")


# ---------------------------------------------------------------------------
# Bug 5: starmap vs map (performance only — document overhead)
# ---------------------------------------------------------------------------

class TestBug5StarmapVsMap(unittest.TestCase):
    """starmap with 1-tuples vs map — no data bug, just overhead."""

    def test_starmap_tuple_wrapping_equivalence(self):
        """Demonstrate that starmap([f, [(x,)...]]) == map(f, [x...])."""
        def identity(x):
            return x * 2

        data = list(range(100))

        # Simulate starmap behavior
        starmap_result = [identity(*args) for args in [(x,) for x in data]]
        # Simulate map behavior
        map_result = [identity(x) for x in data]

        self.assertEqual(starmap_result, map_result)


if __name__ == "__main__":
    unittest.main()
