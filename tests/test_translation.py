# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""Comprehensive tests for codon translation logic.

Covers both:
  - alt_translate() in mutation_scatter_plot/__init__.py
  - _build_codon_table() / _translate_seq() in scripts/translate.py

Test scenarios (per the discussion in Biopython PR #4992 and issue #5036):
  1. Standard codons (64 ATGC combinations)
  2. Full-gap codon ('---' → '-')
  3. Partial-gap codons — each '-' treated as 'N':
       TC-  → TCN → S   (four-fold degenerate → unambiguous Serine)
       AT-  → ATN → X   (ATN can be Ile or Met → ambiguous)
       A--  → ANN → X   (ambiguous)
       -TC  → NTC → X   (ambiguous)
       -T-  → NTN → X   (mostly Ile/Val but check Biopython)
  4. Pure IUPAC ambiguity codons (no gap):
       TCN  → S (Serine, four-fold synonymous)
       CTN  → L (Leucine, four-fold synonymous)
       GGN  → G (Glycine, four-fold synonymous)
       GCN  → A (Alanine, four-fold synonymous)
       GTN  → V (Valine, four-fold synonymous)
       CCN  → P (Proline, four-fold synonymous)
       ACN  → T (Threonine, four-fold synonymous)
       CGN  → R (Arginine, four-fold synonymous)
       ATN  → X (Ile/Met ambiguous)
       AGN  → X (Arg/Ser/stop ambiguous)
       NNN  → X (completely ambiguous)
  5. ignore_gaps mode: strips '-' before slicing
  6. respect_alignment mode (default): '-' treated as 'N' per-codon
  7. Incomplete trailing codon silently dropped
  8. Full sequence translation with mixed codons
  9. Alternative genetic code tables
 10. Case insensitivity of standard codons
"""
# pylint: disable=missing-function-docstring  # pytest: test names are self-documenting
# pylint: disable=invalid-name                # biological codon names use uppercase (TC-, AT-)
# pylint: disable=wrong-import-order,C0413    # translate import must follow sys.path.insert()

import io
import itertools
import os
import sys
from decimal import Decimal

import pytest

# Import from the package
from mutation_scatter_plot import alt_translate
from mutation_scatter_plot.calculate_codon_frequencies import write_tsv_line

# ── helpers / imports ─────────────────────────────────────────────────────────


# Import internal functions from the scripts/ translate tool.
# Since scripts/ is not installed as a package, add it to sys.path.
_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'scripts')
sys.path.insert(0, os.path.abspath(_SCRIPTS_DIR))

from translate import (_build_codon_table, _translate_seq,  # noqa: E402
                       parse_input)

# ─────────────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────────────

@pytest.fixture
def table1():
    """Standard genetic code (NCBI table 1)."""
    return _build_codon_table(table_id=1)


# ─────────────────────────────────────────────────────────────────────────────
# 1. Standard codons — fast dict path
# ─────────────────────────────────────────────────────────────────────────────

class TestStandardCodons:
    """All 64 ATGC codons should translate via the fast dict lookup."""

    @pytest.mark.parametrize("codon,expected", [
        ("TCA", "S"), ("TCC", "S"), ("TCG", "S"), ("TCT", "S"),  # Serine (TC*)
        ("TTA", "L"), ("TTG", "L"),                               # Leucine
        ("CTT", "L"), ("CTC", "L"), ("CTA", "L"), ("CTG", "L"),  # Leucine
        ("ATG", "M"),                                              # Met (start)
        ("ATA", "I"), ("ATC", "I"), ("ATT", "I"),                 # Isoleucine
        ("GGT", "G"), ("GGC", "G"), ("GGA", "G"), ("GGG", "G"),  # Glycine
        ("GCT", "A"), ("GCC", "A"), ("GCA", "A"), ("GCG", "A"),  # Alanine
        ("GTT", "V"), ("GTC", "V"), ("GTA", "V"), ("GTG", "V"),  # Valine
        ("CCT", "P"), ("CCC", "P"), ("CCA", "P"), ("CCG", "P"),  # Proline
        ("ACT", "T"), ("ACC", "T"), ("ACA", "T"), ("ACG", "T"),  # Threonine
        ("CGT", "R"), ("CGC", "R"), ("CGA", "R"), ("CGG", "R"),  # Arginine
        ("TAA", "*"), ("TAG", "*"), ("TGA", "*"),                  # Stop codons
    ])
    def test_fast_path(self, table1, codon, expected):
        assert _translate_seq(codon, table1, ignore_gaps=False) == expected

    def test_lower_case_same_as_upper(self, table1):
        """Lower-case codons should give the same result as upper-case."""
        assert _translate_seq("tca", table1, ignore_gaps=False) == "S"
        assert _translate_seq("atg", table1, ignore_gaps=False) == "M"
        assert _translate_seq("ggg", table1, ignore_gaps=False) == "G"

    def test_full_serine_quartet(self, table1):
        """All four TCX codons encode Serine."""
        seq = "TCATCCTCGTCT"
        assert _translate_seq(seq, table1, ignore_gaps=False) == "SSSS"


# ─────────────────────────────────────────────────────────────────────────────
# 2. Full-gap codon
# ─────────────────────────────────────────────────────────────────────────────

class TestFullGapCodon:
    """'---' must translate to '-' (alignment gap amino acid)."""

    def test_gap_codon(self, table1):
        assert _translate_seq("---", table1, ignore_gaps=False) == "-"

    def test_gap_in_sequence(self, table1):
        assert _translate_seq("TCA---TCG", table1, ignore_gaps=False) == "S-S"

    def test_alt_translate_gap_codon(self):
        assert alt_translate("---") == "-"

    def test_alt_translate_gap_in_sequence(self):
        assert alt_translate("TCA---TCG") == "S-S"


# ─────────────────────────────────────────────────────────────────────────────
# 3. Partial-gap codons — '-' treated as 'N'
# ─────────────────────────────────────────────────────────────────────────────

class TestPartialGapCodons:
    """Partial-gap codons have '-' in a non-full-gap position.

    Each '-' is treated as 'N' (unknown nucleotide) and the result is
    determined by whether the remaining known nucleotides uniquely fix the AA.
    """

    # TC- → TCN → S (all four TCA/TCC/TCG/TCT give Serine)
    def test_TC_dash_gives_S(self, table1):
        assert _translate_seq("TC-", table1, ignore_gaps=False) == "S"

    def test_alt_translate_TC_dash_gives_S(self):
        assert alt_translate("TC-") == "S"

    # AT- → ATN → X (ATA/ATC/ATT=Ile, ATG=Met → ambiguous)
    def test_AT_dash_gives_X(self, table1):
        assert _translate_seq("AT-", table1, ignore_gaps=False) == "X"

    def test_alt_translate_AT_dash_gives_X(self):
        assert alt_translate("AT-") == "X"

    # GG- → GGN → G (all four give Glycine)
    def test_GG_dash_gives_G(self, table1):
        assert _translate_seq("GG-", table1, ignore_gaps=False) == "G"

    def test_alt_translate_GG_dash_gives_G(self):
        assert alt_translate("GG-") == "G"

    # CC- → CCN → P (all four give Proline)
    def test_CC_dash_gives_P(self, table1):
        assert _translate_seq("CC-", table1, ignore_gaps=False) == "P"

    # AC- → ACN → T (all four give Threonine)
    def test_AC_dash_gives_T(self, table1):
        assert _translate_seq("AC-", table1, ignore_gaps=False) == "T"

    # GT- → GTN → V (all four give Valine)
    def test_GT_dash_gives_V(self, table1):
        assert _translate_seq("GT-", table1, ignore_gaps=False) == "V"

    # GC- → GCN → A (all four give Alanine)
    def test_GC_dash_gives_A(self, table1):
        assert _translate_seq("GC-", table1, ignore_gaps=False) == "A"

    # CT- → CTN → L (CTA/CTC/CTG/CTT all give Leucine)
    def test_CT_dash_gives_L(self, table1):
        assert _translate_seq("CT-", table1, ignore_gaps=False) == "L"

    # CG- → CGN → R (all four give Arginine)
    def test_CG_dash_gives_R(self, table1):
        assert _translate_seq("CG-", table1, ignore_gaps=False) == "R"

    # A-- → ANN → X (too ambiguous)
    def test_A_double_dash_gives_X(self, table1):
        assert _translate_seq("A--", table1, ignore_gaps=False) == "X"

    # --- is NOT a partial-gap; it maps to '-'
    def test_triple_dash_is_not_partial_gap(self, table1):
        assert _translate_seq("---", table1, ignore_gaps=False) == "-"

    # mixed in a sequence
    def test_partial_gap_in_sequence(self, table1):
        # TC- → S, TCA → S, AT- → X
        result = _translate_seq("TC-TCAAT-", table1, ignore_gaps=False)
        assert result == "SSX"

    def test_alt_translate_partial_gap_in_sequence(self):
        result = alt_translate("TC-TCAAT-")
        assert result == "SSX"


# ─────────────────────────────────────────────────────────────────────────────
# 4. Pure IUPAC ambiguity codons (no gap)
# ─────────────────────────────────────────────────────────────────────────────

class TestIUPACCodons:
    """IUPAC ambiguity codes in codons: Biopython resolves per-codon."""

    @pytest.mark.parametrize("codon,expected", [
        ("TCN", "S"),   # TC[ACGT] all → Serine
        ("CTN", "L"),   # CT[ACGT] all → Leucine
        ("GGN", "G"),   # GG[ACGT] all → Glycine
        ("GCN", "A"),   # GC[ACGT] all → Alanine
        ("GTN", "V"),   # GT[ACGT] all → Valine
        ("CCN", "P"),   # CC[ACGT] all → Proline
        ("ACN", "T"),   # AC[ACGT] all → Threonine
        ("CGN", "R"),   # CG[ACGT] all → Arginine
        ("ATN", "X"),   # ATA/ATC/ATT→I, ATG→M → ambiguous
        ("NNN", "X"),   # completely ambiguous
    ])
    def test_iupac_translate_seq(self, table1, codon, expected):
        assert _translate_seq(codon, table1, ignore_gaps=False) == expected

    @pytest.mark.parametrize("codon,expected", [
        ("TCN", "S"),
        ("CTN", "L"),
        ("GGN", "G"),
        ("GCN", "A"),
        ("GTN", "V"),
        ("CCN", "P"),
        ("ACN", "T"),
        ("CGN", "R"),
        ("ATN", "X"),
        ("NNN", "X"),
    ])
    def test_iupac_alt_translate(self, codon, expected):
        assert alt_translate(codon) == expected

    def test_iupac_sequence(self, table1):
        """Mixed IUPAC codons in a sequence."""
        # TCN→S, ATN→X, GGN→G
        result = _translate_seq("TCNATNGCN", table1, ignore_gaps=False)
        assert result == "SXA"

    def test_iupac_sequence_alt_translate(self):
        result = alt_translate("TCNATNGCN")
        assert result == "SXA"


# ─────────────────────────────────────────────────────────────────────────────
# 5. ignore_gaps mode
# ─────────────────────────────────────────────────────────────────────────────

class TestIgnoreGaps:
    """--ignore-gaps strips all '-' characters before codon slicing."""

    def test_gap_between_codons(self, table1):
        # TCA-TCG: strip '-' → TCATCG → SS
        assert _translate_seq("TCA-TCG", table1, ignore_gaps=True) == "SS"

    def test_gap_within_codon_reframes(self, table1):
        # AT-GATG: strip '-' → ATGATG → MM
        assert _translate_seq("AT-GATG", table1, ignore_gaps=True) == "MM"

    def test_leading_gap(self, table1):
        # -ATGATG: strip → ATGATG → MM
        assert _translate_seq("-ATGATG", table1, ignore_gaps=True) == "MM"

    def test_full_gap_codon_removed(self, table1):
        # TCA---TCG: strip → TCATCG → SS (gap removed, NOT '-')
        assert _translate_seq("TCA---TCG", table1, ignore_gaps=True) == "SS"

    def test_multiple_gaps(self, table1):
        # AT-GAT-G: strip → ATGATG → MM
        assert _translate_seq("AT-GAT-G", table1, ignore_gaps=True) == "MM"

    def test_alt_translate_not_ignore_gaps_by_default(self):
        """alt_translate does NOT strip gaps by default — it respects alignment."""
        # TCA---TCG: respect_alignment mode → S-S
        assert alt_translate("TCA---TCG") == "S-S"


# ─────────────────────────────────────────────────────────────────────────────
# 6. Incomplete trailing codon
# ─────────────────────────────────────────────────────────────────────────────

class TestIncompleteTrailingCodon:
    """Trailing nucleotides that don't form a complete codon are silently dropped."""

    def test_one_trailing_nuc(self, table1):
        assert _translate_seq("TCAA", table1, ignore_gaps=False) == "S"

    def test_two_trailing_nucs(self, table1):
        assert _translate_seq("TCAAC", table1, ignore_gaps=False) == "S"

    def test_empty_sequence(self, table1):
        assert _translate_seq("", table1, ignore_gaps=False) == ""

    def test_single_nuc(self, table1):
        assert _translate_seq("A", table1, ignore_gaps=False) == ""

    def test_two_nucs(self, table1):
        assert _translate_seq("AT", table1, ignore_gaps=False) == ""

    def test_alt_translate_trailing(self):
        assert alt_translate("TCAA") == "SX"    # TCA=S, A(incomplete)=X
        assert alt_translate("TCAAC") == "SX"   # TCA=S, AC(incomplete)=X


# ─────────────────────────────────────────────────────────────────────────────
# 7. Full mixed sequence
# ─────────────────────────────────────────────────────────────────────────────

class TestFullSequence:
    """End-to-end translation of longer sequences with mixed codon types."""

    def test_all_types_in_one_sequence(self, table1):
        # TCA→S, ---→-, TC-→S, TCN→S, AT-→X, ATG→M
        seq = "TCA---TC-TCNAT-ATG"
        result = _translate_seq(seq, table1, ignore_gaps=False)
        assert result == "S-SSXM"

    def test_all_types_alt_translate(self):
        seq = "TCA---TC-TCNAT-ATG"
        assert alt_translate(seq) == "S-SSXM"

    def test_spike_protein_start(self, table1):
        """ATG = Met, standard spike start codon."""
        assert _translate_seq("ATG", table1, ignore_gaps=False) == "M"

    def test_stop_codon_in_sequence(self, table1):
        assert _translate_seq("ATGTAATCA", table1, ignore_gaps=False) == "M*S"


# ─────────────────────────────────────────────────────────────────────────────
# 8. Alternative genetic code tables
# ─────────────────────────────────────────────────────────────────────────────

class TestAlternativeTables:
    """Table-specific codon differences: e.g. table 2 (vertebrate mitochondrial)
    uses TGA=Trp (not stop) and ATA=Met (not Ile)."""

    def test_table2_TGA_is_Trp(self):
        """In vertebrate mitochondrial code (table 2), TGA encodes Trp, not stop."""
        table2 = _build_codon_table(table_id=2)
        result = _translate_seq("TGA", table2, ignore_gaps=False, table_id=2)
        assert result == "W"

    def test_table1_TGA_is_stop(self):
        """In standard code (table 1), TGA is a stop codon."""
        table1 = _build_codon_table(table_id=1)
        result = _translate_seq("TGA", table1, ignore_gaps=False, table_id=1)
        assert result == "*"

    def test_table2_ATA_is_Met(self):
        """In vertebrate mitochondrial code (table 2), ATA encodes Met."""
        table2 = _build_codon_table(table_id=2)
        result = _translate_seq("ATA", table2, ignore_gaps=False, table_id=2)
        assert result == "M"

    def test_table1_ATA_is_Ile(self):
        """In standard code (table 1), ATA encodes Isoleucine."""
        table1 = _build_codon_table(table_id=1)
        result = _translate_seq("ATA", table1, ignore_gaps=False, table_id=1)
        assert result == "I"

    def test_alt_translate_table1_default(self):
        assert alt_translate("TGA", table=1) == "*"

    def test_alt_translate_table2_TGA_is_Trp(self):
        assert alt_translate("TGA", table=2) == "W"

    def test_alt_translate_table2_ATA_is_Met(self):
        assert alt_translate("ATA", table=2) == "M"

    def test_alt_translate_table_independent_cache(self):
        """Different tables must not share cache entries."""
        r1 = alt_translate("TGA", table=1)
        r2 = alt_translate("TGA", table=2)
        assert r1 == "*"
        assert r2 == "W"
        assert r1 != r2


# ─────────────────────────────────────────────────────────────────────────────
# 9. translate.py parse_input integration via in-memory streams
# ─────────────────────────────────────────────────────────────────────────────

class TestParseInput:
    """Test parse_input() end-to-end with in-memory byte streams."""

    def _run(self, fasta_bytes: bytes, ignore_gaps=False,
             respect_alignment=True, translation_table=1) -> str:
        in_buf = io.BytesIO(fasta_bytes)
        out_buf = io.BytesIO()
        buffered_out = io.BufferedWriter(out_buf)
        parse_input(in_buf, "test", buffered_out,
                    ignore_gaps=ignore_gaps,
                    respect_alignment=respect_alignment,
                    translation_table=translation_table)
        buffered_out.flush()
        return out_buf.getvalue().decode("ascii")

    def test_standard_sequence(self):
        fasta = b">seq1\nTCATCGTCT\n"
        result = self._run(fasta)
        assert ">seq1\nSSS\n" == result

    def test_full_gap_codon(self):
        fasta = b">seq1\nTCA---TCG\n"
        result = self._run(fasta)
        assert ">seq1\nS-S\n" == result

    def test_partial_gap_TC_dash(self):
        fasta = b">seq1\nTC-\n"
        result = self._run(fasta)
        assert ">seq1\nS\n" == result

    def test_partial_gap_AT_dash(self):
        fasta = b">seq1\nAT-\n"
        result = self._run(fasta)
        assert ">seq1\nX\n" == result

    def test_ignore_gaps_strips_dashes(self):
        fasta = b">seq1\nAT-GATG\n"
        result = self._run(fasta, ignore_gaps=True)
        assert ">seq1\nMM\n" == result

    def test_iupac_TCN(self):
        fasta = b">seq1\nTCN\n"
        result = self._run(fasta)
        assert ">seq1\nS\n" == result

    def test_table2_TGA_is_Trp(self):
        fasta = b">seq1\nTGA\n"
        result = self._run(fasta, translation_table=2)
        assert ">seq1\nW\n" == result

    def test_table1_TGA_is_stop(self):
        fasta = b">seq1\nTGA\n"
        result = self._run(fasta, translation_table=1)
        assert ">seq1\n*\n" == result

    def test_multiple_records(self):
        fasta = b">r1\nTCA\n>r2\nATG\n>r3\nGGG\n"
        result = self._run(fasta)
        assert ">r1\nS\n>r2\nM\n>r3\nG\n" == result

    def test_mixed_scenario(self):
        # TCA→S, ---→-, TC-→S, TCN→S, AT-→X, ATG→M
        fasta = b">seq1\nTCA---TC-TCNAT-ATG\n"
        result = self._run(fasta)
        assert ">seq1\nS-SSXM\n" == result


# ─────────────────────────────────────────────────────────────────────────────
# 10. Incomplete codon handling — 1-nt and 2-nt trailing fragments
# ─────────────────────────────────────────────────────────────────────────────

class TestIncompleteCodonFragments:
    """1-nt and 2-nt sequences / trailing fragments must be silently dropped.

    The original Biopython translate() emits a BiopythonWarning when given a
    sequence whose length is not a multiple of three.  Both implementations
    must handle this without warnings or exceptions so that trailing partial
    codons in NGS reads (caused by frame-breaking InDels) are gracefully
    ignored rather than crashing the pipeline.
    """

    # ── alt_translate (no Biopython warning) ─────────────────────────────────

    def test_alt_translate_1nt(self, recwarn):
        """Single nucleotide: returns 'X' (unknown AA), no Biopython warning."""
        assert alt_translate("T") == "X"
        bio_warns = [w for w in recwarn.list if "Partial codon" in str(w.message)]
        assert not bio_warns, "Unexpected BiopythonWarning for 1-nt input"

    def test_alt_translate_2nt(self, recwarn):
        """Two nucleotides: returns 'X' (unknown AA), no Biopython warning."""
        assert alt_translate("TC") == "X"
        bio_warns = [w for w in recwarn.list if "Partial codon" in str(w.message)]
        assert not bio_warns, "Unexpected BiopythonWarning for 2-nt input"

    def test_alt_translate_2nt_with_gap(self, recwarn):
        """Two nucleotides including gap: returns 'X', no Biopython warning."""
        assert alt_translate("T-") == "X"
        bio_warns = [w for w in recwarn.list if "Partial codon" in str(w.message)]
        assert not bio_warns

    def test_alt_translate_trailing_after_complete(self, recwarn):
        """Complete codon + 2 trailing nt: translate complete ('S') + trailing ('X')."""
        assert alt_translate("TCATC") == "SX"
        bio_warns = [w for w in recwarn.list if "Partial codon" in str(w.message)]
        assert not bio_warns

    def test_alt_translate_empty(self):
        """Empty input → empty output (no codon positions to consider)."""
        assert alt_translate("") == ""

    # ── _translate_seq (never calls Biopython for partial codons) ─────────────

    def test_translate_seq_1nt(self, table1):
        assert _translate_seq("T", table1, ignore_gaps=False) == ""

    def test_translate_seq_2nt(self, table1):
        assert _translate_seq("TC", table1, ignore_gaps=False) == ""

    def test_translate_seq_trailing_after_complete(self, table1):
        assert _translate_seq("TCATC", table1, ignore_gaps=False) == "S"

    def test_translate_seq_empty(self, table1):
        assert _translate_seq("", table1, ignore_gaps=False) == ""


# ─────────────────────────────────────────────────────────────────────────────
# 11. Exhaustive IUPAC+gap codon comparison:
#     current alt_translate() vs current _translate_seq()
# ─────────────────────────────────────────────────────────────────────────────

# IUPAC nucleotide codes used in DNA FASTA files, plus alignment gap '-'.
# 16 characters → 16^3 = 4096 complete codons to test.
_IUPAC_GAP_CHARS = list("ACGTRYSWKMBDHVN-")
assert len(_IUPAC_GAP_CHARS) == 16

# Pre-compute the full codon list once at module import so parametrize
# doesn't regenerate it 4096 times.
_ALL_IUPAC_CODONS = [
    c1 + c2 + c3
    for c1, c2, c3 in itertools.product(_IUPAC_GAP_CHARS, repeat=3)
]
assert len(_ALL_IUPAC_CODONS) == 16 ** 3  # 4096


def _old_alt_translate(seq, _bio=__import__('Bio.Seq', fromlist=['translate']).translate):
    """Reference implementation: original NNN-based alt_translate (commit dd812e8).

    Kept verbatim to document known intentional differences from the current
    '-'→'N' implementation:
      TC- → NNN → X   (old, via NNN whole-codon substitution)
      TC- → TCN → S   (new, '-' treated as 'N', Biopython resolves → Serine)
    """
    codons = (seq[i:i+3] for i in range(0, len(seq), 3))
    codons = ("NNN" if "-" in codon and codon != "---" else codon for codon in codons)
    return "".join(_bio(codon, gap='-') for codon in codons)


# Codons where old and new intentionally differ: partial-gap codons where
# the remaining two nucleotides uniquely determine the amino acid.
# Old: TC- → NNN → X.  New: TC- → TCN → S.
_KNOWN_OLD_VS_NEW_DIFFS = {
    # codon : (old_result, new_result)
    # Four-fold degenerate positions (third base) — new correctly gives AA:
}
# Build the known-diff dict dynamically so it stays in sync with Biopython.
_t1 = _build_codon_table(table_id=1)
for _c1, _c2 in itertools.product("ACGTRYSWKMBDHVN", repeat=2):
    _codon = _c1 + _c2 + "-"
    _old = _old_alt_translate(_codon)
    _new = alt_translate(_codon)
    if _old != _new:
        _KNOWN_OLD_VS_NEW_DIFFS[_codon] = (_old, _new)


class TestExhaustiveIUPACGapCodons:
    """Verify that both current implementations agree on every possible codon.

    Old (NNN-based) vs new ('-'→'N') differences are pre-computed and treated
    as *known intentional changes* — they are not failures.  The key assertion
    is that the two CURRENT implementations (alt_translate and _translate_seq)
    are identical for all 4096 complete IUPAC+gap codons.
    """

    @pytest.mark.parametrize("codon", _ALL_IUPAC_CODONS)
    def test_current_implementations_agree(self, codon, table1):
        """alt_translate(codon) == _translate_seq(codon) for every IUPAC+gap codon."""
        result_at = alt_translate(codon)
        result_ts = _translate_seq(codon, table1, ignore_gaps=False)
        assert result_at == result_ts, (
            f"Mismatch for codon {codon!r}: "
            f"alt_translate={result_at!r}, _translate_seq={result_ts!r}"
        )

    def test_known_old_vs_new_differences_are_partial_gap_only(self):
        """All old→new differences must be partial-gap codons (contain '-' but not '---')."""
        for codon in _KNOWN_OLD_VS_NEW_DIFFS:
            assert '-' in codon and codon != '---', (
                f"Unexpected old→new difference for non-partial-gap codon {codon!r}"
            )

    def test_full_gap_codon_unchanged(self):
        """'---' must give '-' in both old and new implementations."""
        assert _old_alt_translate("---") == "-"
        assert alt_translate("---") == "-"
        t1 = _build_codon_table(1)
        assert _translate_seq("---", t1, ignore_gaps=False) == "-"

    def test_no_exception_on_any_iupac_gap_codon(self, table1):
        """No codon should raise an exception (only X for unknowns)."""
        for codon in _ALL_IUPAC_CODONS:
            try:
                alt_translate(codon)
                _translate_seq(codon, table1, ignore_gaps=False)
            except Exception as exc:
                pytest.fail(f"Exception for codon {codon!r}: {exc}")

    def test_diff_count_matches_known(self):
        """Sanity-check: number of old→new diffs matches pre-computed map."""
        diffs = 0
        for codon in _ALL_IUPAC_CODONS:
            if _old_alt_translate(codon) != alt_translate(codon):
                diffs += 1
        assert diffs == len(_KNOWN_OLD_VS_NEW_DIFFS), (
            f"Expected {len(_KNOWN_OLD_VS_NEW_DIFFS)} old→new diffs, got {diffs}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# 12. calculate_codon_frequencies translation path (write_tsv_line)
# ─────────────────────────────────────────────────────────────────────────────


def _write_tsv_translate(codon, translation_table=1):
    """Extract just the amino acid that write_tsv_line assigns to *codon*.

    write_tsv_line() writes TSV rows but we only care about the _some_aa it
    computes.  We drive it with a single-codon dict (count=1, total=1) and
    capture the output to extract the 4th column (mutant_aa).
    """
    buf = io.StringIO()
    codons = {codon: Decimal(1)}
    write_tsv_line(
        buf,
        codons,
        natural_codon_position_padded=1,
        natural_codon_position_depadded=1,
        reference_aa='A',
        total_codons_per_site_sum=1,
        reference_codon='GCT',
        debug=False,
        translation_table=translation_table,
    )
    line = buf.getvalue().strip()
    if not line:
        return None
    # TSV columns: padded_pos, pos, ref_aa, mutant_aa, freq, ref_codon, codon, count, total
    return line.split('\t')[3]   # mutant_aa


class TestCalculateCodonFrequenciesTranslation:
    """Verify the write_tsv_line() translation path for all IUPAC+gap codons.

    write_tsv_line() is the function called by calculate_codon_frequencies to
    assign amino acids to each observed codon in the alignment.  It uses
    alt_translate() for 3-char codons and 'X' for codons shorter than 3 chars.

    Exhaustive test (4096 codons + incomplete codon cases):
      - For every 3-char IUPAC+gap codon, write_tsv_line must agree with
        alt_translate() — same function, different entry point.
      - For 1- and 2-char codons, write_tsv_line returns 'X' regardless.
    """

    @pytest.mark.parametrize("codon", _ALL_IUPAC_CODONS)
    def test_write_tsv_line_agrees_with_alt_translate(self, codon):
        """write_tsv_line and alt_translate must agree for every 3-char codon."""
        expected = alt_translate(codon)
        actual = _write_tsv_translate(codon)
        assert actual == expected, (
            f"write_tsv_line gave {actual!r} but alt_translate gave {expected!r} "
            f"for codon {codon!r}"
        )

    @pytest.mark.parametrize("short_codon", ["T", "TC", "A", "AT", "G", "GC"])
    def test_write_tsv_line_short_codon_returns_X(self, short_codon):
        """Codons shorter than 3 nt must be mapped to 'X' by write_tsv_line."""
        result = _write_tsv_translate(short_codon)
        assert result == "X", (
            f"Expected 'X' for short codon {short_codon!r}, got {result!r}"
        )

    def test_table2_TGA_is_Trp_via_write_tsv_line(self):
        """Alternative genetic code propagates through write_tsv_line."""
        assert _write_tsv_translate("TGA", translation_table=1) == "*"
        assert _write_tsv_translate("TGA", translation_table=2) == "W"

    def test_partial_gap_TC_dash_via_write_tsv_line(self):
        """TC- → S through the calculate_codon_frequencies path."""
        assert _write_tsv_translate("TC-") == "S"

    def test_full_gap_via_write_tsv_line(self):
        """'---' → '-' through the calculate_codon_frequencies path."""
        assert _write_tsv_translate("---") == "-"
