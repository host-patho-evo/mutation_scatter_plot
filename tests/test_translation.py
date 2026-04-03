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

import io
import sys
import os
import pytest

# ── helpers / imports ─────────────────────────────────────────────────────────

# Import from the package
from mutation_scatter_plot import alt_translate

# Import internal functions from the scripts/ translate tool.
# Since scripts/ is not installed as a package, add it to sys.path.
# We also temporarily patch sys.argv to an empty list so that translate.py's
# module-level  (myoptions, myargs) = myparser.parse_args()  call does not
# see pytest's own arguments and raise an OptionParser error.
_SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'scripts')
sys.path.insert(0, os.path.abspath(_SCRIPTS_DIR))

_saved_argv = sys.argv
sys.argv = ['translate']
try:
    from translate import _build_codon_table, _translate_seq, parse_input  # noqa: E402
finally:
    sys.argv = _saved_argv


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
        assert alt_translate("TCAA") == "S"
        assert alt_translate("TCAAC") == "S"


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
        import io as _io
        buffered_out = _io.BufferedWriter(out_buf)
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
