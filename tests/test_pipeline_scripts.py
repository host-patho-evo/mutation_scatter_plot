#!/usr/bin/env python3
"""Tests for the FASTA deduplication pipeline helper scripts.

Covers:
- Per-line UTF-8 → Latin-1 fallback decoding (_decode_fasta_line)
- Timestamp-aware output guards (skip when up-to-date, error when stale)
- --overwrite flag bypasses staleness check
- Accidental '=' detection in --outfile-prefix / --mapping-outfile
- Parent-child naming convention inference in summarize_fasta_pipeline
- Cache freshness helpers: _fresh_discarded_txt, _fresh_tsv
"""

import hashlib
import importlib.util
import os
import subprocess
import sys
import tempfile
import time
import unittest

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'scripts')


def _script(name: str) -> str:
    return os.path.join(SCRIPTS_DIR, name)


def _load_module(script_name: str):
    """Import a scripts/ file as a module (only works for side-effect-free modules)."""
    spec = importlib.util.spec_from_file_location(
        script_name.replace('.py', ''),
        _script(script_name),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _make_nnnx_fasta(path: str, entries: list) -> None:
    """Write a minimal NNNNx.sha256 FASTA to path.

    entries: list of (count: int, seq: str) tuples.
    """
    with open(path, 'w', encoding='utf-8') as f:
        for count, seq in entries:
            h = hashlib.sha256(seq.upper().encode()).hexdigest()
            f.write(f">{count}x.{h}\n{seq}\n")


def _make_plain_fasta(path: str, records: list) -> None:
    """Write a plain FASTA (GISAID-style IDs) to path.

    records: list of (id: str, seq: str) tuples.
    """
    with open(path, 'w', encoding='utf-8') as f:
        for rid, seq in records:
            f.write(f">{rid}\n{seq}\n")


# ── replicate _decode_fasta_line locally ──────────────────────────────────────
# This avoids importing the full script (which calls parse_args at top level).

def _decode_fasta_line(raw: bytes) -> str:
    """Per-line UTF-8 → Latin-1 fallback decoder (mirrors scripts implementation)."""
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Per-line encoding decoder
# ═══════════════════════════════════════════════════════════════════════════════

class TestDecodeFastaLine(unittest.TestCase):
    """Unit tests for the per-line UTF-8 → Latin-1 fallback decoder."""

    def test_pure_ascii_header(self):
        raw = b">EPI_ISL_123456 hCoV-19/Poland/Warsaw/2022\n"
        self.assertEqual(_decode_fasta_line(raw), raw.decode("ascii"))

    def test_pure_ascii_sequence(self):
        raw = b"ATGCATGCATGCNNNN\n"
        self.assertEqual(_decode_fasta_line(raw), raw.decode("ascii"))

    def test_valid_utf8_multibyte(self):
        # é encoded as UTF-8 (0xC3 0xA9) — should decode correctly.
        self.assertEqual(_decode_fasta_line("Łódź".encode("utf-8")), "Łódź")

    def test_latin1_fallback_0xed(self):
        # 0xED is not valid as a UTF-8 continuation byte; in Latin-1 it is 'í'.
        raw = b"Szpital Specjalistyczny im. Edmund\xed Biernackiego"
        result = _decode_fasta_line(raw)
        self.assertIsInstance(result, str)
        self.assertIn("í", result)  # 0xED Latin-1 = í

    def test_latin1_fallback_0xe9(self):
        # 0xE9 alone is invalid UTF-8 start; Latin-1 maps it to 'é'.
        raw = b"H\xf4pital universit\xe9"
        result = _decode_fasta_line(raw)
        self.assertIsInstance(result, str)
        self.assertIn("ô", result)
        self.assertIn("é", result)

    def test_id_part_is_ascii_despite_description_encoding(self):
        # The ID (first word) is always ASCII even if description has non-UTF-8 bytes.
        raw = b">EPI_ISL_9876543 Institut P\xe9diatrique de Lyon\n"
        result = _decode_fasta_line(raw)
        self.assertTrue(result.startswith(">EPI_ISL_9876543"))

    def test_nnnx_id_unaffected(self):
        raw = b">914x.abc123 inst. \xedme Sp\xf3lki\n"
        result = _decode_fasta_line(raw)
        self.assertTrue(result.startswith(">914x.abc123"))


# ═══════════════════════════════════════════════════════════════════════════════
# 2. summarize_fasta_pipeline – naming convention & cache freshness
# ═══════════════════════════════════════════════════════════════════════════════

class TestSummarizePipelineHelpers(unittest.TestCase):
    """Test helpers in summarize_fasta_pipeline (side-effect-free to import)."""

    @classmethod
    def setUpClass(cls):
        cls.mod = _load_module("summarize_fasta_pipeline.py")

    # ── _strip_fasta_suffix ───────────────────────────────────────────────────

    def test_strip_fasta(self):
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.counts.fasta"), "prefix.counts")

    def test_strip_fasta_old(self):
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.counts.fasta.old"), "prefix.counts")

    def test_strip_fasta_orig(self):
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.fasta.orig"), "prefix")

    def test_strip_fasta_ori(self):
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.fasta.ori"), "prefix")

    def test_no_suffix_unchanged(self):
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.tsv"), "prefix.tsv")

    def test_longest_suffix_first(self):
        # .fasta.orig must win over .fasta to avoid leaving '.orig'.
        self.assertEqual(self.mod._strip_fasta_suffix("a.fasta.orig"), "a")

    # ── _fresh_discarded_txt ──────────────────────────────────────────────────

    def test_fresh_discarded_txt_absent(self):
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            child  = os.path.join(d, "a.counts.fasta")
            for p in (parent, child):
                open(p, 'w').close()
            self.assertIsNone(self.mod._fresh_discarded_txt(child, parent))

    def test_fresh_discarded_txt_stale(self):
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            child  = os.path.join(d, "a.counts.fasta")
            txt    = os.path.join(d, "a.counts.discarded_original_ids.txt")
            for p in (txt, parent, child):
                open(p, 'w').close()
            time.sleep(0.05)
            os.utime(child, None)   # child is newest → txt is stale
            self.assertIsNone(self.mod._fresh_discarded_txt(child, parent))

    def test_fresh_discarded_txt_up_to_date(self):
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            child  = os.path.join(d, "a.counts.fasta")
            txt    = os.path.join(d, "a.counts.discarded_original_ids.txt")
            for p in (parent, child):
                open(p, 'w').close()
            time.sleep(0.05)
            open(txt, 'w').close()  # txt is newest
            self.assertEqual(self.mod._fresh_discarded_txt(child, parent), txt)

    # ── _fresh_tsv ────────────────────────────────────────────────────────────

    def test_fresh_tsv_absent(self):
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            open(parent, 'w').close()
            self.assertIsNone(self.mod._fresh_tsv(parent))

    def test_fresh_tsv_stale(self):
        with tempfile.TemporaryDirectory() as d:
            tsv    = os.path.join(d, "a.sha256_to_ids.tsv")
            parent = os.path.join(d, "a.fasta")
            for p in (tsv, parent):
                open(p, 'w').close()
            time.sleep(0.05)
            os.utime(parent, None)  # parent is newest → tsv stale
            self.assertIsNone(self.mod._fresh_tsv(parent))

    def test_fresh_tsv_up_to_date(self):
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            tsv    = os.path.join(d, "a.sha256_to_ids.tsv")
            open(parent, 'w').close()
            time.sleep(0.05)
            open(tsv, 'w').close()  # tsv is newest
            self.assertEqual(self.mod._fresh_tsv(parent), tsv)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. count_same_sequences.py CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestCountSameSequencesCLI(unittest.TestCase):
    """CLI smoke tests for count_same_sequences.py."""

    SCRIPT = _script("count_same_sequences.py")

    def _run(self, *args):
        return subprocess.run(
            [sys.executable, self.SCRIPT] + list(args),
            capture_output=True, text=True,
        )

    # ── '=' detection ─────────────────────────────────────────────────────────

    def test_equals_in_outfile_prefix_rejected(self):
        """Typo --outfile-prefix=infilename=foo is caught immediately."""
        result = self._run(
            "--infilename=/dev/null",
            "--outfile-prefix=infilename=something",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("contains '='", result.stderr + result.stdout)

    def test_equals_in_mapping_outfile_rejected(self):
        result = self._run(
            "--infilename=/dev/null",
            "--mapping-outfile=mapping-outfile=foo.tsv",
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("contains '='", result.stderr + result.stdout)

    # ── timestamp guard ───────────────────────────────────────────────────────

    def test_uptodate_exits_zero(self):
        """Outputs newer than input → prints 'up-to-date' and exits 0."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            prefix = os.path.join(d, "in")
            tsv    = os.path.join(d, "in.sha256_to_ids.tsv")
            counts = os.path.join(d, "in.counts.fasta")
            with open(infile, 'w') as f:
                f.write(">ID1\nATGC\n")
            time.sleep(0.05)
            for p in (counts, tsv):  # outputs are newer
                open(p, 'w').close()
            result = self._run(
                f"--infilename={infile}",
                f"--outfile-prefix={prefix}",
                f"--mapping-outfile={tsv}",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("up-to-date", result.stderr)

    def test_stale_output_errors_without_overwrite(self):
        """Output older than input without --overwrite → non-zero exit."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            prefix = os.path.join(d, "in")
            tsv    = os.path.join(d, "in.sha256_to_ids.tsv")
            counts = os.path.join(d, "in.counts.fasta")
            for p in (counts, tsv):  # outputs created first → will be stale
                open(p, 'w').close()
            time.sleep(0.05)
            with open(infile, 'w') as f:  # input is newest
                f.write(">ID1\nATGC\n")
            result = self._run(
                f"--infilename={infile}",
                f"--outfile-prefix={prefix}",
                f"--mapping-outfile={tsv}",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("stale", result.stderr + result.stdout)

    def test_existing_output_without_overwrite_errors(self):
        """Pre-existing output (no staleness info yet different from stale) → error."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            prefix = os.path.join(d, "in")
            tsv    = os.path.join(d, "in.sha256_to_ids.tsv")
            with open(infile, 'w') as f:
                f.write(">ID1\nATGC\n")
            time.sleep(0.05)
            counts = os.path.join(d, "in.counts.fasta")
            open(counts, 'w').close()   # only one output exists → partial run
            result = self._run(
                f"--infilename={infile}",
                f"--outfile-prefix={prefix}",
                f"--mapping-outfile={tsv}",
            )
            self.assertNotEqual(result.returncode, 0)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. create_list_of_discarded_sequences.py CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestCreateListCLI(unittest.TestCase):
    """CLI smoke tests for create_list_of_discarded_sequences.py."""

    SCRIPT = _script("create_list_of_discarded_sequences.py")

    def _run(self, *args):
        return subprocess.run(
            [sys.executable, self.SCRIPT] + list(args),
            capture_output=True, text=True,
        )

    def test_uptodate_exits_zero(self):
        """Output newer than infilename → skip with exit 0."""
        with tempfile.TemporaryDirectory() as d:
            infile  = os.path.join(d, "kept.fasta")
            outfile = os.path.join(d, "kept.discarded_original_ids.txt")
            _make_nnnx_fasta(infile, [(1, "ATGC")])
            time.sleep(0.05)
            open(outfile, 'w').close()  # outfile is newer
            result = self._run(
                f"--infilename={infile}",
                f"--outfile={outfile}",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("up-to-date", result.stderr)

    def test_stale_errors_without_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            infile  = os.path.join(d, "kept.fasta")
            outfile = os.path.join(d, "kept.discarded_original_ids.txt")
            open(outfile, 'w').close()
            time.sleep(0.05)
            _make_nnnx_fasta(infile, [(1, "ATGC")])  # input is newer
            result = self._run(
                f"--infilename={infile}",
                f"--outfile={outfile}",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("stale", result.stderr + result.stdout)

    def test_overwrite_bypasses_staleness(self):
        """--overwrite causes stale check to be skipped."""
        with tempfile.TemporaryDirectory() as d:
            infile  = os.path.join(d, "kept.fasta")
            outfile = os.path.join(d, "kept.discarded_original_ids.txt")
            open(outfile, 'w').close()
            time.sleep(0.05)
            _make_nnnx_fasta(infile, [(1, "ATGC")])
            result = self._run(
                f"--infilename={infile}",
                f"--outfile={outfile}",
                "--overwrite",
            )
            # No staleness error; any other error is unrelated to timestamps.
            combined = result.stderr + result.stdout
            self.assertNotIn("stale", combined)

    def test_devnull_outfile_skips_guard(self):
        """/dev/null as --outfile bypasses the timestamp guard entirely."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "kept.fasta")
            _make_nnnx_fasta(infile, [(1, "ATGC")])
            result = self._run(
                f"--infilename={infile}",
                "--outfile=/dev/null",
            )
            # Should not fail with a staleness error.
            combined = result.stderr + result.stdout
            self.assertNotIn("stale", combined)
            self.assertNotIn("up-to-date", combined)

    def test_mixed_encoding_fasta_does_not_crash(self):
        """FASTA headers with non-UTF-8 bytes (Latin-1) are read without error."""
        with tempfile.TemporaryDirectory() as d:
            infile  = os.path.join(d, "mixed.fasta")
            outfile = os.path.join(d, "mixed.discarded_original_ids.txt")
            # Write a FASTA with a Latin-1 byte (0xED = í) in a description.
            with open(infile, 'wb') as f:
                seq = b"ATGCATGC"
                h = hashlib.sha256(seq).hexdigest()
                f.write(b">1x." + h.encode() + b" Institut Sp\xedtalu\n")
                f.write(seq + b"\n")
            result = self._run(
                f"--infilename={infile}",
                f"--outfile={outfile}",
            )
            # The script should not crash with a UnicodeDecodeError.
            self.assertNotIn("UnicodeDecodeError", result.stderr)
            self.assertNotIn("UnicodeDecodeError", result.stdout)


# ═══════════════════════════════════════════════════════════════════════════════
# 5. kick.py CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestKickCLI(unittest.TestCase):
    """CLI smoke tests for kick.py."""

    SCRIPT = _script("kick.py")

    def _run(self, *args):
        return subprocess.run(
            [sys.executable, self.SCRIPT] + list(args),
            capture_output=True, text=True,
        )

    def test_missing_full_length_errors(self):
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "test.fasta")
            _make_plain_fasta(infile, [("1x.abc", "ATGCATGC")])
            result = self._run(f"--infile={infile}", "--outfile-prefix=test")
            self.assertNotEqual(result.returncode, 0)

    def test_uptodate_exits_zero(self):
        """All outputs newer than input → skip with exit 0."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "test.fasta")
            prefix = os.path.join(d, "test")
            _make_plain_fasta(infile, [("1x.abc", "ATGCATGC")])
            time.sleep(0.05)
            for sfx in ("exactly_4.fasta", "shorter_4.fasta", "longer_4.fasta"):
                open(os.path.join(d, f"test.{sfx}"), 'w').close()
            result = self._run(
                f"--infile={infile}",
                f"--outfile-prefix={prefix}",
                "--full-length=4",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("up-to-date", result.stderr)

    def test_stale_errors_without_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "test.fasta")
            prefix = os.path.join(d, "test")
            for sfx in ("exactly_4.fasta", "shorter_4.fasta", "longer_4.fasta"):
                open(os.path.join(d, f"test.{sfx}"), 'w').close()
            time.sleep(0.05)
            _make_plain_fasta(infile, [("1x.abc", "ATGCATGC")])  # input newer
            result = self._run(
                f"--infile={infile}",
                f"--outfile-prefix={prefix}",
                "--full-length=4",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("stale", result.stderr + result.stdout)

    def test_overwrite_bypasses_staleness(self):
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "test.fasta")
            prefix = os.path.join(d, "test")
            for sfx in ("exactly_4.fasta", "shorter_4.fasta", "longer_4.fasta"):
                open(os.path.join(d, f"test.{sfx}"), 'w').close()
            time.sleep(0.05)
            _make_plain_fasta(infile, [("seq1", "ATGC"), ("seq2", "ATGCATGCATGC")])
            result = self._run(
                f"--infile={infile}",
                f"--outfile-prefix={prefix}",
                "--full-length=4",
                "--overwrite",
            )
            combined = result.stderr + result.stdout
            self.assertNotIn("stale", combined)


if __name__ == '__main__':
    unittest.main()
