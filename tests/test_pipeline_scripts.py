#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
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
SRC_SCRIPTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'src', 'mutation_scatter_plot', 'scripts',
)


def _script(name: str) -> str:
    return os.path.join(SCRIPTS_DIR, name)


def _src_script(name: str) -> str:
    """Return absolute path to a script inside src/mutation_scatter_plot/scripts/."""
    return os.path.join(SRC_SCRIPTS_DIR, name)


def _load_module(script_name: str):
    """Import a scripts/ file as a module (only works for side-effect-free modules)."""
    spec = importlib.util.spec_from_file_location(
        script_name.replace('.py', ''),
        _script(script_name),
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _touch(path: str) -> None:
    """Create an empty file (or update its mtime), safely using a context manager."""
    with open(path, 'w', encoding='utf-8'):
        pass


def _make_nnnx_fasta(path: str, entries: list) -> None:
    """Write a minimal NNNNx.sha256 FASTA to path.

    entries: list of (count: int, seq: str) tuples.
    """
    with open(path, 'w', encoding='utf-8') as f:
        for count, seq in entries:
            h = hashlib.sha256(seq.upper().encode()).hexdigest()
            f.write(f'>{count}x.{h}\n{seq}\n')


def _make_plain_fasta(path: str, records: list) -> None:
    """Write a plain FASTA (GISAID-style IDs) to path.

    records: list of (id: str, seq: str) tuples.
    """
    with open(path, 'w', encoding='utf-8') as f:
        for rid, seq in records:
            f.write(f'>{rid}\n{seq}\n')


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
        """Pure ASCII FASTA header decodes identically via UTF-8 and ASCII."""
        raw = b">EPI_ISL_123456 hCoV-19/Poland/Warsaw/2022\n"
        self.assertEqual(_decode_fasta_line(raw), raw.decode("ascii"))

    def test_pure_ascii_sequence(self):
        """Pure ASCII sequence line decodes identically."""
        raw = b"ATGCATGCATGCNNNN\n"
        self.assertEqual(_decode_fasta_line(raw), raw.decode("ascii"))

    def test_valid_utf8_multibyte(self):
        """Valid multi-byte UTF-8 sequence decodes correctly."""
        # é encoded as UTF-8 (0xC3 0xA9) — should decode correctly.
        self.assertEqual(_decode_fasta_line("Łódź".encode("utf-8")), "Łódź")

    def test_latin1_fallback_0xed(self):
        """0xED is not valid UTF-8 continuation; Latin-1 fallback maps it to í."""
        raw = b"Szpital Specjalistyczny im. Edmund\xed Biernackiego"
        result = _decode_fasta_line(raw)
        self.assertIsInstance(result, str)
        self.assertIn("í", result)  # 0xED Latin-1 = í

    def test_latin1_fallback_0xe9(self):
        """0xE9 alone is invalid UTF-8 start; Latin-1 maps it to é."""
        raw = b"H\xf4pital universit\xe9"
        result = _decode_fasta_line(raw)
        self.assertIsInstance(result, str)
        self.assertIn("ô", result)
        self.assertIn("é", result)

    def test_id_part_is_ascii_despite_description_encoding(self):
        """The ID (first word) is always ASCII even if description has non-UTF-8 bytes."""
        raw = b">EPI_ISL_9876543 Institut P\xe9diatrique de Lyon\n"
        result = _decode_fasta_line(raw)
        self.assertTrue(result.startswith(">EPI_ISL_9876543"))

    def test_nnnx_id_unaffected(self):
        """NNNNx.hash ID prefix is preserved intact regardless of description encoding."""
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
        """Load the summarize_fasta_pipeline module once for all tests."""
        cls.mod = _load_module("summarize_fasta_pipeline.py")

    # ── _strip_fasta_suffix ───────────────────────────────────────────────────

    def test_strip_fasta(self):
        """Standard .fasta suffix is stripped."""
        self.assertEqual(
            self.mod._strip_fasta_suffix("prefix.counts.fasta"), "prefix.counts"  # pylint: disable=protected-access
        )

    def test_strip_fasta_old(self):
        """.fasta.old suffix is stripped leaving the base."""
        self.assertEqual(
            self.mod._strip_fasta_suffix("prefix.counts.fasta.old"), "prefix.counts"  # pylint: disable=protected-access
        )

    def test_strip_fasta_orig(self):
        """.fasta.orig suffix is stripped."""
        self.assertEqual(
            self.mod._strip_fasta_suffix("prefix.fasta.orig"), "prefix"  # pylint: disable=protected-access
        )

    def test_strip_fasta_ori(self):
        """.fasta.ori suffix is stripped."""
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.fasta.ori"), "prefix")  # pylint: disable=protected-access

    def test_no_suffix_unchanged(self):
        """A path with no recognised FASTA suffix is returned unchanged."""
        self.assertEqual(self.mod._strip_fasta_suffix("prefix.tsv"), "prefix.tsv")  # pylint: disable=protected-access

    def test_longest_suffix_first(self):
        """.fasta.orig must match before .fasta to avoid leaving '.orig'."""
        self.assertEqual(self.mod._strip_fasta_suffix("a.fasta.orig"), "a")  # pylint: disable=protected-access

    # ── _fresh_discarded_txt ──────────────────────────────────────────────────

    def test_fresh_discarded_txt_absent(self):
        """Returns None when the .discarded_original_ids.txt file does not exist."""
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            child = os.path.join(d, "a.counts.fasta")
            for p in (parent, child):
                _touch(p)
            self.assertIsNone(self.mod._fresh_discarded_txt(child, parent))  # pylint: disable=protected-access

    def test_fresh_discarded_txt_stale(self):
        """Returns None when the txt file is older than the child FASTA."""
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            child = os.path.join(d, "a.counts.fasta")
            txt = os.path.join(d, "a.counts.discarded_original_ids.txt")
            for p in (txt, parent, child):
                _touch(p)
            time.sleep(0.05)
            os.utime(child, None)   # child is newest → txt is stale
            self.assertIsNone(self.mod._fresh_discarded_txt(child, parent))  # pylint: disable=protected-access

    def test_fresh_discarded_txt_up_to_date(self):
        """Returns the txt path when it is newer than both FASTA files."""
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            child = os.path.join(d, "a.counts.fasta")
            txt = os.path.join(d, "a.counts.discarded_sha256_hashes.txt")
            for p in (parent, child):
                _touch(p)
            time.sleep(0.05)
            _touch(txt)  # txt is newest
            self.assertEqual(self.mod._fresh_discarded_txt(child, parent), txt)  # pylint: disable=protected-access

    # ── _fresh_tsv ────────────────────────────────────────────────────────────

    def test_fresh_tsv_absent(self):
        """Returns None when the .sha256_to_ids.tsv does not exist."""
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            _touch(parent)
            self.assertIsNone(self.mod._fresh_tsv(parent))  # pylint: disable=protected-access

    def test_fresh_tsv_stale(self):
        """Returns None when the TSV is older than the parent FASTA."""
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, "a.sha256_to_ids.tsv")
            parent = os.path.join(d, "a.fasta")
            for p in (tsv, parent):
                _touch(p)
            time.sleep(0.05)
            os.utime(parent, None)  # parent is newest → tsv stale
            self.assertIsNone(self.mod._fresh_tsv(parent))  # pylint: disable=protected-access

    def test_fresh_tsv_up_to_date(self):
        """Returns the TSV path when it is newer than the parent FASTA."""
        with tempfile.TemporaryDirectory() as d:
            parent = os.path.join(d, "a.fasta")
            tsv = os.path.join(d, "a.sha256_to_ids.tsv")
            _touch(parent)
            time.sleep(0.05)
            _touch(tsv)  # tsv is newest
            self.assertEqual(self.mod._fresh_tsv(parent), tsv)  # pylint: disable=protected-access


# ═══════════════════════════════════════════════════════════════════════════════
# 3. count_same_sequences.py CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestCountSameSequencesCLI(unittest.TestCase):
    """CLI smoke tests for count_same_sequences.py."""

    SCRIPT = _script("count_same_sequences.py")

    def _run(self, *args):
        return subprocess.run(
            [sys.executable, self.SCRIPT] + list(args),
            capture_output=True, text=True, check=False,
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
        """Typo --mapping-outfile=mapping-outfile=foo is caught immediately."""
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
            tsv = os.path.join(d, "in.sha256_to_ids.tsv")
            counts = os.path.join(d, "in.counts.fasta")
            with open(infile, 'w', encoding='utf-8') as f:
                f.write(">ID1\nATGC\n")
            time.sleep(0.05)
            for p in (counts, tsv):  # outputs are newer
                _touch(p)
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
            tsv = os.path.join(d, "in.sha256_to_ids.tsv")
            counts = os.path.join(d, "in.counts.fasta")
            for p in (counts, tsv):  # outputs created first → will be stale
                _touch(p)
            time.sleep(0.05)
            with open(infile, 'w', encoding='utf-8') as f:  # input is newest
                f.write(">ID1\nATGC\n")
            result = self._run(
                f"--infilename={infile}",
                f"--outfile-prefix={prefix}",
                f"--mapping-outfile={tsv}",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("stale", result.stderr + result.stdout)

    def test_existing_output_without_overwrite_errors(self):
        """Pre-existing output (partial run) → error without --overwrite."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            prefix = os.path.join(d, "in")
            tsv = os.path.join(d, "in.sha256_to_ids.tsv")
            with open(infile, 'w', encoding='utf-8') as f:
                f.write(">ID1\nATGC\n")
            time.sleep(0.05)
            counts = os.path.join(d, "in.counts.fasta")
            _touch(counts)   # only one output exists → partial run
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
            capture_output=True, text=True, check=False,
        )

    def test_uptodate_exits_zero(self):
        """Output newer than infilename → skip with exit 0."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "kept.fasta")
            outfile = os.path.join(d, "kept.discarded_original_ids.txt")
            _make_nnnx_fasta(infile, [(1, "ATGC")])
            time.sleep(0.05)
            _touch(outfile)  # outfile is newer
            result = self._run(
                f"--infilename={infile}",
                f"--outfile={outfile}",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("up-to-date", result.stderr)

    def test_stale_errors_without_overwrite(self):
        """Output older than input without --overwrite → non-zero exit."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "kept.fasta")
            outfile = os.path.join(d, "kept.discarded_original_ids.txt")
            _touch(outfile)
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
            infile = os.path.join(d, "kept.fasta")
            outfile = os.path.join(d, "kept.discarded_original_ids.txt")
            _touch(outfile)
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
            infile = os.path.join(d, "mixed.fasta")
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
# 5. split_fasta_entries_by_lengths CLI
# ═══════════════════════════════════════════════════════════════════════════════

class TestSplitFastaByLengthsCLI(unittest.TestCase):
    """CLI smoke tests for split_fasta_entries_by_lengths (formerly kick.py)."""

    SCRIPT = _src_script("split_fasta_entries_by_lengths.py")

    def _run(self, *args):
        return subprocess.run(
            [sys.executable, self.SCRIPT] + list(args),
            capture_output=True, text=True, check=False,
        )

    def test_missing_full_length_errors(self):
        """Missing --full-length argument → non-zero exit."""
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
                _touch(os.path.join(d, f"test.{sfx}"))
            result = self._run(
                f"--infile={infile}",
                f"--outfile-prefix={prefix}",
                "--full-length=4",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("up-to-date", result.stderr)

    def test_stale_errors_without_overwrite(self):
        """Stale output (older than input) without --overwrite → non-zero exit."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "test.fasta")
            prefix = os.path.join(d, "test")
            for sfx in ("exactly_4.fasta", "shorter_4.fasta", "longer_4.fasta"):
                _touch(os.path.join(d, f"test.{sfx}"))
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
        """--overwrite skips the staleness check and proceeds normally."""
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "test.fasta")
            prefix = os.path.join(d, "test")
            for sfx in ("exactly_4.fasta", "shorter_4.fasta", "longer_4.fasta"):
                _touch(os.path.join(d, f"test.{sfx}"))
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

    def test_fasta_typo_collision_backslash_t(self):
        """Test that _extract_subset_to_fasta correctly extracts headers with physical
        backslash-t typos by perfectly masking the unescaping phase.
        """
        # Inject test inside existing TestSummarizePipelineHelpers
        with tempfile.TemporaryDirectory() as tmpdir:
            subset_fasta = os.path.join(tmpdir, "subset.fasta")
            root_fasta = os.path.join(tmpdir, "root.fasta")
            root_tsv = os.path.join(tmpdir, "root.sha256_to_descr_lines.tsv")
            target_out = os.path.join(tmpdir, "subset.discarded_original_entries.fasta")
            hashes_txt = os.path.join(tmpdir, "subset.discarded_sha256_hashes.txt")

            # GISAID sequence physically containing literal '\' and 't' characters
            test_header = "Spike|hCoV-19/USA/TX-123/2021|INCARNATE WORD CONVENT - \\t3400 BRADFORD STREET"

            with open(root_fasta, 'wb') as f:
                f.write(f">{test_header}\nATGC\n".encode('utf-8'))

            # Simulated _auto_generate_tsv escape (it escapes nothing since there's no actual 0x09 tab)
            test_sha = "abcd1234abcd1234abcd1234abcd1234abcd1234abcd1234abcd1234abcd1234"
            with open(root_tsv, 'w', encoding='utf-8') as f:
                f.write(f"{test_sha}\t1\t{test_header}\n")

            with open(hashes_txt, 'w', encoding='utf-8') as f:
                f.write(f"{test_sha}\n")

            # Run extraction using dynamically loaded module from class
            # Since TestSummarizePipelineHelpers doesn't exist at the end of the file, we just load it manually
            # if we are sitting in TestSplitFastaByLengthsCLI or general module.
            mod = _load_module("summarize_fasta_pipeline.py")
            mod._extract_subset_to_fasta(  # pylint: disable=protected-access
                subset_fasta,
                root_fasta,
                root_fasta,  # passes to mapping_path to generate root.sha256_to_descr_lines.tsv internally
                tmpdir,
                "discarded"
            )

            # Must have successfully extracted the literal typo header instead of bypassing it
            self.assertTrue(os.path.exists(target_out))
            with open(target_out, 'r', encoding='utf-8') as f:
                extracted_content = f.read()
            self.assertIn("INCARNATE WORD CONVENT - \\t3400", extracted_content)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. _collect_sha256_set – dict return type & NNNNx counts
# ═══════════════════════════════════════════════════════════════════════════════

class TestCollectSha256Set(unittest.TestCase):
    """Tests for _collect_sha256_set returning dict[str, int] with NNNNx counts."""

    @classmethod
    def setUpClass(cls):
        cls.mod = _load_module("summarize_fasta_pipeline.py")

    # ── return type ───────────────────────────────────────────────────────────

    def test_fasta_returns_dict_not_set(self):
        """_collect_sha256_set returns a dict (not a set) from FASTA fallback."""
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "test.fasta")
            _make_nnnx_fasta(fa, [(3, "ACGT"), (7, "TGCA")])
            result, n_legacy = self.mod._collect_sha256_set(fa)  # pylint: disable=protected-access
            self.assertIsInstance(result, dict)
            self.assertEqual(n_legacy, 0)

    def test_tsv_returns_dict_not_set(self):
        """_collect_sha256_set returns a dict from TSV fast path."""
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "test.fasta")
            tsv = os.path.join(d, "test.sha256_to_ids.tsv")
            _make_nnnx_fasta(fa, [(5, "AAAA"), (12, "CCCC")])
            h1 = hashlib.sha256(b"AAAA").hexdigest()
            h2 = hashlib.sha256(b"CCCC").hexdigest()
            with open(tsv, 'w', encoding='utf-8') as f:
                f.write(f"{h1}\t5\tid_a\tid_b\n")
                f.write(f"{h2}\t12\tid_c\n")
            # Ensure TSV is newer than FASTA
            time.sleep(0.05)
            os.utime(tsv, None)
            result, n_legacy = self.mod._collect_sha256_set(fa)  # pylint: disable=protected-access
            self.assertIsInstance(result, dict)
            self.assertEqual(n_legacy, 0)

    # ── NNNNx counts from FASTA ────────────────────────────────────────────────

    def test_fasta_nnnx_counts(self):
        """FASTA fallback extracts NNNNx count from the header prefix."""
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "test.fasta")
            _make_nnnx_fasta(fa, [(42, "GATTACA"), (999, "CATTAG")])
            result, _ = self.mod._collect_sha256_set(fa)  # pylint: disable=protected-access
            h1 = hashlib.sha256(b"GATTACA").hexdigest()
            h2 = hashlib.sha256(b"CATTAG").hexdigest()
            self.assertEqual(result[h1], 42)
            self.assertEqual(result[h2], 999)

    def test_tsv_nnnx_counts(self):
        """TSV fast path reads the count column as NNNNx."""
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "test.fasta")
            tsv = os.path.join(d, "test.sha256_to_ids.tsv")
            _make_nnnx_fasta(fa, [(1, "AAGG")])
            h1 = hashlib.sha256(b"AAGG").hexdigest()
            h2 = hashlib.sha256(b"TTCC").hexdigest()
            with open(tsv, 'w', encoding='utf-8') as f:
                f.write(f"{h1}\t350\tid1\n")
                f.write(f"{h2}\t17\tid2\n")
            time.sleep(0.05)
            os.utime(tsv, None)
            result, _ = self.mod._collect_sha256_set(fa)  # pylint: disable=protected-access
            self.assertEqual(result[h1], 350)
            self.assertEqual(result[h2], 17)

    # ── legacy GISAID file ─────────────────────────────────────────────────────

    def test_legacy_fasta_returns_empty_dict(self):
        """A plain GISAID FASTA with no sha256 IDs returns an empty dict."""
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "gisaid.fasta")
            _make_plain_fasta(fa, [
                ("EPI_ISL_100", "ATGC"),
                ("EPI_ISL_200", "GCTA"),
                ("EPI_ISL_300", "TTAA"),
            ])
            result, n_legacy = self.mod._collect_sha256_set(fa)  # pylint: disable=protected-access
            self.assertIsInstance(result, dict)
            self.assertEqual(len(result), 0)
            self.assertEqual(n_legacy, 3)

    # ── set operations on dict.keys() ──────────────────────────────────────────

    def test_dict_keys_set_difference(self):
        """dict.keys() supports set difference (child − parent)."""
        with tempfile.TemporaryDirectory() as d:
            parent_fa = os.path.join(d, "parent.fasta")
            child_fa = os.path.join(d, "child.fasta")
            # Parent has 3 sequences, child has 2 (one novel, one shared)
            _make_nnnx_fasta(parent_fa, [
                (100, "AAAA"), (200, "CCCC"), (300, "GGGG"),
            ])
            _make_nnnx_fasta(child_fa, [
                (150, "CCCC"), (50, "TTTT"),  # CCCC shared, TTTT novel
            ])
            parent_dict, _ = self.mod._collect_sha256_set(parent_fa)  # pylint: disable=protected-access
            child_dict, _ = self.mod._collect_sha256_set(child_fa)  # pylint: disable=protected-access
            novel_shas = child_dict.keys() - parent_dict.keys()
            self.assertEqual(len(novel_shas), 1)
            # The novel sha256 is TTTT's hash
            novel_sha = hashlib.sha256(b"TTTT").hexdigest()
            self.assertIn(novel_sha, novel_shas)
            # NNNNx sum of novel
            novel_nnnx = sum(child_dict[s] for s in novel_shas)
            self.assertEqual(novel_nnnx, 50)


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Novel + Traceable invariants
# ═══════════════════════════════════════════════════════════════════════════════

class TestNovelTraceableInvariants(unittest.TestCase):
    """Cross-checks: Novel + Traceable must equal Total for entries and NNNNx.

    Uses synthetic data with made-up sequences and counts.
    """

    @classmethod
    def setUpClass(cls):
        cls.mod = _load_module("summarize_fasta_pipeline.py")

    def _compute_novel_traceable(self, parent_entries, child_entries):
        """Build parent/child FASTAs and compute novel/traceable stats.

        parent_entries: list of (nnnx_count, seq_str)
        child_entries:  list of (nnnx_count, seq_str)

        Returns: (n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum)
        """
        with tempfile.TemporaryDirectory() as d:
            parent_fa = os.path.join(d, "parent.fasta")
            child_fa = os.path.join(d, "child.fasta")
            _make_nnnx_fasta(parent_fa, parent_entries)
            _make_nnnx_fasta(child_fa, child_entries)
            parent_dict, _ = self.mod._collect_sha256_set(parent_fa)  # pylint: disable=protected-access
            child_dict, _ = self.mod._collect_sha256_set(child_fa)  # pylint: disable=protected-access

            novel_shas = child_dict.keys() - parent_dict.keys()
            n_novel = len(novel_shas)
            novel_nnnx = sum(child_dict[s] for s in novel_shas)

            n_rec = len(child_entries)
            n_sum = sum(c for c, _ in child_entries)
            n_trace = n_rec - n_novel
            trace_nnnx = n_sum - novel_nnnx

            return n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum

    def test_all_shared_full_traceability(self):
        """When child is a subset of parent, all entries are traceable."""
        parent = [(500, "AAAA"), (200, "CCCC"), (80, "GGGG"), (10, "TTTT")]
        child = [(500, "AAAA"), (80, "GGGG")]  # strict subset
        n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum = \
            self._compute_novel_traceable(parent, child)
        self.assertEqual(n_novel, 0)
        self.assertEqual(novel_nnnx, 0)
        self.assertEqual(n_trace, n_rec)      # 2
        self.assertEqual(trace_nnnx, n_sum)   # 580
        # Invariant: novel + traceable == total
        self.assertEqual(n_novel + n_trace, n_rec)
        self.assertEqual(novel_nnnx + trace_nnnx, n_sum)

    def test_all_novel_zero_traceability(self):
        """When child has entirely new sequences, nothing is traceable."""
        parent = [(100, "AAAA"), (200, "CCCC")]
        child = [(30, "GATTACA"), (70, "CATTAG"), (15, "TAGGAT")]
        n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum = \
            self._compute_novel_traceable(parent, child)
        self.assertEqual(n_novel, 3)
        self.assertEqual(novel_nnnx, 115)   # 30 + 70 + 15
        self.assertEqual(n_trace, 0)
        self.assertEqual(trace_nnnx, 0)
        # Invariant
        self.assertEqual(n_novel + n_trace, n_rec)
        self.assertEqual(novel_nnnx + trace_nnnx, n_sum)

    def test_mixed_novel_and_shared(self):
        """Mix of shared and novel sequences verifies the invariant."""
        parent = [(1000, "AAAA"), (500, "CCCC"), (250, "GGGG")]
        child = [
            (800, "AAAA"),    # shared, different NNNNx
            (120, "TTTT"),    # novel
            (250, "GGGG"),    # shared, same NNNNx
            (45, "GATTACA"),  # novel
        ]
        n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum = \
            self._compute_novel_traceable(parent, child)
        self.assertEqual(n_rec, 4)
        self.assertEqual(n_sum, 800 + 120 + 250 + 45)   # 1215
        self.assertEqual(n_novel, 2)                     # TTTT + GATTACA
        self.assertEqual(novel_nnnx, 120 + 45)           # 165
        self.assertEqual(n_trace, 2)                     # AAAA + GGGG
        self.assertEqual(trace_nnnx, 800 + 250)          # 1050
        # Invariant
        self.assertEqual(n_novel + n_trace, n_rec)
        self.assertEqual(novel_nnnx + trace_nnnx, n_sum)

    def test_single_entry_shared(self):
        """Single child entry that exists in parent → fully traceable."""
        parent = [(77, "ACGT")]
        child = [(33, "ACGT")]
        n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum = \
            self._compute_novel_traceable(parent, child)
        self.assertEqual(n_novel, 0)
        self.assertEqual(n_trace, 1)
        self.assertEqual(trace_nnnx, 33)
        self.assertEqual(n_novel + n_trace, n_rec)
        self.assertEqual(novel_nnnx + trace_nnnx, n_sum)

    def test_single_entry_novel(self):
        """Single child entry NOT in parent → fully novel."""
        parent = [(77, "ACGT")]
        child = [(33, "TGCA")]
        n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum = \
            self._compute_novel_traceable(parent, child)
        self.assertEqual(n_novel, 1)
        self.assertEqual(novel_nnnx, 33)
        self.assertEqual(n_trace, 0)
        self.assertEqual(trace_nnnx, 0)
        self.assertEqual(n_novel + n_trace, n_rec)
        self.assertEqual(novel_nnnx + trace_nnnx, n_sum)

    def test_large_counts_invariant(self):
        """Invariant holds with large NNNNx values (similar to real data scale)."""
        parent = [
            (576521, "AAAA"), (513175, "CCCC"), (459607, "GGGG"),
            (420410, "TTTT"), (313235, "AACCGG"),
        ]
        child = [
            (576521, "AAAA"),   # shared
            (513175, "CCCC"),   # shared
            (1885, "GATTACA"),  # novel
            (420410, "TTTT"),   # shared
            (313235, "AACCGG"),  # shared
            (750, "CATTAG"),    # novel
        ]
        n_novel, novel_nnnx, n_trace, trace_nnnx, n_rec, n_sum = \
            self._compute_novel_traceable(parent, child)
        self.assertEqual(n_rec, 6)
        self.assertEqual(n_novel, 2)
        self.assertEqual(novel_nnnx, 1885 + 750)   # 2635
        self.assertEqual(n_trace, 4)
        expected_trace = 576521 + 513175 + 420410 + 313235  # 1,823,341
        self.assertEqual(trace_nnnx, expected_trace)
        # Core invariant
        self.assertEqual(n_novel + n_trace, n_rec)
        self.assertEqual(novel_nnnx + trace_nnnx, n_sum)

    def test_nonnegative_values(self):
        """All novel/traceable values are non-negative by construction."""
        parent = [(10, "AA"), (20, "CC")]
        child = [(5, "CC"), (3, "GG"), (8, "TT")]
        n_novel, novel_nnnx, n_trace, trace_nnnx, _, _ = \
            self._compute_novel_traceable(parent, child)
        self.assertGreaterEqual(n_novel, 0)
        self.assertGreaterEqual(novel_nnnx, 0)
        self.assertGreaterEqual(n_trace, 0)
        self.assertGreaterEqual(trace_nnnx, 0)


if __name__ == '__main__':
    unittest.main()
