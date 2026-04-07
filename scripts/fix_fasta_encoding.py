#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
r"""Normalize the character encoding of FASTA files in-place.

Background
----------
GISAID FASTA exports are aggregated from submissions made by thousands of
labs worldwide, each using different software, operating systems, and locale
settings.  As a result, a single FASTA file may contain header lines whose
non-ASCII characters were encoded in three distinct ways:

  1. **Valid UTF-8 multi-byte sequences** — the correct modern encoding.
     Example: Polish ``ś`` (U+015B) encoded as the two bytes ``\xC5 \x9B``,
     or ``ą`` (U+0105) as ``\xC4 \x85``.

  2. **Raw Latin-1 (ISO-8859-1) bytes** — single bytes ``\x80``–``\xFF``
     that represent Western European accented characters in the ISO-8859-1
     code page.  Example: French ``é`` as the single byte ``\xe9``, or
     Spanish ``ó`` as ``\xf3``.

  3. **Literal ``\uXXXX`` escape sequences** — the ASCII representation of
     Unicode code points produced by Java's ``native2ascii`` tool or by
     Python's ``repr()``/``unicode_escape`` codec.  These look like the
     six-character text string ``\u00e9`` rather than the real ``é`` glyph.

All three forms can appear *on the same FASTA header line*, because GISAID's
aggregator concatenated fields from labs that each had a different broken
pipeline.

This script normalises all three forms to proper UTF-8 in a single pass.

Encoding strategy
-----------------
A naïve ``try UTF-8 / except → Latin-1`` per-line approach breaks down when
a single line mixes valid UTF-8 multi-byte sequences (e.g. ``\xC5\x9B`` = ś)
with raw Latin-1 bytes (e.g. ``\xe9`` = é).  In that case the entire line
falls back to Latin-1, which corrupts the UTF-8 multi-byte pairs.

This script uses Python's ``surrogateescape`` error handler instead:

  Step 1 — ``raw.decode('utf-8', errors='surrogateescape')``:
      Valid UTF-8 multi-byte sequences are decoded to the correct Unicode
      character.  Every individual byte that is NOT part of a valid UTF-8
      sequence is mapped to a *surrogate code point* U+DC80..U+DCFF rather
      than raising an exception.

  Step 2 — Surrogate → Latin-1:
      Each surrogate U+DC80+b is converted back to ``chr(b)``, giving the
      Latin-1 interpretation of that raw byte.  This is lossless at the
      byte level: every byte of the original file is represented exactly
      once in the output.

  Step 3 — ``\uXXXX`` unescape:
      Literal six-character escape sequences ``\uXXXX`` are replaced with
      the real Unicode codepoint they represent.  A fast-path check
      (``'\\u' not in s``) makes this essentially free for the vast majority
      of lines (sequence lines, clean headers) that contain no escapes.

  Step 4 — C0 control character stripping (and TAB conversion):
      ASCII control characters 0x01–0x1F (excluding ``\\n`` 0x0A,
      ``\\r`` 0x0D) are stripped from every line. In addition, literal
      TAB (``\\t`` 0x09) characters are converted into standard spaces.

      **Real-world example 1 — GISAID record EPI_ISL_14481342:**

      The raw GISAID FASTA file contains literal TABs embedded internally
      inside the description field (likely a typo from the sequencing lab)::

          >...|INCARNATE WORD CONVENT - \t3400 BRADFORD STREET|...

      If left intact, downstream TSV mapping logs (which rely on tabs as column
      delimiters) would forcefully split this single ID into two columns during
      parsing, permanently corrupting sequence traceability in `summarize_fasta_pipeline.py`.
      This script explicitly converts that TAB into a space to maintain TSV integrity.

      **Real-world example 2 — GISAID record EPI_ISL_2016759:**

      The raw GISAID FASTA file contains the *literal six-character ASCII
      text* ``\\u0003`` (bytes 0x5C 0x75 0x30 0x30 0x30 0x33: backslash,
      u, 0, 0, 0, 3) embedded mid-header::

          >...Run20210508\u0003template bulk upload 20210508.xlsx...

      Confirmed by ``od -c`` on the raw file which shows ``\ u 0 0 0 3``
      (six separate characters), **not** ``\\003`` (the ETX byte).

      The encoding chain in this script:

      * Step 3 (``_unescape_unicode``) converts the literal ``\\u0003`` →
        the real ETX codepoint U+0003 (``\\x03``, one byte).
      * Step 4 (``_strip_c0_controls``) then strips U+0003.
      * Result: ``>...Run20210508template bulk upload 20210508.xlsx...``
        — the header is clean and complete on a single line.

      **Why this matters for downstream tools:**

      ``filterbyname.sh`` v39.06+ is a Java/Kotlin tool in BBTools.
      Java resolves ``\\uXXXX`` escape sequences *natively during string
      parsing*: when it reads the six ASCII characters ``\\u0003`` from the
      file, Java's string layer silently converts them to the real U+0003
      (ETX) codepoint.  ETX is then used by filterbyname.sh as a field or
      line delimiter, producing three failure modes:

      1. **Header truncation**: filterbyname.sh sees the ETX as the end of
         the record name, so the name it registers is truncated at
         ``'[Run20210508``.  The names-file search for ``EPI_ISL_2016759``
         (which appears before the ``\\u0003``) may still succeed OR fail
         depending on whether filterbyname.sh uses exact-name or
         pipe-field matching.

      2. **ignorejunk=t bypass**: if filterbyname.sh classifies the
         record as "junk" because of the ETX codepoint, the
         ``ignorejunk=t`` flag causes it to be forwarded to the output
         stream **without** consulting the names list.  The record
         therefore passes through even when it should have been filtered.

      3. **Malformed output FASTA**: when filterbyname.sh writes the
         extracted record, it emits a real ``\\n`` (LF) where the ETX was,
         splitting the header into two lines.  The tail of the original
         header (``template bulk upload...|Belgium``) is concatenated
         directly onto the nucleotide sequence without a line separator,
         producing a "sequence" that starts with non-nucleotide text.
         The SHA-256 of that mangled string does not match the expected
         sha256, so the downstream audit reports one spurious *Novel
         sha256* and the discarded-entry count is off by one.

      Stripping the ``\\u0003`` literal (via Steps 3 + 4) before any
      BBTools invocation resolves all three failure modes at the source.

Why ``unidecode`` failed
------------------------
The ``unidecode`` command-line tool reads its input as pure UTF-8.  When it
encounters a raw Latin-1 byte such as ``\xe9`` (é) — which is not a valid
UTF-8 continuation byte in isolation — it raises::

    Unable to decode input line NNNNNN: invalid continuation byte,
    start: 260, end: 261

Running this script first normalises the file to strict UTF-8, after which
``unidecode`` runs without errors (if ASCII-only output is required).

Character reference — ``\uXXXX`` escapes seen in GISAID data
-------------------------------------------------------------
The following escapes appear in real GISAID FASTA header lines.  All are
converted to their proper Unicode glyphs by this script:

    \u00c0  À   Latin capital A with grave          (French)
    \u00c1  Á   Latin capital A with acute
    \u00c2  Â   Latin capital A with circumflex
    \u00c3  Ã   Latin capital A with tilde
    \u00c4  Ä   Latin capital A with diaeresis       (German)
    \u00c5  Å   Latin capital A with ring above      (Scandinavian)
    \u00c6  Æ   Latin capital AE ligature
    \u00c7  Ç   Latin capital C with cedilla        (French, Portuguese)
    \u00c8  È   Latin capital E with grave
    \u00c9  É   Latin capital E with acute
    \u00ca  Ê   Latin capital E with circumflex
    \u00cb  Ë   Latin capital E with diaeresis
    \u00cc  Ì   Latin capital I with grave
    \u00cd  Í   Latin capital I with acute          (Spanish)
    \u00ce  Î   Latin capital I with circumflex
    \u00cf  Ï   Latin capital I with diaeresis
    \u00d1  Ñ   Latin capital N with tilde          (Spanish)
    \u00d3  Ó   Latin capital O with acute
    \u00d4  Ô   Latin capital O with circumflex
    \u00d5  Õ   Latin capital O with tilde
    \u00d6  Ö   Latin capital O with diaeresis      (German, Swedish)
    \u00d8  Ø   Latin capital O with stroke         (Norwegian)
    \u00d9  Ù   Latin capital U with grave
    \u00da  Ú   Latin capital U with acute
    \u00db  Û   Latin capital U with circumflex
    \u00dc  Ü   Latin capital U with diaeresis      (German)
    \u00df  ß   Latin small letter sharp S          (German)
    \u00e0  à   Latin small a with grave            (French, Italian)
    \u00e1  á   Latin small a with acute
    \u00e2  â   Latin small a with circumflex
    \u00e3  ã   Latin small a with tilde            (Portuguese)
    \u00e4  ä   Latin small a with diaeresis        (German)
    \u00e5  å   Latin small a with ring above       (Scandinavian)
    \u00e6  æ   Latin small ae ligature
    \u00e7  ç   Latin small c with cedilla         (French, Portuguese)
    \u00e8  è   Latin small e with grave            (French, Italian)
    \u00e9  é   Latin small e with acute            (French) ← very common
    \u00ea  ê   Latin small e with circumflex
    \u00eb  ë   Latin small e with diaeresis
    \u00ec  ì   Latin small i with grave
    \u00ed  í   Latin small i with acute            (Spanish)  ← common
    \u00ee  î   Latin small i with circumflex
    \u00ef  ï   Latin small i with diaeresis
    \u00f1  ñ   Latin small n with tilde            (Spanish)
    \u00f3  ó   Latin small o with acute            (Spanish, Polish)
    \u00f4  ô   Latin small o with circumflex
    \u00f5  õ   Latin small o with tilde
    \u00f6  ö   Latin small o with diaeresis        (German, Swedish) ← common
    \u00f8  ø   Latin small o with stroke           (Norwegian)
    \u00f9  ù   Latin small u with grave
    \u00fa  ú   Latin small u with acute
    \u00fb  û   Latin small u with circumflex
    \u00fc  ü   Latin small u with diaeresis        (German, Swedish) ← common
    \u00fd  ý   Latin small y with acute
    \u00ff  ÿ   Latin small y with diaeresis
    \u0105  ą   Latin small a with ogonek           (Polish)
    \u0107  ć   Latin small c with acute            (Polish)
    \u010d  č   Latin small c with caron            (Czech, Slovak)
    \u0111  đ   Latin small d with stroke           (Croatian)
    \u0119  ę   Latin small e with ogonek           (Polish)
    \u011b  ě   Latin small e with caron            (Czech)
    \u0141  Ł   Latin capital L with stroke        (Polish) ← e.g. ZAKŁAD
    \u0142  ł   Latin small l with stroke           (Polish)
    \u0144  ń   Latin small n with acute            (Polish)
    \u015b  ś   Latin small s with acute            (Polish)
    \u015f  ş   Latin small s with cedilla          (Turkish, Romanian)
    \u0159  ř   Latin small r with caron            (Czech) ← e.g. Kovaříková, laboratoř
    \u0160  Š   Latin capital S with caron          (Czech, Slovak)
    \u0161  š   Latin small s with caron            (Czech, Slovak)
    \u0179  Ź   Latin capital Z with acute          (Polish) ← e.g. KOŹLE
    \u017a  ź   Latin small z with acute            (Polish)
    \u017b  Ż   Latin capital Z with dot above      (Polish) ← e.g. Żarach, Żaganiu
    \u017c  ż   Latin small z with dot above        (Polish)
    \u017e  ž   Latin small z with caron            (Czech, Slovak)
    \u0171  ű   Latin small u with double acute     (Hungarian)
    \u0151  ő   Latin small o with double acute     (Hungarian)
    \u00e3  ã   Latin small a with tilde            (Portuguese: São Paulo)
    \u0219  ș   Latin small s with comma below      (Romanian) ← e.g. Babeș, Timișoara
    \u021b  ț   Latin small t with comma below      (Romanian)
    \u00eb  ë   Latin small e with diaeresis        (Albanian/Kosovo) ← e.g. Prishtinë
    \u00a0      Non-breaking space (NBSP)           (appears in some Italian lab names)
    \u201c  "   Left double quotation mark          (Italian hospital names)
    \u201d  "   Right double quotation mark
    \u2018  '   Left single quotation mark
    \u2019  '   Right single quotation mark / apostrophe  ← common in Italian
                    (e.g. dell\u2019Abruzzo, dell\u2019Aquila)
    \u2013  –   En dash                             (Greek hospital names)
    \u05d3  ד   Hebrew letter Dalet                 (Israel)
    \u05d5  ו   Hebrew letter Vav                   (Israel)
    \u05d9  י   Hebrew letter Yod                   (Israel)
    \u05dc  ל   Hebrew letter Lamed                 (Israel)
    \u05de  מ   Hebrew letter Mem                   (Israel)
    \u05df  ן   Hebrew letter Final Nun             (Israel)
    \u05e2  ע   Hebrew letter Ayin                  (Israel)
    \u05ea  ת   Hebrew letter Tav                   (Israel)
    \u0e3a  ฺ   Thai character Phinthu              (Thailand: Bhumibol adulyadej)
    \u00b2  ²   Superscript Two                     (USA: x²|Utah)

Real-world examples (from spikenuc1207.fasta dry-run)
------------------------------------------------------
Each ``~`` line shows the raw header as stored in the FASTA file (with
literal ``\uXXXX`` sequences).  Each ``+`` line shows the normalised result.

French (Québec, Canada)::

    ~ >...Hospital Universitario de Gran Canaria Dr. Negr\u00edn|...
    + >...Hospital Universitario de Gran Canaria Dr. Negrín|...

    ~ >...Laboratoire de sant\u00e9 publique du Qu\u00e9bec|...
    + >...Laboratoire de santé publique du Québec|...

German (Austria)::

    ~ >...Institut f\u00fcr Virologie am Department f\u00fcr Hygiene|...
    + >...Institut für Virologie am Department für Hygiene|...

    ~ >...Ludwig Boltzmann Institut f\u00fcr Experimentelle und Klinische Traumatologie|...
    + >...Ludwig Boltzmann Institut für Experimentelle und Klinische Traumatologie|...

Swedish (Sweden)::

    ~ >...The Public Health Agency of Sweden|Svartstr\u00f6m|...
    + >...The Public Health Agency of Sweden|Svartström|...

Portuguese (Brazil)::

    ~ >...hCoV-19^^S\u00e3o Paulo|Human|Diagn\u00f3sticos da Am\u00e9rica - DASA|...
    + >...hCoV-19^^São Paulo|Human|Diagnósticos da América - DASA|...

    ~ >...Laborat\u00f3rio de Virologia - Instituto de Medicina Tropical - Universidade de S\u00e3o Paulo|...
    + >...Laboratório de Virologia - Instituto de Medicina Tropical - Universidade de São Paulo|...

Italian (Italy)::

    ~ >...AOU Policlinico Umberto I; Sapienza Universit\u00e0 di Roma|...
    + >...AOU Policlinico Umberto I; Sapienza Università di Roma|...

    ~ >...Fondazione Policlinico Universitario \u201cA. Gemelli\u201d IRCCS|...
    + >...Fondazione Policlinico Universitario "A. Gemelli" IRCCS|...

Belgian French::

    ~ >...Institut de Pathologie et G\u00e9n\u00e9tique (IPEquivalent manual pipeline
--------------------------
The following shell commands achieve the same result but require a Java
runtime (``native2ascii``) and may fail on mixed-encoding files::

    # Convert \\uXXXX escapes to real Unicode (Java native2ascii):
    native2ascii -encoding UTF-8 -reverse input.fasta output.fasta

    # Optionally transliterate remaining non-ASCII to ASCII:
    unidecode output.fasta > output.fasta.tmp && mv output.fasta.tmp output.fasta

   Step 4 — C0 control character stripping:
       The raw GISAID FASTA file for EPI_ISL_2016759 contains the *literal
       six-character text* ``\u0003`` (as verified by ``od -c``: ``\ u 0 0 0 3``,
       not the one-byte ETX ``\003``).  ``filterbyname.sh`` v39.06+ (Java)
       resolves ``\uXXXX`` escapes natively, converting ``\u0003`` to real
       ETX, then uses ETX as a field delimiter.  This truncates the header,
       bypasses name-list filtering via ``ignorejunk=t``, and produces a
       malformed output FASTA where the header tail is concatenated onto the
       sequence, causing an erroneous *Novel sha256* in the audit report.
       Stripping the ``\u0003`` literal here prevents all three failure modes.


Usage examples
--------------
::

    # Preview changes without touching any file (safe, recommended first step):
    fix_fasta_encoding.py --dry-run spikenuc1207.fasta

    # Preview with counts only (no per-line diff):
    fix_fasta_encoding.py --stats-only spikenuc1207.fasta

    # Preview with per-line diffs printed to stdout (verbose):
    fix_fasta_encoding.py --dry-run --verbose spikenuc1207.fasta

    # Fix a single file with scp-style progress meter on stderr:
    fix_fasta_encoding.py --progress spikenuc1207.fasta

    # Fix a single file (renames .fasta → .fasta.orig, writes clean UTF-8):
    fix_fasta_encoding.py spikenuc1207.fasta

    # Fix all FASTA files in a directory:
    fix_fasta_encoding.py /data/seqs/*.fasta

    # Fix all FASTA files recursively (with bash globstar):
    fix_fasta_encoding.py /data/seqs/**/*.fasta

    # Re-process a file that was already fixed (overwrites the .orig backup):
    fix_fasta_encoding.py --overwrite spikenuc1207.fasta

    # After fixing, regenerate pipeline artefacts:
    #   rm *.sha256_to_ids.tsv *.sha256_to_descr_lines.tsv
    #   rm *.discarded_sha256_hashes.txt *.discarded_original_ids.txt
    #   summarize_fasta_pipeline.py /data/seqs/ spikenuc1207 --full-fasta-header

    # If ASCII-only output is needed, run unidecode after this script:
    fix_fasta_encoding.py spikenuc1207.fasta
    unidecode spikenuc1207.fasta > spikenuc1207.ascii.fasta

Progress display (--progress)
------------------------------
The ``--progress`` flag (or automatic detection of a TTY stderr) prints a
single-line scp/dd-style progress gauge that updates in-place::

    spikenuc1207.fasta           87%    55.1 GB   15.4 MB/s  0:09:45  1,231,006 changed

The line is hard-clipped to the terminal width (max 120 chars) so it never
wraps and the carriage-return overwrite works correctly.

Log file (always created)
--------------------------
Every invocation appends a full diagnostic report to::

    {stem}.fix_fasta_encoding.{YYYYMMDD_HHMMSS}.log

where ``{stem}`` is the input path with any FASTA suffix stripped
(``.fasta``, ``.fasta.orig``, ``.fasta.ori``, ``.fasta.old``).  The
datetime is fixed once per invocation so all files processed in one run
share the same log timestamp.

Example log for ``spikenuc1207.fasta`` (63.4 GB, processed 2026-04-04)::

    === 2026-04-04T22:23:27.123456  spikenuc1207.fasta ===
      Wrote spikenuc1207.fasta: 1,413,630 line(s) changed, 14.4 MB/s
      Encoding-mix breakdown (per changed line):
        \\uXXXX escapes only                                 : 1,398,889
        Valid UTF-8 multi-byte only                         :   14,646
        Latin-1 raw bytes  +  valid UTF-8 *** MIXED ***     :       95
      \\uXXXX codepoint frequency (top 30, across all changed lines):
        U+00E9  é  :  504,307
        U+00F3  ó  :  294,344
        U+00FC  ü  :  288,911
        U+00ED  í  :  229,251
        U+00F6  ö  :  165,809
        U+00E4  ä  :   96,839
        U+00D6  Ö  :   90,164
        U+00FA  ú  :   75,447
        U+00E1  á  :   53,409
        U+00A0     :   49,139
        U+0107  ć  :   48,264
        U+2013  –  :   41,515
        U+0161  š  :   39,467
        U+00F1  ñ  :   36,954
        U+010D  č  :   36,863
        U+017E  ž  :   29,029
        U+0308  ̈  :   28,134
        U+00E0  à  :   27,439
        U+00E5  å  :   26,204
        U+00E7  ç  :   24,738
        U+00DC  Ü  :   24,457
        U+00E8  è  :   23,634
        U+0301  ́  :   22,184
        U+00F4  ô  :   20,436
        U+0101  ā  :   18,870
        U+2019  '  :   18,491
        U+0144  ń  :   17,946
        U+201C  "  :   15,531
        U+201D  "  :   14,944
        U+0142  ł  :   12,829
      WARNING — sample lines (Valid UTF-8 multi-byte only):
        '>...Universidade Federal de SÃ£o Paulo...'
        '>...Instituto Nacional de Salud- DirecciÃ³n de InvestigaciÃ³n...'
        '>...VeselÄ\u00abvas Centrs 4...'
        '>...MVZ fÃ¼r Laboratoriumsmedizin...'
      WARNING — sample lines (Latin-1 raw bytes  +  valid UTF-8 *** MIXED ***):
        '>...HÃ\u00b4pitaux de Paris, UniversitÃ© Paris-Est...'
      Renamed spikenuc1207.fasta -> spikenuc1207.fasta.orig
      Written clean UTF-8 -> spikenuc1207.fasta

Notes on the warning categories
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``Valid UTF-8 multi-byte only`` (14,646 lines in example above):
    These headers contain **raw UTF-8 bytes** but *no* ``\\uXXXX`` escapes.
    The garbled-looking ``SÃ£o`` / ``fÃ¼r`` patterns are classic Windows-1252
    Mojibake: the bytes of a UTF-8 sequence were stored as if each byte were
    a separate Latin-1 character, then submitted to GISAID as-is.  The
    surrogate-escape decoder re-encodes them to well-formed UTF-8 output.
    The result is still Mojibake (the source data was already corrupt), but
    downstream tools will at least not crash on UnicodeDecodeError.

``Latin-1 raw bytes + valid UTF-8 *** MIXED ***`` (95 lines in example):
    The most pathological case: a header whose fields were assembled from
    sources with entirely different encodings — some using raw ``\\xNN``
    Latin-1 bytes, others submitting valid multi-byte UTF-8 sequences.
    The surrogate-escape strategy handles all bytes losslessly.

``\\uXXXX escapes + valid UTF-8 multi-byte``:
    Headers where GISAID's export added ``\\uXXXX`` escapes for some
    characters while separately including raw UTF-8 for others (typically
    from different metadata fields).
"""

import argparse
import hashlib
import os
import re
import shutil
import sys
import tempfile
import time
from collections import Counter
from datetime import datetime

VERSION = "202604030000"


def _get_git_version() -> str:
    """Return ``git describe --always --dirty --tags`` output, or ``'unknown'``."""
    import subprocess as _sp  # local import — subprocess not needed at module level
    _here = os.path.dirname(os.path.abspath(__file__))
    try:
        result = _sp.run(
            ["git", "describe", "--always", "--dirty", "--tags"],
            capture_output=True, text=True, check=True,
            cwd=_here,
        )
        ver = result.stdout.strip()
        if ver:
            return ver
    except Exception:  # pylint: disable=broad-except
        pass
    env_ver = os.environ.get("GIT_COMMIT", "").strip()
    return env_ver[:12] if env_ver else "unknown"


_GIT_VERSION: str = _get_git_version()

_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
_parser.add_argument(
    "infiles", nargs="+", metavar="FASTA",
    help="One or more FASTA files to normalize in-place.",
)
_parser.add_argument(
    "--dry-run", action="store_true",
    help=(
        "Show what would be done without modifying any files.  "
        "Prints a summary of changed lines per file."
    ),
)
_parser.add_argument(
    "--overwrite", action="store_true",
    help=(
        "Re-process files even when a .orig backup already exists.  "
        "The existing .orig is overwritten."
    ),
)
_parser.add_argument(
    "--stats-only", action="store_true",
    help=(
        "Like --dry-run but suppresses per-line output; only prints the "
        "per-file summary line."
    ),
)
_parser.add_argument(
    "--verbose", "-v", action="store_true",
    help=(
        "Print a per-line diff of every changed header line.  Works in "
        "both dry-run and real-write modes."
    ),
)
_parser.add_argument(
    "--progress", action="store_true",
    help=(
        "Show a scp/dd-style progress line with bytes read, percentage, "
        "throughput, ETA, and number of lines changed.  Enabled "
        "automatically when stderr is a TTY; this flag forces it also "
        "when stderr is redirected (e.g. to a log file)."
    ),
)
_parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")


# ── helpers ──────────────────────────────────────────────────────────────────

_PROGRESS_INTERVAL = 0.2   # seconds between progress line updates


def _fmt_size(n: int) -> str:
    """Return *n* bytes formatted as a human-readable string (B/KB/MB/GB)."""
    for unit, threshold in (('GB', 1 << 30), ('MB', 1 << 20), ('KB', 1 << 10)):
        if n >= threshold:
            return f"{n / threshold:.1f} {unit}"
    return f"{n} B"


def _fmt_speed(bps: float) -> str:
    """Return *bps* bytes/second as a human-readable speed string."""
    return _fmt_size(int(bps)) + "/s"


def _file_sha256_stream(path: str) -> str:
    """Return sha256 hex digest of *path* by reading in 4 MB chunks (streaming)."""
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        while chunk := fh.read(4 << 20):
            h.update(chunk)
    return h.hexdigest()


def _copy_primary_clean(
        dup_path: str,
        primary_path: str,
        dry_run: bool,
        overwrite: bool,
        log_fh,
) -> int:
    """Copy the already-clean output of *primary_path* to *dup_path*.

    Both files have identical original content (confirmed by size + sha256)
    so the cleaned output is identical.  This avoids re-running the full
    Python per-line encoding pass on the duplicate.

    Returns 1 if the copy was performed (or would be in dry-run), 0 if skipped.
    """
    def _emit(msg: str = '') -> None:
        print(msg, file=sys.stderr)
        if log_fh is not None:
            print(msg, file=log_fh)

    dup_orig = dup_path + '.orig'
    if os.path.exists(dup_orig) and not overwrite:
        _emit(
            f"  Skipping {dup_path}: backup {os.path.basename(dup_orig)} already "
            "exists (use --overwrite to force)."
        )
        return 0
    if not os.path.exists(primary_path):
        _emit(f"  Error: primary clean file {primary_path!r} not found — "
              f"cannot copy to {dup_path!r}.")
        return 0
    size_str = _fmt_size(os.path.getsize(primary_path))
    action = 'Would copy' if dry_run else 'Copying'
    _emit(
        f"  {action} clean output: "
        f"{os.path.basename(primary_path)} \u2192 {os.path.basename(dup_path)} ({size_str})"
    )
    if dry_run:
        return 1
    os.rename(dup_path, dup_orig)
    _emit(f"  Renamed {os.path.basename(dup_path)} \u2192 {os.path.basename(dup_orig)}")
    shutil.copyfile(primary_path, dup_path)
    _emit(f"  Copied clean output \u2192 {os.path.basename(dup_path)}")
    return 1


def _fmt_eta(seconds: float) -> str:
    """Return ETA as 'h:mm:ss' or 'mm:ss'.  Returns '--:--' when unknown."""
    if seconds < 0 or seconds > 86400 * 2:
        return '--:--'
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


def _progress_line(path: str, bytes_done: int, total_bytes: int,
                   changed: int, elapsed: float) -> str:
    """Build a single-line progress string in the style of scp / dd.

    Example output (fits in 80 columns)::

        spikenuc1207.fasta   87%  112 MB   14.2 MB/s   00:08  1,234 changed
    """
    pct = int(bytes_done / total_bytes * 100) if total_bytes else 0
    speed = bytes_done / elapsed if elapsed > 1e-6 else 0
    eta = (total_bytes - bytes_done) / speed if speed > 0 else 0
    cols = min(shutil.get_terminal_size(fallback=(80, 24)).columns, 120)
    name = os.path.basename(path)
    # Cap filename column; the fixed stats portion needs ~52 chars.
    max_name = max(10, min(40, cols - 52))
    if len(name) > max_name:
        name = name[:max_name - 1] + '…'
    line = (
        f"{name:<{max_name}}  {pct:3d}%  {_fmt_size(bytes_done):>8}  "
        f"{_fmt_speed(speed):>11}  {_fmt_eta(eta):>7}  {changed:,} changed"
    )
    # Hard-clip to terminal width so the line never wraps.
    return line[:cols]


def _unescape_unicode(s: str) -> str:
    r"""Convert literal \uXXXX escape sequences to real Unicode characters.

    Fast path: if '\\u' is not in the string the function returns immediately
    with negligible overhead.
    """
    if "\\u" not in s:
        return s
    return re.sub(
        r"\\u([0-9a-fA-F]{4})",
        lambda m: chr(int(m.group(1), 16)),
        s,
    )


# C0 control characters to strip: 0x01–0x1F minus TAB (0x09), LF (0x0A), CR (0x0D).
_C0_STRIP = re.compile(r'[\x01-\x08\x0b\x0c\x0e-\x1f]')


def _strip_c0_controls(s: str) -> str:
    """Remove ASCII C0 control bytes (0x01–0x1f, excluding \\n/\\r) from *s*.
    Converts \\t (TAB) to a space.

    These bytes have no valid meaning in a FASTA file and can break downstream
    parsers that use them as implicit record delimiters (e.g. BBTools
    filterbyname.sh treats ETX \\x03 as a header terminator).
    """
    # Specifically convert literal TABs to spaces to prevent TSV mapping splits downstream
    s = s.replace('\t', ' ')
    if not _C0_STRIP.search(s):
        return s
    return _C0_STRIP.sub('', s)


def _decode_line(raw: bytes) -> str:
    r"""Decode one raw FASTA byte line to a clean Unicode str.

    Handles GISAID FASTA files that mix UTF-8 and Latin-1 *within the same
    line* (e.g. Polish ``ś`` as UTF-8 \\xC5\\x9B and French ``é`` as raw
    Latin-1 \\xe9 in the same header):

      1. ``decode('utf-8', errors='surrogateescape')``: valid multi-byte UTF-8
         sequences decode normally; every individual byte that is not part of
         a valid UTF-8 sequence becomes a surrogate code point U+DC80..U+DCFF.
      2. Each surrogate U+DC80+b is mapped back to chr(b), giving the Latin-1
         interpretation of that byte.  This is lossless — all bytes end up
         correctly decoded.
      3. Literal ``\\uXXXX`` escape sequences (produced by some GISAID export
         pipelines) are then converted to real Unicode codepoints.

    This is why ``unidecode`` previously failed: it received a file where
    some bytes were valid UTF-8 multi-byte sequences and others were raw
    Latin-1 bytes, so any attempt to decode the whole line as pure UTF-8
    raised UnicodeDecodeError.
    """
    s = raw.decode("utf-8", errors="surrogateescape")
    if any(0xDC80 <= ord(ch) <= 0xDCFF for ch in s):
        s = "".join(
            chr(ord(ch) - 0xDC00) if 0xDC80 <= ord(ch) <= 0xDCFF else ch
            for ch in s
        )
    s = _unescape_unicode(s)
    return _strip_c0_controls(s)


def _make_log_path(path: str, dt_str: str) -> str:
    """Return the log file path for *path*.

    Strips known FASTA suffixes (including .orig/.ori) so that
    ``spikenuc1207.fasta`` and ``spikenuc1207.fasta.orig`` both produce
    ``spikenuc1207.fix_fasta_encoding.<dt>.log``.
    """
    stem = path
    for sfx in ('.fasta.orig', '.fasta.ori', '.fasta.old', '.fasta'):
        if stem.endswith(sfx):
            stem = stem[:-len(sfx)]
            break
    return f"{stem}.fix_fasta_encoding.{dt_str}.log"


def _process_file(path: str, dry_run: bool, overwrite: bool,
                  stats_only: bool, verbose: bool, progress: bool,
                  log_fh) -> int:
    """Normalize *path* in-place.  Returns the number of lines changed."""

    def _emit(msg: str = "", **kwargs) -> None:
        """Print *msg* to stderr and tee to *log_fh* (no progress animation)."""
        print(msg, file=sys.stderr, **kwargs)
        if log_fh is not None:
            # Strip the keyword 'end' if present; always use newline in log.
            print(msg, file=log_fh)

    orig_path = path + ".orig"

    if os.path.exists(orig_path) and not overwrite:
        _emit(
            f"  Skipping {path}: backup {os.path.basename(orig_path)} already "
            "exists (use --overwrite to force)."
        )
        return 0

    # Show progress when stderr is a TTY or when --progress forces it.
    show_progress = progress or sys.stderr.isatty()

    total_bytes = os.path.getsize(path)
    bytes_read = 0
    lines_read = 0
    changed = 0
    # Stream output to a temp file immediately — avoids holding the whole
    # file in RAM (input files can be tens of GB).
    # Create a uniquely-named temp file in the same directory so that
    # os.rename() is guaranteed to be atomic (same filesystem).
    tmp_fd, tmp_path = tempfile.mkstemp(
        dir=os.path.dirname(os.path.abspath(path)),
        prefix=".fix_enc_",
        suffix=".tmp",
    )
    # Encoding-mix statistics: counts how many changed lines had each
    # combination of the three encoding forms.
    #   'esc'  — literal \uXXXX escape sequences
    #   'lat1' — raw Latin-1 bytes (not valid UTF-8)
    #   'utf8' — valid UTF-8 multi-byte sequences
    mix_counter: Counter = Counter()
    # Per-codepoint frequency: how many changed lines contained each \uXXXX escape.
    cp_counter: Counter = Counter()
    # Sample lines for each non-trivial encoding category (max 50 per type).
    _max_samples = 50
    warn_samples: dict = {}   # frozenset -> list[str]
    warn_extras: dict = {}    # frozenset -> list[str] (for fallback)
    t_start = time.monotonic()
    t_last = t_start

    with open(path, "rb") as fh, \
            os.fdopen(tmp_fd, "w", encoding="utf-8") as fh_tmp:
        for raw in fh:
            bytes_read += len(raw)
            lines_read += 1
            # Preserve the original line ending (strip nothing here).
            original_text = raw.decode("latin-1")   # lossless reference decode
            clean = _decode_line(raw)
            clean_line = clean if clean.endswith("\n") else clean.rstrip("\r\n") + "\n"
            fh_tmp.write(clean_line)
            if clean.rstrip("\r\n") != original_text.rstrip("\r\n"):
                changed += 1
                # Classify which encoding forms are present on this line.
                decoded_surr = raw.decode('utf-8', errors='surrogateescape')
                enc_types: set[str] = set()
                if re.search(r'\\u[0-9a-fA-F]{4}', original_text):
                    enc_types.add('esc')
                if any(0xDC80 <= ord(ch) <= 0xDCFF for ch in decoded_surr):
                    enc_types.add('lat1')
                if any(ord(ch) > 0x7F and not (0xDC80 <= ord(ch) <= 0xDCFF)
                       for ch in decoded_surr):
                    enc_types.add('utf8')
                if any(0x01 <= b <= 0x1F and b not in (0x09, 0x0A, 0x0D)
                       for b in raw):
                    enc_types.add('ctrl')
                mix_counter[frozenset(enc_types)] += 1
                # Count individual \uXXXX codepoints on this line.
                for cp_hex in re.findall(r'\\u([0-9a-fA-F]{4})', original_text):
                    cp_counter[int(cp_hex, 16)] += 1
                # Collect a sample for any non-pure-esc (warning) category.
                key = frozenset(enc_types)
                if key != frozenset({'esc'}) and key != frozenset():
                    samples = warn_samples.setdefault(key, [])
                    extras = warn_extras.setdefault(key, [])
                    candidate = original_text.rstrip()
                    if len(samples) < _max_samples:
                        # Force visual diversity: reject candidate if it shares a 35-char prefix with an existing sample
                        if not any(candidate[:35] == s[:35] for s in samples):
                            samples.append(candidate)
                        elif len(extras) < _max_samples:
                            extras.append(candidate)

                if verbose or (not stats_only and dry_run):
                    # Print diff line; clear the progress line first if active.
                    if show_progress:
                        print("", file=sys.stderr)  # move off progress row
                        show_progress = False        # don't overwrite the diff
                    print(
                        f"  ~ {original_text.rstrip()!r}\n"
                        f"  + {clean.rstrip()!r}",
                    )

            # Emit progress update every _PROGRESS_INTERVAL seconds.
            now = time.monotonic()
            if show_progress and (now - t_last) >= _PROGRESS_INTERVAL:
                elapsed = now - t_start
                line = _progress_line(path, bytes_read, total_bytes,
                                      changed, elapsed)
                print(f"\r{line}", end="", flush=True, file=sys.stderr)
                t_last = now

    # Final progress update / line clear.
    if show_progress:
        elapsed = time.monotonic() - t_start
        if changed:
            line = _progress_line(path, total_bytes, total_bytes,
                                  changed, elapsed)
            print(f"\r{line}", flush=True, file=sys.stderr)
        else:
            # Clear the progress line if nothing changed.
            width = shutil.get_terminal_size(fallback=(80, 24)).columns
            print(f"\r{' ' * width}\r", end="", flush=True, file=sys.stderr)

    # If nothing changed or dry-run, remove the temp file.
    if changed == 0:
        os.unlink(tmp_path)
        _emit(f"  {path}: no changes needed.")
        return 0

    elapsed = time.monotonic() - t_start
    speed_str = _fmt_speed(total_bytes / elapsed) if elapsed > 1e-6 else ""
    action = "Would write" if dry_run else "Wrote"
    _emit(
        f"  {action} {path}: {changed:,} line(s) changed"
        + (f", {speed_str}" if speed_str else "")
    )

    # Encoding-mix summary — always printed when there are changes.
    _mix_labels = [
        (frozenset({'esc'}),              "\\uXXXX escapes only"),
        (frozenset({'lat1'}),             "Latin-1 raw bytes only"),
        (frozenset({'utf8'}),             "Valid UTF-8 multi-byte only"),
        (frozenset({'ctrl'}),             "C0 control characters only (e.g. ETX \\x03)"),
        (frozenset({'esc', 'utf8'}),      "\\uXXXX escapes  +  valid UTF-8 multi-byte"),
        (frozenset({'esc', 'lat1'}),      "\\uXXXX escapes  +  Latin-1 raw bytes"),
        (frozenset({'lat1', 'utf8'}),     "Latin-1 raw bytes  +  valid UTF-8 *** MIXED ***"),
        (frozenset({'ctrl', 'esc'}),      "C0 control  +  \\uXXXX escapes"),
        (frozenset({'ctrl', 'lat1'}),     "C0 control  +  Latin-1 raw bytes"),
        (frozenset({'ctrl', 'utf8'}),     "C0 control  +  valid UTF-8 multi-byte"),
        (frozenset({'esc', 'lat1', 'utf8'}),
         "All three (\\uXXXX + Latin-1 + UTF-8) *** MIXED ***"),
    ]
    if mix_counter:
        _emit("  Encoding-mix breakdown (per changed line):")
        for key, label in _mix_labels:
            cnt = mix_counter.get(key, 0)
            if cnt:
                _emit(f"    {label:<52s}: {cnt:>8,}")
        # Print any unexpected combinations not in the label table.
        for key, cnt in mix_counter.items():
            if not any(key == k for k, _ in _mix_labels):
                _emit(f"    (other: {sorted(key)}): {cnt:>8,}")

    # Codepoint frequency table (most common first, top 30).
    if cp_counter:
        _emit("  \\uXXXX codepoint frequency (top 30, across all changed lines):")
        for cp, cnt in cp_counter.most_common(30):
            ch = chr(cp)
            _emit(f"    U+{cp:04X}  {ch}  : {cnt:>8,}")

    # Warning samples for non-pure-esc lines (mixed / double-encoded).
    _warn_keys = [
        frozenset({'utf8'}),
        frozenset({'lat1'}),
        frozenset({'ctrl'}),
        frozenset({'esc', 'utf8'}),
        frozenset({'esc', 'lat1'}),
        frozenset({'ctrl', 'esc'}),
        frozenset({'ctrl', 'lat1'}),
        frozenset({'ctrl', 'utf8'}),
        frozenset({'lat1', 'utf8'}),
        frozenset({'esc', 'lat1', 'utf8'}),
    ]
    for key in _warn_keys:
        samples = warn_samples.get(key, [])
        extras = warn_extras.get(key, [])
        # If we couldn't find 5 completely divergent examples, fill the rest from the redundant buffer.
        while len(samples) < _max_samples and extras:
            samples.append(extras.pop(0))

        if not samples:
            continue
        label = next((lb for k, lb in _mix_labels if k == key), str(sorted(key)))
        _emit(f"  WARNING \u2014 sample lines ({label}):")
        for s in samples:
            _emit(f"    {s!r}")

    if dry_run:
        os.unlink(tmp_path)
        return changed

    # Temp file is already written; just rename into place.
    # The temp file lives in the same directory so both os.rename() calls
    # are on the same filesystem (guaranteed atomic on POSIX).
    # Backup original (or overwrite existing backup if --overwrite).
    os.rename(path, orig_path)
    _emit(f"  Renamed {os.path.basename(path)} -> {os.path.basename(orig_path)}")
    os.rename(tmp_path, path)
    _emit(f"  Written clean UTF-8 -> {os.path.basename(path)}")
    return changed


# ── entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    """Normalise FASTA header encoding for each input file."""
    _start_ts = datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    print(
        f"{_start_ts} fix_fasta_encoding.py"
        f"  version {VERSION}  git:{_GIT_VERSION}"
        f"  invoked: {' '.join(sys.argv)}",
        file=sys.stderr,
    )
    opts = _parser.parse_args()
    dt_str = datetime.now().strftime("%Y%m%d_%H%M%S")

    # ── detect identical inputs (size + sha256) ────────────────────────────
    # Identical files produce identical clean output.  We process the first
    # occurrence (primary) normally; for identical duplicates (twins) we just
    # rename the original to .orig and copy the primary's clean output.
    # This avoids re-running the full Python per-line pass on every twin.
    valid_paths = [p for p in opts.infiles if os.path.isfile(p)]
    invalid_paths = [p for p in opts.infiles if not os.path.isfile(p)]

    # Only consider files that will actually be processed (no .orig yet, or --overwrite).
    to_check = [p for p in valid_paths
                if not os.path.exists(p + '.orig') or opts.overwrite]

    by_size: dict[int, list[str]] = {}
    for p in to_check:
        by_size.setdefault(os.path.getsize(p), []).append(p)

    copy_from: dict[str, str] = {}   # twin_path → primary_path
    for sz, size_group in by_size.items():
        if len(size_group) < 2 or sz == 0:
            continue
        print(
            f"Info: {len(size_group)} file(s) share size {_fmt_size(sz)} — "
            "computing sha256 to detect identical inputs…",
            file=sys.stderr,
        )
        sha_to_primary: dict[str, str] = {}
        for p in size_group:
            sha = _file_sha256_stream(p)
            if sha in sha_to_primary:
                twin_of = sha_to_primary[sha]
                copy_from[p] = twin_of
                print(
                    f"  {os.path.basename(p)} ≡ {os.path.basename(twin_of)} "
                    f"(sha256={sha[:16]}…) — will copy clean output instead of reprocessing.",
                    file=sys.stderr,
                )
            else:
                sha_to_primary[sha] = p

    # Primaries first (ensures clean output exists before twins are processed).
    primaries = [p for p in valid_paths if p not in copy_from]
    twins = [p for p in valid_paths if p in copy_from]
    ordered = primaries + twins

    total_files = 0
    total_changed = 0

    for path in invalid_paths:
        print(f"Warning: not a file, skipping: {path}", file=sys.stderr)

    for path in ordered:
        log_path = _make_log_path(path, dt_str)
        print(f"{path}:  (log -> {log_path})", file=sys.stderr)
        with open(log_path, "a", encoding="utf-8") as log_fh:
            log_fh.write(f"\n=== {datetime.now().isoformat()}  {path} ===\n")
            if path in copy_from:
                n = _copy_primary_clean(
                    path,
                    primary_path=copy_from[path],
                    dry_run=opts.dry_run,
                    overwrite=opts.overwrite,
                    log_fh=log_fh,
                )
            else:
                n = _process_file(
                    path,
                    dry_run=opts.dry_run,
                    overwrite=opts.overwrite,
                    stats_only=opts.stats_only,
                    verbose=opts.verbose,
                    progress=opts.progress,
                    log_fh=log_fh,
                )
        total_files += 1
        total_changed += n

    print(
        f"\nDone: {total_files} file(s) processed, "
        f"{total_changed:,} total line(s) changed.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
