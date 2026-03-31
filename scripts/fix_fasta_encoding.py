#!/usr/bin/env python3
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

    ~ >...Institut de Pathologie et G\u00e9n\u00e9tique (IPG)|...
    + >...Institut de Pathologie et Génétique (IPG)|...

French (France)::

    ~ >...Cerballiance Montlh\u00e9ry|...
    + >...Cerballiance Montlhéry|...

Equivalent manual pipeline
--------------------------
The following shell commands achieve the same result but require a Java
runtime (``native2ascii``) and may fail on mixed-encoding files::

    # Convert \uXXXX escapes to real Unicode (Java native2ascii):
    native2ascii -encoding UTF-8 -reverse input.fasta output.fasta

    # Optionally transliterate remaining non-ASCII to ASCII:
    unidecode output.fasta > output.fasta.tmp && mv output.fasta.tmp output.fasta

This script only performs the lossless normalisation step (equivalent to
``native2ascii -reverse`` but in pure Python with no Java dependency), and
additionally handles raw Latin-1 bytes mixed into the same line.  The lossy
ASCII transliteration via ``unidecode`` is not applied here; run it manually
on the output if ASCII-only downstream tools require it.

Usage examples
--------------
::

    # Preview changes without touching any file (safe, recommended first step):
    fix_fasta_encoding.py --dry-run spikenuc1207.fasta

    # Preview with counts only (no per-line diff):
    fix_fasta_encoding.py --stats-only spikenuc1207.fasta

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
"""

import argparse
import os
import re
import shutil
import sys
import tempfile
import time
from collections import Counter


VERSION = "202603312000"

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
    return _unescape_unicode(s)


def _process_file(path: str, dry_run: bool, overwrite: bool,
                  stats_only: bool, verbose: bool, progress: bool) -> int:
    """Normalize *path* in-place.  Returns the number of lines changed."""
    orig_path = path + ".orig"

    if os.path.exists(orig_path) and not overwrite:
        print(
            f"  Skipping {path}: backup {os.path.basename(orig_path)} already "
            "exists (use --overwrite to force).",
            file=sys.stderr,
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
                mix_counter[frozenset(enc_types)] += 1

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
        print(f"  {path}: no changes needed.", file=sys.stderr)
        return 0

    elapsed = time.monotonic() - t_start
    speed_str = _fmt_speed(total_bytes / elapsed) if elapsed > 1e-6 else ""
    action = "Would write" if dry_run else "Wrote"
    print(
        f"  {action} {path}: {changed:,} line(s) changed"
        + (f", {speed_str}" if speed_str else ""),
        file=sys.stderr,
    )

    # Encoding-mix summary — always printed when there are changes.
    _mix_labels = [
        (frozenset({'esc'}),            "\\uXXXX escapes only"),
        (frozenset({'lat1'}),           "Latin-1 raw bytes only"),
        (frozenset({'utf8'}),           "Valid UTF-8 multi-byte only"),
        (frozenset({'esc', 'utf8'}),    "\\uXXXX escapes  +  valid UTF-8 multi-byte"),
        (frozenset({'esc', 'lat1'}),    "\\uXXXX escapes  +  Latin-1 raw bytes"),
        (frozenset({'lat1', 'utf8'}),   "Latin-1 raw bytes  +  valid UTF-8 *** MIXED ***"),
        (frozenset({'esc','lat1','utf8'}),
                                        "All three (\\uXXXX + Latin-1 + UTF-8) *** MIXED ***"),
    ]
    if mix_counter:
        print("  Encoding-mix breakdown (per changed line):", file=sys.stderr)
        for key, label in _mix_labels:
            cnt = mix_counter.get(key, 0)
            if cnt:
                print(f"    {label:<52s}: {cnt:>8,}", file=sys.stderr)
        # Print any unexpected combinations not in the label table.
        for key, cnt in mix_counter.items():
            if not any(key == k for k, _ in _mix_labels):
                print(f"    (other: {sorted(key)}): {cnt:>8,}", file=sys.stderr)

    if dry_run:
        os.unlink(tmp_path)
        return changed

    # Temp file is already written; just rename into place.
    # The temp file lives in the same directory so both os.rename() calls
    # are on the same filesystem (guaranteed atomic on POSIX).
    # Backup original (or overwrite existing backup if --overwrite).
    os.rename(path, orig_path)
    print(f"  Renamed {os.path.basename(path)} -> {os.path.basename(orig_path)}",
          file=sys.stderr)
    os.rename(tmp_path, path)
    print(f"  Written clean UTF-8 -> {os.path.basename(path)}", file=sys.stderr)
    return changed


# ── entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    """Normalise FASTA header encoding for each input file."""
    opts = _parser.parse_args()

    total_files = 0
    total_changed = 0
    for path in opts.infiles:
        if not os.path.isfile(path):
            print(f"Warning: not a file, skipping: {path}", file=sys.stderr)
            continue
        print(f"{path}:", file=sys.stderr)
        n = _process_file(
            path,
            dry_run=opts.dry_run,
            overwrite=opts.overwrite,
            stats_only=opts.stats_only,
            verbose=opts.verbose,
            progress=opts.progress,
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
